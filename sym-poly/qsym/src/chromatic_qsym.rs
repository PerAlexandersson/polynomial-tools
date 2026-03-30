//! Chromatic quasisymmetric functions.
//!
//! X_G = Σ_{proper κ: V→ℤ₊} x_{κ(v_0)} ... x_{κ(v_{n-1})}
//!
//! In the QSym monomial basis, the coefficient of M_α counts proper
//! ordered set partitions of V into independent sets with sizes α.

use std::collections::BTreeMap;

use sym_poly_core::{Composition, Ring};

use crate::basis::QSymBasis;
use crate::qsym_function::QSymFunction;

/// Chromatic quasisymmetric function X_G in the monomial QSym basis.
///
/// Takes raw graph data (number of vertices, edge list).
pub fn chromatic_qsym<C: Ring>(
    n: usize,
    edges: &[(usize, usize)],
) -> QSymFunction<C> {
    if n == 0 {
        return QSymFunction::basis_element(QSymBasis::Monomial, Composition::empty());
    }

    let mut adj: Vec<Vec<bool>> = vec![vec![false; n]; n];
    for &(u, v) in edges {
        adj[u][v] = true;
        adj[v][u] = true;
    }

    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();

    // Enumerate set partitions using standard RGS (rgs[0]=0)
    // For each proper partition, contribute all permutations of class orderings.
    let mut rgs = vec![0usize; n];
    let mut max_class = vec![0usize; n]; // max_class[i] = max(rgs[0..=i])

    loop {
        process_partition(n, &rgs, &adj, &mut terms);

        // Advance to next RGS
        if !next_rgs(&mut rgs, &mut max_class, n) {
            break;
        }
    }

    QSymFunction::from_terms(QSymBasis::Monomial, terms)
}

/// Check if coloring is proper. If so, record contributions for all
/// permutations of class orderings.
fn process_partition<C: Ring>(
    n: usize,
    rgs: &[usize],
    adj: &[Vec<bool>],
    terms: &mut BTreeMap<Composition, C>,
) {
    // Check properness
    for u in 0..n {
        for v in u + 1..n {
            if adj[u][v] && rgs[u] == rgs[v] {
                return;
            }
        }
    }

    // Compute class sizes in RGS order (= order of first occurrence)
    let num_classes = rgs.iter().max().unwrap_or(&0) + 1;
    let mut sizes = vec![0u32; num_classes];
    for &c in rgs {
        sizes[c] += 1;
    }

    // For each of the k! orderings of classes, record the size composition.
    // Classes are distinguishable (by their vertex sets), so even classes with
    // equal sizes give distinct orderings that contribute independently to M_α.
    let k = sizes.len();
    let class_perms = all_index_permutations(k);
    for perm in &class_perms {
        let comp_parts: Vec<u32> = perm.iter().map(|&i| sizes[i]).collect();
        let comp = Composition::new(comp_parts);
        let entry = terms.entry(comp).or_insert_with(C::zero);
        *entry = entry.clone() + C::one();
    }
}

/// Generate all permutations of {0, ..., k-1}.
fn all_index_permutations(k: usize) -> Vec<Vec<usize>> {
    if k == 0 {
        return vec![vec![]];
    }
    let mut result = Vec::new();
    let mut perm: Vec<usize> = (0..k).collect();
    loop {
        result.push(perm.clone());
        // Next permutation in lexicographic order
        let mut i = k - 1;
        while i > 0 && perm[i - 1] >= perm[i] {
            i -= 1;
        }
        if i == 0 {
            break;
        }
        let mut j = k - 1;
        while perm[j] <= perm[i - 1] {
            j -= 1;
        }
        perm.swap(i - 1, j);
        perm[i..].reverse();
    }
    result
}

/// Advance to next standard restricted growth string.
/// Returns false when exhausted.
fn next_rgs(rgs: &mut [usize], max_class: &mut [usize], n: usize) -> bool {
    // Find rightmost position that can be incremented.
    // rgs[0] is always 0 (standard RGS constraint), so never increment it.
    let mut i = n - 1;
    loop {
        if i == 0 {
            return false; // rgs[0] must stay 0
        }
        let max_prev = max_class[i - 1];
        if rgs[i] < max_prev + 1 {
            rgs[i] += 1;
            max_class[i] = max_prev.max(rgs[i]);
            // Reset positions to the right
            for j in i + 1..n {
                rgs[j] = 0;
                max_class[j] = max_class[j - 1];
            }
            return true;
        }
        i -= 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_graph_1v() {
        let f: QSymFunction<i64> = chromatic_qsym(1, &[]);
        assert_eq!(f.coefficient(&Composition::new(vec![1])), 1);
    }

    #[test]
    fn test_edge_graph() {
        // K_2: one set partition {{0},{1}} with sizes (1,1).
        // Permutations of (1,1): just (1,1) (both parts equal → 1 permutation).
        // Coefficient of M_(1,1) = 1.
        // Note: X_{K_2} = 2·m_{(1,1)} in Sym, and m_{(1,1)} = M_{(1,1)}
        // (since (1,1) has only one rearrangement). So coefficient = 1? Or 2?
        //
        // Actually: m_{(1,1)} = M_{(1,1)} since sort(1,1) = (1,1) and there's
        // only one composition sorting to (1,1). So X_{K_2} = 2·M_{(1,1)} in QSym.
        // But our set partition approach gives 1 partition × 1 permutation = 1.
        //
        // The discrepancy: the Sym convention counts ALL proper colorings
        // κ: V → ℤ₊. For K_2, there are 2 proper colorings per pair (i<j):
        // κ(0)=i,κ(1)=j and κ(0)=j,κ(1)=i. But M_{(1,1)} = Σ_{i<j} x_i x_j
        // captures BOTH since x_i·x_j = x_j·x_i. So each (i,j) pair contributes
        // once to M_{(1,1)} via Σ_κ x_{κ(0)}·x_{κ(1)}, giving coefficient 2.
        //
        // Our set partition approach: 1 partition with ordering (1,1).
        // We need to multiply by the number of surjections from the classes
        // to actual colors, which for k classes is 1 (implicit in M_α).
        // But different orderings of the SAME set partition give different
        // colorings. For {{0},{1}} there are 2!/1! = 2 orderings... wait no,
        // both classes have size 1, so there's 2!/1!1! = 1 distinct permutation
        // (1,1) = (1,1).
        //
        // The issue is subtler. Each set partition into k classes yields
        // exactly ONE M_α contribution (summed over all c_1<...<c_k),
        // BUT different orderings of the classes give different κ values
        // that contribute to DIFFERENT terms of M_α.
        //
        // For K_2 with partition {{0},{1}}:
        // - Ordering ({0},{1}): for any c_1<c_2, κ(0)=c_1, κ(1)=c_2 → x_{c_1} x_{c_2}
        // - Ordering ({1},{0}): for any c_1<c_2, κ(1)=c_1, κ(0)=c_2 → x_{c_2} x_{c_1} = x_{c_1} x_{c_2}
        // Both orderings contribute to the SAME monomial! So the coefficient is 2.
        //
        // But with distinct permutations of sizes: (1,1) has only 1 permutation.
        // We're missing a factor. The factor is k!/aut, where aut accounts for
        // equal-sized classes that produce the same composition under permutation.
        //
        // For classes of sizes (1,1), there are 2! = 2 orderings, but only
        // 2!/1! = 1 DISTINCT composition (since both give (1,1)).
        // The number of orderings giving (1,1) is 2! / 1!1! = 1... wait no.
        // There are 2 orderings total: ({0},{1}) and ({1},{0}). Both give
        // sizes (1,1). So the coefficient should be 2, not 1.
        //
        // The fix: for each set partition, the coefficient of M_α is the
        // number of orderings of classes that produce sizes α.
        // For sizes (s_1,...,s_k) and a composition α, the count is
        // k! / (product of factorials of multiplicity of each distinct size).
        // Wait, that's the total number of permutations. For composition (1,1)
        // from sizes (1,1), ALL 2! = 2 orderings give (1,1). So coefficient = 2.
        let f: QSymFunction<i64> = chromatic_qsym(2, &[(0, 1)]);
        assert_eq!(f.coefficient(&Composition::new(vec![1, 1])), 2);
        assert_eq!(f.terms().len(), 1);
    }

    #[test]
    fn test_p3() {
        // Path P3. Proper set partitions:
        // {{0,2},{1}}: sizes (2,1). Orderings: (2,1) and (1,2). Both give distinct comps.
        //   2 orderings total, coefficient of M_(2,1) = 1, M_(1,2) = 1.
        // {{0},{1},{2}}: sizes (1,1,1). 3! = 6 orderings, all give (1,1,1).
        //   Coefficient of M_(1,1,1) = 6.
        let f: QSymFunction<i64> = chromatic_qsym(3, &[(0, 1), (1, 2)]);
        assert_eq!(f.coefficient(&Composition::new(vec![2, 1])), 1);
        assert_eq!(f.coefficient(&Composition::new(vec![1, 2])), 1);
        assert_eq!(f.coefficient(&Composition::new(vec![1, 1, 1])), 6);
        assert_eq!(f.terms().len(), 3);
    }

    #[test]
    fn test_k3() {
        // K3: only partition {{0},{1},{2}}, sizes (1,1,1), 3!=6 orderings.
        let f: QSymFunction<i64> = chromatic_qsym(3, &[(0, 1), (0, 2), (1, 2)]);
        assert_eq!(f.coefficient(&Composition::new(vec![1, 1, 1])), 6);
        assert_eq!(f.terms().len(), 1);
    }

    #[test]
    fn test_is_symmetric() {
        // X_G is symmetric: rearrangements of α should have equal coefficients.
        let f: QSymFunction<i64> = chromatic_qsym(3, &[(0, 1), (1, 2)]);
        assert_eq!(
            f.coefficient(&Composition::new(vec![2, 1])),
            f.coefficient(&Composition::new(vec![1, 2])),
        );
    }
}
