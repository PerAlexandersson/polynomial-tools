//! Stanley's (P,w)-partitions and their generating functions in QSym.
//!
//! Given a poset P on {0,...,n-1} (as cover relations) and a labeling w,
//! the (P,w)-partition generating function is:
//!
//!   Γ(P,w) = Σ_{σ ∈ L(P)} F_{Des(σ)}
//!
//! where L(P) is the set of linear extensions and Des(σ) is the descent
//! composition of σ.
//!
//! This module takes raw poset data (cover relations) so it doesn't depend
//! on the Poset type in combinatoric-core. A Poset-aware wrapper can be
//! added there.

use std::collections::{BTreeMap, BTreeSet};

use sym_poly_core::{Composition, Ring};

use crate::basis::QSymBasis;
use crate::qsym_function::QSymFunction;

/// Compute the (P,w)-partition generating function in the fundamental basis.
///
/// # Arguments
/// * `n` - number of elements (labeled 0..n-1)
/// * `covers` - cover relations (i, j) meaning i <_P j
///
/// For a naturally labeled poset (labels = identity), returns:
///   Γ(P) = Σ_{σ ∈ L(P)} F_{Des(σ)}
pub fn p_partition_generating_function<C: Ring>(
    n: usize,
    covers: &[(usize, usize)],
) -> QSymFunction<C> {
    if n == 0 {
        return QSymFunction::basis_element(QSymBasis::Fundamental, Composition::empty());
    }

    // Build adjacency lists
    let mut children: Vec<Vec<usize>> = vec![vec![]; n];
    let mut parents: Vec<Vec<usize>> = vec![vec![]; n];
    for &(a, b) in covers {
        children[a].push(b);
        parents[b].push(a);
    }

    // Enumerate linear extensions via backtracking
    let extensions = enumerate_linear_extensions(n, &children, &parents);

    // For each linear extension, compute descent composition and accumulate
    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();
    for sigma in &extensions {
        let des_comp = descent_composition(sigma, n);
        let entry = terms.entry(des_comp).or_insert_with(C::zero);
        *entry = entry.clone() + C::one();
    }

    QSymFunction::from_terms(QSymBasis::Fundamental, terms)
}

/// Enumerate all linear extensions of a poset.
///
/// A linear extension is a permutation σ of {0,...,n-1} such that
/// σ(i) appears before σ(j) whenever i <_P j. We return σ as a
/// sequence [σ(0), σ(1), ...] giving the order of elements.
fn enumerate_linear_extensions(
    n: usize,
    children: &[Vec<usize>],
    parents: &[Vec<usize>],
) -> Vec<Vec<usize>> {
    // Find minimal elements (those with no parents)
    let mut in_degree: Vec<usize> = parents.iter().map(|p| p.len()).collect();
    let mut available: BTreeSet<usize> = BTreeSet::new();
    for i in 0..n {
        if in_degree[i] == 0 {
            available.insert(i);
        }
    }

    let mut results = Vec::new();
    let mut current = Vec::new();
    backtrack_extensions(
        n,
        children,
        &mut in_degree,
        &mut available,
        &mut current,
        &mut results,
    );
    results
}

fn backtrack_extensions(
    n: usize,
    children: &[Vec<usize>],
    in_degree: &mut Vec<usize>,
    available: &mut BTreeSet<usize>,
    current: &mut Vec<usize>,
    results: &mut Vec<Vec<usize>>,
) {
    if current.len() == n {
        results.push(current.clone());
        return;
    }

    let choices: Vec<usize> = available.iter().copied().collect();
    for v in choices {
        available.remove(&v);
        current.push(v);

        // Update in-degrees for children of v
        let mut newly_available = Vec::new();
        for &child in &children[v] {
            in_degree[child] -= 1;
            if in_degree[child] == 0 {
                available.insert(child);
                newly_available.push(child);
            }
        }

        backtrack_extensions(n, children, in_degree, available, current, results);

        // Undo
        current.pop();
        available.insert(v);
        for &child in &children[v] {
            in_degree[child] += 1;
        }
        for child in &newly_available {
            available.remove(child);
        }
    }
}

/// Compute the descent composition of a permutation σ.
///
/// The descent set is {i : σ[i] > σ[i+1]} (0-indexed positions).
/// The descent composition records the lengths of ascending runs.
fn descent_composition(sigma: &[usize], n: usize) -> Composition {
    if n <= 1 {
        return Composition::new(vec![n as u32]);
    }
    let mut des_set = BTreeSet::new();
    for i in 0..n - 1 {
        if sigma[i] > sigma[i + 1] {
            des_set.insert((i + 1) as u32);
        }
    }
    Composition::from_descent_set(&des_set, n as u32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chain_3() {
        // Chain 0 < 1 < 2: only one linear extension [0, 1, 2], no descents
        // Γ = F_(3)
        let f: QSymFunction<i64> = p_partition_generating_function(3, &[(0, 1), (1, 2)]);
        assert_eq!(f.coefficient(&Composition::new(vec![3])), 1);
        assert_eq!(f.terms().len(), 1);
    }

    #[test]
    fn test_antichain_2() {
        // Antichain on 2 elements: linear extensions are [0,1] and [1,0]
        // [0,1]: no descents → F_(2)
        // [1,0]: descent at position 1 → F_(1,1)
        // Γ = F_(2) + F_(1,1)
        let f: QSymFunction<i64> = p_partition_generating_function(2, &[]);
        assert_eq!(f.coefficient(&Composition::new(vec![2])), 1);
        assert_eq!(f.coefficient(&Composition::new(vec![1, 1])), 1);
    }

    #[test]
    fn test_antichain_3() {
        // Antichain on 3: 3! = 6 linear extensions
        // Γ should have total coefficient sum = 6
        let f: QSymFunction<i64> = p_partition_generating_function(3, &[]);
        let total: i64 = f.terms().values().sum();
        assert_eq!(total, 6);
    }

    #[test]
    fn test_v_poset() {
        // V-poset: 0 < 2, 1 < 2 (two elements below one)
        // Linear extensions: [0,1,2] and [1,0,2]
        // [0,1,2]: no descents → F_(3)
        // [1,0,2]: descent at 1 → F_(1,2)
        let f: QSymFunction<i64> = p_partition_generating_function(3, &[(0, 2), (1, 2)]);
        assert_eq!(f.coefficient(&Composition::new(vec![3])), 1);
        assert_eq!(f.coefficient(&Composition::new(vec![1, 2])), 1);
        assert_eq!(f.terms().len(), 2);
    }

    #[test]
    fn test_p_partition_is_qsym() {
        // Verify the result is homogeneous of the right degree
        let f: QSymFunction<i64> = p_partition_generating_function(4, &[(0, 2), (1, 3)]);
        for (comp, _) in f.terms() {
            assert_eq!(comp.size(), 4, "all compositions should have degree 4");
        }
    }

    #[test]
    fn test_fence_3() {
        // Fence: 0 < 1, 1 > 2 (= 2 < 1)
        // Cover relations: (0,1), (2,1)
        // Linear extensions: [0,2,1] and [2,0,1]
        // [0,2,1]: 0<2, 2>1 → descent at 2 → F_(2,1)
        // [2,0,1]: 2>0, 0<1 → descent at 1 → F_(1,2)
        let f: QSymFunction<i64> = p_partition_generating_function(3, &[(0, 1), (2, 1)]);
        assert_eq!(f.coefficient(&Composition::new(vec![2, 1])), 1);
        assert_eq!(f.coefficient(&Composition::new(vec![1, 2])), 1);
        assert_eq!(f.terms().len(), 2);
    }
}
