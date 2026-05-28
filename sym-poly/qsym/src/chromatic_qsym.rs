//! Chromatic quasisymmetric functions.
//!
//! X_G = Σ_{proper κ: V→ℤ₊} x_{κ(v_0)} ... x_{κ(v_{n-1})}
//!
//! In the QSym monomial basis, the coefficient of M_α counts proper
//! ordered set partitions of V into independent sets with sizes α.

use std::collections::BTreeMap;

use combinatoric_core::{ordered_set_partitions, Graph};
use sym_poly_core::{Composition, Ring, UnivariatePolynomial};

use crate::basis::QSymBasis;
use crate::qsym_function::QSymFunction;

/// Chromatic quasisymmetric function X_G in the monomial QSym basis.
///
/// Takes raw graph data (number of vertices, edge list).
pub fn chromatic_qsym<C: Ring>(n: usize, edges: &[(usize, usize)]) -> QSymFunction<C> {
    if n == 0 {
        return QSymFunction::basis_element(QSymBasis::Monomial, Composition::empty());
    }

    let mut adj: Vec<Vec<bool>> = vec![vec![false; n]; n];
    for &(u, v) in edges {
        adj[u][v] = true;
        adj[v][u] = true;
    }

    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();

    for partition in ordered_set_partitions(n) {
        process_ordered_partition(&partition.into_blocks(), &adj, &mut terms);
    }

    QSymFunction::from_terms(QSymBasis::Monomial, terms)
}

/// Shareshian--Wachs style asc-weighted chromatic quasisymmetric function.
///
/// The graph vertices are ordered by their labels `0 < 1 < ... < n-1`.  For a
/// proper coloring encoded as an ordered set partition `(B_1, ..., B_k)`, the
/// ascent statistic is
///
/// ```text
/// #{ {i,j} in E : i < j and block(i) < block(j) }.
/// ```
///
/// The coefficient of `M_alpha` is therefore a polynomial in `q`, where the
/// coefficient of `q^a` counts proper ordered set partitions of block-size
/// composition `alpha` with exactly `a` ascents.
pub fn chromatic_qsym_asc(
    n: usize,
    edges: &[(usize, usize)],
) -> QSymFunction<UnivariatePolynomial<i64>> {
    if n == 0 {
        return QSymFunction::basis_element(QSymBasis::Monomial, Composition::empty());
    }

    let normalized_edges = normalize_edges(n, edges);
    let mut adj: Vec<Vec<bool>> = vec![vec![false; n]; n];
    for &(u, v) in &normalized_edges {
        adj[u][v] = true;
        adj[v][u] = true;
    }

    let mut terms: BTreeMap<Composition, UnivariatePolynomial<i64>> = BTreeMap::new();

    for partition in ordered_set_partitions(n) {
        process_ordered_partition_with_ascents(
            &partition.into_blocks(),
            &adj,
            &normalized_edges,
            &mut terms,
        );
    }

    QSymFunction::from_terms(QSymBasis::Monomial, terms)
}

/// LLT-style all-coloring quasisymmetric function with ascent weights.
///
/// This uses the same ascent statistic as [`chromatic_qsym_asc`], but removes
/// the proper-coloring condition. For unit-interval graphs, the symmetric
/// projection is the graph/unicellular LLT polynomial with the same edge set.
pub fn coloring_qsym_asc(
    n: usize,
    edges: &[(usize, usize)],
) -> QSymFunction<UnivariatePolynomial<i64>> {
    if n == 0 {
        return QSymFunction::basis_element(QSymBasis::Monomial, Composition::empty());
    }

    let normalized_edges = normalize_edges(n, edges);
    let mut terms: BTreeMap<Composition, UnivariatePolynomial<i64>> = BTreeMap::new();

    for partition in ordered_set_partitions(n) {
        process_ordered_partition_all_colorings_with_ascents(
            &partition.into_blocks(),
            &normalized_edges,
            &mut terms,
        );
    }

    QSymFunction::from_terms(QSymBasis::Monomial, terms)
}

/// LLT-style all-coloring quasisymmetric function with directed ascent edges.
///
/// The edge `(u, v)` contributes one ascent when the color of `u` is smaller
/// than the color of `v`. This is the circular analogue of
/// [`coloring_qsym_asc`], where wrap edges must keep their circular
/// orientation.
pub fn coloring_qsym_asc_with_ascent_edges(
    n: usize,
    ascent_edges: &[(usize, usize)],
) -> QSymFunction<UnivariatePolynomial<i64>> {
    if n == 0 {
        return QSymFunction::basis_element(QSymBasis::Monomial, Composition::empty());
    }

    let ascent_edges = normalize_directed_edges(n, ascent_edges);
    let mut terms: BTreeMap<Composition, UnivariatePolynomial<i64>> = BTreeMap::new();

    for partition in ordered_set_partitions(n) {
        process_ordered_partition_all_colorings_with_ascents(
            &partition.into_blocks(),
            &ascent_edges,
            &mut terms,
        );
    }

    QSymFunction::from_terms(QSymBasis::Monomial, terms)
}

/// Circular LLT-style all-coloring quasisymmetric function.
pub fn circular_coloring_qsym_asc(area: &[u8]) -> Option<QSymFunction<UnivariatePolynomial<i64>>> {
    let ascent_edges = Graph::circular_unit_interval_directed_edges(area)?;
    Some(coloring_qsym_asc_with_ascent_edges(
        area.len(),
        &ascent_edges,
    ))
}

/// Check if the ordered partition is proper. If so, record its block-size
/// composition.
fn process_ordered_partition<C: Ring>(
    blocks: &[Vec<usize>],
    adj: &[Vec<bool>],
    terms: &mut BTreeMap<Composition, C>,
) {
    for block in blocks {
        for i in 0..block.len() {
            for j in i + 1..block.len() {
                if adj[block[i]][block[j]] {
                    return;
                }
            }
        }
    }

    let comp = Composition::new(blocks.iter().map(|block| block.len() as u32).collect());
    let entry = terms.entry(comp).or_insert_with(C::zero);
    *entry = entry.clone() + C::one();
}

fn process_ordered_partition_with_ascents(
    blocks: &[Vec<usize>],
    adj: &[Vec<bool>],
    edges: &[(usize, usize)],
    terms: &mut BTreeMap<Composition, UnivariatePolynomial<i64>>,
) {
    for block in blocks {
        for i in 0..block.len() {
            for j in i + 1..block.len() {
                if adj[block[i]][block[j]] {
                    return;
                }
            }
        }
    }

    let mut block_of = vec![0usize; adj.len()];
    for (block_idx, block) in blocks.iter().enumerate() {
        for &vertex in block {
            block_of[vertex] = block_idx;
        }
    }

    let ascents = edges
        .iter()
        .filter(|&&(u, v)| block_of[u] < block_of[v])
        .count();
    let comp = Composition::new(blocks.iter().map(|block| block.len() as u32).collect());
    let entry = terms.entry(comp).or_insert_with(UnivariatePolynomial::zero);
    *entry = entry.clone() + UnivariatePolynomial::monomial(ascents, 1);
}

fn process_ordered_partition_all_colorings_with_ascents(
    blocks: &[Vec<usize>],
    edges: &[(usize, usize)],
    terms: &mut BTreeMap<Composition, UnivariatePolynomial<i64>>,
) {
    let mut block_of = vec![0usize; blocks.iter().map(Vec::len).sum()];
    for (block_idx, block) in blocks.iter().enumerate() {
        for &vertex in block {
            block_of[vertex] = block_idx;
        }
    }

    let ascents = edges
        .iter()
        .filter(|&&(u, v)| block_of[u] < block_of[v])
        .count();
    let comp = Composition::new(blocks.iter().map(|block| block.len() as u32).collect());
    let entry = terms.entry(comp).or_insert_with(UnivariatePolynomial::zero);
    *entry = entry.clone() + UnivariatePolynomial::monomial(ascents, 1);
}

fn normalize_edges(n: usize, edges: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut normalized_edges = Vec::new();
    for &(u, v) in edges {
        assert!(u < n && v < n, "edge endpoint out of range");
        if u == v {
            continue;
        }
        normalized_edges.push(if u < v { (u, v) } else { (v, u) });
    }
    normalized_edges.sort_unstable();
    normalized_edges.dedup();
    normalized_edges
}

fn normalize_directed_edges(n: usize, edges: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut normalized_edges = Vec::new();
    for &(u, v) in edges {
        assert!(u < n && v < n, "edge endpoint out of range");
        if u != v {
            normalized_edges.push((u, v));
        }
    }
    normalized_edges.sort_unstable();
    normalized_edges.dedup();
    normalized_edges
}

#[cfg(test)]
mod tests {
    use super::*;
    use combinatoric_core::Graph;
    use std::collections::BTreeMap;
    use sym_poly_sym::{
        chromatic_symmetric, circular_unicellular_llt, hessenberg_area_dot_frobenius_target,
        unicellular_llt,
    };

    use crate::sym_qsym::symmetric_qsym_to_sym;

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

    #[test]
    fn test_chromatic_qsym_projects_to_chromatic_symmetric() {
        let g = Graph::path(3);
        let q = chromatic_qsym::<i64>(g.num_vertices(), g.edges());
        let from_q = symmetric_qsym_to_sym(&q).unwrap();
        let from_sym = chromatic_symmetric::<i64>(&g).to_monomial_basis();
        assert_eq!(from_q, from_sym);
    }

    fn evaluate_q_at_one(f: &QSymFunction<UnivariatePolynomial<i64>>) -> QSymFunction<i64> {
        let terms = f
            .terms()
            .iter()
            .map(|(comp, poly)| {
                let value = poly.coeffs().iter().copied().sum::<i64>();
                (comp.clone(), value)
            })
            .collect();
        QSymFunction::from_terms(f.basis(), terms)
    }

    fn omega_symmetric_projection_by_q_degree(
        f: &QSymFunction<UnivariatePolynomial<i64>>,
    ) -> BTreeMap<u32, sym_poly_sym::SymmetricFunction<i64>> {
        let max_degree = f
            .terms()
            .values()
            .map(|poly| poly.coeffs().len().saturating_sub(1))
            .max()
            .unwrap_or(0);
        let mut result = BTreeMap::new();

        for degree in 0..=max_degree {
            let terms = f
                .terms()
                .iter()
                .filter_map(|(composition, poly)| {
                    let coeff = poly.coeffs().get(degree).copied().unwrap_or(0);
                    if coeff == 0 {
                        None
                    } else {
                        Some((composition.clone(), coeff))
                    }
                })
                .collect();
            let qsym_degree = QSymFunction::from_terms(f.basis(), terms);
            if qsym_degree.is_zero() {
                continue;
            }
            let sym_degree = symmetric_qsym_to_sym(&qsym_degree)
                .expect("chromatic QSym should be symmetric for unit interval graphs");
            result.insert(degree as u32, sym_degree.omega_involution());
        }

        result
    }

    #[test]
    fn test_chromatic_qsym_asc_edge() {
        let f = chromatic_qsym_asc(2, &[(0, 1)]);
        let coeff = f.coefficient(&Composition::new(vec![1, 1]));
        assert_eq!(coeff, UnivariatePolynomial::new(vec![1, 1]));
        assert_eq!(f.terms().len(), 1);
    }

    #[test]
    fn test_coloring_qsym_asc_edge() {
        let f = coloring_qsym_asc(2, &[(0, 1)]);

        assert_eq!(
            f.coefficient(&Composition::new(vec![2])),
            UnivariatePolynomial::new(vec![1])
        );
        assert_eq!(
            f.coefficient(&Composition::new(vec![1, 1])),
            UnivariatePolynomial::new(vec![1, 1])
        );
        assert_eq!(f.terms().len(), 2);
    }

    #[test]
    fn test_chromatic_qsym_asc_specializes_at_one() {
        let g = Graph::path(3);
        let weighted = chromatic_qsym_asc(g.num_vertices(), g.edges());
        let at_one = evaluate_q_at_one(&weighted);
        let unweighted = chromatic_qsym::<i64>(g.num_vertices(), g.edges());
        assert_eq!(at_one, unweighted);
    }

    #[test]
    fn test_hessenberg_target_matches_omega_chromatic_qsym_asc() {
        let area = [0, 1, 1];
        let graph = Graph::unit_interval(&area);
        let qsym = chromatic_qsym_asc(graph.num_vertices(), graph.edges());
        let from_qsym = omega_symmetric_projection_by_q_degree(&qsym);
        let target = hessenberg_area_dot_frobenius_target(&area).unwrap();

        assert_eq!(from_qsym, target);
    }

    #[test]
    fn test_unit_interval_coloring_qsym_projects_to_unicellular_llt() {
        let area = [0, 1, 1];
        let graph = Graph::unit_interval(&area);
        let qsym = coloring_qsym_asc(graph.num_vertices(), graph.edges());
        let from_qsym = symmetric_qsym_to_sym(&qsym).unwrap();

        assert_eq!(from_qsym, unicellular_llt(&area));
    }

    #[test]
    fn test_circular_coloring_qsym_projects_to_circular_unicellular_llt() {
        let area = [1, 1, 1];
        let qsym = circular_coloring_qsym_asc(&area).unwrap();
        let from_qsym = symmetric_qsym_to_sym(&qsym).unwrap();

        assert_eq!(from_qsym, circular_unicellular_llt(&area).unwrap());
    }
}
