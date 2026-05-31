//! Chromatic symmetric functions of graphs.
//!
//! The chromatic symmetric function X_G of a graph G on n vertices is the
//! symmetric function
//!
//! ```text
//! X_G = Σ_{proper colorings κ} x_{κ(1)} x_{κ(2)} ··· x_{κ(n)}
//! ```
//!
//! where the sum is over all proper colorings κ: V(G) → ℤ₊ (no two adjacent
//! vertices receive the same color). This was introduced by Stanley (1995).
//!
//! The **q-chromatic quasisymmetric function** (Shareshian-Wachs, 2016)
//! refines this with a parameter q tracking the number of ascents:
//!
//! ```text
//! X_G(q) = Σ_{proper κ} q^{asc(κ)} x_{κ(1)} x_{κ(2)} ··· x_{κ(n)}
//! ```
//!
//! where asc(κ) = |{(i,j) ∈ E(G) : i < j, κ(i) < κ(j)}|.
//!
//! # Examples
//!
//! ```
//! use combinatoric_core::graph::Graph;
//! use combinatoric_core::partition::Partition;
//! use sym_poly_sym::chromatic_symmetric;
//!
//! // X_{P3} = m_{2,1} + 6 m_{1,1,1}
//! let p3 = Graph::path(3);
//! let x = chromatic_symmetric::<i64>(&p3);
//! let m = x.to_monomial_basis();
//! assert_eq!(m.coefficient(&Partition::new(vec![2, 1])), 1);
//! assert_eq!(m.coefficient(&Partition::new(vec![1, 1, 1])), 6);
//! ```

use std::collections::BTreeMap;

use combinatoric_core::{Graph, Partition, Ring};

use crate::{Basis, SymmetricFunction};

// ---------------------------------------------------------------------------
// Chromatic symmetric function
// ---------------------------------------------------------------------------

/// Compute the chromatic symmetric function X_G in the monomial basis.
///
/// For a graph G on n vertices, X_G = Σ c_λ m_λ where c_λ counts the number
/// of proper colorings of G with color type λ (i.e., the sorted tuple of
/// color multiplicities equals λ), multiplied by the number of distinct
/// monomial assignments for that type.
///
/// The coefficient of m_λ equals the number of proper colorings κ: V → ℤ₊
/// whose sorted color-frequency vector is λ, divided by the number of
/// permutations of equal parts (the automorphism factor), then multiplied
/// back — but since we enumerate all labeled colorings with n colors, the
/// coefficient is simply the count of proper colorings with type λ times
/// the number of ways to assign the actual color values (which is
/// n! / (m_1! m_2! ··· aut(λ)) — but we take the direct approach instead.
pub fn chromatic_symmetric<C: Ring>(g: &Graph) -> SymmetricFunction<C> {
    let n = g.num_vertices();

    // For each partition λ of n, generate all distinct permutations of a
    // coloring vector with color frequencies λ, and count proper ones.
    // This is much faster than iterating over all k^n colorings.
    let partitions = Partition::all_of_size(n as u32);
    let mut terms: BTreeMap<Partition, C> = BTreeMap::new();

    for lambda in &partitions {
        let count = count_proper_colorings_of_type(g, lambda);
        if count != 0 {
            terms.insert(lambda.clone(), C::from_i64(count as i64));
        }
    }

    SymmetricFunction::from_terms(Basis::Monomial, terms)
}

/// Compute the ordered non-proper coloring symmetric function `f_(G, Ω)`.
///
/// Fix a total order `Ω` on the edges of `G`, supplied as `edge_order`.
/// For each non-proper coloring `κ`, let `e = ij` be the first monochromatic
/// edge in that order and set `ν(κ) = κ(i) = κ(j)`. Then
///
/// ```text
/// f_(G, Ω)(x) = Σ_κ (1 / x_{ν(κ)}) ∏_{v ∈ V} x_{κ(v)}
/// ```
///
/// where the sum is over all non-proper colorings.
///
/// This implementation uses the equivalent decomposition
///
/// ```text
/// f_(G, Ω) = Σ_{e ∈ E} X_{G_e},
/// ```
///
/// where `G_e` is obtained by taking the graph of edges earlier than `e` in
/// `Ω` and then identifying the endpoints of `e`.
pub fn first_bad_edge_symmetric<C: Ring>(
    g: &Graph,
    edge_order: &[(usize, usize)],
) -> SymmetricFunction<C> {
    let mut normalized_order: Vec<_> = edge_order.iter().copied().map(normalize_edge).collect();
    normalized_order.sort_unstable();
    normalized_order.dedup();
    assert_eq!(
        normalized_order.len(),
        edge_order.len(),
        "edge_order must not contain duplicates"
    );
    let mut graph_edges = g.edges().to_vec();
    graph_edges.sort_unstable();
    assert_eq!(
        normalized_order.as_slice(),
        graph_edges.as_slice(),
        "edge_order must list each edge of the graph exactly once"
    );

    let mut terms: BTreeMap<Partition, C> = BTreeMap::new();

    for idx in 0..edge_order.len() {
        let prefix_edges: Vec<_> = edge_order[..idx]
            .iter()
            .copied()
            .map(normalize_edge)
            .collect();
        let (u, v) = normalize_edge(edge_order[idx]);
        let prefix_graph = Graph::new(g.num_vertices(), &prefix_edges);
        let ge = prefix_graph.identify_vertices(u, v);
        let contribution = chromatic_symmetric::<C>(&ge);
        for (partition, coeff) in contribution.terms() {
            let entry = terms.entry(partition.clone()).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone();
        }
    }

    SymmetricFunction::from_terms(Basis::Monomial, terms)
}

fn normalize_edge((u, v): (usize, usize)) -> (usize, usize) {
    if u < v {
        (u, v)
    } else {
        (v, u)
    }
}

/// Count proper colorings whose sorted color-frequency vector equals λ,
/// divided by the number of monomials of type λ (giving the m_λ coefficient).
///
/// Strategy: build a coloring vector with frequencies λ (e.g., [1,1,2,3,3,3]
/// for λ = (3,2,1)), generate all distinct permutations, and count those that
/// are proper colorings of g. The count per monomial is this divided by the
/// number of distinct rearrangements of λ-with-zeros, but since we use exactly
/// ℓ(λ) colors, each valid coloring maps to exactly one monomial, so the
/// coefficient = (number of proper permutations) / (number of distinct
/// rearrangements of λ padded to ℓ(λ) colors).
///
/// Actually, with ℓ(λ) colors (no unused colors), every rearrangement of the
/// frequency vector uses all colors, so the number of rearrangements equals
/// ℓ(λ)! / ∏ m_i! where m_i is the multiplicity of part i in λ. But we need
/// to account for ALL possible number of colors (not just ℓ(λ)), since a
/// coloring with unused colors also contributes.
///
/// The simplest correct approach: use the multiset permutation enumeration.
/// For a fixed composition α (like (3,2,1)), the coloring vector [1,1,1,2,2,3]
/// has n!/∏(α_i!) distinct permutations. We test each for properness. The
/// coefficient of m_λ is the count of proper ones divided by ℓ(λ)!/∏ m_i!
/// (the number of relabelings of colors that preserve the type).
fn count_proper_colorings_of_type(g: &Graph, lambda: &Partition) -> i64 {
    let n = g.num_vertices();
    let parts = lambda.parts();

    // Build the canonical coloring vector: color 0 appears λ_1 times,
    // color 1 appears λ_2 times, etc.
    let mut base_coloring = Vec::with_capacity(n);
    for (color, &freq) in parts.iter().enumerate() {
        for _ in 0..freq {
            base_coloring.push(color);
        }
    }
    assert_eq!(base_coloring.len(), n);

    // Generate all distinct permutations (multiset permutations) and count
    // those that are proper colorings. Each proper permutation maps to
    // exactly one monomial x₀^{λ₁} x₁^{λ₂} ···, and the coefficient of
    // m_λ is this count (no division needed — we fix one composition).
    let mut count = 0i64;
    let mut perm = base_coloring;
    perm.sort();
    loop {
        let proper = g.edges().iter().all(|&(u, v)| perm[u] != perm[v]);
        if proper {
            count += 1;
        }
        if !next_multiset_perm(&mut perm) {
            break;
        }
    }

    count
}

/// Next lexicographic permutation of a slice (handles duplicates correctly).
fn next_multiset_perm(perm: &mut [usize]) -> bool {
    let n = perm.len();
    if n <= 1 {
        return false;
    }
    // Find largest i such that perm[i] < perm[i+1]
    let mut i = n - 2;
    loop {
        if perm[i] < perm[i + 1] {
            break;
        }
        if i == 0 {
            return false;
        }
        i -= 1;
    }
    // Find largest j such that perm[i] < perm[j]
    let mut j = n - 1;
    while perm[j] <= perm[i] {
        j -= 1;
    }
    perm.swap(i, j);
    perm[i + 1..].reverse();
    true
}

// ---------------------------------------------------------------------------
// q-Chromatic quasisymmetric function
// ---------------------------------------------------------------------------

/// Compute the q-chromatic symmetric function X_G(q) in the monomial basis.
///
/// The coefficient of q^a · m_λ counts proper colorings with type λ and
/// exactly `a` ascending edges (edges (i,j) with i < j and κ(i) < κ(j)).
///
/// Returns a map from (partition, ascent_count) to coefficient.
/// This is represented as a `BTreeMap<Partition, Vec<i64>>` where the Vec
/// is indexed by ascent count.
pub fn q_chromatic_symmetric(g: &Graph) -> BTreeMap<Partition, Vec<i64>> {
    q_chromatic_symmetric_with_ascent_edges(g.num_vertices(), g.edges(), g.edges())
}

/// Compute the q-chromatic symmetric function with a separate ascent
/// orientation.
///
/// The `proper_edges` determine which colorings are proper, as an undirected
/// graph.  The `ascent_edges` are directed edges `(u, v)` contributing one
/// ascent when `color(u) < color(v)`.  This covers circular unit arc digraphs,
/// where the proper-coloring graph is the underlying simple graph but the
/// ascent statistic uses the circular orientation.
pub fn q_chromatic_symmetric_with_ascent_edges(
    n: usize,
    proper_edges: &[(usize, usize)],
    ascent_edges: &[(usize, usize)],
) -> BTreeMap<Partition, Vec<i64>> {
    assert!(
        ascent_edges.iter().all(|&(u, v)| u < n && v < n && u != v),
        "ascent edge endpoint out of range"
    );

    let graph = Graph::new(n, proper_edges);
    let partitions = Partition::all_of_size(n as u32);
    let mut result: BTreeMap<Partition, Vec<i64>> = BTreeMap::new();

    for lambda in &partitions {
        let q_counts =
            count_proper_colorings_of_type_q_with_ascent_edges(&graph, ascent_edges, lambda);
        if q_counts.iter().any(|&c| c != 0) {
            result.insert(lambda.clone(), q_counts);
        }
    }

    result
}

/// Circular unit interval Frobenius target for the naive dot-action test.
///
/// This returns the coefficient of each `q^d` in
/// `ω X_{\Gamma_a}(x;q)`, where proper colorings use the underlying circular
/// unit interval graph and ascents use the circular edge orientation.
pub fn circular_area_dot_frobenius_target(
    area: &[u8],
) -> Option<BTreeMap<u32, SymmetricFunction<i64>>> {
    let directed_edges = Graph::circular_unit_interval_directed_edges(area)?;
    let graph = Graph::circular_unit_interval(area)?;
    let q_chromatic = q_chromatic_symmetric_with_ascent_edges(
        graph.num_vertices(),
        graph.edges(),
        &directed_edges,
    );
    Some(omega_by_q_degree(q_chromatic))
}

/// Shareshian--Wachs Frobenius target for a Hessenberg area sequence.
///
/// For the unit-interval graph `G` encoded by `area`, this returns the
/// coefficient of each `q^d` in `ω X_G(q)`. Under the Shareshian--Wachs
/// theorem, these are the Frobenius characteristics of the graded dot-action
/// pieces of the corresponding regular semisimple Hessenberg variety.
pub fn hessenberg_area_dot_frobenius_target(
    area: &[u8],
) -> Option<BTreeMap<u32, SymmetricFunction<i64>>> {
    if !is_area_sequence(area) {
        return None;
    }

    let graph = Graph::unit_interval(area);
    let q_chromatic = q_chromatic_symmetric(&graph);
    Some(omega_by_q_degree(q_chromatic))
}

fn omega_by_q_degree(
    q_chromatic: BTreeMap<Partition, Vec<i64>>,
) -> BTreeMap<u32, SymmetricFunction<i64>> {
    let mut by_degree: BTreeMap<u32, BTreeMap<Partition, i64>> = BTreeMap::new();
    for (partition, coeffs) in q_chromatic {
        for (degree, coeff) in coeffs.into_iter().enumerate() {
            if coeff == 0 {
                continue;
            }
            by_degree
                .entry(degree as u32)
                .or_default()
                .insert(partition.clone(), coeff);
        }
    }

    by_degree
        .into_iter()
        .map(|(degree, terms)| {
            let chromatic_degree = SymmetricFunction::from_terms(Basis::Monomial, terms);
            (degree, chromatic_degree.omega_involution())
        })
        .collect()
}

fn is_area_sequence(area: &[u8]) -> bool {
    area.iter().enumerate().all(|(i, &v)| v as usize <= i)
        && area
            .windows(2)
            .all(|w| usize::from(w[1]) <= usize::from(w[0]) + 1)
}

fn count_proper_colorings_of_type_q_with_ascent_edges(
    g: &Graph,
    ascent_edges: &[(usize, usize)],
    lambda: &Partition,
) -> Vec<i64> {
    let n = g.num_vertices();
    let parts = lambda.parts();
    let max_ascents = ascent_edges.len();

    let mut base_coloring = Vec::with_capacity(n);
    for (color, &freq) in parts.iter().enumerate() {
        for _ in 0..freq {
            base_coloring.push(color);
        }
    }

    let mut counts = vec![0i64; max_ascents + 1];
    let mut perm = base_coloring;
    perm.sort();

    loop {
        let proper = g.edges().iter().all(|&(u, v)| perm[u] != perm[v]);
        if proper {
            let ascents = ascent_edges
                .iter()
                .filter(|&&(u, v)| perm[u] < perm[v])
                .count();
            counts[ascents] += 1;
        }
        if !next_multiset_perm(&mut perm) {
            break;
        }
    }

    while counts.len() > 1 && *counts.last().unwrap() == 0 {
        counts.pop();
    }

    counts
}

/// Format the q-chromatic symmetric function for display.
pub fn format_q_chromatic(qcs: &BTreeMap<Partition, Vec<i64>>) -> String {
    let mut parts: Vec<String> = Vec::new();
    for (lambda, coeffs) in qcs {
        let poly_parts: Vec<String> = coeffs
            .iter()
            .enumerate()
            .filter(|(_, &c)| c != 0)
            .map(|(a, &c)| {
                if a == 0 {
                    format!("{}", c)
                } else if a == 1 {
                    if c == 1 {
                        "q".to_string()
                    } else {
                        format!("{}q", c)
                    }
                } else if c == 1 {
                    format!("q^{}", a)
                } else {
                    format!("{}q^{}", c, a)
                }
            })
            .collect();
        if !poly_parts.is_empty() {
            parts.push(format!(
                "({}) m_{{{}}}",
                poly_parts.join(" + "),
                lambda.display()
            ));
        }
    }
    parts.join(" + ")
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- Chromatic symmetric function tests (Mathematica-verified) --

    #[test]
    fn test_edge_graph() {
        // K2: colorings (0,1) and (1,0), both → x₁x₂. coeff of m_{1,1} = 2.
        let g = Graph::new(2, &[(0, 1)]);
        let x = chromatic_symmetric::<i64>(&g);
        let m = x.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![1, 1])), 2);
    }

    #[test]
    fn test_p3_monomial() {
        // X_{P3}: direct polynomial = x₁²x₂ + ... + 6x₁x₂x₃
        // = 1·m_{2,1} + 6·m_{1,1,1}
        let x = chromatic_symmetric::<i64>(&Graph::path(3));
        let m = x.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![2, 1])), 1);
        assert_eq!(m.coefficient(&Partition::new(vec![1, 1, 1])), 6);
        assert_eq!(m.coefficient(&Partition::new(vec![3])), 0);
    }

    #[test]
    fn test_p3_trivial_spec() {
        // X_{P3}(1,1,1) should equal chi_{P3}(3) = 3·2·2 = 12
        let x = chromatic_symmetric::<i64>(&Graph::path(3));
        assert_eq!(x.trivial_specialization(3), 12);
    }

    #[test]
    fn test_k3_monomial() {
        // K3: 6 proper 3-colorings, all type {1,1,1}, one monomial x₁x₂x₃.
        // coeff = 6.
        let x = chromatic_symmetric::<i64>(&Graph::complete(3));
        let m = x.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![1, 1, 1])), 6);
        assert_eq!(m.coefficient(&Partition::new(vec![2, 1])), 0);
    }

    #[test]
    fn test_k4_monomial() {
        // K4: 4! = 24 colorings, all type {1,1,1,1}. coeff = 24.
        let x = chromatic_symmetric::<i64>(&Graph::complete(4));
        let m = x.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![1, 1, 1, 1])), 24);
    }

    #[test]
    fn test_p4_trivial_spec() {
        // chi_{P4}(k) = k(k-1)^3
        // chi(4) = 4·27 = 108
        let x = chromatic_symmetric::<i64>(&Graph::path(4));
        assert_eq!(x.trivial_specialization(4), 108);
        assert_eq!(x.trivial_specialization(3), 3 * 8); // 3·2³ = 24
    }

    #[test]
    fn test_p4_schur_nonneg() {
        // Paths are e-positive (Stanley-Stembridge for unit interval graphs),
        // so all Schur coefficients should be non-negative.
        let x = chromatic_symmetric::<i64>(&Graph::path(4));
        let s = x.to_schur_basis();
        for (_lam, &c) in s.terms() {
            assert!(c >= 0, "Negative Schur coefficient in X_{{P4}}");
        }
    }

    #[test]
    fn test_empty_graph_chromatic() {
        // Empty graph on 3: all 3^3 = 27 colorings are proper (k=3 colors).
        // X(1,1,1) = chi(3) = 27.
        let g = Graph::empty(3);
        let x = chromatic_symmetric::<i64>(&g);
        assert_eq!(x.trivial_specialization(3), 27);
    }

    #[test]
    fn test_first_bad_edge_single_edge() {
        let g = Graph::new(2, &[(0, 1)]);
        let f = first_bad_edge_symmetric::<i64>(&g, &[(0, 1)]);
        let m = f.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![1])), 1);
        assert_eq!(m.terms().len(), 1);
    }

    #[test]
    fn test_first_bad_edge_path3() {
        let g = Graph::path(3);
        let f = first_bad_edge_symmetric::<i64>(&g, &[(0, 1), (1, 2)]);
        let m = f.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![2])), 1);
        assert_eq!(m.coefficient(&Partition::new(vec![1, 1])), 4);
        assert_eq!(m.terms().len(), 2);
    }

    #[test]
    fn test_first_bad_edge_accepts_unsorted_graph_edges() {
        let g = Graph::new(3, &[(1, 2), (0, 1)]);
        let f = first_bad_edge_symmetric::<i64>(&g, &[(0, 1), (1, 2)]);
        let m = f.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![2])), 1);
        assert_eq!(m.coefficient(&Partition::new(vec![1, 1])), 4);
        assert_eq!(m.terms().len(), 2);
    }

    // -- q-Chromatic tests (Mathematica-verified) --

    #[test]
    fn test_q_chromatic_p3() {
        // P3 type {2,1}: all 6 colorings have 1 ascent, rearr=6 → q·m_{2,1}
        // P3 type {1,1,1}: 6 colorings split as 3 with 0 asc + 3 with 1 asc,
        // rearr=1 → (3+3q)·m_{1,1,1} ... wait no, rearr for {1,1,1} with k=3
        // is 3!/(3!) = 1, so we DON'T divide. 6 colorings total.
        // Let me enumerate type {1,1,1}:
        // (0,1,2): asc(0,1)=0<1 yes, asc(1,2)=1<2 yes → 2 asc
        // (0,2,1): asc(0,1)=0<2 yes, asc(1,2)=2>1 no → 1 asc
        // (1,0,2): asc(0,1)=1>0 no, asc(1,2)=0<2 yes → 1 asc
        // (1,2,0): asc(0,1)=1<2 yes, asc(1,2)=2>0 no → 1 asc
        // (2,0,1): asc(0,1)=2>0 no, asc(1,2)=0<1 yes → 1 asc
        // (2,1,0): asc(0,1)=2>1 no, asc(1,2)=1>0 no → 0 asc
        // Counts: asc=0→1, asc=1→4, asc=2→1. Divided by rearr=1: same.
        let qcs = q_chromatic_symmetric(&Graph::path(3));

        let m21 = qcs.get(&Partition::new(vec![2, 1])).unwrap();
        assert_eq!(m21.len(), 2);
        assert_eq!(m21[0], 0);
        assert_eq!(m21[1], 1);

        let m111 = qcs.get(&Partition::new(vec![1, 1, 1])).unwrap();
        assert_eq!(m111[0], 1); // 0 ascents
        assert_eq!(m111[1], 4); // 1 ascent
        assert_eq!(m111[2], 1); // 2 ascents
    }

    #[test]
    fn test_q_chromatic_p3_sum() {
        // Setting q=1 in qX should give the same as X.
        // m_{2,1}: 0+1 = 1. m_{1,1,1}: 1+4+1 = 6. Matches X_{P3}.
        let qcs = q_chromatic_symmetric(&Graph::path(3));
        let m21_total: i64 = qcs[&Partition::new(vec![2, 1])].iter().sum();
        let m111_total: i64 = qcs[&Partition::new(vec![1, 1, 1])].iter().sum();
        assert_eq!(m21_total, 1);
        assert_eq!(m111_total, 6);
    }

    #[test]
    fn test_hessenberg_area_dot_frobenius_target_edgeless() {
        let target = hessenberg_area_dot_frobenius_target(&[0, 0, 0]).unwrap();

        assert_eq!(target.keys().copied().collect::<Vec<_>>(), vec![0]);
        let schur = target[&0].to_schur_basis();
        assert_eq!(schur.coefficient(&Partition::new(vec![3])), 1);
        assert_eq!(schur.coefficient(&Partition::new(vec![2, 1])), 2);
        assert_eq!(schur.coefficient(&Partition::new(vec![1, 1, 1])), 1);
        assert_eq!(schur.terms().len(), 3);
    }

    #[test]
    fn test_hessenberg_area_dot_frobenius_target_complete_graph_k3() {
        let target = hessenberg_area_dot_frobenius_target(&[0, 1, 2]).unwrap();

        assert_eq!(target.keys().copied().collect::<Vec<_>>(), vec![0, 1, 2, 3]);
        let degree_zero = target[&0].to_schur_basis();
        let degree_one = target[&1].to_schur_basis();
        let degree_two = target[&2].to_schur_basis();
        let degree_three = target[&3].to_schur_basis();

        assert_eq!(degree_zero.coefficient(&Partition::new(vec![3])), 1);
        assert_eq!(degree_zero.terms().len(), 1);
        assert_eq!(degree_one.coefficient(&Partition::new(vec![3])), 2);
        assert_eq!(degree_one.terms().len(), 1);
        assert_eq!(degree_two.coefficient(&Partition::new(vec![3])), 2);
        assert_eq!(degree_two.terms().len(), 1);
        assert_eq!(degree_three.coefficient(&Partition::new(vec![3])), 1);
        assert_eq!(degree_three.terms().len(), 1);
    }

    #[test]
    fn test_hessenberg_area_dot_frobenius_target_rejects_invalid_area() {
        assert!(hessenberg_area_dot_frobenius_target(&[0, 2]).is_none());
    }

    #[test]
    fn test_circular_area_target_extends_hessenberg_target() {
        let area = [0, 1, 1];
        assert_eq!(
            circular_area_dot_frobenius_target(&area).unwrap(),
            hessenberg_area_dot_frobenius_target(&area).unwrap()
        );
    }

    #[test]
    fn test_circular_area_target_directed_cycle_s3() {
        let target = circular_area_dot_frobenius_target(&[1, 1, 1]).unwrap();

        assert_eq!(target.keys().copied().collect::<Vec<_>>(), vec![1, 2]);
        let degree_one = target[&1].to_schur_basis();
        let degree_two = target[&2].to_schur_basis();
        assert_eq!(degree_one.coefficient(&Partition::new(vec![3])), 3);
        assert_eq!(degree_one.terms().len(), 1);
        assert_eq!(degree_two.coefficient(&Partition::new(vec![3])), 3);
        assert_eq!(degree_two.terms().len(), 1);
    }

    // -- Edge cases --

    #[test]
    fn test_single_vertex_chromatic() {
        let x = chromatic_symmetric::<i64>(&Graph::empty(1));
        let m = x.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::new(vec![1])), 1);
    }

    #[test]
    fn test_zero_vertex_chromatic_is_one() {
        let x = chromatic_symmetric::<i64>(&Graph::empty(0));
        let m = x.to_monomial_basis();
        assert_eq!(m.coefficient(&Partition::empty()), 1);
        assert_eq!(m.terms().len(), 1);
        assert_eq!(m.trivial_specialization(5), 1);
    }

    #[test]
    fn test_zero_vertex_q_chromatic_is_one() {
        let qcs = q_chromatic_symmetric(&Graph::empty(0));

        assert_eq!(qcs, BTreeMap::from([(Partition::empty(), vec![1])]));
    }
}
