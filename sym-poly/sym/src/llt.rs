//! Graph LLT polynomials and unit-interval specializations.
//!
//! This module currently implements the discrete graph-coloring model in the
//! monomial basis. Coefficients are polynomials in `q`.

use std::collections::BTreeMap;

use sym_poly_core::{Partition, Ring, UnivariatePolynomial};

use crate::kostka::sn_character;
use crate::{Basis, SymmetricFunction};

/// Compute the graph LLT polynomial in the monomial basis.
///
/// The coefficient of `q^a m_\lambda` counts colorings of type `\lambda`
/// with exactly `a` ascents along `q`-edges, subject to the optional
/// `strict_edges` (must be ascents) and `weak_edges` (must not be ascents).
pub fn graph_llt_symmetric(
    n: usize,
    attacking_edges: &[(usize, usize)],
    strict_edges: &[(usize, usize)],
    weak_edges: &[(usize, usize)],
) -> SymmetricFunction<UnivariatePolynomial<i64>> {
    if n == 0 {
        return SymmetricFunction::basis_element(Basis::Monomial, Partition::empty());
    }

    let attacking = normalize_edges(attacking_edges);
    let strict = normalize_edges(strict_edges);
    let weak = normalize_edges(weak_edges);
    let q_edges: Vec<_> = attacking
        .iter()
        .copied()
        .filter(|edge| !strict.contains(edge))
        .collect();

    let mut result_terms = BTreeMap::new();
    for lambda in Partition::all_of_size(n as u32) {
        let coeff = count_llt_colorings_of_type(n, &lambda, &q_edges, &strict, &weak);
        if !coeff.is_zero() {
            result_terms.insert(lambda, coeff);
        }
    }

    SymmetricFunction::from_terms(Basis::Monomial, result_terms)
}

/// The unicellular LLT polynomial attached to a unit-interval area sequence.
pub fn unicellular_llt(area: &[u8]) -> SymmetricFunction<UnivariatePolynomial<i64>> {
    let edges = unit_interval_edges(area);
    graph_llt_symmetric(area.len(), &edges, &[], &[])
}

/// Degree-wise Schur expansion of the unicellular LLT Frobenius target.
///
/// The output maps a `q`-degree to the corresponding Schur-positive symmetric
/// function. It is the representation-theoretic target for the LLT/twin
/// manifold analogue of the Hessenberg dot-action story.
pub fn unicellular_llt_frobenius_target(
    area: &[u8],
) -> Option<BTreeMap<u32, SymmetricFunction<i64>>> {
    if !is_area_sequence(area) {
        return None;
    }

    let schur = unicellular_llt(area).to_schur_basis();
    let mut by_degree: BTreeMap<u32, BTreeMap<Partition, i64>> = BTreeMap::new();
    for (partition, coefficient) in schur.terms() {
        for (degree, &multiplicity) in coefficient.coeffs().iter().enumerate() {
            if multiplicity == 0 {
                continue;
            }
            by_degree
                .entry(degree as u32)
                .or_default()
                .insert(partition.clone(), multiplicity);
        }
    }

    Some(
        by_degree
            .into_iter()
            .map(|(degree, terms)| (degree, SymmetricFunction::from_terms(Basis::Schur, terms)))
            .collect(),
    )
}

/// Graded character values of the abstract LLT representation target.
///
/// This realizes the Schur expansion of [`unicellular_llt_frobenius_target`]
/// as a direct sum of irreducible `S_n` characters. It does not construct the
/// geometric twin-manifold action matrices; rather, it gives the character
/// table that such matrices must have.
pub fn unicellular_llt_character_values_by_degree(
    area: &[u8],
) -> Option<BTreeMap<u32, BTreeMap<Partition, i64>>> {
    let frobenius = unicellular_llt_frobenius_target(area)?;
    let n = area.len() as u32;
    let cycle_types = Partition::all_of_size(n);
    let mut result = BTreeMap::new();

    for (degree, schur_function) in frobenius {
        let mut character_values = BTreeMap::new();
        for cycle_type in &cycle_types {
            let value = schur_function
                .terms()
                .iter()
                .map(|(lambda, &multiplicity)| multiplicity * sn_character(lambda, cycle_type))
                .sum();
            if value != 0 {
                character_values.insert(cycle_type.clone(), value);
            }
        }
        result.insert(degree, character_values);
    }

    Some(result)
}

/// Unit-interval graph edges from an area sequence.
///
/// Vertex `j` is adjacent to `j - gap` for `1 <= gap <= area[j]`.
pub fn unit_interval_edges(area: &[u8]) -> Vec<(usize, usize)> {
    let mut edges = Vec::new();
    for j in 0..area.len() {
        let a = area[j] as usize;
        for gap in 1..=a {
            if gap <= j {
                edges.push((j - gap, j));
            }
        }
    }
    edges.sort_unstable();
    edges
}

fn is_area_sequence(area: &[u8]) -> bool {
    area.iter().enumerate().all(|(i, &v)| v as usize <= i)
        && area
            .windows(2)
            .all(|w| usize::from(w[1]) <= usize::from(w[0]) + 1)
}

fn count_llt_colorings_of_type(
    n: usize,
    lambda: &Partition,
    q_edges: &[(usize, usize)],
    strict_edges: &[(usize, usize)],
    weak_edges: &[(usize, usize)],
) -> UnivariatePolynomial<i64> {
    let parts = lambda.parts();
    let max_ascents = q_edges.len();

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
        if satisfies_strict_weak_constraints(&perm, strict_edges, weak_edges) {
            let ascents = q_edges.iter().filter(|&&(u, v)| perm[u] < perm[v]).count();
            counts[ascents] += 1;
        }
        if !next_multiset_perm(&mut perm) {
            break;
        }
    }

    UnivariatePolynomial::new(counts)
}

fn normalize_edges(edges: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut normalized: Vec<_> = edges
        .iter()
        .map(|&(u, v)| if u < v { (u, v) } else { (v, u) })
        .collect();
    normalized.sort_unstable();
    normalized.dedup();
    normalized
}

fn satisfies_strict_weak_constraints(
    coloring: &[usize],
    strict_edges: &[(usize, usize)],
    weak_edges: &[(usize, usize)],
) -> bool {
    strict_edges.iter().all(|&(u, v)| coloring[u] < coloring[v])
        && weak_edges.iter().all(|&(u, v)| coloring[u] >= coloring[v])
}

fn next_multiset_perm(perm: &mut [usize]) -> bool {
    if perm.len() < 2 {
        return false;
    }

    let mut i = perm.len() - 2;
    while perm[i] >= perm[i + 1] {
        if i == 0 {
            return false;
        }
        i -= 1;
    }

    let mut j = perm.len() - 1;
    while perm[j] <= perm[i] {
        j -= 1;
    }
    perm.swap(i, j);
    perm[i + 1..].reverse();
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frobenius::graded_frobenius_from_character_values;
    use num_rational::Ratio;

    fn p(parts: &[u32]) -> Partition {
        Partition::new(parts.to_vec())
    }

    fn q(n: i64) -> Ratio<i64> {
        Ratio::from_integer(n)
    }

    #[test]
    fn test_unit_interval_edges() {
        assert_eq!(
            unit_interval_edges(&[0, 1, 2]),
            vec![(0, 1), (0, 2), (1, 2)]
        );
    }

    #[test]
    fn test_graph_llt_single_edge() {
        let f = graph_llt_symmetric(2, &[(0, 1)], &[], &[]);
        assert_eq!(
            f.coefficient(&Partition::new(vec![2])),
            UnivariatePolynomial::new(vec![1])
        );
        assert_eq!(
            f.coefficient(&Partition::new(vec![1, 1])),
            UnivariatePolynomial::new(vec![1, 1])
        );
    }

    #[test]
    fn test_graph_llt_strict_edge() {
        let f = graph_llt_symmetric(2, &[(0, 1)], &[(0, 1)], &[]);
        assert_eq!(
            f.coefficient(&Partition::new(vec![2])),
            UnivariatePolynomial::zero()
        );
        assert_eq!(
            f.coefficient(&Partition::new(vec![1, 1])),
            UnivariatePolynomial::new(vec![1])
        );
    }

    #[test]
    fn test_unicellular_llt_matches_complete_graph_k3() {
        let f = unicellular_llt(&[0, 1, 2]);
        assert_eq!(
            f.coefficient(&Partition::new(vec![3])),
            UnivariatePolynomial::new(vec![1])
        );
        assert_eq!(
            f.coefficient(&Partition::new(vec![2, 1])),
            UnivariatePolynomial::new(vec![1, 1, 1])
        );
        assert_eq!(
            f.coefficient(&Partition::new(vec![1, 1, 1])),
            UnivariatePolynomial::new(vec![1, 2, 2, 1])
        );
    }

    #[test]
    fn test_unicellular_llt_frobenius_target_edgeless_s3() {
        let target = unicellular_llt_frobenius_target(&[0, 0, 0]).unwrap();

        assert_eq!(target.keys().copied().collect::<Vec<_>>(), vec![0]);
        assert_eq!(target[&0].coefficient(&p(&[3])), 1);
        assert_eq!(target[&0].coefficient(&p(&[2, 1])), 2);
        assert_eq!(target[&0].coefficient(&p(&[1, 1, 1])), 1);
        assert_eq!(target[&0].terms().len(), 3);
    }

    #[test]
    fn test_unicellular_llt_frobenius_target_rejects_invalid_area() {
        assert!(unicellular_llt_frobenius_target(&[0, 2]).is_none());
        assert!(unicellular_llt_character_values_by_degree(&[0, 2]).is_none());
    }

    #[test]
    fn test_unicellular_llt_characters_reconstruct_frobenius() {
        let area = [0, 1, 1];
        let target = unicellular_llt_frobenius_target(&area).unwrap();
        let character_values = unicellular_llt_character_values_by_degree(&area).unwrap();
        let rational_character_values: BTreeMap<_, BTreeMap<_, _>> = character_values
            .into_iter()
            .map(|(degree, values)| {
                (
                    degree,
                    values
                        .into_iter()
                        .map(|(cycle_type, value)| (cycle_type, q(value)))
                        .collect(),
                )
            })
            .collect();
        let reconstructed = graded_frobenius_from_character_values(&rational_character_values);

        for (&degree, target_degree) in &target {
            let reconstructed_schur = reconstructed[&degree].to_schur_basis();
            let target_schur = target_degree.to_schur_basis();
            for partition in Partition::all_of_size(area.len() as u32) {
                assert_eq!(
                    reconstructed_schur.coefficient(&partition),
                    q(target_schur.coefficient(&partition))
                );
            }
        }
    }

    #[test]
    fn test_unicellular_llt_complete_graph_s3_character_values() {
        let values = unicellular_llt_character_values_by_degree(&[0, 1, 2]).unwrap();

        assert_eq!(values[&0][&p(&[1, 1, 1])], 1);
        assert_eq!(values[&0][&p(&[2, 1])], 1);
        assert_eq!(values[&0][&p(&[3])], 1);
        assert_eq!(values[&3][&p(&[1, 1, 1])], 1);
        assert_eq!(values[&3][&p(&[2, 1])], -1);
        assert_eq!(values[&3][&p(&[3])], 1);
    }
}
