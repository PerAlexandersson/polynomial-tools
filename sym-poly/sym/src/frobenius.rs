//! Frobenius characteristics from symmetric-group character data.
//!
//! The Frobenius characteristic of an `S_n`-module with character `chi` is
//!
//! ```text
//! ch(V) = sum_{mu |- n} chi(mu) / z_mu p_mu.
//! ```
//!
//! This module only constructs the symmetric function. Basis conversion is
//! handled by the existing [`SymmetricFunction`] methods.

use std::collections::BTreeMap;

use sym_poly_core::linear_algebra::{matrix_trace, Matrix};
use sym_poly_core::{Partition, Ring};

use crate::{z_coefficient_i64, Basis, SymmetricFunction};

/// Build the Frobenius characteristic from character values by cycle type.
///
/// Missing conjugacy classes are treated as zero. All supplied cycle types must
/// have the same size.
pub fn frobenius_from_character_values<C: Ring>(
    character_values: &BTreeMap<Partition, C>,
) -> SymmetricFunction<C> {
    let Some(n) = character_values.keys().next().map(Partition::size) else {
        return SymmetricFunction::zero(Basis::PowerSum);
    };
    assert!(
        character_values
            .keys()
            .all(|partition| partition.size() == n),
        "all cycle types must have the same size"
    );

    let mut terms = BTreeMap::new();
    for cycle_type in Partition::all_of_size(n) {
        let value = character_values
            .get(&cycle_type)
            .cloned()
            .unwrap_or_else(C::zero);
        if value.is_zero() {
            continue;
        }
        let coeff = value.exact_div_i64(z_coefficient_i64(&cycle_type));
        if !coeff.is_zero() {
            terms.insert(cycle_type, coeff);
        }
    }

    SymmetricFunction::from_terms(Basis::PowerSum, terms)
}

/// Build the Frobenius characteristic from action matrices keyed by cycle type.
pub fn frobenius_from_trace_matrices<C: Ring>(
    action_matrices: &BTreeMap<Partition, Matrix<C>>,
) -> SymmetricFunction<C> {
    let character_values = action_matrices
        .iter()
        .map(|(cycle_type, matrix)| (cycle_type.clone(), matrix_trace(matrix)))
        .collect();
    frobenius_from_character_values(&character_values)
}

/// Build graded Frobenius characteristics from character values by degree.
pub fn graded_frobenius_from_character_values<C: Ring>(
    graded_character_values: &BTreeMap<u32, BTreeMap<Partition, C>>,
) -> BTreeMap<u32, SymmetricFunction<C>> {
    graded_character_values
        .iter()
        .map(|(&degree, character_values)| {
            (degree, frobenius_from_character_values(character_values))
        })
        .collect()
}

/// Build graded Frobenius characteristics from action matrices by degree.
pub fn graded_frobenius_from_trace_matrices<C: Ring>(
    graded_action_matrices: &BTreeMap<u32, BTreeMap<Partition, Matrix<C>>>,
) -> BTreeMap<u32, SymmetricFunction<C>> {
    graded_action_matrices
        .iter()
        .map(|(&degree, action_matrices)| (degree, frobenius_from_trace_matrices(action_matrices)))
        .collect()
}

/// Build multigraded Frobenius characteristics from character values.
pub fn multigraded_frobenius_from_character_values<C: Ring>(
    multigraded_character_values: &BTreeMap<Vec<u32>, BTreeMap<Partition, C>>,
) -> BTreeMap<Vec<u32>, SymmetricFunction<C>> {
    multigraded_character_values
        .iter()
        .map(|(degree, character_values)| {
            (
                degree.clone(),
                frobenius_from_character_values(character_values),
            )
        })
        .collect()
}

/// Build multigraded Frobenius characteristics from action matrices.
pub fn multigraded_frobenius_from_trace_matrices<C: Ring>(
    multigraded_action_matrices: &BTreeMap<Vec<u32>, BTreeMap<Partition, Matrix<C>>>,
) -> BTreeMap<Vec<u32>, SymmetricFunction<C>> {
    multigraded_action_matrices
        .iter()
        .map(|(degree, action_matrices)| {
            (
                degree.clone(),
                frobenius_from_trace_matrices(action_matrices),
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Ratio;
    use sym_poly_core::{
        left_permutation_basis_action, ungraded_permutation_basis_module, FiniteSnModule,
    };

    type Q = Ratio<i64>;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn p(parts: &[u32]) -> Partition {
        Partition::new(parts.to_vec())
    }

    fn character_table(values: &[(&[u32], i64)]) -> BTreeMap<Partition, Q> {
        values
            .iter()
            .map(|(parts, value)| (p(parts), q(*value)))
            .collect()
    }

    fn permutation_sign(permutation: &[usize]) -> i64 {
        let mut inversions = 0usize;
        for i in 0..permutation.len() {
            for j in i + 1..permutation.len() {
                if permutation[i] > permutation[j] {
                    inversions += 1;
                }
            }
        }
        if inversions % 2 == 0 {
            1
        } else {
            -1
        }
    }

    #[test]
    fn test_trivial_representation_s3() {
        let values = character_table(&[(&[3], 1), (&[2, 1], 1), (&[1, 1, 1], 1)]);

        let frob = frobenius_from_character_values(&values);
        assert_eq!(frob.basis(), Basis::PowerSum);
        assert_eq!(frob.coefficient(&p(&[3])), Ratio::new(1, 3));
        assert_eq!(frob.coefficient(&p(&[2, 1])), Ratio::new(1, 2));
        assert_eq!(frob.coefficient(&p(&[1, 1, 1])), Ratio::new(1, 6));

        let schur = frob.to_schur_basis();
        assert_eq!(schur.coefficient(&p(&[3])), q(1));
        assert_eq!(schur.terms().len(), 1);
    }

    #[test]
    fn test_sign_representation_s3() {
        let values = character_table(&[(&[3], 1), (&[2, 1], -1), (&[1, 1, 1], 1)]);

        let schur = frobenius_from_character_values(&values).to_schur_basis();
        assert_eq!(schur.coefficient(&p(&[1, 1, 1])), q(1));
        assert_eq!(schur.terms().len(), 1);
    }

    #[test]
    fn test_standard_representation_s3() {
        let values = character_table(&[(&[3], -1), (&[2, 1], 0), (&[1, 1, 1], 2)]);

        let schur = frobenius_from_character_values(&values).to_schur_basis();
        assert_eq!(schur.coefficient(&p(&[2, 1])), q(1));
        assert_eq!(schur.terms().len(), 1);
    }

    #[test]
    fn test_frobenius_from_trace_matrices_s2_sign() {
        let matrices =
            BTreeMap::from([(p(&[2]), vec![vec![q(-1)]]), (p(&[1, 1]), vec![vec![q(1)]])]);

        let schur = frobenius_from_trace_matrices(&matrices).to_schur_basis();
        assert_eq!(schur.coefficient(&p(&[1, 1])), q(1));
        assert_eq!(schur.terms().len(), 1);
    }

    #[test]
    fn test_graded_frobenius_from_character_values() {
        let graded_values = BTreeMap::from([
            (0, character_table(&[(&[2], 1), (&[1, 1], 1)])),
            (1, character_table(&[(&[2], -1), (&[1, 1], 1)])),
        ]);

        let graded = graded_frobenius_from_character_values(&graded_values);
        let degree_zero = graded[&0].to_schur_basis();
        let degree_one = graded[&1].to_schur_basis();

        assert_eq!(degree_zero.coefficient(&p(&[2])), q(1));
        assert_eq!(degree_zero.terms().len(), 1);
        assert_eq!(degree_one.coefficient(&p(&[1, 1])), q(1));
        assert_eq!(degree_one.terms().len(), 1);
    }

    #[test]
    fn test_multigraded_frobenius_from_character_values() {
        let multigraded_values = BTreeMap::from([
            (vec![0, 0], character_table(&[(&[2], 1), (&[1, 1], 1)])),
            (vec![1, 0], character_table(&[(&[2], -1), (&[1, 1], 1)])),
        ]);

        let multigraded = multigraded_frobenius_from_character_values(&multigraded_values);
        let degree_zero = multigraded[&vec![0, 0]].to_schur_basis();
        let degree_one = multigraded[&vec![1, 0]].to_schur_basis();

        assert_eq!(degree_zero.coefficient(&p(&[2])), q(1));
        assert_eq!(degree_zero.terms().len(), 1);
        assert_eq!(degree_one.coefficient(&p(&[1, 1])), q(1));
        assert_eq!(degree_one.terms().len(), 1);
    }

    #[test]
    fn test_finite_sn_module_to_multigraded_frobenius() {
        let module =
            FiniteSnModule::<&str, Q>::new(2, vec!["trivial", "sign"], vec![vec![0], vec![1]]);
        let matrices = module.action_matrices_by_degree_and_cycle_type(|permutation, &basis| {
            let coefficient = match basis {
                "trivial" => q(1),
                "sign" => q(permutation_sign(permutation)),
                _ => unreachable!(),
            };
            vec![(basis, coefficient)]
        });

        let frobenius = multigraded_frobenius_from_trace_matrices(&matrices);
        let degree_zero = frobenius[&vec![0]].to_schur_basis();
        let degree_one = frobenius[&vec![1]].to_schur_basis();

        assert_eq!(degree_zero.coefficient(&p(&[2])), q(1));
        assert_eq!(degree_zero.terms().len(), 1);
        assert_eq!(degree_one.coefficient(&p(&[1, 1])), q(1));
        assert_eq!(degree_one.terms().len(), 1);
    }

    #[test]
    fn test_permutation_basis_module_regular_frobenius_s3() {
        let module = ungraded_permutation_basis_module::<Q>(3);
        let matrices = module.action_matrices_by_cycle_type(|permutation, basis| {
            left_permutation_basis_action::<Q>(permutation, basis)
        });

        let schur = frobenius_from_trace_matrices(&matrices).to_schur_basis();

        assert_eq!(schur.coefficient(&p(&[3])), q(1));
        assert_eq!(schur.coefficient(&p(&[2, 1])), q(2));
        assert_eq!(schur.coefficient(&p(&[1, 1, 1])), q(1));
        assert_eq!(schur.terms().len(), 3);
    }
}
