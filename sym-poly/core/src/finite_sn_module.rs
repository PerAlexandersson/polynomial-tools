//! Finite graded `S_n` modules.
//!
//! This module provides a small common representation layer for finite
//! `S_n`-modules. It is intentionally basis-oriented: callers provide a finite
//! basis and an action on basis elements, and the module builds matrices by
//! permutation, by conjugacy class, and by multidegree block.

use std::collections::{BTreeMap, BTreeSet};
use std::marker::PhantomData;

use crate::linear_algebra::{matrix_trace, zero_matrix, Matrix};
use crate::sn_action::{
    action_matrix_from_basis, assert_permutation, compose_permutations,
    conjugacy_class_representatives, inverse_permutation,
};
use crate::{Partition, Ring};

pub type PermutationBasis = Vec<usize>;

/// A finite `S_n` module with a chosen homogeneous basis.
///
/// Basis elements are caller-defined objects. Multidegrees are vectors so the
/// same type covers ordinary grading, bigrading, and ungraded modules.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteSnModule<B, C>
where
    B: Clone + Ord,
    C: Ring,
{
    symmetric_group_degree: usize,
    basis: Vec<B>,
    multidegrees: Vec<Vec<u32>>,
    coefficient: PhantomData<C>,
}

impl<B, C> FiniteSnModule<B, C>
where
    B: Clone + Ord,
    C: Ring,
{
    /// Build a finite `S_n` module from a basis and matching multidegrees.
    pub fn new(symmetric_group_degree: usize, basis: Vec<B>, multidegrees: Vec<Vec<u32>>) -> Self {
        assert_eq!(
            basis.len(),
            multidegrees.len(),
            "basis and multidegree lists have different lengths"
        );
        assert_unique_basis(&basis);
        Self {
            symmetric_group_degree,
            basis,
            multidegrees,
            coefficient: PhantomData,
        }
    }

    /// Build an ungraded module. The unique degree is the empty vector.
    pub fn ungraded(symmetric_group_degree: usize, basis: Vec<B>) -> Self {
        let multidegrees = vec![Vec::new(); basis.len()];
        Self::new(symmetric_group_degree, basis, multidegrees)
    }

    pub fn symmetric_group_degree(&self) -> usize {
        self.symmetric_group_degree
    }

    pub fn basis(&self) -> &[B] {
        &self.basis
    }

    pub fn multidegrees(&self) -> &[Vec<u32>] {
        &self.multidegrees
    }

    pub fn dimension(&self) -> usize {
        self.basis.len()
    }

    pub fn degree_blocks(&self) -> BTreeMap<Vec<u32>, Vec<usize>> {
        let mut blocks = BTreeMap::new();
        for (i, degree) in self.multidegrees.iter().enumerate() {
            blocks
                .entry(degree.clone())
                .or_insert_with(Vec::new)
                .push(i);
        }
        blocks
    }

    /// Build the action matrix of one permutation.
    ///
    /// The action closure receives the permutation and a basis element, and
    /// returns a linear combination in the same basis. Matrices use the column
    /// convention: column `j` is the image of basis vector `j`.
    pub fn action_matrix_by_permutation<F>(&self, permutation: &[usize], mut action: F) -> Matrix<C>
    where
        F: FnMut(&[usize], &B) -> Vec<(B, C)>,
    {
        assert_permutation(permutation);
        assert_eq!(
            permutation.len(),
            self.symmetric_group_degree,
            "permutation has the wrong size for this S_n module"
        );
        action_matrix_from_basis::<C, B, _>(&self.basis, |basis_element| {
            action(permutation, basis_element)
        })
    }

    /// Build action matrices for canonical conjugacy-class representatives.
    pub fn action_matrices_by_cycle_type<F>(&self, mut action: F) -> BTreeMap<Partition, Matrix<C>>
    where
        F: FnMut(&[usize], &B) -> Vec<(B, C)>,
    {
        conjugacy_class_representatives(self.symmetric_group_degree)
            .into_iter()
            .map(|(cycle_type, representative)| {
                let matrix = self
                    .action_matrix_by_permutation(&representative, |permutation, basis| {
                        action(permutation, basis)
                    });
                (cycle_type, matrix)
            })
            .collect()
    }

    /// Character values on canonical conjugacy-class representatives.
    pub fn character_values_by_cycle_type<F>(&self, action: F) -> BTreeMap<Partition, C>
    where
        F: FnMut(&[usize], &B) -> Vec<(B, C)>,
    {
        self.action_matrices_by_cycle_type(action)
            .into_iter()
            .map(|(cycle_type, matrix)| (cycle_type, matrix_trace(&matrix)))
            .collect()
    }

    /// Restrict an action matrix to multidegree blocks.
    pub fn action_matrix_blocks_by_degree(
        &self,
        action_matrix: &[Vec<C>],
    ) -> BTreeMap<Vec<u32>, Matrix<C>> {
        assert_square_matrix(action_matrix, self.dimension());
        assert!(
            self.action_matrix_preserves_multidegrees(action_matrix),
            "action matrix is not multidegree-preserving"
        );
        self.degree_blocks()
            .into_iter()
            .map(|(degree, indices)| (degree, restrict_matrix_to_indices(action_matrix, &indices)))
            .collect()
    }

    /// Check whether an action matrix preserves the chosen multigrading.
    pub fn action_matrix_preserves_multidegrees(&self, action_matrix: &[Vec<C>]) -> bool {
        assert_square_matrix(action_matrix, self.dimension());
        for row in 0..self.dimension() {
            for col in 0..self.dimension() {
                if self.multidegrees[row] != self.multidegrees[col]
                    && !action_matrix[row][col].is_zero()
                {
                    return false;
                }
            }
        }
        true
    }

    /// Build conjugacy-class action matrices and then restrict each to
    /// multidegree blocks.
    pub fn action_matrices_by_degree_and_cycle_type<F>(
        &self,
        action: F,
    ) -> BTreeMap<Vec<u32>, BTreeMap<Partition, Matrix<C>>>
    where
        F: FnMut(&[usize], &B) -> Vec<(B, C)>,
    {
        let matrices = self.action_matrices_by_cycle_type(action);
        let mut result: BTreeMap<Vec<u32>, BTreeMap<Partition, Matrix<C>>> = BTreeMap::new();
        for (cycle_type, matrix) in matrices {
            for (degree, block) in self.action_matrix_blocks_by_degree(&matrix) {
                result
                    .entry(degree)
                    .or_default()
                    .insert(cycle_type.clone(), block);
            }
        }
        result
    }

    /// Multigraded character values on canonical conjugacy-class representatives.
    pub fn character_values_by_degree_and_cycle_type<F>(
        &self,
        action: F,
    ) -> BTreeMap<Vec<u32>, BTreeMap<Partition, C>>
    where
        F: FnMut(&[usize], &B) -> Vec<(B, C)>,
    {
        self.action_matrices_by_degree_and_cycle_type(action)
            .into_iter()
            .map(|(degree, matrices)| {
                let values = matrices
                    .into_iter()
                    .map(|(cycle_type, matrix)| (cycle_type, matrix_trace(&matrix)))
                    .collect();
                (degree, values)
            })
            .collect()
    }
}

/// The fixed-point basis indexed by permutations in `S_n`.
pub fn symmetric_group_permutation_basis(n: usize) -> Vec<PermutationBasis> {
    combinatoric_core::all_permutations_zero_indexed(n)
}

/// A finite `S_n` module whose basis is indexed by permutations in `S_n`.
///
/// This is the natural fixed-point label set for Hessenberg/GKM experiments.
/// The caller supplies a multidegree for each permutation label.
pub fn permutation_basis_module<C, F>(
    symmetric_group_degree: usize,
    mut multidegree: F,
) -> FiniteSnModule<PermutationBasis, C>
where
    C: Ring,
    F: FnMut(&[usize]) -> Vec<u32>,
{
    let basis = symmetric_group_permutation_basis(symmetric_group_degree);
    let multidegrees = basis
        .iter()
        .map(|permutation| multidegree(permutation))
        .collect();
    FiniteSnModule::new(symmetric_group_degree, basis, multidegrees)
}

/// Ungraded permutation-label module on the basis `S_n`.
pub fn ungraded_permutation_basis_module<C: Ring>(
    symmetric_group_degree: usize,
) -> FiniteSnModule<PermutationBasis, C> {
    FiniteSnModule::ungraded(
        symmetric_group_degree,
        symmetric_group_permutation_basis(symmetric_group_degree),
    )
}

/// Left action on permutation labels: `sigma . w = sigma ∘ w`.
pub fn left_permutation_basis_action<C: Ring>(
    permutation: &[usize],
    basis_element: &[usize],
) -> Vec<(PermutationBasis, C)> {
    assert_permutation(permutation);
    assert_permutation(basis_element);
    assert_eq!(
        permutation.len(),
        basis_element.len(),
        "permutation and basis label have different sizes"
    );
    vec![(compose_permutations(permutation, basis_element), C::one())]
}

/// Right-regular left action on permutation labels: `sigma . w = w ∘ sigma^{-1}`.
///
/// The inverse makes this a left `S_n` action while retaining the usual right
/// multiplication convention on labels.
pub fn right_permutation_basis_action<C: Ring>(
    permutation: &[usize],
    basis_element: &[usize],
) -> Vec<(PermutationBasis, C)> {
    assert_permutation(permutation);
    assert_permutation(basis_element);
    assert_eq!(
        permutation.len(),
        basis_element.len(),
        "permutation and basis label have different sizes"
    );
    vec![(
        compose_permutations(basis_element, &inverse_permutation(permutation)),
        C::one(),
    )]
}

fn assert_unique_basis<B: Clone + Ord>(basis: &[B]) {
    let mut seen = BTreeSet::new();
    for basis_element in basis {
        assert!(
            seen.insert(basis_element),
            "basis contains duplicate elements"
        );
    }
}

fn assert_square_matrix<C: Ring>(matrix: &[Vec<C>], expected_dimension: usize) {
    assert_eq!(
        matrix.len(),
        expected_dimension,
        "matrix has the wrong number of rows"
    );
    assert!(
        matrix.iter().all(|row| row.len() == expected_dimension),
        "matrix has the wrong number of columns"
    );
}

fn restrict_matrix_to_indices<C: Ring>(matrix: &[Vec<C>], indices: &[usize]) -> Matrix<C> {
    let dimension = matrix.len();
    assert!(
        indices.iter().all(|&index| index < dimension),
        "matrix index out of range"
    );
    let mut block = zero_matrix::<C>(indices.len(), indices.len());
    for (block_row, &row) in indices.iter().enumerate() {
        for (block_col, &col) in indices.iter().enumerate() {
            block[block_row][block_col] = matrix[row][col].clone();
        }
    }
    block
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linear_algebra::matrix_multiply;
    use crate::sn_action::{cycle_type, identity_permutation, simple_transposition};
    use num_rational::Ratio;

    type Q = Ratio<i64>;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn p(parts: &[u32]) -> Partition {
        Partition::new(parts.to_vec())
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
    fn test_ungraded_permutation_representation_s3_characters() {
        let module = FiniteSnModule::<usize, Q>::ungraded(3, vec![0, 1, 2]);
        let values = module
            .character_values_by_cycle_type(|permutation, &basis| vec![(permutation[basis], q(1))]);

        assert_eq!(values[&p(&[3])], q(0));
        assert_eq!(values[&p(&[2, 1])], q(1));
        assert_eq!(values[&p(&[1, 1, 1])], q(3));
    }

    #[test]
    fn test_permutation_basis_left_regular_characters() {
        let module = ungraded_permutation_basis_module::<Q>(3);
        let values = module.character_values_by_cycle_type(|permutation, basis| {
            left_permutation_basis_action::<Q>(permutation, basis)
        });

        assert_eq!(module.dimension(), 6);
        assert_eq!(values[&p(&[3])], q(0));
        assert_eq!(values[&p(&[2, 1])], q(0));
        assert_eq!(values[&p(&[1, 1, 1])], q(6));
    }

    #[test]
    fn test_permutation_basis_right_regular_characters() {
        let module = ungraded_permutation_basis_module::<Q>(3);
        let values = module.character_values_by_cycle_type(|permutation, basis| {
            right_permutation_basis_action::<Q>(permutation, basis)
        });

        assert_eq!(values[&p(&[3])], q(0));
        assert_eq!(values[&p(&[2, 1])], q(0));
        assert_eq!(values[&p(&[1, 1, 1])], q(6));
    }

    #[test]
    fn test_permutation_basis_actions_are_left_representations() {
        let module = ungraded_permutation_basis_module::<Q>(3);
        let s0 = simple_transposition(3, 0);
        let s1 = simple_transposition(3, 1);
        let product = crate::sn_action::compose_permutations(&s0, &s1);

        let left_s0 = module.action_matrix_by_permutation(&s0, |permutation, basis| {
            left_permutation_basis_action::<Q>(permutation, basis)
        });
        let left_s1 = module.action_matrix_by_permutation(&s1, |permutation, basis| {
            left_permutation_basis_action::<Q>(permutation, basis)
        });
        let left_product = module.action_matrix_by_permutation(&product, |permutation, basis| {
            left_permutation_basis_action::<Q>(permutation, basis)
        });
        assert_eq!(matrix_multiply(&left_s0, &left_s1), left_product);

        let right_s0 = module.action_matrix_by_permutation(&s0, |permutation, basis| {
            right_permutation_basis_action::<Q>(permutation, basis)
        });
        let right_s1 = module.action_matrix_by_permutation(&s1, |permutation, basis| {
            right_permutation_basis_action::<Q>(permutation, basis)
        });
        let right_product = module.action_matrix_by_permutation(&product, |permutation, basis| {
            right_permutation_basis_action::<Q>(permutation, basis)
        });
        assert_eq!(matrix_multiply(&right_s0, &right_s1), right_product);
    }

    #[test]
    fn test_permutation_basis_module_accepts_custom_multidegrees() {
        let module = permutation_basis_module::<Q, _>(3, |permutation| {
            vec![permutation.iter().filter(|&&value| value == 0).count() as u32]
        });

        assert_eq!(module.dimension(), 6);
        assert_eq!(
            module.degree_blocks(),
            BTreeMap::from([(vec![1], vec![0, 1, 2, 3, 4, 5])])
        );
    }

    #[test]
    fn test_graded_trivial_plus_sign_module_s2() {
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

        assert_eq!(module.dimension(), 2);
        assert_eq!(
            module.degree_blocks(),
            BTreeMap::from([(vec![0], vec![0]), (vec![1], vec![1])])
        );
        assert_eq!(matrices[&vec![0]][&p(&[2])], vec![vec![q(1)]]);
        assert_eq!(matrices[&vec![0]][&p(&[1, 1])], vec![vec![q(1)]]);
        assert_eq!(matrices[&vec![1]][&p(&[2])], vec![vec![q(-1)]]);
        assert_eq!(matrices[&vec![1]][&p(&[1, 1])], vec![vec![q(1)]]);
    }

    #[test]
    fn test_detects_action_matrix_that_mixes_degrees() {
        let module = FiniteSnModule::<usize, Q>::new(2, vec![0, 1], vec![vec![0], vec![1]]);
        let swap = crate::sn_action::permutation_matrix::<Q>(&simple_transposition(2, 0));

        assert!(!module.action_matrix_preserves_multidegrees(&swap));
    }

    #[test]
    fn test_action_matrix_by_permutation_checks_size() {
        let module = FiniteSnModule::<usize, Q>::ungraded(3, vec![0, 1, 2]);
        let identity = identity_permutation(3);
        let matrix = module.action_matrix_by_permutation(&identity, |permutation, &basis| {
            vec![(permutation[basis], q(1))]
        });

        assert_eq!(matrix, crate::sn_action::permutation_matrix::<Q>(&identity));
        assert_eq!(cycle_type(&simple_transposition(3, 0)), p(&[2, 1]));
    }
}
