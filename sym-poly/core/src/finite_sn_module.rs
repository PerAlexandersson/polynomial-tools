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
    action_matrix_from_basis, assert_permutation, conjugacy_class_representatives,
};
use crate::{Partition, Ring};

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
