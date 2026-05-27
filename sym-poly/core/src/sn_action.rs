//! Basic `S_n` actions for exact representation computations.
//!
//! Permutations are zero-indexed one-line notation: `perm[i]` is the image of
//! `i`. Matrices use the column convention: column `j` is the coordinate vector
//! of the image of basis vector `j`.

use std::collections::BTreeMap;

use crate::linear_algebra::{zero_matrix, Matrix};
use crate::{Partition, Ring};

/// Check whether `perm` is a permutation of `{0, ..., n-1}`.
pub fn is_permutation(perm: &[usize]) -> bool {
    let n = perm.len();
    let mut seen = vec![false; n];
    for &value in perm {
        if value >= n || seen[value] {
            return false;
        }
        seen[value] = true;
    }
    true
}

/// Panic unless `perm` is a valid zero-indexed permutation.
pub fn assert_permutation(perm: &[usize]) {
    assert!(
        is_permutation(perm),
        "expected a zero-indexed permutation of 0..{}",
        perm.len()
    );
}

pub fn identity_permutation(n: usize) -> Vec<usize> {
    (0..n).collect()
}

/// The simple transposition `s_i = (i, i+1)`.
pub fn simple_transposition(n: usize, i: usize) -> Vec<usize> {
    assert!(i + 1 < n, "simple transposition index out of range");
    let mut perm = identity_permutation(n);
    perm.swap(i, i + 1);
    perm
}

/// Compose zero-indexed permutations as `a ∘ b`.
pub fn compose_permutations(a: &[usize], b: &[usize]) -> Vec<usize> {
    assert_permutation(a);
    assert_permutation(b);
    assert_eq!(a.len(), b.len(), "permutations have different sizes");
    b.iter().map(|&value| a[value]).collect()
}

pub fn inverse_permutation(perm: &[usize]) -> Vec<usize> {
    assert_permutation(perm);
    let mut inverse = vec![0usize; perm.len()];
    for (i, &value) in perm.iter().enumerate() {
        inverse[value] = i;
    }
    inverse
}

/// Cycle type of a zero-indexed permutation.
pub fn cycle_type(perm: &[usize]) -> Partition {
    assert_permutation(perm);
    let n = perm.len();
    let mut seen = vec![false; n];
    let mut parts = Vec::new();

    for start in 0..n {
        if seen[start] {
            continue;
        }
        let mut length = 0u32;
        let mut current = start;
        while !seen[current] {
            seen[current] = true;
            length += 1;
            current = perm[current];
        }
        parts.push(length);
    }

    Partition::new(parts)
}

/// A canonical zero-indexed representative with the requested cycle type.
pub fn conjugacy_class_representative(cycle_type: &Partition) -> Vec<usize> {
    let n = cycle_type.size() as usize;
    let mut perm = identity_permutation(n);
    let mut offset = 0usize;

    for &part in cycle_type.parts() {
        let length = part as usize;
        if length == 0 {
            continue;
        }
        for i in 0..length {
            perm[offset + i] = offset + ((i + 1) % length);
        }
        offset += length;
    }

    perm
}

/// Canonical conjugacy-class representatives for `S_n`, keyed by cycle type.
pub fn conjugacy_class_representatives(n: usize) -> Vec<(Partition, Vec<usize>)> {
    Partition::all_of_size(n as u32)
        .into_iter()
        .map(|partition| {
            let representative = conjugacy_class_representative(&partition);
            (partition, representative)
        })
        .collect()
}

/// Permutation matrix for the natural permutation representation.
///
/// The action convention is `sigma . e_i = e_{sigma(i)}`.
pub fn permutation_matrix<C: Ring>(perm: &[usize]) -> Matrix<C> {
    assert_permutation(perm);
    let n = perm.len();
    let mut matrix = zero_matrix::<C>(n, n);
    for (source, &target) in perm.iter().enumerate() {
        matrix[target][source] = C::one();
    }
    matrix
}

/// Build an action matrix from a named basis and an action on basis elements.
///
/// The action closure returns a linear combination of basis elements. The
/// output matrix uses the column convention.
pub fn action_matrix_from_basis<C, B, F>(basis: &[B], mut action: F) -> Matrix<C>
where
    C: Ring,
    B: Clone + Ord,
    F: FnMut(&B) -> Vec<(B, C)>,
{
    let mut index = BTreeMap::new();
    for (i, basis_element) in basis.iter().enumerate() {
        let old = index.insert(basis_element.clone(), i);
        assert!(old.is_none(), "basis contains duplicate elements");
    }

    let mut matrix = zero_matrix::<C>(basis.len(), basis.len());
    for (col, basis_element) in basis.iter().enumerate() {
        for (image_basis, coeff) in action(basis_element) {
            if coeff.is_zero() {
                continue;
            }
            let Some(&row) = index.get(&image_basis) else {
                panic!("action image contains an element outside the basis");
            };
            matrix[row][col] = matrix[row][col].clone() + coeff;
        }
    }
    matrix
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linear_algebra::{matrix_trace, quotient_action_matrix, QuotientSpace};
    use num_rational::Ratio;

    type Q = Ratio<i64>;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    #[test]
    fn test_basic_permutation_operations() {
        let s0 = simple_transposition(3, 0);
        let s1 = simple_transposition(3, 1);

        assert_eq!(s0, vec![1, 0, 2]);
        assert_eq!(s1, vec![0, 2, 1]);
        assert_eq!(compose_permutations(&s0, &s1), vec![1, 2, 0]);
        assert_eq!(inverse_permutation(&vec![1, 2, 0]), vec![2, 0, 1]);
    }

    #[test]
    fn test_cycle_type() {
        assert_eq!(cycle_type(&[1, 2, 0, 4, 3]), Partition::new(vec![3, 2]));
        assert_eq!(
            cycle_type(&identity_permutation(4)),
            Partition::new(vec![1, 1, 1, 1])
        );
    }

    #[test]
    fn test_conjugacy_class_representative() {
        let partition = Partition::new(vec![3, 2]);
        let representative = conjugacy_class_representative(&partition);

        assert_eq!(cycle_type(&representative), partition);
        assert_eq!(representative, vec![1, 2, 0, 4, 3]);
    }

    #[test]
    fn test_conjugacy_class_representatives_s3() {
        let reps = conjugacy_class_representatives(3);
        let cycle_types: Vec<_> = reps
            .iter()
            .map(|(partition, _)| partition.clone())
            .collect();

        assert_eq!(
            cycle_types,
            vec![
                Partition::new(vec![3]),
                Partition::new(vec![2, 1]),
                Partition::new(vec![1, 1, 1]),
            ]
        );
        for (partition, representative) in reps {
            assert_eq!(cycle_type(&representative), partition);
        }
    }

    #[test]
    fn test_permutation_matrix_and_basis_action_conventions() {
        let swap = simple_transposition(2, 0);
        let matrix = permutation_matrix::<Q>(&swap);
        let basis = vec![0usize, 1usize];
        let from_basis = action_matrix_from_basis::<Q, _, _>(&basis, |&i| vec![(swap[i], q(1))]);

        assert_eq!(matrix, vec![vec![q(0), q(1)], vec![q(1), q(0)]]);
        assert_eq!(from_basis, matrix);
        assert_eq!(matrix_trace(&matrix), q(0));
    }

    #[test]
    fn test_quotient_action_for_sign_representation_of_s2() {
        let quotient = QuotientSpace::from_relations(2, &[vec![q(1), q(1)]]);
        let swap_matrix = permutation_matrix::<Q>(&simple_transposition(2, 0));
        let quotient_matrix = quotient_action_matrix(&quotient, &swap_matrix);

        assert_eq!(quotient_matrix, vec![vec![q(-1)]]);
        assert_eq!(matrix_trace(&quotient_matrix), q(-1));
    }
}
