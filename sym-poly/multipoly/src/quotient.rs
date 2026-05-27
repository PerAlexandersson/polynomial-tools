//! Finite quotient-basis helpers for Groebner bases.

use std::collections::BTreeMap;

use sym_poly_core::linear_algebra::{zero_matrix, Matrix};
use sym_poly_core::sn_action::assert_permutation;
use sym_poly_core::Field;

use crate::groebner::GroebnerBasis;
use crate::monomial_order::{leading_term, monomial_divides};
use crate::MultiPoly;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QuotientBasis {
    pub num_vars: usize,
    pub monomials: Vec<Vec<u32>>,
}

impl QuotientBasis {
    pub fn dimension(&self) -> usize {
        self.monomials.len()
    }
}

/// Find pure-power bounds from leading monomials.
///
/// If the leading ideal contains `x_i^b` for every variable, then every
/// standard monomial has exponent `< b` in variable `i`. The returned bounds
/// are minimal among pure powers found in the supplied leading monomials.
pub fn pure_power_bounds(num_vars: usize, leading_monomials: &[Vec<u32>]) -> Option<Vec<u32>> {
    let mut bounds: Vec<Option<u32>> = vec![None; num_vars];

    for monomial in leading_monomials {
        assert_eq!(monomial.len(), num_vars, "monomial has wrong length");
        let nonzero_positions: Vec<_> = monomial
            .iter()
            .enumerate()
            .filter_map(|(i, &exp)| (exp > 0).then_some((i, exp)))
            .collect();
        if let [(i, exp)] = nonzero_positions.as_slice() {
            bounds[*i] = Some(bounds[*i].map_or(*exp, |old| old.min(*exp)));
        }
    }

    bounds.into_iter().collect()
}

pub fn standard_monomials_from_leading_monomials(
    num_vars: usize,
    leading_monomials: &[Vec<u32>],
) -> Option<Vec<Vec<u32>>> {
    let bounds = pure_power_bounds(num_vars, leading_monomials)?;
    let mut monomials = Vec::new();
    let mut current = vec![0u32; num_vars];
    enumerate_box(&bounds, 0, &mut current, &mut monomials);
    monomials.retain(|monomial| {
        !leading_monomials
            .iter()
            .any(|leading| monomial_divides(leading, monomial))
    });
    monomials.sort_by_key(|monomial| (monomial.iter().sum::<u32>(), monomial.clone()));
    Some(monomials)
}

pub fn quotient_basis<C: Field>(gb: &GroebnerBasis<C>) -> Option<QuotientBasis> {
    let leading_monomials: Vec<_> = gb
        .generators
        .iter()
        .filter_map(|polynomial| leading_term(polynomial, gb.order).map(|lt| lt.exponents))
        .collect();
    let monomials = standard_monomials_from_leading_monomials(gb.num_vars, &leading_monomials)?;
    Some(QuotientBasis {
        num_vars: gb.num_vars,
        monomials,
    })
}

pub fn quotient_coordinates<C: Field>(
    polynomial: &MultiPoly<C>,
    gb: &GroebnerBasis<C>,
    basis: &QuotientBasis,
) -> Option<Vec<C>> {
    assert_eq!(
        polynomial.num_vars(),
        basis.num_vars,
        "polynomial has wrong number of variables"
    );
    let normal = gb.normal_form(polynomial);
    let index: BTreeMap<_, _> = basis
        .monomials
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, monomial)| (monomial, i))
        .collect();
    let mut coordinates = vec![C::zero(); basis.dimension()];
    for (monomial, coeff) in normal.terms() {
        let &i = index.get(monomial)?;
        coordinates[i] = coeff.clone();
    }
    Some(coordinates)
}

pub fn normal_form_in_basis<C: Field>(
    polynomial: &MultiPoly<C>,
    gb: &GroebnerBasis<C>,
    basis: &QuotientBasis,
) -> Option<MultiPoly<C>> {
    let coordinates = quotient_coordinates(polynomial, gb, basis)?;
    let terms = basis
        .monomials
        .iter()
        .cloned()
        .zip(coordinates)
        .filter_map(|(monomial, coeff)| (!coeff.is_zero()).then_some((monomial, coeff)))
        .collect();
    Some(MultiPoly::from_terms(basis.num_vars, terms))
}

/// Apply a zero-indexed variable permutation to a polynomial.
///
/// The convention is `sigma . x_i = x_{sigma(i)}`.
pub fn permute_variables<C: Field>(
    polynomial: &MultiPoly<C>,
    permutation: &[usize],
) -> MultiPoly<C> {
    assert_permutation(permutation);
    assert_eq!(
        polynomial.num_vars(),
        permutation.len(),
        "permutation has wrong size for polynomial"
    );

    let terms = polynomial
        .terms()
        .iter()
        .map(|(exponents, coeff)| {
            let mut new_exponents = vec![0u32; permutation.len()];
            for (source, &target) in permutation.iter().enumerate() {
                new_exponents[target] = exponents[source];
            }
            (new_exponents, coeff.clone())
        })
        .collect();
    MultiPoly::from_terms(polynomial.num_vars(), terms)
}

/// Matrix induced by a variable permutation on a finite quotient basis.
///
/// Matrices use the column convention: column `j` is the coordinate vector of
/// the image of the `j`th standard monomial.
pub fn quotient_action_matrix_by_permutation<C: Field>(
    gb: &GroebnerBasis<C>,
    basis: &QuotientBasis,
    permutation: &[usize],
) -> Option<Matrix<C>> {
    assert_eq!(
        gb.num_vars, basis.num_vars,
        "Groebner basis and quotient basis have incompatible variable counts"
    );
    assert_permutation(permutation);
    assert_eq!(
        permutation.len(),
        basis.num_vars,
        "permutation has wrong size for quotient basis"
    );

    let dim = basis.dimension();
    let mut matrix = zero_matrix::<C>(dim, dim);
    for (col, monomial) in basis.monomials.iter().enumerate() {
        let basis_element = MultiPoly::x_power(basis.num_vars, monomial.clone());
        let image = permute_variables(&basis_element, permutation);
        let coords = quotient_coordinates(&image, gb, basis)?;
        for row in 0..dim {
            matrix[row][col] = coords[row].clone();
        }
    }
    Some(matrix)
}

pub fn quotient_basis_degrees(basis: &QuotientBasis) -> BTreeMap<u32, Vec<usize>> {
    let mut by_degree = BTreeMap::new();
    for (i, monomial) in basis.monomials.iter().enumerate() {
        by_degree
            .entry(monomial.iter().sum::<u32>())
            .or_insert_with(Vec::new)
            .push(i);
    }
    by_degree
}

pub fn restrict_matrix_to_indices<C: Field>(matrix: &[Vec<C>], indices: &[usize]) -> Matrix<C> {
    let dim = matrix.len();
    assert!(
        matrix.iter().all(|row| row.len() == dim),
        "matrix must be square"
    );
    assert!(
        indices.iter().all(|&i| i < dim),
        "matrix index out of range"
    );

    indices
        .iter()
        .map(|&row| {
            indices
                .iter()
                .map(|&col| matrix[row][col].clone())
                .collect()
        })
        .collect()
}

pub fn quotient_action_matrix_degree_blocks<C: Field>(
    basis: &QuotientBasis,
    action_matrix: &[Vec<C>],
) -> BTreeMap<u32, Matrix<C>> {
    assert_eq!(
        action_matrix.len(),
        basis.dimension(),
        "action matrix has the wrong number of rows"
    );
    assert!(
        action_matrix
            .iter()
            .all(|row| row.len() == basis.dimension()),
        "action matrix has the wrong number of columns"
    );

    quotient_basis_degrees(basis)
        .into_iter()
        .map(|(degree, indices)| (degree, restrict_matrix_to_indices(action_matrix, &indices)))
        .collect()
}

pub fn quotient_action_matrices_by_permutation_and_degree<C: Field>(
    gb: &GroebnerBasis<C>,
    basis: &QuotientBasis,
    permutation: &[usize],
) -> Option<BTreeMap<u32, Matrix<C>>> {
    let action = quotient_action_matrix_by_permutation(gb, basis, permutation)?;
    Some(quotient_action_matrix_degree_blocks(basis, &action))
}

pub fn is_degree_preserving_action_matrix<C: Field>(
    basis: &QuotientBasis,
    action_matrix: &[Vec<C>],
) -> bool {
    if action_matrix.len() != basis.dimension()
        || action_matrix
            .iter()
            .any(|row| row.len() != basis.dimension())
    {
        return false;
    }

    let degrees: Vec<u32> = basis
        .monomials
        .iter()
        .map(|monomial| monomial.iter().sum())
        .collect();
    for row in 0..basis.dimension() {
        for col in 0..basis.dimension() {
            if degrees[row] != degrees[col] && !action_matrix[row][col].is_zero() {
                return false;
            }
        }
    }
    true
}

fn enumerate_box(bounds: &[u32], index: usize, current: &mut [u32], monomials: &mut Vec<Vec<u32>>) {
    if index == bounds.len() {
        monomials.push(current.to_vec());
        return;
    }
    for exp in 0..bounds[index] {
        current[index] = exp;
        enumerate_box(bounds, index + 1, current, monomials);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::MonomialOrder;
    use num_rational::Ratio;
    use sym_poly_core::linear_algebra::matrix_trace;

    type Q = Ratio<i64>;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn mono(exponents: &[u32], coefficient: i64) -> MultiPoly<Q> {
        MultiPoly::monomial(2, exponents.to_vec(), q(coefficient))
    }

    fn constant(value: i64) -> MultiPoly<Q> {
        MultiPoly::constant(2, q(value))
    }

    #[test]
    fn test_pure_power_bounds() {
        let leading = vec![vec![1, 0], vec![0, 2], vec![2, 1]];
        assert_eq!(pure_power_bounds(2, &leading), Some(vec![1, 2]));
        assert_eq!(pure_power_bounds(2, &[vec![1, 0]]), None);
    }

    #[test]
    fn test_standard_monomials_from_leading_monomials() {
        let leading = vec![vec![1, 0], vec![0, 2]];
        assert_eq!(
            standard_monomials_from_leading_monomials(2, &leading),
            Some(vec![vec![0, 0], vec![0, 1]])
        );
    }

    #[test]
    fn test_artin_coinvariant_s2_basis() {
        let e1 = mono(&[1, 0], 1) + mono(&[0, 1], 1);
        let e2 = mono(&[1, 1], 1);
        let gb = GroebnerBasis::new(vec![e1, e2], MonomialOrder::Lex);
        let basis = quotient_basis(&gb).expect("Artin quotient should be finite");

        assert_eq!(basis.monomials, vec![vec![0, 0], vec![0, 1]]);
        assert_eq!(basis.dimension(), 2);
        assert_eq!(gb.normal_form(&mono(&[1, 0], 1)), mono(&[0, 1], -1));
    }

    #[test]
    fn test_quotient_coordinates_artin_s2() {
        let e1 = mono(&[1, 0], 1) + mono(&[0, 1], 1);
        let e2 = mono(&[1, 1], 1);
        let gb = GroebnerBasis::new(vec![e1, e2], MonomialOrder::Lex);
        let basis = quotient_basis(&gb).unwrap();

        assert_eq!(
            quotient_coordinates(&(constant(3) + mono(&[1, 0], 2)), &gb, &basis),
            Some(vec![q(3), q(-2)])
        );
        assert_eq!(
            normal_form_in_basis(&(constant(3) + mono(&[1, 0], 2)), &gb, &basis),
            Some(constant(3) + mono(&[0, 1], -2))
        );
    }

    #[test]
    fn test_permute_variables() {
        let f = mono(&[2, 0], 1) + mono(&[1, 1], 3);
        let swapped = permute_variables(&f, &[1, 0]);

        assert_eq!(swapped, mono(&[0, 2], 1) + mono(&[1, 1], 3));
    }

    #[test]
    fn test_quotient_action_matrix_artin_s2() {
        let e1 = mono(&[1, 0], 1) + mono(&[0, 1], 1);
        let e2 = mono(&[1, 1], 1);
        let gb = GroebnerBasis::new(vec![e1, e2], MonomialOrder::Lex);
        let basis = quotient_basis(&gb).unwrap();
        let action = quotient_action_matrix_by_permutation(&gb, &basis, &[1, 0]).unwrap();

        assert_eq!(basis.monomials, vec![vec![0, 0], vec![0, 1]]);
        assert_eq!(action, vec![vec![q(1), q(0)], vec![q(0), q(-1)]]);
        assert_eq!(matrix_trace(&action), q(0));
        assert_eq!(
            quotient_basis_degrees(&basis),
            BTreeMap::from([(0, vec![0]), (1, vec![1])])
        );
        assert!(is_degree_preserving_action_matrix(&basis, &action));
        assert_eq!(
            quotient_action_matrix_degree_blocks(&basis, &action),
            BTreeMap::from([(0, vec![vec![q(1)]]), (1, vec![vec![q(-1)]])])
        );
        assert_eq!(
            quotient_action_matrices_by_permutation_and_degree(&gb, &basis, &[1, 0]).unwrap(),
            BTreeMap::from([(0, vec![vec![q(1)]]), (1, vec![vec![q(-1)]])])
        );
    }
}
