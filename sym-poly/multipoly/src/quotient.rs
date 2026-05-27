//! Finite quotient-basis helpers for Groebner bases.

use std::collections::BTreeMap;

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
}
