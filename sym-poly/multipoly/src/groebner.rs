//! Basic Buchberger algorithm for small exact computations.
//!
//! This is intentionally simple. It is meant to support small research
//! quotient-ring examples before adding Buchberger criteria or specialized
//! symmetric-ideal shortcuts.

use sym_poly_core::{Field, Ring};

use crate::division::{multiply_by_monomial, normal_form};
use crate::monomial_order::{leading_term, monomial_divides, MonomialOrder};
use crate::MultiPoly;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroebnerBasis<C: Ring> {
    pub num_vars: usize,
    pub order: MonomialOrder,
    pub generators: Vec<MultiPoly<C>>,
}

impl<C: Field> GroebnerBasis<C> {
    pub fn new(generators: Vec<MultiPoly<C>>, order: MonomialOrder) -> Self {
        let generators = reduced_groebner_basis(&generators, order);
        let num_vars = generators.first().map(MultiPoly::num_vars).unwrap_or(0);
        Self {
            num_vars,
            order,
            generators,
        }
    }

    pub fn normal_form(&self, polynomial: &MultiPoly<C>) -> MultiPoly<C> {
        normal_form(polynomial, &self.generators, self.order)
    }
}

pub fn s_polynomial<C: Field>(
    f: &MultiPoly<C>,
    g: &MultiPoly<C>,
    order: MonomialOrder,
) -> MultiPoly<C> {
    assert_eq!(
        f.num_vars(),
        g.num_vars(),
        "polynomials have different rings"
    );
    let lt_f = leading_term(f, order).expect("first polynomial is zero");
    let lt_g = leading_term(g, order).expect("second polynomial is zero");
    let lcm = monomial_lcm(&lt_f.exponents, &lt_g.exponents);
    let exp_f = monomial_difference(&lcm, &lt_f.exponents);
    let exp_g = monomial_difference(&lcm, &lt_g.exponents);

    multiply_by_monomial(f, &exp_f, C::one() / lt_f.coefficient)
        - multiply_by_monomial(g, &exp_g, C::one() / lt_g.coefficient)
}

pub fn buchberger_basis<C: Field>(
    generators: &[MultiPoly<C>],
    order: MonomialOrder,
) -> Vec<MultiPoly<C>> {
    if generators.is_empty() {
        return Vec::new();
    }
    let num_vars = generators[0].num_vars();
    assert!(
        generators
            .iter()
            .all(|polynomial| polynomial.num_vars() == num_vars),
        "all generators must have the same number of variables"
    );

    let mut basis: Vec<_> = generators
        .iter()
        .filter(|polynomial| !polynomial.is_zero())
        .map(|polynomial| make_monic(polynomial, order))
        .collect();
    let mut pairs = all_pairs(basis.len());

    while let Some((i, j)) = pairs.pop() {
        let s = s_polynomial(&basis[i], &basis[j], order);
        let remainder = normal_form(&s, &basis, order);
        if remainder.is_zero() {
            continue;
        }
        let new_index = basis.len();
        basis.push(make_monic(&remainder, order));
        for i in 0..new_index {
            pairs.push((i, new_index));
        }
    }

    basis
}

pub fn reduced_groebner_basis<C: Field>(
    generators: &[MultiPoly<C>],
    order: MonomialOrder,
) -> Vec<MultiPoly<C>> {
    let basis = buchberger_basis(generators, order);
    let mut reduced = Vec::new();

    for i in 0..basis.len() {
        let others: Vec<_> = basis
            .iter()
            .enumerate()
            .filter_map(|(j, polynomial)| (i != j).then_some(polynomial.clone()))
            .collect();
        let remainder = normal_form(&basis[i], &others, order);
        if !remainder.is_zero() {
            reduced.push(make_monic(&remainder, order));
        }
    }

    remove_leading_monomial_multiples(&mut reduced, order);
    reduced.sort_by(|a, b| {
        let a_lt = leading_term(a, order).expect("nonzero polynomial");
        let b_lt = leading_term(b, order).expect("nonzero polynomial");
        order.compare(&b_lt.exponents, &a_lt.exponents)
    });
    reduced
}

pub fn is_groebner_basis<C: Field>(basis: &[MultiPoly<C>], order: MonomialOrder) -> bool {
    for (i, j) in all_pairs(basis.len()) {
        let s = s_polynomial(&basis[i], &basis[j], order);
        if !normal_form(&s, basis, order).is_zero() {
            return false;
        }
    }
    true
}

pub fn make_monic<C: Field>(polynomial: &MultiPoly<C>, order: MonomialOrder) -> MultiPoly<C> {
    let lt = leading_term(polynomial, order).expect("zero polynomial has no monic form");
    polynomial.scale(&(C::one() / lt.coefficient))
}

fn monomial_lcm(a: &[u32], b: &[u32]) -> Vec<u32> {
    assert_eq!(a.len(), b.len(), "monomials have different lengths");
    a.iter().zip(b.iter()).map(|(&x, &y)| x.max(y)).collect()
}

fn monomial_difference(a: &[u32], b: &[u32]) -> Vec<u32> {
    assert_eq!(a.len(), b.len(), "monomials have different lengths");
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| {
            assert!(x >= y, "monomial difference would have negative exponent");
            x - y
        })
        .collect()
}

fn all_pairs(n: usize) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            pairs.push((i, j));
        }
    }
    pairs
}

fn remove_leading_monomial_multiples<C: Ring>(basis: &mut Vec<MultiPoly<C>>, order: MonomialOrder) {
    let leading_terms: Vec<_> = basis
        .iter()
        .map(|polynomial| leading_term(polynomial, order).expect("nonzero polynomial"))
        .collect();
    let keep: Vec<_> = leading_terms
        .iter()
        .enumerate()
        .map(|(i, lt_i)| {
            !leading_terms
                .iter()
                .enumerate()
                .any(|(j, lt_j)| i != j && monomial_divides(&lt_j.exponents, &lt_i.exponents))
        })
        .collect();

    let old_basis = std::mem::take(basis);
    *basis = old_basis
        .into_iter()
        .enumerate()
        .filter_map(|(i, polynomial)| keep[i].then_some(polynomial))
        .collect();
}

#[cfg(test)]
mod tests {
    use super::*;
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
    fn test_s_polynomial_cancels_leading_terms() {
        let f = mono(&[1, 0], 1) + mono(&[0, 1], 1);
        let g = mono(&[1, 0], 1) + constant(-1);
        let s = s_polynomial(&f, &g, MonomialOrder::Lex);

        assert_eq!(
            leading_term(&s, MonomialOrder::Lex).unwrap().exponents,
            vec![0, 1]
        );
    }

    #[test]
    fn test_buchberger_for_already_groebner_basis() {
        let g1 = mono(&[1, 0], 1) + mono(&[0, 1], 1);
        let g2 = mono(&[0, 2], 1) + constant(-1);
        let basis = buchberger_basis(&[g1, g2], MonomialOrder::Lex);

        assert!(is_groebner_basis(&basis, MonomialOrder::Lex));
    }

    #[test]
    fn test_reduced_groebner_basis_for_two_linear_relations() {
        let g1 = mono(&[1, 0], 1) + mono(&[0, 1], 1);
        let g2 = mono(&[0, 1], 1) + constant(-1);
        let basis = reduced_groebner_basis(&[g1, g2], MonomialOrder::Lex);

        assert_eq!(
            basis,
            vec![
                mono(&[1, 0], 1) + constant(1),
                mono(&[0, 1], 1) + constant(-1),
            ]
        );
        assert!(is_groebner_basis(&basis, MonomialOrder::Lex));
    }

    #[test]
    fn test_groebner_basis_normal_form() {
        let g1 = mono(&[1, 0], 1) + mono(&[0, 1], 1);
        let g2 = mono(&[0, 1], 1) + constant(-1);
        let gb = GroebnerBasis::new(vec![g1, g2], MonomialOrder::Lex);
        let f = mono(&[2, 0], 1) + mono(&[0, 2], 1);

        assert_eq!(gb.normal_form(&f), constant(2));
    }

    #[test]
    fn test_buchberger_adds_missing_s_polynomial() {
        let f = mono(&[2, 0], 1) - mono(&[0, 1], 1);
        let g = mono(&[1, 1], 1) - constant(1);
        let basis = buchberger_basis(&[f, g], MonomialOrder::Lex);

        assert!(basis.len() > 2);
        assert!(is_groebner_basis(&basis, MonomialOrder::Lex));
    }
}
