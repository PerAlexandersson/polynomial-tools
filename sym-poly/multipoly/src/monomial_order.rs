//! Monomial orders for sparse multivariate polynomials.

use std::cmp::Ordering;

use sym_poly_core::Ring;

use crate::MultiPoly;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MonomialOrder {
    /// Lexicographic order with `x_1 > x_2 > ... > x_n`.
    Lex,
    /// Total degree, then lexicographic order.
    GrLex,
    /// Total degree, then graded reverse lexicographic order.
    GrevLex,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeadingTerm<C> {
    pub exponents: Vec<u32>,
    pub coefficient: C,
}

impl MonomialOrder {
    pub fn compare(&self, a: &[u32], b: &[u32]) -> Ordering {
        assert_eq!(a.len(), b.len(), "monomials have different lengths");
        match self {
            MonomialOrder::Lex => compare_lex(a, b),
            MonomialOrder::GrLex => compare_total_degree(a, b).then_with(|| compare_lex(a, b)),
            MonomialOrder::GrevLex => {
                compare_total_degree(a, b).then_with(|| compare_grevlex_same_degree(a, b))
            }
        }
    }
}

pub fn leading_term<C: Ring>(
    polynomial: &MultiPoly<C>,
    order: MonomialOrder,
) -> Option<LeadingTerm<C>> {
    let term = match order {
        MonomialOrder::Lex => polynomial.terms().iter().next_back(),
        MonomialOrder::GrLex | MonomialOrder::GrevLex => polynomial
            .terms()
            .iter()
            .max_by(|(a, _), (b, _)| order.compare(a, b)),
    };
    term.map(|(exponents, coefficient)| LeadingTerm {
        exponents: exponents.clone(),
        coefficient: coefficient.clone(),
    })
}

pub fn monomial_divides(divisor: &[u32], dividend: &[u32]) -> bool {
    assert_eq!(
        divisor.len(),
        dividend.len(),
        "monomials have different lengths"
    );
    divisor.iter().zip(dividend.iter()).all(|(&a, &b)| a <= b)
}

pub fn monomial_quotient(dividend: &[u32], divisor: &[u32]) -> Option<Vec<u32>> {
    if !monomial_divides(divisor, dividend) {
        return None;
    }
    Some(
        dividend
            .iter()
            .zip(divisor.iter())
            .map(|(&a, &b)| a - b)
            .collect(),
    )
}

fn compare_total_degree(a: &[u32], b: &[u32]) -> Ordering {
    a.iter().sum::<u32>().cmp(&b.iter().sum::<u32>())
}

fn compare_lex(a: &[u32], b: &[u32]) -> Ordering {
    for (&x, &y) in a.iter().zip(b.iter()) {
        match x.cmp(&y) {
            Ordering::Equal => continue,
            ordering => return ordering,
        }
    }
    Ordering::Equal
}

fn compare_grevlex_same_degree(a: &[u32], b: &[u32]) -> Ordering {
    for (&x, &y) in a.iter().zip(b.iter()).rev() {
        match x.cmp(&y) {
            Ordering::Equal => continue,
            Ordering::Less => return Ordering::Greater,
            Ordering::Greater => return Ordering::Less,
        }
    }
    Ordering::Equal
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_monomial_order_lex() {
        let order = MonomialOrder::Lex;
        assert_eq!(order.compare(&[2, 0], &[1, 10]), Ordering::Greater);
        assert_eq!(order.compare(&[1, 1], &[1, 2]), Ordering::Less);
    }

    #[test]
    fn test_monomial_order_grlex() {
        let order = MonomialOrder::GrLex;
        assert_eq!(order.compare(&[1, 1], &[2, 0]), Ordering::Less);
        assert_eq!(order.compare(&[1, 2], &[2, 0]), Ordering::Greater);
    }

    #[test]
    fn test_monomial_order_grevlex() {
        let order = MonomialOrder::GrevLex;
        assert_eq!(order.compare(&[2, 0], &[1, 1]), Ordering::Greater);
        assert_eq!(order.compare(&[0, 2], &[1, 1]), Ordering::Less);
        assert_eq!(order.compare(&[1, 2], &[2, 0]), Ordering::Greater);
    }

    #[test]
    fn test_leading_term() {
        let f = MultiPoly::<i64>::x_power(2, vec![1, 4]) + MultiPoly::<i64>::x_power(2, vec![2, 0]);

        assert_eq!(
            leading_term(&f, MonomialOrder::Lex).unwrap().exponents,
            vec![2, 0]
        );
        assert_eq!(
            leading_term(&f, MonomialOrder::GrLex).unwrap().exponents,
            vec![1, 4]
        );
    }

    #[test]
    fn test_monomial_divisibility_and_quotient() {
        assert!(monomial_divides(&[1, 0, 2], &[3, 0, 2]));
        assert!(!monomial_divides(&[1, 1, 2], &[3, 0, 2]));
        assert_eq!(
            monomial_quotient(&[3, 0, 2], &[1, 0, 2]),
            Some(vec![2, 0, 0])
        );
    }
}
