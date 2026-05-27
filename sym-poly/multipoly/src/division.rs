//! Multivariate polynomial division and normal forms.

use sym_poly_core::{Field, Ring};

use crate::monomial_order::{leading_term, monomial_quotient, MonomialOrder};
use crate::MultiPoly;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DivisionResult<C: Ring> {
    pub quotients: Vec<MultiPoly<C>>,
    pub remainder: MultiPoly<C>,
}

pub fn multiply_by_monomial<C: Ring>(
    polynomial: &MultiPoly<C>,
    exponents: &[u32],
    coefficient: C,
) -> MultiPoly<C> {
    assert_eq!(
        polynomial.num_vars(),
        exponents.len(),
        "monomial has wrong number of variables"
    );
    if polynomial.is_zero() || coefficient.is_zero() {
        return MultiPoly::zero(polynomial.num_vars());
    }

    let terms = polynomial
        .terms()
        .iter()
        .map(|(term_exp, term_coeff)| {
            let new_exp = term_exp
                .iter()
                .zip(exponents.iter())
                .map(|(&a, &b)| a + b)
                .collect();
            (new_exp, term_coeff.clone() * coefficient.clone())
        })
        .collect();
    MultiPoly::from_terms(polynomial.num_vars(), terms)
}

pub fn divide_by_polynomials<C: Field>(
    dividend: &MultiPoly<C>,
    divisors: &[MultiPoly<C>],
    order: MonomialOrder,
) -> DivisionResult<C> {
    assert!(
        divisors
            .iter()
            .all(|divisor| divisor.num_vars() == dividend.num_vars()),
        "all divisors must have the same number of variables as the dividend"
    );
    assert!(
        divisors.iter().all(|divisor| !divisor.is_zero()),
        "division by a zero polynomial"
    );

    let num_vars = dividend.num_vars();
    let mut quotients = vec![MultiPoly::zero(num_vars); divisors.len()];
    let mut remainder = MultiPoly::zero(num_vars);
    let mut current = dividend.clone();

    while let Some(current_lt) = leading_term(&current, order) {
        let mut reduced = false;

        for (i, divisor) in divisors.iter().enumerate() {
            let divisor_lt = leading_term(divisor, order).expect("nonzero divisor");
            let Some(exp_quotient) =
                monomial_quotient(&current_lt.exponents, &divisor_lt.exponents)
            else {
                continue;
            };
            let coeff_quotient = current_lt.coefficient.clone() / divisor_lt.coefficient;
            let quotient_term =
                MultiPoly::monomial(num_vars, exp_quotient.clone(), coeff_quotient.clone());

            quotients[i] = quotients[i].clone() + quotient_term.clone();
            current = current - multiply_by_monomial(divisor, &exp_quotient, coeff_quotient);
            reduced = true;
            break;
        }

        if !reduced {
            let remainder_term =
                MultiPoly::monomial(num_vars, current_lt.exponents, current_lt.coefficient);
            remainder = remainder + remainder_term.clone();
            current = current - remainder_term;
        }
    }

    DivisionResult {
        quotients,
        remainder,
    }
}

pub fn normal_form<C: Field>(
    polynomial: &MultiPoly<C>,
    divisors: &[MultiPoly<C>],
    order: MonomialOrder,
) -> MultiPoly<C> {
    divide_by_polynomials(polynomial, divisors, order).remainder
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
    fn test_multiply_by_monomial() {
        let f = mono(&[1, 0], 2) + mono(&[0, 1], 3);
        let shifted = multiply_by_monomial(&f, &[2, 1], q(-1));

        assert_eq!(shifted, mono(&[3, 1], -2) + mono(&[2, 2], -3));
    }

    #[test]
    fn test_division_exact_single_divisor() {
        let f = mono(&[2, 0], 1) + mono(&[1, 1], 1);
        let g = mono(&[1, 0], 1) + mono(&[0, 1], 1);
        let result = divide_by_polynomials(&f, &[g.clone()], MonomialOrder::Lex);

        assert_eq!(result.quotients, vec![mono(&[1, 0], 1)]);
        assert!(result.remainder.is_zero());
        assert_eq!(result.quotients[0].clone() * g + result.remainder, f);
    }

    #[test]
    fn test_division_with_remainder() {
        let f = mono(&[2, 0], 1) + mono(&[0, 1], 1);
        let g = mono(&[1, 0], 1) + constant(1);
        let result = divide_by_polynomials(&f, &[g.clone()], MonomialOrder::Lex);

        assert_eq!(result.quotients, vec![mono(&[1, 0], 1) - constant(1)]);
        assert_eq!(result.remainder, mono(&[0, 1], 1) + constant(1));
        assert_eq!(result.quotients[0].clone() * g + result.remainder, f);
    }

    #[test]
    fn test_division_by_ordered_list() {
        let f = mono(&[2, 0], 1) + mono(&[0, 2], 1);
        let g1 = mono(&[1, 0], 1) - mono(&[0, 1], 1);
        let g2 = mono(&[0, 1], 1) - constant(1);
        let result = divide_by_polynomials(&f, &[g1.clone(), g2.clone()], MonomialOrder::Lex);

        let reconstruction = result.quotients[0].clone() * g1
            + result.quotients[1].clone() * g2
            + result.remainder.clone();
        assert_eq!(reconstruction, f);
        assert_eq!(result.remainder, constant(2));
    }

    #[test]
    fn test_normal_form() {
        let f = mono(&[2, 0], 1) + mono(&[0, 1], 1);
        let g = mono(&[1, 0], 1) + constant(1);

        assert_eq!(
            normal_form(&f, &[g], MonomialOrder::Lex),
            mono(&[0, 1], 1) + constant(1)
        );
    }
}
