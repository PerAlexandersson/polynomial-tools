//! Multivariate polynomial division and normal forms.

use std::collections::{BTreeMap, BTreeSet, HashSet};

use sym_poly_core::{Field, Ring};

use crate::monomial_order::{leading_term, monomial_quotient, LeadingTerm, MonomialOrder};
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
                .map(|(&a, &b)| {
                    a.checked_add(b)
                        .expect("monomial exponent overflow during multiplication")
                })
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

pub fn normal_form_with_leading_terms<C: Field>(
    polynomial: &MultiPoly<C>,
    divisors: &[MultiPoly<C>],
    divisor_leading_terms: &[LeadingTerm<C>],
    order: MonomialOrder,
) -> MultiPoly<C> {
    assert_eq!(
        divisors.len(),
        divisor_leading_terms.len(),
        "divisor and leading-term lists have different lengths"
    );
    assert!(
        divisors
            .iter()
            .all(|divisor| divisor.num_vars() == polynomial.num_vars()),
        "all divisors must have the same number of variables as the polynomial"
    );

    let num_vars = polynomial.num_vars();
    let mut remainder = MultiPoly::zero(num_vars);
    let mut current = polynomial.clone();

    while let Some(current_lt) = leading_term(&current, order) {
        let mut reduced = false;

        for (divisor, divisor_lt) in divisors.iter().zip(divisor_leading_terms.iter()) {
            let Some(exp_quotient) =
                monomial_quotient(&current_lt.exponents, &divisor_lt.exponents)
            else {
                continue;
            };
            let coeff_quotient = current_lt.coefficient.clone() / divisor_lt.coefficient.clone();
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

    remainder
}

/// Compute several normal forms using one F4-style symbolic preprocessing and
/// matrix reduction step.
///
/// This constructs the reducer multiples forced by the current batch of
/// polynomials, stores them as rows with fixed symbolic pivots, and reduces
/// every input row through the same pivot list. For Buchberger computations
/// this gives a shared preprocessing step for a batch of S-polynomials, while
/// scalar reduction remains available for comparisons through
/// [`normal_form_with_leading_terms`].
pub fn matrix_normal_forms_with_leading_terms<C: Field>(
    polynomials: &[MultiPoly<C>],
    divisors: &[MultiPoly<C>],
    divisor_leading_terms: &[LeadingTerm<C>],
    order: MonomialOrder,
) -> Vec<MultiPoly<C>> {
    if polynomials.is_empty() {
        return Vec::new();
    }
    assert_eq!(
        divisors.len(),
        divisor_leading_terms.len(),
        "divisor and leading-term lists have different lengths"
    );
    let num_vars = polynomials[0].num_vars();
    assert!(
        polynomials
            .iter()
            .all(|polynomial| polynomial.num_vars() == num_vars),
        "all polynomials must have the same number of variables"
    );
    assert!(
        divisors
            .iter()
            .all(|divisor| divisor.num_vars() == num_vars),
        "all divisors must have the same number of variables as the polynomials"
    );
    if divisors.is_empty() {
        return polynomials.to_vec();
    }

    let (monomials, reducer_rows) =
        symbolic_reduction_rows(polynomials, divisors, divisor_leading_terms, order);
    if monomials.is_empty() || reducer_rows.is_empty() {
        return polynomials.to_vec();
    }

    let column_index: BTreeMap<_, _> = monomials
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, monomial)| (monomial, i))
        .collect();
    let reducer_matrix: Vec<_> = reducer_rows
        .iter()
        .map(|(_, row)| polynomial_to_row(row, &column_index))
        .collect();
    let reducer_pivots = reducer_rows
        .iter()
        .map(|(pivot, _)| {
            *column_index
                .get(pivot)
                .expect("reducer pivot is outside the matrix columns")
        })
        .collect::<Vec<_>>();

    polynomials
        .iter()
        .map(|polynomial| {
            let mut row = polynomial_to_row(polynomial, &column_index);
            for (reducer_row, &pivot_col) in reducer_matrix.iter().zip(reducer_pivots.iter()) {
                let pivot_value = row[pivot_col].clone();
                if pivot_value.is_zero() {
                    continue;
                }
                for col in pivot_col..row.len() {
                    row[col] = row[col].clone() - pivot_value.clone() * reducer_row[col].clone();
                }
            }
            row_to_polynomial(num_vars, &monomials, &row)
        })
        .collect()
}

fn symbolic_reduction_rows<C: Field>(
    polynomials: &[MultiPoly<C>],
    divisors: &[MultiPoly<C>],
    divisor_leading_terms: &[LeadingTerm<C>],
    order: MonomialOrder,
) -> (Vec<Vec<u32>>, Vec<(Vec<u32>, MultiPoly<C>)>) {
    let mut monomial_set = BTreeSet::new();
    let mut pending = Vec::new();
    for polynomial in polynomials {
        for monomial in polynomial.terms().keys() {
            if monomial_set.insert(monomial.clone()) {
                pending.push(monomial.clone());
            }
        }
    }

    let mut reducer_pivots = HashSet::new();
    let mut reducer_rows = Vec::new();
    while let Some(monomial) = pending.pop() {
        let Some((divisor_index, exponent_quotient)) =
            first_dividing_leading_term(&monomial, divisor_leading_terms)
        else {
            continue;
        };
        if !reducer_pivots.insert(monomial.clone()) {
            continue;
        }
        let divisor_lt = &divisor_leading_terms[divisor_index];
        let row = multiply_by_monomial(
            &divisors[divisor_index],
            &exponent_quotient,
            C::one() / divisor_lt.coefficient.clone(),
        );
        for new_monomial in row.terms().keys() {
            if monomial_set.insert(new_monomial.clone()) {
                pending.push(new_monomial.clone());
            }
        }
        reducer_rows.push((monomial, row));
    }

    let mut monomials: Vec<_> = monomial_set.into_iter().collect();
    monomials.sort_by(|a, b| order.compare(b, a));
    reducer_rows.sort_by(|(a, _), (b, _)| order.compare(b, a));
    (monomials, reducer_rows)
}

fn first_dividing_leading_term<C>(
    monomial: &[u32],
    divisor_leading_terms: &[LeadingTerm<C>],
) -> Option<(usize, Vec<u32>)> {
    divisor_leading_terms
        .iter()
        .enumerate()
        .find_map(|(i, divisor_lt)| {
            monomial_quotient(monomial, &divisor_lt.exponents).map(|quotient| (i, quotient))
        })
}

fn polynomial_to_row<C: Ring>(
    polynomial: &MultiPoly<C>,
    column_index: &BTreeMap<Vec<u32>, usize>,
) -> Vec<C> {
    let mut row = vec![C::zero(); column_index.len()];
    for (monomial, coefficient) in polynomial.terms() {
        let &column = column_index
            .get(monomial)
            .expect("polynomial contains a monomial outside the matrix columns");
        row[column] = coefficient.clone();
    }
    row
}

fn row_to_polynomial<C: Ring>(num_vars: usize, monomials: &[Vec<u32>], row: &[C]) -> MultiPoly<C> {
    let terms = monomials
        .iter()
        .cloned()
        .zip(row.iter().cloned())
        .filter_map(|(monomial, coefficient)| {
            (!coefficient.is_zero()).then_some((monomial, coefficient))
        })
        .collect();
    MultiPoly::from_terms(num_vars, terms)
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
    #[should_panic(expected = "monomial exponent overflow during multiplication")]
    fn test_multiply_by_monomial_rejects_exponent_overflow() {
        let f = MultiPoly::monomial(1, vec![u32::MAX], q(1));
        let _ = multiply_by_monomial(&f, &[1], q(1));
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

    #[test]
    fn test_normal_form_with_cached_leading_terms() {
        let f = mono(&[2, 0], 1) + mono(&[0, 2], 1);
        let divisors = vec![
            mono(&[1, 0], 1) - mono(&[0, 1], 1),
            mono(&[0, 1], 1) - constant(1),
        ];
        let leading_terms = divisors
            .iter()
            .map(|divisor| leading_term(divisor, MonomialOrder::Lex).unwrap())
            .collect::<Vec<_>>();

        assert_eq!(
            normal_form_with_leading_terms(&f, &divisors, &leading_terms, MonomialOrder::Lex),
            normal_form(&f, &divisors, MonomialOrder::Lex)
        );
    }

    #[test]
    fn test_matrix_normal_forms_match_scalar_normal_forms() {
        let polynomials = vec![
            mono(&[2, 0], 1) + mono(&[0, 2], 1),
            mono(&[2, 1], 1) - mono(&[0, 3], 1) + constant(2),
        ];
        let divisors = vec![
            mono(&[1, 0], 1) - mono(&[0, 1], 1),
            mono(&[0, 2], 1) - constant(1),
        ];
        let leading_terms = divisors
            .iter()
            .map(|divisor| leading_term(divisor, MonomialOrder::Lex).unwrap())
            .collect::<Vec<_>>();

        let matrix_forms = matrix_normal_forms_with_leading_terms(
            &polynomials,
            &divisors,
            &leading_terms,
            MonomialOrder::Lex,
        );
        let scalar_forms = polynomials
            .iter()
            .map(|polynomial| {
                normal_form_with_leading_terms(
                    polynomial,
                    &divisors,
                    &leading_terms,
                    MonomialOrder::Lex,
                )
            })
            .collect::<Vec<_>>();

        assert_eq!(matrix_forms, scalar_forms);
    }
}
