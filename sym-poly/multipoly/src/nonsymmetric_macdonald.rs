//! Nonsymmetric Macdonald polynomials and specializations.
//!
//! This module exposes the full `q,t` filling formula for permuted-basement
//! nonsymmetric Macdonald polynomials, with coefficients in rational functions
//! over `C[q,t]`.
//!
//! The operator-side recursion for the full family also needs an affine
//! `q`-shift operator. That layer is not implemented yet, so the operator
//! interface currently covers the important `q = 0` specialization:
//!
//! - `E_alpha(x; 0, t)`, the nonsymmetric Hall-Littlewood polynomial
//!
//! The tests compare this operator specialization against the filling formula
//! at `q = 0`.

use std::collections::BTreeMap;

use sym_poly_core::{
    qt_constant, qt_unit_monomial, QtPolynomial, QtRationalFunction, RationalFunction, Ring, Ssaf,
};

use crate::key_polynomial::{dominant_rearrangement, sorting_reduced_word};
use crate::multipoly::MultiPoly;
use crate::operators::ttheta_word;

/// Compute the nonsymmetric Hall-Littlewood polynomial `E_alpha(x; 0, t)`.
pub fn nonsymmetric_hall_littlewood<C: Ring>(alpha: &[u32], t: &C) -> MultiPoly<C> {
    let n = alpha.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    let lambda = dominant_rearrangement(alpha);
    let word = sorting_reduced_word(alpha);
    let init = MultiPoly::x_power(n, lambda);
    ttheta_word(&init, &word.iter().rev().copied().collect::<Vec<_>>(), t)
}

/// Compute `E_alpha(x; 0, t)` using Macdonald-family naming.
pub fn nonsymmetric_macdonald_q0<C: Ring>(alpha: &[u32], t: &C) -> MultiPoly<C> {
    nonsymmetric_hall_littlewood(alpha, t)
}

/// Compute the permuted-basement nonsymmetric Macdonald polynomial by the
/// non-attacking filling formula.
///
/// The output is a polynomial in the variables `x_1, ..., x_n` with
/// coefficients in `C(q,t)`.  For a filling `F`, the contribution is
///
/// `q^maj(F) t^coinv(F) x^wt(F) Π (1 - t)/(1 - q^(leg+1)t^(arm+1))`,
///
/// where the product is over the cells returned by
/// [`Ssaf::macdonald_factor_cells`].
pub fn permuted_basement_macdonald_filling_formula<C: Ring>(
    shape: &[u32],
    basement: &[u32],
) -> MultiPoly<QtRationalFunction<C>> {
    assert_eq!(
        shape.len(),
        basement.len(),
        "shape and basement must have the same length"
    );
    let num_vars = shape.len();
    let mut terms: BTreeMap<Vec<u32>, QtRationalFunction<C>> = BTreeMap::new();

    for filling in Ssaf::non_attacking_fillings(shape, basement) {
        let weight = filling.weight_vector();
        let coefficient = macdonald_filling_coefficient(&filling);
        let entry = terms
            .entry(weight)
            .or_insert_with(<QtRationalFunction<C> as Ring>::zero);
        *entry = entry.clone() + coefficient;
    }

    MultiPoly::from_terms(num_vars, terms)
}

/// Identity-basement nonsymmetric Macdonald polynomial by the filling formula.
pub fn nonsymmetric_macdonald_filling_formula<C: Ring>(
    alpha: &[u32],
) -> MultiPoly<QtRationalFunction<C>> {
    let n = alpha.len();
    let basement: Vec<u32> = (1..=n as u32).collect();
    permuted_basement_macdonald_filling_formula(alpha, &basement)
}

fn macdonald_filling_coefficient<C: Ring>(filling: &Ssaf) -> QtRationalFunction<C> {
    let mut coefficient = QtRationalFunction::from_polynomial(qt_unit_monomial(
        filling.major_index() as usize,
        filling.coinversions() as usize,
    ));

    for cell in filling.macdonald_factor_cells() {
        let numerator = one_minus_t();
        let denominator = qt_constant::<C>(1)
            - qt_unit_monomial((cell.leg + 1) as usize, (cell.arm + 1) as usize);
        coefficient = coefficient * RationalFunction::new(numerator, denominator);
    }

    coefficient
}

fn one_minus_t<C: Ring>() -> QtPolynomial<C> {
    qt_constant::<C>(1) - qt_unit_monomial(0, 1)
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;

    use crate::atom_polynomial::{atom_polynomial, t_atom_polynomial};
    use num_rational::Ratio;
    use sym_poly_core::qt_coefficient;

    type Q = Ratio<i64>;
    type Qt = QtPolynomial<Q>;
    type K = QtRationalFunction<Q>;

    fn pow_i64(base: i64, exp: u32) -> i64 {
        (0..exp).fold(1, |acc, _| acc * base)
    }

    fn q0_filling_formula(alpha: &[u32], t: i64) -> MultiPoly<i64> {
        let n = alpha.len();
        let basement: Vec<u32> = (1..=n as u32).collect();
        let mut terms = BTreeMap::new();

        for filling in Ssaf::non_attacking_fillings(alpha, &basement) {
            if filling.major_index() != 0 {
                continue;
            }

            let coeff =
                pow_i64(t, filling.coinversions()) * pow_i64(1 - t, filling.horizontal_descents());
            *terms.entry(filling.weight_vector()).or_insert(0) += coeff;
        }

        MultiPoly::from_terms(n, terms)
    }

    #[test]
    fn test_q0_macdonald_matches_t_atom() {
        let q0: MultiPoly<i64> = nonsymmetric_macdonald_q0(&[1, 0, 2], &2);
        let atom: MultiPoly<i64> = t_atom_polynomial(&[1, 0, 2], &2);
        assert_eq!(q0, atom);
    }

    #[test]
    fn test_q0_macdonald_specializes_to_atom_at_t_zero() {
        let q0: MultiPoly<i64> = nonsymmetric_macdonald_q0(&[0, 2], &0);
        let atom: MultiPoly<i64> = atom_polynomial(&[0, 2]);
        assert_eq!(q0, atom);
    }

    #[test]
    fn test_q0_macdonald_matches_filling_formula() {
        let test_cases: Vec<Vec<u32>> = vec![vec![0, 2], vec![2, 1], vec![1, 0, 2], vec![0, 1, 2]];
        let t_values = [0, 2, -1];

        for alpha in &test_cases {
            for &t in &t_values {
                let operator_side: MultiPoly<i64> = nonsymmetric_macdonald_q0(alpha, &t);
                let filling_side = q0_filling_formula(alpha, t);
                assert_eq!(
                    operator_side, filling_side,
                    "operator/filling mismatch for alpha={alpha:?}, t={t}"
                );
            }
        }
    }

    #[test]
    fn test_q0_macdonald_dominant_is_monomial() {
        let q0: MultiPoly<i64> = nonsymmetric_macdonald_q0(&[3, 2, 1], &5);
        assert_eq!(q0.coefficient(&[3, 2, 1]), 1);
        assert_eq!(q0.terms().len(), 1);
    }

    fn q(value: i64) -> Q {
        Q::from_integer(value)
    }

    fn eval_qt_polynomial(poly: &Qt, q_value: Q, t_value: Q) -> Q {
        let mut result = Q::zero();
        for t_degree in (0..poly.coeffs().len()).rev() {
            let q_poly = poly.coeff(t_degree);
            let mut q_result = Q::zero();
            for q_degree in (0..q_poly.coeffs().len()).rev() {
                q_result = q_result * q_value.clone() + q_poly.coeff(q_degree);
            }
            result = result * t_value.clone() + q_result;
        }
        result
    }

    fn eval_qt_rational_function(f: &K, q_value: Q, t_value: Q) -> Q {
        eval_qt_polynomial(f.numerator(), q_value.clone(), t_value.clone())
            / eval_qt_polynomial(f.denominator(), q_value, t_value)
    }

    #[test]
    fn test_permuted_basement_filling_formula_moura_mandelshtam_example() {
        let poly = permuted_basement_macdonald_filling_formula::<Q>(&[1, 1, 0, 1], &[2, 4, 1, 3]);

        assert_eq!(poly.terms().len(), 3);
        assert!(!poly.coefficient(&[1, 0, 1, 1]).is_zero());
        assert!(!poly.coefficient(&[1, 1, 1, 0]).is_zero());
        assert!(!poly.coefficient(&[0, 1, 1, 1]).is_zero());

        let at_q0_t0 = eval_qt_rational_function(&poly.coefficient(&[0, 1, 1, 1]), q(0), q(0));
        assert_eq!(at_q0_t0, q(1));

        let at_q0_t0 = eval_qt_rational_function(&poly.coefficient(&[1, 1, 1, 0]), q(0), q(0));
        assert_eq!(at_q0_t0, q(0));
    }

    #[test]
    fn test_permuted_basement_filling_formula_q0_matches_operator_side() {
        for alpha in &[vec![0, 2], vec![1, 0, 2], vec![0, 1, 2]] {
            let filling_formula = nonsymmetric_macdonald_filling_formula::<Q>(alpha);
            let operator_side: MultiPoly<i64> = nonsymmetric_macdonald_q0(alpha, &2);

            for (weight, coeff) in filling_formula.terms() {
                let evaluated = eval_qt_rational_function(coeff, q(0), q(2));
                assert_eq!(
                    evaluated,
                    q(operator_side.coefficient(weight)),
                    "q=0 filling/operator mismatch for alpha={alpha:?}, weight={weight:?}"
                );
            }
            for weight in operator_side.terms().keys() {
                assert!(
                    filling_formula.terms().contains_key(weight),
                    "operator side has extra weight {weight:?} for alpha={alpha:?}"
                );
            }
        }
    }

    #[test]
    fn test_permuted_basement_filling_formula_dominant_identity_basement() {
        let poly = nonsymmetric_macdonald_filling_formula::<Q>(&[2, 1]);
        let operator_side: MultiPoly<i64> = nonsymmetric_macdonald_q0(&[2, 1], &2);

        assert_eq!(operator_side.terms().len(), 1);
        assert_eq!(operator_side.coefficient(&[2, 1]), 1);
        assert_eq!(
            eval_qt_rational_function(&poly.coefficient(&[2, 1]), q(0), q(2)),
            q(1)
        );
        for (weight, coeff) in poly.terms() {
            let evaluated = eval_qt_rational_function(coeff, q(0), q(2));
            assert_eq!(
                evaluated,
                q(operator_side.coefficient(weight)),
                "dominant q=0 specialization mismatch at weight {weight:?}"
            );
        }
        assert_eq!(
            qt_coefficient(poly.coefficient(&[2, 1]).numerator(), 0, 0),
            q(1)
        );
    }
}
