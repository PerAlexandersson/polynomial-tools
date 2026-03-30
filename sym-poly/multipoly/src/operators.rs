//! Divided difference operators for multivariate polynomials.
//!
//! - `partial_i(f)` = (f - s_i(f)) / (x_i - x_{i+1})  (simple divided difference)
//! - `pi_i(f)` = partial_i(x_i * f)  (isobaric/Demazure operator)
//!
//! Both always produce polynomials when applied to polynomials.

use sym_poly_core::Ring;
use crate::multipoly::MultiPoly;

/// Simple divided difference ∂_i: (f - s_i(f)) / (x_i - x_{i+1}).
///
/// The i-th simple divided difference (0-indexed). Satisfies:
/// - ∂_i^2 = 0
/// - ∂_i ∂_j = ∂_j ∂_i when |i-j| >= 2
/// - ∂_i ∂_{i+1} ∂_i = ∂_{i+1} ∂_i ∂_{i+1} (braid relation)
pub fn partial_i<C: Ring>(f: &MultiPoly<C>, i: usize) -> MultiPoly<C> {
    // Compute f - s_i(f), then divide by (x_i - x_{i+1})
    let si_f = f.swap_vars(i);
    let numerator = f.clone() - si_f;

    // Division by (x_i - x_{i+1}):
    // For each term c * x^e in numerator, we need to express it as
    // (x_i - x_{i+1}) * something.
    //
    // Key identity: x_i^a * x_{i+1}^b - x_i^b * x_{i+1}^a = (x_i - x_{i+1}) * S
    // where S = Σ_{k=0}^{a-b-1} x_i^{a-1-k} * x_{i+1}^{b+k} (for a > b).
    //
    // Since f - s_i(f) is antisymmetric in x_i, x_{i+1}, it's always divisible.
    //
    // Algorithm: process terms greedily. For each monomial where exp[i] > exp[i+1],
    // contribute (exp[i] - exp[i+1]) terms to the quotient.
    divide_by_xi_minus_xi1(&numerator, i)
}

/// Isobaric (Demazure) divided difference π_i = ∂_i(x_i * f).
///
/// Equivalently: π_i(f) = (x_i * f - x_{i+1} * s_i(f)) / (x_i - x_{i+1}).
///
/// Satisfies:
/// - π_i^2 = π_i (idempotent)
/// - π_i π_{i+1} π_i = π_{i+1} π_i π_{i+1} (braid relation)
pub fn pi_i<C: Ring>(f: &MultiPoly<C>, i: usize) -> MultiPoly<C> {
    let xi_f = f.mul_var(i);
    partial_i(&xi_f, i)
}

/// Divide a polynomial by (x_i - x_{i+1}).
///
/// The polynomial must be antisymmetric in x_i, x_{i+1} (guaranteed when
/// this is called from partial_i on f - s_i(f)).
fn divide_by_xi_minus_xi1<C: Ring>(f: &MultiPoly<C>, i: usize) -> MultiPoly<C> {
    let n = f.num_vars();
    let mut result = MultiPoly::zero(n);

    // For each term with exp[i] > exp[i+1], expand the geometric series
    for (exp, coeff) in f.terms() {
        let a = exp[i];
        let b = exp[i + 1];
        if a <= b {
            // Terms with a <= b should cancel in pairs due to antisymmetry.
            // We only process a > b terms.
            continue;
        }
        // x_i^a * x_{i+1}^b / (x_i - x_{i+1})
        // = x_i^{a-1} * x_{i+1}^b + x_i^{a-2} * x_{i+1}^{b+1} + ... + x_i^b * x_{i+1}^{a-1}
        for k in 0..(a - b) {
            let mut new_exp = exp.clone();
            new_exp[i] = a - 1 - k;
            new_exp[i + 1] = b + k;
            let entry = result.terms.entry(new_exp).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone();
        }
    }

    result.strip_zeros();
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_partial_i_squared_is_zero() {
        // ∂_0^2 = 0. Test on x_1^3 * x_2.
        let f: MultiPoly<i64> = MultiPoly::x_power(3, vec![3, 1, 0]);
        let d1 = partial_i(&f, 0);
        let d2 = partial_i(&d1, 0);
        assert!(d2.is_zero(), "∂_0^2 should be 0, got {:?}", d2);
    }

    #[test]
    fn test_partial_i_on_xi() {
        // ∂_0(x_1) = (x_1 - x_2) / (x_1 - x_2) = 1
        let x1: MultiPoly<i64> = MultiPoly::var(2, 0);
        let result = partial_i(&x1, 0);
        assert_eq!(result.coefficient(&[0, 0]), 1);
        assert_eq!(result.terms().len(), 1);
    }

    #[test]
    fn test_partial_i_on_x1_squared() {
        // ∂_0(x_1^2) = (x_1^2 - x_2^2)/(x_1 - x_2) = x_1 + x_2
        let f: MultiPoly<i64> = MultiPoly::x_power(2, vec![2, 0]);
        let result = partial_i(&f, 0);
        assert_eq!(result.coefficient(&[1, 0]), 1);
        assert_eq!(result.coefficient(&[0, 1]), 1);
        assert_eq!(result.terms().len(), 2);
    }

    #[test]
    fn test_pi_idempotent() {
        // π_0^2 = π_0. Test on x_1^2 * x_2.
        let f: MultiPoly<i64> = MultiPoly::x_power(3, vec![2, 1, 0]);
        let p1 = pi_i(&f, 0);
        let p2 = pi_i(&p1, 0);
        assert_eq!(p1, p2, "π_0 should be idempotent");
    }

    #[test]
    fn test_pi_on_x1() {
        // π_0(x_1) = ∂_0(x_1^2) = x_1 + x_2
        let x1: MultiPoly<i64> = MultiPoly::var(2, 0);
        let result = pi_i(&x1, 0);
        assert_eq!(result.coefficient(&[1, 0]), 1);
        assert_eq!(result.coefficient(&[0, 1]), 1);
    }

    #[test]
    fn test_braid_relation_partial() {
        // ∂_0 ∂_1 ∂_0 = ∂_1 ∂_0 ∂_1 on x_1^3 * x_2^2 * x_3
        let f: MultiPoly<i64> = MultiPoly::x_power(3, vec![3, 2, 1]);
        let lhs = partial_i(&partial_i(&partial_i(&f, 0), 1), 0);
        let rhs = partial_i(&partial_i(&partial_i(&f, 1), 0), 1);
        assert_eq!(lhs, rhs, "braid relation for ∂ failed");
    }

    #[test]
    fn test_braid_relation_pi() {
        // π_0 π_1 π_0 = π_1 π_0 π_1 on x_1^3 * x_2
        let f: MultiPoly<i64> = MultiPoly::x_power(3, vec![3, 1, 0]);
        let lhs = pi_i(&pi_i(&pi_i(&f, 0), 1), 0);
        let rhs = pi_i(&pi_i(&pi_i(&f, 1), 0), 1);
        assert_eq!(lhs, rhs, "braid relation for π failed");
    }

    #[test]
    fn test_commutation_partial() {
        // ∂_0 ∂_2 = ∂_2 ∂_0 (|0-2| >= 2)
        let f: MultiPoly<i64> = MultiPoly::x_power(4, vec![3, 1, 2, 1]);
        let lhs = partial_i(&partial_i(&f, 0), 2);
        let rhs = partial_i(&partial_i(&f, 2), 0);
        assert_eq!(lhs, rhs, "commutation for ∂ failed");
    }
}
