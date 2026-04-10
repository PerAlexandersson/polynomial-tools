//! Fundamental slide polynomials.

use sym_poly_core::{Composition, Ring};

use crate::multipoly::MultiPoly;

/// Compute the monomial slide polynomial M_α for a weak composition α.
///
/// M_α = Σ x^b where b dominance-above α and flat(b) = flat(α).
/// See Assaf-Searles (2016), Definition 3.5.
pub fn monomial_slide_polynomial<C: Ring>(alpha: &[u32]) -> MultiPoly<C> {
    let n = alpha.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    let total: u32 = alpha.iter().sum();
    let flat_alpha: Vec<u32> = alpha.iter().copied().filter(|&x| x > 0).collect();
    let mut result = MultiPoly::zero(n);

    for beta in Composition::weak_integer_compositions(total, n) {
        let beta_parts = beta.parts();
        if !prefix_dominates(beta_parts, alpha) {
            continue;
        }

        let flat_beta: Vec<u32> = beta_parts.iter().copied().filter(|&x| x > 0).collect();
        if flat_beta != flat_alpha {
            continue;
        }

        let monomial = MultiPoly::monomial(n, beta_parts.to_vec(), C::one());
        result = result + monomial;
    }

    result
}

/// Compute the fundamental slide polynomial F_α for a weak composition α.
pub fn fundamental_slide_polynomial<C: Ring>(alpha: &[u32]) -> MultiPoly<C> {
    let n = alpha.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    let total: u32 = alpha.iter().sum();
    let flat_alpha: Vec<u32> = alpha.iter().copied().filter(|&x| x > 0).collect();
    let mut result = MultiPoly::zero(n);

    for beta in Composition::weak_integer_compositions(total, n) {
        let beta_parts = beta.parts();
        if !prefix_dominates(beta_parts, alpha) {
            continue;
        }

        let flat_beta: Vec<u32> = beta_parts.iter().copied().filter(|&x| x > 0).collect();
        if !is_refinement(&flat_beta, &flat_alpha) {
            continue;
        }

        let monomial = MultiPoly::monomial(n, beta_parts.to_vec(), C::one());
        result = result + monomial;
    }

    result
}

fn prefix_dominates(beta: &[u32], alpha: &[u32]) -> bool {
    let mut a_sum = 0u32;
    let mut b_sum = 0u32;
    for i in 0..alpha.len() {
        a_sum += alpha[i];
        b_sum += beta[i];
        if b_sum < a_sum {
            return false;
        }
    }
    true
}

fn is_refinement(flat_beta: &[u32], flat_alpha: &[u32]) -> bool {
    if flat_alpha.is_empty() {
        return flat_beta.is_empty();
    }

    let mut i = 0usize;
    let mut j = 0usize;
    let mut running = 0u32;

    while i < flat_beta.len() && j < flat_alpha.len() {
        running += flat_beta[i];
        if running == flat_alpha[j] {
            j += 1;
            running = 0;
        } else if running > flat_alpha[j] {
            return false;
        }
        i += 1;
    }

    i == flat_beta.len() && j == flat_alpha.len() && running == 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_monomial_slide_zero() {
        let m: MultiPoly<i64> = monomial_slide_polynomial(&[0, 0, 0]);
        assert_eq!(m.coefficient(&[0, 0, 0]), 1);
        assert_eq!(m.terms().len(), 1);
    }

    #[test]
    fn test_monomial_slide_assaf_example() {
        // From Assaf-Searles 2016, Eq. 3.7:
        // M_{(0,2,0,3)} = x1^2*x2^3 + x1^2*x3^3 + x1^2*x4^3 + x2^2*x3^3 + x2^2*x4^3
        let m: MultiPoly<i64> = monomial_slide_polynomial(&[0, 2, 0, 3]);
        assert_eq!(m.terms().len(), 5);
        assert_eq!(m.coefficient(&[2, 3, 0, 0]), 1);
        assert_eq!(m.coefficient(&[2, 0, 3, 0]), 1);
        assert_eq!(m.coefficient(&[2, 0, 0, 3]), 1);
        assert_eq!(m.coefficient(&[0, 2, 3, 0]), 1);
        assert_eq!(m.coefficient(&[0, 2, 0, 3]), 1);
    }

    #[test]
    fn test_monomial_slide_dominant_is_monomial() {
        // For dominant α (weakly decreasing), M_α = x^α
        let m: MultiPoly<i64> = monomial_slide_polynomial(&[3, 2, 1]);
        assert_eq!(m.terms().len(), 1);
        assert_eq!(m.coefficient(&[3, 2, 1]), 1);
    }

    #[test]
    fn test_fundamental_slide_zero() {
        let f: MultiPoly<i64> = fundamental_slide_polynomial(&[0, 0, 0]);
        assert_eq!(f.coefficient(&[0, 0, 0]), 1);
        assert_eq!(f.terms().len(), 1);
    }

    #[test]
    fn test_fundamental_slide_102() {
        // F_(1,0,2) = x1*x3^2 + x1*x2*x3 + x1*x2^2
        let f: MultiPoly<i64> = fundamental_slide_polynomial(&[1, 0, 2]);
        assert_eq!(f.coefficient(&[1, 0, 2]), 1);
        assert_eq!(f.coefficient(&[1, 1, 1]), 1);
        assert_eq!(f.coefficient(&[1, 2, 0]), 1);
        assert_eq!(f.terms().len(), 3);
    }

    #[test]
    fn test_fundamental_slide_11() {
        // F_(1,1) = x1*x2
        let f: MultiPoly<i64> = fundamental_slide_polynomial(&[1, 1]);
        assert_eq!(f.coefficient(&[1, 1]), 1);
        assert_eq!(f.terms().len(), 1);
    }
}
