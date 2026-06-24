//! Fundamental slide polynomials.

use std::collections::{BTreeMap, BTreeSet};

use combinatoric_core::WeakComposition;
use sym_poly_core::{Composition, Ring};

use crate::basis::MultiPolyBasis;
use crate::multipoly::{checked_total_degree, MultiPoly};
use crate::multipoly_function::MultiPolyFunction;

/// Compute the monomial slide polynomial M_α for a weak composition α.
///
/// M_α = Σ x^b where b dominance-above α and flat(b) = flat(α).
/// See Assaf-Searles (2016), Definition 3.5.
pub fn monomial_slide_polynomial<C: Ring>(alpha: &[u32]) -> MultiPoly<C> {
    let n = alpha.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    let total = checked_total_degree(alpha);
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

    let total = checked_total_degree(alpha);
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

/// Expand a sparse polynomial in the fundamental slide basis.
///
/// This uses the unitriangularity of fundamental slides: every monomial in
/// `F_alpha` prefix-dominates `alpha`, and `x^alpha` appears with coefficient
/// one.  Processing weak compositions by increasing prefix sums therefore
/// gives a support-driven conversion without building a dense transition
/// matrix.
pub fn fundamental_slide_expansion<C: Ring>(poly: &MultiPoly<C>) -> MultiPolyFunction<C> {
    let num_vars = poly.num_vars();
    if poly.is_zero() {
        return MultiPolyFunction::zero(MultiPolyBasis::FundSlide, num_vars);
    }

    let mut degrees = BTreeSet::new();
    for exp in poly.terms().keys() {
        degrees.insert(checked_total_degree(exp));
    }

    let mut result_terms = BTreeMap::new();
    for degree in degrees {
        let mut current: BTreeMap<Vec<u32>, C> = poly
            .terms()
            .iter()
            .filter(|(exp, _)| checked_total_degree(exp) == degree)
            .map(|(exp, coeff)| (exp.clone(), coeff.clone()))
            .collect();

        let mut weak_compositions = WeakComposition::all_weak_compositions(degree, num_vars);
        weak_compositions.sort_by_key(|alpha| prefix_sums(alpha.parts()));

        for alpha in weak_compositions {
            let coeff = current.get(alpha.parts()).cloned().unwrap_or_else(C::zero);
            if coeff.is_zero() {
                continue;
            }
            result_terms.insert(alpha.clone(), coeff.clone());

            let slide = fundamental_slide_polynomial::<C>(alpha.parts());
            for (beta, beta_coeff) in slide.terms() {
                let new_coeff = current.get(beta).cloned().unwrap_or_else(C::zero)
                    - coeff.clone() * beta_coeff.clone();
                if new_coeff.is_zero() {
                    current.remove(beta);
                } else {
                    current.insert(beta.clone(), new_coeff);
                }
            }
        }

        assert!(
            current.is_empty(),
            "fundamental slide expansion left a nonzero remainder"
        );
    }

    MultiPolyFunction::from_terms(MultiPolyBasis::FundSlide, num_vars, result_terms)
}

fn prefix_dominates(beta: &[u32], alpha: &[u32]) -> bool {
    let mut a_sum = 0u32;
    let mut b_sum = 0u32;
    for i in 0..alpha.len() {
        a_sum = a_sum
            .checked_add(alpha[i])
            .expect("composition weight overflow");
        b_sum = b_sum
            .checked_add(beta[i])
            .expect("composition weight overflow");
        if b_sum < a_sum {
            return false;
        }
    }
    true
}

fn prefix_sums(alpha: &[u32]) -> Vec<u32> {
    let mut running = 0u32;
    alpha
        .iter()
        .map(|&part| {
            running = running
                .checked_add(part)
                .expect("composition weight overflow");
            running
        })
        .collect()
}

fn is_refinement(flat_beta: &[u32], flat_alpha: &[u32]) -> bool {
    if flat_alpha.is_empty() {
        return flat_beta.is_empty();
    }

    let mut i = 0usize;
    let mut j = 0usize;
    let mut running = 0u32;

    while i < flat_beta.len() && j < flat_alpha.len() {
        running = running
            .checked_add(flat_beta[i])
            .expect("composition weight overflow");
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

    #[test]
    fn test_fundamental_slide_expansion_round_trip() {
        let p = fundamental_slide_polynomial::<i64>(&[1, 0, 2])
            + fundamental_slide_polynomial::<i64>(&[0, 3, 0]).scale(&2);
        let expansion = fundamental_slide_expansion(&p);
        assert_eq!(
            expansion.coefficient(&WeakComposition::from_slice(&[1, 0, 2])),
            1
        );
        assert_eq!(
            expansion.coefficient(&WeakComposition::from_slice(&[0, 3, 0])),
            2
        );
        assert_eq!(expansion.terms().len(), 2);
        assert_eq!(expansion.to_multipoly(), p);
    }

    #[test]
    fn test_fundamental_slide_expansion_matches_dense_transition() {
        let p = MultiPoly::monomial(3, vec![1, 0, 2], 1)
            + MultiPoly::monomial(3, vec![0, 3, 0], 2)
            + MultiPoly::monomial(3, vec![2, 1, 0], 3);
        let fast = fundamental_slide_expansion(&p);
        let dense = MultiPolyFunction::from_multipoly(&p).to_fund_slide_basis();
        assert_eq!(fast, dense);
    }

    #[test]
    #[should_panic(expected = "total degree overflow")]
    fn test_slide_rejects_weight_overflow() {
        let _: MultiPoly<i64> = monomial_slide_polynomial(&[u32::MAX, 1]);
    }
}
