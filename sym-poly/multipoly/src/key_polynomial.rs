//! Key polynomials (Demazure characters) via divided difference operators.
//!
//! The key polynomial κ_α for a weak composition α is defined by:
//!   κ_α = π_{w}(x^{sort(α)})
//! where w is the sorting permutation (the permutation that sorts α into
//! weakly decreasing order) and π_w = π_{i_1} ... π_{i_k} for any
//! reduced word i_1 ... i_k of w.

use crate::multipoly::MultiPoly;
use crate::operators::{pi_i, tpi_word};
use sym_poly_core::Ring;

/// Compute the key polynomial κ_α for a weak composition α = (α_1, ..., α_n).
///
/// The algorithm:
/// 1. Compute the dominant rearrangement λ = sort(α) (weakly decreasing)
/// 2. Find a reduced word for the sorting permutation w (so that w(λ) = α)
/// 3. Apply π_{i_1} ... π_{i_k} to x^λ
pub fn key_polynomial<C: Ring>(alpha: &[u32]) -> MultiPoly<C> {
    let n = alpha.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    // Step 1: Sort α into weakly decreasing order to get λ
    let lambda = dominant_rearrangement(alpha);

    // Step 2: Find the reduced word for the permutation w such that
    // applying w to λ gives α. We use bubble sort to find adjacent transpositions.
    let reduced_word = sorting_reduced_word(alpha);

    // Step 3: Start with x^λ and apply π operators
    let mut result = MultiPoly::x_power(n, lambda);

    // Apply operators right-to-left: π_{i_k} first, then π_{i_{k-1}}, etc.
    // The reduced word gives w = s_{i_1} ... s_{i_k}, and
    // κ_α = π_w(x^λ) = π_{i_1} ... π_{i_k}(x^λ)
    for &i in reduced_word.iter().rev() {
        result = pi_i(&result, i);
    }

    result
}

/// t-deformed operator key polynomial.
pub fn t_key_polynomial<C: Ring>(alpha: &[u32], t: &C) -> MultiPoly<C> {
    let n = alpha.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    let lambda = dominant_rearrangement(alpha);
    let reduced_word = sorting_reduced_word(alpha);
    let init = MultiPoly::x_power(n, lambda);
    tpi_word(
        &init,
        &reduced_word.iter().rev().copied().collect::<Vec<_>>(),
        t,
    )
}

pub(crate) fn dominant_rearrangement(alpha: &[u32]) -> Vec<u32> {
    let mut lambda = alpha.to_vec();
    lambda.sort_unstable_by(|a, b| b.cmp(a));
    lambda
}

/// Find a reduced word for the permutation that sorts α into weakly
/// decreasing order, expressed as a sequence of simple transpositions.
///
/// Uses bubble sort: each swap s_i that moves a larger value left contributes
/// to the reduced word.
pub(crate) fn sorting_reduced_word(alpha: &[u32]) -> Vec<usize> {
    let mut perm: Vec<u32> = alpha.to_vec();
    let mut word = Vec::new();
    let n = perm.len();

    if n < 2 {
        return word;
    }

    // Bubble sort into weakly decreasing order
    loop {
        let mut swapped = false;
        for i in 0..n - 1 {
            if perm[i] < perm[i + 1] {
                perm.swap(i, i + 1);
                word.push(i);
                swapped = true;
            }
        }
        if !swapped {
            break;
        }
    }

    word
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_key_dominant() {
        // κ_{(3,2,1)} = x_1^3 * x_2^2 * x_3 (already dominant, identity permutation)
        let kp: MultiPoly<i64> = key_polynomial(&[3, 2, 1]);
        assert_eq!(kp.coefficient(&[3, 2, 1]), 1);
        assert_eq!(kp.terms().len(), 1);
    }

    #[test]
    fn test_key_single_swap() {
        // κ_{(2,3,1)}: λ = (3,2,1), need s_0 to get (2,3,1)
        // κ_{(2,3,1)} = π_0(x_1^3 x_2^2 x_3) = x_1^3 x_2^2 x_3 ... no wait
        // Let me think: sorted is (3,2,1). To get (2,3,1) from (3,2,1),
        // we swap positions 0,1: s_0(3,2,1) = (2,3,1). So w = s_0.
        // κ = π_0(x^{3,2,1})
        //   = ∂_0(x_1 * x_1^3 * x_2^2 * x_3)
        //   = ∂_0(x_1^4 * x_2^2 * x_3)
        //   = (x_1^4 x_2^2 x_3 - x_1^2 x_2^4 x_3) / (x_1 - x_2)
        //   = x_3 * (x_1^4 x_2^2 - x_1^2 x_2^4) / (x_1 - x_2)
        //   = x_3 * x_1^2 x_2^2 * (x_1^2 - x_2^2) / (x_1 - x_2)
        //   = x_3 * x_1^2 x_2^2 * (x_1 + x_2)
        //   = x_1^3 x_2^2 x_3 + x_1^2 x_2^3 x_3
        let kp: MultiPoly<i64> = key_polynomial(&[2, 3, 1]);
        assert_eq!(kp.coefficient(&[3, 2, 1]), 1);
        assert_eq!(kp.coefficient(&[2, 3, 1]), 1);
        assert_eq!(kp.terms().len(), 2);
    }

    #[test]
    fn test_key_identity_perm() {
        // For any dominant composition (weakly decreasing), κ_α = x^α
        let kp: MultiPoly<i64> = key_polynomial(&[5, 3, 1]);
        assert_eq!(kp.coefficient(&[5, 3, 1]), 1);
        assert_eq!(kp.terms().len(), 1);
    }

    #[test]
    fn test_key_w0_n2() {
        // κ_{(0,2)} for S_2: λ = (2,0), w = s_0 (longest element)
        // κ = π_0(x_1^2) = x_1^2 + x_1*x_2 + x_2^2
        // Wait: π_0(x_1^2) = ∂_0(x_1 * x_1^2) = ∂_0(x_1^3)
        // = (x_1^3 - x_2^3)/(x_1 - x_2) = x_1^2 + x_1 x_2 + x_2^2
        let kp: MultiPoly<i64> = key_polynomial(&[0, 2]);
        assert_eq!(kp.coefficient(&[2, 0]), 1);
        assert_eq!(kp.coefficient(&[1, 1]), 1);
        assert_eq!(kp.coefficient(&[0, 2]), 1);
        assert_eq!(kp.terms().len(), 3);
    }

    #[test]
    fn test_key_positive_coefficients() {
        // Key polynomials always have positive coefficients
        let test_cases: Vec<Vec<u32>> = vec![
            vec![1, 2],
            vec![0, 3],
            vec![1, 0, 2],
            vec![2, 0, 1],
            vec![0, 1, 2],
        ];
        for alpha in &test_cases {
            let kp: MultiPoly<i64> = key_polynomial(alpha);
            for (_, c) in kp.terms() {
                assert!(*c > 0, "negative coefficient in κ_{:?}", alpha);
            }
        }
    }

    #[test]
    fn test_sorting_reduced_word() {
        // (2,3,1) -> (3,2,1): swap at position 0
        assert_eq!(sorting_reduced_word(&[2, 3, 1]), vec![0]);

        // (3,2,1) -> already sorted
        assert_eq!(sorting_reduced_word(&[3, 2, 1]), vec![]);

        // (1,2,3) -> (3,2,1): w_0 = s_0 s_1 s_0
        let word = sorting_reduced_word(&[1, 2, 3]);
        assert_eq!(word.len(), 3); // length of w_0 in S_3
    }

    #[test]
    fn test_t_key_specialization_at_zero() {
        let key: MultiPoly<i64> = key_polynomial(&[2, 3, 1]);
        let t_key: MultiPoly<i64> = t_key_polynomial(&[2, 3, 1], &0);
        assert_eq!(key, t_key);
    }
}
