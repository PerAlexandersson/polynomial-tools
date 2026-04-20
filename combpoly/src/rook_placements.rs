//! Rook placements on Ferrers boards.
//!
//! Two generating polynomials for placements of non-attacking rooks
//! (at most one per row and column) on a Ferrers board μ:
//!
//! ## Standard rook polynomial
//!
//! `rook_polynomial(μ)`: counts all non-attacking rook placements by number of rooks.
//! Recursion (processing the last row, width μ_ℓ):
//!
//!   R(μ₁, …, μ_ℓ) = R(μ₁, …, μ_{ℓ-1}) + μ_ℓ · t · R(μ₁−1, …, μ_{ℓ-1}−1)
//!
//! Placing a rook in the last row has μ_ℓ column choices; each remaining row
//! loses one available column, so widths decrease by 1.
//!
//! ## Non-nesting rook polynomial
//!
//! `non_nesting_rook_polynomial(μ)`: counts non-attacking rook placements
//! where no rook is strictly north-west of another (equivalently, reading
//! left to right by column, row indices are strictly decreasing).
//! Recursion (processing the last row):
//!
//!   R(μ₁, …, μ_ℓ) = R(μ₁, …, μ_{ℓ-1})
//!                   + t · Σ_{c=1}^{μ_ℓ} R(trim(max(0, μ₁−c), …, max(0, μ_{ℓ-1}−c)))
//!
//! Placing a rook in column c forces all remaining rooks to use columns > c,
//! so the available width of row i becomes max(0, μ_i − c).

use std::collections::HashMap;

/// Compute the standard rook polynomial for a Ferrers board.
///
/// `mu` is a partition in weakly decreasing order (μ₁ ≥ μ₂ ≥ … ≥ μ_ℓ > 0).
/// Returns coefficient vector where `coeffs[k]` counts non-attacking rook
/// placements with exactly `k` rooks (nesting allowed).
pub fn rook_polynomial(mu: &[usize]) -> Vec<i64> {
    let mut cache: HashMap<Vec<usize>, Vec<i64>> = HashMap::new();
    rook_rec(mu, &mut cache)
}

fn rook_rec(mu: &[usize], cache: &mut HashMap<Vec<usize>, Vec<i64>>) -> Vec<i64> {
    let end = mu.iter().rposition(|&x| x > 0).map_or(0, |i| i + 1);
    let mu = &mu[..end];

    if mu.is_empty() {
        return vec![1];
    }

    if let Some(cached) = cache.get(mu) {
        return cached.clone();
    }

    let ell = mu.len();
    let last = mu[ell - 1];
    let mp = &mu[..ell - 1];

    // Last row empty
    let mut result = rook_rec(mp, cache);

    // Rook in last row: last choices, remaining board is (μ₁-1, …, μ_{ℓ-1}-1)
    let sub: Vec<usize> = mp
        .iter()
        .map(|&x| x.saturating_sub(1))
        .filter(|&x| x > 0)
        .collect();
    let sub_poly = rook_rec(&sub, cache);
    let needed = sub_poly.len() + 1;
    if result.len() < needed {
        result.resize(needed, 0);
    }
    for (i, &coeff) in sub_poly.iter().enumerate() {
        result[i + 1] += last as i64 * coeff;
    }

    cache.insert(mu.to_vec(), result.clone());
    result
}

/// Compute the non-nesting rook placement polynomial for a Ferrers board.
///
/// `mu` is a partition in weakly decreasing order (μ₁ ≥ μ₂ ≥ … ≥ μ_ℓ > 0).
/// Returns coefficient vector where `coeffs[k]` counts non-nesting rook
/// placements with exactly `k` rooks.
pub fn non_nesting_rook_polynomial(mu: &[usize]) -> Vec<i64> {
    let mut cache: HashMap<Vec<usize>, Vec<i64>> = HashMap::new();
    nn_rook_rec(mu, &mut cache)
}

fn nn_rook_rec(mu: &[usize], cache: &mut HashMap<Vec<usize>, Vec<i64>>) -> Vec<i64> {
    // Normalize: remove trailing zeros
    let end = mu.iter().rposition(|&x| x > 0).map_or(0, |i| i + 1);
    let mu = &mu[..end];

    if mu.is_empty() {
        return vec![1];
    }

    if let Some(cached) = cache.get(mu) {
        return cached.clone();
    }

    let ell = mu.len();
    let last = mu[ell - 1];
    let mp = &mu[..ell - 1]; // Most[lam]

    // Last row empty: recurse on all-but-last
    let mut result = nn_rook_rec(mp, cache);

    // Rook in column c of last row (c = 1, ..., last)
    for c in 1..=last {
        let sub: Vec<usize> = mp
            .iter()
            .map(|&x| x.saturating_sub(c))
            .filter(|&x| x > 0)
            .collect();
        let sub_poly = nn_rook_rec(&sub, cache);
        // Add t · sub_poly to result
        let needed = sub_poly.len() + 1;
        if result.len() < needed {
            result.resize(needed, 0);
        }
        for (i, &coeff) in sub_poly.iter().enumerate() {
            result[i + 1] += coeff;
        }
    }

    cache.insert(mu.to_vec(), result.clone());
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- Standard rook polynomial tests ---

    #[test]
    fn test_rook_empty() {
        assert_eq!(rook_polynomial(&[]), vec![1]);
    }

    #[test]
    fn test_rook_single_row() {
        assert_eq!(rook_polynomial(&[3]), vec![1, 3]);
    }

    #[test]
    fn test_rook_22() {
        // 2x2 square: 1 + 4t + 2t²
        assert_eq!(rook_polynomial(&[2, 2]), vec![1, 4, 2]);
    }

    #[test]
    fn test_rook_square_33() {
        // 3x3: r_k = C(3,k)² · k!  → [1, 9, 18, 6]
        assert_eq!(rook_polynomial(&[3, 3, 3]), vec![1, 9, 18, 6]);
    }

    #[test]
    fn test_rook_staircase_21() {
        // (2,1): 1 + 3t + t²
        assert_eq!(rook_polynomial(&[2, 1]), vec![1, 3, 1]);
    }

    #[test]
    fn test_rook_staircase_321() {
        // (3,2,1): 1 + 6t + 7t² + t³
        assert_eq!(rook_polynomial(&[3, 2, 1]), vec![1, 6, 7, 1]);
    }

    #[test]
    fn test_rook_staircase_4321() {
        // (4,3,2,1): coefficients are reversed Stirling 2nd kind row for n=5
        // S(5,k) = [1, 15, 25, 10, 1] → reversed [1, 10, 25, 15, 1]
        assert_eq!(rook_polynomial(&[4, 3, 2, 1]), vec![1, 10, 25, 15, 1]);
    }

    // --- Non-nesting rook polynomial tests ---

    #[test]
    fn test_empty() {
        assert_eq!(non_nesting_rook_polynomial(&[]), vec![1]);
    }

    #[test]
    fn test_single_row() {
        // R(n) = 1 + n·t
        assert_eq!(non_nesting_rook_polynomial(&[1]), vec![1, 1]);
        assert_eq!(non_nesting_rook_polynomial(&[3]), vec![1, 3]);
        assert_eq!(non_nesting_rook_polynomial(&[5]), vec![1, 5]);
    }

    #[test]
    fn test_21() {
        // R(2,1) = 1 + 3t + t²
        assert_eq!(non_nesting_rook_polynomial(&[2, 1]), vec![1, 3, 1]);
    }

    #[test]
    fn test_22() {
        // R(2,2) = 1 + 4t + t²
        assert_eq!(non_nesting_rook_polynomial(&[2, 2]), vec![1, 4, 1]);
    }

    #[test]
    fn test_staircase_321() {
        // R(3,2,1) = 1 + 6t + 6t² + t³  (palindromic!)
        assert_eq!(non_nesting_rook_polynomial(&[3, 2, 1]), vec![1, 6, 6, 1]);
    }

    #[test]
    fn test_staircase_4321() {
        let p = non_nesting_rook_polynomial(&[4, 3, 2, 1]);
        // Should be palindromic for the staircase
        let n = p.len();
        for i in 0..n {
            assert_eq!(p[i], p[n - 1 - i], "not palindromic at position {i}");
        }
    }

    #[test]
    fn test_rectangle_33() {
        // 3x2 rectangle
        let p = non_nesting_rook_polynomial(&[3, 3]);
        // By hand: empty → R(3) = 1+3t; c=1: sub={2}, R(2)=1+2t → t+2t²;
        // c=2: sub={1}, R(1)=1+t → t+t²; c=3: sub={}, R()=1 → t
        // Total: (1+3t) + (t+2t²) + (t+t²) + t = 1 + 6t + 3t²
        assert_eq!(p, vec![1, 6, 3]);
    }
}
