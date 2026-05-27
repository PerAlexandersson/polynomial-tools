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

use polynomial_tools::{check_weak_interlacing, is_real_rooted};

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

/// Trim trailing zero coefficients, keeping at least one entry.
pub fn trim_coefficient_vector(mut p: Vec<i64>) -> Vec<i64> {
    while p.len() > 1 && p.last().copied() == Some(0) {
        p.pop();
    }
    p
}

/// Enumerate integer partitions of `n` in weakly decreasing order.
pub fn integer_partitions(n: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut buf = Vec::new();
    integer_partitions_rec(n, n, &mut buf, &mut result);
    result
}

fn integer_partitions_rec(
    n: usize,
    max_part: usize,
    buf: &mut Vec<usize>,
    result: &mut Vec<Vec<usize>>,
) {
    if n == 0 {
        result.push(buf.clone());
        return;
    }
    for k in (1..=n.min(max_part)).rev() {
        buf.push(k);
        integer_partitions_rec(n - k, k, buf, result);
        buf.pop();
    }
}

/// Strip the first `c` columns from a Ferrers shape.
pub fn strip_ferrers_columns(mu: &[usize], c: usize) -> Vec<usize> {
    mu.iter()
        .map(|&x| x.saturating_sub(c))
        .filter(|&x| x > 0)
        .collect()
}

fn binom_i128(n: i64, k: i64) -> i128 {
    if k < 0 || n < 0 || k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut num = 1i128;
    let mut den = 1i128;
    for i in 0..k {
        num *= (n - i) as i128;
        den *= (i + 1) as i128;
    }
    num / den
}

fn det_i128(mat: &[Vec<i128>]) -> i128 {
    let n = mat.len();
    if n == 0 {
        return 1;
    }
    if n == 1 {
        return mat[0][0];
    }
    let mut total = 0;
    for j in 0..n {
        if mat[0][j] == 0 {
            continue;
        }
        let mut sub = Vec::with_capacity(n - 1);
        for row in mat.iter().skip(1) {
            let mut sub_row = Vec::with_capacity(n - 1);
            for (col, &entry) in row.iter().enumerate() {
                if col != j {
                    sub_row.push(entry);
                }
            }
            sub.push(sub_row);
        }
        let sign = if j % 2 == 0 { 1 } else { -1 };
        total += sign * mat[0][j] * det_i128(&sub);
    }
    total
}

/// Count non-nesting placements using exactly the rows in `rows_mask`, after
/// stripping the first `strip` columns from a Ferrers shape.
///
/// This is the fixed-row LGV determinant used by the adjacent-packet
/// diagnostics. Rows are indexed from top to bottom, starting at zero.
pub fn fixed_row_non_nesting_count_lgv(eta: &[usize], strip: usize, rows_mask: u64) -> i64 {
    assert!(
        eta.len() <= u64::BITS as usize,
        "row masks support at most 64 rows"
    );
    let rows: Vec<usize> = (0..eta.len())
        .filter(|&r| (rows_mask & (1u64 << r)) != 0)
        .collect();
    let k = rows.len();
    if k == 0 {
        return 1;
    }
    let widths: Vec<i64> = rows
        .iter()
        .map(|&r| eta[r].saturating_sub(strip) as i64)
        .collect();
    if widths.iter().any(|&w| w <= 0) {
        return 0;
    }

    let mut mat = vec![vec![0i128; k]; k];
    for a in 1..=k {
        for b in 1..=k {
            let width = widths[k - b];
            mat[a - 1][b - 1] = binom_i128(width - b as i64 + 1, a as i64 - b as i64 + 1);
        }
    }
    let det = det_i128(&mat);
    i64::try_from(det).expect("LGV determinant overflowed i64")
}

/// Number of set bits in a row mask.
pub fn row_mask_size(mask: u64) -> usize {
    mask.count_ones() as usize
}

/// Coefficient vector for `coeff * t^deg`.
pub fn monomial_coefficient_vector(deg: usize, coeff: i64) -> Vec<i64> {
    if coeff == 0 {
        return vec![0];
    }
    let mut out = vec![0; deg + 1];
    out[deg] = coeff;
    out
}

/// Fixed-row LGV count for the quotient term `(F_c-F_{c+1})/t`, indexed by
/// the remaining rows after deleting the forced boundary rook.
pub fn marked_delta_non_nesting_count_lgv(eta: &[usize], remaining_mask: u64, c: usize) -> i64 {
    assert!(
        eta.len() <= u64::BITS as usize,
        "row masks support at most 64 rows"
    );
    let highest_remaining_row = (0..eta.len())
        .filter(|&r| (remaining_mask & (1u64 << r)) != 0)
        .max();
    let first_forced_row = highest_remaining_row.map_or(0, |r| r + 1);
    let eligible_forced_rows = (first_forced_row..eta.len())
        .filter(|&r| eta[r] > c)
        .count() as i64;

    eligible_forced_rows * fixed_row_non_nesting_count_lgv(eta, c + 1, remaining_mask)
}

/// Add two coefficient vectors and trim trailing zeros.
pub fn add_coefficient_vectors(a: &[i64], b: &[i64]) -> Vec<i64> {
    let n = a.len().max(b.len());
    let mut r = vec![0; n];
    for (i, &c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        r[i] += c;
    }
    trim_coefficient_vector(r)
}

/// Subtract two coefficient vectors and trim trailing zeros.
pub fn subtract_coefficient_vectors(a: &[i64], b: &[i64]) -> Vec<i64> {
    let n = a.len().max(b.len());
    let mut r = vec![0; n];
    for (i, &c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        r[i] -= c;
    }
    trim_coefficient_vector(r)
}

/// Multiply a coefficient vector by `t`.
pub fn multiply_by_t(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0; p.len() + 1];
    for (i, &c) in p.iter().enumerate() {
        r[i + 1] = c;
    }
    trim_coefficient_vector(r)
}

/// Divide a coefficient vector by `t`, returning `None` if the constant term is nonzero.
pub fn divide_by_t(p: &[i64]) -> Option<Vec<i64>> {
    if p.first().copied().unwrap_or(0) != 0 {
        return None;
    }
    if p.len() <= 1 {
        return Some(vec![0]);
    }
    Some(trim_coefficient_vector(p[1..].to_vec()))
}

/// The adjacent packet indexed by remaining old rows after a boundary rook is
/// removed. Returns `(A_J, U_J, L_J)`.
pub fn fixed_remaining_non_nesting_packet_lgv(
    eta: &[usize],
    remaining_mask: u64,
    c: usize,
    p: usize,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    let k = row_mask_size(remaining_mask);
    let d_c1 = fixed_row_non_nesting_count_lgv(eta, c + 1, remaining_mask);
    let d_p1 = fixed_row_non_nesting_count_lgv(eta, p + 1, remaining_mask);
    let tail_sum: i64 = ((c + 2)..=(p + 1))
        .map(|j| fixed_row_non_nesting_count_lgv(eta, j, remaining_mask))
        .sum();
    let marked = marked_delta_non_nesting_count_lgv(eta, remaining_mask, c);

    let a_poly = add_coefficient_vectors(
        &monomial_coefficient_vector(k, d_c1),
        &monomial_coefficient_vector(k + 1, tail_sum),
    );
    let u_poly = monomial_coefficient_vector(k, marked + d_c1 + d_p1);
    let l_poly = monomial_coefficient_vector(k, marked + d_c1 - d_p1);

    (
        trim_coefficient_vector(a_poly),
        trim_coefficient_vector(u_poly),
        trim_coefficient_vector(l_poly),
    )
}

/// Global adjacent packet `(A,U,L)` for a Ferrers top shape `eta`.
pub fn adjacent_non_nesting_packet_global(
    eta: &[usize],
    c: usize,
    p: usize,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    assert!(!eta.is_empty(), "eta must have at least one row");
    assert!(c <= p, "expected c <= p");
    assert!(
        p < eta[0],
        "expected p to be smaller than the first row width"
    );

    let f: Vec<Vec<i64>> = (0..=eta[0])
        .map(|i| non_nesting_rook_polynomial(&strip_ferrers_columns(eta, i)))
        .collect();
    let delta =
        divide_by_t(&subtract_coefficient_vectors(&f[c], &f[c + 1])).expect("constant terms agree");
    let mut a_poly = f[c + 1].clone();
    for item in f.iter().take(p + 2).skip(c + 2) {
        a_poly = add_coefficient_vectors(&a_poly, &multiply_by_t(item));
    }
    let u_poly = add_coefficient_vectors(&add_coefficient_vectors(&delta, &f[c + 1]), &f[p + 1]);
    let l_poly =
        subtract_coefficient_vectors(&add_coefficient_vectors(&delta, &f[c + 1]), &f[p + 1]);
    (a_poly, u_poly, l_poly)
}

/// Degree with the zero polynomial convention used by the rook diagnostics:
/// `degree([0]) == 0`.
pub fn polynomial_degree_zero_convention(p: &[i64]) -> usize {
    p.iter().rposition(|&c| c != 0).unwrap_or(0)
}

/// Directed weak interlacing `f << g`, with zero-polynomial conventions used
/// by the non-nesting-rook diagnostics.
pub fn directed_weak_interlaces(f: &[i64], g: &[i64]) -> Option<bool> {
    let df = polynomial_degree_zero_convention(f);
    let dg = polynomial_degree_zero_convention(g);
    if df == 0 && f.first().copied().unwrap_or(0) == 0 {
        return Some(is_real_rooted(g));
    }
    if dg == 0 && g.first().copied().unwrap_or(0) == 0 {
        return Some(df == 0 && f.first().copied().unwrap_or(0) == 0);
    }
    if df == 0 && dg == 0 {
        return Some(true);
    }
    if df > dg + 1 || dg > df + 1 {
        return None;
    }
    if df <= dg {
        check_weak_interlacing(f, g)
    } else {
        check_weak_interlacing(g, f)
    }
}

/// Whether every coefficient is nonnegative.
pub fn has_nonnegative_coefficients(p: &[i64]) -> bool {
    p.iter().all(|&c| c >= 0)
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

    #[test]
    fn test_integer_partitions_and_column_stripping() {
        assert_eq!(
            integer_partitions(5),
            vec![
                vec![5],
                vec![4, 1],
                vec![3, 2],
                vec![3, 1, 1],
                vec![2, 2, 1],
                vec![2, 1, 1, 1],
                vec![1, 1, 1, 1, 1],
            ]
        );
        assert_eq!(strip_ferrers_columns(&[5, 3, 3, 1], 2), vec![3, 1, 1]);
        assert_eq!(strip_ferrers_columns(&[2, 1], 2), Vec::<usize>::new());
    }

    #[test]
    fn test_lgv_fixed_row_counts() {
        let eta = [3, 2, 1];
        assert_eq!(fixed_row_non_nesting_count_lgv(&eta, 0, 0), 1);
        assert_eq!(fixed_row_non_nesting_count_lgv(&eta, 0, 0b001), 3);
        assert_eq!(fixed_row_non_nesting_count_lgv(&eta, 0, 0b111), 1);
        assert_eq!(fixed_row_non_nesting_count_lgv(&eta, 2, 0b011), 0);
    }

    #[test]
    fn test_lgv_adjacent_packets() {
        let eta = [3, 2, 1];
        assert_eq!(marked_delta_non_nesting_count_lgv(&eta, 0b011, 0), 1);
        assert_eq!(
            fixed_remaining_non_nesting_packet_lgv(&eta, 0b011, 0, 1),
            (vec![0, 0, 1], vec![0, 0, 2], vec![0, 0, 2])
        );
        assert_eq!(
            adjacent_non_nesting_packet_global(&eta, 0, 1),
            (vec![1, 4, 2], vec![5, 9, 2], vec![3, 7, 2])
        );
    }

    #[test]
    fn test_rook_diagnostic_polynomial_helpers() {
        assert_eq!(add_coefficient_vectors(&[1, 2], &[0, 3, 0]), vec![1, 5]);
        assert_eq!(subtract_coefficient_vectors(&[1, 2], &[3, 1]), vec![-2, 1]);
        assert_eq!(multiply_by_t(&[1, 2]), vec![0, 1, 2]);
        assert_eq!(divide_by_t(&[0, 1, 2]), Some(vec![1, 2]));
        assert_eq!(divide_by_t(&[1, 2]), None);
        assert_eq!(polynomial_degree_zero_convention(&[0]), 0);
        assert!(has_nonnegative_coefficients(&[0, 1, 2]));
        assert!(!has_nonnegative_coefficients(&[1, -1]));
    }

    #[test]
    fn test_directed_weak_interlacing_zero_conventions() {
        assert_eq!(directed_weak_interlaces(&[0], &[1, 2, 1]), Some(true));
        assert_eq!(directed_weak_interlaces(&[1, 2, 1], &[0]), Some(false));
        assert_eq!(directed_weak_interlaces(&[1], &[2]), Some(true));
    }
}
