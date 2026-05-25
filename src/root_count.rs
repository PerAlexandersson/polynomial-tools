//! Exact real-root counting over the integers.
//!
//! This module is a lightweight integer alternative to the rational Sturm code
//! and the Bézout-matrix real-rootedness test.  It uses a primitive
//! pseudo-remainder Sturm sequence: every Euclidean remainder is computed by
//! pseudo-division and then made primitive.  That keeps intermediate coefficient
//! growth much smaller than a naive rational PRS, while avoiding the large exact
//! PSD matrices used by the Bézout criterion.
//!
//! The positive-coefficient path is intended for combinatorial polynomials.  If
//! all coefficients have one sign, roots can only be non-positive; after
//! removing powers of `t`, real-rootedness of `f(t)` is equivalent to all roots
//! of `f(-t)` being positive.
//!
//! Public functions in this module count distinct roots of the square-free
//! part.  Real-rootedness tests compare that count with the square-free degree,
//! so repeated real roots are handled without needing multiplicities from the
//! Sturm sequence.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

fn trim(mut p: Vec<BigInt>) -> Vec<BigInt> {
    while p.last().is_some_and(|c| c.is_zero()) {
        p.pop();
    }
    p
}

fn trim_slice(p: &[BigInt]) -> Vec<BigInt> {
    trim(p.to_vec())
}

fn degree(p: &[BigInt]) -> Option<usize> {
    p.iter().rposition(|c| !c.is_zero())
}

fn leading(p: &[BigInt]) -> BigInt {
    degree(p).map(|d| p[d].clone()).unwrap_or_else(BigInt::zero)
}

fn sign_i8(x: &BigInt) -> i8 {
    if x.is_positive() {
        1
    } else if x.is_negative() {
        -1
    } else {
        0
    }
}

fn content(p: &[BigInt]) -> BigInt {
    let mut g = BigInt::zero();
    for c in p {
        if !c.is_zero() {
            let a = c.abs();
            g = if g.is_zero() { a } else { g.gcd(&a) };
        }
    }
    g
}

fn primitive_keep_sign(p: Vec<BigInt>) -> Vec<BigInt> {
    let p = trim(p);
    let g = content(&p);
    if g.is_zero() || g.is_one() {
        p
    } else {
        trim(p.into_iter().map(|c| c / &g).collect())
    }
}

fn primitive_positive(p: Vec<BigInt>) -> Vec<BigInt> {
    let mut p = primitive_keep_sign(p);
    if leading(&p).is_negative() {
        for c in &mut p {
            *c = -c.clone();
        }
    }
    p
}

fn poly_neg(p: &[BigInt]) -> Vec<BigInt> {
    p.iter().map(|c| -c).collect()
}

fn poly_derivative(p: &[BigInt]) -> Vec<BigInt> {
    if p.len() <= 1 {
        return vec![];
    }
    trim(
        p.iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| c * BigInt::from(i))
            .collect(),
    )
}

fn poly_scale(p: &[BigInt], a: &BigInt) -> Vec<BigInt> {
    if a.is_zero() || p.is_empty() {
        return vec![];
    }
    trim(p.iter().map(|c| c * a).collect())
}

fn poly_sub_shifted_scaled(
    a: &[BigInt],
    b: &[BigInt],
    shift: usize,
    scale: &BigInt,
) -> Vec<BigInt> {
    let n = a.len().max(b.len() + shift);
    let mut r = vec![BigInt::zero(); n];
    for (i, c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, c) in b.iter().enumerate() {
        r[i + shift] -= c * scale;
    }
    trim(r)
}

/// Pseudo-remainder with a positive multiplier.
///
/// If `b_pos` is `b` multiplied by a sign so that its leading coefficient is
/// positive, this returns `lc(b_pos)^k * rem(a, b_pos)` for a nonnegative `k`.
/// The scalar is positive, so it does not change Sturm sign variations.
fn pseudo_remainder_positive_multiplier(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let da = match degree(a) {
        Some(d) => d,
        None => return vec![],
    };
    let db = degree(b).expect("pseudo_remainder: division by zero polynomial");
    if da < db {
        return trim_slice(a);
    }

    let mut b_pos = trim_slice(b);
    if leading(&b_pos).is_negative() {
        for c in &mut b_pos {
            *c = -c.clone();
        }
    }
    let lc = leading(&b_pos);
    debug_assert!(lc.is_positive());

    let mut r = trim_slice(a);
    let mut e = da - db + 1;
    while let Some(dr) = degree(&r) {
        if dr < db {
            break;
        }
        let c = leading(&r);
        let shift = dr - db;
        let scaled_r = poly_scale(&r, &lc);
        r = poly_sub_shifted_scaled(&scaled_r, &b_pos, shift, &c);
        e -= 1;
    }
    for _ in 0..e {
        r = poly_scale(&r, &lc);
    }
    trim(r)
}

fn primitive_gcd(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let mut r0 = primitive_positive(trim_slice(a));
    let mut r1 = primitive_positive(trim_slice(b));
    if r0.is_empty() {
        return r1;
    }
    if r1.is_empty() {
        return r0;
    }

    while !r1.is_empty() {
        let rem = primitive_positive(pseudo_remainder_positive_multiplier(&r0, &r1));
        r0 = r1;
        r1 = rem;
    }
    primitive_positive(r0)
}

fn exact_div(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let da = match degree(a) {
        Some(d) => d,
        None => return vec![],
    };
    let db = degree(b).expect("exact_div: division by zero polynomial");
    if da < db {
        return vec![];
    }

    let mut rem = trim_slice(a);
    let mut q = vec![BigInt::zero(); da - db + 1];
    let lc_b = leading(b);
    for k in (0..=da - db).rev() {
        let coeff = rem[db + k].clone();
        if coeff.is_zero() {
            continue;
        }
        assert!(
            (&coeff % &lc_b).is_zero(),
            "exact_div: non-exact leading quotient"
        );
        let qk = coeff / &lc_b;
        q[k] = qk.clone();
        for (j, bj) in b.iter().enumerate().take(db + 1) {
            rem[j + k] -= &qk * bj;
        }
    }
    debug_assert!(trim(rem).is_empty(), "exact_div: nonzero remainder");
    trim(q)
}

fn strip_initial_zeros_bigint(coeffs: &[BigInt]) -> &[BigInt] {
    match coeffs.iter().position(|c| !c.is_zero()) {
        Some(i) => &coeffs[i..],
        None => &[],
    }
}

fn squarefree_part_bigint(coeffs: &[BigInt]) -> Vec<BigInt> {
    let p = primitive_positive(trim_slice(coeffs));
    if degree(&p).unwrap_or(0) == 0 {
        return p;
    }
    let dp = poly_derivative(&p);
    if dp.is_empty() {
        return p;
    }
    let g = primitive_gcd(&p, &dp);
    if degree(&g).unwrap_or(0) == 0 {
        p
    } else {
        primitive_positive(exact_div(&p, &g))
    }
}

fn sturm_chain_squarefree(coeffs: &[BigInt]) -> Vec<Vec<BigInt>> {
    let p0 = squarefree_part_bigint(coeffs);
    if p0.is_empty() || degree(&p0).unwrap_or(0) == 0 {
        return vec![p0];
    }
    let p1 = primitive_keep_sign(poly_derivative(&p0));
    let mut chain = vec![p0, p1];

    loop {
        let n = chain.len();
        let prem = pseudo_remainder_positive_multiplier(&chain[n - 2], &chain[n - 1]);
        if prem.is_empty() {
            break;
        }
        let next = primitive_keep_sign(poly_neg(&primitive_keep_sign(prem)));
        if next.is_empty() {
            break;
        }
        chain.push(next);
    }
    chain
}

fn sign_at_pos_infinity(p: &[BigInt]) -> i8 {
    sign_i8(&leading(p))
}

fn sign_at_neg_infinity(p: &[BigInt]) -> i8 {
    let d = match degree(p) {
        Some(d) => d,
        None => return 0,
    };
    let s = sign_i8(&leading(p));
    if d % 2 == 0 {
        s
    } else {
        -s
    }
}

fn sign_at_zero_plus(p: &[BigInt]) -> i8 {
    for c in p {
        let s = sign_i8(c);
        if s != 0 {
            return s;
        }
    }
    0
}

fn sign_variations<I>(signs: I) -> usize
where
    I: IntoIterator<Item = i8>,
{
    let mut prev = 0i8;
    let mut changes = 0usize;
    for s in signs {
        if s == 0 {
            continue;
        }
        if prev != 0 && prev != s {
            changes += 1;
        }
        prev = s;
    }
    changes
}

fn variations_at_neg_infinity(chain: &[Vec<BigInt>]) -> usize {
    sign_variations(chain.iter().map(|p| sign_at_neg_infinity(p)))
}

fn variations_at_pos_infinity(chain: &[Vec<BigInt>]) -> usize {
    sign_variations(chain.iter().map(|p| sign_at_pos_infinity(p)))
}

fn variations_at_zero_plus(chain: &[Vec<BigInt>]) -> usize {
    sign_variations(chain.iter().map(|p| sign_at_zero_plus(p)))
}

fn alternating_neg_argument(coeffs: &[BigInt]) -> Vec<BigInt> {
    trim(
        coeffs
            .iter()
            .enumerate()
            .map(|(i, c)| if i % 2 == 0 { c.clone() } else { -c })
            .collect(),
    )
}

fn all_one_sign(coeffs: &[BigInt]) -> bool {
    let mut seen_pos = false;
    let mut seen_neg = false;
    for c in coeffs {
        if c.is_positive() {
            seen_pos = true;
        } else if c.is_negative() {
            seen_neg = true;
        }
    }
    !(seen_pos && seen_neg)
}

fn make_nonnegative(mut coeffs: Vec<BigInt>) -> Vec<BigInt> {
    if coeffs.iter().any(|c| c.is_negative()) {
        for c in &mut coeffs {
            *c = -c.clone();
        }
    }
    coeffs
}

/// Degree of the square-free part of a `BigInt` polynomial.
///
/// The zero polynomial and nonzero constants both return `0`.
pub fn squarefree_degree_bigint_coeffs(coeffs: &[BigInt]) -> usize {
    degree(&squarefree_part_bigint(coeffs)).unwrap_or(0)
}

/// Count distinct real roots using a primitive pseudo-remainder Sturm sequence.
///
/// The input may have repeated roots; internally we replace it by its
/// square-free part before counting.
pub fn count_real_roots_prs_bigint_coeffs(coeffs: &[BigInt]) -> usize {
    let chain = sturm_chain_squarefree(coeffs);
    if chain.is_empty() {
        return 0;
    }
    variations_at_neg_infinity(&chain).saturating_sub(variations_at_pos_infinity(&chain))
}

/// Count distinct positive roots using a primitive pseudo-remainder Sturm sequence.
///
/// The count is over the open interval `(0, +infinity)`.
pub fn count_positive_roots_prs_bigint_coeffs(coeffs: &[BigInt]) -> usize {
    let chain = sturm_chain_squarefree(coeffs);
    if chain.is_empty() {
        return 0;
    }
    variations_at_zero_plus(&chain).saturating_sub(variations_at_pos_infinity(&chain))
}

/// Check Newton's inequalities for a nonnegative coefficient sequence.
///
/// This is a cheap necessary condition for real-rootedness.  It is exact and
/// works with arbitrary-size coefficients.
pub fn satisfies_newton_inequalities_bigint(coeffs: &[BigInt]) -> bool {
    let p = trim_slice(strip_initial_zeros_bigint(coeffs));
    let d = match degree(&p) {
        Some(d) if d >= 2 => d,
        _ => return true,
    };
    if p.iter().any(|c| c.is_negative()) {
        return false;
    }

    let mut binom = vec![BigInt::one(); d + 1];
    for k in 1..=d {
        binom[k] = &binom[k - 1] * BigInt::from(d - k + 1) / BigInt::from(k);
    }

    for k in 1..d {
        let lhs = p[k].pow(2) * &binom[k - 1] * &binom[k + 1];
        let rhs = &p[k - 1] * &p[k + 1] * &binom[k] * &binom[k];
        if lhs < rhs {
            return false;
        }
    }
    true
}

/// Exact real-rootedness via primitive integer Sturm/PRS root counting.
///
/// This is a general method: it works for mixed-sign coefficients as well as
/// one-signed combinatorial polynomials.  It is often useful as a fallback when
/// Bézout matrices become large.
pub fn is_real_rooted_prs_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let p = trim_slice(coeffs);
    let d = match degree(&p) {
        Some(d) => d,
        None => return true,
    };
    if d <= 1 {
        return true;
    }
    count_real_roots_prs_bigint_coeffs(&p) == squarefree_degree_bigint_coeffs(&p)
}

/// Exact real-rootedness optimized for one-signed coefficient polynomials.
///
/// Returns `None` if the nonzero coefficients do not all have the same sign.
/// Powers of `t` are removed first, since zero is already a real root.
pub fn is_real_rooted_one_signed_bigint_coeffs(coeffs: &[BigInt]) -> Option<bool> {
    let p = trim_slice(strip_initial_zeros_bigint(coeffs));
    let d = match degree(&p) {
        Some(d) => d,
        None => return Some(true),
    };
    if d <= 1 {
        return Some(true);
    }
    if !all_one_sign(&p) {
        return None;
    }

    let p = make_nonnegative(p);
    if !satisfies_newton_inequalities_bigint(&p) {
        return Some(false);
    }

    let transformed = alternating_neg_argument(&p);
    let positive_roots = count_positive_roots_prs_bigint_coeffs(&transformed);
    let sf_degree = squarefree_degree_bigint_coeffs(&p);
    Some(positive_roots == sf_degree)
}

/// Exact real-rootedness with a fast one-signed path and PRS fallback.
///
/// This avoids constructing Bézout/PSD matrices entirely.  The public
/// `real_rootedness::is_real_rooted_bigint_coeffs` function currently uses the
/// one-signed part of this path and keeps Bézout as the mixed-sign fallback,
/// because Bézout also supports the interlacing routines.
pub fn is_real_rooted_fast_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    if let Some(rr) = is_real_rooted_one_signed_bigint_coeffs(coeffs) {
        rr
    } else {
        is_real_rooted_prs_bigint_coeffs(coeffs)
    }
}

/// Convenience wrapper for `i64` coefficient vectors.
pub fn is_real_rooted_fast_i64(coeffs: &[i64]) -> bool {
    let coeffs: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_real_rooted_fast_bigint_coeffs(&coeffs)
}

/// The largest coefficient bit-size seen in a primitive Sturm/PRS chain.
///
/// This is useful for benchmarking and diagnosing coefficient swell.
pub fn primitive_sturm_max_coefficient_bits(coeffs: &[BigInt]) -> u64 {
    sturm_chain_squarefree(coeffs)
        .iter()
        .flat_map(|p| p.iter())
        .map(|c| c.bits())
        .max()
        .unwrap_or(0)
}

/// Try to reduce a BigInt vector to i64 coefficients.
pub fn bigint_coeffs_to_i64(coeffs: &[BigInt]) -> Option<Vec<i64>> {
    coeffs.iter().map(ToPrimitive::to_i64).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b(v: &[i64]) -> Vec<BigInt> {
        v.iter().map(|&x| BigInt::from(x)).collect()
    }

    fn mul_i64(a: &[i64], b: &[i64]) -> Vec<i64> {
        let mut r = vec![0i64; a.len() + b.len() - 1];
        for (i, &ca) in a.iter().enumerate() {
            for (j, &cb) in b.iter().enumerate() {
                r[i + j] += ca * cb;
            }
        }
        while r.last() == Some(&0) {
            r.pop();
        }
        r
    }

    #[test]
    fn test_pseudo_remainder_basic() {
        // (x^3 - 1) mod (x - 1) = 0.
        let r = pseudo_remainder_positive_multiplier(&b(&[-1, 0, 0, 1]), &b(&[-1, 1]));
        assert!(r.is_empty());

        // (x^2 + 1) mod (x - 2) = 5.  The divisor is monic, so prem = rem.
        let r = pseudo_remainder_positive_multiplier(&b(&[1, 0, 1]), &b(&[-2, 1]));
        assert_eq!(r, b(&[5]));
    }

    #[test]
    fn test_squarefree_part_and_degree() {
        // (1+x)^3 has square-free part 1+x.
        assert_eq!(squarefree_degree_bigint_coeffs(&b(&[1, 3, 3, 1])), 1);

        // (x-1)^2 (x-2) = -2 + 5x - 4x^2 + x^3, square-free degree 2.
        assert_eq!(squarefree_degree_bigint_coeffs(&b(&[-2, 5, -4, 1])), 2);
    }

    #[test]
    fn test_count_real_roots() {
        assert_eq!(count_real_roots_prs_bigint_coeffs(&b(&[1, 0, 1])), 0);
        assert_eq!(count_real_roots_prs_bigint_coeffs(&b(&[-6, 11, -6, 1])), 3);
        assert_eq!(count_real_roots_prs_bigint_coeffs(&b(&[1, 3, 3, 1])), 1);
        assert_eq!(count_real_roots_prs_bigint_coeffs(&b(&[0, -1, 0, 1])), 3);
    }

    #[test]
    fn test_count_positive_roots() {
        // (x-1)(x-2)(x+3) = 6 - 7x + x^3 has two positive roots.
        assert_eq!(
            count_positive_roots_prs_bigint_coeffs(&b(&[6, -7, 0, 1])),
            2
        );
        assert_eq!(count_positive_roots_prs_bigint_coeffs(&b(&[1, 2, 1])), 0);
    }

    #[test]
    fn test_one_signed_real_rootedness() {
        assert_eq!(
            is_real_rooted_one_signed_bigint_coeffs(&b(&[1, 4, 6, 4, 1])),
            Some(true)
        );
        assert_eq!(
            is_real_rooted_one_signed_bigint_coeffs(&b(&[1, 43, 196, 168, 23, 1])),
            Some(false)
        );
        assert_eq!(
            is_real_rooted_one_signed_bigint_coeffs(&b(&[1, -2, 1])),
            None
        );
    }

    #[test]
    fn test_newton_filter() {
        assert!(satisfies_newton_inequalities_bigint(&b(&[1, 4, 6, 4, 1])));
        assert!(!satisfies_newton_inequalities_bigint(&b(&[1, 1, 10, 1])));
    }

    #[test]
    fn test_prs_real_rootedness_general() {
        assert!(is_real_rooted_prs_bigint_coeffs(&b(&[-6, 11, -6, 1])));
        assert!(!is_real_rooted_prs_bigint_coeffs(&b(&[1, 0, 1])));
        assert!(is_real_rooted_prs_bigint_coeffs(&b(&[0, 0, 1, 2, 1]))); // x^2(1+x)^2
    }

    #[test]
    fn test_agrees_with_existing_examples() {
        let cases = [
            (vec![1, 2, 1], true),
            (vec![1, 11, 11, 1], true),
            (vec![1, 0, 1], false),
            (vec![1, 43, 196, 168, 23, 1], false),
            (vec![-15, 23, -9, 1], true),
        ];
        for (coeffs, expected) in cases {
            assert_eq!(
                is_real_rooted_fast_bigint_coeffs(&b(&coeffs)),
                expected,
                "coeffs={coeffs:?}"
            );
        }
    }

    #[test]
    fn test_products_with_known_real_roots() {
        let mut p = vec![1i64];
        for a in 1..=8 {
            p = mul_i64(&p, &[a, 1]); // (x+a)
            assert!(is_real_rooted_fast_bigint_coeffs(&b(&p)));
        }
        let complex_factor = mul_i64(&p, &[1, 0, 1]);
        assert!(!is_real_rooted_fast_bigint_coeffs(&b(&complex_factor)));
    }

    #[test]
    fn test_prs_agrees_with_rational_sturm_on_small_grid() {
        for degree in 0usize..=4 {
            let total = 5usize.pow((degree + 1) as u32);
            for mut mask in 0..total {
                let mut coeffs = Vec::with_capacity(degree + 1);
                for _ in 0..=degree {
                    coeffs.push((mask % 5) as i64 - 2);
                    mask /= 5;
                }

                let bigint = b(&coeffs);
                assert_eq!(
                    is_real_rooted_prs_bigint_coeffs(&bigint),
                    crate::real_rootedness::is_real_rooted_sturm_bigint_coeffs(&bigint),
                    "coeffs={coeffs:?}"
                );
            }
        }
    }
}
