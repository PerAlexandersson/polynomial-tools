//! Linear algebra over exact rationals (BigInt/BigRational).
//!
//! Provides fraction-free Bareiss elimination, Gaussian elimination, positive
//! definiteness checks, determinants, linear system solving, and total
//! non-negativity (TNN) checking.
//!
//! # Total non-negativity
//!
//! A matrix is *totally non-negative* (TNN) if every minor (determinant of every
//! square submatrix) is non-negative. The brute-force check
//! ([`check_total_positivity`]) enumerates all submatrices up to a given size,
//! which is exponential. The Neville elimination method ([`check_tnn_neville`])
//! checks TNN for *all* minors in O(n²m) time by performing adjacent-row
//! subtraction and verifying that no entry becomes negative.

#![allow(clippy::needless_range_loop)]

use crate::Polynomial;
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, ToPrimitive, Zero};

type Q = Ratio<BigInt>;

/// Matrix dimension at which the default positive-definiteness check switches
/// from BigInt Bareiss to modular CRT reconstruction.
pub const MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD: usize = 30;

/// Check if a symmetric BigInt matrix is positive definite.
///
/// Uses Sylvester's criterion.  Small matrices use fraction-free BigInt
/// Bareiss elimination; matrices of size
/// [`MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD`] or larger try the modular
/// CRT path first and fall back to Bareiss if the CRT bounds do not certify in
/// time. The explicit [`is_positive_definite_bareiss`] and
/// [`is_positive_definite_modular`] entry points are kept for benchmarking.
pub fn is_positive_definite(mat: &[Vec<BigInt>]) -> bool {
    if !is_symmetric_bigint_matrix(mat) {
        return false;
    }
    if mat.len() >= MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD {
        is_positive_definite_modular(mat).unwrap_or_else(|| is_positive_definite_bareiss(mat))
    } else {
        is_positive_definite_bareiss(mat)
    }
}

/// Check positive definiteness via BigInt Bareiss leading principal minors.
pub fn is_positive_definite_bareiss(mat: &[Vec<BigInt>]) -> bool {
    if !is_symmetric_bigint_matrix(mat) {
        return false;
    }
    bareiss_leading_principal_minors_bigint(mat)
        .map(|minors| minors.into_iter().all(|minor| minor > BigInt::zero()))
        .unwrap_or(false)
}

/// Try to check positive definiteness via modular CRT.
///
/// This is a second implementation of Sylvester's criterion.  It reconstructs
/// the leading principal minors from determinants computed over prime fields,
/// using Hadamard bounds to certify when CRT reconstruction is exact.
///
/// Returns `None` when the matrix is not square or when the accumulated CRT
/// modulus does not certify all leading minors within the prime budget.
pub fn is_positive_definite_modular(mat: &[Vec<BigInt>]) -> Option<bool> {
    if !is_square_bigint_matrix(mat) {
        return None;
    }
    if !is_symmetric_bigint_matrix(mat) {
        return Some(false);
    }
    modular_leading_principal_minors_bigint(mat)
        .map(|minors| minors.into_iter().all(|minor| minor > BigInt::zero()))
}

/// Check if a symmetric BigInt matrix is positive semi-definite.
///
/// A symmetric matrix is positive semi-definite iff all diagonal pivots
/// are non-negative, and zero pivots have zero entries in their column below.
pub fn is_positive_semidefinite(mat: &[Vec<BigInt>]) -> bool {
    let n = mat.len();
    if !is_symmetric_bigint_matrix(mat) {
        return false;
    }
    let qmat = bigint_to_q(mat);
    let mut a = qmat;

    for k in 0..n {
        if a[k][k] < Q::zero() {
            return false;
        }
        if a[k][k].is_zero() {
            // All entries below must also be zero
            for i in (k + 1)..n {
                if !a[i][k].is_zero() {
                    return false;
                }
            }
            continue;
        }
        let pivot = a[k][k].clone();
        for i in (k + 1)..n {
            if a[i][k].is_zero() {
                continue;
            }
            let factor = a[i][k].clone() / pivot.clone();
            for j in (k + 1)..n {
                let sub = factor.clone() * a[k][j].clone();
                a[i][j] -= sub;
            }
            a[i][k] = Q::zero();
        }
    }
    true
}

/// Compute the determinant of a BigInt matrix using fraction-free elimination.
pub fn determinant(mat: &[Vec<BigInt>]) -> BigInt {
    assert!(
        is_square_bigint_matrix(mat),
        "determinant requires a square matrix"
    );
    let n = mat.len();
    if n == 0 {
        return BigInt::one();
    }

    let mut a = mat.to_vec();
    let mut denom = BigInt::one();
    let mut sign_is_negative = false;

    for k in 0..(n - 1) {
        if a[k][k].is_zero() {
            let Some(pivot_row) = ((k + 1)..n).find(|&i| !a[i][k].is_zero()) else {
                return BigInt::zero();
            };
            a.swap(k, pivot_row);
            sign_is_negative = !sign_is_negative;
        }

        let pivot = a[k][k].clone();
        let mut next_entries = Vec::with_capacity((n - k - 1) * (n - k - 1));
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let numerator = &a[i][j] * &pivot - &a[i][k] * &a[k][j];
                debug_assert!(
                    (&numerator % &denom).is_zero(),
                    "Bareiss division should be exact"
                );
                next_entries.push((i, j, numerator / &denom));
            }
        }
        for (i, j, value) in next_entries {
            a[i][j] = value;
        }
        denom = pivot;
    }

    if sign_is_negative {
        -a[n - 1][n - 1].clone()
    } else {
        a[n - 1][n - 1].clone()
    }
}

fn is_square_bigint_matrix(mat: &[Vec<BigInt>]) -> bool {
    let n = mat.len();
    mat.iter().all(|row| row.len() == n)
}

fn is_symmetric_bigint_matrix(mat: &[Vec<BigInt>]) -> bool {
    let n = mat.len();
    if !is_square_bigint_matrix(mat) {
        return false;
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if mat[i][j] != mat[j][i] {
                return false;
            }
        }
    }
    true
}

/// Compute the leading principal determinants by fraction-free Bareiss
/// elimination over `BigInt`, without row swaps.
///
/// If every Bareiss pivot is nonzero, the `k`th returned value is the
/// determinant of the leading `(k+1) x (k+1)` principal submatrix.  This is
/// the exact integer version of the fixed-order elimination used in
/// Sylvester's criterion.
///
/// Returns `None` if the matrix is not square, a nonfinal required pivot is
/// zero, or an exact Bareiss division unexpectedly fails.  A nonfinal zero
/// pivot does not imply that the whole matrix is singular; it only means this
/// no-pivot leading-principal computation cannot continue.
pub fn bareiss_leading_principal_minors_bigint(mat: &[Vec<BigInt>]) -> Option<Vec<BigInt>> {
    let n = mat.len();
    if mat.iter().any(|row| row.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(Vec::new());
    }

    let mut a = mat.to_vec();
    let mut denom = BigInt::one();
    let mut minors = Vec::with_capacity(n);

    for k in 0..n {
        let pivot = a[k][k].clone();
        if pivot.is_zero() {
            if k == n - 1 {
                minors.push(pivot);
                break;
            }
            return None;
        }
        minors.push(pivot.clone());
        if k == n - 1 {
            break;
        }

        let mut next_entries = Vec::with_capacity((n - k - 1) * (n - k - 1));
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let numerator = &a[i][j] * &pivot - &a[i][k] * &a[k][j];
                if &numerator % &denom != BigInt::zero() {
                    return None;
                }
                next_entries.push((i, j, numerator / &denom));
            }
        }
        for (i, j, value) in next_entries {
            a[i][j] = value;
        }
        denom = pivot;
    }

    Some(minors)
}

/// Compute a determinant by no-pivot fraction-free Bareiss elimination.
///
/// This returns `None` in the same cases as
/// [`bareiss_leading_principal_minors_bigint`].  Use [`determinant`] when a
/// row-swapping determinant routine is needed.
pub fn bareiss_determinant_bigint(mat: &[Vec<BigInt>]) -> Option<BigInt> {
    if mat.is_empty() {
        return Some(BigInt::one());
    }
    bareiss_leading_principal_minors_bigint(mat).and_then(|mut minors| minors.pop())
}

/// Compute leading principal determinants by modular determinant computation
/// and Chinese remaindering.
///
/// For a square integer matrix, this returns the determinant of each leading
/// principal submatrix.  Each determinant is computed modulo a sequence of
/// large primes and reconstructed by CRT once the accumulated modulus exceeds
/// twice a Hadamard bound for every leading minor.
///
/// Unlike [`bareiss_leading_principal_minors_bigint`], this method can handle
/// zero Bareiss pivots because each leading determinant is computed separately
/// with modular pivoting.  It is intended as a comparison/large-entry path, not
/// as a replacement for the simpler Bareiss implementation.
pub fn modular_leading_principal_minors_bigint(mat: &[Vec<BigInt>]) -> Option<Vec<BigInt>> {
    const MAX_MODULAR_PRIMES: usize = 256;

    let n = mat.len();
    if mat.iter().any(|row| row.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(Vec::new());
    }

    let hadamard_squared_bounds = leading_principal_hadamard_squared_bounds(mat)?;
    let mut residues = vec![BigInt::zero(); n];
    let mut modulus = BigInt::one();
    let mut prime_search_start = (1u64 << 61) - 1;

    for _ in 0..MAX_MODULAR_PRIMES {
        let prime = previous_prime_at_or_below(prime_search_start)?;
        prime_search_start = prime.saturating_sub(2);
        let modular_minors = leading_principal_minors_mod_prime(mat, prime);
        crt_update_residues(&mut residues, &mut modulus, &modular_minors, prime)?;

        if crt_modulus_certifies_bounds(&modulus, &hadamard_squared_bounds) {
            return Some(
                residues
                    .into_iter()
                    .map(|r| symmetric_residue(r, &modulus))
                    .collect(),
            );
        }
    }

    None
}

fn leading_principal_hadamard_squared_bounds(mat: &[Vec<BigInt>]) -> Option<Vec<BigInt>> {
    let n = mat.len();
    let mut bounds = Vec::with_capacity(n);
    for k in 1..=n {
        let mut product = BigInt::one();
        for i in 0..k {
            let mut row_square_sum = BigInt::zero();
            for j in 0..k {
                row_square_sum += &mat[i][j] * &mat[i][j];
            }
            product *= row_square_sum;
        }
        bounds.push(product);
    }
    Some(bounds)
}

fn crt_modulus_certifies_bounds(modulus: &BigInt, squared_bounds: &[BigInt]) -> bool {
    let modulus_squared = modulus * modulus;
    squared_bounds.iter().all(|bound| {
        let threshold = BigInt::from(4) * bound;
        modulus_squared > threshold
    })
}

fn symmetric_residue(residue: BigInt, modulus: &BigInt) -> BigInt {
    if &residue * BigInt::from(2) > *modulus {
        residue - modulus
    } else {
        residue
    }
}

fn crt_update_residues(
    residues: &mut [BigInt],
    modulus: &mut BigInt,
    modular_values: &[u64],
    prime: u64,
) -> Option<()> {
    let modulus_mod_prime = bigint_mod_u64(modulus, prime);
    let inverse = mod_inverse_prime(modulus_mod_prime, prime)?;
    let old_modulus = modulus.clone();

    for (residue, &value_mod_prime) in residues.iter_mut().zip(modular_values.iter()) {
        let residue_mod_prime = bigint_mod_u64(residue, prime);
        let correction = if value_mod_prime >= residue_mod_prime {
            value_mod_prime - residue_mod_prime
        } else {
            prime - (residue_mod_prime - value_mod_prime)
        };
        let t = mul_mod_u64(correction, inverse, prime);
        *residue += &old_modulus * BigInt::from(t);
    }

    *modulus = old_modulus * BigInt::from(prime);
    Some(())
}

fn bigint_mod_u64(value: &BigInt, modulus: u64) -> u64 {
    let modulus_big = BigInt::from(modulus);
    let mut reduced = value % &modulus_big;
    if reduced.is_negative() {
        reduced += &modulus_big;
    }
    reduced.to_u64().expect("reduced residue should fit in u64")
}

fn leading_principal_minors_mod_prime(mat: &[Vec<BigInt>], prime: u64) -> Vec<u64> {
    let reduced: Vec<Vec<u64>> = mat
        .iter()
        .map(|row| {
            row.iter()
                .map(|entry| bigint_mod_u64(entry, prime))
                .collect()
        })
        .collect();
    if let Some(minors) = leading_principal_minors_mod_prime_bareiss(&reduced, prime) {
        return minors;
    }

    let n = mat.len();
    (1..=n)
        .map(|prefix_size| determinant_mod_prime_prefix(&reduced, prefix_size, prime))
        .collect()
}

fn leading_principal_minors_mod_prime_bareiss(mat: &[Vec<u64>], prime: u64) -> Option<Vec<u64>> {
    let n = mat.len();
    let mut a = mat.to_vec();
    let mut denom = 1u64;
    let mut minors = Vec::with_capacity(n);

    for k in 0..n {
        let pivot = a[k][k];
        if pivot == 0 {
            if k == n - 1 {
                minors.push(0);
                break;
            }
            return None;
        }
        minors.push(pivot);
        if k == n - 1 {
            break;
        }

        let denom_inverse = mod_inverse_prime(denom, prime)?;
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let left = mul_mod_u64(a[i][j], pivot, prime);
                let right = mul_mod_u64(a[i][k], a[k][j], prime);
                let numerator = sub_mod_u64(left, right, prime);
                a[i][j] = mul_mod_u64(numerator, denom_inverse, prime);
            }
        }
        denom = pivot;
    }

    Some(minors)
}

fn determinant_mod_prime_prefix(mat: &[Vec<u64>], size: usize, prime: u64) -> u64 {
    let mut a: Vec<Vec<u64>> = mat
        .iter()
        .take(size)
        .map(|row| row.iter().take(size).copied().collect())
        .collect();
    let mut det = 1u64;
    let mut sign_is_negative = false;

    for k in 0..size {
        let Some(pivot_row) = (k..size).find(|&i| a[i][k] != 0) else {
            return 0;
        };
        if pivot_row != k {
            a.swap(k, pivot_row);
            sign_is_negative = !sign_is_negative;
        }
        let pivot = a[k][k];
        det = mul_mod_u64(det, pivot, prime);
        let inverse = mod_inverse_prime(pivot, prime).expect("nonzero pivot modulo prime");
        for i in (k + 1)..size {
            if a[i][k] == 0 {
                continue;
            }
            let factor = mul_mod_u64(a[i][k], inverse, prime);
            for j in (k + 1)..size {
                let sub = mul_mod_u64(factor, a[k][j], prime);
                a[i][j] = sub_mod_u64(a[i][j], sub, prime);
            }
        }
    }

    if sign_is_negative && det != 0 {
        prime - det
    } else {
        det
    }
}

fn mul_mod_u64(a: u64, b: u64, modulus: u64) -> u64 {
    ((a as u128 * b as u128) % modulus as u128) as u64
}

fn sub_mod_u64(a: u64, b: u64, modulus: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        modulus - (b - a)
    }
}

fn pow_mod_u64(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1u64;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = mul_mod_u64(result, base, modulus);
        }
        base = mul_mod_u64(base, base, modulus);
        exponent >>= 1;
    }
    result
}

fn mod_inverse_prime(value: u64, prime: u64) -> Option<u64> {
    (value != 0).then(|| pow_mod_u64(value, prime - 2, prime))
}

fn previous_prime_at_or_below(mut n: u64) -> Option<u64> {
    if n < 2 {
        return None;
    }
    if n == 2 {
        return Some(2);
    }
    if n.is_multiple_of(2) {
        n -= 1;
    }
    while n >= 3 {
        if is_prime_u64(n) {
            return Some(n);
        }
        n = n.saturating_sub(2);
    }
    None
}

fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n == p {
            return true;
        }
        if n.is_multiple_of(p) {
            return false;
        }
    }

    let mut d = n - 1;
    let mut s = 0;
    while d.is_multiple_of(2) {
        d /= 2;
        s += 1;
    }

    for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if a >= n {
            continue;
        }
        let mut x = pow_mod_u64(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        let mut probably_prime_for_base = false;
        for _ in 1..s {
            x = mul_mod_u64(x, x, n);
            if x == n - 1 {
                probably_prime_for_base = true;
                break;
            }
        }
        if !probably_prime_for_base {
            return false;
        }
    }
    true
}

/// Compute the leading principal determinants by fraction-free Bareiss
/// elimination over matrices with entries in `Z[t]`.
///
/// The entries are [`Polynomial<BigInt>`] in ascending degree order.  All
/// divisions are checked as exact polynomial divisions over `Z[t]`.  This is
/// useful when a Bézout or Sylvester matrix depends polynomially on a parameter
/// and one wants exact leading principal minors without passing through
/// rational functions.
///
/// Returns `None` if the matrix is not square, a nonfinal required pivot is
/// zero, or an exact polynomial division fails.
pub fn bareiss_leading_principal_minors_polynomial_bigint(
    mat: &[Vec<Polynomial<BigInt>>],
) -> Option<Vec<Polynomial<BigInt>>> {
    let n = mat.len();
    if mat.iter().any(|row| row.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(Vec::new());
    }

    let mut a = mat.to_vec();
    let mut denom = Polynomial::<BigInt>::one();
    let mut minors = Vec::with_capacity(n);

    for k in 0..n {
        let pivot = a[k][k].clone();
        if pivot.is_zero() {
            if k == n - 1 {
                minors.push(pivot);
                break;
            }
            return None;
        }
        minors.push(pivot.clone());
        if k == n - 1 {
            break;
        }

        let mut next_entries = Vec::with_capacity((n - k - 1) * (n - k - 1));
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let numerator = a[i][j].clone() * pivot.clone() - a[i][k].clone() * a[k][j].clone();
                let value = exact_div_polynomial_bigint(&numerator, &denom)?;
                next_entries.push((i, j, value));
            }
        }
        for (i, j, value) in next_entries {
            a[i][j] = value;
        }
        denom = pivot;
    }

    Some(minors)
}

/// Compute a determinant in `Z[t]` by no-pivot fraction-free Bareiss
/// elimination.
///
/// Use this when the fixed pivot order is part of the certificate.  The return
/// value is `None` under the same conditions as
/// [`bareiss_leading_principal_minors_polynomial_bigint`].
pub fn bareiss_determinant_polynomial_bigint(
    mat: &[Vec<Polynomial<BigInt>>],
) -> Option<Polynomial<BigInt>> {
    if mat.is_empty() {
        return Some(Polynomial::one());
    }
    bareiss_leading_principal_minors_polynomial_bigint(mat).and_then(|mut minors| minors.pop())
}

fn exact_div_polynomial_bigint(
    numerator: &Polynomial<BigInt>,
    denominator: &Polynomial<BigInt>,
) -> Option<Polynomial<BigInt>> {
    if denominator.is_zero() {
        return None;
    }
    if numerator.is_zero() {
        return Some(Polynomial::zero());
    }

    let den = denominator.coeffs();
    let den_degree = den.len() - 1;
    let den_lc = den.last().expect("nonzero denominator");
    let mut rem = numerator.coeffs().to_vec();

    if rem.len() < den.len() {
        return if rem.iter().all(BigInt::is_zero) {
            Some(Polynomial::zero())
        } else {
            None
        };
    }

    let mut quotient = vec![BigInt::zero(); rem.len() - den_degree];
    trim_bigint_coeffs(&mut rem);
    while rem.len() >= den.len() && !rem.is_empty() {
        let shift = rem.len() - den.len();
        let rem_lc = rem.last().expect("nonempty remainder").clone();
        if (&rem_lc % den_lc) != BigInt::zero() {
            return None;
        }
        let q = rem_lc / den_lc;
        quotient[shift] = q.clone();
        for (i, c) in den.iter().enumerate() {
            rem[shift + i] -= &q * c;
        }
        trim_bigint_coeffs(&mut rem);
    }

    if rem.is_empty() {
        Some(Polynomial::new(quotient))
    } else {
        None
    }
}

/// Solve Ax = b via Gaussian elimination with full pivoting over ℚ.
///
/// Returns `None` if the system is inconsistent.
/// For underdetermined systems, free variables are set to zero.
pub fn solve_linear_system(a: &[Vec<Q>], b: &[Q]) -> Option<Vec<Q>> {
    let num_rows = a.len();
    if num_rows == 0 {
        return Some(vec![]);
    }
    let num_cols = a[0].len();

    // Build augmented matrix [A | b].
    let mut aug: Vec<Vec<Q>> = a
        .iter()
        .zip(b.iter())
        .map(|(row, bi)| {
            let mut r = row.clone();
            r.push(bi.clone());
            r
        })
        .collect();

    // Forward + upward elimination (reduced row-echelon form).
    let mut pivot_cols: Vec<(usize, usize)> = Vec::new();
    let mut pivot_row = 0;

    for col in 0..num_cols {
        let Some(pr) = (pivot_row..num_rows).find(|&r| !aug[r][col].is_zero()) else {
            continue;
        };

        aug.swap(pivot_row, pr);
        pivot_cols.push((pivot_row, col));

        let pivot_val = aug[pivot_row][col].clone();

        // Eliminate all other rows in this column.
        for row in 0..num_rows {
            if row == pivot_row || aug[row][col].is_zero() {
                continue;
            }
            let factor = aug[row][col].clone() / pivot_val.clone();
            let pivot_snapshot: Vec<_> = aug[pivot_row][col..=num_cols].to_vec();
            for (aug_j, pivot_j) in aug[row][col..=num_cols]
                .iter_mut()
                .zip(pivot_snapshot.iter())
            {
                let sub = pivot_j.clone() * &factor;
                *aug_j -= sub;
            }
        }

        pivot_row += 1;
    }

    // Consistency check: any row [0 ... 0 | nonzero] means no solution.
    for aug_row in &aug[pivot_row..num_rows] {
        if !aug_row[num_cols].is_zero() {
            return None;
        }
    }

    // Extract solution (free variables = 0).
    let mut x = vec![Q::zero(); num_cols];
    for &(pr, pc) in &pivot_cols {
        x[pc] = aug[pr][num_cols].clone() / aug[pr][pc].clone();
    }
    Some(x)
}

// ---------------------------------------------------------------------------
// Total positivity
// ---------------------------------------------------------------------------

/// Check if an integer matrix is totally positive (all minors strictly positive,
/// except those that are trivially zero due to zero rows/columns).
///
/// More precisely, checks that all minors are **non-negative** (totally non-negative, TNN).
/// Set `strict` to require all non-trivially-zero minors to be strictly positive.
///
/// `max_minor_size` limits the largest minor checked. For a full check, pass
/// `min(nrows, ncols)`. Returns `Ok(())` if all minors pass, or `Err(msg)`
/// describing the first violation found.
///
/// # Example
/// ```
/// use polynomial_tools::linalg::check_total_positivity;
/// // Pascal's triangle is totally positive
/// let mat = vec![
///     vec![1, 0, 0],
///     vec![1, 1, 0],
///     vec![1, 2, 1],
///     vec![1, 3, 3],
///     vec![1, 4, 6],
/// ];
/// assert!(check_total_positivity(&mat, 3, false).is_ok());
/// ```
pub fn check_total_positivity(
    mat: &[Vec<i64>],
    max_minor_size: usize,
    strict: bool,
) -> Result<(), String> {
    let nrows = mat.len();
    if nrows == 0 {
        return Ok(());
    }
    let ncols = mat[0].len();
    let max_k = max_minor_size.min(nrows).min(ncols);

    for k in 1..=max_k {
        let row_combos = combinations_usize(nrows, k);
        let col_combos = combinations_usize(ncols, k);
        for rows in &row_combos {
            for cols in &col_combos {
                let sub = extract_submatrix_i64(mat, rows, cols);
                let det = determinant(&sub);
                if det < BigInt::zero() {
                    return Err(format!(
                        "{}x{} minor rows {:?} cols {:?} has det = {} < 0",
                        k, k, rows, cols, det
                    ));
                }
                if strict && det.is_zero() {
                    // Check if this is a non-trivial zero
                    let all_zero = rows.iter().all(|&r| cols.iter().all(|&c| mat[r][c] == 0));
                    if !all_zero {
                        return Err(format!(
                            "{}x{} minor rows {:?} cols {:?} has det = 0 (strict mode)",
                            k, k, rows, cols
                        ));
                    }
                }
            }
        }
    }
    Ok(())
}

/// Check total non-negativity: all minors >= 0.
/// Convenience wrapper around [`check_total_positivity`].
pub fn is_totally_nonnegative(mat: &[Vec<i64>], max_minor_size: usize) -> bool {
    check_total_positivity(mat, max_minor_size, false).is_ok()
}

/// Check total non-negativity via Neville elimination (adjacent-row subtraction).
///
/// This is *much* faster than [`check_total_positivity`]: O(n²m) instead of
/// exponential in `min(n,m)`. It checks ALL minors, not just up to a given size.
///
/// A matrix is TNN iff Neville elimination completes with all multipliers ≥ 0
/// and no negative entries arise during the process.
///
/// Returns `Ok(())` if TNN, or `Err(msg)` describing the failure.
///
/// # Example
/// ```
/// use polynomial_tools::linalg::check_tnn_neville;
/// let pascal = vec![
///     vec![1, 0, 0],
///     vec![1, 1, 0],
///     vec![1, 2, 1],
///     vec![1, 3, 3],
///     vec![1, 4, 6],
/// ];
/// assert!(check_tnn_neville(&pascal).is_ok());
/// ```
pub fn check_tnn_neville(mat: &[Vec<i64>]) -> Result<(), String> {
    let nrows = mat.len();
    if nrows == 0 {
        return Ok(());
    }
    let ncols = mat[0].len();
    if ncols == 0 {
        return Ok(());
    }

    // Convert to rationals for exact arithmetic
    let mut a: Vec<Vec<Q>> = mat
        .iter()
        .map(|row| {
            row.iter()
                .map(|&v| Q::from_integer(BigInt::from(v)))
                .collect()
        })
        .collect();

    // Check all entries are non-negative
    for i in 0..nrows {
        for j in 0..ncols {
            if a[i][j] < Q::zero() {
                return Err(format!("Entry [{},{}] = {} < 0", i, j, mat[i][j]));
            }
        }
    }

    // Neville elimination: for each column k, eliminate bottom-up using adjacent rows
    let min_dim = nrows.min(ncols);
    for k in 0..min_dim {
        for i in (k + 1..nrows).rev() {
            if a[i][k].is_zero() {
                continue;
            }
            if a[i - 1][k].is_zero() {
                // Zero pivot with positive entry below: swap rows
                // Valid for TNN since both rows have zeros in columns < k
                a.swap(i - 1, i);
                continue;
            }
            // Both a[i-1][k] > 0 and a[i][k] > 0: compute multiplier
            let p = a[i][k].clone() / a[i - 1][k].clone();
            // Multiplier must be non-negative (guaranteed since both entries ≥ 0)
            debug_assert!(p >= Q::zero());

            // Subtract: row[i] -= p * row[i-1]
            for j in k..ncols {
                let sub = p.clone() * a[i - 1][j].clone();
                a[i][j] -= sub;
            }
            a[i][k] = Q::zero(); // exact zero in the eliminated position

            // Check all remaining entries in this row are non-negative
            for j in (k + 1)..ncols {
                if a[i][j] < Q::zero() {
                    return Err(format!(
                        "Neville elimination: entry [{},{}] became negative \
                         (multiplier at column {}, rows [{},{}])",
                        i,
                        j,
                        k,
                        i - 1,
                        i
                    ));
                }
            }
        }
    }

    Ok(())
}

/// Fast TNN check via Neville elimination. Returns bool.
/// Equivalent to `check_tnn_neville(mat).is_ok()`.
pub fn is_tnn(mat: &[Vec<i64>]) -> bool {
    check_tnn_neville(mat).is_ok()
}

/// Check TNN via Neville elimination on a BigInt matrix.
/// Use this when entries may overflow i64.
pub fn check_tnn_neville_bigint(mat: &[Vec<BigInt>]) -> Result<(), String> {
    let nrows = mat.len();
    if nrows == 0 {
        return Ok(());
    }
    let ncols = mat[0].len();
    if ncols == 0 {
        return Ok(());
    }

    let mut a: Vec<Vec<Q>> = mat
        .iter()
        .map(|row| row.iter().map(|v| Q::from_integer(v.clone())).collect())
        .collect();

    for i in 0..nrows {
        for j in 0..ncols {
            if a[i][j] < Q::zero() {
                return Err(format!("Entry [{},{}] < 0", i, j));
            }
        }
    }

    let min_dim = nrows.min(ncols);
    for k in 0..min_dim {
        for i in (k + 1..nrows).rev() {
            if a[i][k].is_zero() {
                continue;
            }
            if a[i - 1][k].is_zero() {
                a.swap(i - 1, i);
                continue;
            }
            let p = a[i][k].clone() / a[i - 1][k].clone();
            for j in k..ncols {
                let sub = p.clone() * a[i - 1][j].clone();
                a[i][j] -= sub;
            }
            a[i][k] = Q::zero();
            for j in (k + 1)..ncols {
                if a[i][j] < Q::zero() {
                    return Err(format!(
                        "Neville: entry [{},{}] negative (col {}, rows [{},{}])",
                        i,
                        j,
                        k,
                        i - 1,
                        i
                    ));
                }
            }
        }
    }
    Ok(())
}

fn extract_submatrix_i64(mat: &[Vec<i64>], rows: &[usize], cols: &[usize]) -> Vec<Vec<BigInt>> {
    rows.iter()
        .map(|&r| cols.iter().map(|&c| BigInt::from(mat[r][c])).collect())
        .collect()
}

fn combinations_usize(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut combo = vec![0usize; k];
    fn gen(
        pos: usize,
        start: usize,
        n: usize,
        k: usize,
        combo: &mut Vec<usize>,
        result: &mut Vec<Vec<usize>>,
    ) {
        if pos == k {
            result.push(combo.clone());
            return;
        }
        for i in start..n {
            combo[pos] = i;
            gen(pos + 1, i + 1, n, k, combo, result);
        }
    }
    gen(0, 0, n, k, &mut combo, &mut result);
    result
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn bigint_to_q(mat: &[Vec<BigInt>]) -> Vec<Vec<Q>> {
    mat.iter()
        .map(|row| row.iter().map(|v| Q::from_integer(v.clone())).collect())
        .collect()
}

fn trim_bigint_coeffs(coeffs: &mut Vec<BigInt>) {
    while coeffs.last().is_some_and(BigInt::is_zero) {
        coeffs.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bi(v: i64) -> BigInt {
        BigInt::from(v)
    }

    fn bi_mat(rows: &[&[i64]]) -> Vec<Vec<BigInt>> {
        rows.iter()
            .map(|row| row.iter().map(|&v| bi(v)).collect())
            .collect()
    }

    fn poly(coeffs: &[i64]) -> Polynomial<BigInt> {
        Polynomial::from_i64_coeffs(coeffs)
    }

    // -----------------------------------------------------------------------
    // Determinant
    // -----------------------------------------------------------------------

    #[test]
    fn test_determinant_empty() {
        assert_eq!(determinant(&[]), BigInt::one());
    }

    #[test]
    fn test_determinant_1x1() {
        assert_eq!(determinant(&bi_mat(&[&[7]])), bi(7));
    }

    #[test]
    fn test_determinant_2x2() {
        // [[1, 2], [3, 4]] -> det = 1*4 - 2*3 = -2
        let m = bi_mat(&[&[1, 2], &[3, 4]]);
        assert_eq!(determinant(&m), bi(-2));
    }

    #[test]
    fn test_determinant_3x3() {
        // [[1, 2, 3], [4, 5, 6], [7, 8, 10]] -> det = 1*(50-48) - 2*(40-42) + 3*(32-35) = 2+4-9 = -3
        let m = bi_mat(&[&[1, 2, 3], &[4, 5, 6], &[7, 8, 10]]);
        assert_eq!(determinant(&m), bi(-3));
    }

    #[test]
    fn test_determinant_row_pivot_after_first_step() {
        let m = bi_mat(&[&[1, 1, 0], &[1, 1, 1], &[0, 1, 1]]);
        assert_eq!(determinant(&m), bi(-1));
    }

    #[test]
    fn test_determinant_singular() {
        // [[1, 2], [2, 4]] -> det = 0
        let m = bi_mat(&[&[1, 2], &[2, 4]]);
        assert_eq!(determinant(&m), BigInt::zero());
    }

    #[test]
    #[should_panic(expected = "determinant requires a square matrix")]
    fn test_determinant_rejects_nonsquare_matrix() {
        let m = bi_mat(&[&[1, 2, 3], &[4, 5, 6]]);

        let _ = determinant(&m);
    }

    #[test]
    fn test_determinant_identity() {
        let m = bi_mat(&[&[1, 0, 0], &[0, 1, 0], &[0, 0, 1]]);
        assert_eq!(determinant(&m), bi(1));
    }

    #[test]
    fn test_bareiss_leading_principal_minors_bigint() {
        let m = bi_mat(&[&[2, 1], &[1, 2]]);
        assert_eq!(
            bareiss_leading_principal_minors_bigint(&m),
            Some(vec![bi(2), bi(3)])
        );
        assert_eq!(bareiss_determinant_bigint(&m), Some(bi(3)));
    }

    #[test]
    fn test_bareiss_bigint_no_pivot_zero_pivot() {
        let m = bi_mat(&[&[0, 1], &[1, 0]]);
        assert_eq!(determinant(&m), bi(-1));
        assert_eq!(bareiss_leading_principal_minors_bigint(&m), None);
    }

    #[test]
    fn test_bareiss_bigint_final_zero_pivot() {
        let m = bi_mat(&[&[1, 1], &[1, 1]]);
        assert_eq!(
            bareiss_leading_principal_minors_bigint(&m),
            Some(vec![bi(1), bi(0)])
        );
        assert_eq!(bareiss_determinant_bigint(&m), Some(bi(0)));
    }

    #[test]
    fn test_modular_leading_principal_minors_match_bareiss() {
        let m = bi_mat(&[&[4, 2, 1], &[2, 5, 3], &[1, 3, 6]]);
        let bareiss = bareiss_leading_principal_minors_bigint(&m).unwrap();
        let modular = modular_leading_principal_minors_bigint(&m).unwrap();
        assert_eq!(modular, bareiss);
    }

    #[test]
    fn test_modular_leading_principal_minors_handles_zero_bareiss_pivot() {
        let m = bi_mat(&[&[0, 1], &[1, 0]]);
        assert_eq!(bareiss_leading_principal_minors_bigint(&m), None);
        assert_eq!(
            modular_leading_principal_minors_bigint(&m),
            Some(vec![bi(0), bi(-1)])
        );
    }

    fn shifted_hilbert_like_matrix(size: usize) -> Vec<Vec<BigInt>> {
        let mut mat = vec![vec![BigInt::zero(); size]; size];
        for i in 0..size {
            for j in 0..size {
                mat[i][j] = BigInt::from((i + j + 1) as i64);
            }
            mat[i][i] += BigInt::from((size * size) as i64);
        }
        mat
    }

    fn large_positive_diagonal_matrix(size: usize) -> Vec<Vec<BigInt>> {
        let diagonal = BigInt::one() << 700usize;
        let mut mat = vec![vec![BigInt::zero(); size]; size];
        for (i, row) in mat.iter_mut().enumerate() {
            row[i] = diagonal.clone();
        }
        mat
    }

    #[test]
    fn test_bareiss_polynomial_bigint_2x2() {
        // [[1+z, z], [z, 1+z]] has leading minors 1+z and 1+2z.
        let m = vec![
            vec![poly(&[1, 1]), poly(&[0, 1])],
            vec![poly(&[0, 1]), poly(&[1, 1])],
        ];
        assert_eq!(
            bareiss_leading_principal_minors_polynomial_bigint(&m),
            Some(vec![poly(&[1, 1]), poly(&[1, 2])])
        );
        assert_eq!(
            bareiss_determinant_polynomial_bigint(&m),
            Some(poly(&[1, 2]))
        );
    }

    #[test]
    fn test_bareiss_polynomial_bigint_final_zero_pivot() {
        let m = vec![
            vec![poly(&[1, 1]), poly(&[1, 1])],
            vec![poly(&[1, 1]), poly(&[1, 1])],
        ];
        assert_eq!(
            bareiss_leading_principal_minors_polynomial_bigint(&m),
            Some(vec![poly(&[1, 1]), Polynomial::zero()])
        );
        assert_eq!(
            bareiss_determinant_polynomial_bigint(&m),
            Some(Polynomial::zero())
        );
    }

    #[test]
    fn test_bareiss_polynomial_bigint_nonconstant_exact_division() {
        // The final stage divides by the nonconstant first pivot 1+z.
        let m = vec![
            vec![poly(&[1, 1]), poly(&[0, 1]), Polynomial::zero()],
            vec![poly(&[0, 1]), poly(&[1, 1]), Polynomial::zero()],
            vec![Polynomial::zero(), Polynomial::zero(), poly(&[1, 1])],
        ];
        assert_eq!(
            bareiss_leading_principal_minors_polynomial_bigint(&m),
            Some(vec![poly(&[1, 1]), poly(&[1, 2]), poly(&[1, 3, 2])])
        );
        assert_eq!(
            bareiss_determinant_polynomial_bigint(&m),
            Some(poly(&[1, 3, 2]))
        );
    }

    // -----------------------------------------------------------------------
    // Positive definiteness
    // -----------------------------------------------------------------------

    #[test]
    fn test_pd_identity() {
        let m = bi_mat(&[&[1, 0, 0], &[0, 1, 0], &[0, 0, 1]]);
        assert!(is_positive_definite(&m));
    }

    #[test]
    fn test_pd_2x2() {
        // [[2, 1], [1, 2]] -> pivots: 2, 2 - 1/2 = 3/2 > 0
        let m = bi_mat(&[&[2, 1], &[1, 2]]);
        assert!(is_positive_definite(&m));
        assert!(is_positive_definite_bareiss(&m));
        assert_eq!(is_positive_definite_modular(&m), Some(true));
    }

    #[test]
    fn test_pd_rejects_nonsymmetric_matrix() {
        let m = bi_mat(&[&[1, 2], &[0, 1]]);

        assert!(!is_positive_definite(&m));
        assert!(!is_positive_definite_bareiss(&m));
        assert_eq!(is_positive_definite_modular(&m), Some(false));
    }

    #[test]
    fn test_pd_default_agrees_with_paths_across_threshold() {
        for size in [
            MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD - 1,
            MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD,
        ] {
            let m = shifted_hilbert_like_matrix(size);
            let bareiss = is_positive_definite_bareiss(&m);
            let modular = is_positive_definite_modular(&m).unwrap();
            assert_eq!(bareiss, modular);
            assert_eq!(is_positive_definite(&m), modular);
        }
    }

    #[test]
    fn test_pd_default_falls_back_when_modular_bounds_do_not_certify() {
        let m = large_positive_diagonal_matrix(MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD);

        assert!(modular_leading_principal_minors_bigint(&m).is_none());
        assert_eq!(is_positive_definite_modular(&m), None);
        assert!(is_positive_definite_bareiss(&m));
        assert!(is_positive_definite(&m));
    }

    #[test]
    fn test_pd_negative_definite() {
        // [[-1, 0], [0, -1]] -> first pivot = -1 < 0
        let m = bi_mat(&[&[-1, 0], &[0, -1]]);
        assert!(!is_positive_definite(&m));
        assert_eq!(is_positive_definite_modular(&m), Some(false));
    }

    #[test]
    fn test_pd_zero_matrix() {
        let m = bi_mat(&[&[0, 0], &[0, 0]]);
        assert!(!is_positive_definite(&m)); // PD requires strict positivity
    }

    #[test]
    fn test_pd_indefinite() {
        // [[1, 0], [0, -1]]
        let m = bi_mat(&[&[1, 0], &[0, -1]]);
        assert!(!is_positive_definite(&m));
        assert_eq!(is_positive_definite_modular(&m), Some(false));
    }

    // -----------------------------------------------------------------------
    // Positive semi-definiteness
    // -----------------------------------------------------------------------

    #[test]
    fn test_psd_identity() {
        let m = bi_mat(&[&[1, 0], &[0, 1]]);
        assert!(is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_zero_matrix() {
        let m = bi_mat(&[&[0, 0], &[0, 0]]);
        assert!(is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_rank_deficient() {
        // [[1, 1], [1, 1]] -> eigenvalues 0, 2 -> PSD
        let m = bi_mat(&[&[1, 1], &[1, 1]]);
        assert!(is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_indefinite() {
        // [[1, 2], [2, 1]] -> eigenvalues 3, -1 -> not PSD
        let m = bi_mat(&[&[1, 2], &[2, 1]]);
        assert!(!is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_rejects_nonsquare_matrix() {
        let m = bi_mat(&[&[1, 0, 0], &[0, 1, 0]]);

        assert!(!is_positive_semidefinite(&m));
    }

    #[test]
    fn test_psd_rejects_nonsymmetric_matrix() {
        let m = bi_mat(&[&[1, 0], &[2, 1]]);

        assert!(!is_positive_semidefinite(&m));
    }

    // -----------------------------------------------------------------------
    // Solve linear system
    // -----------------------------------------------------------------------

    #[test]
    fn test_solve_identity() {
        let a = vec![vec![Q::one(), Q::zero()], vec![Q::zero(), Q::one()]];
        let b = vec![Q::from_integer(bi(3)), Q::from_integer(bi(5))];
        let x = solve_linear_system(&a, &b).unwrap();
        assert_eq!(x[0], Q::from_integer(bi(3)));
        assert_eq!(x[1], Q::from_integer(bi(5)));
    }

    #[test]
    fn test_solve_2x2() {
        // x + 2y = 5, 3x + 4y = 11 -> x = 1, y = 2
        let a = vec![
            vec![Q::from_integer(bi(1)), Q::from_integer(bi(2))],
            vec![Q::from_integer(bi(3)), Q::from_integer(bi(4))],
        ];
        let b = vec![Q::from_integer(bi(5)), Q::from_integer(bi(11))];
        let x = solve_linear_system(&a, &b).unwrap();
        assert_eq!(x[0], Q::from_integer(bi(1)));
        assert_eq!(x[1], Q::from_integer(bi(2)));
    }

    #[test]
    fn test_solve_inconsistent() {
        // x + y = 1, x + y = 2 -> no solution
        let a = vec![vec![Q::one(), Q::one()], vec![Q::one(), Q::one()]];
        let b = vec![Q::one(), Q::from_integer(bi(2))];
        assert!(solve_linear_system(&a, &b).is_none());
    }

    // -----------------------------------------------------------------------
    // Total non-negativity
    // -----------------------------------------------------------------------

    #[test]
    fn test_tnn_pascal() {
        let pascal = vec![
            vec![1, 0, 0],
            vec![1, 1, 0],
            vec![1, 2, 1],
            vec![1, 3, 3],
            vec![1, 4, 6],
        ];
        assert!(is_tnn(&pascal));
        assert!(check_tnn_neville(&pascal).is_ok());
        assert!(check_total_positivity(&pascal, 3, false).is_ok());
    }

    #[test]
    fn test_tnn_identity() {
        let id = vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]];
        assert!(is_tnn(&id));
    }

    #[test]
    fn test_tnn_negative_entry() {
        let m = vec![vec![1, -1], vec![0, 1]];
        assert!(!is_tnn(&m));
    }

    #[test]
    fn test_tnn_not_tnn_positive_entries() {
        // All entries ≥ 0 but 2x2 minor det = 1*1 - 2*2 = -3 < 0
        let m = vec![vec![1, 2], vec![2, 1]];
        assert!(!is_tnn(&m));
    }

    #[test]
    fn test_tnn_path_matrix() {
        // Path matrix for P3: entry (i,j) = 1 if i <= j, 0 otherwise
        // These are always TNN
        let m = vec![vec![1, 1, 1], vec![0, 1, 1], vec![0, 0, 1]];
        assert!(is_tnn(&m));
    }

    #[test]
    fn test_tnn_neville_vs_brute() {
        // Verify Neville and brute-force agree on a 3x3 matrix
        let m = vec![vec![1, 1, 0], vec![0, 1, 1], vec![0, 0, 1]];
        let neville = check_tnn_neville(&m).is_ok();
        let brute = check_total_positivity(&m, 3, false).is_ok();
        assert_eq!(neville, brute);
    }

    #[test]
    fn test_tnn_bigint() {
        let m = vec![vec![bi(1), bi(1)], vec![bi(0), bi(1)]];
        assert!(check_tnn_neville_bigint(&m).is_ok());
    }
}
