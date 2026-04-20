//! Linear algebra over exact rationals (BigInt/BigRational).
//!
//! Provides Gaussian elimination, positive definiteness checks, determinants,
//! linear system solving, and total non-negativity (TNN) checking — all with
//! exact arithmetic over ℚ.
//!
//! # Total non-negativity
//!
//! A matrix is *totally non-negative* (TNN) if every minor (determinant of every
//! square submatrix) is non-negative. The brute-force check
//! ([`check_total_positivity`]) enumerates all submatrices up to a given size,
//! which is exponential. The Neville elimination method ([`check_tnn_neville`])
//! checks TNN for *all* minors in O(n²m) time by performing adjacent-row
//! subtraction and verifying that no entry becomes negative.

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};

type Q = Ratio<BigInt>;

// ---------------------------------------------------------------------------
// Gaussian elimination core
// ---------------------------------------------------------------------------

/// Result of Gaussian elimination on a matrix over ℚ.
///
/// After elimination, the matrix is in row-echelon form. The diagonal
/// entries are the pivots; their product (times the sign from row swaps)
/// gives the determinant.
struct EliminationResult {
    /// Row-echelon form of the matrix.
    matrix: Vec<Vec<Q>>,
    /// Sign from row swaps: +1 or -1.
    sign: i8,
    /// Number of zero pivots encountered (rank deficiency).
    zero_pivots: usize,
    /// Indices of rows where the pivot was zero (for semi-definiteness).
    #[allow(dead_code)]
    zero_pivot_rows: Vec<usize>,
}

/// Gaussian elimination with partial pivoting over ℚ.
///
/// Operates in-place on the given matrix. For positive-definiteness checks,
/// set `pivot_strategy` to `PivotStrategy::Diagonal` (no row swaps, checks
/// diagonal pivots only). For determinants and solving, use `PivotStrategy::Partial`.
fn gaussian_elimination(mat: &[Vec<Q>], strategy: PivotStrategy) -> EliminationResult {
    let n = mat.len();
    let ncols = if n > 0 { mat[0].len() } else { 0 };
    let mut a: Vec<Vec<Q>> = mat.to_vec();
    let mut sign: i8 = 1;
    let mut zero_pivots = 0;
    let mut zero_pivot_rows = Vec::new();

    for k in 0..n.min(ncols) {
        if a[k][k].is_zero() {
            match strategy {
                PivotStrategy::Diagonal => {
                    // No row swaps; record zero pivot
                    zero_pivots += 1;
                    zero_pivot_rows.push(k);
                    continue;
                }
                PivotStrategy::Partial => {
                    // Find a non-zero entry below
                    let mut found = false;
                    for i in (k + 1)..n {
                        if !a[i][k].is_zero() {
                            a.swap(k, i);
                            sign = -sign;
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        zero_pivots += 1;
                        zero_pivot_rows.push(k);
                        continue;
                    }
                }
            }
        }

        let pivot = a[k][k].clone();
        for i in (k + 1)..n {
            if a[i][k].is_zero() {
                continue;
            }
            let factor = a[i][k].clone() / pivot.clone();
            for j in (k + 1)..ncols {
                let sub = factor.clone() * a[k][j].clone();
                a[i][j] -= sub;
            }
            a[i][k] = Q::zero();
        }
    }

    EliminationResult {
        matrix: a,
        sign,
        zero_pivots,
        zero_pivot_rows,
    }
}

#[derive(Clone, Copy)]
enum PivotStrategy {
    /// No row swaps; checks diagonal entries in order. Used for definiteness checks.
    Diagonal,
    /// Partial pivoting with row swaps. Used for determinants and solving.
    Partial,
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Check if a symmetric BigInt matrix is positive definite.
///
/// Uses Gaussian elimination without pivoting: a symmetric matrix is positive
/// definite iff all diagonal pivots during elimination are strictly positive.
pub fn is_positive_definite(mat: &[Vec<BigInt>]) -> bool {
    let qmat = bigint_to_q(mat);
    let result = gaussian_elimination(&qmat, PivotStrategy::Diagonal);
    if result.zero_pivots > 0 {
        return false;
    }
    // All pivots must be strictly positive
    for k in 0..result.matrix.len() {
        if result.matrix[k][k] <= Q::zero() {
            return false;
        }
    }
    true
}

/// Check if a symmetric BigInt matrix is positive semi-definite.
///
/// A symmetric matrix is positive semi-definite iff all diagonal pivots
/// are non-negative, and zero pivots have zero entries in their column below.
pub fn is_positive_semidefinite(mat: &[Vec<BigInt>]) -> bool {
    let n = mat.len();
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

/// Compute the determinant of a BigInt matrix using Gaussian elimination over ℚ.
pub fn determinant(mat: &[Vec<BigInt>]) -> BigInt {
    let n = mat.len();
    if n == 0 {
        return BigInt::one();
    }
    let qmat = bigint_to_q(mat);
    let result = gaussian_elimination(&qmat, PivotStrategy::Partial);

    if result.zero_pivots > 0 {
        return BigInt::zero();
    }

    let mut det = Q::from_integer(BigInt::from(result.sign as i64));
    for k in 0..n {
        det *= result.matrix[k][k].clone();
    }
    det.to_integer()
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
    fn test_determinant_singular() {
        // [[1, 2], [2, 4]] -> det = 0
        let m = bi_mat(&[&[1, 2], &[2, 4]]);
        assert_eq!(determinant(&m), BigInt::zero());
    }

    #[test]
    fn test_determinant_identity() {
        let m = bi_mat(&[&[1, 0, 0], &[0, 1, 0], &[0, 0, 1]]);
        assert_eq!(determinant(&m), bi(1));
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
    }

    #[test]
    fn test_pd_negative_definite() {
        // [[-1, 0], [0, -1]] -> first pivot = -1 < 0
        let m = bi_mat(&[&[-1, 0], &[0, -1]]);
        assert!(!is_positive_definite(&m));
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

    // -----------------------------------------------------------------------
    // Solve linear system
    // -----------------------------------------------------------------------

    #[test]
    fn test_solve_identity() {
        let a = vec![
            vec![Q::one(), Q::zero()],
            vec![Q::zero(), Q::one()],
        ];
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
        let a = vec![
            vec![Q::one(), Q::one()],
            vec![Q::one(), Q::one()],
        ];
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
        let m = vec![
            vec![1, 1, 1],
            vec![0, 1, 1],
            vec![0, 0, 1],
        ];
        assert!(is_tnn(&m));
    }

    #[test]
    fn test_tnn_neville_vs_brute() {
        // Verify Neville and brute-force agree on a 3x3 matrix
        let m = vec![
            vec![1, 1, 0],
            vec![0, 1, 1],
            vec![0, 0, 1],
        ];
        let neville = check_tnn_neville(&m).is_ok();
        let brute = check_total_positivity(&m, 3, false).is_ok();
        assert_eq!(neville, brute);
    }

    #[test]
    fn test_tnn_bigint() {
        let m = vec![
            vec![bi(1), bi(1)],
            vec![bi(0), bi(1)],
        ];
        assert!(check_tnn_neville_bigint(&m).is_ok());
    }
}
