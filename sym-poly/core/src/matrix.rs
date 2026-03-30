//! Integer matrix utilities for transition matrix computations.
//!
//! These are extracted from combinatoric-core's transition module and
//! provide the linear algebra foundation for basis conversions in Sym, QSym, etc.

/// Identity matrix of size n.
pub fn identity_matrix(n: usize) -> Vec<Vec<i64>> {
    let mut m = vec![vec![0i64; n]; n];
    for i in 0..n {
        m[i][i] = 1;
    }
    m
}

/// Matrix multiplication: C = A * B.
pub fn mat_mul(a: &[Vec<i64>], b: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let rows = a.len();
    if rows == 0 {
        return vec![];
    }
    let inner = a[0].len();
    let cols = if inner > 0 { b[0].len() } else { 0 };
    let mut result = vec![vec![0i64; cols]; rows];
    for i in 0..rows {
        for k in 0..inner {
            if a[i][k] == 0 {
                continue;
            }
            for j in 0..cols {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    result
}

/// Transpose a matrix.
pub fn transpose(m: &[Vec<i64>]) -> Vec<Vec<i64>> {
    if m.is_empty() {
        return vec![];
    }
    let rows = m.len();
    let cols = m[0].len();
    let mut t = vec![vec![0i64; rows]; cols];
    for i in 0..rows {
        for j in 0..cols {
            t[j][i] = m[i][j];
        }
    }
    t
}

/// Invert an integer matrix using exact Gaussian elimination over rationals.
/// Panics if not invertible.
pub fn invert_integer_matrix(m: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let n = m.len();
    if n == 0 {
        return vec![];
    }

    use num_rational::Ratio;
    type Q = Ratio<i64>;

    let mut aug: Vec<Vec<Q>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = Vec::with_capacity(2 * n);
        for j in 0..n {
            row.push(Q::from_integer(m[i][j]));
        }
        for j in 0..n {
            row.push(if i == j {
                Q::from_integer(1)
            } else {
                Q::from_integer(0)
            });
        }
        aug.push(row);
    }

    // Forward elimination
    for col in 0..n {
        let mut pivot = None;
        for row in col..n {
            if aug[row][col] != Q::from_integer(0) {
                pivot = Some(row);
                break;
            }
        }
        let pivot = pivot.expect("matrix is singular");
        aug.swap(col, pivot);

        let diag = aug[col][col].clone();
        for j in 0..2 * n {
            aug[col][j] = aug[col][j].clone() / diag.clone();
        }

        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = aug[row][col].clone();
            for j in 0..2 * n {
                let val = aug[col][j].clone() * factor.clone();
                aug[row][j] = aug[row][j].clone() - val;
            }
        }
    }

    // Extract inverse (should be integer for our use cases)
    let mut inv = vec![vec![0i64; n]; n];
    for i in 0..n {
        for j in 0..n {
            let val = &aug[i][n + j];
            assert!(
                *val.denom() == 1 || *val.denom() == -1,
                "inverse matrix entry ({},{}) is not integer: {}",
                i,
                j,
                val
            );
            inv[i][j] = if *val.denom() == -1 {
                -(*val.numer())
            } else {
                *val.numer()
            };
        }
    }

    inv
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity() {
        let id = identity_matrix(3);
        assert_eq!(id[0], vec![1, 0, 0]);
        assert_eq!(id[1], vec![0, 1, 0]);
        assert_eq!(id[2], vec![0, 0, 1]);
    }

    #[test]
    fn test_mat_mul_identity() {
        let id = identity_matrix(2);
        let a = vec![vec![1, 2], vec![3, 4]];
        assert_eq!(mat_mul(&a, &id), a);
        assert_eq!(mat_mul(&id, &a), a);
    }

    #[test]
    fn test_transpose() {
        let m = vec![vec![1, 2, 3], vec![4, 5, 6]];
        let t = transpose(&m);
        assert_eq!(t, vec![vec![1, 4], vec![2, 5], vec![3, 6]]);
    }

    #[test]
    fn test_invert() {
        // [[2, 1], [1, 1]] has inverse [[1, -1], [-1, 2]]
        let m = vec![vec![2, 1], vec![1, 1]];
        let inv = invert_integer_matrix(&m);
        assert_eq!(inv, vec![vec![1, -1], vec![-1, 2]]);

        // Verify: M * M^{-1} = I
        let product = mat_mul(&m, &inv);
        assert_eq!(product, identity_matrix(2));
    }
}
