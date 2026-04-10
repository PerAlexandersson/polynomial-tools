//! Lightweight polynomial arithmetic on `&[i64]` coefficient vectors.
//!
//! These are convenience functions for the common pattern in research experiments
//! where polynomials are represented as `Vec<i64>` in ascending degree order
//! (`coeffs[i]` = coefficient of t^i). For generic coefficient types, use
//! [`Polynomial<C>`](crate::polynomial::Polynomial) instead.

/// Trim trailing zeros from a coefficient vector.
pub fn trim(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.len() > 1 && v.last() == Some(&0) {
        v.pop();
    }
    if v.is_empty() {
        vec![0]
    } else {
        v
    }
}

/// Whether the polynomial is zero.
pub fn is_zero(p: &[i64]) -> bool {
    p.iter().all(|&c| c == 0)
}

/// Degree of the polynomial, or `None` if zero.
pub fn degree(p: &[i64]) -> Option<usize> {
    for i in (0..p.len()).rev() {
        if p[i] != 0 {
            return Some(i);
        }
    }
    None
}

/// Add two polynomials.
pub fn add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] += v;
    }
    trim(&r)
}

/// Subtract: a - b.
pub fn sub(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] -= v;
    }
    trim(&r)
}

/// Multiply two polynomials.
pub fn mul(a: &[i64], b: &[i64]) -> Vec<i64> {
    if is_zero(a) || is_zero(b) {
        return vec![0];
    }
    let mut r = vec![0i64; a.len() + b.len() - 1];
    for (i, &av) in a.iter().enumerate() {
        if av == 0 {
            continue;
        }
        for (j, &bv) in b.iter().enumerate() {
            r[i + j] += av * bv;
        }
    }
    trim(&r)
}

/// Multiply by t (shift coefficients up by one degree).
pub fn shift(p: &[i64]) -> Vec<i64> {
    if is_zero(p) {
        return vec![0];
    }
    let mut r = vec![0i64; p.len() + 1];
    for (i, &v) in p.iter().enumerate() {
        r[i + 1] = v;
    }
    r
}

/// Get the coefficient of t^k, or 0 if out of range.
pub fn coeff(p: &[i64], k: usize) -> i64 {
    if k < p.len() {
        p[k]
    } else {
        0
    }
}

/// Evaluate polynomial at an integer point.
pub fn evaluate(p: &[i64], x: i64) -> i64 {
    let mut result = 0i64;
    for &c in p.iter().rev() {
        result = result * x + c;
    }
    result
}

/// Multiply all coefficients by a scalar.
pub fn scale(p: &[i64], c: i64) -> Vec<i64> {
    if c == 0 {
        return vec![0];
    }
    trim(&p.iter().map(|&v| v * c).collect::<Vec<_>>())
}

/// Negate a polynomial.
pub fn neg(p: &[i64]) -> Vec<i64> {
    p.iter().map(|&v| -v).collect()
}

/// Permanent of a matrix whose entries are polynomials (stored as `Vec<i64>`).
///
/// The matrix is `mat[row][col]`, where each entry is a polynomial.
/// Uses expansion by minors (Laplace expansion along the first row).
pub fn permanent(mat: &[Vec<Vec<i64>>]) -> Vec<i64> {
    let n = mat.len();
    if n == 0 {
        return vec![1];
    }
    let m = mat[0].len();
    let mut result = vec![0i64];
    for j in 0..m {
        if is_zero(&mat[0][j]) {
            continue;
        }
        let sub: Vec<Vec<Vec<i64>>> = (1..n)
            .map(|i| {
                (0..m)
                    .filter(|&jj| jj != j)
                    .map(|jj| mat[i][jj].clone())
                    .collect()
            })
            .collect();
        let term = mul(&mat[0][j], &permanent(&sub));
        result = add(&result, &term);
    }
    trim(&result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim() {
        assert_eq!(trim(&[1, 2, 0, 0]), vec![1, 2]);
        assert_eq!(trim(&[0, 0, 0]), vec![0]);
        assert_eq!(trim(&[3]), vec![3]);
    }

    #[test]
    fn test_is_zero() {
        assert!(is_zero(&[0]));
        assert!(is_zero(&[0, 0, 0]));
        assert!(!is_zero(&[1]));
        assert!(!is_zero(&[0, 1]));
    }

    #[test]
    fn test_degree() {
        assert_eq!(degree(&[0]), None);
        assert_eq!(degree(&[3]), Some(0));
        assert_eq!(degree(&[1, 2, 3]), Some(2));
        assert_eq!(degree(&[1, 0, 0]), Some(0));
    }

    #[test]
    fn test_add() {
        assert_eq!(add(&[1, 2], &[3, 4, 5]), vec![4, 6, 5]);
        assert_eq!(add(&[1, -1], &[-1, 1]), vec![0]);
    }

    #[test]
    fn test_sub() {
        assert_eq!(sub(&[3, 4, 5], &[1, 2]), vec![2, 2, 5]);
        assert_eq!(sub(&[1, 2], &[1, 2]), vec![0]);
    }

    #[test]
    fn test_mul() {
        // (1 + t)(1 + t) = 1 + 2t + t^2
        assert_eq!(mul(&[1, 1], &[1, 1]), vec![1, 2, 1]);
        // (1 + t)(1 - t) = 1 - t^2
        assert_eq!(mul(&[1, 1], &[1, -1]), vec![1, 0, -1]);
        assert_eq!(mul(&[0], &[1, 2, 3]), vec![0]);
    }

    #[test]
    fn test_shift() {
        assert_eq!(shift(&[1, 2, 3]), vec![0, 1, 2, 3]);
        assert_eq!(shift(&[0]), vec![0]);
    }

    #[test]
    fn test_coeff() {
        assert_eq!(coeff(&[10, 20, 30], 0), 10);
        assert_eq!(coeff(&[10, 20, 30], 2), 30);
        assert_eq!(coeff(&[10, 20, 30], 5), 0);
    }

    #[test]
    fn test_evaluate() {
        // p(t) = 1 + 2t + 3t^2, p(2) = 1 + 4 + 12 = 17
        assert_eq!(evaluate(&[1, 2, 3], 2), 17);
        assert_eq!(evaluate(&[5], 100), 5);
        assert_eq!(evaluate(&[0, 0, 1], 3), 9);
    }

    #[test]
    fn test_scale() {
        assert_eq!(scale(&[1, 2, 3], 2), vec![2, 4, 6]);
        assert_eq!(scale(&[1, 2, 3], 0), vec![0]);
        assert_eq!(scale(&[1, 2, 3], -1), vec![-1, -2, -3]);
    }

    #[test]
    fn test_permanent_1x1() {
        let mat = vec![vec![vec![0, 1]]]; // entry = t
        assert_eq!(permanent(&mat), vec![0, 1]);
    }

    #[test]
    fn test_permanent_2x2_identity() {
        // [[1, 0], [0, 1]] -> perm = 1*1 + 0*0 = 1
        let mat = vec![
            vec![vec![1], vec![0]],
            vec![vec![0], vec![1]],
        ];
        assert_eq!(permanent(&mat), vec![1]);
    }

    #[test]
    fn test_permanent_2x2_all_ones() {
        // [[1, 1], [1, 1]] -> perm = 1*1 + 1*1 = 2
        let mat = vec![
            vec![vec![1], vec![1]],
            vec![vec![1], vec![1]],
        ];
        assert_eq!(permanent(&mat), vec![2]);
    }

    #[test]
    fn test_permanent_2x2_with_t() {
        // [[1, t], [t, 1]] -> perm = 1*1 + t*t = 1 + t^2
        let mat = vec![
            vec![vec![1], vec![0, 1]],
            vec![vec![0, 1], vec![1]],
        ];
        assert_eq!(permanent(&mat), vec![1, 0, 1]);
    }

    #[test]
    fn test_permanent_empty() {
        let mat: Vec<Vec<Vec<i64>>> = vec![];
        assert_eq!(permanent(&mat), vec![1]);
    }
}
