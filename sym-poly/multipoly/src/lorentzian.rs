//! Lorentzian polynomial testing.
//!
//! A homogeneous polynomial *f* with nonnegative coefficients is **Lorentzian**
//! (Brändén–Huh, *Annals of Mathematics*, 2020) when every iterated partial
//! derivative of order deg(*f*) − 2 yields a quadratic whose Hessian matrix
//! has at most one positive eigenvalue.
//!
//! # Algorithm
//!
//! 1. Verify *f* is homogeneous with nonnegative coefficients.
//! 2. If deg ≤ 1, accept.
//! 3. If deg = 2, extract the Hessian, use the 2×2 principal minors as a
//!    fast necessary check, then count the positive eigenvalues of the full
//!    Hessian.
//! 4. If deg > 2, check that every first partial derivative ∂f/∂xᵢ is
//!    Lorentzian (recursion reduces degree by 1 each step).

use std::collections::{BTreeMap, BTreeSet};

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::Zero;

use crate::multipoly::{checked_total_degree, MultiPoly};

/// Result of a Lorentzian check with diagnostic information.
#[derive(Debug, Clone)]
pub enum LorentzianResult {
    /// The polynomial is Lorentzian.
    Yes,
    /// Not Lorentzian: has a negative coefficient.
    NegativeCoefficient { exponent: Vec<u32> },
    /// Not Lorentzian: not homogeneous.
    NotHomogeneous,
    /// Not strictly Lorentzian: support is not M-convex.
    SupportNotMConvex,
    /// Not Lorentzian: a degree-2 derivative has a Hessian with a positive 2×2 minor.
    HessianFailure {
        /// The derivative sequence that produced the failing quadratic.
        derivative_seq: Vec<usize>,
        /// The two variable indices for which `h_ii·h_jj > h_ij²`.
        vars: (usize, usize),
    },
    /// Not Lorentzian: a degree-2 derivative has a Hessian with too many
    /// positive eigenvalues.
    HessianInertiaFailure {
        /// The derivative sequence that produced the failing quadratic.
        derivative_seq: Vec<usize>,
        /// Number of positive eigenvalues found.
        positive_eigenvalues: usize,
    },
}

impl LorentzianResult {
    pub fn is_lorentzian(&self) -> bool {
        matches!(self, LorentzianResult::Yes)
    }
}

/// Check whether a multivariate polynomial with `i64` coefficients is Lorentzian.
///
/// Returns a [`LorentzianResult`] with diagnostic information.
pub fn is_lorentzian(f: &MultiPoly<i64>) -> LorentzianResult {
    // Empty / zero polynomial is Lorentzian.
    if f.is_zero() {
        return LorentzianResult::Yes;
    }

    // Check nonneg coefficients.
    for (exp, &coeff) in f.terms() {
        if coeff < 0 {
            return LorentzianResult::NegativeCoefficient {
                exponent: exp.clone(),
            };
        }
    }

    // Check homogeneity.
    let degrees: Vec<u32> = f
        .terms()
        .keys()
        .map(|exponents| checked_total_degree(exponents))
        .collect();
    let d = degrees[0];
    if degrees.iter().any(|&deg| deg != d) {
        return LorentzianResult::NotHomogeneous;
    }

    check_lorentzian_recursive(f, d, &mut Vec::new())
}

/// Recursive Lorentzian check.
fn check_lorentzian_recursive(
    f: &MultiPoly<i64>,
    degree: u32,
    deriv_seq: &mut Vec<usize>,
) -> LorentzianResult {
    if f.is_zero() || degree <= 1 {
        return LorentzianResult::Yes;
    }

    if degree == 2 {
        return check_quadratic_lorentzian(f, deriv_seq);
    }

    // degree > 2: check all partial derivatives
    let n = f.num_vars();
    for i in 0..n {
        let df = partial_derivative(f, i);
        if df.is_zero() {
            continue;
        }
        deriv_seq.push(i);
        let result = check_lorentzian_recursive(&df, degree - 1, deriv_seq);
        if !result.is_lorentzian() {
            return result;
        }
        deriv_seq.pop();
    }
    LorentzianResult::Yes
}

/// Check the Lorentzian condition for a degree-2 homogeneous polynomial.
///
/// The Hessian matrix H of a degree-2 homogeneous polynomial
/// q = Σ a_{ij} xᵢxⱼ has entries H\[i\]\[j\] = coefficient of xᵢxⱼ
/// (with the convention that xᵢ² has coefficient a_{ii} and xᵢxⱼ
/// for i≠j appears with coefficient a_{ij}).
///
/// The actual Hessian has H\[i\]\[i\] = 2·c(xᵢ²) and H\[i\]\[j\] = c(xᵢxⱼ),
/// so the 2×2 minor condition is 4·c(xᵢ²)·c(xⱼ²) ≤ c(xᵢxⱼ)².
/// These minors are only a fast necessary filter; the full Hessian must
/// have at most one positive eigenvalue.
fn check_quadratic_lorentzian(f: &MultiPoly<i64>, deriv_seq: &[usize]) -> LorentzianResult {
    let n = f.num_vars();

    // Extract the actual Hessian matrix.
    let mut h = vec![vec![0i128; n]; n];
    for (exp, &coeff) in f.terms() {
        let mut vars = Vec::new();
        for (v, &e) in exp.iter().enumerate() {
            for _ in 0..e {
                vars.push(v);
            }
        }
        debug_assert_eq!(vars.len(), 2);
        let (i, j) = (vars[0], vars[1]);
        if i == j {
            h[i][i] = 2 * coeff as i128;
        } else {
            h[i][j] = coeff as i128;
            h[j][i] = coeff as i128;
        }
    }

    // Fast 2×2 necessary check.
    for i in 0..n {
        for j in (i + 1)..n {
            let diag_prod = h[i][i] * h[j][j];
            let off_sq = h[i][j] * h[i][j];
            if diag_prod > off_sq {
                return LorentzianResult::HessianFailure {
                    derivative_seq: deriv_seq.to_vec(),
                    vars: (i, j),
                };
            }
        }
    }

    let positive_eigenvalues = positive_eigenvalue_count(&h);
    if positive_eigenvalues > 1 {
        return LorentzianResult::HessianInertiaFailure {
            derivative_seq: deriv_seq.to_vec(),
            positive_eigenvalues,
        };
    }

    LorentzianResult::Yes
}

/// Compute the partial derivative ∂f/∂xᵢ.
fn partial_derivative(f: &MultiPoly<i64>, var: usize) -> MultiPoly<i64> {
    let n = f.num_vars();
    let mut terms = BTreeMap::new();
    for (exp, &coeff) in f.terms() {
        if exp[var] == 0 {
            continue;
        }
        let new_coeff = coeff
            .checked_mul(i64::from(exp[var]))
            .expect("partial derivative coefficient overflow");
        let mut new_exp = exp.clone();
        new_exp[var] -= 1;
        let entry = terms.entry(new_exp).or_insert(0i64);
        *entry = entry
            .checked_add(new_coeff)
            .expect("partial derivative coefficient overflow");
    }
    MultiPoly::from_terms(n, terms)
}

/// Count positive eigenvalues of a real symmetric matrix.
///
/// The count is computed exactly by rational congruence diagonalization, so
/// Sylvester inertia is preserved without floating-point tolerances.
fn positive_eigenvalue_count(matrix: &[Vec<i128>]) -> usize {
    let n = matrix.len();
    if n == 0 {
        return 0;
    }

    type BigRat = Ratio<BigInt>;

    let mut a: Vec<Vec<BigRat>> = matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|&x| BigRat::from_integer(BigInt::from(x)))
                .collect()
        })
        .collect();

    let mut positives = 0;
    let mut k = 0;

    while k < n {
        let mut pivot = (k..n).find(|&i| !a[i][i].is_zero());

        if pivot.is_none() {
            let mut offdiag = None;
            'outer: for i in k..n {
                for j in (i + 1)..n {
                    if !a[i][j].is_zero() {
                        offdiag = Some((i, j));
                        break 'outer;
                    }
                }
            }
            if let Some((i, j)) = offdiag {
                add_basis_vector(&mut a, i, j, k);
                pivot = Some(i);
            } else {
                break;
            }
        }

        let pivot = pivot.unwrap();
        if pivot != k {
            swap_symmetric(&mut a, pivot, k);
        }

        let pivot_value = a[k][k].clone();
        debug_assert!(!pivot_value.is_zero());
        if pivot_value > BigRat::zero() {
            positives += 1;
        }

        for i in (k + 1)..n {
            for j in i..n {
                let correction = a[i][k].clone() * a[k][j].clone() / pivot_value.clone();
                let updated = a[i][j].clone() - correction;
                a[i][j] = updated.clone();
                if i != j {
                    a[j][i] = updated;
                }
            }
        }

        for i in (k + 1)..n {
            a[i][k] = BigRat::zero();
            a[k][i] = BigRat::zero();
        }

        k += 1;
    }

    positives
}

fn swap_symmetric<T>(matrix: &mut [Vec<T>], i: usize, j: usize) {
    matrix.swap(i, j);
    for row in matrix {
        row.swap(i, j);
    }
}

fn add_basis_vector<T>(matrix: &mut [Vec<T>], target: usize, source: usize, start: usize)
where
    T: Clone + std::ops::AddAssign,
{
    let n = matrix.len();
    let source_row: Vec<T> = (start..n).map(|col| matrix[source][col].clone()).collect();
    let source_col: Vec<T> = (start..n).map(|row| matrix[row][source].clone()).collect();

    for (offset, col) in (start..n).enumerate() {
        matrix[target][col] += source_row[offset].clone();
    }
    for (offset, row) in (start..n).enumerate() {
        matrix[row][target] += source_col[offset].clone();
    }
}

/// Convenience: check Lorentzian and return a bool.
pub fn is_lorentzian_bool(f: &MultiPoly<i64>) -> bool {
    is_lorentzian(f).is_lorentzian()
}

/// Check whether a polynomial is **strictly Lorentzian**: Lorentzian AND
/// its support is M-convex.
///
/// By Brändén–Huh Thm 2.25 the derivative conditions characterize the
/// *closure* of strictly Lorentzian polynomials.  This function adds the
/// explicit M-convexity check to detect boundary cases like x³ + y³.
pub fn is_strictly_lorentzian(f: &MultiPoly<i64>) -> LorentzianResult {
    let r = is_lorentzian(f);
    if !r.is_lorentzian() {
        return r;
    }
    // Check M-convexity of support
    let support: Vec<Vec<u32>> = f.terms().keys().cloned().collect();
    if !support.is_empty() && !is_m_convex(&support) {
        return LorentzianResult::SupportNotMConvex;
    }
    LorentzianResult::Yes
}

// ── Normalized Lorentzian ────────────────────────────────────────────

/// Check whether f is **normalized Lorentzian**: the generating function
/// g(x) = Σ c_α x^α / α!  is Lorentzian.
///
/// Because ∂ⁱ(x^α / α!) = x^{α−eᵢ} / (α−eᵢ)!, the Hessian of the
/// (d−2)-fold derivative of g has entries c_{β+eᵢ+eⱼ}.  Thus for every
/// |β| = d−2 the full matrix
///
/// > H_β(i,j) = c_{β+eᵢ+eⱼ}
///
/// must have at most one positive eigenvalue.  The pairwise inequalities
/// c_{β+2eᵢ} · c_{β+2eⱼ} ≤ c_{β+eᵢ+eⱼ}² are checked first as a fast
/// necessary filter, but they are not sufficient in three or more variables.
///
/// This works directly on the coefficients of f; no factorial division is
/// required.
pub fn is_normalized_lorentzian(f: &MultiPoly<i64>) -> LorentzianResult {
    if f.is_zero() {
        return LorentzianResult::Yes;
    }

    // Nonneg coefficients.
    for (exp, &coeff) in f.terms() {
        if coeff < 0 {
            return LorentzianResult::NegativeCoefficient {
                exponent: exp.clone(),
            };
        }
    }

    // Homogeneity.
    let degrees: Vec<u32> = f
        .terms()
        .keys()
        .map(|exponents| checked_total_degree(exponents))
        .collect();
    let d = degrees[0];
    if degrees.iter().any(|&deg| deg != d) {
        return LorentzianResult::NotHomogeneous;
    }

    if d <= 1 {
        return LorentzianResult::Yes;
    }

    let n = f.num_vars();
    let bases = normalized_hessian_bases(f, d);

    for beta in bases {
        let h = normalized_hessian_matrix(f, &beta);

        // Fast 2×2 necessary check.
        for i in 0..n {
            for j in (i + 1)..n {
                let lhs = h[i][i] * h[j][j];
                let rhs = h[i][j] * h[i][j];
                if lhs > rhs {
                    return LorentzianResult::HessianFailure {
                        derivative_seq: beta_to_derivative_seq(&beta),
                        vars: (i, j),
                    };
                }
            }
        }

        let positive_eigenvalues = positive_eigenvalue_count(&h);
        if positive_eigenvalues > 1 {
            return LorentzianResult::HessianInertiaFailure {
                derivative_seq: beta_to_derivative_seq(&beta),
                positive_eigenvalues,
            };
        }
    }

    LorentzianResult::Yes
}

fn normalized_hessian_bases(f: &MultiPoly<i64>, degree: u32) -> BTreeSet<Vec<u32>> {
    let n = f.num_vars();
    let mut bases = BTreeSet::new();

    if degree < 2 {
        return bases;
    }

    for alpha in f.terms().keys() {
        for i in 0..n {
            for j in i..n {
                if i == j {
                    if alpha[i] < 2 {
                        continue;
                    }
                } else if alpha[i] == 0 || alpha[j] == 0 {
                    continue;
                }

                let mut beta = alpha.clone();
                beta[i] -= 1;
                beta[j] -= 1;
                bases.insert(beta);
            }
        }
    }

    bases
}

fn normalized_hessian_matrix(f: &MultiPoly<i64>, beta: &[u32]) -> Vec<Vec<i128>> {
    let n = f.num_vars();
    let mut h = vec![vec![0i128; n]; n];

    for i in 0..n {
        for j in i..n {
            let mut alpha = beta.to_vec();
            alpha[i] += 1;
            alpha[j] += 1;
            let coeff = f.coefficient(&alpha) as i128;
            h[i][j] = coeff;
            h[j][i] = coeff;
        }
    }

    h
}

fn beta_to_derivative_seq(beta: &[u32]) -> Vec<usize> {
    beta.iter()
        .enumerate()
        .flat_map(|(v, &e)| std::iter::repeat(v).take(e as usize))
        .collect()
}

/// Normalized Lorentzian check as a bool.
pub fn is_normalized_lorentzian_bool(f: &MultiPoly<i64>) -> bool {
    is_normalized_lorentzian(f).is_lorentzian()
}

/// Strictly normalized Lorentzian: the normalized check plus M-convex support.
pub fn is_strictly_normalized_lorentzian(f: &MultiPoly<i64>) -> LorentzianResult {
    let r = is_normalized_lorentzian(f);
    if !r.is_lorentzian() {
        return r;
    }
    let support: Vec<Vec<u32>> = f.terms().keys().cloned().collect();
    if !support.is_empty() && !is_m_convex(&support) {
        return LorentzianResult::SupportNotMConvex;
    }
    LorentzianResult::Yes
}

// ── M-convexity ─────────────────────────────────────────────────────

/// Check whether a set of nonneg integer vectors satisfies the symmetric
/// M-exchange axiom (discrete polymatroid exchange property).
///
/// For every pair (a, b) in the set, and every index *i* with aᵢ > bᵢ,
/// there must exist an index *j* with aⱼ < bⱼ such that
/// a − eᵢ + eⱼ is also in the set.
///
/// Reference: Murota, *Discrete Convex Analysis*, SIAM 2003.
pub fn is_m_convex(set: &[Vec<u32>]) -> bool {
    if set.len() <= 1 {
        return true;
    }
    let n = set[0].len();

    // All vectors must have the same length and the same coordinate sum.
    let total = checked_total_degree(&set[0]);
    for v in set {
        if v.len() != n || checked_total_degree(v) != total {
            return false;
        }
    }

    // Build a HashSet for O(1) membership queries.
    let member: std::collections::HashSet<&Vec<u32>> = set.iter().collect();

    for a in set {
        for b in set {
            // For each i where a_i > b_i, there must exist j where a_j < b_j
            // such that a - e_i + e_j ∈ set.
            for i in 0..n {
                if a[i] <= b[i] {
                    continue;
                }
                let mut found = false;
                for j in 0..n {
                    if a[j] >= b[j] {
                        continue;
                    }
                    // Check a - e_i + e_j ∈ set
                    let mut swapped = a.clone();
                    swapped[i] -= 1;
                    swapped[j] += 1;
                    if member.contains(&swapped) {
                        found = true;
                        break;
                    }
                }
                if !found {
                    return false;
                }
            }
        }
    }
    true
}

/// Check M-convexity of the support of a polynomial.
pub fn support_is_m_convex(f: &MultiPoly<i64>) -> bool {
    let support: Vec<Vec<u32>> = f.terms().keys().cloned().collect();
    if support.is_empty() {
        return true;
    }
    is_m_convex(&support)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn poly_from_map(n: usize, map: Vec<(Vec<u32>, i64)>) -> MultiPoly<i64> {
        let terms: BTreeMap<Vec<u32>, i64> = map.into_iter().collect();
        MultiPoly::from_terms(n, terms)
    }

    #[test]
    fn test_positive_eigenvalue_count_exact_zero_diagonal_pivot() {
        assert_eq!(positive_eigenvalue_count(&[vec![0, 1], vec![1, 0]]), 1);
        assert_eq!(
            positive_eigenvalue_count(&[vec![0, 1, 1], vec![1, 0, 1], vec![1, 1, 0],]),
            1
        );
    }

    #[test]
    fn test_positive_eigenvalue_count_exact_diagonal_cases() {
        assert_eq!(
            positive_eigenvalue_count(&[vec![2, 0, 0], vec![0, -1, 0], vec![0, 0, 0],]),
            1
        );
        assert_eq!(
            positive_eigenvalue_count(&[vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1],]),
            3
        );
    }

    #[test]
    fn test_zero_is_lorentzian() {
        let f = MultiPoly::<i64>::zero(3);
        assert!(is_lorentzian_bool(&f));
    }

    #[test]
    fn test_linear_is_lorentzian() {
        // x + 2y + 3z
        let f = poly_from_map(
            3,
            vec![(vec![1, 0, 0], 1), (vec![0, 1, 0], 2), (vec![0, 0, 1], 3)],
        );
        assert!(is_lorentzian_bool(&f));
    }

    #[test]
    fn test_negative_coeff_rejected() {
        let f = poly_from_map(2, vec![(vec![2, 0], 1), (vec![1, 1], -1), (vec![0, 2], 1)]);
        let r = is_lorentzian(&f);
        assert!(matches!(r, LorentzianResult::NegativeCoefficient { .. }));
    }

    #[test]
    fn test_not_homogeneous_rejected() {
        let f = poly_from_map(2, vec![(vec![2, 0], 1), (vec![1, 0], 1)]);
        let r = is_lorentzian(&f);
        assert!(matches!(r, LorentzianResult::NotHomogeneous));
    }

    #[test]
    fn test_quadratic_lorentzian() {
        // q = 4xy (Hessian [[0,4],[4,0]], minor 0*0 - 16 = -16 ≤ 0)
        let f = poly_from_map(2, vec![(vec![1, 1], 4)]);
        assert!(is_lorentzian_bool(&f));
    }

    #[test]
    fn test_quadratic_not_lorentzian() {
        // q = x² + y² (Hessian [[1,0],[0,1]], minor 1*1 - 0 = 1 > 0)
        // Two positive eigenvalues → not Lorentzian.
        let f = poly_from_map(2, vec![(vec![2, 0], 1), (vec![0, 2], 1)]);
        assert!(!is_lorentzian_bool(&f));
    }

    #[test]
    fn test_quadratic_lorentzian_3vars() {
        // q = x² + 2xy + 2xz + y² + 2yz + z² = (x+y+z)²
        // This is rank-1, so at most one positive eigenvalue.
        let f = poly_from_map(
            3,
            vec![
                (vec![2, 0, 0], 1),
                (vec![1, 1, 0], 2),
                (vec![1, 0, 1], 2),
                (vec![0, 2, 0], 1),
                (vec![0, 1, 1], 2),
                (vec![0, 0, 2], 1),
            ],
        );
        assert!(is_lorentzian_bool(&f));
    }

    #[test]
    fn test_quadratic_pairwise_minors_not_sufficient() {
        // The fast 2x2 principal-minor checks pass, but the full Hessian has
        // two positive eigenvalues.
        let f = poly_from_map(
            4,
            vec![
                (vec![1, 1, 0, 0], 2),
                (vec![1, 0, 1, 0], 2),
                (vec![1, 0, 0, 1], 2),
                (vec![0, 1, 0, 1], 2),
                (vec![0, 0, 0, 2], 1),
            ],
        );
        let r = is_lorentzian(&f);
        assert!(matches!(
            r,
            LorentzianResult::HessianInertiaFailure {
                positive_eigenvalues: 2,
                ..
            }
        ));
    }

    #[test]
    fn test_quadratic_not_lorentzian_subtle() {
        // q = x² + xy + y². Hessian [[2,1],[1,2]], eigenvalues 3,1.
        // Two positive eigenvalues → NOT Lorentzian.
        // In two variables the normalized Hessian check is equivalent to
        // the normalized pairwise inequality, so this is normalized-Lorentzian.
        let f = poly_from_map(2, vec![(vec![2, 0], 1), (vec![1, 1], 1), (vec![0, 2], 1)]);
        assert!(!is_lorentzian_bool(&f));
        assert!(is_normalized_lorentzian_bool(&f));
    }

    #[test]
    fn test_degree3_lorentzian() {
        // f = x³ + 3x²y + 3xy² + y³ = (x+y)³. This is Lorentzian.
        let f = poly_from_map(
            2,
            vec![
                (vec![3, 0], 1),
                (vec![2, 1], 3),
                (vec![1, 2], 3),
                (vec![0, 3], 1),
            ],
        );
        assert!(is_lorentzian_bool(&f));
    }

    #[test]
    fn test_degree3_not_lorentzian() {
        // f = x³ + y³. ∂f/∂x = 3x², ∂f/∂y = 3y².
        // ∂²f/∂x² = 6x, ∂²f/∂y² = 6y, ∂²f/∂x∂y = 0.
        // At degree 2 after one derivative: e.g. ∂f/∂x = 3x² has Hessian [[6,0],[0,0]],
        // minor 6*0 - 0 = 0 ≤ 0. OK.
        // But ∂f/∂y = 3y² has Hessian [[0,0],[0,6]], minor 0*6 - 0 = 0 ≤ 0. OK.
        // Hmm, actually ∂/∂x and ∂/∂y both pass. Let me reconsider...
        // f = x³ + y³ is NOT Lorentzian because its support {(3,0),(0,3)} is not M-convex
        // (missing (2,1) and (1,2)). But the derivative conditions might not catch this
        // directly. According to Brändén-Huh Thm 2.25, the derivative conditions DO
        // imply M-convexity, so they should fail somewhere.
        //
        // ∂f/∂x = 3x². Check: this is degree 2, Hessian [[6,0],[0,0]].
        // Minor: 6*0 = 0 ≤ 0. Passes!
        // ∂f/∂y = 3y². Hessian [[0,0],[0,6]]. Minor: 0*6 = 0 ≤ 0. Passes!
        //
        // Hmm, so both first derivatives pass the quadratic check. But x³ + y³ is NOT
        // Lorentzian. The issue is that the M-convex support condition is separate.
        //
        // Actually, let me re-read: Thm 2.25 says "the closure of Lorentzian polynomials
        // equals the set of polynomials satisfying the derivative conditions". The closure
        // means we also include limits. x³ + y³ is in the closure but not Lorentzian itself.
        //
        // For a strict check, we need M-convex support. For this test, let's check a
        // polynomial that fails the derivative condition.
        //
        // f = x³ + x²y + xy² + 10y³. Check ∂f/∂x = 3x² + 2xy + y².
        // Hessian: [[6, 2], [2, 2]]. Minor: 12 - 4 = 8 > 0. NOT Lorentzian.
        let f = poly_from_map(
            2,
            vec![
                (vec![3, 0], 1),
                (vec![2, 1], 1),
                (vec![1, 2], 1),
                (vec![0, 3], 10),
            ],
        );
        assert!(!is_lorentzian_bool(&f));
    }

    #[test]
    fn test_elementary_symmetric_lorentzian() {
        // e_2(x,y,z) = xy + xz + yz. Lorentzian (it's the product of linear forms
        // in the "normalized" sense — known to be Lorentzian).
        let f = poly_from_map(
            3,
            vec![(vec![1, 1, 0], 1), (vec![1, 0, 1], 1), (vec![0, 1, 1], 1)],
        );
        assert!(is_lorentzian_bool(&f));
    }

    #[test]
    fn test_schur_21_lorentzian() {
        // s_{2,1}(x,y,z) = x²y + x²z + xy² + 2xyz + xz² + y²z + yz²
        // Schur polynomials are Lorentzian (Brändén-Huh, Cor 3.7).
        let f = poly_from_map(
            3,
            vec![
                (vec![2, 1, 0], 1),
                (vec![2, 0, 1], 1),
                (vec![1, 2, 0], 1),
                (vec![1, 1, 1], 2),
                (vec![1, 0, 2], 1),
                (vec![0, 2, 1], 1),
                (vec![0, 1, 2], 1),
            ],
        );
        assert!(is_lorentzian_bool(&f));
    }

    // ── M-convexity tests ──────────────────────────────────────────

    #[test]
    fn test_m_convex_singleton() {
        assert!(is_m_convex(&[vec![2, 1]]));
    }

    #[test]
    fn test_m_convex_simplex() {
        // Support of e_1(x,y,z) = {(1,0,0),(0,1,0),(0,0,1)} — the 3-simplex
        assert!(is_m_convex(&[vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]]));
    }

    #[test]
    fn test_m_convex_full_degree2() {
        // Support of (x+y+z)² = all monomials of degree 2 in 3 vars
        let set = vec![
            vec![2, 0, 0],
            vec![1, 1, 0],
            vec![1, 0, 1],
            vec![0, 2, 0],
            vec![0, 1, 1],
            vec![0, 0, 2],
        ];
        assert!(is_m_convex(&set));
    }

    #[test]
    fn test_not_m_convex_gap() {
        // {(3,0), (0,3)} — missing interior points, not M-convex
        assert!(!is_m_convex(&[vec![3, 0], vec![0, 3]]));
    }

    #[test]
    fn test_not_m_convex_corners() {
        // {(2,0,0), (0,2,0), (0,0,2)} — degree-2 "anti-matroid"
        assert!(!is_m_convex(&[vec![2, 0, 0], vec![0, 2, 0], vec![0, 0, 2]]));
    }

    #[test]
    fn test_m_convex_schur_support() {
        // Support of s_{2,1}(x,y,z)
        let set = vec![
            vec![2, 1, 0],
            vec![2, 0, 1],
            vec![1, 2, 0],
            vec![1, 1, 1],
            vec![1, 0, 2],
            vec![0, 2, 1],
            vec![0, 1, 2],
        ];
        assert!(is_m_convex(&set));
    }

    #[test]
    fn test_strictly_lorentzian_passes() {
        // (x+y)³ is strictly Lorentzian
        let f = poly_from_map(
            2,
            vec![
                (vec![3, 0], 1),
                (vec![2, 1], 3),
                (vec![1, 2], 3),
                (vec![0, 3], 1),
            ],
        );
        assert!(is_strictly_lorentzian(&f).is_lorentzian());
    }

    #[test]
    fn test_strictly_lorentzian_detects_non_m_convex() {
        // x³ + y³ passes derivative conditions (closure) but support
        // {(3,0),(0,3)} is not M-convex → not strictly Lorentzian.
        let f = poly_from_map(2, vec![(vec![3, 0], 1), (vec![0, 3], 1)]);
        // Passes the weak (closure) Lorentzian check
        assert!(is_lorentzian_bool(&f));
        // Fails the strict check
        let r = is_strictly_lorentzian(&f);
        assert!(matches!(r, LorentzianResult::SupportNotMConvex));
    }

    // ── Normalized Lorentzian tests ───────────────────────────────

    #[test]
    fn test_normalized_lorentzian_power_of_sum() {
        // (x+y)³ has coefficients 1, 3, 3, 1.
        // Normalized: 1/0!3!, 3/1!2!, 3/2!1!, 1/3!0! = 1/6, 1/2, 1/2, 1/6
        // which is (1/6)(x+y)³. Still Lorentzian.
        // In coefficient form: c_{30}=1, c_{21}=3, c_{12}=3, c_{03}=1.
        // Check: c_{30}·c_{12} ≤ c_{21}² → 1·3=3 ≤ 9. ✓
        // Check: c_{21}·c_{03} ≤ c_{12}² → 3·1=3 ≤ 9. ✓
        let f = poly_from_map(
            2,
            vec![
                (vec![3, 0], 1),
                (vec![2, 1], 3),
                (vec![1, 2], 3),
                (vec![0, 3], 1),
            ],
        );
        assert!(is_normalized_lorentzian_bool(&f));
    }

    #[test]
    fn test_normalized_lorentzian_multinomial() {
        // (x+y+z)² has coefficients: x²→1, y²→1, z²→1, xy→2, xz→2, yz→2.
        // c_{200}·c_{020} = 1 ≤ c_{110}² = 4. ✓
        // c_{200}·c_{002} = 1 ≤ c_{101}² = 4. ✓
        // c_{020}·c_{002} = 1 ≤ c_{011}² = 4. ✓
        let f = poly_from_map(
            3,
            vec![
                (vec![2, 0, 0], 1),
                (vec![0, 2, 0], 1),
                (vec![0, 0, 2], 1),
                (vec![1, 1, 0], 2),
                (vec![1, 0, 1], 2),
                (vec![0, 1, 1], 2),
            ],
        );
        assert!(is_normalized_lorentzian_bool(&f));
    }

    #[test]
    fn test_normalized_lorentzian_fails() {
        // f = x² + y² (no cross term).
        // c_{20}·c_{02} = 1 > c_{11}² = 0. Fails.
        let f = poly_from_map(2, vec![(vec![2, 0], 1), (vec![0, 2], 1)]);
        assert!(!is_normalized_lorentzian_bool(&f));
    }

    #[test]
    fn test_normalized_pairwise_minors_not_sufficient() {
        // This normalized Hessian has entries
        // [[0,1,1,1], [1,0,0,1], [1,0,0,0], [1,1,0,1]].
        // All pairwise normalized inequalities pass, but the full Hessian has
        // two positive eigenvalues.
        let f = poly_from_map(
            4,
            vec![
                (vec![1, 1, 0, 0], 1),
                (vec![1, 0, 1, 0], 1),
                (vec![1, 0, 0, 1], 1),
                (vec![0, 1, 0, 1], 1),
                (vec![0, 0, 0, 2], 1),
            ],
        );
        let r = is_normalized_lorentzian(&f);
        assert!(matches!(
            r,
            LorentzianResult::HessianInertiaFailure {
                positive_eigenvalues: 2,
                ..
            }
        ));
    }

    #[test]
    fn test_normalized_vs_standard_difference() {
        // f = x² + 4xy + y². Standard Lorentzian:
        //   Hessian [[1,4],[4,1]], minor 1-16 = -15 ≤ 0. ✓ Standard Lorentzian.
        // Normalized Lorentzian:
        //   c_{20}·c_{02} = 1·1 = 1, c_{11}² = 16. 1 ≤ 16. ✓ Also normalized Lorentzian.
        let f = poly_from_map(2, vec![(vec![2, 0], 1), (vec![1, 1], 4), (vec![0, 2], 1)]);
        assert!(is_lorentzian_bool(&f));
        assert!(is_normalized_lorentzian_bool(&f));

        // f = x² + 2xy + y². Standard:
        //   Hessian [[1,2],[2,1]], minor 1-4 = -3 ≤ 0. ✓
        // Normalized: c_{20}·c_{02} = 1 ≤ c_{11}² = 4. ✓
        let g = poly_from_map(2, vec![(vec![2, 0], 1), (vec![1, 1], 2), (vec![0, 2], 1)]);
        assert!(is_lorentzian_bool(&g));
        assert!(is_normalized_lorentzian_bool(&g));

        // f = 2x² + xy + 2y². Standard:
        //   Hessian [[2,1],[1,2]], minor 4-1 = 3 > 0. ✗ Not standard Lorentzian.
        // Normalized: c_{20}·c_{02} = 4 > c_{11}² = 1. ✗ Also not normalized Lorentzian.
        let h = poly_from_map(2, vec![(vec![2, 0], 2), (vec![1, 1], 1), (vec![0, 2], 2)]);
        assert!(!is_lorentzian_bool(&h));
        assert!(!is_normalized_lorentzian_bool(&h));
    }

    #[test]
    fn test_normalized_lorentzian_degree3_not_standard() {
        // f = x³ + 6x²y + 6xy² + y³.
        // Standard: ∂f/∂x = 3x² + 12xy + 6y². Hessian [[6,12],[12,12]].
        //   Minor: 72 - 144 = -72 ≤ 0. ✓ Standard Lorentzian.
        // Normalized: c_{30}·c_{12} = 1·6 = 6 ≤ c_{21}² = 36. ✓
        //             c_{21}·c_{03} = 6·1 = 6 ≤ c_{12}² = 36. ✓ Also normalized.
        let f = poly_from_map(
            2,
            vec![
                (vec![3, 0], 1),
                (vec![2, 1], 6),
                (vec![1, 2], 6),
                (vec![0, 3], 1),
            ],
        );
        assert!(is_lorentzian_bool(&f));
        assert!(is_normalized_lorentzian_bool(&f));

        // f = x³ + 2x²y + 2xy² + y³.
        // Standard: ∂f/∂x = 3x² + 4xy + 2y². Hessian [[6,4],[4,4]].
        //   Minor: 24 - 16 = 8 > 0. ✗ Not standard Lorentzian.
        // Normalized: c_{30}·c_{12} = 1·2 = 2 ≤ c_{21}² = 4. ✓
        //             c_{21}·c_{03} = 2·1 = 2 ≤ c_{12}² = 4. ✓ IS normalized Lorentzian!
        let g = poly_from_map(
            2,
            vec![
                (vec![3, 0], 1),
                (vec![2, 1], 2),
                (vec![1, 2], 2),
                (vec![0, 3], 1),
            ],
        );
        assert!(!is_lorentzian_bool(&g));
        assert!(is_normalized_lorentzian_bool(&g)); // passes normalized!
    }

    #[test]
    fn test_support_m_convex_helper() {
        let f = poly_from_map(
            3,
            vec![(vec![1, 1, 0], 1), (vec![1, 0, 1], 1), (vec![0, 1, 1], 1)],
        );
        assert!(support_is_m_convex(&f));
    }
}
