//! Find linear recurrences among sequences of polynomials.
//!
//! Given polynomials P_1(t), P_2(t), ..., P_m(t), searches for a recurrence:
//!
//!   f(n,t) P_n(t) = sum_{r,d} c_{r,d}(n,t) D^d P_{n-r}(t)  [+ g(n,t)]
//!
//! where c_{r,d}(n,t) are polynomial coefficients in n and t, D^d is the
//! d-th derivative in t, and f(n,t) is an optional LHS denominator.  When
//! `alternating_sign` is enabled, the search also allows terms of the form
//! (-1)^n c_{r,d}(n,t) D^d P_{n-r}(t).
//!
//! This reduces to solving a linear system over the rationals.

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};
use std::fmt;

pub type BigRational = Ratio<BigInt>;

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------

/// Options controlling the search space for `find_polynomial_recurrence`.
#[derive(Debug, Clone)]
pub struct RecurrenceOptions {
    /// Max degree of coefficient polynomials in the variable t.
    pub var_deg: usize,
    /// Max degree of coefficient polynomials in the index n.
    pub idx_deg: usize,
    /// Max order of differentiation (0 = no derivatives).
    pub diff_deg: usize,
    /// How many previous terms the recurrence may use.
    pub rec_len: usize,
    /// If false, allow an additive inhomogeneous polynomial g(n,t).
    pub homogeneous: bool,
    /// Max degree in t of the inhomogeneous term g(n,t).
    ///
    /// Only used when `homogeneous == false`.
    pub inhomo_var_deg: usize,
    /// Max degree in n of the inhomogeneous term g(n,t).
    ///
    /// Only used when `homogeneous == false`.
    pub inhomo_idx_deg: usize,
    /// Degree in t of the LHS factor f(n,t) beyond the implicit constant 1.
    pub denom_var_deg: usize,
    /// Degree in n of the LHS factor f(n,t) beyond the implicit constant 1.
    pub denom_idx_deg: usize,
    /// Also allow recurrence coefficient terms multiplied by (-1)^n.
    pub alternating_sign: bool,
    /// Use modular consistency checks to reject candidates before exact solving.
    ///
    /// This is a probabilistic prefilter: bad primes can make a rationally
    /// solvable system look inconsistent modulo every tested prime. Keep it
    /// disabled for fully exact exhaustive searches.
    pub modular_prefilter: bool,
}

impl Default for RecurrenceOptions {
    fn default() -> Self {
        Self {
            var_deg: 1,
            idx_deg: 1,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            inhomo_var_deg: 1,
            inhomo_idx_deg: 1,
            denom_var_deg: 0,
            denom_idx_deg: 0,
            alternating_sign: false,
            modular_prefilter: false,
        }
    }
}

impl RecurrenceOptions {
    /// Return these options with alternating `(-1)^n` recurrence terms enabled
    /// or disabled.
    pub fn with_alternating_sign(mut self, alternating_sign: bool) -> Self {
        self.alternating_sign = alternating_sign;
        self
    }
}

// ---------------------------------------------------------------------------
// Data types
// ---------------------------------------------------------------------------

/// A polynomial in two variables (n, t).
///
/// Stored as `coeffs[i][j]` = coefficient of n^i t^j.
#[derive(Debug, Clone)]
pub struct BivarPoly {
    pub coeffs: Vec<Vec<BigRational>>,
}

/// Extra sign factor attached to one recurrence term.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum RecurrenceSign {
    /// No extra sign factor.
    #[default]
    None,
    /// Multiply the term by (-1)^n.
    AlternatingN,
}

impl RecurrenceSign {
    fn from_family_index(sign_idx: usize) -> Self {
        match sign_idx {
            0 => Self::None,
            1 => Self::AlternatingN,
            _ => unreachable!("invalid recurrence sign family"),
        }
    }
}

/// One term in a recurrence: sign(n) * coeff(n,t) * D^d P_{n-r}(t).
#[derive(Debug, Clone)]
pub struct RecurrenceTerm {
    /// Recurrence offset r (P_{n-r}).
    pub offset: usize,
    /// Derivative order d.
    pub deriv_order: usize,
    /// Optional multiplicative sign factor.
    pub sign: RecurrenceSign,
    /// Coefficient polynomial c(n,t).
    pub coeff: BivarPoly,
}

impl RecurrenceTerm {
    /// Construct an ordinary non-alternating recurrence term.
    pub fn new(offset: usize, deriv_order: usize, coeff: BivarPoly) -> Self {
        Self {
            offset,
            deriv_order,
            sign: RecurrenceSign::None,
            coeff,
        }
    }

    /// Return this term with a different sign factor.
    pub fn with_sign(mut self, sign: RecurrenceSign) -> Self {
        self.sign = sign;
        self
    }

    /// True when this term is multiplied by `(-1)^n`.
    pub fn is_alternating(&self) -> bool {
        self.sign == RecurrenceSign::AlternatingN
    }
}

/// A polynomial recurrence found by `find_polynomial_recurrence`.
#[derive(Debug, Clone)]
pub struct Recurrence {
    /// RHS terms: each is c(n,t) D^d P_{n-r}.
    pub terms: Vec<RecurrenceTerm>,
    /// LHS factor f(n,t) if non-trivial (None means f = 1).
    pub denominator: Option<BivarPoly>,
    /// Inhomogeneous additive term g(n,t), if present.
    pub inhomogeneous: Option<BivarPoly>,
}

// ---------------------------------------------------------------------------
// Polynomial helpers (univariate, coefficient-vector representation)
// ---------------------------------------------------------------------------

/// Convert an `i64` coefficient vector to exact rational coefficients.
fn i64_poly_to_rational(coeffs: &[i64]) -> Vec<BigRational> {
    coeffs
        .iter()
        .map(|&c| BigRational::from_integer(BigInt::from(c)))
        .collect()
}

fn i64_polys_to_rational(polys: &[Vec<i64>]) -> Vec<Vec<BigRational>> {
    polys.iter().map(|p| i64_poly_to_rational(p)).collect()
}

fn rational_zero_poly() -> Vec<BigRational> {
    vec![BigRational::zero()]
}

/// d-th derivative of a polynomial given as a rational coefficient vector.
fn poly_nth_derivative_rational(coeffs: &[BigRational], d: usize) -> Vec<BigRational> {
    let mut result = coeffs.to_vec();
    for _ in 0..d {
        if result.len() <= 1 {
            return rational_zero_poly();
        }
        result = result[1..]
            .iter()
            .enumerate()
            .map(|(i, c)| c.clone() * BigRational::from_integer(BigInt::from(i + 1)))
            .collect();
    }
    if result.is_empty() {
        rational_zero_poly()
    } else {
        result
    }
}

/// Coefficient of t^k (0 when k is out of range).
fn poly_coeff_rational(coeffs: &[BigRational], k: usize) -> BigRational {
    if k < coeffs.len() {
        coeffs[k].clone()
    } else {
        BigRational::zero()
    }
}

/// Degree of a polynomial (returns 0 for the zero polynomial).
fn poly_degree_rational(coeffs: &[BigRational]) -> usize {
    coeffs.iter().rposition(|c| !c.is_zero()).unwrap_or(0)
}

fn trim_poly_rational(mut coeffs: Vec<BigRational>) -> Vec<BigRational> {
    while coeffs.len() > 1 && coeffs.last().is_some_and(|c| c.is_zero()) {
        coeffs.pop();
    }
    if coeffs.is_empty() {
        rational_zero_poly()
    } else {
        coeffs
    }
}

fn poly_is_zero_rational(coeffs: &[BigRational]) -> bool {
    coeffs.iter().all(|c| c.is_zero())
}

fn poly_equal_rational(lhs: &[BigRational], rhs: &[BigRational]) -> bool {
    let max_len = lhs.len().max(rhs.len());
    (0..max_len).all(|idx| poly_coeff_rational(lhs, idx) == poly_coeff_rational(rhs, idx))
}

fn poly_add_scaled_assign(
    target: &mut Vec<BigRational>,
    source: &[BigRational],
    scale: &BigRational,
) {
    if scale.is_zero() || poly_is_zero_rational(source) {
        return;
    }
    if target.len() < source.len() {
        target.resize(source.len(), BigRational::zero());
    }
    for (idx, coeff) in source.iter().enumerate() {
        if !coeff.is_zero() {
            target[idx] += coeff.clone() * scale.clone();
        }
    }
}

fn poly_mul_rational(lhs: &[BigRational], rhs: &[BigRational]) -> Vec<BigRational> {
    if poly_is_zero_rational(lhs) || poly_is_zero_rational(rhs) {
        return rational_zero_poly();
    }
    let mut result = vec![BigRational::zero(); lhs.len() + rhs.len() - 1];
    for (i, a) in lhs.iter().enumerate() {
        if a.is_zero() {
            continue;
        }
        for (j, b) in rhs.iter().enumerate() {
            if !b.is_zero() {
                result[i + j] += a.clone() * b.clone();
            }
        }
    }
    trim_poly_rational(result)
}

fn poly_div_exact_rational(
    numerator: &[BigRational],
    denominator: &[BigRational],
) -> Result<Vec<BigRational>, RecurrenceEvaluationError> {
    let denominator = trim_poly_rational(denominator.to_vec());
    if poly_is_zero_rational(&denominator) {
        return Err(RecurrenceEvaluationError::ZeroDenominator);
    }

    let mut remainder = trim_poly_rational(numerator.to_vec());
    if poly_is_zero_rational(&remainder) {
        return Ok(rational_zero_poly());
    }

    let denom_deg = poly_degree_rational(&denominator);
    let denom_lc = denominator[denom_deg].clone();
    let mut quotient =
        vec![BigRational::zero(); poly_degree_rational(&remainder).saturating_sub(denom_deg) + 1];

    while !poly_is_zero_rational(&remainder) && poly_degree_rational(&remainder) >= denom_deg {
        let rem_deg = poly_degree_rational(&remainder);
        let shift = rem_deg - denom_deg;
        let factor = remainder[rem_deg].clone() / denom_lc.clone();
        quotient[shift] += factor.clone();
        for (idx, denom_coeff) in denominator.iter().enumerate().take(denom_deg + 1) {
            let target = idx + shift;
            remainder[target] -= factor.clone() * denom_coeff.clone();
        }
        remainder = trim_poly_rational(remainder);
    }

    if poly_is_zero_rational(&remainder) {
        Ok(trim_poly_rational(quotient))
    } else {
        Err(RecurrenceEvaluationError::NonPolynomialQuotient)
    }
}

fn bivar_eval_n(poly: &BivarPoly, n: usize) -> Vec<BigRational> {
    bivar_eval_n_with_powers(poly, &rational_index_powers(n, bivar_n_degree(poly)))
}

fn bivar_eval_n_with_powers(poly: &BivarPoly, n_powers: &[BigRational]) -> Vec<BigRational> {
    let width = poly.coeffs.iter().map(Vec::len).max().unwrap_or(0);
    if width == 0 {
        return rational_zero_poly();
    }
    let mut result = vec![BigRational::zero(); width];
    for (i, row) in poly.coeffs.iter().enumerate() {
        let n_power = n_powers
            .get(i)
            .expect("bivariate evaluation requires enough powers of n");
        for (j, coeff) in row.iter().enumerate() {
            if !coeff.is_zero() {
                result[j] += coeff.clone() * n_power.clone();
            }
        }
    }
    trim_poly_rational(result)
}

fn bivar_n_degree(poly: &BivarPoly) -> usize {
    poly.coeffs.len().saturating_sub(1)
}

fn recurrence_n_degree(rec: &Recurrence) -> usize {
    let mut degree = rec.denominator.as_ref().map_or(0, bivar_n_degree);
    for term in &rec.terms {
        degree = degree.max(bivar_n_degree(&term.coeff));
    }
    if let Some(inhomogeneous) = &rec.inhomogeneous {
        degree = degree.max(bivar_n_degree(inhomogeneous));
    }
    degree
}

use crate::linalg;

const MODULAR_PREFILTER_PRIMES: [i64; 3] = [1_000_000_007, 1_000_000_009, 998_244_353];
const MODULAR_PREFILTER_INCONSISTENT_PRIMES_TO_REJECT: usize = 2;
// Full-rank modular fits first check one held-out row, which cheaply rejects
// many false prefix fits.  A full held-out sweep is reserved for larger spaces.
const MODULAR_PREFILTER_FIRST_TAIL_VERIFY_MIN_UNKNOWNS: usize = 1;
const MODULAR_PREFILTER_TAIL_VERIFY_MIN_UNKNOWNS: usize = 128;

fn mod_norm(value: i64, modulus: i64) -> i64 {
    value.rem_euclid(modulus)
}

fn mod_neg(value: i64, modulus: i64) -> i64 {
    if value == 0 {
        0
    } else {
        modulus - value
    }
}

fn mod_mul(lhs: i64, rhs: i64, modulus: i64) -> i64 {
    if (0..modulus).contains(&lhs) && (0..modulus).contains(&rhs) && modulus <= 3_000_000_000 {
        return (lhs * rhs) % modulus;
    }
    ((lhs as i128 * rhs as i128).rem_euclid(modulus as i128)) as i64
}

fn mod_pow(mut base: i64, mut exp: i64, modulus: i64) -> i64 {
    let mut acc = 1;
    base = mod_norm(base, modulus);
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mod_mul(acc, base, modulus);
        }
        base = mod_mul(base, base, modulus);
        exp >>= 1;
    }
    acc
}

fn mod_index_powers(n: usize, max_power: usize, modulus: i64) -> Vec<i64> {
    let base = mod_norm(n as i64, modulus);
    let mut powers = Vec::with_capacity(max_power + 1);
    let mut current = 1;
    for _ in 0..=max_power {
        powers.push(current);
        current = mod_mul(current, base, modulus);
    }
    powers
}

fn rational_index_powers(n: usize, max_power: usize) -> Vec<BigRational> {
    let base = BigInt::from(n);
    let mut powers = Vec::with_capacity(max_power + 1);
    let mut current = BigInt::one();
    for _ in 0..=max_power {
        powers.push(BigRational::from_integer(current.clone()));
        current *= &base;
    }
    powers
}

fn mod_inv(value: i64, modulus: i64) -> Option<i64> {
    let value = mod_norm(value, modulus);
    (value != 0).then(|| mod_pow(value, modulus - 2, modulus))
}

fn bigint_mod_i64(value: &BigInt, modulus: i64) -> i64 {
    let modulus_big = BigInt::from(modulus);
    let residue = ((value % &modulus_big) + &modulus_big) % &modulus_big;
    residue
        .to_i64()
        .expect("residue modulo i64 prime fits in i64")
}

fn rational_mod_prime(value: &BigRational, modulus: i64) -> Option<i64> {
    let numerator = bigint_mod_i64(value.numer(), modulus);
    let denominator = bigint_mod_i64(value.denom(), modulus);
    let denominator_inv = mod_inv(denominator, modulus)?;
    Some(mod_mul(numerator, denominator_inv, modulus))
}

fn poly_coeff_mod(coeffs: &[Vec<i64>], poly_idx: usize, coeff_idx: usize) -> i64 {
    coeffs
        .get(poly_idx)
        .and_then(|coeffs| coeffs.get(coeff_idx))
        .copied()
        .unwrap_or(0)
}

fn poly_coeff_mod_slice(coeffs: &[i64], coeff_idx: usize) -> i64 {
    coeffs.get(coeff_idx).copied().unwrap_or(0)
}

fn add_modular_solution_term(
    lhs: &mut i64,
    solution: &[u64],
    col: usize,
    coeff: i64,
    modulus: i64,
) {
    let coeff = mod_norm(coeff, modulus);
    if coeff == 0 {
        return;
    }
    let term = mod_mul(coeff, solution[col] as i64, modulus);
    *lhs += term;
    if *lhs >= modulus {
        *lhs -= modulus;
    }
}

fn rational_polys_mod_prime(polys: &[Vec<BigRational>], modulus: i64) -> Option<Vec<Vec<i64>>> {
    polys
        .iter()
        .map(|poly| {
            poly.iter()
                .map(|coeff| rational_mod_prime(coeff, modulus))
                .collect()
        })
        .collect()
}

fn rational_derivs_mod_prime(
    derivs: &[Vec<Vec<BigRational>>],
    modulus: i64,
    max_diff_deg: usize,
) -> Option<Vec<Vec<Vec<i64>>>> {
    derivs
        .iter()
        .map(|poly_derivs| {
            poly_derivs
                .iter()
                .take(max_diff_deg + 1)
                .map(|poly| {
                    poly.iter()
                        .map(|coeff| rational_mod_prime(coeff, modulus))
                        .collect()
                })
                .collect()
        })
        .collect()
}

struct ModularSystemResult {
    consistent: bool,
    full_rank_pivot_rows: Option<Vec<usize>>,
    full_rank_solution: Option<Vec<u64>>,
}

fn recurrence_system_consistent_mod_images(
    polys: &[Vec<BigRational>],
    polys_mod: &[Vec<i64>],
    derivs_mod: &[Vec<Vec<i64>>],
    opts: &RecurrenceOptions,
    modulus: i64,
) -> ModularSystemResult {
    let m = polys.len();
    if m <= opts.rec_len {
        return ModularSystemResult {
            consistent: true,
            full_rank_pivot_rows: None,
            full_rank_solution: None,
        };
    }

    let denom_start: usize = 0;
    let denom_w = opts.denom_var_deg + 1;
    let num_denom_vars = (opts.denom_idx_deg + 1) * denom_w - 1;
    let denom_col = |i: usize, j: usize| -> usize {
        let flat = i * denom_w + j;
        denom_start + flat - 1
    };

    let coeff_start = denom_start + num_denom_vars;
    let vars_per_coeff = (opts.idx_deg + 1) * (opts.var_deg + 1);
    let sign_family_count = if opts.alternating_sign { 2 } else { 1 };
    let num_coeff_vars = opts.rec_len * (opts.diff_deg + 1) * sign_family_count * vars_per_coeff;
    let coeff_col = |r: usize, d: usize, sign_idx: usize, i: usize, j: usize| -> usize {
        coeff_start
            + (((r - 1) * (opts.diff_deg + 1) + d) * sign_family_count + sign_idx) * vars_per_coeff
            + i * (opts.var_deg + 1)
            + j
    };

    let inhomo_start = coeff_start + num_coeff_vars;
    let inhomo_w = opts.inhomo_var_deg + 1;
    let num_inhomo_vars = if opts.homogeneous {
        0
    } else {
        (opts.inhomo_idx_deg + 1) * inhomo_w
    };
    let inhomo_col = |i: usize, j: usize| -> usize { inhomo_start + i * inhomo_w + j };

    let num_vars = inhomo_start + num_inhomo_vars;
    if num_vars == 0 {
        return ModularSystemResult {
            consistent: true,
            full_rank_pivot_rows: None,
            full_rank_solution: None,
        };
    }
    let prime = modulus as u64;

    let max_poly_deg = polys
        .iter()
        .map(|p| poly_degree_rational(p))
        .max()
        .unwrap_or(0);
    let max_j = opts
        .var_deg
        .max(opts.denom_var_deg)
        .max(if opts.homogeneous {
            0
        } else {
            opts.inhomo_var_deg
        });
    let max_t_deg = max_j + max_poly_deg;
    let eqs_per_nn = max_t_deg + 1;
    let mut rows = Vec::with_capacity((m - opts.rec_len) * eqs_per_nn);
    let mut row_indices = Vec::with_capacity((m - opts.rec_len) * eqs_per_nn);
    let max_n_pow = opts
        .idx_deg
        .max(opts.denom_idx_deg)
        .max(if opts.homogeneous {
            0
        } else {
            opts.inhomo_idx_deg
        });

    for (eq_idx, nn) in (opts.rec_len + 1..=m).enumerate() {
        let current_idx = nn - 1;
        let n_powers = mod_index_powers(nn, max_n_pow, modulus);

        for l in 0..=max_t_deg {
            let equation_row = eq_idx * eqs_per_nn + l;
            let cur_l = poly_coeff_mod(polys_mod, current_idx, l);
            let mut row = linalg::SparseModRow::new(mod_neg(cur_l, modulus) as u64, prime);

            if num_denom_vars > 0 {
                for (i, &n_power) in n_powers.iter().enumerate().take(opts.denom_idx_deg + 1) {
                    for j in 0..=opts.denom_var_deg {
                        if i == 0 && j == 0 {
                            continue;
                        }
                        if l < j {
                            continue;
                        }
                        let pc = poly_coeff_mod(polys_mod, current_idx, l - j);
                        if pc == 0 {
                            continue;
                        }
                        let col = denom_col(i, j);
                        row.add_reduced_entry_sorted(
                            col,
                            mod_mul(pc, n_power, modulus) as u64,
                            prime,
                        );
                    }
                }
            }

            for r in 1..=opts.rec_len {
                for (d, ref_poly) in derivs_mod[nn - 1 - r]
                    .iter()
                    .enumerate()
                    .take(opts.diff_deg + 1)
                {
                    for sign_idx in 0..sign_family_count {
                        for (i, &n_power) in n_powers.iter().enumerate().take(opts.idx_deg + 1) {
                            let n_factor = if sign_idx == 1 && nn % 2 == 1 {
                                mod_neg(n_power, modulus)
                            } else {
                                n_power
                            };
                            for j in 0..=opts.var_deg {
                                if l < j {
                                    continue;
                                }
                                let rc = poly_coeff_mod_slice(ref_poly, l - j);
                                if rc == 0 {
                                    continue;
                                }
                                let col = coeff_col(r, d, sign_idx, i, j);
                                row.add_reduced_entry_sorted(
                                    col,
                                    mod_neg(mod_mul(rc, n_factor, modulus), modulus) as u64,
                                    prime,
                                );
                            }
                        }
                    }
                }
            }

            if !opts.homogeneous && l <= opts.inhomo_var_deg {
                for (i, &n_power) in n_powers.iter().enumerate().take(opts.inhomo_idx_deg + 1) {
                    let col = inhomo_col(i, l);
                    row.add_reduced_entry_sorted(col, mod_neg(n_power, modulus) as u64, prime);
                }
            }

            if row.is_contradiction() {
                return ModularSystemResult {
                    consistent: false,
                    full_rank_pivot_rows: None,
                    full_rank_solution: None,
                };
            }
            if !row.is_tautology() {
                row_indices.push(equation_row);
                rows.push(row);
            }
        }
    }

    let row_order = linalg::SparseModRowOrder::IncreasingNonzeros;
    let result =
        linalg::sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode(
            rows,
            num_vars,
            prime,
            row_order,
            linalg::SparseModSolutionMode::FullRank,
        )
        .expect("recurrence modular prefilter builds a well-formed prime-field system");
    let full_rank = result.consistent && result.rank == num_vars;
    let full_rank_pivot_rows = full_rank.then(|| {
        result
            .pivot_rows
            .iter()
            .map(|&row| row_indices[row])
            .collect()
    });
    let full_rank_solution = full_rank.then_some(result.solution).flatten();
    ModularSystemResult {
        consistent: result.consistent,
        full_rank_pivot_rows,
        full_rank_solution,
    }
}

fn recurrence_solution_holds_from_mod_images(
    polys: &[Vec<BigRational>],
    polys_mod: &[Vec<i64>],
    derivs_mod: &[Vec<Vec<i64>>],
    opts: &RecurrenceOptions,
    modulus: i64,
    solution: &[u64],
    start_nn: usize,
) -> bool {
    let m = polys.len();
    let first_nn = start_nn.max(opts.rec_len + 1);
    if first_nn > m {
        return true;
    }

    let denom_start: usize = 0;
    let denom_w = opts.denom_var_deg + 1;
    let num_denom_vars = (opts.denom_idx_deg + 1) * denom_w - 1;
    let denom_col = |i: usize, j: usize| -> usize {
        let flat = i * denom_w + j;
        denom_start + flat - 1
    };

    let coeff_start = denom_start + num_denom_vars;
    let vars_per_coeff = (opts.idx_deg + 1) * (opts.var_deg + 1);
    let sign_family_count = if opts.alternating_sign { 2 } else { 1 };
    let num_coeff_vars = opts.rec_len * (opts.diff_deg + 1) * sign_family_count * vars_per_coeff;
    let coeff_col = |r: usize, d: usize, sign_idx: usize, i: usize, j: usize| -> usize {
        coeff_start
            + (((r - 1) * (opts.diff_deg + 1) + d) * sign_family_count + sign_idx) * vars_per_coeff
            + i * (opts.var_deg + 1)
            + j
    };

    let inhomo_start = coeff_start + num_coeff_vars;
    let inhomo_w = opts.inhomo_var_deg + 1;
    let num_inhomo_vars = if opts.homogeneous {
        0
    } else {
        (opts.inhomo_idx_deg + 1) * inhomo_w
    };
    let inhomo_col = |i: usize, j: usize| -> usize { inhomo_start + i * inhomo_w + j };

    let num_vars = inhomo_start + num_inhomo_vars;
    if solution.len() != num_vars {
        return false;
    }

    let max_poly_deg = polys
        .iter()
        .map(|p| poly_degree_rational(p))
        .max()
        .unwrap_or(0);
    let max_j = opts
        .var_deg
        .max(opts.denom_var_deg)
        .max(if opts.homogeneous {
            0
        } else {
            opts.inhomo_var_deg
        });
    let max_t_deg = max_j + max_poly_deg;
    let max_n_pow = opts
        .idx_deg
        .max(opts.denom_idx_deg)
        .max(if opts.homogeneous {
            0
        } else {
            opts.inhomo_idx_deg
        });

    for nn in first_nn..=m {
        let current_idx = nn - 1;
        let n_powers = mod_index_powers(nn, max_n_pow, modulus);

        for l in 0..=max_t_deg {
            let cur_l = poly_coeff_mod(polys_mod, current_idx, l);
            let rhs = mod_neg(cur_l, modulus);
            let mut lhs = 0;

            if num_denom_vars > 0 {
                for (i, &n_power) in n_powers.iter().enumerate().take(opts.denom_idx_deg + 1) {
                    for j in 0..=opts.denom_var_deg {
                        if i == 0 && j == 0 {
                            continue;
                        }
                        if l < j {
                            continue;
                        }
                        let pc = poly_coeff_mod(polys_mod, current_idx, l - j);
                        if pc == 0 {
                            continue;
                        }
                        add_modular_solution_term(
                            &mut lhs,
                            solution,
                            denom_col(i, j),
                            mod_mul(pc, n_power, modulus),
                            modulus,
                        );
                    }
                }
            }

            for r in 1..=opts.rec_len {
                for (d, ref_poly) in derivs_mod[nn - 1 - r]
                    .iter()
                    .enumerate()
                    .take(opts.diff_deg + 1)
                {
                    for sign_idx in 0..sign_family_count {
                        for (i, &n_power) in n_powers.iter().enumerate().take(opts.idx_deg + 1) {
                            let n_factor = if sign_idx == 1 && nn % 2 == 1 {
                                mod_neg(n_power, modulus)
                            } else {
                                n_power
                            };
                            for j in 0..=opts.var_deg {
                                if l < j {
                                    continue;
                                }
                                let rc = poly_coeff_mod_slice(ref_poly, l - j);
                                if rc == 0 {
                                    continue;
                                }
                                add_modular_solution_term(
                                    &mut lhs,
                                    solution,
                                    coeff_col(r, d, sign_idx, i, j),
                                    mod_neg(mod_mul(rc, n_factor, modulus), modulus),
                                    modulus,
                                );
                            }
                        }
                    }
                }
            }

            if !opts.homogeneous && l <= opts.inhomo_var_deg {
                for (i, &n_power) in n_powers.iter().enumerate().take(opts.inhomo_idx_deg + 1) {
                    add_modular_solution_term(
                        &mut lhs,
                        solution,
                        inhomo_col(i, l),
                        mod_neg(n_power, modulus),
                        modulus,
                    );
                }
            }

            if lhs != rhs {
                return false;
            }
        }
    }

    true
}

fn recurrence_system_consistent_mod_prime(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
    modulus: i64,
) -> Option<bool> {
    let polys_mod = rational_polys_mod_prime(polys, modulus)?;
    let derivs_mod = rational_derivs_mod_prime(derivs, modulus, opts.diff_deg)?;
    Some(
        recurrence_system_consistent_mod_images(polys, &polys_mod, &derivs_mod, opts, modulus)
            .consistent,
    )
}

struct ModularPrefilterPrimeImages {
    modulus: i64,
    polys: Vec<Vec<i64>>,
    derivs: Vec<Vec<Vec<i64>>>,
}

struct ModularPrefilterCache {
    primes: Vec<ModularPrefilterPrimeImages>,
}

impl ModularPrefilterCache {
    fn new(
        polys: &[Vec<BigRational>],
        derivs: &[Vec<Vec<BigRational>>],
        max_diff_deg: usize,
    ) -> Self {
        let primes = MODULAR_PREFILTER_PRIMES
            .iter()
            .filter_map(|&modulus| {
                let polys_mod = rational_polys_mod_prime(polys, modulus)?;
                let derivs_mod = rational_derivs_mod_prime(derivs, modulus, max_diff_deg)?;
                Some(ModularPrefilterPrimeImages {
                    modulus,
                    polys: polys_mod,
                    derivs: derivs_mod,
                })
            })
            .collect();
        Self { primes }
    }
}

fn modular_prefilter_rejects(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
) -> bool {
    if !opts.modular_prefilter {
        return false;
    }

    let mut inconsistent_primes = 0;
    for &prime in &MODULAR_PREFILTER_PRIMES {
        if let Some(consistent) = recurrence_system_consistent_mod_prime(polys, derivs, opts, prime)
        {
            if !consistent {
                inconsistent_primes += 1;
                if inconsistent_primes >= MODULAR_PREFILTER_INCONSISTENT_PRIMES_TO_REJECT {
                    return true;
                }
            }
        }
    }

    false
}

#[derive(Default)]
struct ModularPrefilterResult {
    rejected: bool,
    full_rank_pivot_rows: Option<Vec<usize>>,
}

fn modular_prefilter_with_cache(
    polys: &[Vec<BigRational>],
    fit_len: usize,
    opts: &RecurrenceOptions,
    cache: &ModularPrefilterCache,
) -> ModularPrefilterResult {
    if !opts.modular_prefilter {
        return ModularPrefilterResult::default();
    }

    let mut inconsistent_primes = 0;
    let unknowns = count_unknowns(opts);
    let tail_verify_end =
        if fit_len >= polys.len() || unknowns < MODULAR_PREFILTER_FIRST_TAIL_VERIFY_MIN_UNKNOWNS {
            None
        } else if unknowns >= MODULAR_PREFILTER_TAIL_VERIFY_MIN_UNKNOWNS {
            Some(polys.len())
        } else {
            Some(fit_len + 1)
        };
    for prime_images in &cache.primes {
        if fit_len > polys.len()
            || polys.len() > prime_images.polys.len()
            || polys.len() > prime_images.derivs.len()
        {
            continue;
        }
        if prime_images
            .derivs
            .first()
            .is_some_and(|derivs| derivs.len() <= opts.diff_deg)
        {
            continue;
        }

        let result = recurrence_system_consistent_mod_images(
            &polys[..fit_len],
            &prime_images.polys[..fit_len],
            &prime_images.derivs[..fit_len],
            opts,
            prime_images.modulus,
        );
        let tail_consistent = match (tail_verify_end, result.full_rank_solution.as_ref()) {
            (Some(verify_len), Some(solution)) => recurrence_solution_holds_from_mod_images(
                &polys[..verify_len],
                &prime_images.polys[..verify_len],
                &prime_images.derivs[..verify_len],
                opts,
                prime_images.modulus,
                solution,
                fit_len + 1,
            ),
            _ => true,
        };
        if result.consistent && !tail_consistent {
            // A full-rank modular fit has a unique solution with a nonzero
            // pivot determinant modulo this prime. If a held-out row fails,
            // the exact candidate cannot satisfy the same held-out equation.
            return ModularPrefilterResult {
                rejected: true,
                full_rank_pivot_rows: None,
            };
        }
        if !result.consistent {
            inconsistent_primes += 1;
            if inconsistent_primes >= MODULAR_PREFILTER_INCONSISTENT_PRIMES_TO_REJECT {
                return ModularPrefilterResult {
                    rejected: true,
                    full_rank_pivot_rows: None,
                };
            }
        } else if result.full_rank_pivot_rows.is_none() {
            if let Some(verify_len) = tail_verify_end {
                if verify_len > fit_len {
                    // Rank-deficient prefix fits are the expensive false
                    // positives: they have many modular solutions, so no
                    // unique held-out solution can be checked.  Add the
                    // held-out rows modulo the same prime before paying for an
                    // exact rational solve.
                    let extended = recurrence_system_consistent_mod_images(
                        &polys[..verify_len],
                        &prime_images.polys[..verify_len],
                        &prime_images.derivs[..verify_len],
                        opts,
                        prime_images.modulus,
                    );
                    if !extended.consistent {
                        inconsistent_primes += 1;
                        if inconsistent_primes >= MODULAR_PREFILTER_INCONSISTENT_PRIMES_TO_REJECT {
                            return ModularPrefilterResult {
                                rejected: true,
                                full_rank_pivot_rows: None,
                            };
                        }
                    }
                }
            }
        } else if let Some(pivot_rows) = result.full_rank_pivot_rows {
            return ModularPrefilterResult {
                rejected: false,
                full_rank_pivot_rows: Some(pivot_rows),
            };
        }
    }

    ModularPrefilterResult {
        rejected: false,
        full_rank_pivot_rows: None,
    }
}

// ---------------------------------------------------------------------------
// Main algorithm
// ---------------------------------------------------------------------------

fn rational_derivatives_up_to(
    polys: &[Vec<BigRational>],
    max_diff_deg: usize,
) -> Vec<Vec<Vec<BigRational>>> {
    polys
        .iter()
        .map(|p| {
            (0..=max_diff_deg)
                .map(|d| poly_nth_derivative_rational(p, d))
                .collect()
        })
        .collect()
}

/// Search for a polynomial recurrence satisfied by a sequence of polynomials.
///
/// `polys[i]` is P_{i+1}(t) given as a coefficient vector (index = power of t).
/// Returns `None` if no recurrence is found with the given options.
pub fn find_polynomial_recurrence_rational(
    polys: &[Vec<BigRational>],
    opts: &RecurrenceOptions,
) -> Option<Recurrence> {
    let derivs = rational_derivatives_up_to(polys, opts.diff_deg);
    find_polynomial_recurrence_rational_with_derivs(polys, &derivs, opts)
}

fn find_polynomial_recurrence_rational_with_derivs(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
) -> Option<Recurrence> {
    find_polynomial_recurrence_rational_with_derivs_and_rows(polys, derivs, opts, None)
}

fn find_polynomial_recurrence_rational_with_derivs_and_rows(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
    exact_row_indices: Option<&[usize]>,
) -> Option<Recurrence> {
    let m = polys.len();
    if m <= opts.rec_len {
        return None;
    }

    if modular_prefilter_rejects(polys, derivs, opts) {
        return None;
    }

    // --- Assign column indices to unknowns ---

    // 1) Denominator unknowns d[i][j], skipping (0,0) which is the fixed constant 1.
    let denom_start: usize = 0;
    let denom_w = opts.denom_var_deg + 1; // width in j
                                          // Total slots minus 1 for the fixed constant-1 at (0,0).
    let num_denom_vars = (opts.denom_idx_deg + 1) * denom_w - 1;

    let denom_col = |i: usize, j: usize| -> usize {
        let flat = i * denom_w + j; // (0,0) has flat=0
        denom_start + flat - 1 // skip (0,0)
    };

    // 2) Recurrence coefficient unknowns c[r][d][sign][i][j].
    let coeff_start = denom_start + num_denom_vars;
    let vars_per_coeff = (opts.idx_deg + 1) * (opts.var_deg + 1);
    let sign_family_count = if opts.alternating_sign { 2 } else { 1 };
    let num_coeff_vars = opts.rec_len * (opts.diff_deg + 1) * sign_family_count * vars_per_coeff;

    let coeff_col = |r: usize, d: usize, sign_idx: usize, i: usize, j: usize| -> usize {
        coeff_start
            + (((r - 1) * (opts.diff_deg + 1) + d) * sign_family_count + sign_idx) * vars_per_coeff
            + i * (opts.var_deg + 1)
            + j
    };

    // 3) Inhomogeneous unknowns (if needed).
    let inhomo_start = coeff_start + num_coeff_vars;
    let inhomo_w = opts.inhomo_var_deg + 1;
    let num_inhomo_vars = if opts.homogeneous {
        0
    } else {
        (opts.inhomo_idx_deg + 1) * inhomo_w
    };

    let inhomo_col = |i: usize, j: usize| -> usize { inhomo_start + i * inhomo_w + j };

    let num_vars = inhomo_start + num_inhomo_vars;
    if num_vars == 0 {
        return None;
    }

    // --- Determine max t-degree across all equations ---
    let max_poly_deg = polys
        .iter()
        .map(|p| poly_degree_rational(p))
        .max()
        .unwrap_or(0);
    let max_j = opts
        .var_deg
        .max(opts.denom_var_deg)
        .max(if opts.homogeneous {
            0
        } else {
            opts.inhomo_var_deg
        });
    let max_t_deg = max_j + max_poly_deg;
    let eqs_per_nn = max_t_deg + 1;

    // --- Build linear system Ax = b ---
    let num_nn = m - opts.rec_len; // nn = rec_len+1 ..= m  (1-based)
    let num_rows = num_nn * eqs_per_nn;
    let all_row_indices;
    let row_indices = if let Some(row_indices) = exact_row_indices {
        if row_indices.is_empty() {
            return None;
        }
        row_indices
    } else {
        all_row_indices = (0..num_rows).collect::<Vec<_>>();
        &all_row_indices
    };

    let zero = BigRational::zero();
    let mut matrix: Vec<Vec<BigRational>> = vec![vec![zero.clone(); num_vars]; row_indices.len()];
    let mut rhs: Vec<BigRational> = vec![zero.clone(); row_indices.len()];
    let max_n_pow = opts
        .idx_deg
        .max(opts.denom_idx_deg)
        .max(if opts.homogeneous {
            0
        } else {
            opts.inhomo_idx_deg
        });
    let n_power_table = (opts.rec_len + 1..=m)
        .map(|nn| rational_index_powers(nn, max_n_pow))
        .collect::<Vec<_>>();

    for (row, &equation_row) in row_indices.iter().enumerate() {
        if equation_row >= num_rows {
            return None;
        }
        let eq_idx = equation_row / eqs_per_nn;
        let l = equation_row % eqs_per_nn;
        let nn = opts.rec_len + 1 + eq_idx;
        // nn is 1-based; polys[nn-1] is P_nn(t).
        let current = &polys[nn - 1];
        let n_powers = &n_power_table[eq_idx];

        // RHS = coefficient of t^l in P_nn(t)  (moved to RHS with negation).
        let cur_l = poly_coeff_rational(current, l);
        rhs[row] = -cur_l;

        // Denominator unknowns: d[i][j] * nn^i * (coeff of t^{l-j} in P_nn).
        if num_denom_vars > 0 {
            for i in 0..=opts.denom_idx_deg {
                for j in 0..=opts.denom_var_deg {
                    if i == 0 && j == 0 {
                        continue;
                    }
                    if l < j {
                        continue;
                    }
                    let pc = poly_coeff_rational(current, l - j);
                    if pc.is_zero() {
                        continue;
                    }
                    let val = pc * n_powers[i].clone();
                    matrix[row][denom_col(i, j)] = val;
                }
            }
        }

        // Recurrence coefficients:
        // −sign(nn) * c[r][d][sign][i][j] * nn^i
        //     * (coeff of t^{l-j} in D^d P_{nn-r}).
        for r in 1..=opts.rec_len {
            for d in 0..=opts.diff_deg {
                let ref_poly = &derivs[nn - 1 - r][d];
                for sign_idx in 0..sign_family_count {
                    for i in 0..=opts.idx_deg {
                        let mut ni = n_powers[i].clone();
                        if sign_idx == 1 && nn % 2 == 1 {
                            ni = -ni;
                        }
                        for j in 0..=opts.var_deg {
                            if l < j {
                                continue;
                            }
                            let rc = poly_coeff_rational(ref_poly, l - j);
                            if rc.is_zero() {
                                continue;
                            }
                            let val = -rc * ni.clone();
                            matrix[row][coeff_col(r, d, sign_idx, i, j)] = val;
                        }
                    }
                }
            }
        }

        // Inhomogeneous unknowns: −c_inh[i][j] * nn^i * delta(l,j).
        if !opts.homogeneous {
            for i in 0..=opts.inhomo_idx_deg {
                if l <= opts.inhomo_var_deg {
                    let val = -n_powers[i].clone();
                    matrix[row][inhomo_col(i, l)] += val;
                }
            }
        }
    }

    // --- Solve ---
    let solution = if exact_row_indices.is_some() && matrix.len() == num_vars {
        linalg::solve_full_rank_square_linear_system(&matrix, &rhs)?
    } else {
        linalg::solve_linear_system(&matrix, &rhs)?
    };

    // Check: are all recurrence coefficients zero?  (Trivial / degenerate.)
    let all_zero = (coeff_start..coeff_start + num_coeff_vars).all(|c| solution[c].is_zero());
    if all_zero {
        return None;
    }

    // NOTE: No normalization here. The system is non-homogeneous (the fixed
    // constant 1 in f(n,t) pins the scale), so the solution is uniquely
    // determined — not up to scaling.  Dividing by GCD would lose actual
    // coefficient values (e.g. turning P(n) = 2 P(n-1) into P(n) = P(n-1)).

    // --- Build result ---
    let mut terms = Vec::new();
    for r in 1..=opts.rec_len {
        for d in 0..=opts.diff_deg {
            for sign_idx in 0..sign_family_count {
                let bv = extract_bivar(
                    &solution,
                    |i, j| coeff_col(r, d, sign_idx, i, j),
                    opts.idx_deg,
                    opts.var_deg,
                );
                if !bv.is_zero() {
                    terms.push(RecurrenceTerm {
                        offset: r,
                        deriv_order: d,
                        sign: RecurrenceSign::from_family_index(sign_idx),
                        coeff: bv,
                    });
                }
            }
        }
    }

    let denominator = if num_denom_vars > 0 {
        // Build denominator BivarPoly manually since (0,0) is not a variable
        // (it's the fixed constant 1).
        let mut bv = extract_bivar(
            &solution,
            |i, j| {
                if i == 0 && j == 0 {
                    // Dummy index — will be overwritten below.
                    0
                } else {
                    denom_col(i, j)
                }
            },
            opts.denom_idx_deg,
            opts.denom_var_deg,
        );
        bv.coeffs[0][0] = BigRational::one();
        if bv.is_one() {
            None
        } else {
            Some(bv)
        }
    } else {
        None
    };

    let inhomogeneous = if !opts.homogeneous {
        let bv = extract_bivar(
            &solution,
            inhomo_col,
            opts.inhomo_idx_deg,
            opts.inhomo_var_deg,
        );
        if bv.is_zero() {
            None
        } else {
            Some(bv)
        }
    } else {
        None
    };

    let recurrence = Recurrence {
        terms,
        denominator,
        inhomogeneous,
    };

    if exact_row_indices.is_some()
        && !recurrence_holds_rational_with_derivs(polys, derivs, &recurrence, opts.rec_len)
    {
        return None;
    }

    Some(recurrence)
}

/// Search for a polynomial recurrence satisfied by an integer-coefficient
/// polynomial sequence.
///
/// This is a convenience wrapper around [`find_polynomial_recurrence_rational`].
pub fn find_polynomial_recurrence(
    polys: &[Vec<i64>],
    opts: &RecurrenceOptions,
) -> Option<Recurrence> {
    let rational_polys = i64_polys_to_rational(polys);
    find_polynomial_recurrence_rational(&rational_polys, opts)
}

// ---------------------------------------------------------------------------
// Helpers for building the result
// ---------------------------------------------------------------------------

/// Extract a BivarPoly from a solution vector using a column-index mapping.
fn extract_bivar(
    solution: &[BigRational],
    col_fn: impl Fn(usize, usize) -> usize,
    max_i: usize,
    max_j: usize,
) -> BivarPoly {
    let coeffs: Vec<Vec<BigRational>> = (0..=max_i)
        .map(|i| {
            (0..=max_j)
                .map(|j| solution[col_fn(i, j)].clone())
                .collect()
        })
        .collect();
    // Trim trailing zero rows/cols if desired (we keep them for now).
    BivarPoly { coeffs }
}

// ---------------------------------------------------------------------------
// BivarPoly helpers
// ---------------------------------------------------------------------------

impl BivarPoly {
    /// True when all coefficients are zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs
            .iter()
            .all(|row| row.iter().all(|c| c.is_zero()))
    }

    /// True when the polynomial equals the constant 1.
    pub fn is_one(&self) -> bool {
        for (i, row) in self.coeffs.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                if i == 0 && j == 0 {
                    if !c.is_one() {
                        return false;
                    }
                } else if !c.is_zero() {
                    return false;
                }
            }
        }
        true
    }

    /// True when the polynomial is a single monomial.
    fn is_monomial(&self) -> bool {
        let count = self
            .coeffs
            .iter()
            .flat_map(|row| row.iter())
            .filter(|c| !c.is_zero())
            .count();
        count <= 1
    }

    /// Number of non-zero terms.
    fn num_terms(&self) -> usize {
        self.coeffs
            .iter()
            .flat_map(|row| row.iter())
            .filter(|c| !c.is_zero())
            .count()
    }
}

// ---------------------------------------------------------------------------
// Display
// ---------------------------------------------------------------------------

fn fmt_rational(r: &BigRational) -> String {
    if r.denom() == &BigInt::one() {
        format!("{}", r.numer())
    } else {
        format!("{}/{}", r.numer(), r.denom())
    }
}

/// Format an exact rational coefficient as an integer or `num/den` string.
pub fn format_rational_coeff(r: &BigRational) -> String {
    fmt_rational(r)
}

/// Parse an exact rational coefficient from an integer or `num/den` string.
pub fn parse_rational_coeff(input: &str) -> Result<BigRational, String> {
    let token = input.trim();
    if token.is_empty() {
        return Err("empty rational coefficient".to_string());
    }
    if let Some((num, den)) = token.split_once('/') {
        let numerator = num
            .trim()
            .parse::<BigInt>()
            .map_err(|e| format!("invalid numerator `{}`: {e}", num.trim()))?;
        let denominator = den
            .trim()
            .parse::<BigInt>()
            .map_err(|e| format!("invalid denominator `{}`: {e}", den.trim()))?;
        if denominator.is_zero() {
            return Err(format!("zero denominator in coefficient `{token}`"));
        }
        Ok(BigRational::new(numerator, denominator))
    } else {
        let integer = token
            .parse::<BigInt>()
            .map_err(|e| format!("invalid integer `{token}`: {e}"))?;
        Ok(BigRational::from_integer(integer))
    }
}

fn fmt_monomial(c: &BigRational, n_pow: usize, t_pow: usize) -> String {
    let var = match (n_pow, t_pow) {
        (0, 0) => return fmt_rational(c),
        (1, 0) => "n".into(),
        (i, 0) => format!("n^{i}"),
        (0, 1) => "t".into(),
        (0, j) => format!("t^{j}"),
        (1, 1) => "nt".into(),
        (i, 1) => format!("n^{i}t"),
        (1, j) => format!("nt^{j}"),
        (i, j) => format!("n^{i}t^{j}"),
    };

    if c.is_one() {
        var
    } else if *c == BigRational::from(BigInt::from(-1)) {
        format!("-{var}")
    } else {
        format!("{}{var}", fmt_rational(c))
    }
}

impl fmt::Display for BivarPoly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut terms: Vec<String> = Vec::new();
        // Iterate j (t-power) first, then i (n-power), for natural ordering.
        let max_j = self.coeffs.first().map_or(0, |r| r.len().saturating_sub(1));
        let max_i = self.coeffs.len().saturating_sub(1);
        for j in 0..=max_j {
            for i in 0..=max_i {
                let c = &self.coeffs[i][j];
                if c.is_zero() {
                    continue;
                }
                terms.push(fmt_monomial(c, i, j));
            }
        }

        if terms.is_empty() {
            return write!(f, "0");
        }

        let mut result = terms[0].clone();
        for t in &terms[1..] {
            if let Some(rest) = t.strip_prefix('-') {
                result.push_str(" - ");
                result.push_str(rest);
            } else {
                result.push_str(" + ");
                result.push_str(t);
            }
        }
        write!(f, "{result}")
    }
}

fn fmt_rational_latex(r: &BigRational) -> String {
    if r.denom() == &BigInt::one() {
        format!("{}", r.numer())
    } else {
        format!("\\frac{{{}}}{{{}}}", r.numer(), r.denom())
    }
}

fn fmt_monomial_latex(c: &BigRational, n_pow: usize, t_pow: usize) -> String {
    let var = match (n_pow, t_pow) {
        (0, 0) => return fmt_rational_latex(c),
        (1, 0) => "n".into(),
        (i, 0) => format!("n^{{{i}}}"),
        (0, 1) => "t".into(),
        (0, j) => format!("t^{{{j}}}"),
        (1, 1) => "nt".into(),
        (i, 1) => format!("n^{{{i}}}t"),
        (1, j) => format!("nt^{{{j}}}"),
        (i, j) => format!("n^{{{i}}}t^{{{j}}}"),
    };

    if c.is_one() {
        var
    } else if *c == BigRational::from(BigInt::from(-1)) {
        format!("-{var}")
    } else {
        format!("{}{var}", fmt_rational_latex(c))
    }
}

impl BivarPoly {
    /// Format as a LaTeX expression.
    pub fn to_latex(&self) -> String {
        let mut terms: Vec<String> = Vec::new();
        let max_j = self.coeffs.first().map_or(0, |r| r.len().saturating_sub(1));
        let max_i = self.coeffs.len().saturating_sub(1);
        for j in 0..=max_j {
            for i in 0..=max_i {
                let c = &self.coeffs[i][j];
                if c.is_zero() {
                    continue;
                }
                terms.push(fmt_monomial_latex(c, i, j));
            }
        }
        if terms.is_empty() {
            return "0".into();
        }
        let mut result = terms[0].clone();
        for t in &terms[1..] {
            if let Some(rest) = t.strip_prefix('-') {
                result.push_str(" - ");
                result.push_str(rest);
            } else {
                result.push_str(" + ");
                result.push_str(t);
            }
        }
        result
    }
}

fn fmt_poly_ref_latex(offset: usize, deriv_order: usize) -> String {
    let sub = if offset == 0 {
        "n".into()
    } else {
        format!("n-{offset}")
    };
    match deriv_order {
        0 => format!("P({sub})"),
        1 => format!("P'({sub})"),
        d => format!("P^{{({d})}}({sub})"),
    }
}

fn fmt_poly_ref_latex_with_sign(sign: RecurrenceSign, offset: usize, deriv_order: usize) -> String {
    let pref = fmt_poly_ref_latex(offset, deriv_order);
    match sign {
        RecurrenceSign::None => pref,
        RecurrenceSign::AlternatingN => format!("(-1)^n {pref}"),
    }
}

#[derive(Copy, Clone)]
enum CodeStyle {
    Mathematica,
    Sage,
}

fn recurrence_sign_code(style: CodeStyle, sign: RecurrenceSign) -> Option<&'static str> {
    match (style, sign) {
        (_, RecurrenceSign::None) => None,
        (CodeStyle::Mathematica, RecurrenceSign::AlternatingN) => Some("(-1)^n"),
        (CodeStyle::Sage, RecurrenceSign::AlternatingN) => Some("(-1)**n"),
    }
}

fn fmt_rational_code(style: CodeStyle, r: &BigRational) -> String {
    match style {
        CodeStyle::Mathematica => {
            if r.denom() == &BigInt::one() {
                format!("{}", r.numer())
            } else {
                format!("{}/{}", r.numer(), r.denom())
            }
        }
        CodeStyle::Sage => {
            if r.denom() == &BigInt::one() {
                format!("{}", r.numer())
            } else if r.numer() < &BigInt::zero() {
                format!("-QQ({})/{}", -r.numer().clone(), r.denom())
            } else {
                format!("QQ({})/{}", r.numer(), r.denom())
            }
        }
    }
}

fn fmt_power_code(style: CodeStyle, var: &str, pow: usize) -> Option<String> {
    match pow {
        0 => None,
        1 => Some(var.to_string()),
        _ => Some(match style {
            CodeStyle::Mathematica => format!("{var}^{pow}"),
            CodeStyle::Sage => format!("{var}**{pow}"),
        }),
    }
}

fn fmt_monomial_code_abs(
    style: CodeStyle,
    coeff_abs: &BigRational,
    n_pow: usize,
    t_pow: usize,
) -> String {
    let has_vars = n_pow > 0 || t_pow > 0;
    let mut factors = Vec::new();

    if !(coeff_abs.is_one() && has_vars) {
        let coeff_text = fmt_rational_code(style, coeff_abs);
        let needs_wrap = has_vars && coeff_text.contains('/');
        factors.push(if needs_wrap {
            format!("({coeff_text})")
        } else {
            coeff_text
        });
    }

    if let Some(n_term) = fmt_power_code(style, "n", n_pow) {
        factors.push(n_term);
    }
    if let Some(t_term) = fmt_power_code(style, "t", t_pow) {
        factors.push(t_term);
    }

    if factors.is_empty() {
        "1".to_string()
    } else {
        factors.join("*")
    }
}

fn join_signed_terms(terms: &[(bool, String)]) -> String {
    if terms.is_empty() {
        return "0".to_string();
    }

    let mut result = String::new();
    for (idx, (negative, body)) in terms.iter().enumerate() {
        if idx == 0 {
            if *negative {
                result.push('-');
            }
            result.push_str(body);
        } else if *negative {
            result.push_str(" - ");
            result.push_str(body);
        } else {
            result.push_str(" + ");
            result.push_str(body);
        }
    }
    result
}

fn wrap_if_sum(expr: &str) -> String {
    if expr.contains(" + ") || expr.contains(" - ") {
        format!("({expr})")
    } else {
        expr.to_string()
    }
}

fn fmt_univariate_poly_rational_code(style: CodeStyle, coeffs: &[BigRational]) -> String {
    let mut terms = Vec::new();
    for (pow, coeff) in coeffs.iter().enumerate() {
        if coeff.is_zero() {
            continue;
        }
        let negative = coeff < &BigRational::zero();
        let abs_coeff = if negative {
            -coeff.clone()
        } else {
            coeff.clone()
        };
        let body = fmt_monomial_code_abs(style, &abs_coeff, 0, pow);
        terms.push((negative, body));
    }
    join_signed_terms(&terms)
}

fn fmt_python_rational(r: &BigRational) -> String {
    if r.denom() == &BigInt::one() {
        format!("{}", r.numer())
    } else {
        format!("Fraction({}, {})", r.numer(), r.denom())
    }
}

fn fmt_python_poly_literal(coeffs: &[BigRational]) -> String {
    format!(
        "[{}]",
        coeffs
            .iter()
            .map(fmt_python_rational)
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn fmt_python_bivar_literal(poly: &BivarPoly) -> String {
    format!(
        "[{}]",
        poly.coeffs
            .iter()
            .map(|row| fmt_python_poly_literal(row))
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn python_sign_name(sign: RecurrenceSign) -> &'static str {
    match sign {
        RecurrenceSign::None => "none",
        RecurrenceSign::AlternatingN => "alternating_n",
    }
}

fn python_term_literal(term: &RecurrenceTerm) -> String {
    format!(
        "{{\"offset\": {}, \"deriv_order\": {}, \"sign\": \"{}\", \"coeff\": {}}}",
        term.offset,
        term.deriv_order,
        python_sign_name(term.sign),
        fmt_python_bivar_literal(&term.coeff)
    )
}

impl Recurrence {
    /// Largest recurrence offset used by any term.
    pub fn max_offset(&self) -> usize {
        self.terms.iter().map(|term| term.offset).max().unwrap_or(0)
    }

    fn term_to_mathematica_code(&self, term: &RecurrenceTerm) -> String {
        let coeff = term.coeff.to_mathematica_code();
        let idx = if term.offset == 0 {
            "n".to_string()
        } else {
            format!("n - {}", term.offset)
        };
        let pref = match term.deriv_order {
            0 => format!("P[{idx}, t]"),
            1 => format!("D[P[{idx}, t], t]"),
            d => format!("D[P[{idx}, t], {{t, {d}}}]"),
        };
        let sign = recurrence_sign_code(CodeStyle::Mathematica, term.sign);

        if term.coeff.is_one() {
            if let Some(sign) = sign {
                format!("{sign}*{pref}")
            } else {
                pref
            }
        } else if coeff == "-1" {
            if let Some(sign) = sign {
                format!("-{sign}*{pref}")
            } else {
                format!("-{pref}")
            }
        } else if let Some(sign) = sign {
            format!("{sign}*{}*{pref}", wrap_if_sum(&coeff))
        } else {
            format!("{}*{pref}", wrap_if_sum(&coeff))
        }
    }

    fn term_to_sage_code(&self, term: &RecurrenceTerm) -> String {
        let coeff = term.coeff.to_sage_code();
        let idx = if term.offset == 0 {
            "n".to_string()
        } else {
            format!("n - {}", term.offset)
        };
        let pref = match term.deriv_order {
            0 => format!("P({idx})"),
            1 => format!("P({idx}).derivative(t)"),
            d => format!("P({idx}).derivative(t, {d})"),
        };
        let sign = recurrence_sign_code(CodeStyle::Sage, term.sign);

        if term.coeff.is_one() {
            if let Some(sign) = sign {
                format!("{sign}*{pref}")
            } else {
                pref
            }
        } else if coeff == "-1" {
            if let Some(sign) = sign {
                format!("-{sign}*{pref}")
            } else {
                format!("-{pref}")
            }
        } else if let Some(sign) = sign {
            format!("{sign}*{}*{pref}", wrap_if_sum(&coeff))
        } else {
            format!("{}*{pref}", wrap_if_sum(&coeff))
        }
    }

    fn rhs_to_mathematica_code(&self) -> String {
        let mut pieces: Vec<String> = self
            .terms
            .iter()
            .map(|term| self.term_to_mathematica_code(term))
            .collect();
        if let Some(ref inh) = self.inhomogeneous {
            pieces.push(inh.to_mathematica_code());
        }
        let rhs = if pieces.is_empty() {
            "0".to_string()
        } else {
            pieces.join(" + ")
        };
        if let Some(ref denom) = self.denominator {
            format!(
                "Expand[Together[({})/({})]]",
                rhs,
                denom.to_mathematica_code()
            )
        } else {
            format!("Expand[{rhs}]")
        }
    }

    fn rhs_to_sage_code(&self) -> String {
        let mut pieces: Vec<String> = self
            .terms
            .iter()
            .map(|term| self.term_to_sage_code(term))
            .collect();
        if let Some(ref inh) = self.inhomogeneous {
            pieces.push(inh.to_sage_code());
        }
        let rhs = if pieces.is_empty() {
            "0".to_string()
        } else {
            pieces.join(" + ")
        };
        if let Some(ref denom) = self.denominator {
            format!("R(K({rhs}) / K({}))", denom.to_sage_code())
        } else {
            format!("R({rhs})")
        }
    }

    pub fn to_mathematica_definition_rational(&self, initial_polys: &[Vec<BigRational>]) -> String {
        let base_count = self.max_offset().min(initial_polys.len());
        let start_n = base_count + 1;
        let mut lines = vec!["ClearAll[P];".to_string()];
        for (idx, coeffs) in initial_polys.iter().take(base_count).enumerate() {
            lines.push(format!(
                "P[{}, t_] := {};",
                idx + 1,
                fmt_univariate_poly_rational_code(CodeStyle::Mathematica, coeffs)
            ));
        }
        lines.push(format!(
            "P[n_Integer /; n >= {start_n}, t_] := P[n, t] = {};",
            self.rhs_to_mathematica_code()
        ));
        lines.push(String::new());
        lines.push("(* Example: Table[P[n, t], {n, 1, 10}] *)".to_string());
        lines.join("\n")
    }

    pub fn to_mathematica_definition(&self, initial_polys: &[Vec<i64>]) -> String {
        let rational_polys = i64_polys_to_rational(initial_polys);
        self.to_mathematica_definition_rational(&rational_polys)
    }

    pub fn to_sage_definition_rational(&self, initial_polys: &[Vec<BigRational>]) -> String {
        let base_count = self.max_offset().min(initial_polys.len());
        let mut lines = vec![
            "R.<t> = PolynomialRing(QQ)".to_string(),
            "K = R.fraction_field()".to_string(),
            "_P_cache = {".to_string(),
        ];
        for (idx, coeffs) in initial_polys.iter().take(base_count).enumerate() {
            let coeff_list = coeffs
                .iter()
                .map(|coeff| fmt_rational_code(CodeStyle::Sage, coeff))
                .collect::<Vec<_>>()
                .join(", ");
            lines.push(format!("    {}: R([{}]),", idx + 1, coeff_list));
        }
        lines.push("}".to_string());
        lines.push(String::new());
        lines.push("def P(n):".to_string());
        lines.push("    if n < 1:".to_string());
        lines.push("        raise ValueError(\"n must be a positive integer\")".to_string());
        lines.push("    if n in _P_cache:".to_string());
        lines.push("        return _P_cache[n]".to_string());
        lines.push(format!("    value = {}", self.rhs_to_sage_code()));
        lines.push("    _P_cache[n] = value".to_string());
        lines.push("    return value".to_string());
        lines.push(String::new());
        lines.push("# Example: [P(n) for n in range(1, 11)]".to_string());
        lines.join("\n")
    }

    pub fn to_sage_definition(&self, initial_polys: &[Vec<i64>]) -> String {
        let rational_polys = i64_polys_to_rational(initial_polys);
        self.to_sage_definition_rational(&rational_polys)
    }

    pub fn to_python_definition_rational(&self, initial_polys: &[Vec<BigRational>]) -> String {
        let base_count = self.max_offset().min(initial_polys.len());
        let mut lines = vec![
            "from fractions import Fraction".to_string(),
            String::new(),
            "def _trim(p):".to_string(),
            "    p = list(p)".to_string(),
            "    while len(p) > 1 and p[-1] == 0:".to_string(),
            "        p.pop()".to_string(),
            "    return p or [0]".to_string(),
            String::new(),
            "def _is_zero(p):".to_string(),
            "    return all(coeff == 0 for coeff in p)".to_string(),
            String::new(),
            "def _add(a, b):".to_string(),
            "    length = max(len(a), len(b))".to_string(),
            "    out = [0] * length".to_string(),
            "    for i in range(length):".to_string(),
            "        out[i] = (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)"
                .to_string(),
            "    return _trim(out)".to_string(),
            String::new(),
            "def _scale(p, c):".to_string(),
            "    return _trim([c * coeff for coeff in p])".to_string(),
            String::new(),
            "def _mul(a, b):".to_string(),
            "    if _is_zero(a) or _is_zero(b):".to_string(),
            "        return [0]".to_string(),
            "    out = [0] * (len(a) + len(b) - 1)".to_string(),
            "    for i, ai in enumerate(a):".to_string(),
            "        for j, bj in enumerate(b):".to_string(),
            "            out[i + j] += ai * bj".to_string(),
            "    return _trim(out)".to_string(),
            String::new(),
            "def _derivative(p, order):".to_string(),
            "    out = list(p)".to_string(),
            "    for _ in range(order):".to_string(),
            "        if len(out) <= 1:".to_string(),
            "            return [0]".to_string(),
            "        out = [(i + 1) * out[i + 1] for i in range(len(out) - 1)]".to_string(),
            "    return _trim(out)".to_string(),
            String::new(),
            "def _div_exact(numerator, denominator):".to_string(),
            "    denominator = _trim(denominator)".to_string(),
            "    if _is_zero(denominator):".to_string(),
            "        raise ZeroDivisionError(\"zero recurrence denominator\")".to_string(),
            "    remainder = _trim(numerator)".to_string(),
            "    if _is_zero(remainder):".to_string(),
            "        return [0]".to_string(),
            "    if len(remainder) < len(denominator):".to_string(),
            "        raise ValueError(\"non-exact polynomial quotient\")".to_string(),
            "    quotient = [0] * (len(remainder) - len(denominator) + 1)".to_string(),
            "    denominator_degree = len(denominator) - 1".to_string(),
            "    denominator_lc = denominator[-1]".to_string(),
            "    while not _is_zero(remainder) and len(remainder) - 1 >= denominator_degree:"
                .to_string(),
            "        shift = len(remainder) - 1 - denominator_degree".to_string(),
            "        factor = remainder[-1] / denominator_lc".to_string(),
            "        quotient[shift] += factor".to_string(),
            "        for i, coeff in enumerate(denominator):".to_string(),
            "            remainder[i + shift] -= factor * coeff".to_string(),
            "        remainder = _trim(remainder)".to_string(),
            "    if not _is_zero(remainder):".to_string(),
            "        raise ValueError(\"non-exact polynomial quotient\")".to_string(),
            "    return _trim(quotient)".to_string(),
            String::new(),
            "def _bivar(coeffs, n):".to_string(),
            "    width = max((len(row) for row in coeffs), default=0)".to_string(),
            "    out = [0] * width".to_string(),
            "    n_power = 1".to_string(),
            "    for row in coeffs:".to_string(),
            "        for t_power, coeff in enumerate(row):".to_string(),
            "            out[t_power] += coeff * n_power".to_string(),
            "        n_power *= n".to_string(),
            "    return _trim(out)".to_string(),
            String::new(),
            "_P_cache = {".to_string(),
        ];
        for (idx, coeffs) in initial_polys.iter().take(base_count).enumerate() {
            lines.push(format!(
                "    {}: {},",
                idx + 1,
                fmt_python_poly_literal(coeffs)
            ));
        }
        lines.push("}".to_string());
        lines.push(format!(
            "_terms = [{}]",
            self.terms
                .iter()
                .map(python_term_literal)
                .collect::<Vec<_>>()
                .join(", ")
        ));
        lines.push(format!(
            "_denominator = {}",
            self.denominator
                .as_ref()
                .map(fmt_python_bivar_literal)
                .unwrap_or_else(|| "None".to_string())
        ));
        lines.push(format!(
            "_inhomogeneous = {}",
            self.inhomogeneous
                .as_ref()
                .map(fmt_python_bivar_literal)
                .unwrap_or_else(|| "None".to_string())
        ));
        lines.push(String::new());
        lines.push("def P(n):".to_string());
        lines.push("    if n < 1:".to_string());
        lines.push("        raise ValueError(\"n must be a positive integer\")".to_string());
        lines.push("    if n in _P_cache:".to_string());
        lines.push("        return _P_cache[n]".to_string());
        lines.push("    start = max(_P_cache) + 1 if _P_cache else 1".to_string());
        lines.push("    for k in range(start, n + 1):".to_string());
        lines.push("        rhs = [0]".to_string());
        lines.push("        for term in _terms:".to_string());
        lines.push(
            "            ref = _derivative(P(k - term[\"offset\"]), term[\"deriv_order\"])"
                .to_string(),
        );
        lines.push("            coeff = _bivar(term[\"coeff\"], k)".to_string());
        lines.push("            product = _mul(coeff, ref)".to_string());
        lines
            .push("            if term[\"sign\"] == \"alternating_n\" and k % 2 == 1:".to_string());
        lines.push("                product = _scale(product, -1)".to_string());
        lines.push("            rhs = _add(rhs, product)".to_string());
        lines.push("        if _inhomogeneous is not None:".to_string());
        lines.push("            rhs = _add(rhs, _bivar(_inhomogeneous, k))".to_string());
        lines.push("        if _denominator is not None:".to_string());
        lines.push("            rhs = _div_exact(rhs, _bivar(_denominator, k))".to_string());
        lines.push("        _P_cache[k] = rhs".to_string());
        lines.push("    return _P_cache[n]".to_string());
        lines.push(String::new());
        lines.push("# Example: [P(n) for n in range(1, 11)]".to_string());
        lines.join("\n")
    }

    pub fn to_python_definition(&self, initial_polys: &[Vec<i64>]) -> String {
        let rational_polys = i64_polys_to_rational(initial_polys);
        self.to_python_definition_rational(&rational_polys)
    }

    /// Format the recurrence as a LaTeX expression.
    pub fn to_latex(&self) -> String {
        let mut s = String::new();
        if let Some(ref denom) = self.denominator {
            let dl = denom.to_latex();
            if denom.num_terms() > 1 {
                s.push_str(&format!("\\bigl({dl}\\bigr) "));
            } else {
                s.push_str(&format!("{dl} \\cdot "));
            }
        }
        s.push_str("P(n) = ");

        let mut first = true;
        for term in &self.terms {
            let cl = term.coeff.to_latex();
            let pref = fmt_poly_ref_latex_with_sign(term.sign, term.offset, term.deriv_order);

            if first {
                first = false;
                if cl == "1" {
                    s.push_str(&pref);
                } else if cl == "-1" {
                    s.push_str(&format!("-{pref}"));
                } else if term.coeff.num_terms() > 1 {
                    s.push_str(&format!("\\bigl({cl}\\bigr) {pref}"));
                } else {
                    s.push_str(&format!("{cl} \\cdot {pref}"));
                }
            } else if cl.starts_with('-') && term.coeff.is_monomial() {
                let pos = &cl[1..];
                if pos == "1" {
                    s.push_str(&format!(" - {pref}"));
                } else {
                    s.push_str(&format!(" - {pos} \\cdot {pref}"));
                }
            } else if term.coeff.num_terms() > 1 {
                s.push_str(&format!(" + \\bigl({cl}\\bigr) {pref}"));
            } else if cl == "1" {
                s.push_str(&format!(" + {pref}"));
            } else {
                s.push_str(&format!(" + {cl} \\cdot {pref}"));
            }
        }

        if let Some(ref inh) = self.inhomogeneous {
            if !first {
                s.push_str(" + ");
            }
            s.push_str(&inh.to_latex());
        }

        s
    }
}

impl BivarPoly {
    fn to_code(&self, style: CodeStyle) -> String {
        let mut terms: Vec<(bool, String)> = Vec::new();
        let max_j = self.coeffs.first().map_or(0, |r| r.len().saturating_sub(1));
        let max_i = self.coeffs.len().saturating_sub(1);
        for j in 0..=max_j {
            for i in 0..=max_i {
                let c = &self.coeffs[i][j];
                if c.is_zero() {
                    continue;
                }
                let negative = c < &BigRational::zero();
                let coeff_abs = if negative { -c.clone() } else { c.clone() };
                terms.push((negative, fmt_monomial_code_abs(style, &coeff_abs, i, j)));
            }
        }
        join_signed_terms(&terms)
    }

    pub fn to_mathematica_code(&self) -> String {
        self.to_code(CodeStyle::Mathematica)
    }

    pub fn to_sage_code(&self) -> String {
        self.to_code(CodeStyle::Sage)
    }
}

fn fmt_poly_ref(offset: usize, deriv_order: usize) -> String {
    let sub = if offset == 0 {
        "n".into()
    } else {
        format!("n-{offset}")
    };
    match deriv_order {
        0 => format!("P({sub})"),
        1 => format!("P'({sub})"),
        d => format!("P^({d})({sub})"),
    }
}

fn fmt_poly_ref_with_sign(sign: RecurrenceSign, offset: usize, deriv_order: usize) -> String {
    let pref = fmt_poly_ref(offset, deriv_order);
    match sign {
        RecurrenceSign::None => pref,
        RecurrenceSign::AlternatingN => format!("(-1)^n {pref}"),
    }
}

impl fmt::Display for Recurrence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // LHS
        if let Some(ref denom) = self.denominator {
            write!(f, "({denom}) ")?;
        }
        write!(f, "P(n) = ")?;

        let mut first = true;
        for term in &self.terms {
            let cs = format!("{}", term.coeff);
            let pref = fmt_poly_ref_with_sign(term.sign, term.offset, term.deriv_order);

            if first {
                first = false;
                if cs == "1" {
                    write!(f, "{pref}")?;
                } else if cs == "-1" {
                    write!(f, "-{pref}")?;
                } else if term.coeff.num_terms() > 1 {
                    write!(f, "({cs}) {pref}")?;
                } else {
                    write!(f, "{cs} {pref}")?;
                }
            } else {
                // Determine sign for nice formatting.
                if cs.starts_with('-') && term.coeff.is_monomial() {
                    // Single negative term: pull out the minus sign.
                    let pos = &cs[1..];
                    if pos == "1" {
                        write!(f, " - {pref}")?;
                    } else {
                        write!(f, " - {pos} {pref}")?;
                    }
                } else if term.coeff.num_terms() > 1 {
                    write!(f, " + ({cs}) {pref}")?;
                } else if cs == "1" {
                    write!(f, " + {pref}")?;
                } else {
                    write!(f, " + {cs} {pref}")?;
                }
            }
        }

        if let Some(ref inh) = self.inhomogeneous {
            if !first {
                write!(f, " + ")?;
            }
            write!(f, "{inh}")?;
        }

        Ok(())
    }
}

/// Errors that can occur while replaying a recurrence to generate rows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RecurrenceEvaluationError {
    /// More initial rows are needed before the recurrence can be evaluated.
    InsufficientInitialRows { required: usize, provided: usize },
    /// An offset points before the first available initial/generated row.
    OffsetBeforeInitialRow {
        n: usize,
        offset: usize,
        first_index: usize,
    },
    /// The evaluated LHS factor is the zero polynomial.
    ZeroDenominator,
    /// A denominator recurrence produced a rational function rather than a polynomial.
    NonPolynomialQuotient,
}

impl fmt::Display for RecurrenceEvaluationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InsufficientInitialRows { required, provided } => write!(
                f,
                "insufficient initial rows: need at least {required}, got {provided}"
            ),
            Self::OffsetBeforeInitialRow {
                n,
                offset,
                first_index,
            } => write!(
                f,
                "cannot evaluate P({n}): offset {offset} refers before first index {first_index}"
            ),
            Self::ZeroDenominator => write!(f, "evaluated denominator is the zero polynomial"),
            Self::NonPolynomialQuotient => write!(
                f,
                "recurrence evaluation did not divide exactly to a polynomial"
            ),
        }
    }
}

impl std::error::Error for RecurrenceEvaluationError {}

impl Recurrence {
    /// Evaluate the next row after `rows`, where `rows[0]` is `P(first_index)`.
    pub fn evaluate_next_rational(
        &self,
        rows: &[Vec<BigRational>],
        first_index: usize,
    ) -> Result<Vec<BigRational>, RecurrenceEvaluationError> {
        let required = self.max_offset();
        if rows.len() < required {
            return Err(RecurrenceEvaluationError::InsufficientInitialRows {
                required,
                provided: rows.len(),
            });
        }

        let n = first_index + rows.len();
        let mut rhs = rational_zero_poly();

        for term in &self.terms {
            let Some(reference_index) = n.checked_sub(term.offset) else {
                return Err(RecurrenceEvaluationError::OffsetBeforeInitialRow {
                    n,
                    offset: term.offset,
                    first_index,
                });
            };
            if reference_index < first_index {
                return Err(RecurrenceEvaluationError::OffsetBeforeInitialRow {
                    n,
                    offset: term.offset,
                    first_index,
                });
            }
            let row_index = reference_index - first_index;
            let Some(ref_poly) = rows.get(row_index) else {
                return Err(RecurrenceEvaluationError::OffsetBeforeInitialRow {
                    n,
                    offset: term.offset,
                    first_index,
                });
            };

            let deriv = poly_nth_derivative_rational(ref_poly, term.deriv_order);
            let coeff = bivar_eval_n(&term.coeff, n);
            let product = poly_mul_rational(&coeff, &deriv);
            let sign = match term.sign {
                RecurrenceSign::None => BigRational::one(),
                RecurrenceSign::AlternatingN if n % 2 == 1 => {
                    BigRational::from_integer(BigInt::from(-1))
                }
                RecurrenceSign::AlternatingN => BigRational::one(),
            };
            poly_add_scaled_assign(&mut rhs, &product, &sign);
        }

        if let Some(inhomogeneous) = &self.inhomogeneous {
            let inh = bivar_eval_n(inhomogeneous, n);
            poly_add_scaled_assign(&mut rhs, &inh, &BigRational::one());
        }

        if let Some(denominator) = &self.denominator {
            poly_div_exact_rational(&rhs, &bivar_eval_n(denominator, n))
        } else {
            Ok(trim_poly_rational(rhs))
        }
    }

    /// Generate `total_rows` rows, including the supplied initial rows.
    pub fn generate_rows_rational(
        &self,
        initial_rows: &[Vec<BigRational>],
        first_index: usize,
        total_rows: usize,
    ) -> Result<Vec<Vec<BigRational>>, RecurrenceEvaluationError> {
        if total_rows <= initial_rows.len() {
            return Ok(initial_rows[..total_rows].to_vec());
        }
        let required = self.max_offset();
        if initial_rows.len() < required {
            return Err(RecurrenceEvaluationError::InsufficientInitialRows {
                required,
                provided: initial_rows.len(),
            });
        }

        let mut rows = initial_rows.to_vec();
        while rows.len() < total_rows {
            let next = self.evaluate_next_rational(&rows, first_index)?;
            rows.push(next);
        }
        Ok(rows)
    }
}

// ---------------------------------------------------------------------------
// JSON schema
// ---------------------------------------------------------------------------

/// Current JSON schema tag for serialized recurrences.
pub const RECURRENCE_JSON_SCHEMA: &str = "polynomial-tools.recurrence.v1";

/// Stable JSON representation of a recurrence together with initial rows.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceJson {
    /// Schema identifier. Currently `polynomial-tools.recurrence.v1`.
    pub schema: String,
    /// Index of the first polynomial in `initial_polynomials`.
    pub first_index: usize,
    /// Initial coefficient rows, in ascending powers of `t`.
    pub initial_polynomials: Vec<Vec<String>>,
    /// The recurrence data.
    pub recurrence: RecurrenceJsonData,
    /// Optional metadata from an adaptive search.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub search: Option<RecurrenceJsonSearch>,
}

/// Stable JSON representation of a recurrence without initial rows.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceJsonData {
    pub terms: Vec<RecurrenceTermJson>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub denominator: Option<BivarPolyJson>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub inhomogeneous: Option<BivarPolyJson>,
}

/// Stable JSON representation of one recurrence term.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceTermJson {
    pub offset: usize,
    pub deriv_order: usize,
    pub sign: String,
    pub coeff: BivarPolyJson,
}

/// Stable JSON representation of a bivariate polynomial in `(n,t)`.
///
/// `coeffs[i][j]` is the coefficient of `n^i t^j`, stored as an exact rational
/// string.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct BivarPolyJson {
    pub coeffs: Vec<Vec<String>>,
}

/// Optional metadata describing how an adaptive recurrence was found.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceJsonSearch {
    pub recurrence_text: String,
    pub source_rows: usize,
    pub skip_prefix: usize,
    pub unknowns: usize,
    pub weighted_unknowns: usize,
    pub equations: usize,
    pub fit_polynomials: usize,
    pub verification_polynomials: usize,
    pub candidates_tried: usize,
    pub options: RecurrenceOptionsJson,
}

/// Search-space options for the recurrence that was found.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceOptionsJson {
    pub var_deg: usize,
    pub idx_deg: usize,
    pub diff_deg: usize,
    pub rec_len: usize,
    pub homogeneous: bool,
    pub inhomo_var_deg: usize,
    pub inhomo_idx_deg: usize,
    pub denom_var_deg: usize,
    pub denom_idx_deg: usize,
    pub alternating_sign: bool,
}

impl From<&RecurrenceOptions> for RecurrenceOptionsJson {
    fn from(opts: &RecurrenceOptions) -> Self {
        Self {
            var_deg: opts.var_deg,
            idx_deg: opts.idx_deg,
            diff_deg: opts.diff_deg,
            rec_len: opts.rec_len,
            homogeneous: opts.homogeneous,
            inhomo_var_deg: opts.inhomo_var_deg,
            inhomo_idx_deg: opts.inhomo_idx_deg,
            denom_var_deg: opts.denom_var_deg,
            denom_idx_deg: opts.denom_idx_deg,
            alternating_sign: opts.alternating_sign,
        }
    }
}

impl RecurrenceJson {
    /// Create a JSON record from a recurrence and initial rows.
    pub fn from_recurrence_rational(
        recurrence: &Recurrence,
        first_index: usize,
        initial_polynomials: &[Vec<BigRational>],
        search: Option<RecurrenceJsonSearch>,
    ) -> Self {
        Self {
            schema: RECURRENCE_JSON_SCHEMA.to_string(),
            first_index,
            initial_polynomials: initial_polynomials
                .iter()
                .map(|row| row.iter().map(format_rational_coeff).collect())
                .collect(),
            recurrence: RecurrenceJsonData::from_recurrence(recurrence),
            search,
        }
    }

    /// Convert this JSON record back to a recurrence and initial rows.
    pub fn to_recurrence_parts(
        &self,
    ) -> Result<(Recurrence, usize, Vec<Vec<BigRational>>), String> {
        if self.schema != RECURRENCE_JSON_SCHEMA {
            return Err(format!(
                "unsupported recurrence JSON schema `{}`",
                self.schema
            ));
        }
        let recurrence = self.recurrence.to_recurrence()?;
        let initial_polynomials = self
            .initial_polynomials
            .iter()
            .enumerate()
            .map(|(row_idx, row)| {
                row.iter()
                    .enumerate()
                    .map(|(coeff_idx, coeff)| {
                        parse_rational_coeff(coeff).map_err(|e| {
                            format!(
                                "invalid initial coefficient at row {}, coefficient {}: {e}",
                                row_idx + 1,
                                coeff_idx
                            )
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()
                    .map(trim_poly_rational)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok((recurrence, self.first_index, initial_polynomials))
    }
}

impl RecurrenceJsonData {
    fn from_recurrence(recurrence: &Recurrence) -> Self {
        Self {
            terms: recurrence
                .terms
                .iter()
                .map(RecurrenceTermJson::from_term)
                .collect(),
            denominator: recurrence
                .denominator
                .as_ref()
                .map(BivarPolyJson::from_bivar),
            inhomogeneous: recurrence
                .inhomogeneous
                .as_ref()
                .map(BivarPolyJson::from_bivar),
        }
    }

    fn to_recurrence(&self) -> Result<Recurrence, String> {
        Ok(Recurrence {
            terms: self
                .terms
                .iter()
                .map(RecurrenceTermJson::to_term)
                .collect::<Result<Vec<_>, _>>()?,
            denominator: self
                .denominator
                .as_ref()
                .map(BivarPolyJson::to_bivar)
                .transpose()?,
            inhomogeneous: self
                .inhomogeneous
                .as_ref()
                .map(BivarPolyJson::to_bivar)
                .transpose()?,
        })
    }
}

impl RecurrenceTermJson {
    fn from_term(term: &RecurrenceTerm) -> Self {
        Self {
            offset: term.offset,
            deriv_order: term.deriv_order,
            sign: match term.sign {
                RecurrenceSign::None => "none",
                RecurrenceSign::AlternatingN => "alternating_n",
            }
            .to_string(),
            coeff: BivarPolyJson::from_bivar(&term.coeff),
        }
    }

    fn to_term(&self) -> Result<RecurrenceTerm, String> {
        let sign = match self.sign.as_str() {
            "none" => RecurrenceSign::None,
            "alternating_n" => RecurrenceSign::AlternatingN,
            other => return Err(format!("unsupported recurrence sign `{other}`")),
        };
        Ok(RecurrenceTerm {
            offset: self.offset,
            deriv_order: self.deriv_order,
            sign,
            coeff: self.coeff.to_bivar()?,
        })
    }
}

impl BivarPolyJson {
    fn from_bivar(poly: &BivarPoly) -> Self {
        Self {
            coeffs: poly
                .coeffs
                .iter()
                .map(|row| row.iter().map(format_rational_coeff).collect())
                .collect(),
        }
    }

    fn to_bivar(&self) -> Result<BivarPoly, String> {
        let coeffs = self
            .coeffs
            .iter()
            .enumerate()
            .map(|(i, row)| {
                row.iter()
                    .enumerate()
                    .map(|(j, coeff)| {
                        parse_rational_coeff(coeff).map_err(|e| {
                            format!("invalid bivariate coefficient ({i},{j}) `{coeff}`: {e}")
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(BivarPoly { coeffs })
    }
}

// ---------------------------------------------------------------------------
// Recurrence verification
// ---------------------------------------------------------------------------

/// Check whether a recurrence holds at the 1-based row index `nn`.
///
/// This tests the exact polynomial identity
/// `f(n,t) P_n(t) = RHS(n,t)` without using the linear system that produced
/// the recurrence.
pub fn recurrence_holds_at_rational(
    polys: &[Vec<BigRational>],
    rec: &Recurrence,
    nn: usize,
) -> bool {
    recurrence_holds_at_rational_with_degree(polys, rec, nn, recurrence_n_degree(rec))
}

fn recurrence_holds_at_rational_with_degree(
    polys: &[Vec<BigRational>],
    rec: &Recurrence,
    nn: usize,
    rec_n_degree: usize,
) -> bool {
    if nn == 0 || nn > polys.len() {
        return false;
    }

    let n_powers = rational_index_powers(nn, rec_n_degree);
    let current = &polys[nn - 1];
    let lhs = if let Some(denom) = &rec.denominator {
        poly_mul_rational(&bivar_eval_n_with_powers(denom, &n_powers), current)
    } else {
        trim_poly_rational(current.clone())
    };

    let mut rhs = rational_zero_poly();
    for term in &rec.terms {
        if term.offset == 0 || term.offset >= nn {
            return false;
        }
        let ref_poly = &polys[nn - 1 - term.offset];
        let deriv = poly_nth_derivative_rational(ref_poly, term.deriv_order);
        let coeff = bivar_eval_n_with_powers(&term.coeff, &n_powers);
        let product = poly_mul_rational(&coeff, &deriv);
        let sign = match term.sign {
            RecurrenceSign::None => BigRational::one(),
            RecurrenceSign::AlternatingN if nn % 2 == 1 => {
                BigRational::from_integer(BigInt::from(-1))
            }
            RecurrenceSign::AlternatingN => BigRational::one(),
        };
        poly_add_scaled_assign(&mut rhs, &product, &sign);
    }

    if let Some(inhomogeneous) = &rec.inhomogeneous {
        let inh = bivar_eval_n_with_powers(inhomogeneous, &n_powers);
        poly_add_scaled_assign(&mut rhs, &inh, &BigRational::one());
    }

    poly_equal_rational(&lhs, &rhs)
}

fn recurrence_holds_at_rational_with_derivs_and_degree(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    rec: &Recurrence,
    nn: usize,
    rec_n_degree: usize,
) -> bool {
    if nn == 0 || nn > polys.len() || nn > derivs.len() {
        return false;
    }

    let n_powers = rational_index_powers(nn, rec_n_degree);
    let current = &polys[nn - 1];
    let lhs = if let Some(denom) = &rec.denominator {
        poly_mul_rational(&bivar_eval_n_with_powers(denom, &n_powers), current)
    } else {
        trim_poly_rational(current.clone())
    };

    let mut rhs = rational_zero_poly();
    for term in &rec.terms {
        if term.offset == 0 || term.offset >= nn {
            return false;
        }
        let Some(poly_derivs) = derivs.get(nn - 1 - term.offset) else {
            return false;
        };
        let Some(deriv) = poly_derivs.get(term.deriv_order) else {
            return false;
        };
        let coeff = bivar_eval_n_with_powers(&term.coeff, &n_powers);
        let product = poly_mul_rational(&coeff, deriv);
        let sign = match term.sign {
            RecurrenceSign::None => BigRational::one(),
            RecurrenceSign::AlternatingN if nn % 2 == 1 => {
                BigRational::from_integer(BigInt::from(-1))
            }
            RecurrenceSign::AlternatingN => BigRational::one(),
        };
        poly_add_scaled_assign(&mut rhs, &product, &sign);
    }

    if let Some(inhomogeneous) = &rec.inhomogeneous {
        let inh = bivar_eval_n_with_powers(inhomogeneous, &n_powers);
        poly_add_scaled_assign(&mut rhs, &inh, &BigRational::one());
    }

    poly_equal_rational(&lhs, &rhs)
}

/// Check whether a recurrence holds for all admissible rows from `start_nn`.
///
/// The row index is 1-based. Rows before `rec_len + 1` are skipped because
/// there are not enough previous polynomials for the recurrence.
pub fn recurrence_holds_from_rational(
    polys: &[Vec<BigRational>],
    rec: &Recurrence,
    rec_len: usize,
    start_nn: usize,
) -> bool {
    let first = start_nn.max(rec_len + 1);
    let rec_n_degree = recurrence_n_degree(rec);
    (first..=polys.len())
        .all(|nn| recurrence_holds_at_rational_with_degree(polys, rec, nn, rec_n_degree))
}

fn recurrence_holds_from_rational_with_derivs(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    rec: &Recurrence,
    rec_len: usize,
    start_nn: usize,
) -> bool {
    let first = start_nn.max(rec_len + 1);
    let rec_n_degree = recurrence_n_degree(rec);
    (first..=polys.len()).all(|nn| {
        recurrence_holds_at_rational_with_derivs_and_degree(polys, derivs, rec, nn, rec_n_degree)
    })
}

/// Check whether a recurrence holds on every admissible row.
pub fn recurrence_holds_rational(
    polys: &[Vec<BigRational>],
    rec: &Recurrence,
    rec_len: usize,
) -> bool {
    recurrence_holds_from_rational(polys, rec, rec_len, rec_len + 1)
}

fn recurrence_holds_rational_with_derivs(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    rec: &Recurrence,
    rec_len: usize,
) -> bool {
    recurrence_holds_from_rational_with_derivs(polys, derivs, rec, rec_len, rec_len + 1)
}

// ---------------------------------------------------------------------------
// Parameter counting
// ---------------------------------------------------------------------------

/// Parameter counts for one candidate recurrence search space.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CandidateComplexity {
    /// Actual number of linear unknowns in the system.
    pub raw_unknowns: usize,
    /// Weighted complexity used only for ordering adaptive candidates.
    pub weighted_unknowns: usize,
    /// Non-alternating, non-derivative recurrence coefficient unknowns.
    pub ordinary_unknowns: usize,
    /// Non-alternating derivative recurrence coefficient unknowns.
    pub derivative_unknowns: usize,
    /// Alternating, non-derivative recurrence coefficient unknowns.
    pub alternating_unknowns: usize,
    /// Alternating derivative recurrence coefficient unknowns.
    pub alternating_derivative_unknowns: usize,
    /// LHS denominator unknowns.
    pub denominator_unknowns: usize,
    /// Inhomogeneous-term unknowns.
    pub inhomogeneous_unknowns: usize,
}

/// Count and weight the parameter slots for a candidate recurrence space.
///
/// The raw count is used for the linear-algebra dimension check. The weighted
/// count is only a heuristic search score: derivative parameters count double,
/// denominator parameters count double, alternating parameters are delayed, and
/// inhomogeneous parameters are delayed further.
pub fn candidate_complexity(opts: &RecurrenceOptions) -> CandidateComplexity {
    let denominator_unknowns = (opts.denom_idx_deg + 1) * (opts.denom_var_deg + 1) - 1;
    let vars_per_coeff = (opts.idx_deg + 1) * (opts.var_deg + 1);
    let ordinary_unknowns = opts.rec_len * vars_per_coeff;
    let derivative_unknowns = opts.rec_len * opts.diff_deg * vars_per_coeff;
    let alternating_unknowns = if opts.alternating_sign {
        opts.rec_len * vars_per_coeff
    } else {
        0
    };
    let alternating_derivative_unknowns = if opts.alternating_sign {
        opts.rec_len * opts.diff_deg * vars_per_coeff
    } else {
        0
    };
    let inhomogeneous_unknowns = if opts.homogeneous {
        0
    } else {
        (opts.inhomo_idx_deg + 1) * (opts.inhomo_var_deg + 1)
    };
    let raw_unknowns = ordinary_unknowns
        + derivative_unknowns
        + alternating_unknowns
        + alternating_derivative_unknowns
        + denominator_unknowns
        + inhomogeneous_unknowns;
    let weighted_unknowns = ordinary_unknowns
        + 2 * derivative_unknowns
        + 3 * alternating_unknowns
        + 4 * alternating_derivative_unknowns
        + 2 * denominator_unknowns
        + 4 * inhomogeneous_unknowns;
    CandidateComplexity {
        raw_unknowns,
        weighted_unknowns,
        ordinary_unknowns,
        derivative_unknowns,
        alternating_unknowns,
        alternating_derivative_unknowns,
        denominator_unknowns,
        inhomogeneous_unknowns,
    }
}

/// Count the total number of unknowns for a given set of recurrence options.
pub fn count_unknowns(opts: &RecurrenceOptions) -> usize {
    candidate_complexity(opts).raw_unknowns
}

/// Weighted unknown count used for adaptive candidate ordering.
pub fn count_weighted_unknowns(opts: &RecurrenceOptions) -> usize {
    candidate_complexity(opts).weighted_unknowns
}

/// Count the total number of equations for a given set of recurrence options
/// and rational input polynomials.
pub fn count_equations_rational(polys: &[Vec<BigRational>], opts: &RecurrenceOptions) -> usize {
    let m = polys.len();
    if m <= opts.rec_len {
        return 0;
    }
    let max_poly_deg = polys
        .iter()
        .map(|p| poly_degree_rational(p))
        .max()
        .unwrap_or(0);
    let max_j = opts
        .var_deg
        .max(opts.denom_var_deg)
        .max(if opts.homogeneous {
            0
        } else {
            opts.inhomo_var_deg
        });
    let eqs_per_nn = max_j + max_poly_deg + 1;
    (m - opts.rec_len) * eqs_per_nn
}

/// Count the total number of equations for integer input polynomials.
pub fn count_equations(polys: &[Vec<i64>], opts: &RecurrenceOptions) -> usize {
    let rational_polys = i64_polys_to_rational(polys);
    count_equations_rational(&rational_polys, opts)
}

// ---------------------------------------------------------------------------
// Adaptive search
// ---------------------------------------------------------------------------

/// Options controlling the adaptive search bounds.
#[derive(Debug, Clone)]
pub struct AdaptiveSearchOptions {
    /// Number of initial polynomials to ignore before searching.
    pub skip_prefix: usize,
    /// Minimum recurrence length to try.
    pub min_rec_len: usize,
    /// Maximum recurrence length to try.
    pub max_rec_len: usize,
    /// Minimum degree in t for coefficients.
    pub min_var_deg: usize,
    /// Maximum degree in t for coefficients.
    pub max_var_deg: usize,
    /// Minimum degree in n for coefficients.
    pub min_idx_deg: usize,
    /// Maximum degree in n for coefficients.
    pub max_idx_deg: usize,
    /// Minimum derivative order.
    pub min_diff_deg: usize,
    /// Maximum derivative order.
    pub max_diff_deg: usize,
    /// Also search inhomogeneous recurrences.
    pub try_inhomogeneous: bool,
    /// Minimum degree in t for the inhomogeneous term when enabled.
    pub min_inhomo_var_deg: usize,
    /// Maximum degree in t for the inhomogeneous term when enabled.
    pub max_inhomo_var_deg: usize,
    /// Minimum degree in n for the inhomogeneous term when enabled.
    pub min_inhomo_idx_deg: usize,
    /// Maximum degree in n for the inhomogeneous term when enabled.
    pub max_inhomo_idx_deg: usize,
    /// Also search with LHS denominators.
    pub try_denominator: bool,
    /// Also search recurrence terms multiplied by (-1)^n.
    pub try_alternating_sign: bool,
    /// Maximum denom_var_deg when try_denominator is true.
    pub max_denom_var_deg: usize,
    /// Maximum denom_idx_deg when try_denominator is true.
    pub max_denom_idx_deg: usize,
    /// Minimum surplus: equations - unknowns must be >= this.
    pub min_margin: usize,
    /// Use all input rows for fitting and skip held-out verification.
    pub no_verify: bool,
    /// Extra rows to add after the first prefix that clears `min_margin`.
    pub fit_extra_rows: usize,
    /// Probabilistically reject inconsistent candidates modulo large primes
    /// before exact rational solving.
    pub modular_prefilter: bool,
    /// Print each candidate tried to stderr.
    pub verbose: bool,
}

impl Default for AdaptiveSearchOptions {
    fn default() -> Self {
        Self {
            skip_prefix: 0,
            min_rec_len: 1,
            max_rec_len: 5,
            min_var_deg: 0,
            max_var_deg: 3,
            min_idx_deg: 0,
            max_idx_deg: 3,
            min_diff_deg: 0,
            max_diff_deg: 2,
            try_inhomogeneous: false,
            min_inhomo_var_deg: 0,
            max_inhomo_var_deg: 3,
            min_inhomo_idx_deg: 0,
            max_inhomo_idx_deg: 3,
            try_denominator: false,
            try_alternating_sign: false,
            max_denom_var_deg: 2,
            max_denom_idx_deg: 2,
            min_margin: 1,
            no_verify: false,
            fit_extra_rows: 1,
            modular_prefilter: false,
            verbose: false,
        }
    }
}

impl AdaptiveSearchOptions {
    /// Return these search options with alternating `(-1)^n` terms enabled or
    /// disabled.
    pub fn with_alternating_sign(mut self, try_alternating_sign: bool) -> Self {
        self.try_alternating_sign = try_alternating_sign;
        self
    }
}

/// Result of an adaptive recurrence search, including metadata.
#[derive(Debug, Clone)]
pub struct AdaptiveSearchResult {
    /// The recurrence found.
    pub recurrence: Recurrence,
    /// The options that produced it.
    pub opts: RecurrenceOptions,
    /// Number of unknowns in the winning system.
    pub num_unknowns: usize,
    /// Weighted unknown count used for adaptive ordering.
    pub weighted_unknowns: usize,
    /// Number of equations in the winning system.
    pub num_equations: usize,
    /// Number of polynomials used to fit the winning system after skip_prefix.
    pub fit_polynomials: usize,
    /// Number of held-out polynomials verified after fitting.
    pub verification_polynomials: usize,
    /// Number of candidates actually solved (not just counted).
    pub candidates_tried: usize,
    /// Counters for the adaptive search stages up to the returned result.
    pub diagnostics: AdaptiveSearchDiagnostics,
}

/// Counters describing how adaptive recurrence search reached a result.
#[derive(Debug, Clone, Default)]
pub struct AdaptiveSearchDiagnostics {
    /// Number of parameter candidates generated in the searched candidate lists.
    pub generated_candidates: usize,
    /// Number of candidates inspected before the search returned.
    pub considered_candidates: usize,
    /// Candidates skipped because no prefix had enough rows/equations.
    pub insufficient_fit_rows: usize,
    /// Candidates rejected by the equation-count margin after selecting a fit prefix.
    pub equation_bound_rejections: usize,
    /// Homogeneous candidates rejected by the degree-bound precheck.
    pub degree_bound_rejections: usize,
    /// Candidates rejected by the modular prefilter before exact rational solving.
    pub modular_prefilter_rejections: usize,
    /// Candidates sent to exact rational recurrence solving.
    pub exact_solve_attempts: usize,
    /// Exact rational solves that did not produce a recurrence.
    pub failed_exact_solves: usize,
    /// Prefix fits that failed held-out verification.
    pub heldout_verification_failures: usize,
    /// Whether the automatic denominator escalation pass was entered.
    pub denominator_escalation_entered: bool,
}

/// Generate candidate parameter sets.
///
/// Candidates are sorted by weighted parameter count. The raw parameter count
/// is still used later for the linear-algebra solvability check; this ordering
/// only delays more complicated explanations such as derivatives,
/// alternating-sign terms, denominators, and inhomogeneous terms.
fn generate_candidates(m: usize, search: &AdaptiveSearchOptions) -> Vec<RecurrenceOptions> {
    let min_rl = search.min_rec_len.max(1).min(m.saturating_sub(1));
    let max_rl = search.max_rec_len.min(m.saturating_sub(1));
    let mut candidates = Vec::new();

    if min_rl > max_rl {
        return candidates;
    }

    for rec_len in min_rl..=max_rl {
        for diff_deg in search.min_diff_deg..=search.max_diff_deg {
            for idx_deg in search.min_idx_deg..=search.max_idx_deg {
                for var_deg in search.min_var_deg..=search.max_var_deg {
                    let alternating_choices = if search.try_alternating_sign {
                        vec![false, true]
                    } else {
                        vec![false]
                    };
                    for alternating_sign in alternating_choices {
                        candidates.push(RecurrenceOptions {
                            rec_len,
                            var_deg,
                            idx_deg,
                            diff_deg,
                            homogeneous: true,
                            inhomo_var_deg: 0,
                            inhomo_idx_deg: 0,
                            denom_var_deg: 0,
                            denom_idx_deg: 0,
                            alternating_sign,
                            modular_prefilter: search.modular_prefilter,
                        });

                        if search.try_inhomogeneous {
                            for inhomo_idx_deg in
                                search.min_inhomo_idx_deg..=search.max_inhomo_idx_deg
                            {
                                for inhomo_var_deg in
                                    search.min_inhomo_var_deg..=search.max_inhomo_var_deg
                                {
                                    candidates.push(RecurrenceOptions {
                                        rec_len,
                                        var_deg,
                                        idx_deg,
                                        diff_deg,
                                        homogeneous: false,
                                        inhomo_var_deg,
                                        inhomo_idx_deg,
                                        denom_var_deg: 0,
                                        denom_idx_deg: 0,
                                        alternating_sign,
                                        modular_prefilter: search.modular_prefilter,
                                    });
                                }
                            }
                        }

                        if search.try_denominator {
                            for dvd in 0..=search.max_denom_var_deg {
                                for did in 0..=search.max_denom_idx_deg {
                                    if dvd == 0 && did == 0 {
                                        continue; // already covered above
                                    }
                                    candidates.push(RecurrenceOptions {
                                        rec_len,
                                        var_deg,
                                        idx_deg,
                                        diff_deg,
                                        homogeneous: true,
                                        inhomo_var_deg: 0,
                                        inhomo_idx_deg: 0,
                                        denom_var_deg: dvd,
                                        denom_idx_deg: did,
                                        alternating_sign,
                                        modular_prefilter: search.modular_prefilter,
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    candidates.sort_by_key(|opts| {
        let complexity = candidate_complexity(opts);
        let has_denominator = opts.denom_var_deg > 0 || opts.denom_idx_deg > 0;
        (
            complexity.weighted_unknowns,
            complexity.raw_unknowns,
            opts.alternating_sign,
            has_denominator,
            !opts.homogeneous,
            opts.diff_deg,
            opts.rec_len,
            opts.idx_deg,
            opts.var_deg,
            opts.denom_idx_deg + opts.denom_var_deg,
            opts.inhomo_idx_deg + opts.inhomo_var_deg,
        )
    });
    candidates
}

fn fitting_polynomial_count(
    polys: &[Vec<BigRational>],
    opts: &RecurrenceOptions,
    search: &AdaptiveSearchOptions,
) -> Option<usize> {
    let m = polys.len();
    if m <= opts.rec_len {
        return None;
    }

    if search.no_verify {
        return Some(m);
    }

    let max_fit = m.checked_sub(1)?;
    if max_fit <= opts.rec_len {
        return None;
    }

    let unknowns = count_unknowns(opts);
    for fit_len in opts.rec_len + 1..=max_fit {
        let equations = count_equations_rational(&polys[..fit_len], opts);
        if equations >= unknowns + search.min_margin {
            return Some((fit_len + search.fit_extra_rows).min(max_fit));
        }
    }

    None
}

fn verify_heldout_tail(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    rec: &Recurrence,
    opts: &RecurrenceOptions,
    fit_len: usize,
    search: &AdaptiveSearchOptions,
) -> bool {
    search.no_verify
        || recurrence_holds_from_rational_with_derivs(polys, derivs, rec, opts.rec_len, fit_len + 1)
}

fn homogeneous_degree_bound_rejects(
    polys: &[Vec<BigRational>],
    derivs: &[Vec<Vec<BigRational>>],
    opts: &RecurrenceOptions,
    fit_len: usize,
) -> bool {
    if !opts.homogeneous || opts.denom_var_deg > 0 || opts.denom_idx_deg > 0 {
        return false;
    }

    for nn in opts.rec_len + 1..=fit_len {
        let current = &polys[nn - 1];
        if poly_is_zero_rational(current) {
            continue;
        }
        let lhs_degree = poly_degree_rational(current);
        let mut rhs_degree_bound: Option<usize> = None;
        for r in 1..=opts.rec_len {
            for ref_poly in derivs[nn - 1 - r].iter().take(opts.diff_deg + 1) {
                if !poly_is_zero_rational(ref_poly) {
                    let degree = poly_degree_rational(ref_poly) + opts.var_deg;
                    rhs_degree_bound =
                        Some(rhs_degree_bound.map_or(degree, |bound| bound.max(degree)));
                }
            }
        }
        if rhs_degree_bound.is_none_or(|degree| degree < lhs_degree) {
            return true;
        }
    }

    false
}

/// Search for the simplest polynomial recurrence by trying parameter
/// combinations in order of ascending complexity.
pub fn find_recurrence_adaptive_rational(
    polys: &[Vec<BigRational>],
    search: &AdaptiveSearchOptions,
) -> Option<AdaptiveSearchResult> {
    let polys = polys.get(search.skip_prefix..).unwrap_or(&[]);
    let m = polys.len();
    if m < 2 {
        return None;
    }
    let derivs = rational_derivatives_up_to(polys, search.max_diff_deg);
    let modular_cache = search
        .modular_prefilter
        .then(|| ModularPrefilterCache::new(polys, &derivs, search.max_diff_deg));

    // First pass: use the options as given.
    let candidates = generate_candidates(m, search);
    let mut tried = 0;
    let mut diagnostics = AdaptiveSearchDiagnostics {
        generated_candidates: candidates.len(),
        ..Default::default()
    };

    for opts in &candidates {
        diagnostics.considered_candidates += 1;
        let complexity = candidate_complexity(opts);
        let unknowns = complexity.raw_unknowns;
        let weighted_unknowns = complexity.weighted_unknowns;
        let Some(fit_len) = fitting_polynomial_count(polys, opts, search) else {
            diagnostics.insufficient_fit_rows += 1;
            continue;
        };
        let fit_polys = &polys[..fit_len];
        let equations = count_equations_rational(fit_polys, opts);

        if equations < unknowns + search.min_margin {
            diagnostics.equation_bound_rejections += 1;
            continue;
        }

        if homogeneous_degree_bound_rejects(polys, &derivs, opts, fit_len) {
            diagnostics.degree_bound_rejections += 1;
            continue;
        }

        tried += 1;

        if search.verbose {
            eprintln!(
                "  try #{tried}: rec_len={} var_deg={} idx_deg={} diff_deg={} \
                 alternating={} denom=({},{}) homog={} \
                 inhomo=({},{}) fit_rows={} verify_rows={} \
                 (unknowns={unknowns}, weighted={weighted_unknowns}, \
                 equations={equations}, margin={})",
                opts.rec_len,
                opts.var_deg,
                opts.idx_deg,
                opts.diff_deg,
                opts.alternating_sign,
                opts.denom_var_deg,
                opts.denom_idx_deg,
                opts.homogeneous,
                opts.inhomo_var_deg,
                opts.inhomo_idx_deg,
                fit_len,
                m - fit_len,
                equations - unknowns,
            );
        }

        let fit_derivs = &derivs[..fit_len];
        let mut exact_row_indices = None;
        let cached_solve_opts = if let Some(cache) = &modular_cache {
            let prefilter = modular_prefilter_with_cache(polys, fit_len, opts, cache);
            if prefilter.rejected {
                diagnostics.modular_prefilter_rejections += 1;
                continue;
            }
            exact_row_indices = prefilter.full_rank_pivot_rows;
            let mut solve_opts = opts.clone();
            solve_opts.modular_prefilter = false;
            Some(solve_opts)
        } else {
            None
        };
        let solve_opts = cached_solve_opts.as_ref().unwrap_or(opts);
        diagnostics.exact_solve_attempts += 1;
        if let Some(rec) = find_polynomial_recurrence_rational_with_derivs_and_rows(
            fit_polys,
            fit_derivs,
            solve_opts,
            exact_row_indices.as_deref(),
        ) {
            if !verify_heldout_tail(polys, &derivs, &rec, opts, fit_len, search) {
                diagnostics.heldout_verification_failures += 1;
                if search.verbose {
                    eprintln!("  -> fitted prefix but failed held-out verification");
                }
                continue;
            }
            if search.verbose {
                eprintln!("  -> found!");
            }
            return Some(AdaptiveSearchResult {
                recurrence: rec,
                opts: opts.clone(),
                num_unknowns: unknowns,
                weighted_unknowns,
                num_equations: equations,
                fit_polynomials: fit_len,
                verification_polynomials: m - fit_len,
                candidates_tried: tried,
                diagnostics,
            });
        } else {
            diagnostics.failed_exact_solves += 1;
        }
    }

    // Second pass: if denominator wasn't already tried, automatically
    // escalate to rational coefficients (LHS denominator f(n,t)).
    // Only try if we have enough data to avoid spurious fits.
    if !search.try_denominator && m >= 6 {
        let mut rational_search = search.clone();
        rational_search.try_denominator = true;
        // Use moderate denominator degrees if not already set.
        if rational_search.max_denom_idx_deg == 0 {
            rational_search.max_denom_idx_deg = 1;
        }

        let rational_candidates = generate_candidates(m, &rational_search);
        diagnostics.denominator_escalation_entered = true;
        diagnostics.generated_candidates += rational_candidates
            .iter()
            .filter(|opts| opts.denom_var_deg > 0 || opts.denom_idx_deg > 0)
            .count();
        for opts in &rational_candidates {
            // Skip candidates without a denominator (already tried above).
            if opts.denom_var_deg == 0 && opts.denom_idx_deg == 0 {
                continue;
            }

            diagnostics.considered_candidates += 1;
            let complexity = candidate_complexity(opts);
            let unknowns = complexity.raw_unknowns;
            let weighted_unknowns = complexity.weighted_unknowns;
            let Some(fit_len) = fitting_polynomial_count(polys, opts, &rational_search) else {
                diagnostics.insufficient_fit_rows += 1;
                continue;
            };
            let fit_polys = &polys[..fit_len];
            let equations = count_equations_rational(fit_polys, opts);

            if equations < unknowns + search.min_margin {
                diagnostics.equation_bound_rejections += 1;
                continue;
            }

            if homogeneous_degree_bound_rejects(polys, &derivs, opts, fit_len) {
                diagnostics.degree_bound_rejections += 1;
                continue;
            }

            tried += 1;

            if search.verbose {
                eprintln!(
                    "  try #{tried} (rational): rec_len={} var_deg={} idx_deg={} diff_deg={} \
                     alternating={} denom=({},{}) \
                     fit_rows={} verify_rows={} \
                     (unknowns={unknowns}, weighted={weighted_unknowns}, \
                     equations={equations}, margin={})",
                    opts.rec_len,
                    opts.var_deg,
                    opts.idx_deg,
                    opts.diff_deg,
                    opts.alternating_sign,
                    opts.denom_var_deg,
                    opts.denom_idx_deg,
                    fit_len,
                    m - fit_len,
                    equations - unknowns,
                );
            }

            let fit_derivs = &derivs[..fit_len];
            let mut exact_row_indices = None;
            let cached_solve_opts = if let Some(cache) = &modular_cache {
                let prefilter = modular_prefilter_with_cache(polys, fit_len, opts, cache);
                if prefilter.rejected {
                    diagnostics.modular_prefilter_rejections += 1;
                    continue;
                }
                exact_row_indices = prefilter.full_rank_pivot_rows;
                let mut solve_opts = opts.clone();
                solve_opts.modular_prefilter = false;
                Some(solve_opts)
            } else {
                None
            };
            let solve_opts = cached_solve_opts.as_ref().unwrap_or(opts);
            diagnostics.exact_solve_attempts += 1;
            if let Some(rec) = find_polynomial_recurrence_rational_with_derivs_and_rows(
                fit_polys,
                fit_derivs,
                solve_opts,
                exact_row_indices.as_deref(),
            ) {
                if !verify_heldout_tail(polys, &derivs, &rec, opts, fit_len, search) {
                    diagnostics.heldout_verification_failures += 1;
                    if search.verbose {
                        eprintln!("  -> fitted prefix but failed held-out verification");
                    }
                    continue;
                }
                if search.verbose {
                    eprintln!("  -> found (rational)!");
                }
                return Some(AdaptiveSearchResult {
                    recurrence: rec,
                    opts: opts.clone(),
                    num_unknowns: unknowns,
                    weighted_unknowns,
                    num_equations: equations,
                    fit_polynomials: fit_len,
                    verification_polynomials: m - fit_len,
                    candidates_tried: tried,
                    diagnostics,
                });
            } else {
                diagnostics.failed_exact_solves += 1;
            }
        }
    }

    None
}

/// Search for the simplest polynomial recurrence for an integer-coefficient
/// polynomial sequence.
///
/// This is a convenience wrapper around [`find_recurrence_adaptive_rational`].
pub fn find_recurrence_adaptive(
    polys: &[Vec<i64>],
    search: &AdaptiveSearchOptions,
) -> Option<AdaptiveSearchResult> {
    let rational_polys = i64_polys_to_rational(polys);
    find_recurrence_adaptive_rational(&rational_polys, search)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: assert a recurrence is found and matches expected display.
    fn assert_recurrence(polys: &[Vec<i64>], opts: &RecurrenceOptions, expected: &str) {
        let rec = find_polynomial_recurrence(polys, opts)
            .unwrap_or_else(|| panic!("Expected recurrence, got None"));
        let display = format!("{rec}");
        assert_eq!(display, expected, "\npolys: {:?}\nopts: {:?}", polys, opts);
    }

    fn br(n: i64, d: i64) -> BigRational {
        BigRational::new(BigInt::from(n), BigInt::from(d))
    }

    #[test]
    fn rational_input_geometric() {
        let polys: Vec<Vec<BigRational>> = vec![
            vec![br(1, 2)],
            vec![br(1, 4)],
            vec![br(1, 8)],
            vec![br(1, 16)],
            vec![br(1, 32)],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence_rational(&polys, &opts)
            .expect("should find rational recurrence");
        assert_eq!(format!("{rec}"), "P(n) = 1/2 P(n-1)");
    }

    #[test]
    fn rational_input_large_coefficients() {
        let base = BigInt::one() << 100usize;
        let polys: Vec<Vec<BigRational>> = (0..5)
            .map(|i| vec![BigRational::from_integer(&base << i)])
            .collect();
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence_rational(&polys, &opts)
            .expect("should find recurrence with large coefficients");
        assert_eq!(format!("{rec}"), "P(n) = 2 P(n-1)");
    }

    #[test]
    fn modular_prefilter_rejects_inconsistent_scalar_candidate() {
        let polys = i64_polys_to_rational(&[vec![1], vec![2], vec![5]]);
        let derivs = polys
            .iter()
            .map(|p| vec![poly_nth_derivative_rational(p, 0)])
            .collect::<Vec<_>>();
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            modular_prefilter: true,
            ..Default::default()
        };
        assert!(modular_prefilter_rejects(&polys, &derivs, &opts));
        assert!(find_polynomial_recurrence_rational(&polys, &opts).is_none());
    }

    #[test]
    fn modular_prefilter_preserves_true_recurrence() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            modular_prefilter: true,
            ..Default::default()
        };
        assert_recurrence(&polys, &opts, "P(n) = P(n-1) + P(n-2)");
    }

    #[test]
    fn modular_tail_verification_checks_solution_on_heldout_rows() {
        let polys = i64_polys_to_rational(&[
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
        ]);
        let derivs = rational_derivatives_up_to(&polys, 0);
        let modulus = MODULAR_PREFILTER_PRIMES[0];
        let polys_mod = rational_polys_mod_prime(&polys, modulus).expect("modular polys");
        let derivs_mod =
            rational_derivs_mod_prime(&derivs, modulus, 0).expect("modular derivatives");
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            modular_prefilter: true,
            ..Default::default()
        };

        assert!(recurrence_solution_holds_from_mod_images(
            &polys,
            &polys_mod,
            &derivs_mod,
            &opts,
            modulus,
            &[1, 1],
            4,
        ));
        assert!(!recurrence_solution_holds_from_mod_images(
            &polys,
            &polys_mod,
            &derivs_mod,
            &opts,
            modulus,
            &[1, 2],
            4,
        ));
    }

    #[test]
    fn modular_prefilter_rejects_rank_deficient_prefix_with_bad_heldout_row() {
        let polys = i64_polys_to_rational(&[vec![0], vec![0], vec![0], vec![1]]);
        let derivs = rational_derivatives_up_to(&polys, 0);
        let cache = ModularPrefilterCache::new(&polys, &derivs, 0);
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            modular_prefilter: true,
            ..Default::default()
        };

        let prefilter = modular_prefilter_with_cache(&polys, 3, &opts, &cache);
        assert!(prefilter.rejected);
        assert!(prefilter.full_rank_pivot_rows.is_none());
    }

    #[test]
    fn modular_prefilter_tolerates_single_bad_prime() {
        let bad_prime = MODULAR_PREFILTER_PRIMES[0];
        let polys = i64_polys_to_rational(&[vec![bad_prime * bad_prime], vec![bad_prime], vec![1]]);
        let derivs = polys
            .iter()
            .map(|p| vec![poly_nth_derivative_rational(p, 0)])
            .collect::<Vec<_>>();
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            modular_prefilter: true,
            ..Default::default()
        };
        assert!(!modular_prefilter_rejects(&polys, &derivs, &opts));
        assert!(find_polynomial_recurrence_rational(&polys, &opts).is_some());
    }

    #[test]
    fn fibonacci() {
        // P_n = P_{n-1} + P_{n-2}
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            ..Default::default()
        };
        assert_recurrence(&polys, &opts, "P(n) = P(n-1) + P(n-2)");
    }

    #[test]
    fn binomial_expansion() {
        // P_n = (1+t)^n, so P_n = (1+t) P_{n-1}
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1, 1],
            vec![1, 2, 1],
            vec![1, 3, 3, 1],
            vec![1, 4, 6, 4, 1],
        ];
        let opts = RecurrenceOptions {
            var_deg: 1,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        assert_recurrence(&polys, &opts, "P(n) = (1 + t) P(n-1)");
    }

    #[test]
    fn factorial() {
        // P_n = n! (as constant polynomials): P_n = n P_{n-1}
        // Using 1-based indexing: P_1=1, P_2=2, P_3=6, P_4=24, P_5=120
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![6], vec![24], vec![120]];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 1,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        assert_recurrence(&polys, &opts, "P(n) = n P(n-1)");
    }

    #[test]
    fn chebyshev() {
        // T_n(t) = 2t T_{n-1}(t) - T_{n-2}(t)
        // T_0=1, T_1=t, T_2=2t^2-1, T_3=4t^3-3t, T_4=8t^4-8t^2+1
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![0, 1],
            vec![-1, 0, 2],
            vec![0, -3, 0, 4],
            vec![1, 0, -8, 0, 8],
        ];
        let opts = RecurrenceOptions {
            var_deg: 1,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            ..Default::default()
        };
        assert_recurrence(&polys, &opts, "P(n) = 2t P(n-1) - P(n-2)");
    }

    #[test]
    fn eulerian_with_derivative() {
        // A_n(t) = (1 + (n-1)t) A_{n-1}(t) + t(1-t) A'_{n-1}(t)
        // A_1=1, A_2=1, A_3=1+t, A_4=1+4t+t^2, A_5=1+11t+11t^2+t^3
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![1, 1],
            vec![1, 4, 1],
            vec![1, 11, 11, 1],
            vec![1, 26, 66, 26, 1],
        ];
        let opts = RecurrenceOptions {
            var_deg: 2,
            idx_deg: 1,
            diff_deg: 1,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence(&polys, &opts).expect("should find recurrence");
        let display = format!("{rec}");
        // With 1-based indexing: A_0=P_1, A_1=P_2, etc.
        // So the recurrence becomes P_n = [1 + (n-2)t] P_{n-1} + t(1-t) P'_{n-1}
        // Coefficient of P(n-1): 1 - 2t + nt
        // Coefficient of P'(n-1): t - t^2
        assert_eq!(display, "P(n) = (1 - 2t + nt) P(n-1) + (t - t^2) P'(n-1)");
    }

    #[test]
    fn mathematica_definition_contains_initial_values_and_rule() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence(&polys, &opts).expect("should find recurrence");
        let code = rec.to_mathematica_definition(&polys);
        assert!(code.contains("P[1, t_] := 1;"));
        assert!(code.contains("P[2, t_] := 1;"));
        assert!(code.contains(
            "P[n_Integer /; n >= 3, t_] := P[n, t] = Expand[P[n - 1, t] + P[n - 2, t]];"
        ));
    }

    #[test]
    fn sage_definition_contains_derivative_rule() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![1, 1],
            vec![1, 4, 1],
            vec![1, 11, 11, 1],
            vec![1, 26, 66, 26, 1],
        ];
        let opts = RecurrenceOptions {
            var_deg: 2,
            idx_deg: 1,
            diff_deg: 1,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence(&polys, &opts).expect("should find recurrence");
        let code = rec.to_sage_definition(&polys);
        assert!(code.contains("R.<t> = PolynomialRing(QQ)"));
        assert!(code.contains("_P_cache = {"));
        assert!(code.contains("1: R([1]),"));
        assert!(code.contains("value = R("));
        assert!(code.contains("P(n - 1).derivative(t)"));
    }

    #[test]
    fn python_definition_contains_exact_fraction_rule() {
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![1], vec![2], vec![3], vec![5], vec![8]];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence(&polys, &opts).expect("should find recurrence");
        let code = rec.to_python_definition(&polys);
        assert!(code.contains("from fractions import Fraction"));
        assert!(code.contains("def P(n):"));
        assert!(code.contains("_P_cache = {"));
        assert!(code.contains("1: [1],"));
        assert!(code.contains("\"offset\": 1"));
        assert!(code.contains("_div_exact"));
    }

    #[test]
    fn python_definition_contains_derivative_metadata() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![1, 1],
            vec![1, 4, 1],
            vec![1, 11, 11, 1],
            vec![1, 26, 66, 26, 1],
        ];
        let opts = RecurrenceOptions {
            var_deg: 2,
            idx_deg: 1,
            diff_deg: 1,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence(&polys, &opts).expect("should find recurrence");
        let code = rec.to_python_definition(&polys);
        assert!(code.contains("\"deriv_order\": 1"));
        assert!(code.contains("ref = _derivative(P(k - term[\"offset\"]),"));
        assert!(code.contains("coeff = _bivar(term[\"coeff\"], k)"));
    }

    #[test]
    fn geometric_sequence() {
        // P_n = 2^{n-1}: P_n = 2 P_{n-1}
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![2],
            vec![4],
            vec![8],
            vec![16],
            vec![32],
            vec![64],
            vec![128],
            vec![256],
            vec![512],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        assert_recurrence(&polys, &opts, "P(n) = 2 P(n-1)");
    }

    #[test]
    fn alternating_sign_terms() {
        // P_n = (-1)^n P_{n-1}; with P_1=1 this gives pairs 1,1,-1,-1,...
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![-1],
            vec![-1],
            vec![1],
            vec![1],
            vec![-1],
            vec![-1],
        ];
        let ordinary_opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        assert!(find_polynomial_recurrence(&polys, &ordinary_opts).is_none());

        let alternating_opts = RecurrenceOptions {
            alternating_sign: true,
            ..ordinary_opts
        };
        let rec = find_polynomial_recurrence(&polys, &alternating_opts)
            .expect("should find alternating-sign recurrence");
        assert_eq!(format!("{rec}"), "P(n) = (-1)^n P(n-1)");
        assert_eq!(rec.to_latex(), "P(n) = (-1)^n P(n-1)");
        assert!(rec
            .to_mathematica_definition(&polys)
            .contains("Expand[(-1)^n*P[n - 1, t]]"));
        assert!(rec
            .to_sage_definition(&polys)
            .contains("value = R((-1)**n*P(n - 1))"));
    }

    #[test]
    fn constant_coeffs_with_separate_inhomogeneous_bounds() {
        // P_n = P_{n-1} + n + t^2
        let polys: Vec<Vec<i64>> = vec![
            vec![0],
            vec![2, 0, 1],
            vec![5, 0, 2],
            vec![9, 0, 3],
            vec![14, 0, 4],
            vec![20, 0, 5],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: false,
            inhomo_var_deg: 2,
            inhomo_idx_deg: 1,
            ..Default::default()
        };
        assert_recurrence(&polys, &opts, "P(n) = P(n-1) + n + t^2");
    }

    #[test]
    fn adaptive_geometric() {
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![4], vec![8], vec![16], vec![32]];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = 2 P(n-1)");
        assert!(result.verification_polynomials >= 1);
        assert!(result.fit_polynomials < polys.len());
    }

    #[test]
    fn homogeneous_degree_bound_rejects_impossible_rhs_degree() {
        let polys = i64_polys_to_rational(&[vec![1], vec![1, 1], vec![1, 2, 1]]);
        let derivs = rational_derivatives_up_to(&polys, 0);
        let too_low = RecurrenceOptions {
            rec_len: 1,
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            homogeneous: true,
            ..Default::default()
        };
        let enough_t_degree = RecurrenceOptions {
            var_deg: 1,
            ..too_low.clone()
        };

        assert!(homogeneous_degree_bound_rejects(
            &polys, &derivs, &too_low, 3
        ));
        assert!(!homogeneous_degree_bound_rejects(
            &polys,
            &derivs,
            &enough_t_degree,
            3
        ));
    }

    #[test]
    fn adaptive_no_verify_uses_all_rows() {
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![4], vec![8], vec![16], vec![32]];
        let search = AdaptiveSearchOptions {
            no_verify: true,
            ..Default::default()
        };
        let result = find_recurrence_adaptive(&polys, &search).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = 2 P(n-1)");
        assert_eq!(result.fit_polynomials, polys.len());
        assert_eq!(result.verification_polynomials, 0);
    }

    #[test]
    fn adaptive_rejects_false_prefix_fit_on_heldout_row() {
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![4], vec![8], vec![17]];
        let search = AdaptiveSearchOptions {
            max_rec_len: 1,
            max_var_deg: 0,
            max_idx_deg: 0,
            max_diff_deg: 0,
            ..Default::default()
        };
        assert!(find_recurrence_adaptive(&polys, &search).is_none());
    }

    #[test]
    fn weighted_complexity_delays_alternating_terms() {
        let ordinary = RecurrenceOptions {
            rec_len: 2,
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            homogeneous: true,
            ..Default::default()
        };
        let alternating = RecurrenceOptions {
            rec_len: 1,
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            homogeneous: true,
            alternating_sign: true,
            ..Default::default()
        };
        assert_eq!(count_unknowns(&ordinary), count_unknowns(&alternating));
        assert!(count_weighted_unknowns(&ordinary) < count_weighted_unknowns(&alternating));
    }

    #[test]
    fn adaptive_respects_min_bounds() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
        ];
        let search = AdaptiveSearchOptions {
            min_rec_len: 2,
            max_rec_len: 2,
            min_var_deg: 0,
            max_var_deg: 0,
            min_idx_deg: 0,
            max_idx_deg: 0,
            min_diff_deg: 0,
            max_diff_deg: 0,
            ..Default::default()
        };
        let result = find_recurrence_adaptive(&polys, &search).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = P(n-1) + P(n-2)");
        assert_eq!(result.opts.rec_len, 2);
    }

    #[test]
    fn adaptive_can_skip_prefix() {
        let polys: Vec<Vec<i64>> = vec![
            vec![9],
            vec![1],
            vec![2],
            vec![4],
            vec![8],
            vec![16],
            vec![32],
        ];
        let search = AdaptiveSearchOptions {
            skip_prefix: 1,
            ..Default::default()
        };
        let result = find_recurrence_adaptive(&polys, &search).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = 2 P(n-1)");
    }

    #[test]
    fn no_recurrence() {
        // Random polynomials unlikely to satisfy a short recurrence.
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1, 7],
            vec![3, 1, 5],
            vec![2, 9, 1, 4],
            vec![7, 2, 8, 1, 3],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            ..Default::default()
        };
        assert!(find_polynomial_recurrence(&polys, &opts).is_none());
    }

    #[test]
    fn index_dependent_coefficients() {
        // P_n = (2n-1) P_{n-1} - (n-1)^2 P_{n-2}  (Legendre-related)
        // P_1=1, P_2=1, then P_3 = 3*1 - 1*1 = 2, P_4 = 5*2 - 4*1 = 6,
        // P_5 = 7*6 - 9*2 = 24, P_6 = 9*24 - 16*6 = 120
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![6],
            vec![24],
            vec![120],
            vec![720],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 2,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence(&polys, &opts);
        // This should find SOME recurrence (the simplest might be P_n = n P_{n-1}).
        assert!(rec.is_some());
    }

    // --- Adaptive search tests ---

    #[test]
    fn adaptive_fibonacci() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
        ];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = P(n-1) + P(n-2)");
        assert_eq!(result.opts.rec_len, 2);
        assert_eq!(result.opts.diff_deg, 0);
    }

    #[test]
    fn adaptive_modular_prefilter_finds_fibonacci() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
        ];
        let search = AdaptiveSearchOptions {
            modular_prefilter: true,
            ..Default::default()
        };
        let result = find_recurrence_adaptive(&polys, &search).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = P(n-1) + P(n-2)");
        assert!(result.opts.modular_prefilter);
    }

    #[test]
    fn adaptive_factorial() {
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![6], vec![24], vec![120]];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = n P(n-1)");
        assert_eq!(result.opts.rec_len, 1);
        assert_eq!(result.opts.idx_deg, 1);
    }

    #[test]
    fn adaptive_binomial() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1, 1],
            vec![1, 2, 1],
            vec![1, 3, 3, 1],
            vec![1, 4, 6, 4, 1],
        ];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = (1 + t) P(n-1)");
        assert_eq!(result.opts.rec_len, 1);
    }

    #[test]
    fn adaptive_chebyshev() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![0, 1],
            vec![-1, 0, 2],
            vec![0, -3, 0, 4],
            vec![1, 0, -8, 0, 8],
        ];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        assert_eq!(
            format!("{}", result.recurrence),
            "P(n) = 2t P(n-1) - P(n-2)"
        );
    }

    #[test]
    fn adaptive_eulerian() {
        // With held-out verification, prefix-only non-derivative fits are
        // rejected, and the standard derivative recurrence is found.
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![1, 1],
            vec![1, 4, 1],
            vec![1, 11, 11, 1],
            vec![1, 26, 66, 26, 1],
        ];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        assert_eq!(
            format!("{}", result.recurrence),
            "P(n) = (1 - 2t + nt) P(n-1) + (t - t^2) P'(n-1)"
        );
        assert_eq!(result.opts.diff_deg, 1);
        assert_eq!(result.opts.rec_len, 1);

        // With max_rec_len=1, the same derivative recurrence is still found.
        let search = AdaptiveSearchOptions {
            max_rec_len: 1,
            max_diff_deg: 2,
            ..Default::default()
        };
        let result = find_recurrence_adaptive(&polys, &search).unwrap();
        assert_eq!(
            format!("{}", result.recurrence),
            "P(n) = (1 - 2t + nt) P(n-1) + (t - t^2) P'(n-1)"
        );
        assert_eq!(result.opts.diff_deg, 1);
    }

    #[test]
    fn adaptive_alternating_sign() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![-1],
            vec![-1],
            vec![1],
            vec![1],
            vec![-1],
            vec![-1],
        ];
        let search = AdaptiveSearchOptions {
            try_alternating_sign: true,
            max_rec_len: 1,
            max_var_deg: 0,
            max_idx_deg: 0,
            max_diff_deg: 0,
            ..Default::default()
        };
        let result = find_recurrence_adaptive(&polys, &search).unwrap();
        assert_eq!(format!("{}", result.recurrence), "P(n) = (-1)^n P(n-1)");
        assert!(result.opts.alternating_sign);
    }

    #[test]
    fn adaptive_short_sequence() {
        // m=3 with constant polys: too few equations for any 2-term recurrence.
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![1], vec![2]];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default());
        assert!(result.is_none());
    }

    #[test]
    fn recurrence_json_roundtrip_generates_rows() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
            vec![21],
        ];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        let rational_polys = i64_polys_to_rational(&polys);
        let initial_count = result.recurrence.max_offset();
        let recurrence_json = RecurrenceJson::from_recurrence_rational(
            &result.recurrence,
            1,
            &rational_polys[..initial_count],
            Some(RecurrenceJsonSearch {
                recurrence_text: result.recurrence.to_string(),
                source_rows: polys.len(),
                skip_prefix: 0,
                unknowns: result.num_unknowns,
                weighted_unknowns: result.weighted_unknowns,
                equations: result.num_equations,
                fit_polynomials: result.fit_polynomials,
                verification_polynomials: result.verification_polynomials,
                candidates_tried: result.candidates_tried,
                options: RecurrenceOptionsJson::from(&result.opts),
            }),
        );

        let encoded = serde_json::to_string(&recurrence_json).unwrap();
        let decoded: RecurrenceJson = serde_json::from_str(&encoded).unwrap();
        let (recurrence, first_index, initial_polys) = decoded.to_recurrence_parts().unwrap();
        let generated = recurrence
            .generate_rows_rational(&initial_polys, first_index, polys.len())
            .unwrap();
        assert_eq!(generated, rational_polys);
    }

    #[test]
    fn recurrence_json_generates_denominator_rows() {
        let mut polys: Vec<Vec<BigRational>> = vec![vec![BigRational::one()]];
        for n in 2..=8 {
            let prev = polys.last().unwrap()[0].clone();
            polys.push(vec![prev / BigRational::from_integer(BigInt::from(n + 1))]);
        }
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            denom_var_deg: 0,
            denom_idx_deg: 1,
            ..Default::default()
        };
        let recurrence =
            find_polynomial_recurrence_rational(&polys, &opts).expect("should find denominator");
        let recurrence_json =
            RecurrenceJson::from_recurrence_rational(&recurrence, 1, &polys[..1], None);
        let (recurrence, first_index, initial_polys) =
            recurrence_json.to_recurrence_parts().unwrap();
        let generated = recurrence
            .generate_rows_rational(&initial_polys, first_index, polys.len())
            .unwrap();
        assert_eq!(generated, polys);
    }

    #[test]
    fn test_count_unknowns() {
        let opts = RecurrenceOptions {
            var_deg: 1,
            idx_deg: 1,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            inhomo_var_deg: 0,
            inhomo_idx_deg: 0,
            denom_var_deg: 0,
            denom_idx_deg: 0,
            alternating_sign: false,
            modular_prefilter: false,
        };
        // vars_per_coeff = 2*2 = 4, num_coeff_vars = 2*1*4 = 8, denom = 0
        assert_eq!(count_unknowns(&opts), 8);
    }

    #[test]
    fn test_count_equations() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            inhomo_var_deg: 0,
            inhomo_idx_deg: 0,
            denom_var_deg: 0,
            denom_idx_deg: 0,
            alternating_sign: false,
            modular_prefilter: false,
        };
        // m=7, rec_len=2, num_nn=5, max_poly_deg=0, max_j=0, eqs_per_nn=1
        assert_eq!(count_equations(&polys, &opts), 5);
    }
}
