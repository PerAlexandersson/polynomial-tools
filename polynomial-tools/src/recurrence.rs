//! Find linear recurrences among sequences of polynomials.
//!
//! Given polynomials P_1(t), P_2(t), ..., P_m(t), searches for a recurrence:
//!
//!   f(n,t) P_n(t) = sum_{r,d} c_{r,d}(n,t) D^d P_{n-r}(t)  [+ g(n,t)]
//!
//! where c_{r,d}(n,t) are polynomial coefficients in n and t, D^d is the
//! d-th derivative in t, and f(n,t) is an optional LHS denominator.
//!
//! This reduces to solving a linear system over the rationals.

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};
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
        }
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

/// One term in a recurrence: coeff(n,t) * D^d P_{n-r}(t).
#[derive(Debug, Clone)]
pub struct RecurrenceTerm {
    /// Recurrence offset r (P_{n-r}).
    pub offset: usize,
    /// Derivative order d.
    pub deriv_order: usize,
    /// Coefficient polynomial c(n,t).
    pub coeff: BivarPoly,
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

use crate::linalg;

// ---------------------------------------------------------------------------
// Main algorithm
// ---------------------------------------------------------------------------

/// Search for a polynomial recurrence satisfied by a sequence of polynomials.
///
/// `polys[i]` is P_{i+1}(t) given as a coefficient vector (index = power of t).
/// Returns `None` if no recurrence is found with the given options.
pub fn find_polynomial_recurrence_rational(
    polys: &[Vec<BigRational>],
    opts: &RecurrenceOptions,
) -> Option<Recurrence> {
    let m = polys.len();
    if m <= opts.rec_len {
        return None;
    }

    // Pre-compute all needed derivatives.
    let derivs: Vec<Vec<Vec<BigRational>>> = polys
        .iter()
        .map(|p| {
            (0..=opts.diff_deg)
                .map(|d| poly_nth_derivative_rational(p, d))
                .collect()
        })
        .collect();

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

    // 2) Recurrence coefficient unknowns c[r][d][i][j].
    let coeff_start = denom_start + num_denom_vars;
    let vars_per_coeff = (opts.idx_deg + 1) * (opts.var_deg + 1);
    let num_coeff_vars = opts.rec_len * (opts.diff_deg + 1) * vars_per_coeff;

    let coeff_col = |r: usize, d: usize, i: usize, j: usize| -> usize {
        coeff_start
            + ((r - 1) * (opts.diff_deg + 1) + d) * vars_per_coeff
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

    let zero = BigRational::zero();
    let mut matrix: Vec<Vec<BigRational>> = vec![vec![zero.clone(); num_vars]; num_rows];
    let mut rhs: Vec<BigRational> = vec![zero.clone(); num_rows];

    for (eq_idx, nn) in (opts.rec_len + 1..=m).enumerate() {
        // nn is 1-based; polys[nn-1] is P_nn(t).
        let current = &polys[nn - 1];

        for l in 0..=max_t_deg {
            let row = eq_idx * eqs_per_nn + l;

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
                        let val = pc
                            * BigRational::from_integer(num_traits::pow::pow(BigInt::from(nn), i));
                        matrix[row][denom_col(i, j)] = val;
                    }
                }
            }

            // Recurrence coefficients: −c[r][d][i][j] * nn^i * (coeff of t^{l-j} in D^d P_{nn-r}).
            for r in 1..=opts.rec_len {
                for d in 0..=opts.diff_deg {
                    let ref_poly = &derivs[nn - 1 - r][d];
                    for i in 0..=opts.idx_deg {
                        let ni = num_traits::pow::pow(BigInt::from(nn as i64), i);
                        for j in 0..=opts.var_deg {
                            if l < j {
                                continue;
                            }
                            let rc = poly_coeff_rational(ref_poly, l - j);
                            if rc.is_zero() {
                                continue;
                            }
                            let val = -rc * BigRational::from_integer(ni.clone());
                            matrix[row][coeff_col(r, d, i, j)] = val;
                        }
                    }
                }
            }

            // Inhomogeneous unknowns: −c_inh[i][j] * nn^i * delta(l,j).
            if !opts.homogeneous {
                for i in 0..=opts.inhomo_idx_deg {
                    if l <= opts.inhomo_var_deg {
                        let val =
                            -BigRational::from_integer(num_traits::pow::pow(BigInt::from(nn), i));
                        matrix[row][inhomo_col(i, l)] += val;
                    }
                }
            }
        }
    }

    // --- Solve ---
    let solution = linalg::solve_linear_system(&matrix, &rhs)?;

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
            let bv = extract_bivar(
                &solution,
                |i, j| coeff_col(r, d, i, j),
                opts.idx_deg,
                opts.var_deg,
            );
            if !bv.is_zero() {
                terms.push(RecurrenceTerm {
                    offset: r,
                    deriv_order: d,
                    coeff: bv,
                });
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

    Some(Recurrence {
        terms,
        denominator,
        inhomogeneous,
    })
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

#[derive(Copy, Clone)]
enum CodeStyle {
    Mathematica,
    Sage,
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

impl Recurrence {
    fn max_offset(&self) -> usize {
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

        if term.coeff.is_one() {
            pref
        } else if coeff == "-1" {
            format!("-{pref}")
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

        if term.coeff.is_one() {
            pref
        } else if coeff == "-1" {
            format!("-{pref}")
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
            let pref = fmt_poly_ref_latex(term.offset, term.deriv_order);

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
            let pref = fmt_poly_ref(term.offset, term.deriv_order);

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

// ---------------------------------------------------------------------------
// Parameter counting
// ---------------------------------------------------------------------------

/// Count the total number of unknowns for a given set of recurrence options.
pub fn count_unknowns(opts: &RecurrenceOptions) -> usize {
    let num_denom_vars = (opts.denom_idx_deg + 1) * (opts.denom_var_deg + 1) - 1;
    let vars_per_coeff = (opts.idx_deg + 1) * (opts.var_deg + 1);
    let num_coeff_vars = opts.rec_len * (opts.diff_deg + 1) * vars_per_coeff;
    let num_inhomo_vars = if opts.homogeneous {
        0
    } else {
        (opts.inhomo_idx_deg + 1) * (opts.inhomo_var_deg + 1)
    };
    num_denom_vars + num_coeff_vars + num_inhomo_vars
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
    /// Require every offset 1..=rec_len to appear with a non-zero coefficient.
    pub require_all_offsets: bool,
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
    /// Maximum denom_var_deg when try_denominator is true.
    pub max_denom_var_deg: usize,
    /// Maximum denom_idx_deg when try_denominator is true.
    pub max_denom_idx_deg: usize,
    /// Minimum surplus: equations - unknowns must be >= this.
    pub min_margin: usize,
    /// Print each candidate tried to stderr.
    pub verbose: bool,
}

impl Default for AdaptiveSearchOptions {
    fn default() -> Self {
        Self {
            skip_prefix: 0,
            require_all_offsets: false,
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
            max_denom_var_deg: 2,
            max_denom_idx_deg: 2,
            min_margin: 1,
            verbose: false,
        }
    }
}

fn recurrence_uses_all_offsets(rec: &Recurrence, rec_len: usize) -> bool {
    (1..=rec_len).all(|offset| rec.terms.iter().any(|term| term.offset == offset))
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
    /// Number of equations in the winning system.
    pub num_equations: usize,
    /// Number of candidates actually solved (not just counted).
    pub candidates_tried: usize,
}

/// Generate candidate parameter sets, sorted by (unknowns, diff_deg, rec_len, idx_deg, var_deg).
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
                    });

                    if search.try_inhomogeneous {
                        for inhomo_idx_deg in search.min_inhomo_idx_deg..=search.max_inhomo_idx_deg
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
                                });
                            }
                        }
                    }
                }
            }
        }
    }

    candidates.sort_by_key(|opts| {
        (
            count_unknowns(opts),
            opts.diff_deg,
            opts.rec_len,
            opts.idx_deg,
            opts.var_deg,
        )
    });
    candidates
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

    // First pass: use the options as given.
    let candidates = generate_candidates(m, search);
    let mut tried = 0;

    for opts in &candidates {
        let unknowns = count_unknowns(opts);
        let equations = count_equations_rational(polys, opts);

        if equations < unknowns + search.min_margin {
            continue;
        }

        tried += 1;

        if search.verbose {
            eprintln!(
                "  try #{tried}: rec_len={} var_deg={} idx_deg={} diff_deg={} \
                 denom=({},{}) homog={} \
                 inhomo=({},{}) \
                 (unknowns={unknowns}, equations={equations}, margin={})",
                opts.rec_len,
                opts.var_deg,
                opts.idx_deg,
                opts.diff_deg,
                opts.denom_var_deg,
                opts.denom_idx_deg,
                opts.homogeneous,
                opts.inhomo_var_deg,
                opts.inhomo_idx_deg,
                equations - unknowns,
            );
        }

        if let Some(rec) = find_polynomial_recurrence_rational(polys, opts) {
            if search.require_all_offsets && !recurrence_uses_all_offsets(&rec, opts.rec_len) {
                continue;
            }
            if search.verbose {
                eprintln!("  -> found!");
            }
            return Some(AdaptiveSearchResult {
                recurrence: rec,
                opts: opts.clone(),
                num_unknowns: unknowns,
                num_equations: equations,
                candidates_tried: tried,
            });
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
        for opts in &rational_candidates {
            // Skip candidates without a denominator (already tried above).
            if opts.denom_var_deg == 0 && opts.denom_idx_deg == 0 {
                continue;
            }

            let unknowns = count_unknowns(opts);
            let equations = count_equations_rational(polys, opts);

            if equations < unknowns + search.min_margin {
                continue;
            }

            tried += 1;

            if search.verbose {
                eprintln!(
                    "  try #{tried} (rational): rec_len={} var_deg={} idx_deg={} diff_deg={} \
                     denom=({},{}) \
                     (unknowns={unknowns}, equations={equations}, margin={})",
                    opts.rec_len,
                    opts.var_deg,
                    opts.idx_deg,
                    opts.diff_deg,
                    opts.denom_var_deg,
                    opts.denom_idx_deg,
                    equations - unknowns,
                );
            }

            if let Some(rec) = find_polynomial_recurrence_rational(polys, opts) {
                if search.require_all_offsets && !recurrence_uses_all_offsets(&rec, opts.rec_len) {
                    continue;
                }
                if search.verbose {
                    eprintln!("  -> found (rational)!");
                }
                return Some(AdaptiveSearchResult {
                    recurrence: rec,
                    opts: opts.clone(),
                    num_unknowns: unknowns,
                    num_equations: equations,
                    candidates_tried: tried,
                });
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
            require_all_offsets: true,
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
        // With 6 polys, the adaptive search finds a 3-term non-derivative recurrence
        // (6 unknowns) before the 1-term derivative recurrence (12 unknowns).
        // Both are valid; the search correctly picks the simpler one.
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![1, 1],
            vec![1, 4, 1],
            vec![1, 11, 11, 1],
            vec![1, 26, 66, 26, 1],
        ];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        assert_eq!(result.opts.diff_deg, 0);
        assert_eq!(result.opts.rec_len, 3);

        // With max_diff_deg=0 disabled, force derivative search only:
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
    fn adaptive_short_sequence() {
        // m=3 with constant polys: too few equations for any 2-term recurrence.
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![1], vec![2]];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default());
        assert!(result.is_none());
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
        };
        // m=7, rec_len=2, num_nn=5, max_poly_deg=0, max_j=0, eqs_per_nn=1
        assert_eq!(count_equations(&polys, &opts), 5);
    }
}
