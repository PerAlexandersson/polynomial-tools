//! Dense univariate polynomial toolkit for combinatorial research.
//!
//! Provides [`Polynomial<C>`] with generic coefficients, plus specialized routines
//! for real-rootedness, interlacing, gamma-positivity, resultants, Ehrhart theory,
//! recurrence search, and standard polynomial sequences — all with exact arithmetic.
//!
//! # Quick start
//!
//! Most property-checking functions accept `&[i64]` coefficient vectors in ascending
//! degree order (`coeffs[i]` = coefficient of t^i):
//!
//! ```
//! use polynomial_tools::*;
//!
//! // Eulerian polynomial A_4(t) = 1 + 11t + 11t^2 + t^3
//! assert!(is_real_rooted(&[1, 11, 11, 1]));     // Bézout matrix (default, fast)
//! assert!(is_palindromic(&[1, 11, 11, 1]));
//! assert!(is_gamma_positive(&[1, 11, 11, 1]));   // gamma = [1, 8]
//!
//! // Strict interlacing (degree diff = 1, Bézout-based)
//! assert_eq!(check_interlacing(&[8, -6, 1], &[-15, 23, -9, 1]), Some(true));
//! ```
//!
//! For polynomial arithmetic, use [`Polynomial<C>`]:
//!
//! ```
//! use polynomial_tools::Polynomial;
//!
//! let p = Polynomial::<i64>::from_i64_coeffs(&[1, 1]); // 1 + t
//! let q = p.clone() * p.clone(); // 1 + 2t + t^2
//! assert_eq!(format!("{}", q), "1 + 2t + t^2");
//! assert!(q.is_palindromic());
//! ```
//!
//! # Modules
//!
//! - [`polynomial`] — `Polynomial<C>` with `CoeffRing`/`FieldRing` traits, arithmetic,
//!   derivative, evaluate, shift, reverse, dilate, GCD, division, Lagrange interpolation
//! - [`linalg`] — Exact linear algebra over ℚ: Gaussian elimination, positive
//!   definiteness/semi-definiteness, determinants, linear system solving,
//!   total non-negativity via Neville elimination
//! - [`real_rootedness`] — Bézout matrix (default) and Sturm chain real-rootedness,
//!   strict/weak interlacing (including same-degree via Cauchy bound reduction),
//!   log-concavity, ultra-log-concavity, palindromic check, gamma-positivity,
//!   resultant, discriminant, Ehrhart ↔ h*-vector, display utilities
//! - [`sturm`] — Sturm chains for exact root isolation (used internally)
//! - [`recurrence`] — Adaptive recurrence search for polynomial sequences
//! - [`sequences`] — Standard sequences: Eulerian, Narayana, type B Eulerian,
//!   Chebyshev T/U, Hermite
//! - [`parse`] — Flexible polynomial parsing (comma/space-separated, bracketed,
//!   expanded polynomial notation)

pub mod polynomial;
pub mod linalg;
pub mod sturm;

pub use linalg::{
    check_total_positivity,
    is_totally_nonnegative,
    check_tnn_neville,
    check_tnn_neville_bigint,
    is_tnn,
    is_positive_definite,
    is_positive_semidefinite,
    determinant,
};

pub mod real_rootedness;
pub mod recurrence;
pub mod sequences;
pub mod parse;

pub use polynomial::{Polynomial, CoeffRing, FieldRing};
pub use real_rootedness::{
    // Default (Bézout-based) methods
    is_real_rooted,
    check_interlacing,
    check_weak_interlacing,
    // Sturm-chain methods (slower, but can isolate roots)
    is_real_rooted_sturm,
    check_interlacing_sturm,
    real_roots,
    // Bézout matrix directly
    bezout_matrix,
    // Concavity and symmetry
    is_log_concave,
    is_ultra_log_concave,
    is_palindromic,
    is_gamma_positive,
    gamma_coefficients,
    // Resultant and discriminant
    resultant,
    discriminant,
    sylvester_matrix,
    // Ehrhart polynomial <-> h*-vector
    hstar_to_ehrhart,
    ehrhart_to_hstar,
    ehrhart_to_hstar_with_denom,
    // Display
    format_poly,
    format_poly_var,
};
pub use real_rootedness::{
    is_real_rooted_bezout,
    check_interlacing_bezout,
    check_weak_interlacing_bezout,
};
pub use parse::{parse_polynomial, parse_polynomials};
