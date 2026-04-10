//! Dense univariate polynomial toolkit for combinatorial research.
//!
//! Provides [`Polynomial<C>`] with generic coefficients, plus specialized routines
//! for real-rootedness, interlacing, gamma-positivity, Stapledon decomposition, resultants, Ehrhart theory,
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
//!   Stapledon decomposition, resultant, discriminant, Ehrhart ↔ h*-vector,
//!   display utilities
//! - [`sturm`] — Sturm chains for exact root isolation (used internally)
//! - [`recurrence`] — Adaptive recurrence search for polynomial sequences
//! - [`sequences`] — Standard sequences: Eulerian, Narayana, type B Eulerian,
//!   Chebyshev T/U, Hermite
//! - [`parse`] — Flexible polynomial parsing (comma/space-separated, bracketed,
//!   expanded polynomial notation)

pub mod linalg;
pub mod polynomial;
pub mod sturm;
pub mod vec_poly;

pub use linalg::{
    check_tnn_neville, check_tnn_neville_bigint, check_total_positivity, determinant,
    is_positive_definite, is_positive_semidefinite, is_tnn, is_totally_nonnegative,
};

pub mod parse;
pub mod real_rootedness;
pub mod recurrence;
pub mod sequences;

pub use parse::{parse_polynomial, parse_polynomials};
pub use polynomial::{CoeffRing, FieldRing, Polynomial};
pub use real_rootedness::{
    // Bézout matrix directly
    bezout_matrix,
    check_interlacing,
    check_interlacing_sturm,
    check_weak_interlacing,
    discriminant,
    ehrhart_to_hstar,
    ehrhart_to_hstar_with_denom,
    // Display
    format_poly,
    format_poly_var,
    gamma_coefficients,
    // Ehrhart polynomial <-> h*-vector
    hstar_to_ehrhart,
    is_gamma_positive,
    // Concavity and symmetry
    is_log_concave,
    is_palindromic,
    // Default (Bézout-based) methods
    is_real_rooted,
    // Sturm-chain methods (slower, but can isolate roots)
    is_real_rooted_sturm,
    is_ultra_log_concave,
    real_roots,
    // Resultant and discriminant
    resultant,
    stapledon_decomposition,
    sylvester_matrix,
};
pub use real_rootedness::{
    check_interlacing_bezout, check_weak_interlacing_bezout, is_real_rooted_bezout,
};
