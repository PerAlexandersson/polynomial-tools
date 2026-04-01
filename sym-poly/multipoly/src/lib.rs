//! Multivariate polynomials with divided difference operators and nonsymmetric bases.
//!
//! This crate provides:
//! - [`MultiPoly<C>`] — sparse multivariate polynomials
//! - [`operators`] — simple, Demazure, and t-deformed operators (∂_i, π_i, θ_i)
//! - [`key_polynomial`] — Demazure characters via π operators
//! - [`atom_polynomial`] — Demazure atoms via θ operators
//! - [`nonsymmetric_macdonald`] — operator-side `q = 0` Macdonald / Hall-Littlewood
//! - [`schubert_polynomial`] — Schubert polynomials via divided differences
//! - [`slide_polynomial`] — fundamental slide polynomials

pub mod atom_polynomial;
pub mod key_polynomial;
pub mod multipoly;
pub mod nonsymmetric_macdonald;
pub mod operators;
pub mod schubert_polynomial;
pub mod slide_polynomial;

pub use atom_polynomial::{atom_polynomial, t_atom_polynomial};
pub use key_polynomial::{key_polynomial, t_key_polynomial};
pub use multipoly::MultiPoly;
pub use nonsymmetric_macdonald::{nonsymmetric_hall_littlewood, nonsymmetric_macdonald_q0};
pub use operators::{
    partial_i, partial_word, pi_i, pi_word, theta_i, theta_word, tpi_i, tpi_word, ttheta_i,
    ttheta_word,
};
pub use schubert_polynomial::schubert_polynomial;
pub use slide_polynomial::fundamental_slide_polynomial;
