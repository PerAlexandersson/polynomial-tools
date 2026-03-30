//! Multivariate polynomials with divided difference operators and key polynomials.
//!
//! This crate provides:
//! - [`MultiPoly<C>`] — sparse multivariate polynomials
//! - [`operators`] — simple and isobaric divided difference operators (∂_i, π_i)
//! - [`key_polynomial`] — Demazure characters via π operators

pub mod multipoly;
pub mod operators;
pub mod key_polynomial;

pub use multipoly::MultiPoly;
pub use operators::{partial_i, pi_i};
pub use key_polynomial::key_polynomial;
