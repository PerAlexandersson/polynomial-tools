//! Symmetric functions in the six classical bases with generic coefficients.
//!
//! This crate provides [`SymmetricFunction<C>`] supporting the monomial, elementary,
//! complete homogeneous, power sum, Schur, and forgotten bases. All 30 basis-pair
//! conversions are implemented via cached transition matrices.
//!
//! Built on [`sym_poly_core`] for the Ring trait, Partition type, and matrix utilities.

pub mod basis;
pub mod symmetric_function;
pub mod kostka;
pub mod transition;

pub use basis::Basis;
pub use symmetric_function::SymmetricFunction;
