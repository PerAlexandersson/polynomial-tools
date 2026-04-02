//! Shared foundation for symmetric polynomial algebras.
//!
//! This crate provides the common building blocks used by `sym-poly-sym` (symmetric functions),
//! `sym-poly-qsym` (quasisymmetric functions), and `sym-poly-multipoly` (multivariate polynomials):
//!
//! - [`Ring`] trait for generic coefficients (i64, BigInt, rationals)
//! - [`Partition`] and [`Composition`] types for basis indexing
//! - [`UnivariatePolynomial`] for polynomial-valued coefficients
//! - [`Tableau`] for basic Young tableau combinatorics
//! - [`BasisIndex`] trait abstracting over index types
//! - [`FormalSum`] for generic linear combinations of basis elements
//! - [`matrix`] utilities for transition matrix computation
//! - [`TransitionCache`] for per-algebra cached transition matrices

pub mod composition {
    pub use combinatoric_core::composition::*;
}
pub mod formal_sum;
pub mod index;
pub mod matrix;
pub mod partition {
    pub use combinatoric_core::partition::*;
}
pub mod polynomial;
pub mod ring {
    pub use combinatoric_core::ring::*;
}
pub mod tableau;
pub mod transition_cache;

pub use combinatoric_core::{Composition, Partition, Ring};
pub use formal_sum::FormalSum;
pub use index::BasisIndex;
pub use polynomial::UnivariatePolynomial;
pub use tableau::Tableau;
pub use transition_cache::TransitionCache;
