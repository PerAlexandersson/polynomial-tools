//! Shared foundation for symmetric polynomial algebras.
//!
//! This crate provides the common building blocks used by `sym-poly-sym` (symmetric functions),
//! `sym-poly-qsym` (quasisymmetric functions), and `sym-poly-multipoly` (multivariate polynomials):
//!
//! - [`Ring`] trait for generic coefficients (i64, BigInt, rationals)
//! - [`Partition`] and [`Composition`] types for basis indexing
//! - [`BasisIndex`] trait abstracting over index types
//! - [`FormalSum`] for generic linear combinations of basis elements
//! - [`matrix`] utilities for transition matrix computation
//! - [`TransitionCache`] for per-algebra cached transition matrices

pub mod ring;
pub mod partition;
pub mod composition;
pub mod index;
pub mod formal_sum;
pub mod matrix;
pub mod transition_cache;

pub use ring::Ring;
pub use partition::Partition;
pub use composition::Composition;
pub use index::BasisIndex;
pub use formal_sum::FormalSum;
pub use transition_cache::TransitionCache;
