//! Shared foundation for symmetric polynomial algebras.
//!
//! This crate provides the common building blocks used by `sym-poly-sym` (symmetric functions),
//! `sym-poly-qsym` (quasisymmetric functions), and `sym-poly-multipoly` (multivariate polynomials):
//!
//! - [`Ring`] trait for generic coefficients (i64, BigInt, rationals)
//! - [`Partition`] and [`Composition`] types for basis indexing
//! - [`UnivariatePolynomial`] for polynomial-valued coefficients
//! - [`Tableau`], [`SkewTableau`], and their lazy standard-tableau iterators
//! - [`BasisIndex`] trait abstracting over index types
//! - [`PrimeField`] for modular computations over finite prime fields
//! - [`FormalSum`] for generic linear combinations of basis elements
//! - [`matrix`] utilities for transition matrix computation
//! - [`linear_algebra`] utilities for exact row reduction and quotient spaces
//! - [`sn_action`] utilities for small exact `S_n`-module computations
//! - [`RationalFunction`] and [`QtRationalFunction`] for `Q(q,t)` coefficient fields
//! - [`TransitionCache`] for per-algebra cached transition matrices

pub mod composition {
    pub use combinatoric_core::composition::*;
}
pub mod crt;
pub mod field;
pub mod finite_field;
pub mod formal_sum;
pub mod index;
pub mod linear_algebra;
pub mod matrix;
pub mod partition {
    pub use combinatoric_core::partition::*;
}
pub mod polynomial;
pub mod rational_function;
pub mod ring {
    pub use combinatoric_core::ring::*;
}
pub mod sn_action;
pub mod ssaf;
pub mod tableau;
pub mod transition_cache;

pub use combinatoric_core::{Composition, Partition, Ring, WeakComposition};
pub use crt::{
    chinese_remainder, chinese_remainder_pair, rational_reconstruction, symmetric_residue,
};
pub use field::Field;
pub use finite_field::PrimeField;
pub use formal_sum::FormalSum;
pub use index::BasisIndex;
pub use polynomial::UnivariatePolynomial;
pub use rational_function::{
    qt_coefficient, qt_constant, qt_monomial, qt_polynomial_constant, qt_rational_monomial,
    qt_unit_monomial, QtPolynomial, QtRationalFunction, RationalFunction,
};
pub use ssaf::Ssaf;
pub use tableau::{SkewTableau, StandardSkewTableauxIter, StandardTableauxIter, Tableau};
pub use transition_cache::TransitionCache;
