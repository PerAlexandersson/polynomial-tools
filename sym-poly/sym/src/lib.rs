//! Symmetric functions in the six classical bases with generic coefficients.
//!
//! This crate provides [`SymmetricFunction<C>`] supporting the monomial, elementary,
//! complete homogeneous, power sum, Schur, and forgotten bases. All 30 basis-pair
//! conversions are implemented via cached transition matrices.
//!
//! Built on [`sym_poly_core`] for the Ring trait, Partition type, and matrix utilities.

pub mod basis;
pub mod chromatic;
pub mod frobenius;
pub mod kostka;
pub mod llt;
pub mod shifted_lr;
pub mod symmetric_function;
pub mod transition;

pub use basis::Basis;
pub use chromatic::{chromatic_symmetric, first_bad_edge_symmetric};
pub use frobenius::{frobenius_from_character_values, frobenius_from_trace_matrices};
pub use llt::{graph_llt_symmetric, unicellular_llt, unit_interval_edges};
pub use shifted_lr::{
    shifted_littlewood_richardson_coefficient, shifted_littlewood_richardson_stats,
    shifted_schur_evaluation, ShiftedLrError, ShiftedLrStats,
};
pub use symmetric_function::SymmetricFunction;
