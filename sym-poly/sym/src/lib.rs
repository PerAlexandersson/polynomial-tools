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
pub mod hessenberg_gkm;
pub mod kostka;
pub mod llt;
pub mod macdonald;
pub mod shifted_lr;
pub mod symmetric_function;
pub mod transition;
pub mod twin_gkm;

pub use basis::Basis;
pub use chromatic::{
    chromatic_symmetric, circular_area_dot_frobenius_target, first_bad_edge_symmetric,
    hessenberg_area_dot_frobenius_target, q_chromatic_symmetric_with_ascent_edges,
};
pub use frobenius::{
    frobenius_from_character_values, frobenius_from_trace_matrices,
    graded_frobenius_from_character_values, graded_frobenius_from_trace_matrices,
    multigraded_frobenius_from_character_values, multigraded_frobenius_from_trace_matrices,
};
pub use hessenberg_gkm::{
    affine_shadow_circular_gkm_hilbert, hessenberg_gkm_dot_action_matrices,
    hessenberg_gkm_dot_character_values_by_degree,
    hessenberg_gkm_dot_character_values_by_degree_packed,
    hessenberg_gkm_dot_character_values_packed_crt,
    hessenberg_gkm_dot_character_values_packed_mod_prime, hessenberg_gkm_dot_frobenius,
    hessenberg_gkm_dot_frobenius_packed, naive_circular_gkm_dot_character_values_by_degree,
    naive_circular_gkm_dot_frobenius,
};
pub use llt::{
    circular_unicellular_llt, circular_unicellular_llt_character_values_by_degree,
    circular_unicellular_llt_frobenius_target, circular_unicellular_llt_q_plus_one,
    circular_unicellular_llt_q_plus_one_e_expansion,
    circular_unicellular_llt_q_plus_one_is_e_positive, directed_graph_llt_symmetric,
    graph_llt_symmetric, unicellular_llt, unicellular_llt_character_values_by_degree,
    unicellular_llt_frobenius_target, unit_interval_edges,
};
pub use macdonald::{
    delta_eigenvalue, delta_modified_macdonald, delta_prime_eigenvalue,
    delta_prime_modified_macdonald, macdonald_b_alphabet, macdonald_b_eigenvalue, nabla_eigenvalue,
    nabla_modified_macdonald, qt_coefficient, qt_constant, qt_monomial, ModifiedMacdonaldExpansion,
    QtPolynomial,
};
pub use shifted_lr::{
    shifted_littlewood_richardson_coefficient, shifted_littlewood_richardson_stats,
    shifted_schur_evaluation, ShiftedLrError, ShiftedLrStats,
};
pub use symmetric_function::SymmetricFunction;
pub use twin_gkm::{
    twin_gkm_dagger_action_matrices, twin_gkm_dagger_character_values_by_degree,
    twin_gkm_dagger_frobenius,
};
