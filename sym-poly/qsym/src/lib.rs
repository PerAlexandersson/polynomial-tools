//! Quasisymmetric functions with generic coefficients.
//!
//! This crate provides [`QSymFunction<C>`] supporting six bases indexed
//! by compositions:
//!
//! - **M_α** (monomial) and **F_α** (Gessel fundamental) — integer coefficients
//! - **S_α** (quasisymmetric Schur) and **S*_α** (dual immaculate)
//! - **Ψ_α** and **Φ_α** (type 1 and type 2 power sums) — rational coefficients
//!
//! The power sum bases follow Ballantine--Daugherty--Hicks--Mason--Niese,
//! *Quasisymmetric Power Sums*, JCTA 2020,
//! <https://doi.org/10.1016/j.jcta.2020.105273>.
//!
//! All pairwise basis conversions are supported, along with the omega
//! involution (ω² = id, ω(Ψ_α) = (-1)^{n-ℓ(α)} Ψ_{α^r}) and the `psi`
//! involution on the fundamental basis.
//!
//! Additional features:
//! - [`p_partition`]: Stanley's (P,w)-partition generating functions (Ψ̃-positive for naturally labeled posets)
//! - [`chromatic_qsym`]: chromatic quasisymmetric functions, including an
//!   asc-weighted Shareshian--Wachs style refinement
//! - [`sym_qsym`]: maps between Sym and QSym
//! - [`power_sum`]: normalized Ψ̃ = Ψ/z and Φ̃ = Φ/z utilities

pub mod basis;
pub mod chromatic_qsym;
pub mod p_partition;
pub mod peak;
pub mod power_sum;
pub mod qsym_function;
pub mod schur_qsym;
pub mod sym_qsym;
pub mod transition;

pub use basis::QSymBasis;
pub use chromatic_qsym::{
    chromatic_qsym, chromatic_qsym_asc, circular_coloring_qsym_asc, coloring_qsym_asc,
    coloring_qsym_asc_with_ascent_edges,
};
pub use p_partition::{
    p_partition_generating_function, p_partition_generating_function_with_labels,
    p_partition_linear_extensions, p_partition_linear_extensions_with_labels,
    strict_p_partition_generating_function, PPartitionLinearExtension,
};
pub use peak::{is_peak_set, peak_quasisymmetric};
pub use qsym_function::QSymFunction;
pub use schur_qsym::{
    composition_to_descent_set, descent_set_to_composition, dual_immaculate,
    dual_immaculate_monomial, fundamental_slide, qsym_schur, row_strict_dual_immaculate,
    row_strict_dual_immaculate_fundamental, row_strict_qsym_schur, row_strict_young_qsym_schur,
    young_qsym_schur, CompositionTableau, ImmaculateTableau, ReverseCompositionTableau,
};
pub use sym_qsym::{qsym_to_sym, sym_to_qsym, symmetric_qsym_to_sym};
