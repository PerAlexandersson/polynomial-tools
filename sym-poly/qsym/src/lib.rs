//! Quasisymmetric functions with generic coefficients.
//!
//! This crate provides [`QSymFunction<C>`] supporting four bases indexed
//! by compositions:
//!
//! - **M_α** (monomial) and **F_α** (Gessel fundamental) — integer coefficients
//! - **Ψ_α** and **Φ_α** (type 1 and type 2 power sums) — rational coefficients
//!
//! The power sum bases follow Ballantine--Daugherty--Hicks--Mason--Niese,
//! *Quasisymmetric Power Sums*, JCTA 2020,
//! <https://doi.org/10.1016/j.jcta.2020.105273>.
//!
//! All pairwise basis conversions are supported, along with the omega
//! involution (ω² = id, ω(Ψ_α) = (-1)^{n-ℓ(α)} Ψ_{α^r}).
//!
//! Additional features:
//! - [`p_partition`]: Stanley's (P,w)-partition generating functions (Ψ̃-positive for naturally labeled posets)
//! - [`chromatic_qsym`]: chromatic quasisymmetric functions
//! - [`sym_qsym`]: maps between Sym and QSym
//! - [`power_sum`]: normalized Ψ̃ = Ψ/z and Φ̃ = Φ/z utilities

pub mod basis;
pub mod qsym_function;
pub mod transition;
pub mod power_sum;
pub mod sym_qsym;
pub mod p_partition;
pub mod chromatic_qsym;

pub use basis::QSymBasis;
pub use qsym_function::QSymFunction;
pub use sym_qsym::{sym_to_qsym, qsym_to_sym};
pub use p_partition::p_partition_generating_function;
pub use chromatic_qsym::chromatic_qsym;
