//! Multivariate polynomials with divided difference operators and nonsymmetric bases.
//!
//! This crate provides:
//! - [`MultiPoly<C>`] — sparse multivariate polynomials
//! - [`MultiPolyFunction<C>`] — polynomials in nonsymmetric bases (key, atom, slide)
//! - [`MultiPolyBasis`] — basis enum (Monomial, Key, Atom, MonSlide, FundSlide)
//! - [`operators`] — simple, Demazure, and t-deformed operators (∂_i, π_i, θ_i)
//! - [`key_polynomial`] — Demazure characters via π operators
//! - [`atom_polynomial`] — Demazure atoms via θ operators
//! - [`kohnert`] — Kohnert diagrams and Assaf Yamanouchi tests
//! - [`nonsymmetric_macdonald`] — operator-side `q = 0` Macdonald / Hall-Littlewood
//! - [`schubert_polynomial`] — Schubert polynomials via divided differences
//! - [`slide_polynomial`] — monomial and fundamental slide polynomials

pub mod atom_polynomial;
pub mod basis;
pub mod key_polynomial;
pub mod kohnert;
pub mod lorentzian;
pub mod multipoly;
pub mod multipoly_function;
pub mod nonsymmetric_macdonald;
pub mod operators;
pub mod schubert_polynomial;
pub mod slide_polynomial;
pub mod transition;

pub use atom_polynomial::{atom_polynomial, t_atom_polynomial};
pub use basis::MultiPolyBasis;
pub use key_polynomial::{key_polynomial, t_key_polynomial};
pub use kohnert::{
    canonical_labeling, cells_in_col, column_pairing, diagram_from_labeling, diagram_weight,
    format_diagram, is_yamanouchi, kohnert_diagrams, kohnert_moves, label_pairing, max_col,
    rectify_labeled, rectify_labeled_column_star, rothe_diagram, sorted_rows_in_col,
    yamanouchi_diagrams, Cell, Diagram, Labeling,
};
pub use lorentzian::{
    is_lorentzian, is_lorentzian_bool, is_m_convex, is_normalized_lorentzian,
    is_normalized_lorentzian_bool, is_strictly_lorentzian, is_strictly_normalized_lorentzian,
    support_is_m_convex, LorentzianResult,
};
pub use multipoly::MultiPoly;
pub use multipoly_function::MultiPolyFunction;
pub use nonsymmetric_macdonald::{nonsymmetric_hall_littlewood, nonsymmetric_macdonald_q0};
pub use operators::{
    partial_i, partial_word, pi_i, pi_word, theta_i, theta_word, tpi_i, tpi_word, ttheta_i,
    ttheta_word,
};
pub use schubert_polynomial::{
    schubert_polynomial, schubert_to_atom, schubert_to_fund_slide, schubert_to_key,
    schubert_to_monomial,
};
pub use slide_polynomial::{fundamental_slide_polynomial, monomial_slide_polynomial};
