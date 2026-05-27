//! Multivariate polynomials with divided difference operators and nonsymmetric bases.
//!
//! This crate provides:
//! - [`MultiPoly<C>`] — sparse multivariate polynomials
//! - [`MultiPolyFunction<C>`] — polynomials in nonsymmetric bases (key, atom, slide)
//! - [`MultiPolyBasis`] — basis enum (Monomial, Key, Atom, MonSlide, FundSlide)
//! - [`monomial_order`] — monomial orders and leading terms
//! - [`division`] — multivariate division and normal forms
//! - [`groebner`] — basic Buchberger algorithm for small exact quotients
//! - [`quotient`] — finite standard-monomial quotient bases
//! - [`operators`] — simple, Demazure, and t-deformed operators (∂_i, π_i, θ_i)
//! - [`key_polynomial`] — Demazure characters via π operators
//! - [`atom_polynomial`] — Demazure atoms via θ operators
//! - [`kohnert`] — Kohnert diagrams and Assaf Yamanouchi tests
//! - [`nonsymmetric_macdonald`] — operator-side `q = 0` Macdonald / Hall-Littlewood
//! - [`schubert_polynomial`] — Schubert polynomials via divided differences
//! - [`slide_polynomial`] — monomial and fundamental slide polynomials

pub mod atom_polynomial;
pub mod basis;
pub mod division;
pub mod groebner;
pub mod indexed_variables;
pub mod key_polynomial;
pub mod kohnert;
pub mod lorentzian;
pub mod monomial_order;
pub mod multipoly;
pub mod multipoly_function;
pub mod nonsymmetric_macdonald;
pub mod operators;
pub mod quotient;
pub mod quotient_module;
pub mod schubert_polynomial;
pub mod slide_polynomial;
pub mod symmetric_polynomials;
pub mod transition;

pub use atom_polynomial::{atom_polynomial, t_atom_polynomial};
pub use basis::MultiPolyBasis;
pub use division::{divide_by_polynomials, multiply_by_monomial, normal_form, DivisionResult};
pub use groebner::{
    buchberger_basis, is_groebner_basis, make_monic, reduced_groebner_basis, s_polynomial,
    GroebnerBasis,
};
pub use indexed_variables::{
    ideal_generators_are_invariant_under_index_permutation, ideal_generators_are_sn_invariant,
    is_multidegree_preserving_action_matrix,
    quotient_action_matrices_by_index_permutation_and_multidegree,
    quotient_action_matrices_by_multidegree_and_cycle_type,
    quotient_action_matrix_multidegree_blocks, quotient_basis_multidegrees, IndexedVariables,
};
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
pub use monomial_order::{
    leading_term, monomial_divides, monomial_quotient, LeadingTerm, MonomialOrder,
};
pub use multipoly::MultiPoly;
pub use multipoly_function::MultiPolyFunction;
pub use nonsymmetric_macdonald::{nonsymmetric_hall_littlewood, nonsymmetric_macdonald_q0};
pub use operators::{
    partial_i, partial_word, pi_i, pi_word, theta_i, theta_word, tpi_i, tpi_word, ttheta_i,
    ttheta_word,
};
pub use quotient::{
    is_degree_preserving_action_matrix, normal_form_in_basis, permute_variables, pure_power_bounds,
    quotient_action_matrices_by_permutation_and_degree, quotient_action_matrix_by_permutation,
    quotient_action_matrix_degree_blocks, quotient_basis, quotient_basis_degrees,
    quotient_coordinates, restrict_matrix_to_indices, standard_monomials_from_leading_monomials,
    QuotientBasis,
};
pub use quotient_module::{PolynomialQuotientSnModule, PolynomialQuotientSnModuleError};
pub use schubert_polynomial::{
    schubert_polynomial, schubert_to_atom, schubert_to_fund_slide, schubert_to_key,
    schubert_to_monomial,
};
pub use slide_polynomial::{fundamental_slide_polynomial, monomial_slide_polynomial};
pub use symmetric_polynomials::{elementary_symmetric_generators, elementary_symmetric_polynomial};
