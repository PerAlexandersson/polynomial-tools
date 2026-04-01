//! Shared library for algebraic combinatorics: partitions, compositions,
//! symmetric functions, graphs, posets, chromatic functions, and key polynomials.
//!
//! Core types (Ring, Partition, Composition, SymmetricFunction) are re-exported
//! from the `sym-poly-core` and `sym-poly-sym` crates. Graph, Poset, chromatic,
//! and key polynomial modules remain here.

// Re-export foundation types from sym-poly-core
pub mod ring {
    pub use sym_poly_core::ring::*;
}
pub mod partition {
    pub use sym_poly_core::partition::*;
}
pub mod composition {
    pub use sym_poly_core::composition::*;
}

// Re-export symmetric function types from sym-poly-sym
pub mod symmetric_function {
    pub use sym_poly_sym::basis::*;
    pub use sym_poly_sym::symmetric_function::*;
}
pub mod kostka {
    pub use sym_poly_sym::kostka::*;
}
pub mod transition {
    pub use sym_poly_sym::transition::*;
}

// Local modules (not moved to sym-poly)
pub mod chromatic;
pub mod graph;
pub mod key_polynomial;
pub mod poset;

// Top-level re-exports for convenience
pub use graph::Graph;
pub use sym_poly_core::Composition;
pub use sym_poly_core::Partition;
pub use sym_poly_core::Ring;
pub use sym_poly_sym::Basis;
pub use sym_poly_sym::SymmetricFunction;
