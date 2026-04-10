//! Shared library for foundational combinatorics: partitions, set partitions,
//! compositions, permutations, graphs, posets, and related utilities.

pub mod composition;
pub mod graph;
pub mod key_polynomial;
pub mod partition;
pub mod permutation;
pub mod poset;
pub mod ring;
pub mod set_partition;

// Top-level re-exports for convenience
pub use composition::{Composition, WeakComposition};
pub use graph::Graph;
pub use partition::Partition;
pub use permutation::{
    all_permutations_zero_indexed, compose_permutations, inverse_permutation, longest_permutation,
    next_permutation, permutation_from_simple_transpositions, reduced_word,
};
pub use ring::Ring;
pub use set_partition::{
    ordered_set_partitions, set_partitions, OrderedSetPartition, SetPartition,
};
