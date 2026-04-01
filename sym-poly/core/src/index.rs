//! Generic basis indexing trait for graded algebras.
//!
//! Partitions index symmetric functions (Sym), compositions index
//! quasisymmetric functions (QSym), and set partitions will index
//! noncommutative symmetric functions (NCSym) in the future.

use std::fmt;
use std::hash::Hash;

use crate::composition::Composition;
use crate::partition::Partition;

/// An index type for a graded algebra of symmetric-like functions.
///
/// Each basis element is labeled by an index of this type, and the
/// algebra is graded by the `degree` function. The trait provides
/// enough structure for transition matrix caching and formal sum operations.
pub trait BasisIndex: Clone + Eq + Ord + Hash + fmt::Debug + fmt::Display + Sized {
    /// The degree (grading value) of this index.
    fn degree(&self) -> u32;

    /// Whether this is the index of degree 0 (the "empty" index).
    fn is_empty(&self) -> bool;

    /// The empty (degree-0) index.
    fn empty() -> Self;

    /// All indices of a given degree, in a canonical order.
    ///
    /// For partitions: reverse lexicographic order (matching `Partition::all_of_size`).
    /// For compositions: lexicographic order (matching `Composition::integer_compositions`).
    fn all_of_degree(n: u32) -> Vec<Self>;
}

// ---------------------------------------------------------------------------
// Partition implements BasisIndex
// ---------------------------------------------------------------------------

impl BasisIndex for Partition {
    fn degree(&self) -> u32 {
        self.size()
    }

    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn empty() -> Self {
        Partition::empty()
    }

    fn all_of_degree(n: u32) -> Vec<Self> {
        Partition::all_of_size(n)
    }
}

// ---------------------------------------------------------------------------
// Composition implements BasisIndex
// ---------------------------------------------------------------------------

impl BasisIndex for Composition {
    fn degree(&self) -> u32 {
        self.size()
    }

    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn empty() -> Self {
        Composition::empty()
    }

    fn all_of_degree(n: u32) -> Vec<Self> {
        Composition::integer_compositions(n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_partition_basis_index() {
        let p = Partition::new(vec![3, 2, 1]);
        assert_eq!(p.degree(), 6);
        assert!(!p.is_empty());
        assert!(Partition::empty().is_empty());

        let all = Partition::all_of_degree(4);
        assert_eq!(all.len(), 5); // p(4) = 5
    }

    #[test]
    fn test_composition_basis_index() {
        let c = Composition::new(vec![2, 1, 1]);
        assert_eq!(c.degree(), 4);
        assert!(!c.is_empty());
        assert!(Composition::empty().is_empty());

        let all = Composition::all_of_degree(4);
        assert_eq!(all.len(), 8); // 2^(4-1) = 8
    }
}
