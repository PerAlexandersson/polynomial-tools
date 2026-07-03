//! Lah symmetric functions.
//!
//! The implementation follows Pétréolle--Sokal's forest definition.  An
//! unordered forest of increasing ordered trees is represented first by a
//! parent assignment: each non-root vertex has a parent with smaller label.
//! For each parent assignment, the ordered-children contribution is the product
//! of factorials of the outdegrees.

use std::collections::BTreeMap;

use sym_poly_core::Partition;

use crate::{Basis, SymmetricFunction};

/// Expansion of the Lah symmetric function `L_{n,k}^{(infty)+}` in the
/// elementary basis.
pub fn lah_symmetric_elementary(n: usize, k: usize) -> SymmetricFunction<i64> {
    assert!(k >= 1, "k must be positive");
    assert!(k <= n, "k must be at most n");

    let mut child_counts = vec![0usize; n];
    let mut terms: BTreeMap<Partition, i64> = BTreeMap::new();
    enumerate_parent_assignments(2, n, k, 1, &mut child_counts, &mut terms);

    SymmetricFunction::from_terms(Basis::Elementary, terms)
}

/// Expansion of the Lah symmetric function `L_{n,k}^{(infty)+}` in the
/// monomial basis.
pub fn lah_symmetric_monomial(n: usize, k: usize) -> SymmetricFunction<i64> {
    lah_symmetric_elementary(n, k).to_monomial_basis()
}

fn enumerate_parent_assignments(
    vertex: usize,
    n: usize,
    target_roots: usize,
    roots_so_far: usize,
    child_counts: &mut [usize],
    terms: &mut BTreeMap<Partition, i64>,
) {
    if vertex > n {
        if roots_so_far == target_roots {
            let mut parts: Vec<u32> = child_counts
                .iter()
                .copied()
                .filter(|&count| count > 0)
                .map(|count| u32::try_from(count).expect("outdegree does not fit in u32"))
                .collect();
            parts.sort_unstable_by(|a, b| b.cmp(a));
            let partition = Partition::from_sorted(parts);
            let weight: i64 = child_counts.iter().map(|&count| factorial(count)).product();
            *terms.entry(partition).or_insert(0) += weight;
        }
        return;
    }

    if roots_so_far < target_roots {
        enumerate_parent_assignments(
            vertex + 1,
            n,
            target_roots,
            roots_so_far + 1,
            child_counts,
            terms,
        );
    }

    for parent in 1..vertex {
        child_counts[parent - 1] += 1;
        enumerate_parent_assignments(
            vertex + 1,
            n,
            target_roots,
            roots_so_far,
            child_counts,
            terms,
        );
        child_counts[parent - 1] -= 1;
    }
}

fn factorial(n: usize) -> i64 {
    (1..=n)
        .map(|i| i64::try_from(i).expect("factorial input too large"))
        .product()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn partition(parts: &[u32]) -> Partition {
        Partition::new(parts.to_vec())
    }

    #[test]
    fn test_lah_n4_elementary_expansion() {
        let l41 = lah_symmetric_elementary(4, 1);
        assert_eq!(l41.coefficient(&partition(&[3])), 6);
        assert_eq!(l41.coefficient(&partition(&[2, 1])), 8);
        assert_eq!(l41.coefficient(&partition(&[1, 1, 1])), 1);
        assert_eq!(l41.terms().len(), 3);

        let l42 = lah_symmetric_elementary(4, 2);
        assert_eq!(l42.coefficient(&partition(&[2])), 8);
        assert_eq!(l42.coefficient(&partition(&[1, 1])), 7);
        assert_eq!(l42.terms().len(), 2);

        let l43 = lah_symmetric_elementary(4, 3);
        assert_eq!(l43.coefficient(&partition(&[1])), 6);
        assert_eq!(l43.terms().len(), 1);

        let l44 = lah_symmetric_elementary(4, 4);
        assert_eq!(l44.coefficient(&Partition::empty()), 1);
        assert_eq!(l44.terms().len(), 1);
    }

    #[test]
    fn test_lah_n4_monomial_expansion() {
        let l41 = lah_symmetric_monomial(4, 1);
        assert_eq!(l41.coefficient(&partition(&[3])), 1);
        assert_eq!(l41.coefficient(&partition(&[2, 1])), 11);
        assert_eq!(l41.coefficient(&partition(&[1, 1, 1])), 36);

        let l42 = lah_symmetric_monomial(4, 2);
        assert_eq!(l42.coefficient(&partition(&[2])), 7);
        assert_eq!(l42.coefficient(&partition(&[1, 1])), 22);
    }

    #[test]
    fn test_lah_forest_counts_match_closed_formula_small() {
        let expected = [
            vec![1],
            vec![1, 1],
            vec![3, 3, 1],
            vec![15, 15, 6, 1],
            vec![105, 105, 45, 10, 1],
        ];

        for n in 1..=5 {
            for k in 1..=n {
                let total: i64 = lah_symmetric_elementary(n, k).terms().values().sum();
                assert_eq!(total, expected[n - 1][k - 1]);
            }
        }
    }
}
