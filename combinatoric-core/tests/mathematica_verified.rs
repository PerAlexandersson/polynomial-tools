//! Integration tests verified against Mathematica's SymmetricFunctions.m and CombinatoricTools.m.
//!
//! All expected values in this file were computed using Mathematica and cross-checked.

use combinatoric_core::composition::Composition;
use combinatoric_core::kostka;
use combinatoric_core::partition::Partition;
use combinatoric_core::ring::Ring;
use combinatoric_core::symmetric_function::{Basis, SymmetricFunction};
use num_rational::Ratio;
use std::collections::BTreeSet;

type Q = Ratio<i64>;

fn p(parts: &[u32]) -> Partition {
    Partition::new(parts.to_vec())
}

fn q(n: i64, d: i64) -> Q {
    Ratio::new(n, d)
}

// ===========================================================================
// Partition utilities — verified against CombinatoricTools.m
// ===========================================================================

#[test]
fn test_conjugate_self_conjugate() {
    // (4,3,2,1) is self-conjugate
    assert_eq!(p(&[4, 3, 2, 1]).conjugate_partition(), p(&[4, 3, 2, 1]));
}

#[test]
fn test_conjugate_various() {
    assert_eq!(p(&[5, 3, 1]).conjugate_partition(), p(&[3, 2, 2, 1, 1]));
    assert_eq!(p(&[4]).conjugate_partition(), p(&[1, 1, 1, 1]));
    assert_eq!(p(&[1, 1, 1, 1]).conjugate_partition(), p(&[4]));
    assert_eq!(Partition::empty().conjugate_partition(), Partition::empty());
}

#[test]
fn test_partition_n_mathematica() {
    // PartitionN[{4,3,2,1}] = 0*4 + 1*3 + 2*2 + 3*1 = 10
    assert_eq!(p(&[4, 3, 2, 1]).partition_n(), 10);
    assert_eq!(p(&[5]).partition_n(), 0);
    assert_eq!(p(&[1, 1, 1]).partition_n(), 3);
    assert_eq!(Partition::empty().partition_n(), 0);
}

#[test]
fn test_hook_lengths_mathematica() {
    // HookLengths[{4,3,2,1}] = {{7,5,3,1},{5,3,1},{3,1},{1}}
    assert_eq!(
        p(&[4, 3, 2, 1]).hook_lengths(),
        vec![vec![7, 5, 3, 1], vec![5, 3, 1], vec![3, 1], vec![1]]
    );
    // HookLengths[{3}] = {{3}}... no, wait. (3) means one row of 3 cells.
    // Hooks: arm + leg + 1. Conj = (1,1,1).
    // (0,0): arm=2, leg=0 → 3. (0,1): arm=1, leg=0 → 2. (0,2): arm=0, leg=0 → 1.
    assert_eq!(p(&[3]).hook_lengths(), vec![vec![3, 2, 1]]);
}

#[test]
fn test_arm_leg() {
    let lam = p(&[4, 3, 2, 1]);
    // (0,0): arm = 4-0-1 = 3, leg = conj[0]-0-1 = 3
    assert_eq!(lam.partition_arm(0, 0), Some(3));
    assert_eq!(lam.partition_leg(0, 0), Some(3));
    // (0,3): arm = 4-3-1 = 0, leg = conj[3]-0-1 = 0
    assert_eq!(lam.partition_arm(0, 3), Some(0));
    assert_eq!(lam.partition_leg(0, 3), Some(0));
    // Out of diagram
    assert_eq!(lam.partition_arm(0, 4), None);
    assert_eq!(lam.partition_arm(4, 0), None);
}

#[test]
fn test_diagram_boxes_mathematica() {
    // DiagramBoxes[{4,3,2,1}] = 10 boxes (Mathematica uses 1-indexed)
    let boxes = p(&[4, 3, 2, 1]).diagram_boxes();
    assert_eq!(boxes.len(), 10);
    assert!(boxes.contains(&(0, 0)));
    assert!(boxes.contains(&(0, 3)));
    assert!(boxes.contains(&(3, 0)));
    assert!(!boxes.contains(&(1, 3)));
}

#[test]
fn test_skew_diagram_boxes() {
    let outer = p(&[3, 2]);
    let inner = p(&[1]);
    let boxes = outer.skew_diagram_boxes(&inner);
    // Row 0: cols 1, 2. Row 1: cols 0, 1.
    assert_eq!(boxes.len(), 4);
    assert!(boxes.contains(&(0, 1)));
    assert!(boxes.contains(&(0, 2)));
    assert!(boxes.contains(&(1, 0)));
    assert!(boxes.contains(&(1, 1)));
}

#[test]
fn test_partition_add_box_mathematica() {
    // PartitionAddBox[{4,3,2,1}] = {{5,3,2,1}, {4,4,2,1}, {4,3,3,1}, {4,3,2,2}, {4,3,2,1,1}}
    let added = p(&[4, 3, 2, 1]).partition_add_box();
    assert_eq!(added.len(), 5);
    assert!(added.contains(&p(&[5, 3, 2, 1])));
    assert!(added.contains(&p(&[4, 4, 2, 1])));
    assert!(added.contains(&p(&[4, 3, 3, 1])));
    assert!(added.contains(&p(&[4, 3, 2, 2])));
    assert!(added.contains(&p(&[4, 3, 2, 1, 1])));
}

#[test]
fn test_partition_remove_box_mathematica() {
    // PartitionRemoveBox[{4,3,2,1}] = {{3,3,2,1}, {4,2,2,1}, {4,3,1,1}, {4,3,2}}
    let removed = p(&[4, 3, 2, 1]).partition_remove_box();
    assert_eq!(removed.len(), 4);
    assert!(removed.contains(&p(&[3, 3, 2, 1])));
    assert!(removed.contains(&p(&[4, 2, 2, 1])));
    assert!(removed.contains(&p(&[4, 3, 1, 1])));
    assert!(removed.contains(&p(&[4, 3, 2])));
}

#[test]
fn test_partition_less_equal() {
    assert!(p(&[2, 1]).partition_less_equal(&p(&[3, 2])));
    assert!(p(&[2, 1]).partition_less_equal(&p(&[2, 1])));
    assert!(!p(&[3, 1]).partition_less_equal(&p(&[2, 2])));
    assert!(Partition::empty().partition_less_equal(&p(&[1])));
}

#[test]
fn test_dominance_order_mathematica() {
    // From Mathematica: verified dominance relations on partitions of 5
    assert!(p(&[4, 1]).dominated_by(&p(&[5])));
    assert!(!p(&[5]).dominated_by(&p(&[4, 1])));
    assert!(p(&[3, 2]).dominated_by(&p(&[4, 1])));
    assert!(!p(&[4, 1]).dominated_by(&p(&[3, 2])));
    assert!(p(&[3, 1, 1]).dominated_by(&p(&[4, 1])));
    assert!(p(&[3, 1, 1]).dominated_by(&p(&[3, 2])));
    assert!(p(&[2, 2, 1]).dominated_by(&p(&[3, 2])));
    assert!(p(&[2, 2, 1]).dominated_by(&p(&[3, 1, 1])));
    assert!(p(&[2, 1, 1, 1]).dominated_by(&p(&[2, 2, 1])));
    assert!(p(&[1, 1, 1, 1, 1]).dominated_by(&p(&[2, 1, 1, 1])));
    // Different sizes: never dominates
    assert!(!p(&[3]).dominated_by(&p(&[2, 2])));
}

#[test]
fn test_strictly_dominated() {
    assert!(p(&[2, 2, 1]).strictly_dominated_by(&p(&[3, 1, 1])));
    assert!(!p(&[3, 1, 1]).strictly_dominated_by(&p(&[3, 1, 1])));
}

#[test]
fn test_partition_part_count() {
    // (3,2,2,1): m_1=1, m_2=2, m_3=1
    let counts = p(&[3, 2, 2, 1]).partition_part_count();
    assert_eq!(counts[1], 1);
    assert_eq!(counts[2], 2);
    assert_eq!(counts[3], 1);
}

#[test]
fn test_z_coefficients_mathematica() {
    // Verified against Mathematica ZCoefficient
    assert_eq!(p(&[1]).z_coefficient(), 1);
    assert_eq!(p(&[2]).z_coefficient(), 2);
    assert_eq!(p(&[1, 1]).z_coefficient(), 2);
    assert_eq!(p(&[3]).z_coefficient(), 3);
    assert_eq!(p(&[2, 1]).z_coefficient(), 2);
    assert_eq!(p(&[1, 1, 1]).z_coefficient(), 6);
    assert_eq!(p(&[4]).z_coefficient(), 4);
    assert_eq!(p(&[3, 1]).z_coefficient(), 3);
    assert_eq!(p(&[2, 2]).z_coefficient(), 8);
    assert_eq!(p(&[2, 1, 1]).z_coefficient(), 4);
    assert_eq!(p(&[1, 1, 1, 1]).z_coefficient(), 24);
    assert_eq!(p(&[5]).z_coefficient(), 5);
    assert_eq!(p(&[4, 1]).z_coefficient(), 4);
    assert_eq!(p(&[3, 2]).z_coefficient(), 6);
    assert_eq!(p(&[3, 1, 1]).z_coefficient(), 6);
    assert_eq!(p(&[2, 2, 1]).z_coefficient(), 8);
    assert_eq!(p(&[2, 1, 1, 1]).z_coefficient(), 12);
    assert_eq!(p(&[1, 1, 1, 1, 1]).z_coefficient(), 120);
    assert_eq!(p(&[6]).z_coefficient(), 6);
    assert_eq!(p(&[3, 3]).z_coefficient(), 18);
    assert_eq!(p(&[3, 2, 1]).z_coefficient(), 6);
    assert_eq!(p(&[2, 2, 2]).z_coefficient(), 48);
}

#[test]
fn test_partition_core_mathematica() {
    // PartitionCore[{4,3,2,1}, 2] = {4,3,2,1} (is a 2-core: no even hooks)
    assert_eq!(p(&[4, 3, 2, 1]).partition_core(2), p(&[4, 3, 2, 1]));
    // PartitionCore[{4,3,2,1}, 3] = {1}
    assert_eq!(p(&[4, 3, 2, 1]).partition_core(3), p(&[1]));
    // (3,2,1) is also a 2-core (staircase partitions are 2-cores)
    assert_eq!(p(&[3, 2, 1]).partition_core(2), p(&[3, 2, 1]));
    assert_eq!(Partition::empty().partition_core(2), Partition::empty());
}

#[test]
fn test_partition_quotient() {
    // PartitionQuotient[{4,3,2,1}, 2] should be 2 partitions summing correctly
    let quot = p(&[4, 3, 2, 1]).partition_quotient(2);
    assert_eq!(quot.len(), 2);
    // Core is (4,3,2,1) itself (2-core), so quotient sizes should sum to 0
    let quot_sizes: u32 = quot.iter().map(|q| q.size()).sum();
    assert_eq!(quot_sizes, 0); // (|λ| - |core|) / d = (10-10)/2 = 0

    // PartitionQuotient[{4,3,2,1}, 3]: core = {1}, so quotient sizes sum to (10-1)/3 = 3
    let quot3 = p(&[4, 3, 2, 1]).partition_quotient(3);
    assert_eq!(quot3.len(), 3);
    let quot3_sizes: u32 = quot3.iter().map(|q| q.size()).sum();
    assert_eq!(quot3_sizes, 3);
}

#[test]
fn test_remove_vertical_strip() {
    let lam = p(&[3, 2]);
    let strips = lam.partition_remove_vertical_strip(2);
    for mu in &strips {
        assert_eq!(lam.size() - mu.size(), 2);
        assert!(mu.partition_less_equal(&lam));
    }
    assert!(!strips.is_empty());
}

#[test]
fn test_all_of_size_bounded() {
    // Partitions of 6 with at most 3 parts, each at most 3
    let parts = Partition::all_of_size_bounded(6, 3, 3);
    // Only (3,3), (3,2,1), (2,2,2) fit
    assert_eq!(parts.len(), 3);
    assert!(parts.contains(&p(&[3, 3])));
    assert!(parts.contains(&p(&[3, 2, 1])));
    assert!(parts.contains(&p(&[2, 2, 2])));
}

#[test]
fn test_all_of_size_dominance_order() {
    // In dominance order, (n) is first, (1^n) is last
    let parts = Partition::all_of_size_dominance_order(4);
    assert_eq!(parts.len(), 5);
    assert_eq!(parts[0], p(&[4]));
    assert_eq!(parts[parts.len() - 1], p(&[1, 1, 1, 1]));
    // Check it's a valid linear extension: each partition dominates its successor
    for i in 0..parts.len() - 1 {
        assert!(
            parts[i + 1].dominated_by(&parts[i]) || !parts[i].dominated_by(&parts[i + 1]),
            "dominance order violation at {} vs {}",
            parts[i],
            parts[i + 1]
        );
    }
}

#[test]
fn test_partition_enumeration_counts() {
    // p(0)=1, p(1)=1, p(2)=2, p(3)=3, p(4)=5, p(5)=7, p(6)=11, p(7)=15
    assert_eq!(Partition::all_of_size(0).len(), 1);
    assert_eq!(Partition::all_of_size(1).len(), 1);
    assert_eq!(Partition::all_of_size(2).len(), 2);
    assert_eq!(Partition::all_of_size(3).len(), 3);
    assert_eq!(Partition::all_of_size(4).len(), 5);
    assert_eq!(Partition::all_of_size(5).len(), 7);
    assert_eq!(Partition::all_of_size(6).len(), 11);
    assert_eq!(Partition::all_of_size(7).len(), 15);
    assert_eq!(Partition::all_of_size(10).len(), 42);
}

#[test]
fn test_empty_partition_edge_cases() {
    let e = Partition::empty();
    assert!(e.is_empty());
    assert_eq!(e.size(), 0);
    assert_eq!(e.num_parts(), 0);
    assert_eq!(e.part(0), 0);
    assert_eq!(e.conjugate_partition(), Partition::empty());
    assert_eq!(e.partition_n(), 0);
    assert_eq!(e.hook_lengths(), Vec::<Vec<u32>>::new());
    assert_eq!(e.diagram_boxes(), Vec::<(usize, usize)>::new());
    assert_eq!(e.durfee_square(), 0);
    assert_eq!(e.partition_add_box(), vec![p(&[1])]);
    assert_eq!(e.partition_remove_box(), vec![]);
    assert_eq!(e.partition_core(2), Partition::empty());
}

#[test]
fn test_parse_edge_cases() {
    assert_eq!(Partition::parse("").unwrap(), Partition::empty());
    assert_eq!(Partition::parse("0").unwrap(), Partition::empty());
    assert_eq!(Partition::parse("1").unwrap(), p(&[1]));
    assert_eq!(Partition::parse("3,2,1").unwrap(), p(&[3, 2, 1]));
    assert!(Partition::parse("abc").is_err());
}

// ===========================================================================
// Composition tests
// ===========================================================================

#[test]
fn test_composition_enumeration_counts() {
    // |C(n)| = 2^{n-1} for n >= 1
    assert_eq!(Composition::integer_compositions(0).len(), 1);
    assert_eq!(Composition::integer_compositions(1).len(), 1);
    assert_eq!(Composition::integer_compositions(2).len(), 2);
    assert_eq!(Composition::integer_compositions(3).len(), 4);
    assert_eq!(Composition::integer_compositions(4).len(), 8);
    assert_eq!(Composition::integer_compositions(5).len(), 16);
    assert_eq!(Composition::integer_compositions(6).len(), 32);
}

#[test]
fn test_composition_k_counts() {
    // C(n,k) = binom(n-1, k-1)
    assert_eq!(Composition::integer_compositions_k(5, 1).len(), 1);
    assert_eq!(Composition::integer_compositions_k(5, 2).len(), 4);
    assert_eq!(Composition::integer_compositions_k(5, 3).len(), 6);
    assert_eq!(Composition::integer_compositions_k(5, 4).len(), 4);
    assert_eq!(Composition::integer_compositions_k(5, 5).len(), 1);
    assert_eq!(Composition::integer_compositions_k(5, 6).len(), 0);
    assert_eq!(Composition::integer_compositions_k(0, 0).len(), 1);
    assert_eq!(Composition::integer_compositions_k(3, 0).len(), 0);
}

#[test]
fn test_weak_composition_counts() {
    // Weak compositions of n into k parts: binom(n+k-1, k-1)
    assert_eq!(Composition::weak_integer_compositions(3, 3).len(), 10); // binom(5,2)
    assert_eq!(Composition::weak_integer_compositions(0, 3).len(), 1);
    assert_eq!(Composition::weak_integer_compositions(2, 1).len(), 1);
}

#[test]
fn test_composition_descent_set_mathematica() {
    // CompositionToDescentSet[{2,1,3}] = {2, 3}
    let c = Composition::new(vec![2, 1, 3]);
    let des = c.composition_to_descent_set();
    let expected: BTreeSet<u32> = [2, 3].into_iter().collect();
    assert_eq!(des, expected);
}

#[test]
fn test_composition_ribbon_mathematica() {
    // CompositionToRibbon[{2,1,3}] = ({4,2,2}, {1,1})
    let c = Composition::new(vec![2, 1, 3]);
    let (outer, inner) = c.composition_to_ribbon();
    // Check size: |outer| - |inner| = |composition|
    assert_eq!(outer.size() - inner.size(), 6);
    // Inner must be contained in outer
    assert!(inner.partition_less_equal(&outer));
}

#[test]
fn test_composition_refinements_mathematica() {
    // CompositionRefinements[{4}] has 8 refinements (= compositions of 4)
    let refs = Composition::new(vec![4]).composition_refinements();
    assert_eq!(refs.len(), 8);

    // CompositionRefinements[{2,1}] = all ways to split 2 then 1
    // Splits of 2: (2), (1,1). Splits of 1: (1). So 2 * 1 = 2 refinements.
    let refs2 = Composition::new(vec![2, 1]).composition_refinements();
    assert_eq!(refs2.len(), 2);
}

#[test]
fn test_composition_to_partition() {
    assert_eq!(
        Composition::new(vec![2, 3, 1]).to_partition(),
        p(&[3, 2, 1])
    );
}

#[test]
fn test_empty_composition() {
    let c = Composition::empty();
    assert!(c.is_empty());
    assert_eq!(c.size(), 0);
    assert_eq!(c.num_parts(), 0);
    assert_eq!(c.composition_to_descent_set(), BTreeSet::new());
}

// ===========================================================================
// Ring trait tests
// ===========================================================================

#[test]
fn test_ring_i64() {
    assert_eq!(i64::zero(), 0);
    assert_eq!(i64::one(), 1);
    assert_eq!(i64::minus_one(), -1);
    assert!(i64::zero().is_zero());
    assert!(!i64::one().is_zero());
    assert_eq!(i64::from_i64(-42), -42);
    assert_eq!(42i64.exact_div_i64(7), 6);
    assert_eq!(0i64.exact_div_i64(5), 0);
    assert_eq!((-12i64).exact_div_i64(3), -4);
}

#[test]
#[should_panic(expected = "inexact division")]
fn test_ring_i64_inexact_div_panics() {
    7i64.exact_div_i64(3);
}

#[test]
fn test_ring_rational() {
    let half = q(1, 2);
    let third = q(1, 3);
    assert_eq!(half.clone() + third.clone(), q(5, 6));
    assert_eq!(half.clone() * third.clone(), q(1, 6));
    assert!(!half.is_zero());
    assert!(Q::zero().is_zero());
    // Exact division always works for rationals
    assert_eq!(half.exact_div_i64(3), q(1, 6));
    assert_eq!(Q::from_i64(10).exact_div_i64(3), q(10, 3));
}

#[test]
fn test_ring_bigint() {
    use num_bigint::BigInt;
    let a = BigInt::from_i64(100);
    let b = BigInt::from_i64(-50);
    assert_eq!(a.clone() + b.clone(), BigInt::from(50));
    assert_eq!(a.clone() * b.clone(), BigInt::from(-5000));
    assert_eq!(a.exact_div_i64(25), BigInt::from(4));
}

// ===========================================================================
// Kostka coefficients — verified against Mathematica
// ===========================================================================

#[test]
fn test_kostka_matrix_n4_mathematica() {
    // Full K matrix for n=4, verified against Mathematica
    let expected: Vec<(&[u32], &[u32], i64)> = vec![
        (&[4], &[4], 1),
        (&[4], &[3, 1], 1),
        (&[4], &[2, 2], 1),
        (&[4], &[2, 1, 1], 1),
        (&[4], &[1, 1, 1, 1], 1),
        (&[3, 1], &[4], 0),
        (&[3, 1], &[3, 1], 1),
        (&[3, 1], &[2, 2], 1),
        (&[3, 1], &[2, 1, 1], 2),
        (&[3, 1], &[1, 1, 1, 1], 3),
        (&[2, 2], &[4], 0),
        (&[2, 2], &[3, 1], 0),
        (&[2, 2], &[2, 2], 1),
        (&[2, 2], &[2, 1, 1], 1),
        (&[2, 2], &[1, 1, 1, 1], 2),
        (&[2, 1, 1], &[4], 0),
        (&[2, 1, 1], &[3, 1], 0),
        (&[2, 1, 1], &[2, 2], 0),
        (&[2, 1, 1], &[2, 1, 1], 1),
        (&[2, 1, 1], &[1, 1, 1, 1], 3),
        (&[1, 1, 1, 1], &[4], 0),
        (&[1, 1, 1, 1], &[3, 1], 0),
        (&[1, 1, 1, 1], &[2, 2], 0),
        (&[1, 1, 1, 1], &[2, 1, 1], 0),
        (&[1, 1, 1, 1], &[1, 1, 1, 1], 1),
    ];
    for (lam, mu, expected_k) in expected {
        let got = kostka::kostka_coefficient(&p(lam), &p(mu));
        assert_eq!(got, expected_k, "K({:?}, {:?})", lam, mu);
    }
}

#[test]
fn test_kostka_matrix_n5_mathematica() {
    // Spot-check several entries from n=5 Mathematica output
    assert_eq!(kostka::kostka_coefficient(&p(&[5]), &p(&[1, 1, 1, 1, 1])), 1);
    assert_eq!(kostka::kostka_coefficient(&p(&[4, 1]), &p(&[1, 1, 1, 1, 1])), 4);
    assert_eq!(kostka::kostka_coefficient(&p(&[3, 2]), &p(&[1, 1, 1, 1, 1])), 5);
    assert_eq!(kostka::kostka_coefficient(&p(&[3, 1, 1]), &p(&[1, 1, 1, 1, 1])), 6);
    assert_eq!(kostka::kostka_coefficient(&p(&[2, 2, 1]), &p(&[1, 1, 1, 1, 1])), 5);
    assert_eq!(kostka::kostka_coefficient(&p(&[2, 1, 1, 1]), &p(&[1, 1, 1, 1, 1])), 4);
    assert_eq!(kostka::kostka_coefficient(&p(&[4, 1]), &p(&[3, 1, 1])), 2);
    assert_eq!(kostka::kostka_coefficient(&p(&[3, 2]), &p(&[3, 1, 1])), 1);
    assert_eq!(kostka::kostka_coefficient(&p(&[3, 2]), &p(&[2, 2, 1])), 2);
    assert_eq!(kostka::kostka_coefficient(&p(&[2, 2, 1]), &p(&[2, 1, 1, 1])), 2);
}

#[test]
fn test_kostka_edge_cases() {
    // Size mismatch → 0
    assert_eq!(kostka::kostka_coefficient(&p(&[3]), &p(&[2])), 0);
    // Empty
    assert_eq!(
        kostka::kostka_coefficient(&Partition::empty(), &Partition::empty()),
        1
    );
    // K(λ, λ) = 1 for all λ (diagonal)
    for n in 0..=6 {
        for lam in Partition::all_of_size(n) {
            assert_eq!(
                kostka::kostka_coefficient(&lam, &lam),
                1,
                "K({}, {}) should be 1",
                lam,
                lam
            );
        }
    }
}

#[test]
fn test_inverse_kostka_product_identity_n4() {
    let (parts, inv) = kostka::inverse_kostka_matrix(4);
    let k = parts.len();
    let mut kmat = vec![vec![0i64; k]; k];
    for i in 0..k {
        for j in 0..k {
            kmat[i][j] = kostka::kostka_coefficient(&parts[i], &parts[j]);
        }
    }
    for i in 0..k {
        for j in 0..k {
            let mut sum = 0i64;
            for l in 0..k {
                sum += kmat[i][l] * inv[l][j];
            }
            let expected = if i == j { 1 } else { 0 };
            assert_eq!(sum, expected, "K*K^-1 != I at ({},{})", i, j);
        }
    }
}

// ===========================================================================
// Character table — verified against Mathematica
// ===========================================================================

#[test]
fn test_character_table_s4_mathematica() {
    // Full S_4 character table from Mathematica
    // Rows: {4}, {3,1}, {2,2}, {2,1,1}, {1^4}
    // Cols: {4}, {3,1}, {2,2}, {2,1,1}, {1^4}
    let expected: Vec<Vec<i64>> = vec![
        vec![1, 1, 1, 1, 1],
        vec![-1, 0, -1, 1, 3],
        vec![0, -1, 2, 0, 2],
        vec![1, 0, -1, -1, 3],
        vec![-1, 1, 1, -1, 1],
    ];
    let partitions = Partition::all_of_size(4);
    for (i, lam) in partitions.iter().enumerate() {
        for (j, mu) in partitions.iter().enumerate() {
            let got = kostka::sn_character(lam, mu);
            assert_eq!(
                got, expected[i][j],
                "χ^{}({}) = {} but expected {}",
                lam, mu, got, expected[i][j]
            );
        }
    }
}

#[test]
fn test_character_table_s5_mathematica() {
    // Full S_5 character table from Mathematica
    let expected: Vec<Vec<i64>> = vec![
        vec![1, 1, 1, 1, 1, 1, 1],
        vec![-1, 0, -1, 1, 0, 2, 4],
        vec![0, -1, 1, -1, 1, 1, 5],
        vec![1, 0, 0, 0, -2, 0, 6],
        vec![0, 1, -1, -1, 1, -1, 5],
        vec![-1, 0, 1, 1, 0, -2, 4],
        vec![1, -1, -1, 1, 1, -1, 1],
    ];
    let partitions = Partition::all_of_size(5);
    for (i, lam) in partitions.iter().enumerate() {
        for (j, mu) in partitions.iter().enumerate() {
            let got = kostka::sn_character(lam, mu);
            assert_eq!(
                got, expected[i][j],
                "χ^{}({}) = {} but expected {}",
                lam, mu, got, expected[i][j]
            );
        }
    }
}

#[test]
fn test_character_orthogonality_n5() {
    // Row orthogonality: Σ_μ χ^λ(μ) χ^ν(μ) / z_μ = δ_{λν}
    let (parts, table) = kostka::character_table(5);
    let k = parts.len();
    for i in 0..k {
        for j in 0..k {
            let mut sum = 0f64;
            for l in 0..k {
                let z = parts[l].z_coefficient() as f64;
                sum += (table[i][l] as f64) * (table[j][l] as f64) / z;
            }
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (sum - expected).abs() < 1e-10,
                "n=5 orthogonality failed for ({}, {}): got {}",
                i,
                j,
                sum
            );
        }
    }

    // Column orthogonality: Σ_λ χ^λ(μ) χ^λ(ν) = δ_{μν} * z_μ
    for i in 0..k {
        for j in 0..k {
            let mut sum = 0i64;
            for l in 0..k {
                sum += table[l][i] * table[l][j];
            }
            let expected = if i == j {
                parts[i].z_coefficient() as i64
            } else {
                0
            };
            assert_eq!(
                sum, expected,
                "n=5 column orthogonality failed for ({}, {})",
                i, j
            );
        }
    }
}

// ===========================================================================
// Basis conversions — verified against Mathematica SymmetricFunctions.m
// ===========================================================================

#[test]
fn test_schur_to_monomial_degree3_mathematica() {
    // s[2,1] = m[2,1] + 2*m[1,1,1]
    let s: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(p(&[2, 1]));
    let m = s.to_monomial_basis();
    assert_eq!(m.coefficient(&p(&[2, 1])), 1);
    assert_eq!(m.coefficient(&p(&[1, 1, 1])), 2);
    assert_eq!(m.coefficient(&p(&[3])), 0);
    assert_eq!(m.terms().len(), 2);
}

#[test]
fn test_schur_to_monomial_degree4_mathematica() {
    // s[3,1] = m[3,1] + m[2,2] + 2*m[2,1,1] + 3*m[1,1,1,1]
    let s: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(p(&[3, 1]));
    let m = s.to_monomial_basis();
    assert_eq!(m.coefficient(&p(&[3, 1])), 1);
    assert_eq!(m.coefficient(&p(&[2, 2])), 1);
    assert_eq!(m.coefficient(&p(&[2, 1, 1])), 2);
    assert_eq!(m.coefficient(&p(&[1, 1, 1, 1])), 3);
    assert_eq!(m.coefficient(&p(&[4])), 0);

    // s[2,2] = m[2,2] + m[2,1,1] + 2*m[1,1,1,1]
    let s22: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(p(&[2, 2]));
    let m22 = s22.to_monomial_basis();
    assert_eq!(m22.coefficient(&p(&[2, 2])), 1);
    assert_eq!(m22.coefficient(&p(&[2, 1, 1])), 1);
    assert_eq!(m22.coefficient(&p(&[1, 1, 1, 1])), 2);
    assert_eq!(m22.terms().len(), 3);

    // s[2,1,1] = m[2,1,1] + 3*m[1,1,1,1]
    let s211: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(p(&[2, 1, 1]));
    let m211 = s211.to_monomial_basis();
    assert_eq!(m211.coefficient(&p(&[2, 1, 1])), 1);
    assert_eq!(m211.coefficient(&p(&[1, 1, 1, 1])), 3);
    assert_eq!(m211.terms().len(), 2);
}

#[test]
fn test_elementary_to_schur_mathematica() {
    // e[2,1] = s[2,1] + s[1,1,1]
    let e: SymmetricFunction<i64> = SymmetricFunction::elementary_symmetric(p(&[2, 1]));
    let s = e.to_schur_basis();
    assert_eq!(s.coefficient(&p(&[2, 1])), 1);
    assert_eq!(s.coefficient(&p(&[1, 1, 1])), 1);
    assert_eq!(s.terms().len(), 2);

    // e[4] = s[1,1,1,1]
    let e4: SymmetricFunction<i64> = SymmetricFunction::elementary_symmetric(p(&[4]));
    let s4 = e4.to_schur_basis();
    assert_eq!(s4.coefficient(&p(&[1, 1, 1, 1])), 1);
    assert_eq!(s4.terms().len(), 1);
}

#[test]
fn test_elementary_to_monomial_mathematica() {
    // e[3] = m[1,1,1]
    let e: SymmetricFunction<i64> = SymmetricFunction::elementary_symmetric(p(&[3]));
    let m = e.to_monomial_basis();
    assert_eq!(m.coefficient(&p(&[1, 1, 1])), 1);
    assert_eq!(m.terms().len(), 1);

    // e[2,1] = m[2,1] + 3*m[1,1,1]
    let e21: SymmetricFunction<i64> = SymmetricFunction::elementary_symmetric(p(&[2, 1]));
    let m21 = e21.to_monomial_basis();
    assert_eq!(m21.coefficient(&p(&[2, 1])), 1);
    assert_eq!(m21.coefficient(&p(&[1, 1, 1])), 3);
}

#[test]
fn test_complete_h_to_schur_mathematica() {
    // h[2,1] = s[3] + s[2,1]
    let h: SymmetricFunction<i64> = SymmetricFunction::complete_h_symmetric(p(&[2, 1]));
    let s = h.to_schur_basis();
    assert_eq!(s.coefficient(&p(&[3])), 1);
    assert_eq!(s.coefficient(&p(&[2, 1])), 1);
    assert_eq!(s.terms().len(), 2);

    // h[4] = s[4]
    let h4: SymmetricFunction<i64> = SymmetricFunction::complete_h_symmetric(p(&[4]));
    let s4 = h4.to_schur_basis();
    assert_eq!(s4.coefficient(&p(&[4])), 1);
    assert_eq!(s4.terms().len(), 1);
}

#[test]
fn test_complete_h_to_monomial_mathematica() {
    // h[3] = m[3] + m[2,1] + m[1,1,1]
    let h: SymmetricFunction<i64> = SymmetricFunction::complete_h_symmetric(p(&[3]));
    let m = h.to_monomial_basis();
    assert_eq!(m.coefficient(&p(&[3])), 1);
    assert_eq!(m.coefficient(&p(&[2, 1])), 1);
    assert_eq!(m.coefficient(&p(&[1, 1, 1])), 1);

    // h[2,1] = m[3] + 2*m[2,1] + 3*m[1,1,1]
    let h21: SymmetricFunction<i64> = SymmetricFunction::complete_h_symmetric(p(&[2, 1]));
    let m21 = h21.to_monomial_basis();
    assert_eq!(m21.coefficient(&p(&[3])), 1);
    assert_eq!(m21.coefficient(&p(&[2, 1])), 2);
    assert_eq!(m21.coefficient(&p(&[1, 1, 1])), 3);
}

#[test]
fn test_forgotten_to_monomial_degree3_mathematica() {
    // f[3] = m[3]
    let f: SymmetricFunction<i64> = SymmetricFunction::forgotten_symmetric(p(&[3]));
    let m = f.to_monomial_basis();
    assert_eq!(m.coefficient(&p(&[3])), 1);
    assert_eq!(m.terms().len(), 1);

    // f[2,1] = -2*m[3] - m[2,1]
    let f21: SymmetricFunction<i64> = SymmetricFunction::forgotten_symmetric(p(&[2, 1]));
    let m21 = f21.to_monomial_basis();
    assert_eq!(m21.coefficient(&p(&[3])), -2);
    assert_eq!(m21.coefficient(&p(&[2, 1])), -1);

    // f[1,1,1] = m[3] + m[2,1] + m[1,1,1]
    let f111: SymmetricFunction<i64> = SymmetricFunction::forgotten_symmetric(p(&[1, 1, 1]));
    let m111 = f111.to_monomial_basis();
    assert_eq!(m111.coefficient(&p(&[3])), 1);
    assert_eq!(m111.coefficient(&p(&[2, 1])), 1);
    assert_eq!(m111.coefficient(&p(&[1, 1, 1])), 1);
}

#[test]
fn test_forgotten_to_monomial_degree4_mathematica() {
    // f[4] = -m[4]
    let f4: SymmetricFunction<i64> = SymmetricFunction::forgotten_symmetric(p(&[4]));
    let m = f4.to_monomial_basis();
    assert_eq!(m.coefficient(&p(&[4])), -1);
    assert_eq!(m.terms().len(), 1);

    // f[3,1] = 2*m[4] + m[3,1]
    let f31: SymmetricFunction<i64> = SymmetricFunction::forgotten_symmetric(p(&[3, 1]));
    let m31 = f31.to_monomial_basis();
    assert_eq!(m31.coefficient(&p(&[4])), 2);
    assert_eq!(m31.coefficient(&p(&[3, 1])), 1);

    // f[2,1,1] = -3*m[4] - 2*m[3,1] - 2*m[2,2] - m[2,1,1]
    let f211: SymmetricFunction<i64> = SymmetricFunction::forgotten_symmetric(p(&[2, 1, 1]));
    let m211 = f211.to_monomial_basis();
    assert_eq!(m211.coefficient(&p(&[4])), -3);
    assert_eq!(m211.coefficient(&p(&[3, 1])), -2);
    assert_eq!(m211.coefficient(&p(&[2, 2])), -2);
    assert_eq!(m211.coefficient(&p(&[2, 1, 1])), -1);
}

#[test]
fn test_power_sum_to_schur_mathematica() {
    // p[2,1] = s[3] - s[1,1,1]
    let ps: SymmetricFunction<i64> = SymmetricFunction::power_sum_symmetric(p(&[2, 1]));
    let s = ps.to_schur_basis();
    assert_eq!(s.coefficient(&p(&[3])), 1);
    assert_eq!(s.coefficient(&p(&[1, 1, 1])), -1);
    assert_eq!(s.coefficient(&p(&[2, 1])), 0);

    // p[4] = s[4] - s[3,1] + s[2,1,1] - s[1,1,1,1]
    let p4: SymmetricFunction<i64> = SymmetricFunction::power_sum_symmetric(p(&[4]));
    let s4 = p4.to_schur_basis();
    assert_eq!(s4.coefficient(&p(&[4])), 1);
    assert_eq!(s4.coefficient(&p(&[3, 1])), -1);
    assert_eq!(s4.coefficient(&p(&[2, 2])), 0);
    assert_eq!(s4.coefficient(&p(&[2, 1, 1])), 1);
    assert_eq!(s4.coefficient(&p(&[1, 1, 1, 1])), -1);
}

#[test]
fn test_schur_to_power_sum_degree3_mathematica() {
    // s[3] = (1/3)*p[3] + (1/2)*p[2,1] + (1/6)*p[1,1,1]
    let s3: SymmetricFunction<Q> = SymmetricFunction::schur_symmetric(p(&[3]));
    let ps = s3.to_power_sum_basis();
    assert_eq!(ps.coefficient(&p(&[3])), q(1, 3));
    assert_eq!(ps.coefficient(&p(&[2, 1])), q(1, 2));
    assert_eq!(ps.coefficient(&p(&[1, 1, 1])), q(1, 6));

    // s[2,1] = -(1/3)*p[3] + (1/3)*p[1,1,1]
    let s21: SymmetricFunction<Q> = SymmetricFunction::schur_symmetric(p(&[2, 1]));
    let ps21 = s21.to_power_sum_basis();
    assert_eq!(ps21.coefficient(&p(&[3])), q(-1, 3));
    assert_eq!(ps21.coefficient(&p(&[2, 1])), q(0, 1));
    assert_eq!(ps21.coefficient(&p(&[1, 1, 1])), q(1, 3));

    // s[1,1,1] = (1/3)*p[3] - (1/2)*p[2,1] + (1/6)*p[1,1,1]
    let s111: SymmetricFunction<Q> = SymmetricFunction::schur_symmetric(p(&[1, 1, 1]));
    let ps111 = s111.to_power_sum_basis();
    assert_eq!(ps111.coefficient(&p(&[3])), q(1, 3));
    assert_eq!(ps111.coefficient(&p(&[2, 1])), q(-1, 2));
    assert_eq!(ps111.coefficient(&p(&[1, 1, 1])), q(1, 6));
}

#[test]
fn test_schur_to_power_sum_degree4_mathematica() {
    // s[3,1] = -(1/4)*p[4] - (1/8)*p[2,2] + (1/4)*p[2,1,1] + (1/8)*p[1^4]
    let s31: SymmetricFunction<Q> = SymmetricFunction::schur_symmetric(p(&[3, 1]));
    let ps = s31.to_power_sum_basis();
    assert_eq!(ps.coefficient(&p(&[4])), q(-1, 4));
    assert_eq!(ps.coefficient(&p(&[3, 1])), q(0, 1));
    assert_eq!(ps.coefficient(&p(&[2, 2])), q(-1, 8));
    assert_eq!(ps.coefficient(&p(&[2, 1, 1])), q(1, 4));
    assert_eq!(ps.coefficient(&p(&[1, 1, 1, 1])), q(1, 8));
}

// ===========================================================================
// Roundtrip tests: convert through all basis pairs
// ===========================================================================

#[test]
fn test_all_basis_roundtrips_degree4() {
    // For each basis pair, convert there and back, verify identity
    let bases = [
        Basis::Monomial,
        Basis::Elementary,
        Basis::CompleteH,
        Basis::Schur,
        Basis::Forgotten,
    ];

    // Start with s[3,1] + 2*s[2,2] (a non-trivial degree-4 function)
    let original: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(p(&[3, 1]))
        + SymmetricFunction::scaled_basis_element(Basis::Schur, p(&[2, 2]), 2);

    for &basis in &bases {
        let converted = original.to_basis(basis);
        let back = converted.to_schur_basis();
        assert_eq!(
            back.coefficient(&p(&[3, 1])),
            1,
            "roundtrip through {:?} failed for s[3,1]",
            basis
        );
        assert_eq!(
            back.coefficient(&p(&[2, 2])),
            2,
            "roundtrip through {:?} failed for s[2,2]",
            basis
        );
        assert_eq!(
            back.coefficient(&p(&[4])),
            0,
            "roundtrip through {:?} introduced s[4]",
            basis
        );
    }
}

#[test]
fn test_power_sum_roundtrips_degree3() {
    // Roundtrip through power-sum (requires rational)
    let bases = [
        Basis::Monomial,
        Basis::Elementary,
        Basis::CompleteH,
        Basis::Schur,
    ];

    let original: SymmetricFunction<Q> = SymmetricFunction::schur_symmetric(p(&[2, 1]));
    for &basis in &bases {
        let via_ps = original.to_power_sum_basis().to_basis(basis);
        let direct = original.to_basis(basis);
        for lam in Partition::all_of_size(3) {
            assert_eq!(
                via_ps.coefficient(&lam),
                direct.coefficient(&lam),
                "power-sum roundtrip failed via {:?} at {}",
                basis,
                lam
            );
        }
    }
}

// ===========================================================================
// Symmetric function operations
// ===========================================================================

#[test]
fn test_multiply_schur_via_conversion() {
    // s[2] * s[1] should give s[3] + s[2,1] (by Pieri's rule)
    let s2: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(p(&[2]));
    let s1: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(p(&[1]));
    let prod = s2.multiply(&s1);
    assert_eq!(prod.basis(), Basis::Schur);
    assert_eq!(prod.coefficient(&p(&[3])), 1);
    assert_eq!(prod.coefficient(&p(&[2, 1])), 1);
    assert_eq!(prod.terms().len(), 2);
}

#[test]
fn test_multiply_schur_degree4() {
    // s[2] * s[2] = s[4] + s[3,1] + s[2,2]
    let s2: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(p(&[2]));
    let prod = s2.multiply(&s2.clone());
    assert_eq!(prod.coefficient(&p(&[4])), 1);
    assert_eq!(prod.coefficient(&p(&[3, 1])), 1);
    assert_eq!(prod.coefficient(&p(&[2, 2])), 1);
    assert_eq!(prod.coefficient(&p(&[2, 1, 1])), 0);
    assert_eq!(prod.terms().len(), 3);
}

#[test]
fn test_multiply_elementary() {
    // e[2] * e[1] = e[2,1] (multiplicative basis)
    let e2: SymmetricFunction<i64> = SymmetricFunction::elementary_symmetric(p(&[2]));
    let e1: SymmetricFunction<i64> = SymmetricFunction::elementary_symmetric(p(&[1]));
    let prod = e2.multiply(&e1);
    assert_eq!(prod.coefficient(&p(&[2, 1])), 1);
    assert_eq!(prod.terms().len(), 1);
}

#[test]
fn test_scale() {
    let s: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(p(&[2, 1]));
    let scaled = s.scale(&3);
    assert_eq!(scaled.coefficient(&p(&[2, 1])), 3);

    let zero_scaled = s.scale(&0);
    assert!(zero_scaled.is_zero());
}

#[test]
fn test_sub_and_neg() {
    let s1: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(p(&[2, 1]));
    let s2: SymmetricFunction<i64> =
        SymmetricFunction::scaled_basis_element(Basis::Schur, p(&[2, 1]), 3);
    let diff = s2 - s1;
    assert_eq!(diff.coefficient(&p(&[2, 1])), 2);

    let s: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(p(&[2, 1]));
    let neg = -s;
    assert_eq!(neg.coefficient(&p(&[2, 1])), -1);
}

#[test]
fn test_degree() {
    let s: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(p(&[2, 1]));
    assert_eq!(s.degree(), Some(3));

    let mixed: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(p(&[2, 1]))
        + SymmetricFunction::schur_symmetric(p(&[3, 1]));
    assert_eq!(mixed.degree(), None); // mixed degree

    let z: SymmetricFunction<i64> = SymmetricFunction::zero(Basis::Schur);
    assert_eq!(z.degree(), None);
}

#[test]
fn test_omega_involution_mathematica() {
    // ω(e_n) = h_n: check via Schur conversion
    // e_3 = s[1,1,1], h_3 = s[3], ω(s[1,1,1]) = s[3]
    let e3: SymmetricFunction<i64> = SymmetricFunction::elementary_symmetric(p(&[3]));
    let omega_e3 = e3.omega_involution();
    assert_eq!(omega_e3.basis(), Basis::CompleteH);
    let omega_e3_schur = omega_e3.to_schur_basis();
    // h_3 in Schur = s[3]
    assert_eq!(omega_e3_schur.coefficient(&p(&[3])), 1);
    assert_eq!(omega_e3_schur.terms().len(), 1);

    // ω(p_{2,1}) = (-1)^{3-2} p_{2,1} = -p_{2,1}
    let ps: SymmetricFunction<i64> = SymmetricFunction::power_sum_symmetric(p(&[2, 1]));
    let omega_ps = ps.omega_involution();
    assert_eq!(omega_ps.coefficient(&p(&[2, 1])), -1);

    // ω(p_{1,1,1}) = (-1)^{3-3} p_{1,1,1} = p_{1,1,1}
    let p111: SymmetricFunction<i64> = SymmetricFunction::power_sum_symmetric(p(&[1, 1, 1]));
    let omega_p111 = p111.omega_involution();
    assert_eq!(omega_p111.coefficient(&p(&[1, 1, 1])), 1);

    // ω is involutory on every basis
    let f: SymmetricFunction<i64> = SymmetricFunction::forgotten_symmetric(p(&[2, 1]));
    let omega2 = f.omega_involution().omega_involution();
    assert_eq!(omega2.basis(), Basis::Forgotten);
    assert_eq!(omega2.coefficient(&p(&[2, 1])), 1);
}

#[test]
fn test_hall_inner_product_mathematica() {
    // <s[2,1], s[2,1]> = 1 (Schur orthonormality)
    let s21: SymmetricFunction<Q> = SymmetricFunction::schur_symmetric(p(&[2, 1]));
    assert_eq!(s21.hall_inner_product(&s21), Q::one());

    // <s[3], s[2,1]> = 0
    let s3: SymmetricFunction<Q> = SymmetricFunction::schur_symmetric(p(&[3]));
    assert_eq!(s3.hall_inner_product(&s21), Q::zero());

    // <s[3], s[3]> = 1
    assert_eq!(s3.hall_inner_product(&s3), Q::one());

    // <p[2,1], p[2,1]> = z_{2,1} = 2
    let p21: SymmetricFunction<Q> = SymmetricFunction::power_sum_symmetric(p(&[2, 1]));
    assert_eq!(p21.hall_inner_product(&p21), Q::from_i64(2));

    // <p[3], p[3]> = z_3 = 3
    let p3: SymmetricFunction<Q> = SymmetricFunction::power_sum_symmetric(p(&[3]));
    assert_eq!(p3.hall_inner_product(&p3), Q::from_i64(3));

    // <p[3], p[2,1]> = 0
    assert_eq!(p3.hall_inner_product(&p21), Q::zero());

    // <h[2,1], e[2,1]> = 1
    let h21: SymmetricFunction<Q> = SymmetricFunction::complete_h_symmetric(p(&[2, 1]));
    let e21: SymmetricFunction<Q> = SymmetricFunction::elementary_symmetric(p(&[2, 1]));
    assert_eq!(h21.hall_inner_product(&e21), Q::one());

    // <h[2,1], e[1,1,1]> = 3
    let e111: SymmetricFunction<Q> = SymmetricFunction::elementary_symmetric(p(&[1, 1, 1]));
    assert_eq!(h21.hall_inner_product(&e111), Q::from_i64(3));
}

#[test]
fn test_positive_coefficients() {
    let s: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(p(&[2, 1]));
    let m = s.to_monomial_basis();
    assert!(m.positive_coefficients()); // Kostka numbers are non-negative

    let f: SymmetricFunction<i64> = SymmetricFunction::forgotten_symmetric(p(&[2, 1]));
    let fm = f.to_monomial_basis();
    assert!(!fm.positive_coefficients()); // f[2,1] = -2*m[3] - m[2,1]
}

// ===========================================================================
// Conversion edge cases
// ===========================================================================

#[test]
fn test_conversion_degree_0() {
    // The empty partition represents degree 0
    let s: SymmetricFunction<i64> = SymmetricFunction::schur_symmetric(Partition::empty());
    let m = s.to_monomial_basis();
    assert_eq!(m.coefficient(&Partition::empty()), 1);
}

#[test]
fn test_conversion_degree_1() {
    // All bases agree at degree 1: s[1] = m[1] = e[1] = h[1] = p[1]
    let bases = [
        Basis::Monomial,
        Basis::Elementary,
        Basis::CompleteH,
        Basis::Schur,
    ];
    for &from in &bases {
        for &to in &bases {
            let f: SymmetricFunction<i64> = SymmetricFunction::basis_element(from, p(&[1]));
            let converted = f.to_basis(to);
            assert_eq!(
                converted.coefficient(&p(&[1])),
                1,
                "{:?} -> {:?} at degree 1 failed",
                from,
                to
            );
            assert_eq!(converted.terms().len(), 1);
        }
    }
}

#[test]
fn test_conversion_zero_function() {
    let z: SymmetricFunction<i64> = SymmetricFunction::zero(Basis::Schur);
    let m = z.to_monomial_basis();
    assert!(m.is_zero());
    assert_eq!(m.basis(), Basis::Monomial);
}

#[test]
fn test_all_30_basis_pairs_degree3() {
    // Comprehensive: for each of the 30 ordered pairs of distinct bases,
    // convert s[2,1] and verify the roundtrip.
    let bases = [
        Basis::Monomial,
        Basis::Elementary,
        Basis::CompleteH,
        Basis::PowerSum,
        Basis::Schur,
        Basis::Forgotten,
    ];

    let original: SymmetricFunction<Q> = SymmetricFunction::schur_symmetric(p(&[2, 1]));

    for &b1 in &bases {
        let in_b1 = original.to_basis(b1);
        for &b2 in &bases {
            let in_b2 = in_b1.to_basis(b2);
            let back = in_b2.to_schur_basis();
            assert_eq!(
                back.coefficient(&p(&[2, 1])),
                Q::one(),
                "roundtrip {:?} -> {:?} -> Schur failed for coeff of [2,1]",
                b1,
                b2,
            );
            assert_eq!(
                back.coefficient(&p(&[3])),
                Q::zero(),
                "roundtrip {:?} -> {:?} -> Schur introduced [3]",
                b1,
                b2,
            );
        }
    }
}

// ── Kostka-to-LR reduction ──────────────────────────────────────────────

#[test]
fn test_kostka_to_lr_reduction() {
    // Verify K_{λ,μ} = c^ρ_{λ,σ} for straight shapes (ν = ∅).
    // Compute both sides: K via DP, c^ρ_{λ,σ} via skew Schur expansion.
    let cases: Vec<(Vec<u32>, Vec<u32>)> = vec![
        (vec![3, 2, 1], vec![1, 1, 1, 1, 1, 1]),
        (vec![4, 2], vec![2, 2, 1, 1]),
        (vec![4, 3], vec![2, 2, 2, 1]),
        (vec![3, 3, 3], vec![1, 1, 1, 1, 1, 1, 1, 1, 1]),
        (vec![4, 3, 2, 1], vec![2, 2, 2, 2, 1, 1]),
    ];

    for (lam_v, mu_v) in &cases {
        let lambda = p(lam_v);
        let mu = p(mu_v);

        // Compute K_{λ, μ} via Kostka DP
        let k_val = kostka::kostka_coefficient(&lambda, &mu);

        // Compute c^ρ_{λ, σ} via the LR reduction
        let (rho, _lam, sigma) = Partition::kostka_to_lr(&lambda, &Partition::empty(), &mu);
        let skew: SymmetricFunction<i64> =
            SymmetricFunction::skew_schur_function(&rho, &sigma);
        let in_schur = skew.to_schur_basis();
        let c_val = in_schur.coefficient(&lambda);

        assert_eq!(
            k_val, c_val,
            "K_{{({}), ({})}} = {} but c^{{({})}}_{{{}, ({})}} = {}",
            lambda, mu, k_val, rho, lambda, sigma, c_val
        );
    }
}
