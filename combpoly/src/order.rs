//! Bruhat order and weak order ideals of permutations.

use crate::statistics::{self, Stat};
use std::collections::{BTreeSet, HashSet, VecDeque};

/// Compute the Bruhat lower order ideal of a permutation.
/// Returns all permutations sigma such that sigma <= perm in Bruhat order.
pub fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> {
    let n = perm.len();
    let mut visited: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut queue: BTreeSet<Vec<u8>> = BTreeSet::new();
    queue.insert(perm.to_vec());

    while let Some(current) = queue.pop_last() {
        for i in 0..n {
            for j in i + 1..n {
                if current[i] > current[j] {
                    let mut child = current.clone();
                    child.swap(i, j);
                    if !visited.contains(&child) {
                        queue.insert(child);
                    }
                }
            }
        }
        visited.insert(current);
    }

    visited.into_iter().collect()
}

/// Compute the weak lower order ideal of a permutation.
/// Returns all permutations sigma such that sigma <= perm in weak order.
pub fn weak_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> {
    let n = perm.len();
    if n == 0 {
        return vec![vec![]];
    }
    let mut visited: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut queue: BTreeSet<Vec<u8>> = BTreeSet::new();
    queue.insert(perm.to_vec());

    while let Some(current) = queue.pop_last() {
        for i in 0..n - 1 {
            if current[i] > current[i + 1] {
                let mut child = current.clone();
                child.swap(i, i + 1);
                if !visited.contains(&child) {
                    queue.insert(child);
                }
            }
        }
        visited.insert(current);
    }

    visited.into_iter().collect()
}

/// Compute the polynomial Σ t^{stat(σ)} over the weak lower ideal of `perm`,
/// without materializing the ideal.
///
/// Uses a HashSet for O(1) amortized membership checks (vs O(log n) for BTreeSet)
/// and a VecDeque for BFS. This is much faster for large ideals.
pub fn weak_ideal_poly(perm: &[u8], stat: Stat) -> Vec<i64> {
    let n = perm.len();
    if n == 0 {
        return vec![1];
    }
    let mut coeffs: Vec<i64> = Vec::new();
    let mut visited: HashSet<Vec<u8>> = HashSet::new();
    let mut queue: VecDeque<Vec<u8>> = VecDeque::new();

    visited.insert(perm.to_vec());
    queue.push_back(perm.to_vec());

    while let Some(current) = queue.pop_front() {
        // Accumulate statistic directly
        let val = statistics::compute(&current, stat);
        while coeffs.len() <= val {
            coeffs.push(0);
        }
        coeffs[val] += 1;

        // Generate children (adjacent transpositions reducing inversions)
        for i in 0..n - 1 {
            if current[i] > current[i + 1] {
                let mut child = current.clone();
                child.swap(i, i + 1);
                if visited.insert(child.clone()) {
                    queue.push_back(child);
                }
            }
        }
    }

    if coeffs.is_empty() {
        coeffs.push(0);
    }
    coeffs
}

/// Compute the polynomial Σ t^{stat(σ)} over the Bruhat lower ideal of `perm`,
/// without materializing the ideal.
pub fn bruhat_ideal_poly(perm: &[u8], stat: Stat) -> Vec<i64> {
    let n = perm.len();
    if n == 0 {
        return vec![1];
    }
    let mut coeffs: Vec<i64> = Vec::new();
    let mut visited: HashSet<Vec<u8>> = HashSet::new();
    let mut queue: VecDeque<Vec<u8>> = VecDeque::new();

    visited.insert(perm.to_vec());
    queue.push_back(perm.to_vec());

    while let Some(current) = queue.pop_front() {
        let val = statistics::compute(&current, stat);
        while coeffs.len() <= val {
            coeffs.push(0);
        }
        coeffs[val] += 1;

        for i in 0..n {
            for j in i + 1..n {
                if current[i] > current[j] {
                    let mut child = current.clone();
                    child.swap(i, j);
                    if visited.insert(child.clone()) {
                        queue.push_back(child);
                    }
                }
            }
        }
    }

    if coeffs.is_empty() {
        coeffs.push(0);
    }
    coeffs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bruhat_ideal_identity() {
        assert_eq!(bruhat_lower_ideal(&[1, 2, 3]).len(), 1);
    }

    #[test]
    fn test_bruhat_ideal_w0() {
        assert_eq!(bruhat_lower_ideal(&[3, 2, 1]).len(), 6);
        assert_eq!(bruhat_lower_ideal(&[4, 3, 2, 1]).len(), 24);
    }

    #[test]
    fn test_weak_ideal_w0() {
        assert_eq!(weak_lower_ideal(&[3, 2, 1]).len(), 6);
    }

    #[test]
    fn test_weak_subset_of_bruhat() {
        let perm = vec![3, 1, 4, 2];
        let bruhat: BTreeSet<_> = bruhat_lower_ideal(&perm).into_iter().collect();
        let weak: BTreeSet<_> = weak_lower_ideal(&perm).into_iter().collect();
        assert!(weak.is_subset(&bruhat));
    }
}
