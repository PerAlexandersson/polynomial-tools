//! Utilities for the fixed-descent insertion framework.
//!
//! This module packages the combinatorics that appears repeatedly in the
//! "insert the new maximum" recurrence for permutations with a prescribed
//! descent set.
//!
//! The central setup is:
//! - fix a target descent set `S ⊆ [n-1]`, encoded as a bitmask,
//! - choose an insertion position `p ∈ [n]`,
//! - insert the value `n` at position `p` into a permutation of `[n-1]`.
//!
//! The target descent set determines which insertion positions are valid,
//! which source descent sets can occur, and how the long-swaps statistic
//! changes under insertion.
//!
//! The helper names in this module are intentionally explicit so that the
//! corresponding mathematical objects are easy to identify in code and in
//! accompanying notes.

use std::collections::BTreeMap;

use crate::permutation::all_permutations;
use crate::statistics::{compute, descent_set_bitmask, Stat};

/// Data describing the local source-set relation for consecutive valid
/// insertion positions.
///
/// If `left_insertion_position < right_insertion_position` are consecutive
/// valid insertion positions for a target descent set `S`, then the source
/// descent sets differ in a highly localized way:
/// - the base source sets differ by exchanging at most one endpoint of the
///   descent run that starts at `left_insertion_position`,
/// - the augmented source sets do the same, together with the deterministic
///   boundary change at `left_insertion_position - 1`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConsecutiveValidInsertionPositionSourceDescentSetData {
    /// The size `n` of the target permutations in `S_n`.
    pub target_permutation_size: u8,
    /// The left element `p_a` of the consecutive valid pair.
    pub left_insertion_position: u8,
    /// The right element `p_b` of the consecutive valid pair.
    pub right_insertion_position: u8,
    /// The last descent position in the descent run of `S` that begins at
    /// `left_insertion_position`.
    pub descent_run_end_position: u8,
    /// The shared part of the two base source descent sets outside the local
    /// exchange.
    pub common_fixed_part_source_descent_set_mask: u64,
    /// The base source descent set `S'_{p_a}`.
    pub left_base_source_descent_set_mask: u64,
    /// The base source descent set `S'_{p_b}`.
    pub right_base_source_descent_set_mask: u64,
    /// The augmented source descent set `S''_{p_a}`, when it exists.
    pub left_augmented_source_descent_set_mask: Option<u64>,
    /// The augmented source descent set `S''_{p_b}`, when it exists.
    pub right_augmented_source_descent_set_mask: Option<u64>,
}

/// A `q`-refined modified source distribution for a fixed target descent set
/// and insertion position.
///
/// The source permutations live in `S_{n-1}` and are refined by
/// `q = pos(n-1)`, the 1-indexed position of the previous maximum in the
/// source permutation.
///
/// The key count stored here is
/// `N_p(s, q) = #{ π in source(p) : modified_source_swaps(π, p) = s, pos(n-1)=q }`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QRefinedModifiedSourceDistributionData {
    /// The size `n` of the target permutations in `S_n`.
    pub target_permutation_size: u8,
    /// The target descent set `S ⊆ [n-1]`.
    pub target_descent_set_mask: u64,
    /// The insertion position `p ∈ [n]`.
    pub insertion_position: u8,
    /// The base source descent set `S'_p`.
    pub base_source_descent_set_mask: u64,
    /// The augmented source descent set `S''_p`, when it exists.
    pub augmented_source_descent_set_mask: Option<u64>,
    /// The counts `N_p(s, q)`, keyed by `(modified_source_swaps, q)`.
    pub counts_by_modified_source_swaps_and_previous_maximum_position: BTreeMap<(usize, u8), usize>,
}

impl QRefinedModifiedSourceDistributionData {
    /// Build the polynomial `F̃_q(t)` for a fixed previous-maximum position `q`.
    pub fn modified_source_polynomial_for_previous_maximum_position(
        &self,
        previous_maximum_position: u8,
    ) -> Vec<i64> {
        let mut polynomial = Vec::<i64>::new();
        for (&(modified_source_swaps, q), &count) in
            &self.counts_by_modified_source_swaps_and_previous_maximum_position
        {
            if q != previous_maximum_position {
                continue;
            }
            if polynomial.len() <= modified_source_swaps {
                polynomial.resize(modified_source_swaps + 1, 0);
            }
            polynomial[modified_source_swaps] += count as i64;
        }
        if polynomial.is_empty() {
            return vec![0];
        }
        while polynomial.len() > 1 && *polynomial.last().unwrap() == 0 {
            polynomial.pop();
        }
        polynomial
    }

    /// Build the output polynomial `L^{(p)}(t)` obtained after applying the
    /// `ε₂` staircase shift to the refined source distribution.
    pub fn output_long_swaps_polynomial_after_insertion(&self) -> Vec<i64> {
        let mut polynomial = Vec::<i64>::new();
        for (&(modified_source_swaps, previous_maximum_position), &count) in
            &self.counts_by_modified_source_swaps_and_previous_maximum_position
        {
            let output_long_swaps = modified_source_swaps
                + usize::from(insertion_creates_new_long_swap_with_previous_maximum(
                    self.insertion_position,
                    previous_maximum_position,
                ));
            if polynomial.len() <= output_long_swaps {
                polynomial.resize(output_long_swaps + 1, 0);
            }
            polynomial[output_long_swaps] += count as i64;
        }
        if polynomial.is_empty() {
            return vec![0];
        }
        while polynomial.len() > 1 && *polynomial.last().unwrap() == 0 {
            polynomial.pop();
        }
        polynomial
    }
}

impl ConsecutiveValidInsertionPositionSourceDescentSetData {
    /// The positions appearing only in the left base source descent set.
    pub fn left_only_base_source_positions(&self) -> Vec<u8> {
        descent_positions_from_bitmask(
            self.left_base_source_descent_set_mask & !self.right_base_source_descent_set_mask,
            self.target_permutation_size - 1,
        )
    }

    /// The positions appearing only in the right base source descent set.
    pub fn right_only_base_source_positions(&self) -> Vec<u8> {
        descent_positions_from_bitmask(
            self.right_base_source_descent_set_mask & !self.left_base_source_descent_set_mask,
            self.target_permutation_size - 1,
        )
    }
}

/// Return whether the descent set bitmask contains the given 1-indexed
/// descent position.
pub fn target_descent_set_contains_position(target_descent_set_mask: u64, position: u8) -> bool {
    position >= 1 && ((target_descent_set_mask >> (position - 1)) & 1) == 1
}

/// Convert a descent-set bitmask into its sorted list of 1-indexed positions.
///
/// The parameter `n` refers to permutations in `S_n`, so valid descent
/// positions lie in `[1, n-1]`.
pub fn descent_positions_from_bitmask(descent_set_mask: u64, n: u8) -> Vec<u8> {
    let mut positions = Vec::new();
    for position in 1..n {
        if target_descent_set_contains_position(descent_set_mask, position) {
            positions.push(position);
        }
    }
    positions
}

/// Return the valid insertion positions for a target descent set.
///
/// A position `p ∈ [n]` is valid exactly when either:
/// - `p < n`, `p ∈ S`, and `p-1 ∉ S`, so `p` is the left boundary of a
///   descent run, or
/// - `p = n` and `n-1 ∉ S`.
pub fn valid_insertion_positions_for_target_descent_set(
    target_descent_set_mask: u64,
    n: u8,
) -> Vec<u8> {
    let mut positions = Vec::new();
    for insertion_position in 1..n {
        if target_descent_set_contains_position(target_descent_set_mask, insertion_position)
            && (insertion_position < 2
                || !target_descent_set_contains_position(
                    target_descent_set_mask,
                    insertion_position - 1,
                ))
        {
            positions.push(insertion_position);
        }
    }
    if n >= 1 && (n == 1 || !target_descent_set_contains_position(target_descent_set_mask, n - 1)) {
        positions.push(n);
    }
    positions
}

/// Return whether the given insertion position is valid for the target
/// descent set.
pub fn is_valid_insertion_position_for_target_descent_set(
    target_descent_set_mask: u64,
    insertion_position: u8,
    n: u8,
) -> bool {
    valid_insertion_positions_for_target_descent_set(target_descent_set_mask, n)
        .contains(&insertion_position)
}

/// Compute the base source descent set `S'_p`.
///
/// This is the descent set of a source permutation in `S_{n-1}` away from the
/// unconstrained position `p-1`.
///
/// # Panics
///
/// Panics if `insertion_position` is not in `[1, n]`.
pub fn base_source_descent_set_for_target_descent_set_and_insertion_position(
    target_descent_set_mask: u64,
    insertion_position: u8,
    n: u8,
) -> u64 {
    assert!(
        (1..=n).contains(&insertion_position),
        "insertion position must lie in [1, n]"
    );

    if n <= 2 {
        return 0;
    }
    if insertion_position == n {
        return target_descent_set_mask;
    }

    let mut source_descent_set_mask = 0u64;

    if insertion_position == 1 {
        for target_descent_position in 2..n {
            if target_descent_set_contains_position(
                target_descent_set_mask,
                target_descent_position,
            ) {
                source_descent_set_mask |= 1 << (target_descent_position - 2);
            }
        }
        return source_descent_set_mask;
    }

    for source_descent_position in 1..=(insertion_position.saturating_sub(2)) {
        if target_descent_set_contains_position(target_descent_set_mask, source_descent_position) {
            source_descent_set_mask |= 1 << (source_descent_position - 1);
        }
    }

    for target_descent_position in (insertion_position + 1)..n {
        if target_descent_set_contains_position(target_descent_set_mask, target_descent_position) {
            source_descent_set_mask |= 1 << (target_descent_position - 2);
        }
    }

    source_descent_set_mask
}

/// Compute the augmented source descent set `S''_p = S'_p ∪ {p-1}` when it
/// exists.
///
/// The augmented source exists exactly for `1 < p < n`.
pub fn augmented_source_descent_set_for_target_descent_set_and_insertion_position(
    target_descent_set_mask: u64,
    insertion_position: u8,
    n: u8,
) -> Option<u64> {
    if insertion_position <= 1 || insertion_position >= n {
        return None;
    }
    Some(
        base_source_descent_set_for_target_descent_set_and_insertion_position(
            target_descent_set_mask,
            insertion_position,
            n,
        ) | (1 << (insertion_position - 2)),
    )
}

/// Return the last descent position in the descent run that starts at a valid
/// insertion position.
///
/// If `p < n` is a valid insertion position, then `p` is the left boundary of
/// a descent run in the target descent set. This function returns the right
/// endpoint of that run.
pub fn descent_run_end_position_for_valid_insertion_position(
    target_descent_set_mask: u64,
    insertion_position: u8,
    n: u8,
) -> Option<u8> {
    if insertion_position >= n
        || !is_valid_insertion_position_for_target_descent_set(
            target_descent_set_mask,
            insertion_position,
            n,
        )
        || !target_descent_set_contains_position(target_descent_set_mask, insertion_position)
    {
        return None;
    }

    let mut descent_run_end_position = insertion_position;
    while descent_run_end_position + 1 < n
        && target_descent_set_contains_position(
            target_descent_set_mask,
            descent_run_end_position + 1,
        )
    {
        descent_run_end_position += 1;
    }
    Some(descent_run_end_position)
}

/// Describe the local source-set relation for a consecutive valid pair
/// `p_a < p_b`.
///
/// Returns `None` unless the two insertion positions are consecutive elements
/// of `P(S)`.
pub fn consecutive_valid_insertion_position_source_descent_set_data(
    target_descent_set_mask: u64,
    left_insertion_position: u8,
    right_insertion_position: u8,
    n: u8,
) -> Option<ConsecutiveValidInsertionPositionSourceDescentSetData> {
    if left_insertion_position >= right_insertion_position {
        return None;
    }

    let valid_positions =
        valid_insertion_positions_for_target_descent_set(target_descent_set_mask, n);
    let left_index = valid_positions
        .iter()
        .position(|&position| position == left_insertion_position)?;
    if valid_positions.get(left_index + 1).copied()? != right_insertion_position {
        return None;
    }

    let descent_run_end_position = descent_run_end_position_for_valid_insertion_position(
        target_descent_set_mask,
        left_insertion_position,
        n,
    )?;

    let mut common_fixed_part_source_descent_set_mask = 0u64;
    for descent_position in 1..left_insertion_position.saturating_sub(1) {
        if target_descent_set_contains_position(target_descent_set_mask, descent_position) {
            common_fixed_part_source_descent_set_mask |= 1 << (descent_position - 1);
        }
    }
    for target_descent_position in (right_insertion_position + 1)..n {
        if target_descent_set_contains_position(target_descent_set_mask, target_descent_position) {
            common_fixed_part_source_descent_set_mask |= 1 << (target_descent_position - 2);
        }
    }

    Some(ConsecutiveValidInsertionPositionSourceDescentSetData {
        target_permutation_size: n,
        left_insertion_position,
        right_insertion_position,
        descent_run_end_position,
        common_fixed_part_source_descent_set_mask,
        left_base_source_descent_set_mask:
            base_source_descent_set_for_target_descent_set_and_insertion_position(
                target_descent_set_mask,
                left_insertion_position,
                n,
            ),
        right_base_source_descent_set_mask:
            base_source_descent_set_for_target_descent_set_and_insertion_position(
                target_descent_set_mask,
                right_insertion_position,
                n,
            ),
        left_augmented_source_descent_set_mask:
            augmented_source_descent_set_for_target_descent_set_and_insertion_position(
                target_descent_set_mask,
                left_insertion_position,
                n,
            ),
        right_augmented_source_descent_set_mask:
            augmented_source_descent_set_for_target_descent_set_and_insertion_position(
                target_descent_set_mask,
                right_insertion_position,
                n,
            ),
    })
}

/// Return whether insertion breaks a consecutive ascending pair at the
/// insertion boundary.
///
/// This is the quantity `ε₁` in the fixed-descent insertion recurrence.
pub fn insertion_breaks_consecutive_ascending_pair_at_boundary(
    source_permutation: &[u8],
    insertion_position: u8,
) -> bool {
    let n = source_permutation.len() as u8 + 1;
    if insertion_position <= 1 || insertion_position >= n {
        return false;
    }
    source_permutation[(insertion_position - 2) as usize] + 1
        == source_permutation[(insertion_position - 1) as usize]
}

/// Return whether inserting the new maximum creates the new long swap
/// `(n-1, n)`.
///
/// The parameter `position_of_n_minus_1_in_source_permutation` is 1-indexed.
/// This is the quantity `ε₂`.
pub fn insertion_creates_new_long_swap_with_previous_maximum(
    insertion_position: u8,
    position_of_n_minus_1_in_source_permutation: u8,
) -> bool {
    insertion_position >= 2 && position_of_n_minus_1_in_source_permutation <= insertion_position - 2
}

/// Return the 1-indexed position of a value in a permutation.
pub fn one_indexed_position_of_value(permutation: &[u8], value: u8) -> Option<u8> {
    permutation
        .iter()
        .position(|&entry| entry == value)
        .map(|index| index as u8 + 1)
}

/// Compute the modified source statistic
/// `swaps(pi) + ε₁(pi,p)`.
pub fn modified_source_swaps_for_insertion_position(
    source_permutation: &[u8],
    insertion_position: u8,
) -> usize {
    compute(source_permutation, Stat::Swaps)
        + usize::from(insertion_breaks_consecutive_ascending_pair_at_boundary(
            source_permutation,
            insertion_position,
        ))
}

/// Compute the long-swaps statistic after inserting the new maximum.
///
/// This uses the insertion formula
/// `swaps(sigma) = swaps(pi) + ε₁ + ε₂`.
pub fn long_swaps_after_inserting_new_maximum(
    source_permutation: &[u8],
    insertion_position: u8,
) -> usize {
    let new_maximum_value = source_permutation.len() as u8;
    let position_of_previous_maximum =
        one_indexed_position_of_value(source_permutation, new_maximum_value)
            .expect("source permutation must contain n-1");
    modified_source_swaps_for_insertion_position(source_permutation, insertion_position)
        + usize::from(insertion_creates_new_long_swap_with_previous_maximum(
            insertion_position,
            position_of_previous_maximum,
        ))
}

/// Group all permutations of `[1..source_permutation_size]` by descent-set
/// bitmask.
///
/// This is useful when repeatedly computing source distributions at a fixed
/// size `n-1` for many target descent sets and insertion positions.
pub fn permutations_grouped_by_descent_set_bitmask(
    source_permutation_size: u8,
) -> BTreeMap<u64, Vec<Vec<u8>>> {
    let mut permutations_by_descent_set_bitmask: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
    for permutation in all_permutations(source_permutation_size) {
        permutations_by_descent_set_bitmask
            .entry(descent_set_bitmask(&permutation))
            .or_default()
            .push(permutation);
    }
    permutations_by_descent_set_bitmask
}

/// Build the `q`-refined modified source distribution from a pre-grouped
/// source-permutation table.
///
/// The input map should group permutations of `[1..n-1]` by their descent-set
/// bitmask, typically obtained from
/// [`permutations_grouped_by_descent_set_bitmask`].
pub fn q_refined_modified_source_distribution_from_grouped_source_permutations(
    target_descent_set_mask: u64,
    insertion_position: u8,
    n: u8,
    source_permutations_grouped_by_descent_set_bitmask: &BTreeMap<u64, Vec<Vec<u8>>>,
) -> QRefinedModifiedSourceDistributionData {
    let base_source_descent_set_mask =
        base_source_descent_set_for_target_descent_set_and_insertion_position(
            target_descent_set_mask,
            insertion_position,
            n,
        );
    let augmented_source_descent_set_mask =
        augmented_source_descent_set_for_target_descent_set_and_insertion_position(
            target_descent_set_mask,
            insertion_position,
            n,
        );

    let mut counts_by_modified_source_swaps_and_previous_maximum_position = BTreeMap::new();

    if let Some(source_class) =
        source_permutations_grouped_by_descent_set_bitmask.get(&base_source_descent_set_mask)
    {
        for source_permutation in source_class {
            let previous_maximum_position =
                one_indexed_position_of_value(source_permutation, n - 1)
                    .expect("source permutation must contain n-1");
            let modified_source_swaps = modified_source_swaps_for_insertion_position(
                source_permutation,
                insertion_position,
            );
            *counts_by_modified_source_swaps_and_previous_maximum_position
                .entry((modified_source_swaps, previous_maximum_position))
                .or_default() += 1;
        }
    }

    if let Some(augmented_source_descent_set_mask) = augmented_source_descent_set_mask {
        if let Some(source_class) = source_permutations_grouped_by_descent_set_bitmask
            .get(&augmented_source_descent_set_mask)
        {
            for source_permutation in source_class {
                let previous_maximum_position =
                    one_indexed_position_of_value(source_permutation, n - 1)
                        .expect("source permutation must contain n-1");
                let modified_source_swaps = compute(source_permutation, Stat::Swaps);
                *counts_by_modified_source_swaps_and_previous_maximum_position
                    .entry((modified_source_swaps, previous_maximum_position))
                    .or_default() += 1;
            }
        }
    }

    QRefinedModifiedSourceDistributionData {
        target_permutation_size: n,
        target_descent_set_mask,
        insertion_position,
        base_source_descent_set_mask,
        augmented_source_descent_set_mask,
        counts_by_modified_source_swaps_and_previous_maximum_position,
    }
}

/// Build the `q`-refined modified source distribution for `(n, S, p)` by
/// enumerating and grouping all source permutations of size `n-1`.
pub fn q_refined_modified_source_distribution_for_target_descent_set_and_insertion_position(
    target_descent_set_mask: u64,
    insertion_position: u8,
    n: u8,
) -> QRefinedModifiedSourceDistributionData {
    let source_permutations_grouped_by_descent_set_bitmask =
        permutations_grouped_by_descent_set_bitmask(n.saturating_sub(1));
    q_refined_modified_source_distribution_from_grouped_source_permutations(
        target_descent_set_mask,
        insertion_position,
        n,
        &source_permutations_grouped_by_descent_set_bitmask,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::permutation::all_permutations;
    use crate::statistics::descent_set_bitmask;

    fn positions_to_mask(positions: &[u8]) -> u64 {
        let mut mask = 0u64;
        for &position in positions {
            mask |= 1 << (position - 1);
        }
        mask
    }

    #[test]
    fn valid_insertion_positions_match_examples() {
        assert_eq!(
            valid_insertion_positions_for_target_descent_set(positions_to_mask(&[2, 5]), 6),
            vec![2, 5]
        );
        assert_eq!(
            valid_insertion_positions_for_target_descent_set(positions_to_mask(&[2, 3, 6]), 7),
            vec![2, 6]
        );
        assert_eq!(
            valid_insertion_positions_for_target_descent_set(positions_to_mask(&[]), 4),
            vec![4]
        );
    }

    #[test]
    fn source_descent_sets_match_paper_examples() {
        let target_descent_set_mask = positions_to_mask(&[2, 5]);
        assert_eq!(
            descent_positions_from_bitmask(
                base_source_descent_set_for_target_descent_set_and_insertion_position(
                    target_descent_set_mask,
                    2,
                    6
                ),
                5
            ),
            vec![4]
        );
        assert_eq!(
            descent_positions_from_bitmask(
                base_source_descent_set_for_target_descent_set_and_insertion_position(
                    target_descent_set_mask,
                    5,
                    6
                ),
                5
            ),
            vec![2]
        );
        assert_eq!(
            descent_positions_from_bitmask(
                augmented_source_descent_set_for_target_descent_set_and_insertion_position(
                    target_descent_set_mask,
                    2,
                    6
                )
                .unwrap(),
                5
            ),
            vec![1, 4]
        );
    }

    #[test]
    fn consecutive_valid_pair_local_structure_matches_theory() {
        let target_descent_set_mask = positions_to_mask(&[2, 3, 5]);
        let data = consecutive_valid_insertion_position_source_descent_set_data(
            target_descent_set_mask,
            2,
            5,
            6,
        )
        .unwrap();

        assert_eq!(data.descent_run_end_position, 3);
        assert_eq!(
            descent_positions_from_bitmask(data.left_base_source_descent_set_mask, 5),
            vec![2, 4]
        );
        assert_eq!(
            descent_positions_from_bitmask(data.right_base_source_descent_set_mask, 5),
            vec![2, 3]
        );
        assert_eq!(data.left_only_base_source_positions(), vec![4]);
        assert_eq!(data.right_only_base_source_positions(), vec![3]);
    }

    #[test]
    fn insertion_formula_matches_direct_computation() {
        for n in 2..=7 {
            for source_permutation in all_permutations(n - 1) {
                for insertion_position in 1..=n {
                    let mut target_permutation = Vec::with_capacity(n as usize);
                    for &value in source_permutation
                        .iter()
                        .take((insertion_position - 1) as usize)
                    {
                        target_permutation.push(value);
                    }
                    target_permutation.push(n);
                    for &value in source_permutation
                        .iter()
                        .skip((insertion_position - 1) as usize)
                    {
                        target_permutation.push(value);
                    }

                    assert_eq!(
                        long_swaps_after_inserting_new_maximum(
                            &source_permutation,
                            insertion_position
                        ),
                        compute(&target_permutation, Stat::Swaps),
                    );
                }
            }
        }
    }

    #[test]
    fn source_descent_sets_match_direct_insertion_filtering() {
        for n in 3..=7 {
            for target_descent_set_mask in 0u64..(1 << (n - 1)) {
                let valid_positions =
                    valid_insertion_positions_for_target_descent_set(target_descent_set_mask, n);
                for insertion_position in valid_positions {
                    let mut observed_source_descent_sets = Vec::new();
                    for source_permutation in all_permutations(n - 1) {
                        let mut target_permutation = Vec::with_capacity(n as usize);
                        for &value in source_permutation
                            .iter()
                            .take((insertion_position - 1) as usize)
                        {
                            target_permutation.push(value);
                        }
                        target_permutation.push(n);
                        for &value in source_permutation
                            .iter()
                            .skip((insertion_position - 1) as usize)
                        {
                            target_permutation.push(value);
                        }

                        if descent_set_bitmask(&target_permutation) == target_descent_set_mask {
                            observed_source_descent_sets
                                .push(descent_set_bitmask(&source_permutation));
                        }
                    }

                    observed_source_descent_sets.sort_unstable();
                    observed_source_descent_sets.dedup();

                    let mut expected_source_descent_sets = vec![
                        base_source_descent_set_for_target_descent_set_and_insertion_position(
                            target_descent_set_mask,
                            insertion_position,
                            n,
                        ),
                    ];
                    if let Some(augmented_source_descent_set_mask) =
                        augmented_source_descent_set_for_target_descent_set_and_insertion_position(
                            target_descent_set_mask,
                            insertion_position,
                            n,
                        )
                    {
                        expected_source_descent_sets.push(augmented_source_descent_set_mask);
                    }
                    expected_source_descent_sets.sort_unstable();
                    expected_source_descent_sets.dedup();

                    assert_eq!(observed_source_descent_sets, expected_source_descent_sets);
                }
            }
        }
    }

    #[test]
    fn grouped_source_permutations_match_direct_grouping() {
        let grouped = permutations_grouped_by_descent_set_bitmask(4);
        let mut direct: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for permutation in all_permutations(4) {
            direct
                .entry(descent_set_bitmask(&permutation))
                .or_default()
                .push(permutation);
        }
        assert_eq!(grouped, direct);
    }

    #[test]
    fn q_refined_distribution_reconstructs_output_polynomial() {
        for n in 3..=7 {
            let grouped_source_permutations = permutations_grouped_by_descent_set_bitmask(n - 1);
            for target_descent_set_mask in 0u64..(1 << (n - 1)) {
                for insertion_position in
                    valid_insertion_positions_for_target_descent_set(target_descent_set_mask, n)
                {
                    let distribution =
                        q_refined_modified_source_distribution_from_grouped_source_permutations(
                            target_descent_set_mask,
                            insertion_position,
                            n,
                            &grouped_source_permutations,
                        );

                    let mut direct_polynomial = Vec::<i64>::new();
                    for source_permutation in all_permutations(n - 1) {
                        let mut target_permutation = Vec::with_capacity(n as usize);
                        for &value in source_permutation
                            .iter()
                            .take((insertion_position - 1) as usize)
                        {
                            target_permutation.push(value);
                        }
                        target_permutation.push(n);
                        for &value in source_permutation
                            .iter()
                            .skip((insertion_position - 1) as usize)
                        {
                            target_permutation.push(value);
                        }

                        if descent_set_bitmask(&target_permutation) != target_descent_set_mask {
                            continue;
                        }

                        let output_long_swaps = compute(&target_permutation, Stat::Swaps);
                        if direct_polynomial.len() <= output_long_swaps {
                            direct_polynomial.resize(output_long_swaps + 1, 0);
                        }
                        direct_polynomial[output_long_swaps] += 1;
                    }

                    if direct_polynomial.is_empty() {
                        direct_polynomial.push(0);
                    }
                    while direct_polynomial.len() > 1 && *direct_polynomial.last().unwrap() == 0 {
                        direct_polynomial.pop();
                    }

                    assert_eq!(
                        distribution.output_long_swaps_polynomial_after_insertion(),
                        direct_polynomial,
                    );
                }
            }
        }
    }

    #[test]
    fn direct_and_grouped_q_refined_distribution_builders_agree() {
        for n in 3..=6 {
            for target_descent_set_mask in 0u64..(1 << (n - 1)) {
                for insertion_position in
                    valid_insertion_positions_for_target_descent_set(target_descent_set_mask, n)
                {
                    let direct =
                        q_refined_modified_source_distribution_for_target_descent_set_and_insertion_position(
                            target_descent_set_mask,
                            insertion_position,
                            n,
                        );
                    let grouped =
                        q_refined_modified_source_distribution_from_grouped_source_permutations(
                            target_descent_set_mask,
                            insertion_position,
                            n,
                            &permutations_grouped_by_descent_set_bitmask(n - 1),
                        );
                    assert_eq!(direct, grouped);
                }
            }
        }
    }
}
