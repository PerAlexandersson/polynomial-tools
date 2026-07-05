//! Basic permutation utilities shared across combinatorial crates.

/// Return whether `perm` is a permutation of `1..=n` in one-line notation.
pub fn is_one_indexed_permutation(perm: &[usize]) -> bool {
    let n = perm.len();
    let mut seen = vec![false; n];
    for &value in perm {
        if !(1..=n).contains(&value) || seen[value - 1] {
            return false;
        }
        seen[value - 1] = true;
    }
    true
}

/// Panic unless `perm` is a permutation of `1..=n` in one-line notation.
pub fn assert_one_indexed_permutation(perm: &[usize]) {
    assert!(
        is_one_indexed_permutation(perm),
        "expected a one-indexed permutation of 1..={}",
        perm.len()
    );
}

/// Advance `perm` to the next lexicographic permutation.
///
/// Returns `false` if `perm` is already the final permutation in lexicographic
/// order, leaving it unchanged.
pub fn next_permutation<T: Ord>(perm: &mut [T]) -> bool {
    let n = perm.len();
    if n <= 1 {
        return false;
    }

    let mut i = n - 2;
    while perm[i] >= perm[i + 1] {
        if i == 0 {
            return false;
        }
        i -= 1;
    }

    let mut j = n - 1;
    while perm[j] <= perm[i] {
        j -= 1;
    }

    perm.swap(i, j);
    perm[i + 1..].reverse();
    true
}

/// Return all permutations of `{0, ..., n-1}` in lexicographic order.
pub fn all_permutations_zero_indexed(n: usize) -> Vec<Vec<usize>> {
    if n == 0 {
        return vec![vec![]];
    }

    let mut result = Vec::new();
    let mut perm: Vec<usize> = (0..n).collect();
    loop {
        result.push(perm.clone());
        if !next_permutation(&mut perm) {
            break;
        }
    }
    result
}

/// Return all permutations of `{1, ..., n}` in lexicographic order.
pub fn all_permutations_one_indexed(n: usize) -> Vec<Vec<usize>> {
    if n == 0 {
        return vec![vec![]];
    }

    let mut result = Vec::new();
    let mut perm: Vec<usize> = (1..=n).collect();
    loop {
        result.push(perm.clone());
        if !next_permutation(&mut perm) {
            break;
        }
    }
    result
}

/// Return the longest permutation in `S_n` in one-line notation.
pub fn longest_permutation(n: usize) -> Vec<usize> {
    (1..=n).rev().collect()
}

/// Return the inverse of a permutation in one-line notation.
///
/// Panics if `perm` is not a permutation of `1..=n`.
pub fn inverse_permutation(perm: &[usize]) -> Vec<usize> {
    let n = perm.len();
    let mut inv = vec![0usize; n];
    let mut seen = vec![false; n];
    for (i, &v) in perm.iter().enumerate() {
        assert!(
            (1..=n).contains(&v),
            "permutation entry {v} is out of range for size {n}"
        );
        assert!(!seen[v - 1], "permutation entry {v} occurs more than once");
        seen[v - 1] = true;
        inv[v - 1] = i + 1;
    }
    inv
}

/// Compose permutations in one-line notation: return `a ∘ b`.
pub fn compose_permutations(a: &[usize], b: &[usize]) -> Vec<usize> {
    assert_eq!(a.len(), b.len());
    (0..a.len()).map(|i| a[b[i] - 1]).collect()
}

/// Compute the permutation represented by a word of simple transpositions.
///
/// The word `[a_1, ..., a_k]` represents `s_{a_1} ... s_{a_k}` in the standard
/// 1-indexed Coxeter generators of `S_n`, applied right-to-left.
pub fn permutation_from_simple_transpositions(n: usize, word: &[usize]) -> Vec<usize> {
    let mut perm: Vec<usize> = (1..=n).collect();
    for &s in word.iter().rev() {
        assert!(
            (1..n).contains(&s),
            "simple transposition s_{s} is out of range for S_{n}"
        );
        perm.swap(s - 1, s);
    }
    perm
}

/// Compute a reduced word using the leftmost-descent bubble-sort algorithm.
///
/// The result is 0-indexed: entry `i` means the simple transposition `s_i`
/// swapping positions `i` and `i+1`.
pub fn reduced_word(perm: &[usize]) -> Vec<usize> {
    let mut pi = perm.to_vec();
    let mut word = Vec::new();
    loop {
        let mut found = false;
        for i in 0..pi.len().saturating_sub(1) {
            if pi[i] > pi[i + 1] {
                word.push(i);
                pi.swap(i, i + 1);
                found = true;
                break;
            }
        }
        if !found {
            break;
        }
    }
    word
}

fn optimist_sort_step_unchecked(
    perm: &mut [usize],
    positions: &mut Vec<usize>,
    values: &mut Vec<usize>,
) -> bool {
    positions.clear();
    positions.extend((0..perm.len()).filter(|&i| perm[i] != i + 1));
    if positions.is_empty() {
        return false;
    }

    let target_position = perm[positions[0]] - 1;
    let shift = positions
        .iter()
        .position(|&position| position == target_position)
        .expect("target position of an unfixed value should be unfixed");

    values.clear();
    values.extend(positions.iter().map(|&position| perm[position]));
    let len = positions.len();
    for (i, &position) in positions.iter().enumerate() {
        perm[position] = values[(i + len - shift) % len];
    }
    true
}

/// Apply one optimist-sorting step to a one-indexed permutation.
///
/// The operation follows the OEIS A345453 convention.  Look at the unfixed
/// positions, take the first one, and rotate the values in all unfixed
/// positions so that the value at the first unfixed position moves to its
/// target position.
///
/// Returns `true` if the permutation changed, and `false` if it was already
/// sorted.
pub fn optimist_sort_step(perm: &mut [usize]) -> bool {
    assert_one_indexed_permutation(perm);
    let mut positions = Vec::with_capacity(perm.len());
    let mut values = Vec::with_capacity(perm.len());
    optimist_sort_step_unchecked(perm, &mut positions, &mut values)
}

fn unfixed_standardization_unchecked(
    perm: &[usize],
    standardized: &mut Vec<usize>,
    positions: &mut Vec<usize>,
    rank_by_value: &mut Vec<usize>,
) {
    positions.clear();
    positions.extend((0..perm.len()).filter(|&i| perm[i] != i + 1));

    standardized.clear();
    if positions.is_empty() {
        return;
    }

    rank_by_value.clear();
    rank_by_value.resize(perm.len() + 1, 0);
    for (rank, &position) in positions.iter().enumerate() {
        rank_by_value[position + 1] = rank + 1;
    }
    standardized.extend(
        positions
            .iter()
            .map(|&position| rank_by_value[perm[position]]),
    );
}

/// Delete fixed points from a one-indexed permutation and standardize the rest.
///
/// The remaining positions and values are both the same subset of `1..=n`.
/// This returns the relative one-indexed permutation on that subset.  It is the
/// reduced state seen by the optimist-sorting process after ignoring positions
/// that are already correct.
pub fn unfixed_standardization(perm: &[usize]) -> Vec<usize> {
    assert_one_indexed_permutation(perm);
    let mut standardized = Vec::with_capacity(perm.len());
    let mut positions = Vec::with_capacity(perm.len());
    let mut rank_by_value = Vec::with_capacity(perm.len() + 1);
    unfixed_standardization_unchecked(perm, &mut standardized, &mut positions, &mut rank_by_value);
    standardized
}

fn optimist_sort_steps_reduced_unchecked_with_scratch(
    perm: &[usize],
    state: &mut Vec<usize>,
    work: &mut Vec<usize>,
    next_state: &mut Vec<usize>,
    positions: &mut Vec<usize>,
    values: &mut Vec<usize>,
    rank_by_value: &mut Vec<usize>,
) -> usize {
    unfixed_standardization_unchecked(perm, state, positions, rank_by_value);

    let mut steps = 0;
    while !state.is_empty() {
        work.clear();
        work.extend_from_slice(state);
        optimist_sort_step_unchecked(work, positions, values);
        unfixed_standardization_unchecked(work, next_state, positions, rank_by_value);
        std::mem::swap(state, next_state);
        steps += 1;
    }
    steps
}

/// Return the number of optimist-sorting steps needed to sort a permutation.
///
/// The input is a one-indexed permutation in one-line notation.
pub fn optimist_sort_steps(perm: &[usize]) -> usize {
    assert_one_indexed_permutation(perm);
    let mut state = Vec::with_capacity(perm.len());
    let mut work = Vec::with_capacity(perm.len());
    let mut next_state = Vec::with_capacity(perm.len());
    let mut positions = Vec::with_capacity(perm.len());
    let mut values = Vec::with_capacity(perm.len());
    let mut rank_by_value = Vec::with_capacity(perm.len() + 1);
    optimist_sort_steps_reduced_unchecked_with_scratch(
        perm,
        &mut state,
        &mut work,
        &mut next_state,
        &mut positions,
        &mut values,
        &mut rank_by_value,
    )
}

fn binomial_u128(n: usize, k: usize) -> u128 {
    let k = k.min(n - k);
    let mut result = 1u128;
    for i in 1..=k {
        result = result * (n - k + i) as u128 / i as u128;
    }
    result
}

fn derangement_distribution_backtrack(
    position: usize,
    current: &mut [usize],
    used: &mut [bool],
    counts: &mut [u128],
    state: &mut Vec<usize>,
    work: &mut Vec<usize>,
    next_state: &mut Vec<usize>,
    positions: &mut Vec<usize>,
    values: &mut Vec<usize>,
    rank_by_value: &mut Vec<usize>,
) {
    let n = current.len();
    if position == n {
        let steps = optimist_sort_steps_reduced_unchecked_with_scratch(
            current,
            state,
            work,
            next_state,
            positions,
            values,
            rank_by_value,
        );
        counts[steps] += 1;
        return;
    }

    for value in 1..=n {
        if value == position + 1 || used[value] {
            continue;
        }
        current[position] = value;
        used[value] = true;
        derangement_distribution_backtrack(
            position + 1,
            current,
            used,
            counts,
            state,
            work,
            next_state,
            positions,
            values,
            rank_by_value,
        );
        used[value] = false;
    }
}

/// Count derangements of `1..=n` by the number of optimist-sorting steps.
///
/// This is useful because a general permutation first reduces to the
/// derangement induced by its unfixed positions.
pub fn optimist_sort_derangement_step_distribution(n: usize) -> Vec<u128> {
    if n == 0 {
        return vec![1];
    }

    let mut counts = vec![0u128; n];
    let mut current = vec![0usize; n];
    let mut used = vec![false; n + 1];
    let mut state = Vec::with_capacity(n);
    let mut work = Vec::with_capacity(n);
    let mut next_state = Vec::with_capacity(n);
    let mut positions = Vec::with_capacity(n);
    let mut values = Vec::with_capacity(n);
    let mut rank_by_value = Vec::with_capacity(n + 1);

    derangement_distribution_backtrack(
        0,
        &mut current,
        &mut used,
        &mut counts,
        &mut state,
        &mut work,
        &mut next_state,
        &mut positions,
        &mut values,
        &mut rank_by_value,
    );
    counts
}

/// Count permutations of `1..=n` by optimist-sorting steps via derangements.
///
/// If a permutation has exactly `r` unfixed positions, its reduced state is a
/// derangement of size `r`, and the choice of the unfixed positions contributes
/// the factor `binomial(n,r)`.
pub fn optimist_sort_step_distribution_via_derangements(n: usize) -> Vec<u128> {
    let mut counts = vec![0u128; n.max(1)];
    for r in 0..=n {
        let scale = binomial_u128(n, r);
        for (steps, derangement_count) in optimist_sort_derangement_step_distribution(r)
            .into_iter()
            .enumerate()
        {
            counts[steps] += scale * derangement_count;
        }
    }
    counts
}

/// Count permutations of `1..=n` by the number of optimist-sorting steps.
///
/// The returned vector has length `max(n, 1)`, so entry `k` is the number of
/// permutations requiring exactly `k` steps.  This is the row convention used
/// by OEIS A345453 for `n >= 1`.
pub fn optimist_sort_step_distribution(n: usize) -> Vec<u128> {
    if n == 0 {
        return vec![1];
    }

    let mut counts = vec![0u128; n];
    let mut perm: Vec<usize> = (1..=n).collect();
    let mut state = Vec::with_capacity(n);
    let mut work = Vec::with_capacity(n);
    let mut next_state = Vec::with_capacity(n);
    let mut positions = Vec::with_capacity(n);
    let mut values = Vec::with_capacity(n);
    let mut rank_by_value = Vec::with_capacity(n + 1);
    loop {
        let steps = optimist_sort_steps_reduced_unchecked_with_scratch(
            &perm,
            &mut state,
            &mut work,
            &mut next_state,
            &mut positions,
            &mut values,
            &mut rank_by_value,
        );
        counts[steps] += 1;
        if !next_permutation(&mut perm) {
            break;
        }
    }
    counts
}

/// Stable standardization of a word.
///
/// The result is a one-indexed permutation: equal letters are ordered by their
/// original positions.  For example, `["b", "a", "b"]` standardizes to
/// `[2, 1, 3]`.
pub fn stable_standardization<T: Ord>(word: &[T]) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..word.len()).collect();
    indices.sort_by(|&i, &j| word[i].cmp(&word[j]).then_with(|| i.cmp(&j)));

    let mut standardized = vec![0usize; word.len()];
    for (rank, index) in indices.into_iter().enumerate() {
        standardized[index] = rank + 1;
    }
    standardized
}

/// Return the number of optimist-sorting steps for a word.
///
/// This uses stable standardization, so repeated letters are distinguished by
/// their original left-to-right order before applying the permutation sorting
/// rule.
pub fn optimist_sort_steps_word<T: Ord>(word: &[T]) -> usize {
    optimist_sort_steps(&stable_standardization(word))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_next_permutation() {
        let mut perm = vec![0, 1, 2];
        assert!(next_permutation(&mut perm));
        assert_eq!(perm, vec![0, 2, 1]);
    }

    #[test]
    fn test_next_permutation_last() {
        let mut perm = vec![2, 1, 0];
        assert!(!next_permutation(&mut perm));
        assert_eq!(perm, vec![2, 1, 0]);
    }

    #[test]
    fn test_all_permutations_zero_indexed() {
        assert_eq!(all_permutations_zero_indexed(0), vec![Vec::<usize>::new()]);
        assert_eq!(
            all_permutations_zero_indexed(3),
            vec![
                vec![0, 1, 2],
                vec![0, 2, 1],
                vec![1, 0, 2],
                vec![1, 2, 0],
                vec![2, 0, 1],
                vec![2, 1, 0],
            ]
        );
    }

    #[test]
    fn test_all_permutations_one_indexed() {
        assert_eq!(all_permutations_one_indexed(0), vec![Vec::<usize>::new()]);
        assert_eq!(
            all_permutations_one_indexed(3),
            vec![
                vec![1, 2, 3],
                vec![1, 3, 2],
                vec![2, 1, 3],
                vec![2, 3, 1],
                vec![3, 1, 2],
                vec![3, 2, 1],
            ]
        );
    }

    #[test]
    fn test_longest_permutation() {
        assert_eq!(longest_permutation(0), Vec::<usize>::new());
        assert_eq!(longest_permutation(4), vec![4, 3, 2, 1]);
    }

    #[test]
    fn test_inverse_permutation() {
        assert!(is_one_indexed_permutation(&[3, 1, 4, 2]));
        assert!(!is_one_indexed_permutation(&[3, 1, 4, 4]));
        assert_eq!(inverse_permutation(&[3, 1, 4, 2]), vec![2, 4, 1, 3]);
    }

    #[test]
    fn test_compose_permutations() {
        let a = vec![2, 3, 1];
        let b = vec![3, 1, 2];
        assert_eq!(compose_permutations(&a, &b), vec![1, 2, 3]);
    }

    #[test]
    fn test_permutation_from_simple_transpositions() {
        assert_eq!(
            permutation_from_simple_transpositions(4, &[1, 3, 2]),
            vec![3, 1, 4, 2]
        );
    }

    #[test]
    fn test_reduced_word() {
        assert_eq!(reduced_word(&[1, 2, 3]), Vec::<usize>::new());
        assert_eq!(reduced_word(&[3, 2, 1]), vec![0, 1, 0]);
    }

    #[test]
    fn test_optimist_sort_step() {
        let mut perm = vec![2, 3, 1, 4];
        assert!(optimist_sort_step(&mut perm));
        assert_eq!(perm, vec![1, 2, 3, 4]);

        assert!(!optimist_sort_step(&mut perm));
        assert_eq!(perm, vec![1, 2, 3, 4]);
    }

    #[test]
    fn test_optimist_sort_steps() {
        assert_eq!(optimist_sort_steps(&[1, 2, 3, 4]), 0);
        assert_eq!(optimist_sort_steps(&[2, 3, 1, 4]), 1);
        assert_eq!(optimist_sort_steps(&[4, 3, 2, 1]), 2);
    }

    #[test]
    fn test_unfixed_standardization() {
        assert_eq!(unfixed_standardization(&[1, 3, 2, 4]), vec![2, 1]);
        assert_eq!(unfixed_standardization(&[2, 1, 4, 3]), vec![2, 1, 4, 3]);
        assert_eq!(unfixed_standardization(&[1, 2, 3]), Vec::<usize>::new());
    }

    #[test]
    fn test_optimist_sort_step_distribution_matches_a345453() {
        let expected = [
            vec![1],
            vec![1, 1],
            vec![1, 5, 0],
            vec![1, 17, 6, 0],
            vec![1, 49, 64, 6, 0],
            vec![1, 129, 432, 150, 8, 0],
            vec![1, 321, 2356, 2016, 336, 10, 0],
            vec![1, 769, 11340, 19868, 7564, 764, 14, 0],
        ];

        for (n, row) in expected.iter().enumerate() {
            let expected_row: Vec<u128> = row.iter().map(|&value| value as u128).collect();
            assert_eq!(optimist_sort_step_distribution(n + 1), expected_row);
            assert_eq!(
                optimist_sort_step_distribution_via_derangements(n + 1),
                expected_row
            );
        }
    }

    #[test]
    fn test_stable_standardization_and_word_steps() {
        assert_eq!(stable_standardization(&["b", "a", "b"]), vec![2, 1, 3]);
        assert_eq!(optimist_sort_steps_word(&["b", "a", "b"]), 1);
        assert_eq!(optimist_sort_steps_word(&["a", "b", "b"]), 0);
    }
}
