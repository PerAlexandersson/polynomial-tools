//! Basic permutation utilities shared across combinatorial crates.

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
    fn test_longest_permutation() {
        assert_eq!(longest_permutation(0), Vec::<usize>::new());
        assert_eq!(longest_permutation(4), vec![4, 3, 2, 1]);
    }

    #[test]
    fn test_inverse_permutation() {
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
}
