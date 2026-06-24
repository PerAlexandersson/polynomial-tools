//! Closed meanders and rooted meandric permutations.
//!
//! A closed meandric system of order `n` is represented here by an ordered pair
//! of noncrossing perfect matchings on the zero-indexed ground set
//! `{0, ..., 2n - 1}`.  The pair is connected, hence a meander, exactly when
//! alternating upper and lower arcs gives a single cycle.
//!
//! The emitted meandric permutations use the standard one-indexed road labels
//! `{1, ..., 2n}`.  We root the cyclic order at label `1` and take the first
//! step along an upper arc.  We do not quotient by reversal, reflection, or
//! cyclic relabeling.  With this convention the counts are OEIS A005315.

/// A noncrossing perfect matching stored as a partner table.
///
/// If `matching[i] == j`, then `matching[j] == i`.  The labels are zero-indexed.
pub type NoncrossingMatching = Vec<usize>;

/// A rooted meandric permutation in one-line cyclic-word notation.
///
/// The entries are one-indexed road labels.  The empty vector is the unique
/// order-zero meander.
pub type MeandricPermutation = Vec<usize>;

/// Initial values of OEIS A005315, the closed meandric numbers.
pub const CLOSED_MEANDRIC_NUMBERS_INITIAL: &[usize] = &[
    1, 1, 2, 8, 42, 262, 1828, 13820, 110954, 933458, 8152860, 73424650,
];

/// Return all noncrossing perfect matchings on `{0, ..., 2n - 1}`.
pub fn noncrossing_perfect_matchings(n: usize) -> Vec<NoncrossingMatching> {
    let mut memo = vec![None; n + 1];
    noncrossing_perfect_matchings_cached(n, &mut memo)
}

fn noncrossing_perfect_matchings_cached(
    n: usize,
    memo: &mut [Option<Vec<NoncrossingMatching>>],
) -> Vec<NoncrossingMatching> {
    if let Some(matchings) = &memo[n] {
        return matchings.clone();
    }

    let matchings = if n == 0 {
        vec![Vec::new()]
    } else {
        let size = 2 * n;
        let mut result = Vec::new();

        for partner_zero in (1..size).step_by(2) {
            let left_arcs = (partner_zero - 1) / 2;
            let right_arcs = n - 1 - left_arcs;
            let right_shift = partner_zero + 1;

            let left_matchings = noncrossing_perfect_matchings_cached(left_arcs, memo);
            let right_matchings = noncrossing_perfect_matchings_cached(right_arcs, memo);

            for left in &left_matchings {
                for right in &right_matchings {
                    let mut matching = vec![usize::MAX; size];
                    matching[0] = partner_zero;
                    matching[partner_zero] = 0;

                    for (index, &partner) in left.iter().enumerate() {
                        matching[1 + index] = 1 + partner;
                    }
                    for (index, &partner) in right.iter().enumerate() {
                        matching[right_shift + index] = right_shift + partner;
                    }

                    result.push(matching);
                }
            }
        }

        result
    };

    memo[n] = Some(matchings.clone());
    matchings
}

/// Return the rooted cyclic order of a connected ordered arch pair.
///
/// The input matchings use zero-indexed partner tables.  The returned
/// permutation uses one-indexed road labels, starts at label `1`, and first
/// follows an upper arc.  Disconnected arch pairs return `None`.
pub fn rooted_meandric_permutation_from_arch_pair(
    upper: &[usize],
    lower: &[usize],
) -> Option<MeandricPermutation> {
    assert_eq!(
        upper.len(),
        lower.len(),
        "upper and lower matchings must have the same size"
    );
    if upper.is_empty() {
        return Some(Vec::new());
    }

    let size = upper.len();
    let mut vertex = 0usize;
    let mut use_upper = true;
    let mut seen = vec![false; size];
    let mut word = Vec::with_capacity(size);

    for _ in 0..size {
        if seen[vertex] {
            return None;
        }
        seen[vertex] = true;
        word.push(vertex + 1);
        vertex = if use_upper {
            upper[vertex]
        } else {
            lower[vertex]
        };
        use_upper = !use_upper;
    }

    (vertex == 0).then_some(word)
}

/// Check whether an ordered arch pair is connected.
pub fn is_connected_arch_pair(upper: &[usize], lower: &[usize]) -> bool {
    rooted_meandric_permutation_from_arch_pair(upper, lower).is_some()
}

/// Iterator over rooted meandric permutations of order `n`.
#[derive(Debug, Clone)]
pub struct RootedMeandricPermutations {
    matchings: Vec<NoncrossingMatching>,
    upper_index: usize,
    lower_index: usize,
}

impl RootedMeandricPermutations {
    /// Create an iterator over rooted meandric permutations of order `n`.
    pub fn new(n: usize) -> Self {
        Self {
            matchings: noncrossing_perfect_matchings(n),
            upper_index: 0,
            lower_index: 0,
        }
    }

    /// Number of ordered arch pairs scanned by this generator.
    pub fn arch_pair_count(&self) -> usize {
        self.matchings.len() * self.matchings.len()
    }

    fn advance_pair(&mut self) -> Option<(usize, usize)> {
        if self.upper_index >= self.matchings.len() {
            return None;
        }

        let pair = (self.upper_index, self.lower_index);
        self.lower_index += 1;
        if self.lower_index == self.matchings.len() {
            self.lower_index = 0;
            self.upper_index += 1;
        }
        Some(pair)
    }
}

impl Iterator for RootedMeandricPermutations {
    type Item = MeandricPermutation;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((upper_index, lower_index)) = self.advance_pair() {
            let upper = &self.matchings[upper_index];
            let lower = &self.matchings[lower_index];
            if let Some(permutation) = rooted_meandric_permutation_from_arch_pair(upper, lower) {
                return Some(permutation);
            }
        }
        None
    }
}

/// Return an iterator over rooted meandric permutations of order `n`.
pub fn rooted_meandric_permutations(n: usize) -> RootedMeandricPermutations {
    RootedMeandricPermutations::new(n)
}

/// Count rooted meandric permutations of order `n`.
pub fn rooted_meandric_permutation_count(n: usize) -> usize {
    rooted_meandric_permutations(n).count()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn noncrossing_matching_counts_are_catalan() {
        let catalan = [1usize, 1, 2, 5, 14, 42, 132, 429, 1430];
        for (n, &expected) in catalan.iter().enumerate() {
            assert_eq!(noncrossing_perfect_matchings(n).len(), expected);
        }
    }

    #[test]
    fn order_two_rooted_meandric_permutations() {
        let permutations: Vec<_> = rooted_meandric_permutations(2).collect();
        assert_eq!(permutations, vec![vec![1, 2, 3, 4], vec![1, 4, 3, 2]]);
    }

    #[test]
    fn disconnected_arch_pair_returns_none() {
        let upper = vec![1, 0, 3, 2];
        let lower = vec![1, 0, 3, 2];
        assert_eq!(
            rooted_meandric_permutation_from_arch_pair(&upper, &lower),
            None
        );
        assert!(!is_connected_arch_pair(&upper, &lower));
    }

    #[test]
    fn rooted_meandric_counts_match_a005315() {
        for (n, &expected) in CLOSED_MEANDRIC_NUMBERS_INITIAL.iter().take(9).enumerate() {
            assert_eq!(rooted_meandric_permutation_count(n), expected, "n={n}");
        }
    }

    #[test]
    fn iterator_reports_arch_pair_count() {
        let iterator = RootedMeandricPermutations::new(5);
        assert_eq!(iterator.arch_pair_count(), 42 * 42);
    }
}
