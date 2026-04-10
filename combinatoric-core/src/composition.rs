use crate::partition::Partition;
use std::collections::BTreeSet;
use std::fmt;
use std::hash::Hash;

// =========================================================================
// WeakComposition
// =========================================================================

/// A weak integer composition: an ordered tuple of non-negative integers
/// with a fixed number of parts (including trailing zeros).
///
/// Unlike [`Composition`], trailing zeros are **not** stripped — the length
/// encodes the number of variables, which is structural for polynomial bases
/// (key, atom, slide, etc.).
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct WeakComposition(Vec<u32>);

impl WeakComposition {
    /// Create a weak composition from parts. Length is preserved exactly.
    pub fn new(parts: Vec<u32>) -> Self {
        WeakComposition(parts)
    }

    /// Create from a slice.
    pub fn from_slice(s: &[u32]) -> Self {
        WeakComposition(s.to_vec())
    }

    /// The zero weak composition with `k` parts (all zeros).
    pub fn zero(k: usize) -> Self {
        WeakComposition(vec![0; k])
    }

    /// The parts as a slice.
    pub fn parts(&self) -> &[u32] {
        &self.0
    }

    /// Number of parts (= number of variables).
    pub fn num_vars(&self) -> usize {
        self.0.len()
    }

    /// Sum of all parts.
    pub fn degree(&self) -> u32 {
        self.0.iter().sum()
    }

    /// Whether all parts are zero.
    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|&p| p == 0)
    }

    /// Whether the parts are weakly decreasing (dominant composition).
    pub fn is_dominant(&self) -> bool {
        self.0.windows(2).all(|w| w[0] >= w[1])
    }

    /// Remove zero parts to get a (strong) composition.
    pub fn flat(&self) -> Composition {
        let parts: Vec<u32> = self.0.iter().copied().filter(|&p| p > 0).collect();
        Composition::new(parts)
    }

    /// Sort non-zero parts into a partition (the "content" of the composition).
    pub fn content(&self) -> Partition {
        let mut parts: Vec<u32> = self.0.iter().copied().filter(|&p| p > 0).collect();
        parts.sort_unstable_by(|a, b| b.cmp(a));
        Partition::from_sorted(parts)
    }

    /// Convert to a `Composition` by stripping trailing zeros.
    pub fn to_composition(&self) -> Composition {
        Composition::new(self.0.clone())
    }

    /// Create from a `Composition` by padding with zeros to `num_vars` parts.
    /// Panics if `num_vars < c.num_parts()`.
    pub fn from_composition(c: &Composition, num_vars: usize) -> Self {
        assert!(
            num_vars >= c.num_parts(),
            "num_vars ({}) < composition length ({})",
            num_vars,
            c.num_parts()
        );
        let mut parts = c.parts().to_vec();
        parts.resize(num_vars, 0);
        WeakComposition(parts)
    }

    /// All weak compositions of `d` into exactly `k` parts (parts >= 0).
    pub fn all_weak_compositions(d: u32, k: usize) -> Vec<Self> {
        if k == 0 {
            return if d == 0 {
                vec![WeakComposition(vec![])]
            } else {
                vec![]
            };
        }
        let mut results = Vec::new();
        Self::enumerate_helper(d, k, &mut Vec::with_capacity(k), &mut results);
        results
    }

    fn enumerate_helper(
        remaining: u32,
        parts_left: usize,
        current: &mut Vec<u32>,
        results: &mut Vec<WeakComposition>,
    ) {
        if parts_left == 1 {
            let mut c = current.clone();
            c.push(remaining);
            results.push(WeakComposition(c));
            return;
        }
        for part in 0..=remaining {
            current.push(part);
            Self::enumerate_helper(remaining - part, parts_left - 1, current, results);
            current.pop();
        }
    }
}

impl fmt::Display for WeakComposition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s: Vec<String> = self.0.iter().map(|p| p.to_string()).collect();
        write!(f, "({})", s.join(","))
    }
}

/// An integer composition: an ordered tuple of positive integers.
///
/// Unlike partitions, the order of parts matters.
/// Matches `IntegerCompositions` etc. from CombinatoricTools.m.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Composition(Vec<u32>);

impl Composition {
    /// Create a composition from parts. Strips trailing zeros.
    pub fn new(parts: Vec<u32>) -> Self {
        let mut p = parts;
        while p.last() == Some(&0) {
            p.pop();
        }
        Composition(p)
    }

    /// The empty composition.
    pub fn empty() -> Self {
        Composition(vec![])
    }

    /// The parts as a slice, in the given order.
    pub fn parts(&self) -> &[u32] {
        &self.0
    }

    /// Number of parts.
    pub fn num_parts(&self) -> usize {
        self.0.len()
    }

    /// Sum of all parts.
    pub fn size(&self) -> u32 {
        self.0.iter().sum()
    }

    /// Whether this is the empty composition.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Sort into a partition.
    pub fn to_partition(&self) -> Partition {
        Partition::new(self.0.clone())
    }

    // -----------------------------------------------------------------------
    // Enumeration — IntegerCompositions from CombinatoricTools.m
    // -----------------------------------------------------------------------

    /// All compositions of n (with all positive parts).
    /// `IntegerCompositions[n]` in Mathematica.
    pub fn integer_compositions(n: u32) -> Vec<Composition> {
        if n == 0 {
            return vec![Composition::empty()];
        }
        let mut results = Vec::new();
        Self::compositions_helper(n, &mut vec![], &mut results);
        results
    }

    fn compositions_helper(remaining: u32, current: &mut Vec<u32>, results: &mut Vec<Composition>) {
        if remaining == 0 {
            results.push(Composition(current.clone()));
            return;
        }
        for part in 1..=remaining {
            current.push(part);
            Self::compositions_helper(remaining - part, current, results);
            current.pop();
        }
    }

    /// All compositions of n with exactly k parts.
    /// `IntegerCompositions[n, k]` in Mathematica.
    pub fn integer_compositions_k(n: u32, k: usize) -> Vec<Composition> {
        if k == 0 {
            return if n == 0 {
                vec![Composition::empty()]
            } else {
                vec![]
            };
        }
        let mut results = Vec::new();
        Self::compositions_k_helper(n, k, &mut vec![], &mut results);
        results
    }

    fn compositions_k_helper(
        remaining: u32,
        parts_left: usize,
        current: &mut Vec<u32>,
        results: &mut Vec<Composition>,
    ) {
        if parts_left == 0 {
            if remaining == 0 {
                results.push(Composition(current.clone()));
            }
            return;
        }
        let max_part = if parts_left > 1 {
            remaining - (parts_left as u32 - 1)
        } else {
            remaining
        };
        for part in 1..=max_part {
            current.push(part);
            Self::compositions_k_helper(remaining - part, parts_left - 1, current, results);
            current.pop();
        }
    }

    /// All weak compositions of n with exactly k parts (parts >= 0).
    /// `WeakIntegerCompositions[n, k]` in Mathematica.
    pub fn weak_integer_compositions(n: u32, k: usize) -> Vec<Composition> {
        if k == 0 {
            return if n == 0 {
                vec![Composition::empty()]
            } else {
                vec![]
            };
        }
        let mut results = Vec::new();
        Self::weak_compositions_helper(n, k, &mut vec![], &mut results);
        results
    }

    fn weak_compositions_helper(
        remaining: u32,
        parts_left: usize,
        current: &mut Vec<u32>,
        results: &mut Vec<Composition>,
    ) {
        if parts_left == 1 {
            let mut c = current.clone();
            c.push(remaining);
            results.push(Composition(c));
            return;
        }
        for part in 0..=remaining {
            current.push(part);
            Self::weak_compositions_helper(remaining - part, parts_left - 1, current, results);
            current.pop();
        }
    }

    // -----------------------------------------------------------------------
    // Descent set / ribbon conversions
    // -----------------------------------------------------------------------

    /// Convert composition to descent set: partial sums {α_1, α_1+α_2, ...} excluding the total.
    /// `CompositionToDescentSet` in Mathematica.
    pub fn composition_to_descent_set(&self) -> BTreeSet<u32> {
        let mut set = BTreeSet::new();
        let mut sum = 0u32;
        for (i, &part) in self.0.iter().enumerate() {
            sum += part;
            if i < self.0.len() - 1 {
                set.insert(sum);
            }
        }
        set
    }

    /// Reconstruct composition from descent set and total n.
    /// `DescentSetToComposition` in Mathematica.
    pub fn from_descent_set(descent_set: &BTreeSet<u32>, n: u32) -> Composition {
        if n == 0 {
            return Composition::empty();
        }
        let mut parts = Vec::new();
        let mut prev = 0u32;
        for &d in descent_set {
            parts.push(d - prev);
            prev = d;
        }
        parts.push(n - prev);
        Composition(parts)
    }

    /// Convert composition to ribbon (border strip) as a skew shape (outer, inner).
    /// `CompositionToRibbon` in Mathematica.
    pub fn composition_to_ribbon(&self) -> (Partition, Partition) {
        if self.is_empty() {
            return (Partition::empty(), Partition::empty());
        }
        let k = self.num_parts();
        let mut outer_parts = Vec::new();
        let mut inner_parts = Vec::new();

        let reversed: Vec<u32> = self.0.iter().rev().copied().collect();
        let mut offset = 0u32;

        for i in 0..k {
            let width = reversed[i];
            let right = offset + width;
            outer_parts.push(right);
            inner_parts.push(offset);
            if i < k - 1 {
                offset = right - 1;
            }
        }

        outer_parts.reverse();
        inner_parts.reverse();

        (
            Partition::from_sorted(outer_parts),
            Partition::from_sorted(inner_parts),
        )
    }

    /// All refinements of a composition (split each part into a composition).
    /// `CompositionRefinements` in Mathematica.
    pub fn composition_refinements(&self) -> Vec<Composition> {
        if self.is_empty() {
            return vec![Composition::empty()];
        }
        let part_compositions: Vec<Vec<Vec<u32>>> = self
            .0
            .iter()
            .map(|&p| {
                Composition::integer_compositions(p)
                    .into_iter()
                    .map(|c| c.0)
                    .collect()
            })
            .collect();

        let mut results = vec![vec![]];
        for choices in &part_compositions {
            let mut new_results = Vec::new();
            for existing in &results {
                for choice in choices {
                    let mut combined = existing.clone();
                    combined.extend(choice);
                    new_results.push(combined);
                }
            }
            results = new_results;
        }

        results.into_iter().map(Composition).collect()
    }

    /// Binary word encoding of composition: 1 at part boundaries, 0 elsewhere.
    /// `CompositionWord` in Mathematica.
    pub fn composition_word(&self) -> Vec<u8> {
        if self.is_empty() {
            return vec![];
        }
        let n = self.size() as usize;
        let mut word = vec![0u8; n - 1];
        let mut pos = 0usize;
        for (i, &part) in self.0.iter().enumerate() {
            pos += part as usize;
            if i < self.0.len() - 1 && pos > 0 && pos <= n - 1 {
                word[pos - 1] = 1;
            }
        }
        word
    }

    /// Reconstruct composition from binary word.
    /// `WordComposition` in Mathematica.
    pub fn from_composition_word(word: &[u8]) -> Composition {
        if word.is_empty() {
            return Composition(vec![1]);
        }
        let mut parts = Vec::new();
        let mut current = 1u32;
        for &bit in word {
            if bit == 1 {
                parts.push(current);
                current = 1;
            } else {
                current += 1;
            }
        }
        parts.push(current);
        Composition(parts)
    }

    /// Apply the slinky rule (Egge-Loehr-Warrington) to straighten s_α into ±s_λ.
    ///
    /// Returns `Some((partition, sign))` where sign is ±1, or `None` if s_α = 0.
    pub fn composition_slinky(&self) -> Option<(Partition, i64)> {
        if self.is_empty() {
            return Some((Partition::empty(), 1));
        }

        let mut parts: Vec<i64> = self.0.iter().map(|&p| p as i64).collect();
        let mut sign: i64 = 1;

        loop {
            let bad = parts.windows(2).position(|w| w[0] < w[1]);
            match bad {
                None => {
                    let p: Vec<u32> = parts.iter().map(|&x| x as u32).collect();
                    return Some((Partition::from_sorted(p), sign));
                }
                Some(i) => {
                    let a = parts[i];
                    let b = parts[i + 1];
                    if b == a + 1 {
                        return None;
                    }
                    parts[i] = b - 1;
                    parts[i + 1] = a + 1;
                    sign = -sign;

                    if parts[i] <= 0 || parts[i + 1] <= 0 {
                        return None;
                    }
                }
            }
        }
    }
}

impl fmt::Display for Composition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            write!(f, "∅")
        } else {
            let s: Vec<String> = self.0.iter().map(|p| p.to_string()).collect();
            write!(f, "({})", s.join(","))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compositions_of_4() {
        let comps = Composition::integer_compositions(4);
        assert_eq!(comps.len(), 8);
    }

    #[test]
    fn test_compositions_k() {
        let comps = Composition::integer_compositions_k(5, 3);
        assert_eq!(comps.len(), 6);
        for c in &comps {
            assert_eq!(c.size(), 5);
            assert_eq!(c.num_parts(), 3);
        }
    }

    #[test]
    fn test_weak_compositions() {
        let comps = Composition::weak_integer_compositions(3, 2);
        assert_eq!(comps.len(), 4);
    }

    #[test]
    fn test_descent_set_roundtrip() {
        let c = Composition::new(vec![2, 3, 1]);
        let des = c.composition_to_descent_set();
        let n = c.size();
        let c2 = Composition::from_descent_set(&des, n);
        assert_eq!(c, c2);
    }

    #[test]
    fn test_ribbon() {
        let c = Composition::new(vec![2, 3]);
        let (outer, inner) = c.composition_to_ribbon();
        assert_eq!(outer.size() - inner.size(), c.size());
    }

    #[test]
    fn test_refinements() {
        let c = Composition::new(vec![3]);
        let refs = c.composition_refinements();
        assert_eq!(refs.len(), 4);
    }

    #[test]
    fn test_composition_word_roundtrip() {
        let c = Composition::new(vec![2, 1, 3]);
        let word = c.composition_word();
        let c2 = Composition::from_composition_word(&word);
        assert_eq!(c, c2);
    }

    #[test]
    fn test_slinky_partition() {
        let c = Composition::new(vec![3, 2, 1]);
        let (lam, sign) = c.composition_slinky().unwrap();
        assert_eq!(lam, Partition::new(vec![3, 2, 1]));
        assert_eq!(sign, 1);
    }

    #[test]
    fn test_slinky_zero_12() {
        let c = Composition::new(vec![1, 2]);
        assert!(c.composition_slinky().is_none());
    }

    // WeakComposition tests

    #[test]
    fn test_weak_composition_preserves_trailing_zeros() {
        let wc = WeakComposition::new(vec![2, 0, 1, 0]);
        assert_eq!(wc.num_vars(), 4);
        assert_eq!(wc.parts(), &[2, 0, 1, 0]);
    }

    #[test]
    fn test_weak_composition_degree() {
        let wc = WeakComposition::new(vec![2, 0, 3]);
        assert_eq!(wc.degree(), 5);
    }

    #[test]
    fn test_weak_composition_is_dominant() {
        assert!(WeakComposition::new(vec![3, 2, 0]).is_dominant());
        assert!(!WeakComposition::new(vec![1, 2, 0]).is_dominant());
    }

    #[test]
    fn test_weak_composition_flat() {
        let wc = WeakComposition::new(vec![2, 0, 3, 0, 1]);
        assert_eq!(wc.flat(), Composition::new(vec![2, 3, 1]));
    }

    #[test]
    fn test_weak_composition_content() {
        let wc = WeakComposition::new(vec![1, 0, 3, 2]);
        assert_eq!(wc.content(), Partition::new(vec![3, 2, 1]));
    }

    #[test]
    fn test_weak_composition_enumeration_count() {
        // C(d+k-1, k-1) = C(3+3-1, 3-1) = C(5,2) = 10
        assert_eq!(WeakComposition::all_weak_compositions(3, 3).len(), 10);
        // C(0+3-1, 3-1) = C(2,2) = 1
        assert_eq!(WeakComposition::all_weak_compositions(0, 3).len(), 1);
        // C(2+1-1, 1-1) = C(2,0) = 1
        assert_eq!(WeakComposition::all_weak_compositions(2, 1).len(), 1);
        // k=0, d=0 -> 1
        assert_eq!(WeakComposition::all_weak_compositions(0, 0).len(), 1);
        // k=0, d>0 -> 0
        assert_eq!(WeakComposition::all_weak_compositions(3, 0).len(), 0);
    }

    #[test]
    fn test_weak_composition_enumeration_length() {
        for wc in WeakComposition::all_weak_compositions(4, 3) {
            assert_eq!(wc.num_vars(), 3);
            assert_eq!(wc.degree(), 4);
        }
    }

    #[test]
    fn test_weak_composition_from_composition() {
        let c = Composition::new(vec![2, 3, 1]);
        let wc = WeakComposition::from_composition(&c, 5);
        assert_eq!(wc.parts(), &[2, 3, 1, 0, 0]);
    }

    #[test]
    fn test_weak_composition_ordering() {
        let a = WeakComposition::new(vec![0, 3]);
        let b = WeakComposition::new(vec![1, 2]);
        assert!(a < b); // lexicographic
    }
}
