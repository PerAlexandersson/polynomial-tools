use num_bigint::BigUint;
use num_traits::One;
use std::fmt;

/// A partition of a non-negative integer: a weakly decreasing sequence of positive integers.
///
/// Stored as `Vec<u32>` in weakly decreasing order with no trailing zeros.
/// Following Mathematica conventions from CombinatoricTools.m.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Partition(Vec<u32>);

impl Partition {
    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    /// Create a partition from parts. Sorts descending and strips zeros.
    pub fn new(mut parts: Vec<u32>) -> Self {
        parts.sort_unstable_by(|a, b| b.cmp(a));
        parts.retain(|&x| x > 0);
        Partition(parts)
    }

    /// Create a partition from parts that are already sorted descending.
    /// Strips trailing zeros but does NOT re-sort. Caller must guarantee order.
    pub fn from_sorted(mut parts: Vec<u32>) -> Self {
        parts.retain(|&x| x > 0);
        Partition(parts)
    }

    /// The empty partition (partition of 0).
    pub fn empty() -> Self {
        Partition(vec![])
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /// The parts as a slice.
    pub fn parts(&self) -> &[u32] {
        &self.0
    }

    /// Number of (nonzero) parts, i.e. the length ℓ(λ).
    pub fn num_parts(&self) -> usize {
        self.0.len()
    }

    /// The i-th part (0-indexed). Returns 0 if i >= num_parts.
    pub fn part(&self, i: usize) -> u32 {
        if i < self.0.len() {
            self.0[i]
        } else {
            0
        }
    }

    /// The size |λ| = sum of all parts.
    pub fn size(&self) -> u32 {
        checked_part_sum(&self.0)
    }

    /// Whether this is the empty partition.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    // -----------------------------------------------------------------------
    // Partition combinatorics — matching CombinatoricTools.m naming
    // -----------------------------------------------------------------------

    /// Conjugate (transpose) partition.
    /// `ConjugatePartition` in Mathematica.
    pub fn conjugate_partition(&self) -> Partition {
        if self.is_empty() {
            return Partition::empty();
        }
        let max_part = self.0[0] as usize;
        let mut conj = vec![0u32; max_part];
        for &p in &self.0 {
            for j in 0..(p as usize) {
                conj[j] += 1;
            }
        }
        Partition(conj)
    }

    /// Test Young diagram containment: self ⊆ other (entrywise).
    /// `PartitionLessEqualQ` in Mathematica.
    pub fn partition_less_equal(&self, other: &Partition) -> bool {
        let n = self.num_parts().max(other.num_parts());
        for i in 0..n {
            if self.part(i) > other.part(i) {
                return false;
            }
        }
        true
    }

    /// Dominance order: self ≤_dom other iff partial sums of self ≤ partial sums of other.
    /// `PartitionDominatesQ[other, self]` in Mathematica (note: Mathematica asks "does p1 dominate p2").
    /// Here we ask: does `other` dominate `self`?
    pub fn dominated_by(&self, other: &Partition) -> bool {
        if self.size() != other.size() {
            return false;
        }
        let n = self.num_parts().max(other.num_parts());
        let mut sum_self: u32 = 0;
        let mut sum_other: u32 = 0;
        for i in 0..n {
            sum_self = sum_self
                .checked_add(self.part(i))
                .expect("partition partial sum overflow");
            sum_other = sum_other
                .checked_add(other.part(i))
                .expect("partition partial sum overflow");
            if sum_self > sum_other {
                return false;
            }
        }
        true
    }

    /// Strict dominance: `other` strictly dominates `self`.
    pub fn strictly_dominated_by(&self, other: &Partition) -> bool {
        self.dominated_by(other) && self != other
    }

    /// All partitions obtained by adding one box.
    /// `PartitionAddBox` in Mathematica.
    pub fn partition_add_box(&self) -> Vec<Partition> {
        let mut results = Vec::new();
        let n = self.num_parts();
        // Can add box to each row where it doesn't violate decreasing order
        for i in 0..=n {
            let prev = if i == 0 { u32::MAX } else { self.part(i - 1) };
            let curr = self.part(i);
            if curr < prev {
                let mut new_parts = self.0.clone();
                if i < n {
                    new_parts[i] = new_parts[i]
                        .checked_add(1)
                        .expect("partition part overflow");
                } else {
                    new_parts.push(1);
                }
                results.push(Partition(new_parts));
            }
        }
        results
    }

    /// All partitions obtained by removing one box.
    /// `PartitionRemoveBox` in Mathematica.
    pub fn partition_remove_box(&self) -> Vec<Partition> {
        let mut results = Vec::new();
        let n = self.num_parts();
        for i in 0..n {
            let curr = self.0[i];
            let next = self.part(i + 1);
            if curr > next {
                let mut new_parts = self.0.clone();
                new_parts[i] -= 1;
                // Remove trailing zeros
                new_parts.retain(|&x| x > 0);
                results.push(Partition(new_parts));
            }
        }
        results
    }

    /// Arm length at box (r, c) (0-indexed): λ_r - c - 1.
    /// `PartitionArm` in Mathematica (but Mathematica uses 1-indexed).
    /// Returns None if the box is not in the diagram.
    pub fn partition_arm(&self, row: usize, col: usize) -> Option<u32> {
        let part_r = self.part(row);
        if col < part_r as usize {
            Some(part_r - col as u32 - 1)
        } else {
            None
        }
    }

    /// Leg length at box (r, c) (0-indexed): λ'_c - r - 1.
    /// `PartitionLeg` in Mathematica.
    pub fn partition_leg(&self, row: usize, col: usize) -> Option<u32> {
        let conj = self.conjugate_partition();
        let part_c = conj.part(col);
        if row < part_c as usize {
            Some(part_c - row as u32 - 1)
        } else {
            None
        }
    }

    /// Hook length at box (r, c): arm + leg + 1.
    pub fn hook_length(&self, row: usize, col: usize) -> Option<u32> {
        match (self.partition_arm(row, col), self.partition_leg(row, col)) {
            (Some(a), Some(l)) => Some(
                a.checked_add(l)
                    .and_then(|sum| sum.checked_add(1))
                    .expect("hook length overflow"),
            ),
            _ => None,
        }
    }

    /// Table of all hook lengths.
    /// `HookLengths` in Mathematica.
    pub fn hook_lengths(&self) -> Vec<Vec<u32>> {
        let conj = self.conjugate_partition();
        self.0
            .iter()
            .enumerate()
            .map(|(r, &part_r)| {
                (0..part_r as usize)
                    .map(|c| {
                        let arm = part_r - c as u32 - 1;
                        let leg = conj.part(c) - r as u32 - 1;
                        arm.checked_add(leg)
                            .and_then(|sum| sum.checked_add(1))
                            .expect("hook length overflow")
                    })
                    .collect()
            })
            .collect()
    }

    /// Count standard Young tableaux of shape λ via the hook-length formula:
    /// f^λ = |λ|! / ∏ hook(i,j).
    /// `NumberOfStandardYoungTableaux` in Mathematica.
    pub fn count_syt(&self) -> BigUint {
        let n = self.size() as u64;
        let mut numerator = BigUint::one();
        for i in 2..=n {
            numerator *= i;
        }
        for row in &self.hook_lengths() {
            for &h in row {
                numerator /= h as u64;
            }
        }
        numerator
    }

    /// All (row, col) coordinates in the diagram (0-indexed).
    /// `DiagramBoxes` in Mathematica.
    pub fn diagram_boxes(&self) -> Vec<(usize, usize)> {
        let mut boxes = Vec::new();
        for (r, &part_r) in self.0.iter().enumerate() {
            for c in 0..part_r as usize {
                boxes.push((r, c));
            }
        }
        boxes
    }

    /// All (row, col) in skew shape self / inner (0-indexed).
    pub fn skew_diagram_boxes(&self, inner: &Partition) -> Vec<(usize, usize)> {
        let mut boxes = Vec::new();
        for (r, &part_r) in self.0.iter().enumerate() {
            let inner_r = inner.part(r);
            for c in inner_r as usize..part_r as usize {
                boxes.push((r, c));
            }
        }
        boxes
    }

    /// n(λ) = Σ (i) * λ_i (0-indexed i), equivalently Σ binom(λ'_j, 2).
    /// `PartitionN` in Mathematica.
    pub fn partition_n(&self) -> u32 {
        self.0
            .iter()
            .enumerate()
            .try_fold(0u32, |total, (i, &p)| {
                let index = u32::try_from(i).expect("partition index overflow");
                let term = index.checked_mul(p).expect("partition n overflow");
                total.checked_add(term)
            })
            .expect("partition n overflow")
    }

    /// Part multiplicity vector: m_i = number of parts equal to i.
    /// `PartitionPartCount` in Mathematica.
    pub fn partition_part_count(&self) -> Vec<u32> {
        if self.is_empty() {
            return vec![];
        }
        let max_part = self.0[0] as usize;
        let mut counts = vec![0u32; max_part + 1];
        for &p in &self.0 {
            counts[p as usize] += 1;
        }
        counts
    }

    /// z_λ = Π i^{m_i} * m_i! where m_i is the multiplicity of i in λ.
    /// `ZCoefficient` in Mathematica. Used for power-sum inner product normalization.
    pub fn z_coefficient(&self) -> u64 {
        let counts = self.partition_part_count();
        let mut z: u64 = 1;
        for (i, &m) in counts.iter().enumerate() {
            if i == 0 || m == 0 {
                continue;
            }
            // i^m * m!
            for _ in 0..m {
                z = z.checked_mul(i as u64).expect("z coefficient overflow");
            }
            for k in 1..=m as u64 {
                z = z.checked_mul(k).expect("z coefficient overflow");
            }
        }
        z
    }

    /// Durfee square size: largest d such that λ_d ≥ d (1-indexed, so we check part(d-1) >= d).
    /// `Durfee` in Mathematica.
    pub fn durfee_square(&self) -> usize {
        let mut d = 0;
        while d < self.num_parts() && self.0[d] as usize > d {
            d += 1;
        }
        d
    }

    // -----------------------------------------------------------------------
    // Cores and quotients (abacus)
    // -----------------------------------------------------------------------

    /// Convert partition to beta-set (first-column hook lengths):
    /// β_i = λ_i + (ℓ(λ) - 1 - i) for i = 0..ℓ(λ)-1.
    fn to_beta_set(&self) -> Vec<u32> {
        let n = self.num_parts();
        (0..n)
            .map(|i| {
                let shift = u32::try_from(n - 1 - i).expect("partition index overflow");
                self.0[i].checked_add(shift).expect("beta-set overflow")
            })
            .collect()
    }

    /// Reconstruct partition from beta-set (sorted descending).
    fn from_beta_set(beta: &[u32]) -> Partition {
        let n = beta.len();
        if n == 0 {
            return Partition::empty();
        }
        let parts: Vec<u32> = (0..n)
            .map(|i| {
                if beta[i] >= (n - 1 - i) as u32 {
                    beta[i] - (n - 1 - i) as u32
                } else {
                    0
                }
            })
            .collect();
        Partition::from_sorted(parts)
    }

    /// d-core of partition (remove all d-hooks).
    /// `PartitionCore` in Mathematica.
    pub fn partition_core(&self, d: u32) -> Partition {
        assert!(d > 0, "core modulus must be positive");
        // Use beta-set approach: reduce each beta number mod d, then reconstruct
        let beta = self.to_beta_set();
        // Sort residues: for each residue class, collect and reassign minimally
        let mut residue_classes: Vec<Vec<u32>> = vec![vec![]; d as usize];
        for &b in &beta {
            residue_classes[(b % d) as usize].push(b);
        }
        // For each class, replace values with the minimal set: r, r+d, r+2d, ...
        let mut new_beta = Vec::new();
        for (r, class) in residue_classes.iter().enumerate() {
            for k in 0..class.len() {
                let residue = u32::try_from(r).expect("partition residue overflow");
                let multiple = u32::try_from(k)
                    .expect("partition beta-set overflow")
                    .checked_mul(d)
                    .expect("partition beta-set overflow");
                new_beta.push(
                    residue
                        .checked_add(multiple)
                        .expect("partition beta-set overflow"),
                );
            }
        }
        new_beta.sort_unstable_by(|a, b| b.cmp(a));
        Partition::from_beta_set(&new_beta)
    }

    /// d-quotient of partition: a tuple of d partitions.
    /// `PartitionQuotient` in Mathematica.
    pub fn partition_quotient(&self, d: u32) -> Vec<Partition> {
        assert!(d > 0, "quotient modulus must be positive");
        let beta = self.to_beta_set();
        let mut residue_classes: Vec<Vec<u32>> = vec![vec![]; d as usize];
        for &b in &beta {
            residue_classes[(b % d) as usize].push(b);
        }
        residue_classes
            .iter()
            .map(|class| {
                // Divide by d to get the quotient partition's beta-set
                let mut quot_beta: Vec<u32> = class.iter().map(|&b| b / d).collect();
                quot_beta.sort_unstable_by(|a, b| b.cmp(a));
                Partition::from_beta_set(&quot_beta)
            })
            .collect()
    }

    // -----------------------------------------------------------------------
    // Horizontal and vertical strips
    // -----------------------------------------------------------------------

    /// All partitions μ such that self/μ is a horizontal strip of size k.
    /// (At most one box removed per column.)
    /// `PartitionRemoveHorizontalStrip` in Mathematica.
    pub fn partition_remove_horizontal_strip(&self, k: u32) -> Vec<Partition> {
        if k == 0 {
            return vec![self.clone()];
        }
        if k > self.size() {
            return vec![];
        }
        let outer = self;
        let conj = outer.conjugate_partition();
        let mut results = Vec::new();
        Self::horizontal_strip_helper(outer, &conj, k, 0, vec![], &mut results);
        results
    }

    fn horizontal_strip_helper(
        outer: &Partition,
        conj_outer: &Partition,
        remaining: u32,
        row: usize,
        current: Vec<u32>,
        results: &mut Vec<Partition>,
    ) {
        let _ = conj_outer; // kept for API compatibility
        if remaining == 0 {
            // Fill in remaining rows unchanged
            let mut parts = current;
            for i in row..outer.num_parts() {
                parts.push(outer.0[i]);
            }
            results.push(Partition::from_sorted(parts));
            return;
        }
        if row > outer.num_parts() {
            return;
        }

        let outer_r = outer.part(row);
        let lower = outer.part(row + 1); // mu[r] >= outer[r+1] for horizontal strip
        let upper = outer_r;
        let max_remove = (upper - lower).min(remaining);

        for remove in 0..=max_remove {
            let new_r = upper - remove;
            if !current.is_empty() && new_r > *current.last().unwrap() {
                continue;
            }
            let mut next = current.clone();
            next.push(new_r);
            Self::horizontal_strip_helper(
                outer,
                conj_outer,
                remaining - remove,
                row + 1,
                next,
                results,
            );
        }
    }

    /// All partitions μ such that self/μ is a vertical strip of size k.
    /// `PartitionRemoveVerticalStrip` in Mathematica.
    pub fn partition_remove_vertical_strip(&self, k: u32) -> Vec<Partition> {
        // Vertical strip of λ/μ iff horizontal strip of λ'/μ'
        let conj = self.conjugate_partition();
        conj.partition_remove_horizontal_strip(k)
            .into_iter()
            .map(|mu| mu.conjugate_partition())
            .collect()
    }

    // -----------------------------------------------------------------------
    // Partition join
    // -----------------------------------------------------------------------

    /// Union of parts (sort combined parts into a partition).
    /// `PartitionJoin` in Mathematica. Used for multiplicative bases (e, h, p).
    pub fn partition_join(&self, other: &Partition) -> Partition {
        let mut parts: Vec<u32> = self.0.iter().chain(other.0.iter()).copied().collect();
        parts.sort_unstable_by(|a, b| b.cmp(a));
        Partition(parts)
    }

    // -----------------------------------------------------------------------
    // Disjoint union of skew shapes
    // -----------------------------------------------------------------------

    /// Construct a skew shape ρ/σ that is the disjoint union of the given
    /// skew shapes (no shared rows or columns between any two shapes).
    pub fn disjoint_union_skew_shapes(shapes: &[(Partition, Partition)]) -> (Partition, Partition) {
        if shapes.is_empty() {
            return (Partition::empty(), Partition::empty());
        }
        if shapes.len() == 1 {
            return shapes[0].clone();
        }

        let widths: Vec<u32> = shapes.iter().map(|(outer, _)| outer.part(0)).collect();
        let heights: Vec<usize> = shapes.iter().map(|(outer, _)| outer.num_parts()).collect();

        let mut col_offsets = vec![0u32; shapes.len()];
        for k in (0..shapes.len() - 1).rev() {
            col_offsets[k] = col_offsets[k + 1]
                .checked_add(widths[k + 1])
                .expect("skew-shape column offset overflow");
        }

        let mut row_offsets = vec![0usize; shapes.len()];
        for k in 1..shapes.len() {
            row_offsets[k] = row_offsets[k - 1]
                .checked_add(heights[k - 1])
                .expect("skew-shape row offset overflow");
        }

        let total_rows = heights
            .iter()
            .try_fold(0usize, |total, &height| total.checked_add(height))
            .expect("skew-shape row offset overflow");
        let mut rho_parts = vec![0u32; total_rows];
        let mut sigma_parts = vec![0u32; total_rows];

        for (k, (outer, inner)) in shapes.iter().enumerate() {
            for i in 0..heights[k] {
                rho_parts[row_offsets[k] + i] = col_offsets[k]
                    .checked_add(outer.part(i))
                    .expect("skew-shape part overflow");
                sigma_parts[row_offsets[k] + i] = col_offsets[k]
                    .checked_add(inner.part(i))
                    .expect("skew-shape part overflow");
            }
        }

        (
            Partition::from_sorted(rho_parts),
            Partition::from_sorted(sigma_parts),
        )
    }

    /// Reduce a Kostka coefficient to a Littlewood-Richardson coefficient.
    pub fn kostka_to_lr(
        lambda: &Partition,
        nu: &Partition,
        mu: &Partition,
    ) -> (Partition, Partition, Partition) {
        let mut shapes: Vec<(Partition, Partition)> = Vec::new();

        if !nu.is_empty() {
            shapes.push((nu.clone(), Partition::empty()));
        }

        for &part in mu.parts() {
            shapes.push((Partition::new(vec![part]), Partition::empty()));
        }

        let (rho, sigma) = Self::disjoint_union_skew_shapes(&shapes);
        (rho, lambda.clone(), sigma)
    }

    // -----------------------------------------------------------------------
    // Enumeration
    // -----------------------------------------------------------------------

    /// All partitions of n (in reverse lexicographic order).
    pub fn all_of_size(n: u32) -> Vec<Partition> {
        Self::all_of_size_bounded(n, n as usize, n)
    }

    /// All partitions of n with at most `max_parts` parts, each at most `max_part`.
    pub fn all_of_size_bounded(n: u32, max_parts: usize, max_part: u32) -> Vec<Partition> {
        let mut result = Vec::new();
        Self::enumerate_helper(n, max_parts, max_part, &mut vec![], &mut result);
        result
    }

    fn enumerate_helper(
        remaining: u32,
        max_parts: usize,
        max_part: u32,
        current: &mut Vec<u32>,
        results: &mut Vec<Partition>,
    ) {
        if remaining == 0 {
            results.push(Partition(current.clone()));
            return;
        }
        if max_parts == 0 {
            return;
        }
        let upper = remaining.min(max_part);
        let lower = if max_parts > 0 {
            1
        } else {
            return;
        };
        for part in (lower..=upper).rev() {
            current.push(part);
            Self::enumerate_helper(
                remaining - part,
                max_parts - 1,
                part, // next part <= current
                current,
                results,
            );
            current.pop();
        }
    }

    /// All partitions of n in dominance order (largest first).
    pub fn all_of_size_dominance_order(n: u32) -> Vec<Partition> {
        let mut parts = Self::all_of_size(n);
        parts.sort_by(|a, b| {
            let len = a.num_parts().max(b.num_parts());
            let mut sa = 0u32;
            let mut sb = 0u32;
            for i in 0..len {
                sa = sa
                    .checked_add(a.part(i))
                    .expect("partition partial sum overflow");
                sb = sb
                    .checked_add(b.part(i))
                    .expect("partition partial sum overflow");
                if sa != sb {
                    return sb.cmp(&sa);
                }
            }
            std::cmp::Ordering::Equal
        });
        parts
    }

    /// All sub-partitions μ ⊆ λ (component-wise: μ_i ≤ λ_i for all i).
    ///
    /// Returns partitions in lexicographic order. Includes both ∅ and λ itself.
    pub fn sub_partitions(&self) -> Vec<Partition> {
        let n = self.num_parts();
        if n == 0 {
            return vec![Partition::empty()];
        }
        let mut results = Vec::new();
        let mut mu = vec![0u32; n];
        Self::sub_partitions_helper(self.parts(), &mut mu, 0, u32::MAX, &mut results);
        results
    }

    fn sub_partitions_helper(
        lambda: &[u32],
        mu: &mut Vec<u32>,
        pos: usize,
        max_val: u32,
        results: &mut Vec<Partition>,
    ) {
        if pos == lambda.len() {
            results.push(Partition::new(mu.clone()));
            return;
        }
        let upper = lambda[pos].min(max_val);
        for v in 0..=upper {
            mu[pos] = v;
            Self::sub_partitions_helper(lambda, mu, pos + 1, v, results);
        }
    }

    // -----------------------------------------------------------------------
    // Parsing
    // -----------------------------------------------------------------------

    /// Parse "5,3,1" or "5.3.1" into a Partition.
    pub fn parse(s: &str) -> Result<Partition, String> {
        if s.is_empty() || s == "0" {
            return Ok(Partition::empty());
        }
        let sep = if s.contains(',') { ',' } else { '.' };
        let parts: Result<Vec<u32>, _> = s.split(sep).map(|x| x.trim().parse::<u32>()).collect();
        match parts {
            Ok(p) => Ok(Partition::new(p)),
            Err(e) => Err(format!("Failed to parse partition '{}': {}", s, e)),
        }
    }

    /// Display as comma-separated parts, or "∅" for empty.
    pub fn display(&self) -> String {
        if self.is_empty() {
            "∅".to_string()
        } else {
            self.0
                .iter()
                .map(|p| p.to_string())
                .collect::<Vec<_>>()
                .join(",")
        }
    }
}

fn checked_part_sum(parts: &[u32]) -> u32 {
    parts
        .iter()
        .try_fold(0u32, |total, &part| total.checked_add(part))
        .expect("partition size overflow")
}

impl fmt::Display for Partition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.display())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic() {
        let p = Partition::new(vec![3, 5, 1]);
        assert_eq!(p.parts(), &[5, 3, 1]);
        assert_eq!(p.size(), 9);
        assert_eq!(p.num_parts(), 3);
    }

    #[test]
    #[should_panic(expected = "partition size overflow")]
    fn test_size_rejects_overflow() {
        let p = Partition::from_sorted(vec![u32::MAX, 1]);

        let _ = p.size();
    }

    #[test]
    fn test_conjugate() {
        let p = Partition::new(vec![4, 2, 1]);
        assert_eq!(p.conjugate_partition(), Partition::new(vec![3, 2, 1, 1]));
        assert_eq!(p.conjugate_partition().conjugate_partition(), p);
    }

    #[test]
    fn test_dominance() {
        let p = Partition::new(vec![3, 1, 1]);
        let q = Partition::new(vec![2, 2, 1]);
        assert!(q.dominated_by(&p));
        assert!(!p.dominated_by(&q));
    }

    #[test]
    fn test_add_remove_box() {
        let p = Partition::new(vec![2, 1]);
        let added = p.partition_add_box();
        assert_eq!(added.len(), 3);

        let p2 = Partition::new(vec![3, 1]);
        let removed = p2.partition_remove_box();
        assert_eq!(removed.len(), 2);
    }

    #[test]
    fn test_hook_lengths() {
        let p = Partition::new(vec![3, 2]);
        let hooks = p.hook_lengths();
        assert_eq!(hooks, vec![vec![4, 3, 1], vec![2, 1]]);
    }

    #[test]
    fn test_count_syt() {
        let p = Partition::new(vec![3, 2, 1]);
        assert_eq!(p.count_syt(), BigUint::from(16u32));
        let q = Partition::new(vec![4, 3, 2, 1]);
        assert_eq!(q.count_syt(), BigUint::from(768u32));
        assert_eq!(Partition::new(vec![1]).count_syt(), BigUint::one());
        assert_eq!(Partition::empty().count_syt(), BigUint::one());
    }

    #[test]
    fn test_partition_n() {
        let p = Partition::new(vec![3, 2, 1]);
        assert_eq!(p.partition_n(), 4);
    }

    #[test]
    #[should_panic(expected = "partition n overflow")]
    fn test_partition_n_rejects_overflow() {
        let p = Partition::from_sorted(vec![u32::MAX, u32::MAX, u32::MAX]);

        let _ = p.partition_n();
    }

    #[test]
    fn test_z_coefficient() {
        let p = Partition::new(vec![2, 2, 1]);
        assert_eq!(p.z_coefficient(), 8);
    }

    #[test]
    #[should_panic(expected = "z coefficient overflow")]
    fn test_z_coefficient_rejects_overflow() {
        let p = Partition::new(vec![1; 65]);

        let _ = p.z_coefficient();
    }

    #[test]
    fn test_durfee() {
        let p = Partition::new(vec![4, 3, 2, 1]);
        assert_eq!(p.durfee_square(), 2);

        let q = Partition::new(vec![3, 3, 3]);
        assert_eq!(q.durfee_square(), 3);
    }

    #[test]
    fn test_enumeration() {
        let parts = Partition::all_of_size(5);
        assert_eq!(parts.len(), 7);
    }

    #[test]
    fn test_core_quotient() {
        let p = Partition::new(vec![3, 2, 1]);
        let core = p.partition_core(2);
        let hooks = core.hook_lengths();
        for row in &hooks {
            for &h in row {
                assert!(h % 2 != 0, "2-core should have no even hooks, got {}", h);
            }
        }
    }

    #[test]
    #[should_panic(expected = "core modulus must be positive")]
    fn test_partition_core_rejects_zero_modulus() {
        let p = Partition::new(vec![3, 2, 1]);

        let _ = p.partition_core(0);
    }

    #[test]
    #[should_panic(expected = "quotient modulus must be positive")]
    fn test_partition_quotient_rejects_zero_modulus() {
        let p = Partition::new(vec![3, 2, 1]);

        let _ = p.partition_quotient(0);
    }

    #[test]
    #[should_panic(expected = "beta-set overflow")]
    fn test_beta_set_rejects_overflow() {
        let p = Partition::from_sorted(vec![u32::MAX, 1]);

        let _ = p.partition_core(2);
    }

    #[test]
    fn test_horizontal_strip() {
        let p = Partition::new(vec![3, 2]);
        let strips = p.partition_remove_horizontal_strip(2);
        for mu in &strips {
            assert_eq!(p.size() - mu.size(), 2);
            assert!(mu.partition_less_equal(&p));
        }
    }

    #[test]
    fn test_partition_join() {
        let a = Partition::new(vec![3, 1]);
        let b = Partition::new(vec![2, 2]);
        assert_eq!(a.partition_join(&b), Partition::new(vec![3, 2, 2, 1]));
    }

    #[test]
    fn test_parse() {
        assert_eq!(
            Partition::parse("5,3,1").unwrap(),
            Partition::new(vec![5, 3, 1])
        );
        assert_eq!(Partition::parse("0").unwrap(), Partition::empty());
    }

    #[test]
    fn test_disjoint_union_two_rows() {
        let shapes = vec![
            (Partition::new(vec![3]), Partition::empty()),
            (Partition::new(vec![2]), Partition::empty()),
        ];
        let (rho, sigma) = Partition::disjoint_union_skew_shapes(&shapes);
        assert_eq!(rho, Partition::new(vec![5, 2]));
        assert_eq!(sigma, Partition::new(vec![2]));
    }

    #[test]
    #[should_panic(expected = "skew-shape part overflow")]
    fn test_disjoint_union_rejects_part_overflow() {
        let wide = (Partition::from_sorted(vec![u32::MAX]), Partition::empty());
        let shift = (Partition::from_sorted(vec![1]), Partition::empty());

        let _ = Partition::disjoint_union_skew_shapes(&[wide, shift]);
    }

    #[test]
    fn test_sub_partitions_empty() {
        assert_eq!(
            Partition::empty().sub_partitions(),
            vec![Partition::empty()]
        );
    }

    #[test]
    fn test_sub_partitions_single_row() {
        // Sub-partitions of (3): (0), (1), (2), (3) = 4
        let subs = Partition::new(vec![3]).sub_partitions();
        assert_eq!(subs.len(), 4);
        assert!(subs.contains(&Partition::empty()));
        assert!(subs.contains(&Partition::new(vec![3])));
    }

    #[test]
    fn test_sub_partitions_21() {
        // Sub-partitions of (2,1): mu_1 <= 2, mu_2 <= 1, mu_1 >= mu_2
        // (0,0), (1,0), (1,1), (2,0), (2,1) = 5
        let subs = Partition::new(vec![2, 1]).sub_partitions();
        assert_eq!(subs.len(), 5);
        for mu in &subs {
            assert!(mu.partition_less_equal(&Partition::new(vec![2, 1])));
        }
    }

    #[test]
    fn test_sub_partitions_contains_extremes() {
        let lam = Partition::new(vec![3, 2, 1]);
        let subs = lam.sub_partitions();
        assert!(subs.contains(&Partition::empty()));
        assert!(subs.contains(&lam));
    }
}
