//! Set partitions and ordered set partitions of `{0, 1, ..., n-1}`.
//!
//! The canonical encoding for an ordinary set partition uses standard
//! restricted growth strings (RGS). Ordered set partitions are obtained by
//! taking all permutations of the blocks of an ordinary set partition.

use crate::permutation::next_permutation;

/// An ordinary set partition of `{0, 1, ..., n-1}`.
///
/// Blocks are stored with entries in increasing order, and the blocks
/// themselves are stored in canonical order by increasing least element.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SetPartition {
    blocks: Vec<Vec<usize>>,
}

impl SetPartition {
    /// Construct a set partition from blocks.
    ///
    /// The input blocks may be unordered; each block is sorted increasingly and
    /// the blocks are reordered canonically by their least element. Panics if
    /// the blocks do not form a partition of `{0, 1, ..., n-1}` for some `n`.
    pub fn new(mut blocks: Vec<Vec<usize>>) -> Self {
        normalize_blocks(&mut blocks);
        blocks.sort_unstable_by_key(|block| block[0]);
        validate_ground_set(&blocks);
        Self { blocks }
    }

    /// Construct a set partition from a standard restricted growth string.
    ///
    /// Panics if `rgs` is not a valid standard restricted growth string.
    pub fn from_restricted_growth_string(rgs: &[usize]) -> Self {
        assert!(
            is_standard_restricted_growth_string(rgs),
            "invalid standard restricted growth string: {rgs:?}"
        );
        if rgs.is_empty() {
            return Self { blocks: Vec::new() };
        }

        let num_blocks = rgs.iter().copied().max().unwrap_or(0) + 1;
        let mut blocks = vec![Vec::new(); num_blocks];
        for (i, &block) in rgs.iter().enumerate() {
            blocks[block].push(i);
        }
        Self { blocks }
    }

    /// Return the blocks of the partition.
    pub fn blocks(&self) -> &[Vec<usize>] {
        &self.blocks
    }

    /// Consume the partition and return its blocks.
    pub fn into_blocks(self) -> Vec<Vec<usize>> {
        self.blocks
    }

    /// Number of blocks.
    pub fn num_blocks(&self) -> usize {
        self.blocks.len()
    }

    /// Total size of the ground set.
    pub fn ground_set_size(&self) -> usize {
        self.blocks.iter().map(Vec::len).sum()
    }

    /// Whether the partition is empty.
    pub fn is_empty(&self) -> bool {
        self.blocks.is_empty()
    }

    /// The block-size sequence.
    pub fn block_sizes(&self) -> Vec<usize> {
        self.blocks.iter().map(Vec::len).collect()
    }

    /// The standard restricted growth string of the partition.
    pub fn to_restricted_growth_string(&self) -> Vec<usize> {
        let n = self.ground_set_size();
        let mut rgs = vec![0usize; n];
        for (block_idx, block) in self.blocks.iter().enumerate() {
            for &x in block {
                rgs[x] = block_idx;
            }
        }
        rgs
    }
}

/// An ordered set partition of `{0, 1, ..., n-1}`.
///
/// Blocks are stored in their given order, with each individual block sorted
/// increasingly.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct OrderedSetPartition {
    blocks: Vec<Vec<usize>>,
}

impl OrderedSetPartition {
    /// Construct an ordered set partition from ordered blocks.
    ///
    /// Each block is sorted increasingly, but the block order is preserved.
    /// Panics if the blocks do not form a partition of `{0, 1, ..., n-1}` for
    /// some `n`.
    pub fn new(mut blocks: Vec<Vec<usize>>) -> Self {
        normalize_blocks(&mut blocks);
        validate_ground_set(&blocks);
        Self { blocks }
    }

    /// Return the blocks of the ordered partition.
    pub fn blocks(&self) -> &[Vec<usize>] {
        &self.blocks
    }

    /// Consume the ordered partition and return its blocks.
    pub fn into_blocks(self) -> Vec<Vec<usize>> {
        self.blocks
    }

    /// Number of blocks.
    pub fn num_blocks(&self) -> usize {
        self.blocks.len()
    }

    /// Total size of the ground set.
    pub fn ground_set_size(&self) -> usize {
        self.blocks.iter().map(Vec::len).sum()
    }

    /// Whether the partition is empty.
    pub fn is_empty(&self) -> bool {
        self.blocks.is_empty()
    }

    /// The block-size sequence.
    pub fn block_sizes(&self) -> Vec<usize> {
        self.blocks.iter().map(Vec::len).collect()
    }

    /// Forget the block order and return the underlying ordinary set partition.
    pub fn as_set_partition(&self) -> SetPartition {
        SetPartition::new(self.blocks.clone())
    }
}

/// Iterator over ordinary set partitions of `{0, 1, ..., n-1}`.
#[derive(Debug, Clone)]
pub struct SetPartitions {
    rgs: Vec<usize>,
    max_class: Vec<usize>,
    first: bool,
    done: bool,
}

impl SetPartitions {
    /// Create the iterator over all set partitions of `{0, 1, ..., n-1}`.
    pub fn new(n: usize) -> Self {
        Self {
            rgs: vec![0; n],
            max_class: vec![0; n],
            first: true,
            done: false,
        }
    }
}

impl Iterator for SetPartitions {
    type Item = SetPartition;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        if self.first {
            self.first = false;
            return Some(SetPartition::from_restricted_growth_string(&self.rgs));
        }
        if !next_standard_restricted_growth_string(&mut self.rgs, &mut self.max_class) {
            self.done = true;
            return None;
        }
        Some(SetPartition::from_restricted_growth_string(&self.rgs))
    }
}

/// Iterator over ordered set partitions of `{0, 1, ..., n-1}`.
#[derive(Debug, Clone)]
pub struct OrderedSetPartitions {
    set_partitions: SetPartitions,
    current_blocks: Option<Vec<Vec<usize>>>,
    perm: Vec<usize>,
    first_perm: bool,
}

impl OrderedSetPartitions {
    /// Create the iterator over all ordered set partitions of `{0, 1, ..., n-1}`.
    pub fn new(n: usize) -> Self {
        Self {
            set_partitions: SetPartitions::new(n),
            current_blocks: None,
            perm: Vec::new(),
            first_perm: false,
        }
    }
}

impl Iterator for OrderedSetPartitions {
    type Item = OrderedSetPartition;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(blocks) = &self.current_blocks {
                if self.first_perm {
                    self.first_perm = false;
                    return Some(permuted_ordered_partition(blocks, &self.perm));
                }
                if next_permutation(&mut self.perm) {
                    return Some(permuted_ordered_partition(blocks, &self.perm));
                }
            }

            let next_partition = self.set_partitions.next()?;
            let blocks = next_partition.into_blocks();
            self.perm = (0..blocks.len()).collect();
            self.current_blocks = Some(blocks);
            self.first_perm = true;
        }
    }
}

/// Return an iterator over all set partitions of `{0, 1, ..., n-1}`.
pub fn set_partitions(n: usize) -> SetPartitions {
    SetPartitions::new(n)
}

/// Return an iterator over all ordered set partitions of `{0, 1, ..., n-1}`.
pub fn ordered_set_partitions(n: usize) -> OrderedSetPartitions {
    OrderedSetPartitions::new(n)
}

fn normalize_blocks(blocks: &mut Vec<Vec<usize>>) {
    for block in blocks.iter_mut() {
        block.sort_unstable();
        assert!(
            !block.is_empty(),
            "set partitions cannot contain empty blocks"
        );
        for i in 1..block.len() {
            assert!(
                block[i - 1] != block[i],
                "duplicate element {} inside a block",
                block[i]
            );
        }
    }
}

fn validate_ground_set(blocks: &[Vec<usize>]) {
    let n: usize = blocks.iter().map(Vec::len).sum();
    if n == 0 {
        assert!(blocks.is_empty(), "the empty partition must have no blocks");
        return;
    }

    let mut seen = vec![false; n];
    for block in blocks {
        for &x in block {
            assert!(
                x < n,
                "element {x} is out of range for a partition of {{0, ..., {}}}",
                n - 1
            );
            assert!(!seen[x], "element {x} appears in more than one block");
            seen[x] = true;
        }
    }
    assert!(
        seen.into_iter().all(|flag| flag),
        "blocks do not cover the full ground set {{0, ..., {}}}",
        n - 1
    );
}

fn permuted_ordered_partition(blocks: &[Vec<usize>], perm: &[usize]) -> OrderedSetPartition {
    let ordered_blocks = perm.iter().map(|&i| blocks[i].clone()).collect();
    OrderedSetPartition {
        blocks: ordered_blocks,
    }
}

fn is_standard_restricted_growth_string(rgs: &[usize]) -> bool {
    if rgs.is_empty() {
        return true;
    }
    if rgs[0] != 0 {
        return false;
    }
    let mut max_seen = 0usize;
    for &x in &rgs[1..] {
        if x > max_seen + 1 {
            return false;
        }
        max_seen = max_seen.max(x);
    }
    true
}

fn next_standard_restricted_growth_string(rgs: &mut [usize], max_class: &mut [usize]) -> bool {
    if rgs.is_empty() {
        return false;
    }

    let mut i = rgs.len() - 1;
    loop {
        if i == 0 {
            return false;
        }
        let max_prev = max_class[i - 1];
        if rgs[i] < max_prev + 1 {
            rgs[i] += 1;
            max_class[i] = max_prev.max(rgs[i]);
            for j in i + 1..rgs.len() {
                rgs[j] = 0;
                max_class[j] = max_class[j - 1];
            }
            return true;
        }
        i -= 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_set_partition_rgs_roundtrip() {
        let p = SetPartition::new(vec![vec![2, 0], vec![1, 4], vec![3]]);
        assert_eq!(p.blocks(), &[vec![0, 2], vec![1, 4], vec![3]]);
        assert_eq!(p.to_restricted_growth_string(), vec![0, 1, 0, 2, 1]);
        assert_eq!(
            SetPartition::from_restricted_growth_string(&[0, 1, 0, 2, 1]),
            p
        );
    }

    #[test]
    fn test_ordered_set_partition_forget_order() {
        let p = OrderedSetPartition::new(vec![vec![2, 0], vec![3], vec![1, 4]]);
        assert_eq!(p.blocks(), &[vec![0, 2], vec![3], vec![1, 4]]);
        assert_eq!(
            p.as_set_partition(),
            SetPartition::new(vec![vec![0, 2], vec![1, 4], vec![3]])
        );
    }

    #[test]
    fn test_set_partition_counts_small() {
        let bell = [1usize, 1, 2, 5, 15, 52];
        for (n, &count) in bell.iter().enumerate() {
            assert_eq!(set_partitions(n).count(), count);
        }
    }

    #[test]
    fn test_ordered_set_partition_counts_small() {
        let fubini = [1usize, 1, 3, 13, 75, 541];
        for (n, &count) in fubini.iter().enumerate() {
            assert_eq!(ordered_set_partitions(n).count(), count);
        }
    }

    #[test]
    fn test_ordered_set_partitions_of_three() {
        let parts: Vec<Vec<Vec<usize>>> = ordered_set_partitions(3)
            .map(OrderedSetPartition::into_blocks)
            .collect();
        assert!(parts.contains(&vec![vec![0, 1, 2]]));
        assert!(parts.contains(&vec![vec![0], vec![1, 2]]));
        assert!(parts.contains(&vec![vec![1, 2], vec![0]]));
        assert!(parts.contains(&vec![vec![0], vec![1], vec![2]]));
        assert!(parts.contains(&vec![vec![2], vec![1], vec![0]]));
        assert_eq!(parts.len(), 13);
    }
}
