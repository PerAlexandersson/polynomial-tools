//! Verify the staircase-rook <-> set-partition bijection and the induced
//! "big block" statistic on rook placements.
//!
//! A standard rook placement on the staircase board delta_{n-1} can be viewed
//! as a set of pairs (i,j) with 1 <= i < j <= n, no repeated i, and no repeated j.
//! Interpreting each rook as a directed edge i -> j gives a directed graph on [n]
//! with indegree/outdegree at most 1 and all edges increasing, hence a disjoint
//! union of directed paths. The vertex sets of those path components form a set
//! partition of [n].

use std::collections::{BTreeMap, BTreeSet};

type Placement = Vec<(usize, usize)>;
type Partition = Vec<Vec<usize>>;

fn staircase_placements(n: usize) -> Vec<Placement> {
    let mut used_cols = vec![false; n + 1];
    let mut current = Vec::new();
    let mut out = Vec::new();
    rec_placements(1, n, &mut used_cols, &mut current, &mut out);
    out
}

fn rec_placements(
    row: usize,
    n: usize,
    used_cols: &mut [bool],
    current: &mut Placement,
    out: &mut Vec<Placement>,
) {
    if row >= n {
        out.push(current.clone());
        return;
    }

    rec_placements(row + 1, n, used_cols, current, out);

    for col in (row + 1)..=n {
        if !used_cols[col] {
            used_cols[col] = true;
            current.push((row, col));
            rec_placements(row + 1, n, used_cols, current, out);
            current.pop();
            used_cols[col] = false;
        }
    }
}

fn placement_to_partition(n: usize, placement: &Placement) -> Partition {
    let mut next = vec![None; n + 1];
    let mut prev = vec![None; n + 1];
    for &(i, j) in placement {
        next[i] = Some(j);
        prev[j] = Some(i);
    }

    let mut blocks = Vec::new();
    for start in 1..=n {
        if prev[start].is_none() {
            let mut block = Vec::new();
            let mut v = start;
            block.push(v);
            while let Some(w) = next[v] {
                v = w;
                block.push(v);
            }
            blocks.push(block);
        }
    }
    blocks
}

fn partition_to_placement(partition: &Partition) -> Placement {
    let mut placement = Vec::new();
    for block in partition {
        for w in block.windows(2) {
            placement.push((w[0], w[1]));
        }
    }
    placement.sort();
    placement
}

fn set_partitions(n: usize) -> Vec<Partition> {
    if n == 0 {
        return vec![Vec::new()];
    }

    let mut rgs = vec![0usize; n];
    let mut max_class = vec![0usize; n];
    let mut out = Vec::new();

    loop {
        out.push(rgs_to_partition(&rgs));
        if !next_rgs(&mut rgs, &mut max_class) {
            break;
        }
    }

    out
}

fn rgs_to_partition(rgs: &[usize]) -> Partition {
    let num_classes = rgs.iter().copied().max().unwrap_or(0) + 1;
    let mut blocks = vec![Vec::new(); num_classes];
    for (i, &c) in rgs.iter().enumerate() {
        blocks[c].push(i + 1);
    }
    blocks
}

fn next_rgs(rgs: &mut [usize], max_class: &mut [usize]) -> bool {
    let n = rgs.len();
    let mut i = n - 1;
    loop {
        if i == 0 {
            return false;
        }
        let max_prev = max_class[i - 1];
        if rgs[i] < max_prev + 1 {
            rgs[i] += 1;
            max_class[i] = max_prev.max(rgs[i]);
            for j in (i + 1)..n {
                rgs[j] = 0;
                max_class[j] = max_class[j - 1];
            }
            return true;
        }
        i -= 1;
    }
}

fn big_blocks(partition: &Partition, j: usize) -> usize {
    partition.iter().filter(|block| block.len() >= j).count()
}

fn long_path_components(n: usize, placement: &Placement, j: usize) -> usize {
    placement_to_partition(n, placement)
        .iter()
        .filter(|block| block.len() >= j)
        .count()
}

fn main() {
    let max_n = 8;

    println!("=== Staircase rook <-> set partition bijection ===\n");

    for n in 1..=max_n {
        let placements = staircase_placements(n);
        let partitions = set_partitions(n);

        let placement_partitions: BTreeSet<Partition> = placements
            .iter()
            .map(|p| placement_to_partition(n, p))
            .collect();
        let partition_set: BTreeSet<Partition> = partitions.into_iter().collect();

        let image_ok = placement_partitions == partition_set;
        let inverse_ok = placements
            .iter()
            .all(|p| partition_to_placement(&placement_to_partition(n, p)) == {
                let mut q = p.clone();
                q.sort();
                q
            });

        println!(
            "n={}: placements={} partitions={} image_ok={} inverse_ok={}",
            n,
            placements.len(),
            partition_set.len(),
            image_ok,
            inverse_ok
        );
    }

    println!("\n=== Examples ===\n");
    let examples: Vec<Placement> = vec![
        vec![],
        vec![(1, 2)],
        vec![(1, 3), (3, 4)],
        vec![(1, 2), (2, 4)],
        vec![(1, 3), (2, 4)],
    ];
    for p in &examples {
        let n = p
            .iter()
            .flat_map(|&(i, j)| [i, j])
            .max()
            .unwrap_or(1)
            .max(4);
        let part = placement_to_partition(n, p);
        println!("placement={:?} -> partition={:?}", p, part);
    }

    println!("\n=== Big-block statistic on staircase rook placements ===\n");
    for j in 2..=5 {
        println!("j={j}:");
        for n in 1..=max_n {
            let placements = staircase_placements(n);
            let mut coeffs: BTreeMap<usize, usize> = BTreeMap::new();
            for p in &placements {
                let bb = long_path_components(n, p, j);
                *coeffs.entry(bb).or_insert(0) += 1;
            }
            let degree = coeffs.keys().copied().max().unwrap_or(0);
            let poly: Vec<usize> = (0..=degree).map(|k| coeffs.get(&k).copied().unwrap_or(0)).collect();
            println!("  n={}: {:?}", n, poly);
        }
        println!();
    }

    println!("=== Direct verification against set partitions ===\n");
    for j in 2..=5 {
        let mut all_ok = true;
        for n in 1..=max_n {
            let placements = staircase_placements(n);
            let partitions = set_partitions(n);
            let mut rook_counts: BTreeMap<usize, usize> = BTreeMap::new();
            let mut part_counts: BTreeMap<usize, usize> = BTreeMap::new();
            for p in &placements {
                *rook_counts.entry(long_path_components(n, p, j)).or_insert(0) += 1;
            }
            for pi in &partitions {
                *part_counts.entry(big_blocks(pi, j)).or_insert(0) += 1;
            }
            if rook_counts != part_counts {
                all_ok = false;
                println!(
                    "  FAIL: j={} n={} rook={:?} part={:?}",
                    j, n, rook_counts, part_counts
                );
            }
        }
        println!("j={}: {}", j, if all_ok { "all match ✓" } else { "FAIL" });
    }
}
