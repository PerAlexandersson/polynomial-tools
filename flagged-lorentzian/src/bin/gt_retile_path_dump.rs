use std::collections::{BTreeMap, VecDeque};

use flagged_lorentzian::{
    add_patterns, elementary_row_exchange_neighbors, enumerate_tableaux, gt, pair_envelope,
    sharp_flag, RowFlaggedSkewShape, SkewGtPattern, SkewShape, TableauRecord,
};

type Content = Vec<u32>;
type GtSum = Vec<Vec<u32>>;
type ActiveRow = Vec<u32>;
type Envelope = Vec<u32>;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct FullKey {
    active_row: ActiveRow,
    envelope: Envelope,
    gt_sum: GtSum,
}

#[derive(Debug, Clone)]
struct SingleData {
    gt: SkewGtPattern,
    sharp_flag: Envelope,
}

#[derive(Debug, Clone)]
struct PairData {
    left_pos: usize,
    right_pos: usize,
    left_idx: usize,
    right_idx: usize,
    active_row: ActiveRow,
    envelope: Envelope,
    gt_sum: GtSum,
}

impl PairData {
    fn full_key(&self) -> FullKey {
        FullKey {
            active_row: self.active_row.clone(),
            envelope: self.envelope.clone(),
            gt_sum: self.gt_sum.clone(),
        }
    }
}

fn main() {
    let shape = SkewShape::from_parts(vec![4, 4, 3, 2], vec![3, 2, 2]);
    let alphabet = 5;
    let lower_label = 4;
    let beta = vec![1, 2, 1, 0, 0];
    let max_depth = 3;
    let flagged = RowFlaggedSkewShape::ordinary(shape.clone(), alphabet);
    let tableaux = enumerate_tableaux(&flagged, None).unwrap();
    let single_data: Vec<_> = tableaux
        .iter()
        .map(|tableau| SingleData {
            gt: SkewGtPattern::from_tableau(&shape, &tableau.values, alphabet),
            sharp_flag: sharp_flag(&shape, &tableau.values),
        })
        .collect();

    let lower = lower_label as usize - 1;
    let upper = lower + 1;
    let ii = add_units(&beta, lower, 2);
    let mixed = add_units(&add_units(&beta, lower, 1), upper, 1);
    let jj = add_units(&beta, upper, 2);
    let ii_indices = indices_with_content(&tableaux, &ii);
    let mixed_indices = indices_with_content(&tableaux, &mixed);
    let jj_indices = indices_with_content(&tableaux, &jj);

    let negative_pairs = pair_data(&single_data, &ii_indices, &jj_indices, lower_label);
    let positive_pairs = pair_data(&single_data, &mixed_indices, &mixed_indices, lower_label);
    let mut negative_by_key = BTreeMap::<FullKey, Vec<PairData>>::new();
    let mut positive_by_key = BTreeMap::<FullKey, Vec<PairData>>::new();
    for pair in negative_pairs {
        negative_by_key
            .entry(pair.full_key())
            .or_default()
            .push(pair);
    }
    for pair in positive_pairs.clone() {
        positive_by_key
            .entry(pair.full_key())
            .or_default()
            .push(pair);
    }

    for (key, negative_preimages) in negative_by_key {
        let positive_count = positive_by_key.get(&key).map_or(0, Vec::len);
        if negative_preimages.len() <= positive_count {
            continue;
        }

        let mut candidate_targets: Vec<_> = positive_pairs
            .iter()
            .filter(|pair| pair.active_row == key.active_row && pair.envelope == key.envelope)
            .filter(|pair| pair.gt_sum != key.gt_sum)
            .cloned()
            .collect();
        candidate_targets.sort_by_key(|pair| {
            let distance = distance_between(&key.gt_sum, &pair.gt_sum);
            (
                distance.l1,
                distance.rows_changed,
                pair.left_pos,
                pair.right_pos,
            )
        });
        candidate_targets.dedup_by_key(|pair| pair.gt_sum.clone());

        let Some(target) = candidate_targets.first() else {
            continue;
        };
        let distance = distance_between(&key.gt_sum, &target.gt_sum);
        if distance.l1 < 4 {
            continue;
        }

        let Some(path) =
            exchange_path(&key.gt_sum, &target.gt_sum, lower_label as usize, max_depth)
        else {
            println!("No path found within depth {max_depth}");
            return;
        };

        println!("depth-{} retile witness", path.len().saturating_sub(1));
        println!("shape=(4,4,3,2)/(3,2,2), beta={beta:?}, active labels=4,5");
        println!("active={:?}, envelope={:?}", key.active_row, key.envelope);
        println!(
            "negative_count={}, positive_count_same_full_sum={positive_count}",
            negative_preimages.len()
        );
        let negative = &negative_preimages[0];
        println!(
            "negative pair: (A{}, C{}) = {} / {}",
            negative.left_pos,
            negative.right_pos,
            format_tableau(&shape, &tableaux[negative.left_idx].values),
            format_tableau(&shape, &tableaux[negative.right_idx].values)
        );
        println!(
            "target pair: (M{}, M{}) = {} / {}",
            target.left_pos,
            target.right_pos,
            format_tableau(&shape, &tableaux[target.left_idx].values),
            format_tableau(&shape, &tableaux[target.right_idx].values)
        );
        println!(
            "difference negative-target: {}",
            summarize_diff(&distance.diff)
        );
        println!();
        println!("exchange path:");
        for (step, array) in path.iter().enumerate() {
            let mixed_splits = positive_pairs
                .iter()
                .filter(|pair| pair.active_row == key.active_row && pair.envelope == key.envelope)
                .filter(|pair| &pair.gt_sum == array)
                .count();
            println!("step {step}: mixed_splits={mixed_splits}");
            print_rows(array);
        }
        return;
    }

    println!("No L1>=4 deficit found.");
}

#[derive(Debug, Clone)]
struct Distance {
    l1: u32,
    rows_changed: usize,
    diff: Vec<Vec<i32>>,
}

fn distance_between(left: &[Vec<u32>], right: &[Vec<u32>]) -> Distance {
    let diff = gt::subtract_pattern_sums(left, right);
    Distance {
        l1: diff
            .iter()
            .flatten()
            .map(|entry| entry.unsigned_abs())
            .sum(),
        rows_changed: diff
            .iter()
            .filter(|row| row.iter().any(|&entry| entry != 0))
            .count(),
        diff,
    }
}

fn exchange_path(
    start: &[Vec<u32>],
    target: &[Vec<u32>],
    fixed_level: usize,
    max_depth: usize,
) -> Option<Vec<GtSum>> {
    let start = start.to_vec();
    let target = target.to_vec();
    let mut parent = BTreeMap::<GtSum, Option<GtSum>>::from([(start.clone(), None)]);
    let mut queue = VecDeque::from([(start.clone(), 0usize)]);

    while let Some((array, depth)) = queue.pop_front() {
        if array == target {
            let mut path = Vec::new();
            let mut current = Some(array);
            while let Some(node) = current {
                current = parent.get(&node).cloned().flatten();
                path.push(node);
            }
            path.reverse();
            return Some(path);
        }
        if depth == max_depth {
            continue;
        }
        for neighbor in elementary_row_exchange_neighbors(&array, fixed_level) {
            if parent.contains_key(&neighbor) {
                continue;
            }
            parent.insert(neighbor.clone(), Some(array.clone()));
            queue.push_back((neighbor, depth + 1));
        }
    }

    None
}

fn pair_data(
    single_data: &[SingleData],
    left_indices: &[usize],
    right_indices: &[usize],
    lower_label: u32,
) -> Vec<PairData> {
    let mut out = Vec::new();
    for (left_pos, &left_idx) in left_indices.iter().enumerate() {
        for (right_pos, &right_idx) in right_indices.iter().enumerate() {
            let left = &single_data[left_idx];
            let right = &single_data[right_idx];
            let gt_sum = add_patterns(&left.gt, &right.gt);
            out.push(PairData {
                left_pos,
                right_pos,
                left_idx,
                right_idx,
                active_row: gt_sum[lower_label as usize].clone(),
                envelope: pair_envelope(&left.sharp_flag, &right.sharp_flag),
                gt_sum,
            });
        }
    }
    out
}

fn indices_with_content(tableaux: &[TableauRecord], content: &[u32]) -> Vec<usize> {
    tableaux
        .iter()
        .enumerate()
        .filter_map(|(idx, tableau)| (tableau.content.as_slice() == content).then_some(idx))
        .collect()
}

fn add_units(content: &[u32], index: usize, amount: u32) -> Content {
    let mut out = content.to_vec();
    out[index] += amount;
    out
}

fn print_rows(rows: &[Vec<u32>]) {
    for (level, row) in rows.iter().enumerate() {
        println!("  r{level}: {row:?}");
    }
}

fn summarize_diff(diff: &[Vec<i32>]) -> String {
    let parts: Vec<_> = diff
        .iter()
        .enumerate()
        .filter(|(_, row)| row.iter().any(|&entry| entry != 0))
        .map(|(level, row)| format!("r{level}={row:?}"))
        .collect();
    if parts.is_empty() {
        "0".to_string()
    } else {
        parts.join("; ")
    }
}

fn format_tableau(shape: &SkewShape, values: &[u32]) -> String {
    let rows = shape
        .cells()
        .iter()
        .map(|cell| cell.row)
        .max()
        .map_or(0, |row| row + 1);
    let cols = shape
        .cells()
        .iter()
        .map(|cell| cell.col)
        .max()
        .map_or(0, |col| col + 1);

    let mut out = Vec::new();
    for row in 0..rows {
        let mut row_values = Vec::new();
        for col in 0..cols {
            let value = shape
                .cell_index(flagged_lorentzian::Cell { row, col })
                .map(|idx| values[idx].to_string())
                .unwrap_or_else(|| ".".to_string());
            row_values.push(value);
        }
        out.push(row_values.join(" "));
    }
    format!("[{}]", out.join(" / "))
}
