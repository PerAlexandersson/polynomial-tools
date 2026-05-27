use std::collections::{BTreeMap, VecDeque};

use flagged_lorentzian::{
    add_patterns, descent_data_for_values, elementary_row_exchange_neighbors, enumerate_tableaux,
    pair_envelope, sharp_flag, subtract_pattern_sums, DescentData, DescentStatistic,
    RowFlaggedSkewShape, SkewGtPattern, SkewShape, TableauRecord,
};

type Content = Vec<u32>;
type GtSum = Vec<Vec<u32>>;
type GtRow = Vec<u32>;
type Envelope = Vec<u32>;

#[derive(Debug, Clone)]
struct SingleData {
    gt: SkewGtPattern,
    sharp_flag: Envelope,
    descent: DescentData,
}

#[derive(Debug, Clone)]
struct PairData {
    left_pos: usize,
    right_pos: usize,
    left_idx: usize,
    right_idx: usize,
    active_row: GtRow,
    envelope: Envelope,
    descent_pair: (DescentData, DescentData),
    gt_sum: GtSum,
}

fn main() {
    let shape = SkewShape::from_parts(vec![4, 4, 3, 2], vec![2, 2, 1]);
    let alphabet = 5;
    let lower_label = 4;
    let beta = vec![2, 3, 1, 0, 0];
    let bad_active = vec![8, 8, 6, 2];
    let bad_envelope = vec![3, 4, 4, 5];

    let flagged = RowFlaggedSkewShape::ordinary(shape.clone(), alphabet);
    let tableaux = enumerate_tableaux(&flagged, None).unwrap();
    let reading_orders = DescentStatistic::Global.reading_orders(&shape);
    let single_data: Vec<_> = tableaux
        .iter()
        .map(|tableau| SingleData {
            gt: SkewGtPattern::from_tableau(&shape, &tableau.values, alphabet),
            sharp_flag: sharp_flag(&shape, &tableau.values),
            descent: descent_data_for_values(&tableau.values, &reading_orders),
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

    let negative = pair_data(&single_data, &ii_indices, &jj_indices, lower_label);
    let positive = pair_data(&single_data, &mixed_indices, &mixed_indices, lower_label);

    let bad_negative: Vec<_> = negative
        .iter()
        .filter(|pair| pair.active_row == bad_active && pair.envelope == bad_envelope)
        .cloned()
        .collect();
    let exact_positive_count = positive
        .iter()
        .filter(|pair| pair.active_row == bad_active && pair.envelope == bad_envelope)
        .count();

    println!("first exact-envelope failure requiring envelope drop");
    println!("shape=(4,4,3,2)/(2,2,1), beta={beta:?}, active labels=4,5");
    println!("active={bad_active:?}, exact envelope={bad_envelope:?}");
    println!(
        "negative_count={}, exact_positive_count={exact_positive_count}",
        bad_negative.len()
    );

    let negative = bad_negative
        .first()
        .expect("expected the hardcoded exact-envelope witness");
    print_pair("negative", negative, &shape, &tableaux, &single_data);

    let mut candidates: Vec<_> = positive
        .iter()
        .filter(|pair| pair.active_row == negative.active_row)
        .filter(|pair| less_equal(&pair.envelope, &negative.envelope))
        .filter(|pair| pair.descent_pair == negative.descent_pair)
        .filter_map(|pair| {
            let path = exchange_path(&negative.gt_sum, &pair.gt_sum, lower_label as usize, 1)?;
            Some((pair.clone(), path))
        })
        .collect();
    candidates.sort_by_key(|(pair, _)| {
        let diff = subtract_pattern_sums(&negative.gt_sum, &pair.gt_sum);
        (
            pair.envelope.clone(),
            l1_distance(&diff),
            rows_changed(&diff),
            pair.gt_sum.clone(),
            pair.left_pos,
            pair.right_pos,
        )
    });

    println!();
    println!(
        "nonincrease candidates with same global descent pair and exchange depth <= 1: {}",
        candidates.len()
    );
    for (idx, (target, path)) in candidates.iter().take(8).enumerate() {
        let diff = subtract_pattern_sums(&negative.gt_sum, &target.gt_sum);
        println!();
        println!("candidate {idx}");
        print_pair("target", target, &shape, &tableaux, &single_data);
        println!(
            "envelope drop: {:?} -> {:?}",
            negative.envelope, target.envelope
        );
        println!("difference negative-target: {}", summarize_diff(&diff));
        println!("exchange path:");
        for (step, rows) in path.iter().enumerate() {
            println!("  step {step}");
            print_rows(rows);
        }
    }
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
                descent_pair: unordered_descent_pair(&left.descent, &right.descent),
                gt_sum,
            });
        }
    }
    out
}

fn print_pair(
    label: &str,
    pair: &PairData,
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    single_data: &[SingleData],
) {
    println!(
        "{label}: ({}, {}) = {} / {}",
        pair.left_pos,
        pair.right_pos,
        format_tableau(shape, &tableaux[pair.left_idx].values),
        format_tableau(shape, &tableaux[pair.right_idx].values)
    );
    println!(
        "  sharp flags: {:?}, {:?}; envelope={:?}; descents=({}, {})",
        single_data[pair.left_idx].sharp_flag,
        single_data[pair.right_idx].sharp_flag,
        pair.envelope,
        pair.descent_pair.0,
        pair.descent_pair.1
    );
    println!("  pair-sum rows:");
    print_rows(&pair.gt_sum);
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

fn unordered_descent_pair(left: &DescentData, right: &DescentData) -> (DescentData, DescentData) {
    if left <= right {
        (left.clone(), right.clone())
    } else {
        (right.clone(), left.clone())
    }
}

fn less_equal(left: &[u32], right: &[u32]) -> bool {
    left.iter().zip(right).all(|(&left, &right)| left <= right)
}

fn l1_distance(rows: &[Vec<i32>]) -> u32 {
    rows.iter()
        .flatten()
        .map(|entry| entry.unsigned_abs())
        .sum()
}

fn rows_changed(rows: &[Vec<i32>]) -> usize {
    rows.iter()
        .filter(|row| row.iter().any(|&entry| entry != 0))
        .count()
}

fn summarize_diff(diff: &[Vec<i32>]) -> String {
    let parts: Vec<_> = diff
        .iter()
        .enumerate()
        .filter(|(_, row)| row.iter().any(|&entry| entry != 0))
        .map(|(level, row)| format!("r{level}={row:?}"))
        .collect();
    if parts.is_empty() {
        "all rows preserved".to_string()
    } else {
        parts.join(", ")
    }
}

fn print_rows(rows: &[Vec<u32>]) {
    for (level, row) in rows.iter().enumerate() {
        println!("    r{level}: {row:?}");
    }
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
