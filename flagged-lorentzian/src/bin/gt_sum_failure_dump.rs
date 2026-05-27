use flagged_lorentzian::{
    add_patterns, add_rows, enumerate_tableaux, subtract_pattern_sums, RowFlaggedSkewShape,
    SkewGtPattern, SkewShape, TableauRecord,
};

type Content = Vec<u32>;

fn main() {
    let shape = SkewShape::from_parts(vec![4, 2, 1], vec![2]);
    let alphabet = 5;
    let lower_label = 4;
    let beta = vec![1, 1, 1, 0, 0];
    let witness_gt_sum = vec![
        vec![4, 0, 0],
        vec![5, 1, 0],
        vec![6, 1, 1],
        vec![6, 3, 1],
        vec![6, 4, 2],
        vec![8, 4, 2],
    ];

    let flagged = RowFlaggedSkewShape::ordinary(shape.clone(), alphabet);
    let tableaux = enumerate_tableaux(&flagged, None).unwrap();
    let lower = lower_label as usize - 1;
    let upper = lower + 1;
    let ii = add_units(&beta, lower, 2);
    let mixed = add_units(&add_units(&beta, lower, 1), upper, 1);
    let jj = add_units(&beta, upper, 2);

    let ii_indices = indices_with_content(&tableaux, &ii);
    let mixed_indices = indices_with_content(&tableaux, &mixed);
    let jj_indices = indices_with_content(&tableaux, &jj);

    println!("first full-GT-sum failure");
    println!("shape=(4,2,1)/(2), beta={beta:?}, active labels=4,5");
    println!(
        "ii={}, mixed={}, jj={}",
        ii_indices.len(),
        mixed_indices.len(),
        jj_indices.len()
    );
    println!("witness pair GT sum:");
    print_unsigned_rows(&witness_gt_sum);
    println!();

    let mut witness_negative_pairs = Vec::new();
    for (left_pos, &left_idx) in ii_indices.iter().enumerate() {
        for (right_pos, &right_idx) in jj_indices.iter().enumerate() {
            let pair_sum = pair_gt_sum(&shape, &tableaux[left_idx], &tableaux[right_idx], alphabet);
            if pair_sum == witness_gt_sum {
                witness_negative_pairs.push((left_pos, right_pos, left_idx, right_idx));
            }
        }
    }

    let mut witness_positive_pairs = Vec::new();
    for (left_pos, &left_idx) in mixed_indices.iter().enumerate() {
        for (right_pos, &right_idx) in mixed_indices.iter().enumerate() {
            let pair_sum = pair_gt_sum(&shape, &tableaux[left_idx], &tableaux[right_idx], alphabet);
            if pair_sum == witness_gt_sum {
                witness_positive_pairs.push((left_pos, right_pos, left_idx, right_idx));
            }
        }
    }

    println!("negative pairs with witness full GT sum:");
    for (left_pos, right_pos, left_idx, right_idx) in &witness_negative_pairs {
        println!(
            "(A{left_pos}, C{right_pos}) = {}, {}",
            format_tableau(&shape, &tableaux[*left_idx].values),
            format_tableau(&shape, &tableaux[*right_idx].values)
        );
    }
    println!("positive pairs with witness full GT sum:");
    if witness_positive_pairs.is_empty() {
        println!("  none");
    } else {
        for (left_pos, right_pos, left_idx, right_idx) in &witness_positive_pairs {
            println!(
                "(M{left_pos}, M{right_pos}) = {}, {}",
                format_tableau(&shape, &tableaux[*left_idx].values),
                format_tableau(&shape, &tableaux[*right_idx].values)
            );
        }
    }

    let witness_active = witness_gt_sum[lower_label as usize].clone();
    println!();
    println!("closest active-row-compatible mixed targets:");
    let mut target_infos = Vec::new();
    for (left_pos, &left_idx) in mixed_indices.iter().enumerate() {
        for (right_pos, &right_idx) in mixed_indices.iter().enumerate() {
            let left_gt = SkewGtPattern::from_tableau(&shape, &tableaux[left_idx].values, alphabet);
            let right_gt =
                SkewGtPattern::from_tableau(&shape, &tableaux[right_idx].values, alphabet);
            let active_sum = add_rows(left_gt.row(lower_label), right_gt.row(lower_label));
            if active_sum != witness_active {
                continue;
            }
            let positive_sum = add_patterns(&left_gt, &right_gt);
            let difference = subtract_pattern_sums(&witness_gt_sum, &positive_sum);
            target_infos.push(TargetInfo {
                left_pos,
                right_pos,
                left_idx,
                right_idx,
                l1_distance: l1_distance(&difference),
                rows_changed: rows_changed(&difference),
                summary: summarize_signed_rows(&difference),
            });
        }
    }
    target_infos.sort_by_key(|info| {
        (
            info.l1_distance,
            info.rows_changed,
            info.left_pos,
            info.right_pos,
        )
    });
    for info in target_infos.iter().take(16) {
        println!(
            "(M{}, M{}) l1={} rows={} :: {} / {} :: {}",
            info.left_pos,
            info.right_pos,
            info.l1_distance,
            info.rows_changed,
            format_tableau(&shape, &tableaux[info.left_idx].values),
            format_tableau(&shape, &tableaux[info.right_idx].values),
            info.summary
        );
    }
}

#[derive(Debug)]
struct TargetInfo {
    left_pos: usize,
    right_pos: usize,
    left_idx: usize,
    right_idx: usize,
    l1_distance: u32,
    rows_changed: usize,
    summary: String,
}

fn pair_gt_sum(
    shape: &SkewShape,
    left: &TableauRecord,
    right: &TableauRecord,
    alphabet: usize,
) -> Vec<Vec<u32>> {
    let left_gt = SkewGtPattern::from_tableau(shape, &left.values, alphabet);
    let right_gt = SkewGtPattern::from_tableau(shape, &right.values, alphabet);
    add_patterns(&left_gt, &right_gt)
}

fn print_unsigned_rows(rows: &[Vec<u32>]) {
    for (level, row) in rows.iter().enumerate() {
        println!("  r{level}: {row:?}");
    }
}

fn summarize_signed_rows(rows: &[Vec<i32>]) -> String {
    let parts: Vec<_> = rows
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
