use flagged_lorentzian::{
    add_patterns, add_rows, descent_data_for_values, enumerate_tableaux, pair_envelope, sharp_flag,
    subtract_pattern_sums, DescentStatistic, RowFlaggedSkewShape, SkewGtPattern, SkewShape,
    TableauRecord,
};

type Content = Vec<u32>;

#[derive(Debug, Clone)]
struct SingleData {
    gt: SkewGtPattern,
    sharp_flag: Vec<u32>,
    descent: flagged_lorentzian::DescentData,
}

#[derive(Debug, Clone)]
struct PairData {
    left_pos: usize,
    right_pos: usize,
    left_idx: usize,
    right_idx: usize,
    active_row: Vec<u32>,
    envelope: Vec<u32>,
    gt_sum: Vec<Vec<u32>>,
    descent_pair: (
        flagged_lorentzian::DescentData,
        flagged_lorentzian::DescentData,
    ),
}

fn main() {
    let shape = SkewShape::from_parts(vec![2, 2, 2, 2, 1], vec![2, 2, 1]);
    let row_flags = vec![3, 4, 5, 6, 6];
    let alphabet = 6;
    let lower_label = 4;
    let beta = vec![0, 0, 1, 0, 0, 1];

    let flagged = RowFlaggedSkewShape::new(shape.clone(), row_flags, alphabet);
    let tableaux = enumerate_tableaux(&flagged, None).unwrap();
    let reading_orders = DescentStatistic::Componentwise.reading_orders(&shape);
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

    println!("first lower-label active-row failure");
    println!("shape=(2,2,2,2,1)/(2,2,1), flags=(3,4,5,6,6), beta={beta:?}, active labels=4,5");
    println!(
        "ii={}, mixed={}, jj={}, negative_pairs={}, positive_pairs={}",
        ii_indices.len(),
        mixed_indices.len(),
        jj_indices.len(),
        negative.len(),
        positive.len()
    );

    let witness = negative
        .iter()
        .find(|pair| {
            !positive
                .iter()
                .any(|target| target.active_row == pair.active_row)
        })
        .expect("expected active-row witness");
    print_pair("negative witness", witness, &shape, &tableaux, &single_data);

    let mut candidates = positive.clone();
    candidates.sort_by_key(|target| {
        let active_diff = signed_row_diff(&witness.active_row, &target.active_row);
        let full_diff = subtract_pattern_sums(&witness.gt_sum, &target.gt_sum);
        (
            l1_row(&active_diff),
            target.envelope.clone(),
            l1_rows(&full_diff),
            rows_changed(&full_diff),
            target.left_pos,
            target.right_pos,
        )
    });

    println!();
    println!("closest mixed targets by active-row distance:");
    for target in candidates.iter().take(12) {
        let active_diff = signed_row_diff(&witness.active_row, &target.active_row);
        let full_diff = subtract_pattern_sums(&witness.gt_sum, &target.gt_sum);
        print_pair("target", target, &shape, &tableaux, &single_data);
        println!("  active difference witness-target: {active_diff:?}");
        println!("  full difference: {}", summarize_diff(&full_diff));
        println!();
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
            let descent_pair = if left.descent <= right.descent {
                (left.descent.clone(), right.descent.clone())
            } else {
                (right.descent.clone(), left.descent.clone())
            };
            out.push(PairData {
                left_pos,
                right_pos,
                left_idx,
                right_idx,
                active_row: add_rows(left.gt.row(lower_label), right.gt.row(lower_label)),
                envelope: pair_envelope(&left.sharp_flag, &right.sharp_flag),
                gt_sum,
                descent_pair,
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
        "  active={:?}, envelope={:?}, sharp flags={:?}, {:?}, descents=({}, {})",
        pair.active_row,
        pair.envelope,
        single_data[pair.left_idx].sharp_flag,
        single_data[pair.right_idx].sharp_flag,
        pair.descent_pair.0,
        pair.descent_pair.1
    );
    println!("  pair-sum rows:");
    for (level, row) in pair.gt_sum.iter().enumerate() {
        println!("    r{level}: {row:?}");
    }
}

fn signed_row_diff(left: &[u32], right: &[u32]) -> Vec<i32> {
    left.iter()
        .zip(right)
        .map(|(&left, &right)| left as i32 - right as i32)
        .collect()
}

fn l1_row(row: &[i32]) -> u32 {
    row.iter().map(|entry| entry.unsigned_abs()).sum()
}

fn l1_rows(rows: &[Vec<i32>]) -> u32 {
    rows.iter().map(|row| l1_row(row)).sum()
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
