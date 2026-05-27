use flagged_lorentzian::{
    add_patterns, enumerate_tableaux, pair_envelope, sharp_flag, RowFlaggedSkewShape,
    SkewGtPattern, SkewShape, TableauRecord,
};

fn main() {
    let shape = SkewShape::from_parts(vec![6, 4], vec![2]);
    let alphabet = 5;
    let lower_label = 4;
    let beta = vec![0, 2, 4, 0, 0];
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

    let left_idx = ii_indices[2];
    let right_idx = jj_indices[5];
    let left_gt = SkewGtPattern::from_tableau(&shape, &tableaux[left_idx].values, alphabet);
    let right_gt = SkewGtPattern::from_tableau(&shape, &tableaux[right_idx].values, alphabet);
    let witness_sum = add_patterns(&left_gt, &right_gt);
    let witness_envelope = pair_envelope(
        &sharp_flag(&shape, &tableaux[left_idx].values),
        &sharp_flag(&shape, &tableaux[right_idx].values),
    );

    println!("two-row exact pair-sum witness");
    println!("shape=(6,4)/(2), beta={beta:?}, active labels=4,5");
    println!(
        "negative: A0={} / C2={}",
        format_tableau(&shape, &tableaux[left_idx].values),
        format_tableau(&shape, &tableaux[right_idx].values)
    );
    println!("negative envelope={witness_envelope:?}");
    println!("negative pair-sum:");
    print_rows(&witness_sum);
    println!();

    println!("mixed targets with same pair-sum and same envelope:");
    let mut count = 0usize;
    for (left_pos, &mixed_left_idx) in mixed_indices.iter().enumerate() {
        for (right_pos, &mixed_right_idx) in mixed_indices.iter().enumerate() {
            let mixed_left_gt =
                SkewGtPattern::from_tableau(&shape, &tableaux[mixed_left_idx].values, alphabet);
            let mixed_right_gt =
                SkewGtPattern::from_tableau(&shape, &tableaux[mixed_right_idx].values, alphabet);
            let sum = add_patterns(&mixed_left_gt, &mixed_right_gt);
            if sum != witness_sum {
                continue;
            }
            let envelope = pair_envelope(
                &sharp_flag(&shape, &tableaux[mixed_left_idx].values),
                &sharp_flag(&shape, &tableaux[mixed_right_idx].values),
            );
            if envelope != witness_envelope {
                continue;
            }
            count += 1;
            println!(
                "M{left_pos}/M{right_pos}: {} / {}",
                format_tableau(&shape, &tableaux[mixed_left_idx].values),
                format_tableau(&shape, &tableaux[mixed_right_idx].values)
            );
            println!("  left diff:");
            print_signed_rows(&subtract_rows(left_gt.rows(), mixed_left_gt.rows()));
            println!("  right diff:");
            print_signed_rows(&subtract_rows(right_gt.rows(), mixed_right_gt.rows()));
        }
    }
    println!("target_count={count}");
}

fn indices_with_content(tableaux: &[TableauRecord], content: &[u32]) -> Vec<usize> {
    tableaux
        .iter()
        .enumerate()
        .filter_map(|(idx, tableau)| (tableau.content.as_slice() == content).then_some(idx))
        .collect()
}

fn add_units(content: &[u32], index: usize, amount: u32) -> Vec<u32> {
    let mut out = content.to_vec();
    out[index] += amount;
    out
}

fn subtract_rows(left: &[Vec<u32>], right: &[Vec<u32>]) -> Vec<Vec<i32>> {
    left.iter()
        .zip(right)
        .map(|(left, right)| {
            left.iter()
                .zip(right)
                .map(|(&left, &right)| left as i32 - right as i32)
                .collect()
        })
        .collect()
}

fn print_rows(rows: &[Vec<u32>]) {
    for (level, row) in rows.iter().enumerate() {
        println!("  r{level}: {row:?}");
    }
}

fn print_signed_rows(rows: &[Vec<i32>]) {
    for (level, row) in rows.iter().enumerate() {
        if row.iter().any(|&entry| entry != 0) {
            println!("    r{level}: {row:?}");
        }
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
