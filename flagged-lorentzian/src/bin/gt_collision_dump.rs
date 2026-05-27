use flagged_lorentzian::{
    add_patterns, add_rows, enumerate_tableaux, pair_envelope, sharp_flag, subtract_pattern_sums,
    RowFlaggedSkewShape, SkewGtPattern, SkewShape, TableauRecord,
};

type Content = Vec<u32>;

fn main() {
    let shape = SkewShape::from_parts(vec![6, 2], vec![2]);
    let alphabet = 5;
    let flagged = RowFlaggedSkewShape::ordinary(shape.clone(), alphabet);
    let tableaux = enumerate_tableaux(&flagged, None).unwrap();

    let beta = vec![0, 0, 2, 0, 2];
    let lower_label = 4;
    let lower = lower_label as usize - 1;
    let upper = lower + 1;
    let ii = add_units(&beta, lower, 2);
    let mixed = add_units(&add_units(&beta, lower, 1), upper, 1);
    let jj = add_units(&beta, upper, 2);

    let ii_indices = indices_with_content(&tableaux, &ii);
    let mixed_indices = indices_with_content(&tableaux, &mixed);
    let jj_indices = indices_with_content(&tableaux, &jj);

    println!("shape=(6,2)/(2), beta={beta:?}, active labels=4,5");
    println!(
        "ii={}, mixed={}, jj={}",
        ii_indices.len(),
        mixed_indices.len(),
        jj_indices.len()
    );
    println!();

    print_family("A", &shape, &tableaux, &ii_indices, alphabet, lower_label);
    print_family("C", &shape, &tableaux, &jj_indices, alphabet, lower_label);
    print_family(
        "M",
        &shape,
        &tableaux,
        &mixed_indices,
        alphabet,
        lower_label,
    );

    let a2 = &tableaux[ii_indices[2]];
    let c2 = &tableaux[jj_indices[2]];
    let m4 = &tableaux[mixed_indices[4]];

    println!();
    println!("selected missing active-row target:");
    println!("negative: (A2, C2)");
    println!("positive: (M4, M4)");
    println!();
    print_selected_pattern("A2", &shape, a2, alphabet);
    print_selected_pattern("C2", &shape, c2, alphabet);
    print_selected_pattern("M4", &shape, m4, alphabet);

    let a2_gt = SkewGtPattern::from_tableau(&shape, &a2.values, alphabet);
    let c2_gt = SkewGtPattern::from_tableau(&shape, &c2.values, alphabet);
    let m4_gt = SkewGtPattern::from_tableau(&shape, &m4.values, alphabet);
    let negative_sum = add_patterns(&a2_gt, &c2_gt);
    let positive_sum = add_patterns(&m4_gt, &m4_gt);
    let difference = subtract_pattern_sums(&negative_sum, &positive_sum);

    println!("pair GT row sums, negative (A2+C2):");
    print_unsigned_rows(&negative_sum);
    println!("pair GT row sums, positive (M4+M4):");
    print_unsigned_rows(&positive_sum);
    println!("difference negative - positive:");
    print_signed_rows(&difference);

    println!();
    println!("all active-row-compatible mixed targets for (A2, C2):");
    let negative_active = add_rows(a2_gt.row(lower_label), c2_gt.row(lower_label));
    for (left_pos, &left_idx) in mixed_indices.iter().enumerate() {
        for (right_pos, &right_idx) in mixed_indices.iter().enumerate() {
            let left_gt = SkewGtPattern::from_tableau(&shape, &tableaux[left_idx].values, alphabet);
            let right_gt =
                SkewGtPattern::from_tableau(&shape, &tableaux[right_idx].values, alphabet);
            let positive_active = add_rows(left_gt.row(lower_label), right_gt.row(lower_label));
            if positive_active != negative_active {
                continue;
            }

            let positive_sum = add_patterns(&left_gt, &right_gt);
            let difference = subtract_pattern_sums(&negative_sum, &positive_sum);
            println!(
                "(M{left_pos}, M{right_pos}) lower-difference-summary: {}",
                summarize_signed_rows(&difference)
            );
        }
    }
}

fn print_family(
    name: &str,
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    indices: &[usize],
    alphabet: usize,
    lower_label: u32,
) {
    println!("{name} tableaux:");
    for (pos, &idx) in indices.iter().enumerate() {
        let tableau = &tableaux[idx];
        let gt = SkewGtPattern::from_tableau(shape, &tableau.values, alphabet);
        println!(
            "{name}{pos}: {}  active={:?}  sharp={:?}",
            format_tableau(shape, &tableau.values),
            gt.row(lower_label),
            sharp_flag(shape, &tableau.values)
        );
    }
    println!();
}

fn print_selected_pattern(name: &str, shape: &SkewShape, tableau: &TableauRecord, alphabet: usize) {
    let gt = SkewGtPattern::from_tableau(shape, &tableau.values, alphabet);
    println!("{name}: {}", format_tableau(shape, &tableau.values));
    println!("  sharp={:?}", sharp_flag(shape, &tableau.values));
    println!(
        "  envelope with itself={:?}",
        pair_envelope(
            &sharp_flag(shape, &tableau.values),
            &sharp_flag(shape, &tableau.values),
        )
    );
    print_unsigned_rows(gt.rows());
}

fn print_unsigned_rows(rows: &[Vec<u32>]) {
    for (level, row) in rows.iter().enumerate() {
        println!("  r{level}: {row:?}");
    }
}

fn print_signed_rows(rows: &[Vec<i32>]) {
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
