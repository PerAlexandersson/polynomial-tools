use std::collections::BTreeMap;

use flagged_lorentzian::{
    enumerate_tableaux, is_gt_array, RowFlaggedSkewShape, SkewGtPattern, SkewShape, TableauRecord,
};

type GtRows = Vec<Vec<u32>>;

fn main() {
    let shape = SkewShape::from_parts(vec![6, 4], vec![2]);
    let alphabet = 5;
    let lower_label = 4usize;
    let beta = vec![0, 0, 1, 2, 3];
    let flagged = RowFlaggedSkewShape::ordinary(shape.clone(), alphabet);
    let tableaux = enumerate_tableaux(&flagged, None).unwrap();
    let gt_rows: Vec<_> = tableaux
        .iter()
        .map(|tableau| {
            SkewGtPattern::from_tableau(&shape, &tableau.values, alphabet)
                .rows()
                .to_vec()
        })
        .collect();

    let lower = lower_label - 1;
    let upper = lower + 1;
    let ii = add_units(&beta, lower, 2);
    let mixed = add_units(&add_units(&beta, lower, 1), upper, 1);
    let jj = add_units(&beta, upper, 2);
    let ii_indices = indices_with_content(&tableaux, &ii);
    let mixed_indices = indices_with_content(&tableaux, &mixed);
    let jj_indices = indices_with_content(&tableaux, &jj);
    let mixed_gt_to_pos: BTreeMap<_, _> = mixed_indices
        .iter()
        .enumerate()
        .map(|(pos, &idx)| (gt_rows[idx].clone(), pos))
        .collect();

    println!("carrier collision dump");
    println!("shape=(6,4)/(2), beta={beta:?}, active labels=4,5");
    println!(
        "ii={}, mixed={}, jj={}",
        ii_indices.len(),
        mixed_indices.len(),
        jj_indices.len()
    );
    println!();

    println!("ii tableaux:");
    for (pos, &idx) in ii_indices.iter().enumerate() {
        println!("A{pos}: {}", format_tableau(&shape, &tableaux[idx].values));
    }
    println!("jj tableaux:");
    for (pos, &idx) in jj_indices.iter().enumerate() {
        println!("C{pos}: {}", format_tableau(&shape, &tableaux[idx].values));
    }
    println!("mixed tableaux:");
    for (pos, &idx) in mixed_indices.iter().enumerate() {
        println!("M{pos}: {}", format_tableau(&shape, &tableaux[idx].values));
    }
    println!();

    for (left_pos, right_pos) in [(1usize, 1usize), (2, 0)] {
        let left_idx = ii_indices[left_pos];
        let right_idx = jj_indices[right_pos];
        println!(
            "source (A{left_pos}, C{right_pos}) = {} / {}",
            format_tableau(&shape, &tableaux[left_idx].values),
            format_tableau(&shape, &tableaux[right_idx].values)
        );
        let candidates = carrier_targets(
            &gt_rows[left_idx],
            &gt_rows[right_idx],
            &mixed_gt_to_pos,
            lower_label,
        );
        for candidate in candidates {
            println!(
                "  -> (M{}, M{}) via {}",
                candidate.left_pos,
                candidate.right_pos,
                format_carrier(&candidate.carrier)
            );
        }
    }
}

#[derive(Debug, Clone)]
struct CarrierTarget {
    left_pos: usize,
    right_pos: usize,
    carrier: Vec<[i32; 2]>,
}

fn carrier_targets(
    left: &[Vec<u32>],
    right: &[Vec<u32>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    lower_label: usize,
) -> Vec<CarrierTarget> {
    let mut out = Vec::new();
    for carrier in carrier_differences(left.len(), lower_label) {
        let Some(left_image) = apply_carrier(left, &carrier, -1) else {
            continue;
        };
        let Some(&left_pos) = mixed_gt_to_pos.get(&left_image) else {
            continue;
        };
        let Some(right_image) = apply_carrier(right, &carrier, 1) else {
            continue;
        };
        let Some(&right_pos) = mixed_gt_to_pos.get(&right_image) else {
            continue;
        };
        out.push(CarrierTarget {
            left_pos,
            right_pos,
            carrier,
        });
    }
    out
}

fn carrier_differences(row_count: usize, lower_label: usize) -> Vec<Vec<[i32; 2]>> {
    let mut carriers = Vec::new();
    for start_level in (1..=lower_label).rev() {
        let mut carrier = vec![[0, 0]; row_count];
        extend_carrier(&mut carriers, &mut carrier, start_level, lower_label);
    }
    carriers
}

fn extend_carrier(
    carriers: &mut Vec<Vec<[i32; 2]>>,
    carrier: &mut Vec<[i32; 2]>,
    level: usize,
    lower_label: usize,
) {
    if level == lower_label {
        carrier[level] = [1, 0];
        carriers.push(carrier.clone());
        carrier[level] = [0, 1];
        carriers.push(carrier.clone());
        carrier[level] = [0, 0];
        return;
    }
    carrier[level] = [1, -1];
    extend_carrier(carriers, carrier, level + 1, lower_label);
    carrier[level] = [-1, 1];
    extend_carrier(carriers, carrier, level + 1, lower_label);
    carrier[level] = [0, 0];
}

fn apply_carrier(rows: &[Vec<u32>], carrier: &[[i32; 2]], sign: i32) -> Option<GtRows> {
    let mut out = rows.to_vec();
    for (row, delta) in out.iter_mut().zip(carrier) {
        for col in 0..2 {
            let value = row[col] as i32 + sign * delta[col];
            if value < 0 {
                return None;
            }
            row[col] = value as u32;
        }
    }
    is_gt_array(&out).then_some(out)
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

fn format_carrier(carrier: &[[i32; 2]]) -> String {
    carrier
        .iter()
        .enumerate()
        .filter(|(_, delta)| **delta != [0, 0])
        .map(|(level, delta)| format!("D{level}=({},{})", delta[0], delta[1]))
        .collect::<Vec<_>>()
        .join(", ")
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
