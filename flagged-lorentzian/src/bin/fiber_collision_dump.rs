use std::collections::{BTreeMap, BTreeSet};

use flagged_lorentzian::{
    bender_knuth_e_images, bender_knuth_f_images, enumerate_tableaux, RowFlaggedSkewShape,
    SkewShape, TableauRecord,
};

type Content = Vec<u32>;

fn main() {
    let shape = SkewShape::from_parts(vec![6, 2], vec![2]);
    let flagged = RowFlaggedSkewShape::ordinary(shape.clone(), 5);
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
    let mixed_value_to_pos: BTreeMap<_, _> = mixed_indices
        .iter()
        .enumerate()
        .map(|(pos, &idx)| (tableaux[idx].values.clone(), pos))
        .collect();

    let left_images: Vec<Vec<usize>> = ii_indices
        .iter()
        .map(|&idx| {
            image_positions(
                bender_knuth_f_images(&shape, &tableaux[idx].values, lower_label),
                &mixed_value_to_pos,
            )
        })
        .collect();
    let right_images: Vec<Vec<usize>> = jj_indices
        .iter()
        .map(|&idx| {
            image_positions(
                bender_knuth_e_images(&shape, &tableaux[idx].values, lower_label),
                &mixed_value_to_pos,
            )
        })
        .collect();

    println!("shape=(6,2)/(2), beta={beta:?}, active labels=4,5");
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
        println!("    f-targets: {:?}", left_images[pos]);
    }
    println!();

    println!("jj tableaux:");
    for (pos, &idx) in jj_indices.iter().enumerate() {
        println!("C{pos}: {}", format_tableau(&shape, &tableaux[idx].values));
        println!("    e-targets: {:?}", right_images[pos]);
    }
    println!();

    println!("mixed tableaux:");
    for (pos, &idx) in mixed_indices.iter().enumerate() {
        println!("M{pos}: {}", format_tableau(&shape, &tableaux[idx].values));
    }
    println!();

    let mut target_preimages = BTreeMap::<(usize, usize), Vec<(usize, usize)>>::new();
    let mut negative_edges = BTreeMap::<(usize, usize), Vec<(usize, usize)>>::new();

    for (left_pos, left_targets) in left_images.iter().enumerate() {
        for (right_pos, right_targets) in right_images.iter().enumerate() {
            let mut targets = BTreeSet::new();
            for &left_target in left_targets {
                for &right_target in right_targets {
                    let target = (left_target, right_target);
                    targets.insert(target);
                    target_preimages
                        .entry(target)
                        .or_default()
                        .push((left_pos, right_pos));
                }
            }
            negative_edges.insert((left_pos, right_pos), targets.into_iter().collect());
        }
    }

    println!("negative pair edges:");
    for (negative, targets) in &negative_edges {
        println!(
            "(A{}, C{}) -> {}",
            negative.0,
            negative.1,
            format_targets(targets)
        );
    }
    println!();

    println!("mixed targets with multiple preimages:");
    for (target, preimages) in target_preimages
        .iter()
        .filter(|(_, preimages)| preimages.len() > 1)
    {
        println!(
            "(M{}, M{}) <- {}",
            target.0,
            target.1,
            format_preimages(preimages)
        );
    }
}

fn indices_with_content(tableaux: &[TableauRecord], content: &[u32]) -> Vec<usize> {
    tableaux
        .iter()
        .enumerate()
        .filter_map(|(idx, tableau)| (tableau.content.as_slice() == content).then_some(idx))
        .collect()
}

fn image_positions(
    images: Vec<Vec<u32>>,
    mixed_value_to_pos: &BTreeMap<Vec<u32>, usize>,
) -> Vec<usize> {
    let mut positions: Vec<_> = images
        .into_iter()
        .filter_map(|image| mixed_value_to_pos.get(&image).copied())
        .collect();
    positions.sort_unstable();
    positions.dedup();
    positions
}

fn add_units(content: &[u32], index: usize, amount: u32) -> Content {
    let mut out = content.to_vec();
    out[index] += amount;
    out
}

fn format_targets(targets: &[(usize, usize)]) -> String {
    let parts: Vec<_> = targets
        .iter()
        .map(|(left, right)| format!("(M{left},M{right})"))
        .collect();
    format!("[{}]", parts.join(", "))
}

fn format_preimages(preimages: &[(usize, usize)]) -> String {
    let parts: Vec<_> = preimages
        .iter()
        .map(|(left, right)| format!("(A{left},C{right})"))
        .collect();
    parts.join(", ")
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
