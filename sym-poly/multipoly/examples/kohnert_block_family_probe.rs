use std::collections::BTreeMap;

use sym_poly_multipoly::{
    canonical_labeling, diagram_weight, kohnert_diagrams, Cell, Diagram, Labeling,
};

type ColumnType = Vec<usize>;
type LabeledState = Vec<(usize, usize)>;

struct BlockFamily {
    name: &'static str,
    left: Vec<ColumnType>,
    block: ColumnType,
    right: Vec<ColumnType>,
    max_n: usize,
}

fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    if args.get(1).map(String::as_str) == Some("--search") {
        let max_row = parse_arg(&args, 2, 4);
        let block_len = parse_arg(&args, 3, 4);
        let cell_cap = parse_arg(&args, 4, 12);
        let diagram_cap = parse_arg(&args, 5, 200_000);
        search_small_single_block_counterexamples(max_row, block_len, cell_cap, diagram_cap);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--show-front-drop") {
        show_first_suffix_front_lowering();
        return;
    }

    let families = [
        BlockFamily {
            name: "holey_135_with_right_boundary",
            left: vec![vec![1, 3]],
            block: vec![1, 3, 5],
            right: vec![vec![1, 5]],
            max_n: 5,
        },
        BlockFamily {
            name: "holey_124_with_low_boundary",
            left: vec![vec![1, 2]],
            block: vec![1, 2, 4],
            right: vec![vec![1, 4]],
            max_n: 5,
        },
        BlockFamily {
            name: "two_level_13_with_right_drop",
            left: vec![vec![1]],
            block: vec![1, 3],
            right: vec![vec![3]],
            max_n: 6,
        },
    ];

    for family in families {
        println!("== {} ==", family.name);
        for n in 1..=family.max_n {
            let (diagram, block_start, block_len) = build_family(&family, n);
            println!(
                "N={n}, columns={}, cells={}, southwest={}",
                max_col(&diagram),
                diagram.len(),
                is_southwest(&diagram)
            );
            match kohnert_diagrams(&diagram, 500_000) {
                Ok(diagrams) => {
                    let mut yamanouchi = Vec::new();
                    let mut weight_counts = BTreeMap::<Vec<u32>, usize>::new();
                    for candidate in diagrams {
                        let Some(labeling) = canonical_labeling(&diagram, &candidate) else {
                            continue;
                        };
                        let rectified = sym_poly_multipoly::rectify_labeled(&labeling, 1);
                        if rectified.iter().all(|(cell, label)| cell.row == *label) {
                            let weight = diagram_weight(&candidate);
                            *weight_counts.entry(weight.clone()).or_insert(0) += 1;
                            yamanouchi.push((weight, labeling));
                        }
                    }
                    println!("  AYK count: {}", yamanouchi.len());
                    for (weight, count) in &weight_counts {
                        println!("    {count} * {:?}", weight);
                    }
                    if yamanouchi.len() <= 12 {
                        for (idx, (weight, labeling)) in yamanouchi.iter().enumerate() {
                            let states = block_states(labeling, block_start, block_len);
                            println!(
                                "    Y{} wt={:?}, block states: {}",
                                idx + 1,
                                weight,
                                format_states(&states)
                            );
                        }
                    }
                }
                Err(err) => println!("  skipped: {err}"),
            }
        }
        println!();
    }
}

fn show_first_suffix_front_lowering() {
    let max_row = 4;
    let block_len = 4;
    let diagram_cap = 200_000;
    let cell_cap = 12;
    let all_types = column_types(max_row, true);
    let nonempty_types = column_types(max_row, false);

    for left in &all_types {
        for block in &nonempty_types {
            for right in &all_types {
                let family = BlockFamily {
                    name: "front-drop",
                    left: vec![left.clone()],
                    block: block.clone(),
                    right: vec![right.clone()],
                    max_n: block_len,
                };
                let (diagram, block_start, len) = build_family(&family, block_len);
                if diagram.len() > cell_cap || !is_southwest(&diagram) {
                    continue;
                }
                let Ok(diagrams) = kohnert_diagrams(&diagram, diagram_cap) else {
                    continue;
                };
                for candidate in diagrams {
                    let Some(labeling) = canonical_labeling(&diagram, &candidate) else {
                        continue;
                    };
                    for col in block_start..block_start + len - 1 {
                        let suffix = restrict_labeling(&labeling, col + 1);
                        let rectified = sym_poly_multipoly::rectify_labeled(&suffix, col + 1);
                        for label in block {
                            let Some(original_row) = row_of_label(&labeling, col + 1, *label)
                            else {
                                continue;
                            };
                            let Some(rectified_row) = row_of_label(&rectified, col + 1, *label)
                            else {
                                continue;
                            };
                            if rectified_row < original_row {
                                println!("First suffix-front lowering example:");
                                println!(
                                    "  left={left:?}, block={block:?}, right={right:?}, col={col}, label={label}"
                                );
                                println!(
                                    "  original row={original_row}, rectified row={rectified_row}"
                                );
                                println!(
                                    "  block states: {}",
                                    format_states(&block_states(&labeling, block_start, len))
                                );
                                println!(
                                    "  suffix front states before: {}",
                                    format_column(&labeling, col + 1)
                                );
                                println!(
                                    "  suffix front states after:  {}",
                                    format_column(&rectified, col + 1)
                                );
                                println!("  full labeling: {}", format_labeling(&labeling));
                                println!("  rectified suffix: {}", format_labeling(&rectified));
                                println!("  rectification trace:");
                                trace_rectification(suffix, col + 1);
                                return;
                            }
                        }
                    }
                }
            }
        }
    }
}

fn trace_rectification(mut labeling: Labeling, min_col: usize) {
    let mut step = 0usize;
    loop {
        let max_column = max_col_labeling(&labeling);
        let mut changed = false;
        for col in (min_col..max_column).rev() {
            if sym_poly_multipoly::rectify_labeled_column_star(&mut labeling, col) {
                step += 1;
                println!(
                    "    step {step}, star at col {col}: {}",
                    format_labeling(&labeling)
                );
                changed = true;
                break;
            }
        }
        if !changed {
            break;
        }
    }
}

fn parse_arg(args: &[String], index: usize, default: usize) -> usize {
    args.get(index)
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}

fn search_small_single_block_counterexamples(
    max_row: usize,
    block_len: usize,
    cell_cap: usize,
    diagram_cap: usize,
) {
    let mut tested = 0usize;
    let mut skipped_large = 0usize;
    let mut nonsouthwest = 0usize;
    let mut ayk_total = 0usize;
    let mut bad_examples = 0usize;
    let mut kohnert_total = 0usize;
    let mut nonmonotone_kohnert = 0usize;
    let mut non_rank_labeled_kohnert = 0usize;
    let mut suffix_front_lowering_kohnert = 0usize;
    let mut suffix_threshold_greedy_mismatch = 0usize;
    let all_types = column_types(max_row, true);
    let nonempty_types = column_types(max_row, false);

    for left in &all_types {
        for block in &nonempty_types {
            for right in &all_types {
                let family = BlockFamily {
                    name: "search",
                    left: vec![left.clone()],
                    block: block.clone(),
                    right: vec![right.clone()],
                    max_n: block_len,
                };
                let (diagram, block_start, len) = build_family(&family, block_len);
                if diagram.len() > cell_cap {
                    skipped_large += 1;
                    continue;
                }
                if !is_southwest(&diagram) {
                    nonsouthwest += 1;
                    continue;
                }
                tested += 1;
                let diagrams = match kohnert_diagrams(&diagram, diagram_cap) {
                    Ok(diagrams) => diagrams,
                    Err(_) => {
                        skipped_large += 1;
                        continue;
                    }
                };
                for candidate in diagrams {
                    let Some(labeling) = canonical_labeling(&diagram, &candidate) else {
                        continue;
                    };
                    kohnert_total += 1;
                    let states = block_states(&labeling, block_start, len);
                    if !is_labelwise_drop_monotone(&states) {
                        nonmonotone_kohnert += 1;
                    }
                    if !is_rank_labeled(&states, block) {
                        non_rank_labeled_kohnert += 1;
                    }
                    if has_suffix_front_lowering(&labeling, block_start, len, block) {
                        suffix_front_lowering_kohnert += 1;
                    }
                    if has_suffix_threshold_greedy_mismatch(&labeling, block_start, len, block) {
                        suffix_threshold_greedy_mismatch += 1;
                    }
                    let rectified = sym_poly_multipoly::rectify_labeled(&labeling, 1);
                    if !rectified.iter().all(|(cell, label)| cell.row == *label) {
                        continue;
                    }
                    ayk_total += 1;
                    let interior = &states[1..states.len() - 1];
                    if !interior.windows(2).all(|pair| pair[0] == pair[1]) {
                        bad_examples += 1;
                        println!("Potential counterexample to R=1 interior rigidity:");
                        println!("  left={left:?}, block={block:?}, right={right:?}");
                        println!("  weight={:?}", diagram_weight(&candidate));
                        println!("  block states: {}", format_states(&states));
                        if bad_examples >= 10 {
                            println!("Stopping after 10 potential counterexamples.");
                            return;
                        }
                    }
                }
            }
        }
    }

    println!("Small single-block search complete.");
    println!("  max_row={max_row}, block_len={block_len}, cell_cap={cell_cap}");
    println!("  tested southwest families: {tested}");
    println!("  nonsouthwest families skipped: {nonsouthwest}");
    println!("  large/capped families skipped: {skipped_large}");
    println!("  labeled Kohnert diagrams checked: {kohnert_total}");
    println!("  nonmonotone Kohnert block words: {nonmonotone_kohnert}");
    println!("  non-rank-labeled Kohnert block words: {non_rank_labeled_kohnert}");
    println!("  suffix-front lowering violations: {suffix_front_lowering_kohnert}");
    println!("  suffix-threshold greedy mismatches: {suffix_threshold_greedy_mismatch}");
    println!("  AYK diagrams checked: {ayk_total}");
    println!("  potential counterexamples: {bad_examples}");
}

fn build_family(family: &BlockFamily, n: usize) -> (Diagram, usize, usize) {
    let mut columns = family.left.clone();
    let block_start = columns.len() + 1;
    for _ in 0..n {
        columns.push(family.block.clone());
    }
    let block_len = n;
    columns.extend(family.right.clone());

    let mut diagram = Diagram::new();
    for (col_idx, rows) in columns.iter().enumerate() {
        for row in rows {
            diagram.insert(Cell {
                col: col_idx + 1,
                row: *row,
            });
        }
    }
    (diagram, block_start, block_len)
}

fn block_states(labeling: &Labeling, block_start: usize, block_len: usize) -> Vec<LabeledState> {
    (block_start..block_start + block_len)
        .map(|col| {
            let mut state = labeling
                .iter()
                .filter(|(cell, _)| cell.col == col)
                .map(|(cell, label)| (cell.row, *label))
                .collect::<Vec<_>>();
            state.sort();
            state
        })
        .collect()
}

fn format_states(states: &[LabeledState]) -> String {
    states
        .iter()
        .map(|state| {
            let entries = state
                .iter()
                .map(|(row, label)| format!("({row},{label})"))
                .collect::<Vec<_>>()
                .join(",");
            format!("{{{entries}}}")
        })
        .collect::<Vec<_>>()
        .join(" ")
}

fn format_column(labeling: &Labeling, col: usize) -> String {
    let mut state = labeling
        .iter()
        .filter(|(cell, _)| cell.col == col)
        .map(|(cell, label)| (cell.row, *label))
        .collect::<Vec<_>>();
    state.sort();
    format_states(&[state])
}

fn format_labeling(labeling: &Labeling) -> String {
    let mut entries = labeling
        .iter()
        .map(|(cell, label)| (cell.col, cell.row, *label))
        .collect::<Vec<_>>();
    entries.sort();
    entries
        .into_iter()
        .map(|(col, row, label)| format!("({col},{row};{label})"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn max_col(diagram: &Diagram) -> usize {
    diagram.iter().map(|cell| cell.col).max().unwrap_or(0)
}

fn max_col_labeling(labeling: &Labeling) -> usize {
    labeling.keys().map(|cell| cell.col).max().unwrap_or(0)
}

fn is_southwest(diagram: &Diagram) -> bool {
    for left in diagram {
        for right in diagram {
            if left.col < right.col && right.row < left.row {
                if !diagram.contains(&Cell {
                    col: left.col,
                    row: right.row,
                }) {
                    return false;
                }
            }
        }
    }
    true
}

fn is_labelwise_drop_monotone(states: &[LabeledState]) -> bool {
    for pair in states.windows(2) {
        for (_, label) in &pair[0] {
            let left_row = pair[0]
                .iter()
                .find_map(|(row, current_label)| {
                    if current_label == label {
                        Some(*row)
                    } else {
                        None
                    }
                })
                .expect("label appears in left state");
            let right_row = pair[1]
                .iter()
                .find_map(|(row, current_label)| {
                    if current_label == label {
                        Some(*row)
                    } else {
                        None
                    }
                })
                .expect("same source block has same labels");
            if left_row < right_row {
                return false;
            }
        }
    }
    true
}

fn is_rank_labeled(states: &[LabeledState], block: &[usize]) -> bool {
    let mut labels = block.to_vec();
    labels.sort();
    states.iter().all(|state| {
        let mut by_row = state.clone();
        by_row.sort();
        by_row
            .iter()
            .map(|(_, label)| *label)
            .eq(labels.iter().copied())
    })
}

fn has_suffix_front_lowering(
    labeling: &Labeling,
    block_start: usize,
    block_len: usize,
    block: &[usize],
) -> bool {
    for col in block_start..block_start + block_len - 1 {
        let suffix = restrict_labeling(labeling, col + 1);
        let rectified = sym_poly_multipoly::rectify_labeled(&suffix, col + 1);
        for label in block {
            let Some(original_row) = row_of_label(labeling, col + 1, *label) else {
                continue;
            };
            let Some(rectified_row) = row_of_label(&rectified, col + 1, *label) else {
                continue;
            };
            if rectified_row < original_row {
                return true;
            }
        }
    }
    false
}

fn has_suffix_threshold_greedy_mismatch(
    labeling: &Labeling,
    block_start: usize,
    block_len: usize,
    block: &[usize],
) -> bool {
    for col in block_start..block_start + block_len - 1 {
        let suffix = restrict_labeling(labeling, col + 1);
        let rectified = sym_poly_multipoly::rectify_labeled(&suffix, col + 1);
        let mut rows = labeling
            .iter()
            .filter_map(|(cell, _)| {
                if cell.col == col + 1 {
                    Some(cell.row)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        rows.sort();

        let Some(greedy) = greedy_from_thresholds(&rows, block, &rectified, col + 1) else {
            return true;
        };
        for (row, label) in greedy {
            if row_of_label(labeling, col + 1, label) != Some(row) {
                return true;
            }
        }
    }
    false
}

fn greedy_from_thresholds(
    rows: &[usize],
    block: &[usize],
    threshold_labeling: &Labeling,
    threshold_col: usize,
) -> Option<Vec<(usize, usize)>> {
    let mut labels = block.to_vec();
    let mut available = rows.to_vec();
    let mut result = Vec::new();
    labels.sort();
    for label in labels {
        let required_row = row_of_label(threshold_labeling, threshold_col, label);
        let position = available.iter().position(|row| {
            required_row
                .map(|required| *row >= required)
                .unwrap_or(true)
        })?;
        let row = available.remove(position);
        result.push((row, label));
    }
    Some(result)
}

fn restrict_labeling(labeling: &Labeling, min_col: usize) -> Labeling {
    labeling
        .iter()
        .filter(|(cell, _)| cell.col >= min_col)
        .map(|(cell, label)| (*cell, *label))
        .collect()
}

fn row_of_label(labeling: &Labeling, col: usize, label: usize) -> Option<usize> {
    labeling
        .iter()
        .filter_map(|(cell, current_label)| {
            if cell.col == col && *current_label == label {
                Some(cell.row)
            } else {
                None
            }
        })
        .min()
}

fn column_types(max_row: usize, include_empty: bool) -> Vec<ColumnType> {
    let start = if include_empty { 0 } else { 1 };
    (start..(1usize << max_row))
        .map(|mask| {
            (0..max_row)
                .filter_map(|idx| {
                    if mask & (1usize << idx) == 0 {
                        None
                    } else {
                        Some(idx + 1)
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect()
}
