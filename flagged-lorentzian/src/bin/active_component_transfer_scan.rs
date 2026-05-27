use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use flagged_lorentzian::{
    active_component_crystal_e_images, active_component_crystal_f_images, bender_knuth_e_images,
    bender_knuth_f_images, enumerate_tableaux, families::skew_shapes_of_size, RowFlaggedSkewShape,
    SkewShape, TableauRecord,
};

type Content = Vec<u32>;

#[derive(Debug, Clone)]
struct Args {
    lambda: Option<Vec<u32>>,
    mu: Vec<u32>,
    alphabet: usize,
    lower_label: u32,
    operator: OperatorMode,
    max_skew_size: u32,
    max_outer_extra: u32,
    connected_only: bool,
    tableau_limit: Option<usize>,
    stop_at_first_failure: bool,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            lambda: None,
            mu: Vec::new(),
            alphabet: 5,
            lower_label: 1,
            operator: OperatorMode::ActiveComponentCrystal,
            max_skew_size: 7,
            max_outer_extra: 8,
            connected_only: false,
            tableau_limit: None,
            stop_at_first_failure: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OperatorMode {
    ActiveComponentCrystal,
    BenderKnuthFreeCell,
    RightmostCrackedColumn,
    AnyCrackedColumn,
}

fn main() {
    let args = match Args::parse_from_env() {
        Ok(args) => args,
        Err(message) => {
            eprintln!("{message}");
            eprintln!();
            eprintln!("{}", help_text());
            std::process::exit(2);
        }
    };
    let started = Instant::now();

    if let Some(lambda) = args.lambda.clone() {
        let shape = SkewShape::from_parts(lambda, args.mu.clone());
        let flagged_shape = RowFlaggedSkewShape::ordinary(shape, args.alphabet);
        let outcome = scan_one_shape(&flagged_shape, &args);
        print_single_outcome(&flagged_shape, &outcome, started.elapsed().as_secs_f64());
        std::process::exit(if outcome.limit_exceeded { 2 } else { 0 });
    }

    let mut shapes_checked = 0usize;
    let mut fibers_checked = 0usize;
    let mut total_negative_pairs = 0u128;
    let mut total_matched_pairs = 0u128;
    let mut total_edges = 0u128;
    let mut skipped_by_limit = 0usize;

    for skew_size in 0..=args.max_skew_size {
        let max_outer_size = skew_size + args.max_outer_extra;
        for shape in skew_shapes_of_size(skew_size, max_outer_size) {
            if args.connected_only && !shape.is_connected() {
                continue;
            }

            let flagged_shape = RowFlaggedSkewShape::ordinary(shape, args.alphabet);
            let outcome = scan_one_shape(&flagged_shape, &args);
            if outcome.limit_exceeded {
                skipped_by_limit += 1;
                continue;
            }

            shapes_checked += 1;
            fibers_checked += outcome.fibers_checked;
            total_negative_pairs += outcome.negative_pairs;
            total_matched_pairs += outcome.matched_pairs;
            total_edges += outcome.edge_count;

            if args.stop_at_first_failure && outcome.first_failure.is_some() {
                println!("FAIL");
                println!("outer={:?}", flagged_shape.shape().outer().parts());
                println!("inner={:?}", flagged_shape.shape().inner().parts());
                println!("{}", outcome.first_failure.unwrap());
                println!(
                    "shapes_checked={shapes_checked}, fibers_checked={fibers_checked}, total_negative_pairs={total_negative_pairs}, total_matched_pairs={total_matched_pairs}, total_edges={total_edges}, skipped_by_limit={skipped_by_limit}, elapsed={:.3}s",
                    started.elapsed().as_secs_f64()
                );
                std::process::exit(1);
            }
        }

        println!(
            "skew_size={skew_size}: shapes_checked={shapes_checked}, fibers_checked={fibers_checked}, total_negative_pairs={total_negative_pairs}, total_matched_pairs={total_matched_pairs}, total_edges={total_edges}, skipped_by_limit={skipped_by_limit}, elapsed={:.3}s",
            started.elapsed().as_secs_f64()
        );
    }

    println!("PASS");
    println!(
        "shapes_checked={shapes_checked}, fibers_checked={fibers_checked}, total_negative_pairs={total_negative_pairs}, total_matched_pairs={total_matched_pairs}, total_edges={total_edges}, skipped_by_limit={skipped_by_limit}, elapsed={:.3}s",
        started.elapsed().as_secs_f64()
    );
}

impl Args {
    fn parse_from_env() -> Result<Self, String> {
        let mut args = Args::default();
        let mut iter = std::env::args().skip(1);
        while let Some(flag) = iter.next() {
            match flag.as_str() {
                "--help" | "-h" => {
                    println!("{}", help_text());
                    std::process::exit(0);
                }
                "--lambda" => args.lambda = Some(parse_u32_vec(&take_value(&mut iter, &flag)?)?),
                "--mu" => args.mu = parse_u32_vec(&take_value(&mut iter, &flag)?)?,
                "--alphabet" => {
                    args.alphabet = take_value(&mut iter, &flag)?
                        .parse()
                        .map_err(|err| format!("invalid --alphabet: {err}"))?
                }
                "--lower-label" => {
                    args.lower_label = take_value(&mut iter, &flag)?
                        .parse()
                        .map_err(|err| format!("invalid --lower-label: {err}"))?
                }
                "--operator" => {
                    args.operator = match take_value(&mut iter, &flag)?.as_str() {
                        "active-component-crystal" => OperatorMode::ActiveComponentCrystal,
                        "bender-knuth-free-cell" => OperatorMode::BenderKnuthFreeCell,
                        "rightmost-cracked-column" => OperatorMode::RightmostCrackedColumn,
                        "any-cracked-column" => OperatorMode::AnyCrackedColumn,
                        other => {
                            return Err(format!(
                                "invalid --operator `{other}`; expected active-component-crystal, bender-knuth-free-cell, rightmost-cracked-column, or any-cracked-column"
                            ))
                        }
                    }
                }
                "--max-skew-size" => {
                    args.max_skew_size = take_value(&mut iter, &flag)?
                        .parse()
                        .map_err(|err| format!("invalid --max-skew-size: {err}"))?
                }
                "--max-outer-extra" => {
                    args.max_outer_extra = take_value(&mut iter, &flag)?
                        .parse()
                        .map_err(|err| format!("invalid --max-outer-extra: {err}"))?
                }
                "--connected-only" => args.connected_only = true,
                "--tableau-limit" => {
                    args.tableau_limit = Some(
                        take_value(&mut iter, &flag)?
                            .parse()
                            .map_err(|err| format!("invalid --tableau-limit: {err}"))?,
                    )
                }
                "--stop-at-first-failure" => args.stop_at_first_failure = true,
                other => return Err(format!("unknown argument `{other}`")),
            }
        }
        if args.lower_label as usize >= args.alphabet {
            return Err("--lower-label must be smaller than --alphabet".to_string());
        }
        Ok(args)
    }
}

#[derive(Debug)]
struct ShapeScanOutcome {
    tableaux: usize,
    fibers_checked: usize,
    negative_pairs: u128,
    positive_pairs: u128,
    matched_pairs: u128,
    edge_count: u128,
    first_failure: Option<MatchingFailure>,
    limit_exceeded: bool,
}

#[derive(Debug, Clone)]
struct MatchingFailure {
    beta: Content,
    negative_pairs: u128,
    positive_pairs: u128,
    matched_pairs: u128,
    edge_count: u128,
    isolated_negative_pairs: u128,
    first_isolated_pair: Option<IsolatedPair>,
}

#[derive(Debug, Clone)]
struct IsolatedPair {
    left_values: Vec<u32>,
    right_values: Vec<u32>,
    left_image_count: usize,
    right_image_count: usize,
}

impl std::fmt::Display for MatchingFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "beta={:?}, negative_pairs={}, positive_pairs={}, matched_pairs={}, edge_count={}, isolated_negative_pairs={}",
            self.beta,
            self.negative_pairs,
            self.positive_pairs,
            self.matched_pairs,
            self.edge_count,
            self.isolated_negative_pairs
        )?;
        if let Some(pair) = &self.first_isolated_pair {
            write!(
                f,
                ", first_isolated_left={:?}, first_isolated_right={:?}, left_image_count={}, right_image_count={}",
                pair.left_values,
                pair.right_values,
                pair.left_image_count,
                pair.right_image_count
            )?;
        }
        Ok(())
    }
}

fn scan_one_shape(flagged_shape: &RowFlaggedSkewShape, args: &Args) -> ShapeScanOutcome {
    let tableaux = match enumerate_tableaux(flagged_shape, args.tableau_limit) {
        Ok(tableaux) => tableaux,
        Err(_) => {
            return ShapeScanOutcome {
                tableaux: 0,
                fibers_checked: 0,
                negative_pairs: 0,
                positive_pairs: 0,
                matched_pairs: 0,
                edge_count: 0,
                first_failure: None,
                limit_exceeded: true,
            }
        }
    };

    let mut by_content = BTreeMap::<Content, Vec<usize>>::new();
    for (idx, tableau) in tableaux.iter().enumerate() {
        by_content
            .entry(tableau.content.clone())
            .or_default()
            .push(idx);
    }

    let lower = args.lower_label as usize - 1;
    let upper = lower + 1;
    let mut seen_beta = BTreeSet::new();
    let mut fibers_checked = 0usize;
    let mut negative_pairs = 0u128;
    let mut positive_pairs = 0u128;
    let mut matched_pairs = 0u128;
    let mut edge_count = 0u128;
    let mut first_failure = None;

    for content in by_content.keys() {
        let Some(beta) = subtract_units(content, lower, 2) else {
            continue;
        };
        let ii = add_units(&beta, lower, 2);
        let mixed = add_units(&add_units(&beta, lower, 1), upper, 1);
        let jj = add_units(&beta, upper, 2);
        if content != &ii || !seen_beta.insert(beta.clone()) {
            continue;
        }

        let (Some(ii_indices), Some(mixed_indices), Some(jj_indices)) = (
            by_content.get(&ii),
            by_content.get(&mixed),
            by_content.get(&jj),
        ) else {
            continue;
        };

        fibers_checked += 1;
        let result = scan_fiber(
            flagged_shape.shape(),
            &tableaux,
            &beta,
            ii_indices,
            mixed_indices,
            jj_indices,
            args.lower_label,
            args.operator,
        );
        negative_pairs += result.negative_pairs;
        positive_pairs += result.positive_pairs;
        matched_pairs += result.matched_pairs;
        edge_count += result.edge_count;

        if result.matched_pairs < result.negative_pairs && first_failure.is_none() {
            first_failure = Some(result);
        }
    }

    ShapeScanOutcome {
        tableaux: tableaux.len(),
        fibers_checked,
        negative_pairs,
        positive_pairs,
        matched_pairs,
        edge_count,
        first_failure,
        limit_exceeded: false,
    }
}

fn scan_fiber(
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    beta: &[u32],
    ii_indices: &[usize],
    mixed_indices: &[usize],
    jj_indices: &[usize],
    lower_label: u32,
    operator: OperatorMode,
) -> MatchingFailure {
    let mut mixed_value_to_pos = BTreeMap::new();
    for (pos, &idx) in mixed_indices.iter().enumerate() {
        mixed_value_to_pos.insert(tableaux[idx].values.clone(), pos);
    }

    let left_images: Vec<Vec<usize>> = ii_indices
        .iter()
        .map(|&idx| {
            f_image_positions(
                shape,
                &tableaux[idx].values,
                lower_label,
                operator,
                tableaux,
                mixed_indices,
                &mixed_value_to_pos,
            )
        })
        .collect();
    let right_images: Vec<Vec<usize>> = jj_indices
        .iter()
        .map(|&idx| {
            e_image_positions(
                shape,
                &tableaux[idx].values,
                lower_label,
                operator,
                tableaux,
                mixed_indices,
                &mixed_value_to_pos,
            )
        })
        .collect();

    let right_count = mixed_indices.len() * mixed_indices.len();
    let mut edges = vec![Vec::new(); ii_indices.len() * jj_indices.len()];
    let mut edge_count = 0u128;
    let mut isolated_negative_pairs = 0u128;
    let mut first_isolated_pair = None;

    for (left_pos, left_targets) in left_images.iter().enumerate() {
        for (right_pos, right_targets) in right_images.iter().enumerate() {
            let negative_pos = left_pos * jj_indices.len() + right_pos;
            let mut targets = BTreeSet::new();
            for &left_target in left_targets {
                for &right_target in right_targets {
                    targets.insert(left_target * mixed_indices.len() + right_target);
                }
            }
            if targets.is_empty() {
                isolated_negative_pairs += 1;
                if first_isolated_pair.is_none() {
                    first_isolated_pair = Some(IsolatedPair {
                        left_values: tableaux[ii_indices[left_pos]].values.clone(),
                        right_values: tableaux[jj_indices[right_pos]].values.clone(),
                        left_image_count: left_targets.len(),
                        right_image_count: right_targets.len(),
                    });
                }
            }
            edge_count += targets.len() as u128;
            edges[negative_pos] = targets.into_iter().collect();
        }
    }

    let matched_pairs = max_bipartite_matching(&edges, right_count) as u128;
    MatchingFailure {
        beta: beta.to_vec(),
        negative_pairs: ii_indices.len() as u128 * jj_indices.len() as u128,
        positive_pairs: mixed_indices.len() as u128 * mixed_indices.len() as u128,
        matched_pairs,
        edge_count,
        isolated_negative_pairs,
        first_isolated_pair,
    }
}

fn f_images(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
    operator: OperatorMode,
) -> Vec<Vec<u32>> {
    match operator {
        OperatorMode::ActiveComponentCrystal => {
            active_component_crystal_f_images(shape, values, lower_label)
        }
        OperatorMode::BenderKnuthFreeCell => bender_knuth_f_images(shape, values, lower_label),
        OperatorMode::RightmostCrackedColumn => bender_knuth_f_images(shape, values, lower_label),
        OperatorMode::AnyCrackedColumn => bender_knuth_f_images(shape, values, lower_label),
    }
}

fn e_images(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
    operator: OperatorMode,
) -> Vec<Vec<u32>> {
    match operator {
        OperatorMode::ActiveComponentCrystal => {
            active_component_crystal_e_images(shape, values, lower_label)
        }
        OperatorMode::BenderKnuthFreeCell => bender_knuth_e_images(shape, values, lower_label),
        OperatorMode::RightmostCrackedColumn => bender_knuth_e_images(shape, values, lower_label),
        OperatorMode::AnyCrackedColumn => bender_knuth_e_images(shape, values, lower_label),
    }
}

fn f_image_positions(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
    operator: OperatorMode,
    tableaux: &[TableauRecord],
    mixed_indices: &[usize],
    mixed_value_to_pos: &BTreeMap<Vec<u32>, usize>,
) -> Vec<usize> {
    let mut positions = image_positions(
        f_images(shape, values, lower_label, operator),
        mixed_value_to_pos,
    );
    if operator == OperatorMode::RightmostCrackedColumn {
        positions.extend(cracked_column_f_positions(
            shape,
            values,
            lower_label,
            tableaux,
            mixed_indices,
        ));
        positions.sort_unstable();
        positions.dedup();
    } else if operator == OperatorMode::AnyCrackedColumn {
        positions.extend(any_cracked_column_f_positions(
            shape,
            values,
            lower_label,
            tableaux,
            mixed_indices,
        ));
        positions.sort_unstable();
        positions.dedup();
    }
    positions
}

fn e_image_positions(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
    operator: OperatorMode,
    tableaux: &[TableauRecord],
    mixed_indices: &[usize],
    mixed_value_to_pos: &BTreeMap<Vec<u32>, usize>,
) -> Vec<usize> {
    let mut positions = image_positions(
        e_images(shape, values, lower_label, operator),
        mixed_value_to_pos,
    );
    if operator == OperatorMode::RightmostCrackedColumn {
        positions.extend(cracked_column_e_positions(
            shape,
            values,
            lower_label,
            tableaux,
            mixed_indices,
        ));
        positions.sort_unstable();
        positions.dedup();
    } else if operator == OperatorMode::AnyCrackedColumn {
        positions.extend(any_cracked_column_e_positions(
            shape,
            values,
            lower_label,
            tableaux,
            mixed_indices,
        ));
        positions.sort_unstable();
        positions.dedup();
    }
    positions
}

fn any_cracked_column_f_positions(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
    tableaux: &[TableauRecord],
    mixed_indices: &[usize],
) -> Vec<usize> {
    let upper_label = lower_label + 1;
    let mut positions = Vec::new();
    for (top_idx, cracked_col) in vertical_pairs(shape, values, lower_label) {
        positions.extend(mixed_indices.iter().enumerate().filter_map(
            |(mixed_pos, &tableau_idx)| {
                let candidate = &tableaux[tableau_idx].values;
                cracked_f_candidate(
                    shape,
                    values,
                    candidate,
                    top_idx,
                    cracked_col,
                    lower_label,
                    upper_label,
                )
                .then_some(mixed_pos)
            },
        ));
    }
    positions
}

fn any_cracked_column_e_positions(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
    tableaux: &[TableauRecord],
    mixed_indices: &[usize],
) -> Vec<usize> {
    let upper_label = lower_label + 1;
    let mut positions = Vec::new();
    for (bottom_idx, cracked_col) in vertical_pair_bottoms(shape, values, lower_label) {
        positions.extend(mixed_indices.iter().enumerate().filter_map(
            |(mixed_pos, &tableau_idx)| {
                let candidate = &tableaux[tableau_idx].values;
                cracked_e_candidate(
                    shape,
                    values,
                    candidate,
                    bottom_idx,
                    cracked_col,
                    lower_label,
                    upper_label,
                )
                .then_some(mixed_pos)
            },
        ));
    }
    positions
}

fn cracked_column_f_positions(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
    tableaux: &[TableauRecord],
    mixed_indices: &[usize],
) -> Vec<usize> {
    let Some((top_idx, cracked_col)) = rightmost_vertical_pair(shape, values, lower_label) else {
        return Vec::new();
    };
    let upper_label = lower_label + 1;
    mixed_indices
        .iter()
        .enumerate()
        .filter_map(|(mixed_pos, &tableau_idx)| {
            let candidate = &tableaux[tableau_idx].values;
            cracked_f_candidate(
                shape,
                values,
                candidate,
                top_idx,
                cracked_col,
                lower_label,
                upper_label,
            )
            .then_some(mixed_pos)
        })
        .collect()
}

fn cracked_column_e_positions(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
    tableaux: &[TableauRecord],
    mixed_indices: &[usize],
) -> Vec<usize> {
    let Some((bottom_idx, cracked_col)) = leftmost_vertical_pair_bottom(shape, values, lower_label)
    else {
        return Vec::new();
    };
    let upper_label = lower_label + 1;
    mixed_indices
        .iter()
        .enumerate()
        .filter_map(|(mixed_pos, &tableau_idx)| {
            let candidate = &tableaux[tableau_idx].values;
            cracked_e_candidate(
                shape,
                values,
                candidate,
                bottom_idx,
                cracked_col,
                lower_label,
                upper_label,
            )
            .then_some(mixed_pos)
        })
        .collect()
}

fn cracked_f_candidate(
    shape: &SkewShape,
    values: &[u32],
    candidate: &[u32],
    top_idx: usize,
    cracked_col: usize,
    lower_label: u32,
    upper_label: u32,
) -> bool {
    candidate[top_idx] == upper_label
        && shape.cells().iter().enumerate().all(|(idx, cell)| {
            if cell.col < cracked_col {
                candidate[idx] == values[idx]
            } else if idx != top_idx && values[idx] == lower_label {
                candidate[idx] == lower_label
            } else {
                true
            }
        })
}

fn cracked_e_candidate(
    shape: &SkewShape,
    values: &[u32],
    candidate: &[u32],
    bottom_idx: usize,
    cracked_col: usize,
    lower_label: u32,
    upper_label: u32,
) -> bool {
    candidate[bottom_idx] == lower_label
        && shape.cells().iter().enumerate().all(|(idx, cell)| {
            if cell.col > cracked_col {
                candidate[idx] == values[idx]
            } else if idx != bottom_idx && values[idx] == upper_label {
                candidate[idx] == upper_label
            } else {
                true
            }
        })
}

fn rightmost_vertical_pair(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
) -> Option<(usize, usize)> {
    let upper_label = lower_label + 1;
    shape
        .cells()
        .iter()
        .enumerate()
        .filter_map(|(idx, cell)| {
            shape
                .cell_index(flagged_lorentzian::Cell {
                    row: cell.row + 1,
                    col: cell.col,
                })
                .filter(|&below_idx| values[idx] == lower_label && values[below_idx] == upper_label)
                .map(|_| (idx, cell.col))
        })
        .max_by_key(|(_, col)| *col)
}

fn vertical_pairs(shape: &SkewShape, values: &[u32], lower_label: u32) -> Vec<(usize, usize)> {
    let upper_label = lower_label + 1;
    shape
        .cells()
        .iter()
        .enumerate()
        .filter_map(|(idx, cell)| {
            shape
                .cell_index(flagged_lorentzian::Cell {
                    row: cell.row + 1,
                    col: cell.col,
                })
                .filter(|&below_idx| values[idx] == lower_label && values[below_idx] == upper_label)
                .map(|_| (idx, cell.col))
        })
        .collect()
}

fn leftmost_vertical_pair_bottom(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
) -> Option<(usize, usize)> {
    let upper_label = lower_label + 1;
    shape
        .cells()
        .iter()
        .enumerate()
        .filter_map(|(idx, cell)| {
            shape
                .cell_index(flagged_lorentzian::Cell {
                    row: cell.row + 1,
                    col: cell.col,
                })
                .filter(|&below_idx| values[idx] == lower_label && values[below_idx] == upper_label)
                .map(|below_idx| (below_idx, cell.col))
        })
        .min_by_key(|(_, col)| *col)
}

fn vertical_pair_bottoms(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
) -> Vec<(usize, usize)> {
    let upper_label = lower_label + 1;
    shape
        .cells()
        .iter()
        .enumerate()
        .filter_map(|(idx, cell)| {
            shape
                .cell_index(flagged_lorentzian::Cell {
                    row: cell.row + 1,
                    col: cell.col,
                })
                .filter(|&below_idx| values[idx] == lower_label && values[below_idx] == upper_label)
                .map(|below_idx| (below_idx, cell.col))
        })
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

fn max_bipartite_matching(edges: &[Vec<usize>], right_count: usize) -> usize {
    let mut match_right = vec![usize::MAX; right_count];
    let mut seen = vec![0usize; right_count];
    let mut matched = 0usize;

    for left in 0..edges.len() {
        if augment(left, edges, &mut match_right, &mut seen, left + 1) {
            matched += 1;
        }
    }

    matched
}

fn augment(
    left: usize,
    edges: &[Vec<usize>],
    match_right: &mut [usize],
    seen: &mut [usize],
    stamp: usize,
) -> bool {
    for &right in &edges[left] {
        if seen[right] == stamp {
            continue;
        }
        seen[right] = stamp;
        if match_right[right] == usize::MAX
            || augment(match_right[right], edges, match_right, seen, stamp)
        {
            match_right[right] = left;
            return true;
        }
    }
    false
}

fn subtract_units(content: &[u32], index: usize, amount: u32) -> Option<Content> {
    let mut out = content.to_vec();
    *out.get_mut(index)? = out[index].checked_sub(amount)?;
    Some(out)
}

fn add_units(content: &[u32], index: usize, amount: u32) -> Content {
    let mut out = content.to_vec();
    out[index] += amount;
    out
}

fn print_single_outcome(
    flagged_shape: &RowFlaggedSkewShape,
    outcome: &ShapeScanOutcome,
    elapsed_seconds: f64,
) {
    if outcome.limit_exceeded {
        println!("SKIPPED: tableau limit exceeded");
        return;
    }

    println!("outer={:?}", flagged_shape.shape().outer().parts());
    println!("inner={:?}", flagged_shape.shape().inner().parts());
    println!("alphabet={}", flagged_shape.alphabet());
    println!("tableaux={}", outcome.tableaux);
    println!("fibers_checked={}", outcome.fibers_checked);
    println!("negative_pairs={}", outcome.negative_pairs);
    println!("positive_pairs={}", outcome.positive_pairs);
    println!("matched_pairs={}", outcome.matched_pairs);
    println!("edge_count={}", outcome.edge_count);
    println!("elapsed={elapsed_seconds:.3}s");
    if let Some(failure) = &outcome.first_failure {
        println!("first_failure={failure}");
    } else {
        println!("all_negative_pairs_match_by_active_component_moves");
    }
}

fn take_value(iter: &mut impl Iterator<Item = String>, flag: &str) -> Result<String, String> {
    iter.next()
        .ok_or_else(|| format!("missing value after {flag}"))
}

fn parse_u32_vec(input: &str) -> Result<Vec<u32>, String> {
    if input.trim().is_empty() {
        return Ok(Vec::new());
    }

    input
        .split(',')
        .map(|part| {
            part.trim()
                .parse::<u32>()
                .map_err(|err| format!("invalid integer `{part}`: {err}"))
        })
        .collect()
}

fn help_text() -> &'static str {
    "Check whether local active moves admit a 2x2 matching.

For each negative pair (A,C), the scanner allows a local f_i move on A and a
local e_i move on C.  It then checks whether the induced bipartite graph from
negative pairs to mixed pairs has a matching covering all negative pairs.

USAGE:
  active_component_transfer_scan [OPTIONS]

OPTIONS:
  --lambda PARTS              Outer partition, e.g. 4,3,1. If omitted, scan a family.
  --mu PARTS                  Inner partition, e.g. 3,1. Default: empty.
  --alphabet N                Alphabet size. Default: 5.
  --lower-label I             Use crystal operators on labels I,I+1. Default: 1.
  --operator KIND             active-component-crystal, bender-knuth-free-cell,
                              rightmost-cracked-column, or any-cracked-column.
                              Default: active-component-crystal.
  --max-skew-size N           Maximum skew size for family scans. Default: 7.
  --max-outer-extra N         Family outer sizes go up to skew_size + this. Default: 8.
  --connected-only            Restrict family scans to connected skew shapes.
  --tableau-limit N           Skip any shape whose tableau enumeration exceeds N.
  --stop-at-first-failure     Stop as soon as active-component moves cannot match a fiber.
  --help                      Print this help.

EXAMPLE:
  cargo run -q -p flagged-lorentzian --bin active_component_transfer_scan -- \\
    --lambda 4,2 --mu 2 --alphabet 5"
}
