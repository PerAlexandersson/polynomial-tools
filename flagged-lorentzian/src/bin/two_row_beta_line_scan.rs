use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use flagged_lorentzian::{
    enumerate_tableaux, families::skew_shapes_of_size, is_gt_array, RowFlaggedSkewShape,
    SkewGtPattern, SkewShape, TableauRecord,
};

type Content = Vec<u32>;
type GtRows = Vec<Vec<u32>>;
type FrameKey = Vec<(usize, Vec<u32>)>;

#[derive(Debug, Clone)]
struct Args {
    lambda: Option<Vec<u32>>,
    mu: Vec<u32>,
    row_flags: Option<Vec<u32>>,
    alphabet: usize,
    lower_label: u32,
    max_skew_size: u32,
    max_outer_extra: u32,
    tableau_limit: Option<usize>,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            lambda: None,
            mu: Vec::new(),
            row_flags: None,
            alphabet: 5,
            lower_label: 4,
            max_skew_size: 7,
            max_outer_extra: 8,
            tableau_limit: None,
        }
    }
}

#[derive(Debug, Clone)]
struct ShapeOutcome {
    lines_checked: usize,
    source_elements: usize,
    failure: Option<LineFailure>,
    limit_exceeded: bool,
}

#[derive(Debug, Clone)]
struct LineFailure {
    beta: Content,
    z: GtRows,
    ranks: BTreeMap<u32, usize>,
    message: String,
}

#[derive(Debug, Clone)]
struct Split {
    pos: usize,
    content: Content,
    rows: GtRows,
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
        let flagged_shape = flagged_shape_for_args(shape, &args);
        let outcome = scan_one_shape(&flagged_shape, &args);
        print_single_outcome(&flagged_shape, &outcome, started.elapsed().as_secs_f64());
        std::process::exit(if outcome.failure.is_some() { 1 } else { 0 });
    }

    let mut shapes_checked = 0usize;
    let mut skipped_by_limit = 0usize;
    let mut lines_checked = 0usize;
    let mut source_elements = 0usize;

    for skew_size in 0..=args.max_skew_size {
        let max_outer_size = skew_size + args.max_outer_extra;
        for shape in skew_shapes_of_size(skew_size, max_outer_size) {
            if shape.row_count() != 2 {
                continue;
            }
            let flagged_shape = flagged_shape_for_args(shape, &args);
            let outcome = scan_one_shape(&flagged_shape, &args);
            if outcome.limit_exceeded {
                skipped_by_limit += 1;
                continue;
            }
            shapes_checked += 1;
            lines_checked += outcome.lines_checked;
            source_elements += outcome.source_elements;
            if let Some(failure) = outcome.failure {
                println!("FAIL");
                println!("outer={:?}", flagged_shape.shape().outer().parts());
                println!("inner={:?}", flagged_shape.shape().inner().parts());
                println!("flags={:?}", flagged_shape.row_flags());
                println!(
                    "beta={:?}, {}, ranks={:?}, z={}",
                    failure.beta,
                    failure.message,
                    failure.ranks,
                    summarize_gt(&failure.z)
                );
                println!(
                    "shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, lines_checked={lines_checked}, source_elements={source_elements}, elapsed={:.3}s",
                    started.elapsed().as_secs_f64()
                );
                std::process::exit(1);
            }
        }
        println!(
            "skew_size={skew_size}: shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, lines_checked={lines_checked}, source_elements={source_elements}, elapsed={:.3}s",
            started.elapsed().as_secs_f64()
        );
    }

    println!("PASS");
    println!(
        "shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, lines_checked={lines_checked}, source_elements={source_elements}, elapsed={:.3}s",
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
                "--row-flags" => {
                    args.row_flags = Some(parse_u32_vec(&take_value(&mut iter, &flag)?)?)
                }
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
                "--tableau-limit" => {
                    args.tableau_limit = Some(
                        take_value(&mut iter, &flag)?
                            .parse()
                            .map_err(|err| format!("invalid --tableau-limit: {err}"))?,
                    )
                }
                other => return Err(format!("unknown argument `{other}`")),
            }
        }
        if args.lower_label == 0 || args.lower_label as usize >= args.alphabet {
            return Err("--lower-label must be between 1 and alphabet-1".to_string());
        }
        Ok(args)
    }
}

fn scan_one_shape(flagged_shape: &RowFlaggedSkewShape, args: &Args) -> ShapeOutcome {
    let tableaux = match enumerate_tableaux(flagged_shape, args.tableau_limit) {
        Ok(tableaux) => tableaux,
        Err(_) => {
            return ShapeOutcome {
                lines_checked: 0,
                source_elements: 0,
                failure: None,
                limit_exceeded: true,
            }
        }
    };
    let gt_rows: Vec<_> = tableaux
        .iter()
        .map(|tableau| {
            SkewGtPattern::from_tableau(flagged_shape.shape(), &tableau.values, args.alphabet)
                .rows()
                .to_vec()
        })
        .collect();
    let gt_to_pos: BTreeMap<_, _> = gt_rows
        .iter()
        .enumerate()
        .map(|(pos, rows)| (rows.clone(), pos))
        .collect();

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
    let mut lines_checked = 0usize;
    let mut source_elements = 0usize;

    for content in by_content.keys() {
        let Some(beta) = subtract_units(content, lower, 2) else {
            continue;
        };
        if !seen_beta.insert(beta.clone()) {
            continue;
        }
        let ii = add_units(&beta, lower, 2);
        let mixed = add_units(&add_units(&beta, lower, 1), upper, 1);
        let jj = add_units(&beta, upper, 2);
        let (Some(ii_indices), Some(mixed_indices), Some(jj_indices)) = (
            by_content.get(&ii),
            by_content.get(&mixed),
            by_content.get(&jj),
        ) else {
            continue;
        };

        for &left_idx in ii_indices {
            for &right_idx in jj_indices {
                let z = add_gt_rows(&gt_rows[left_idx], &gt_rows[right_idx]);
                let splits = split_poset(&z, &gt_rows, &gt_to_pos, &tableaux);
                let line: Vec<_> = splits
                    .iter()
                    .filter(|split| outside_content_matches(&split.content, &beta, lower, upper))
                    .collect();
                if line.is_empty() {
                    return ShapeOutcome {
                        lines_checked,
                        source_elements,
                        failure: Some(LineFailure {
                            beta,
                            z,
                            ranks: BTreeMap::new(),
                            message: "empty fixed-outside-content line".to_string(),
                        }),
                        limit_exceeded: false,
                    };
                }
                let total_i = z[args.lower_label as usize][0] + z[args.lower_label as usize][1]
                    - z[args.lower_label as usize - 1][0]
                    - z[args.lower_label as usize - 1][1];
                let mut ranks = BTreeMap::<u32, usize>::new();
                for split in &line {
                    *ranks.entry(split.content[lower]).or_default() += 1;
                }
                lines_checked += 1;
                source_elements += ranks.get(&(beta[lower] + 2)).copied().unwrap_or(0);

                if let Some(message) = rank_failure_message(&ranks, total_i, beta[lower] + 1) {
                    return ShapeOutcome {
                        lines_checked,
                        source_elements,
                        failure: Some(LineFailure {
                            beta,
                            z,
                            ranks,
                            message,
                        }),
                        limit_exceeded: false,
                    };
                }
                if ranks.get(&(beta[lower] + 2)).copied().unwrap_or(0)
                    > ranks.get(&(beta[lower] + 1)).copied().unwrap_or(0)
                {
                    return ShapeOutcome {
                        lines_checked,
                        source_elements,
                        failure: Some(LineFailure {
                            beta,
                            z,
                            ranks,
                            message: "central injection count fails".to_string(),
                        }),
                        limit_exceeded: false,
                    };
                }
                if let Some(message) =
                    adjacent_matching_failure(&line, lower, args.lower_label as usize)
                {
                    return ShapeOutcome {
                        lines_checked,
                        source_elements,
                        failure: Some(LineFailure {
                            beta,
                            z,
                            ranks,
                            message,
                        }),
                        limit_exceeded: false,
                    };
                }
                if let Some(message) =
                    active_rectangle_fiber_failure(&line, args.lower_label as usize, lower)
                {
                    return ShapeOutcome {
                        lines_checked,
                        source_elements,
                        failure: Some(LineFailure {
                            beta,
                            z,
                            ranks,
                            message,
                        }),
                        limit_exceeded: false,
                    };
                }
                if let Some(message) = active_rectangle_scd_rule_failure(
                    &line,
                    args.lower_label as usize,
                    lower,
                    beta[lower] + 2,
                ) {
                    return ShapeOutcome {
                        lines_checked,
                        source_elements,
                        failure: Some(LineFailure {
                            beta,
                            z,
                            ranks,
                            message,
                        }),
                        limit_exceeded: false,
                    };
                }
                let target_pos: BTreeSet<_> = mixed_indices.iter().copied().collect();
                let has_target = splits
                    .iter()
                    .filter(|split| split.content.as_slice() == ii.as_slice())
                    .all(|source| {
                        splits.iter().any(|target| {
                            target_pos.contains(&target.pos)
                                && carrier_delta(
                                    &source.rows,
                                    &target.rows,
                                    args.lower_label as usize,
                                )
                                .is_some()
                        })
                    });
                if !has_target {
                    return ShapeOutcome {
                        lines_checked,
                        source_elements,
                        failure: Some(LineFailure {
                            beta,
                            z,
                            ranks,
                            message: "some source has no carrier target on beta-line".to_string(),
                        }),
                        limit_exceeded: false,
                    };
                }
            }
        }
    }

    ShapeOutcome {
        lines_checked,
        source_elements,
        failure: None,
        limit_exceeded: false,
    }
}

fn rank_failure_message(
    ranks: &BTreeMap<u32, usize>,
    total_i: u32,
    central_rank: u32,
) -> Option<String> {
    for (&rank, &count) in ranks {
        let dual_rank = total_i.checked_sub(rank)?;
        let dual_count = ranks.get(&dual_rank).copied().unwrap_or(0);
        if count != dual_count {
            return Some(format!(
                "rank symmetry fails at rank {rank}: count={count}, dual rank {dual_rank} count={dual_count}, total_i={total_i}"
            ));
        }
    }
    let mut previous_rank = None;
    let mut previous_count = 0usize;
    for (&rank, &count) in ranks.range(..=central_rank) {
        if let Some(previous_rank) = previous_rank {
            if rank != previous_rank + 1 {
                return Some(format!(
                    "rank support has a gap before central rank: previous={previous_rank}, current={rank}"
                ));
            }
        }
        if count < previous_count {
            return Some(format!(
                "rank sequence decreases before center: previous_count={previous_count}, rank={rank}, count={count}"
            ));
        }
        previous_rank = Some(rank);
        previous_count = count;
    }
    None
}

fn active_rectangle_fiber_failure(
    line: &[&Split],
    lower_label: usize,
    active_index: usize,
) -> Option<String> {
    let mut by_frame = BTreeMap::<FrameKey, BTreeSet<(u32, u32, u32)>>::new();
    for split in line {
        let active_row = &split.rows[lower_label];
        by_frame
            .entry(frame_key(&split.rows, lower_label))
            .or_default()
            .insert((active_row[0], active_row[1], split.content[active_index]));
    }

    for (frame, active_rows) in by_frame {
        let top_min = active_rows.iter().map(|row| row.0).min().unwrap();
        let top_max = active_rows.iter().map(|row| row.0).max().unwrap();
        let bottom_min = active_rows.iter().map(|row| row.1).min().unwrap();
        let bottom_max = active_rows.iter().map(|row| row.1).max().unwrap();
        for top in top_min..=top_max {
            for bottom in bottom_min..=bottom_max {
                let lower_row = frame.iter().find(|(level, _)| *level == lower_label - 1)?;
                let rank = top + bottom - lower_row.1[0] - lower_row.1[1];
                if !active_rows.contains(&(top, bottom, rank)) {
                    return Some(format!(
                        "active rows over a fixed frame are not a rectangle: missing ({top},{bottom}) over frame {}",
                        format_frame(&frame)
                    ));
                }
            }
        }
    }
    None
}

fn active_rectangle_scd_rule_failure(
    line: &[&Split],
    lower_label: usize,
    active_index: usize,
    source_rank: u32,
) -> Option<String> {
    let sources: Vec<_> = line
        .iter()
        .enumerate()
        .filter(|(_, split)| split.content[active_index] == source_rank)
        .map(|(line_pos, _)| line_pos)
        .collect();
    let mut used_targets = BTreeSet::new();
    let mut unassigned: BTreeSet<_> = sources.iter().copied().collect();

    let mut boundary_sources: Vec<_> = sources
        .iter()
        .copied()
        .filter(|&source_pos| {
            local_scd_predecessor(line, source_pos, lower_label, active_index).is_none()
        })
        .collect();
    boundary_sources.sort_by_key(|&source_pos| {
        (
            carrier_targets_for(line, source_pos, lower_label, active_index).len(),
            source_pos,
        )
    });
    for source_pos in boundary_sources {
        let target_pos = match choose_recursive_fallback(
            line,
            source_pos,
            lower_label,
            active_index,
            &used_targets,
        ) {
            Ok(target_pos) => target_pos,
            Err(message) => return Some(message),
        };
        if let Err(message) = assign_rule_target(
            line,
            source_pos,
            target_pos,
            lower_label,
            &mut used_targets,
            &mut unassigned,
        ) {
            return Some(message);
        }
    }

    for &source_pos in &sources {
        if !unassigned.contains(&source_pos) {
            continue;
        }
        let carrier_targets = carrier_targets_for(line, source_pos, lower_label, active_index);
        if carrier_targets.len() == 1 {
            if let Err(message) = assign_rule_target(
                line,
                source_pos,
                carrier_targets[0].0,
                lower_label,
                &mut used_targets,
                &mut unassigned,
            ) {
                return Some(message);
            }
        }
    }

    for &source_pos in &sources {
        if !unassigned.contains(&source_pos) {
            continue;
        }
        if let Some(target_pos) = local_scd_predecessor(line, source_pos, lower_label, active_index)
        {
            if !used_targets.contains(&target_pos) {
                if let Err(message) = assign_rule_target(
                    line,
                    source_pos,
                    target_pos,
                    lower_label,
                    &mut used_targets,
                    &mut unassigned,
                ) {
                    return Some(message);
                }
            }
        }
    }

    while !unassigned.is_empty() {
        let before = unassigned.len();
        let current: Vec<_> = unassigned.iter().copied().collect();
        for source_pos in current {
            if !unassigned.contains(&source_pos) {
                continue;
            }
            let available_targets: Vec<_> =
                carrier_targets_for(line, source_pos, lower_label, active_index)
                    .into_iter()
                    .filter(|(target_pos, _)| !used_targets.contains(target_pos))
                    .collect();
            if let Ok(target_pos) = choose_from_available_recursive_fallback(
                line,
                source_pos,
                lower_label,
                active_index,
                available_targets,
            ) {
                if let Err(message) = assign_rule_target(
                    line,
                    source_pos,
                    target_pos,
                    lower_label,
                    &mut used_targets,
                    &mut unassigned,
                ) {
                    return Some(message);
                }
            }
        }
        if unassigned.len() == before {
            let source_pos = *unassigned.iter().next().unwrap();
            return Some(format_rule_source(
                "recursive rule could not force a remaining source",
                line,
                source_pos,
                lower_label,
            ));
        }
    }
    None
}

fn carrier_targets_for(
    line: &[&Split],
    source_pos: usize,
    lower_label: usize,
    active_index: usize,
) -> Vec<(usize, usize)> {
    let source = line[source_pos];
    line.iter()
        .enumerate()
        .filter(|(_, target)| target.content[active_index] + 1 == source.content[active_index])
        .filter_map(|(target_pos, target)| {
            let carrier = carrier_delta(&source.rows, &target.rows, lower_label)?;
            let start = carrier.iter().position(|&row| row != [0, 0])?;
            Some((target_pos, start))
        })
        .collect()
}

fn choose_recursive_fallback(
    line: &[&Split],
    source_pos: usize,
    lower_label: usize,
    active_index: usize,
    used_targets: &BTreeSet<usize>,
) -> Result<usize, String> {
    let available_targets: Vec<_> =
        carrier_targets_for(line, source_pos, lower_label, active_index)
            .into_iter()
            .filter(|(target_pos, _)| !used_targets.contains(target_pos))
            .collect();
    choose_from_available_recursive_fallback(
        line,
        source_pos,
        lower_label,
        active_index,
        available_targets,
    )
}

fn choose_from_available_recursive_fallback(
    line: &[&Split],
    source_pos: usize,
    lower_label: usize,
    active_index: usize,
    available_targets: Vec<(usize, usize)>,
) -> Result<usize, String> {
    let Some(max_start) = available_targets.iter().map(|&(_, start)| start).max() else {
        return Err(format_rule_source(
            "recursive rule has no available carrier target",
            line,
            source_pos,
            lower_label,
        ));
    };
    let shortest_targets: Vec<_> = available_targets
        .into_iter()
        .filter(|&(_, start)| start == max_start)
        .map(|(target_pos, _)| target_pos)
        .collect();
    if shortest_targets.len() == 1 {
        return Ok(shortest_targets[0]);
    }

    let unclaimed_by_local: Vec<_> = shortest_targets
        .iter()
        .copied()
        .filter(|&target_pos| {
            local_scd_successor(line, target_pos, lower_label, active_index).is_none()
        })
        .collect();
    if unclaimed_by_local.len() == 1 {
        Ok(unclaimed_by_local[0])
    } else {
        Err(format!(
            "{}; shortest_targets={}, unclaimed_by_local={}",
            format_rule_source(
                "recursive fallback is not forced",
                line,
                source_pos,
                lower_label
            ),
            shortest_targets.len(),
            unclaimed_by_local.len()
        ))
    }
}

fn assign_rule_target(
    line: &[&Split],
    source_pos: usize,
    target_pos: usize,
    lower_label: usize,
    used_targets: &mut BTreeSet<usize>,
    unassigned: &mut BTreeSet<usize>,
) -> Result<(), String> {
    let source = line[source_pos];
    let target = line[target_pos];
    if carrier_delta(&source.rows, &target.rows, lower_label).is_none() {
        return Err(format_rule_source(
            "recursive rule chose a non-carrier target",
            line,
            source_pos,
            lower_label,
        ));
    }
    if !used_targets.insert(target_pos) {
        return Err(format!(
            "{}; repeated target active row ({},{})",
            format_rule_source(
                "recursive rule target collision",
                line,
                source_pos,
                lower_label
            ),
            target.rows[lower_label][0],
            target.rows[lower_label][1]
        ));
    }
    unassigned.remove(&source_pos);
    Ok(())
}

fn local_scd_predecessor(
    line: &[&Split],
    source_pos: usize,
    lower_label: usize,
    active_index: usize,
) -> Option<usize> {
    let source = line[source_pos];
    let positions: Vec<_> = line
        .iter()
        .enumerate()
        .filter(|(_, split)| same_frame(&split.rows, &source.rows, lower_label))
        .map(|(line_pos, _)| line_pos)
        .collect();
    let top_min = positions
        .iter()
        .map(|&line_pos| line[line_pos].rows[lower_label][0])
        .min()?;
    let top_max = positions
        .iter()
        .map(|&line_pos| line[line_pos].rows[lower_label][0])
        .max()?;
    let bottom_min = positions
        .iter()
        .map(|&line_pos| line[line_pos].rows[lower_label][1])
        .min()?;
    let bottom_max = positions
        .iter()
        .map(|&line_pos| line[line_pos].rows[lower_label][1])
        .max()?;
    let top_span = top_max - top_min;
    let bottom_span = bottom_max - bottom_min;
    let source_row = &source.rows[lower_label];
    let u = source_row[0] - top_min;
    let v = source_row[1] - bottom_min;
    let (pred_u, pred_v) = product_chain_predecessor(u, v, top_span, bottom_span)?;
    let target_row = (top_min + pred_u, bottom_min + pred_v);
    positions.into_iter().find(|&target_pos| {
        let target = line[target_pos];
        target.content[active_index] + 1 == source.content[active_index]
            && target.rows[lower_label][0] == target_row.0
            && target.rows[lower_label][1] == target_row.1
    })
}

fn format_rule_source(
    message: &str,
    line: &[&Split],
    source_pos: usize,
    lower_label: usize,
) -> String {
    let source = line[source_pos];
    format!(
        "{message}: source active row ({},{}), frame {}",
        source.rows[lower_label][0],
        source.rows[lower_label][1],
        format_frame(&frame_key(&source.rows, lower_label))
    )
}

fn product_chain_predecessor(
    u: u32,
    v: u32,
    top_span: u32,
    bottom_span: u32,
) -> Option<(u32, u32)> {
    if top_span <= bottom_span {
        if u + v <= top_span {
            u.checked_sub(1).map(|previous_u| (previous_u, v))
        } else {
            v.checked_sub(1).map(|previous_v| (u, previous_v))
        }
    } else if u + v <= bottom_span {
        v.checked_sub(1).map(|previous_v| (u, previous_v))
    } else {
        u.checked_sub(1).map(|previous_u| (previous_u, v))
    }
}

fn local_scd_successor(
    line: &[&Split],
    target_pos: usize,
    lower_label: usize,
    active_index: usize,
) -> Option<usize> {
    let target = line[target_pos];
    let positions: Vec<_> = line
        .iter()
        .enumerate()
        .filter(|(_, split)| same_frame(&split.rows, &target.rows, lower_label))
        .map(|(line_pos, _)| line_pos)
        .collect();
    let top_min = positions
        .iter()
        .map(|&line_pos| line[line_pos].rows[lower_label][0])
        .min()?;
    let top_max = positions
        .iter()
        .map(|&line_pos| line[line_pos].rows[lower_label][0])
        .max()?;
    let bottom_min = positions
        .iter()
        .map(|&line_pos| line[line_pos].rows[lower_label][1])
        .min()?;
    let bottom_max = positions
        .iter()
        .map(|&line_pos| line[line_pos].rows[lower_label][1])
        .max()?;
    let top_span = top_max - top_min;
    let bottom_span = bottom_max - bottom_min;
    let target_row = (
        target.rows[lower_label][0] - top_min,
        target.rows[lower_label][1] - bottom_min,
    );

    positions
        .into_iter()
        .filter(|&source_pos| {
            line[source_pos].content[active_index] == target.content[active_index] + 1
        })
        .find(|&source_pos| {
            let source_row = &line[source_pos].rows[lower_label];
            let u = source_row[0] - top_min;
            let v = source_row[1] - bottom_min;
            product_chain_predecessor(u, v, top_span, bottom_span) == Some(target_row)
        })
}

fn frame_key(rows: &[Vec<u32>], lower_label: usize) -> FrameKey {
    rows.iter()
        .enumerate()
        .filter(|(level, _)| *level != lower_label)
        .map(|(level, row)| (level, row.clone()))
        .collect()
}

fn same_frame(left: &[Vec<u32>], right: &[Vec<u32>], lower_label: usize) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .enumerate()
            .all(|(level, (left_row, right_row))| level == lower_label || left_row == right_row)
}

fn format_frame(frame: &[(usize, Vec<u32>)]) -> String {
    frame
        .iter()
        .map(|(level, row)| format!("r{level}={row:?}"))
        .collect::<Vec<_>>()
        .join(";")
}

fn adjacent_matching_failure(
    line: &[&Split],
    active_index: usize,
    lower_label: usize,
) -> Option<String> {
    let mut positions_by_rank = BTreeMap::<u32, Vec<usize>>::new();
    for (pos, split) in line.iter().enumerate() {
        positions_by_rank
            .entry(split.content[active_index])
            .or_default()
            .push(pos);
    }

    for (&rank, upper_positions) in &positions_by_rank {
        if rank == 0 {
            continue;
        }
        let Some(lower_positions) = positions_by_rank.get(&(rank - 1)) else {
            continue;
        };
        if upper_positions.len() > lower_positions.len() {
            continue;
        }

        let lower_position_to_target: BTreeMap<_, _> = lower_positions
            .iter()
            .copied()
            .enumerate()
            .map(|(target_pos, line_pos)| (line_pos, target_pos))
            .collect();
        let mut edges = vec![Vec::new(); upper_positions.len()];
        for (source_pos, &upper_line_pos) in upper_positions.iter().enumerate() {
            let source = line[upper_line_pos];
            for &lower_line_pos in lower_positions {
                let target = line[lower_line_pos];
                if carrier_delta(&source.rows, &target.rows, lower_label).is_some() {
                    edges[source_pos].push(lower_position_to_target[&lower_line_pos]);
                }
            }
        }
        let matched = max_bipartite_matching(&edges, lower_positions.len());
        if matched < upper_positions.len() {
            return Some(format!(
                "adjacent carrier matching fails from rank {rank} to {}: matched={matched}, sources={}, targets={}",
                rank - 1,
                upper_positions.len(),
                lower_positions.len()
            ));
        }
    }
    None
}

fn max_bipartite_matching(edges: &[Vec<usize>], right_count: usize) -> usize {
    let mut right_match = vec![None; right_count];
    let mut matched = 0usize;
    for left in 0..edges.len() {
        let mut seen = vec![false; right_count];
        if augment(left, edges, &mut seen, &mut right_match) {
            matched += 1;
        }
    }
    matched
}

fn augment(
    left: usize,
    edges: &[Vec<usize>],
    seen: &mut [bool],
    right_match: &mut [Option<usize>],
) -> bool {
    for &right in &edges[left] {
        if seen[right] {
            continue;
        }
        seen[right] = true;
        if right_match[right].is_none_or(|previous| augment(previous, edges, seen, right_match)) {
            right_match[right] = Some(left);
            return true;
        }
    }
    false
}

fn split_poset(
    z: &[Vec<u32>],
    gt_rows: &[GtRows],
    gt_to_pos: &BTreeMap<GtRows, usize>,
    tableaux: &[TableauRecord],
) -> Vec<Split> {
    let mut splits = Vec::new();
    for (pos, rows) in gt_rows.iter().enumerate() {
        let Some(complement) = subtract_gt_rows(z, rows) else {
            continue;
        };
        if !is_gt_array(&complement) {
            continue;
        }
        if !gt_to_pos.contains_key(&complement) {
            continue;
        }
        splits.push(Split {
            pos,
            content: tableaux[pos].content.clone(),
            rows: rows.clone(),
        });
    }
    splits
}

fn add_gt_rows(left: &[Vec<u32>], right: &[Vec<u32>]) -> GtRows {
    left.iter()
        .zip(right)
        .map(|(left, right)| left.iter().zip(right).map(|(&a, &b)| a + b).collect())
        .collect()
}

fn subtract_gt_rows(left: &[Vec<u32>], right: &[Vec<u32>]) -> Option<GtRows> {
    left.iter()
        .zip(right)
        .map(|(left, right)| {
            left.iter()
                .zip(right)
                .map(|(&a, &b)| a.checked_sub(b))
                .collect::<Option<Vec<_>>>()
        })
        .collect()
}

fn outside_content_matches(content: &[u32], beta: &[u32], lower: usize, upper: usize) -> bool {
    content
        .iter()
        .enumerate()
        .all(|(idx, &count)| idx == lower || idx == upper || beta.get(idx).copied() == Some(count))
}

fn carrier_delta(
    source: &[Vec<u32>],
    target: &[Vec<u32>],
    lower_label: usize,
) -> Option<Vec<[i32; 2]>> {
    if source.len() != target.len() || source.len() <= lower_label {
        return None;
    }
    let d: Vec<[i32; 2]> = source
        .iter()
        .zip(target)
        .map(|(source_row, target_row)| {
            [
                source_row[0] as i32 - target_row[0] as i32,
                source_row[1] as i32 - target_row[1] as i32,
            ]
        })
        .collect();
    if d.iter().skip(lower_label + 1).any(|&row| row != [0, 0]) {
        return None;
    }
    if d[lower_label] != [1, 0] && d[lower_label] != [0, 1] {
        return None;
    }
    let start = d.iter().position(|&row| row != [0, 0])?;
    if start > lower_label {
        return None;
    }
    if d.iter().take(start).any(|&row| row != [0, 0]) {
        return None;
    }
    if d.iter()
        .take(lower_label)
        .skip(start)
        .any(|&row| row != [1, -1] && row != [-1, 1])
    {
        return None;
    }
    Some(d)
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

fn summarize_gt(rows: &[Vec<u32>]) -> String {
    rows.iter()
        .enumerate()
        .filter(|(_, row)| row.iter().any(|&entry| entry != 0))
        .map(|(level, row)| format!("r{level}={row:?}"))
        .collect::<Vec<_>>()
        .join(";")
}

fn flagged_shape_for_args(shape: SkewShape, args: &Args) -> RowFlaggedSkewShape {
    match &args.row_flags {
        Some(row_flags) => RowFlaggedSkewShape::new(shape, row_flags.clone(), args.alphabet),
        None => RowFlaggedSkewShape::ordinary(shape, args.alphabet),
    }
}

fn print_single_outcome(
    flagged_shape: &RowFlaggedSkewShape,
    outcome: &ShapeOutcome,
    elapsed_seconds: f64,
) {
    if outcome.limit_exceeded {
        println!("SKIPPED: tableau limit exceeded");
        return;
    }
    println!("outer={:?}", flagged_shape.shape().outer().parts());
    println!("inner={:?}", flagged_shape.shape().inner().parts());
    println!("flags={:?}", flagged_shape.row_flags());
    if let Some(failure) = &outcome.failure {
        println!("FAIL");
        println!(
            "beta={:?}, {}, ranks={:?}, z={}",
            failure.beta,
            failure.message,
            failure.ranks,
            summarize_gt(&failure.z)
        );
    } else {
        println!("PASS");
    }
    println!(
        "lines_checked={}, source_elements={}, elapsed={elapsed_seconds:.3}s",
        outcome.lines_checked, outcome.source_elements
    );
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
    "Scan fixed-Z, fixed-outside-content beta-lines in two-row carrier fibers.

USAGE:
  two_row_beta_line_scan [OPTIONS]

OPTIONS:
  --lambda PARTS          Outer partition. If omitted, scan a family.
  --mu PARTS              Inner partition. Default: empty.
  --row-flags PARTS       Row upper flags. Omit for ordinary skew tableaux.
  --alphabet N            Alphabet size. Default: 5.
  --lower-label I         Active labels are I,I+1. Default: 4.
  --max-skew-size N       Maximum skew size for family scans. Default: 7.
  --max-outer-extra N     Family outer sizes go up to skew_size + this. Default: 8.
  --tableau-limit N       Skip shapes whose tableau enumeration exceeds N.
  --help                  Print this help."
}
