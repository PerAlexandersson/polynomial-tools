use std::collections::BTreeMap;

use flagged_lorentzian::{
    enumerate_tableaux, is_gt_array, RowFlaggedSkewShape, SkewGtPattern, SkewShape, TableauRecord,
};

type GtRows = Vec<Vec<u32>>;

#[derive(Debug, Clone)]
struct Args {
    lambda: Vec<u32>,
    mu: Vec<u32>,
    row_flags: Option<Vec<u32>>,
    alphabet: usize,
    lower_label: usize,
    beta: Vec<u32>,
    show_z_groups: bool,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            lambda: vec![6, 4],
            mu: vec![2],
            row_flags: None,
            alphabet: 5,
            lower_label: 4,
            beta: vec![0, 0, 1, 2, 3],
            show_z_groups: false,
        }
    }
}

#[derive(Debug, Clone)]
struct CarrierTarget {
    left_pos: usize,
    right_pos: usize,
    carrier: Vec<[i32; 2]>,
}

#[derive(Debug, Clone)]
enum TraceEvent {
    Try {
        depth: usize,
        source: usize,
        target: usize,
    },
    Conflict {
        depth: usize,
        target: usize,
        previous: usize,
    },
    Move {
        depth: usize,
        source: usize,
        old_target: Option<usize>,
        new_target: usize,
    },
}

fn main() {
    let args = match parse_args() {
        Ok(args) => args,
        Err(message) => {
            eprintln!("{message}");
            eprintln!();
            eprintln!("{}", help_text());
            std::process::exit(2);
        }
    };

    let shape = SkewShape::from_parts(args.lambda.clone(), args.mu.clone());
    let flagged = match &args.row_flags {
        Some(row_flags) => {
            RowFlaggedSkewShape::new(shape.clone(), row_flags.clone(), args.alphabet)
        }
        None => RowFlaggedSkewShape::ordinary(shape.clone(), args.alphabet),
    };
    let tableaux = enumerate_tableaux(&flagged, None).unwrap();
    let gt_rows: Vec<_> = tableaux
        .iter()
        .map(|tableau| {
            SkewGtPattern::from_tableau(&shape, &tableau.values, args.alphabet)
                .rows()
                .to_vec()
        })
        .collect();

    let lower = args.lower_label - 1;
    let upper = lower + 1;
    let ii = add_units(&args.beta, lower, 2);
    let mixed = add_units(&add_units(&args.beta, lower, 1), upper, 1);
    let jj = add_units(&args.beta, upper, 2);
    let ii_indices = indices_with_content(&tableaux, &ii);
    let mixed_indices = indices_with_content(&tableaux, &mixed);
    let jj_indices = indices_with_content(&tableaux, &jj);
    let mixed_count = mixed_indices.len();
    let mixed_gt_to_pos: BTreeMap<_, _> = mixed_indices
        .iter()
        .enumerate()
        .map(|(pos, &idx)| (gt_rows[idx].clone(), pos))
        .collect();

    println!("two-row carrier fiber dump");
    println!("outer={:?}", shape.outer().parts());
    println!("inner={:?}", shape.inner().parts());
    println!("flags={:?}", flagged.row_flags());
    println!("alphabet={}", args.alphabet);
    println!(
        "beta={:?}, active labels={},{}",
        args.beta,
        lower + 1,
        upper + 1
    );
    println!(
        "ii={}, mixed={}, jj={}, sources={}, ordered_targets={}",
        ii_indices.len(),
        mixed_indices.len(),
        jj_indices.len(),
        ii_indices.len() * jj_indices.len(),
        mixed_indices.len() * mixed_indices.len()
    );
    println!();

    print_tableaux("A", &shape, &tableaux, &ii_indices);
    print_tableaux("C", &shape, &tableaux, &jj_indices);
    print_tableaux("M", &shape, &tableaux, &mixed_indices);

    let source_pairs: Vec<_> = (0..ii_indices.len())
        .flat_map(|left_pos| (0..jj_indices.len()).map(move |right_pos| (left_pos, right_pos)))
        .collect();
    let mut edges_by_source: Vec<Vec<usize>> = Vec::with_capacity(source_pairs.len());
    let mut carriers_by_source: Vec<BTreeMap<usize, Vec<[i32; 2]>>> =
        Vec::with_capacity(source_pairs.len());
    for &(left_pos, right_pos) in &source_pairs {
        let left_idx = ii_indices[left_pos];
        let right_idx = jj_indices[right_pos];
        let targets = carrier_targets(
            &gt_rows[left_idx],
            &gt_rows[right_idx],
            &mixed_gt_to_pos,
            args.lower_label,
        );
        let mut carrier_by_target = BTreeMap::new();
        for target in targets {
            let target_index = target.left_pos * mixed_count + target.right_pos;
            carrier_by_target
                .entry(target_index)
                .or_insert(target.carrier);
        }
        edges_by_source.push(carrier_by_target.keys().copied().collect());
        carriers_by_source.push(carrier_by_target);
    }

    println!("carrier edges:");
    for (source, &(left_pos, right_pos)) in source_pairs.iter().enumerate() {
        let pieces: Vec<_> = edges_by_source[source]
            .iter()
            .map(|&target| {
                let carrier = &carriers_by_source[source][&target];
                format!(
                    "{} via {}",
                    target_pair_string(target, mixed_count),
                    format_carrier(carrier)
                )
            })
            .collect();
        println!(
            "  S{source}=(A{left_pos},C{right_pos}): {}",
            pieces.join("; ")
        );
    }
    println!();

    if args.show_z_groups {
        print_z_groups(
            &gt_rows,
            &ii_indices,
            &mixed_indices,
            &jj_indices,
            &source_pairs,
            mixed_count,
        );
        println!();
    }

    let (source_match, events) = canonical_matching_with_trace(&edges_by_source, mixed_count);
    println!("canonical augmenting events:");
    for event in events {
        match event {
            TraceEvent::Try {
                depth,
                source,
                target,
            } => println!(
                "  {}try S{source} -> {}",
                "  ".repeat(depth),
                target_pair_string(target, mixed_count)
            ),
            TraceEvent::Conflict {
                depth,
                target,
                previous,
            } => println!(
                "  {}occupied {} by S{previous}",
                "  ".repeat(depth),
                target_pair_string(target, mixed_count)
            ),
            TraceEvent::Move {
                depth,
                source,
                old_target,
                new_target,
            } => {
                let old = old_target
                    .map(|target| target_pair_string(target, mixed_count))
                    .unwrap_or_else(|| "none".to_string());
                println!(
                    "  {}move S{source}: {old} -> {}",
                    "  ".repeat(depth),
                    target_pair_string(new_target, mixed_count)
                );
            }
        }
    }
    println!();

    println!("final matching:");
    for (source, target) in source_match.iter().enumerate() {
        let (left_pos, right_pos) = source_pairs[source];
        match target {
            Some(target) => println!(
                "  S{source}=(A{left_pos},C{right_pos}) -> {}",
                target_pair_string(*target, mixed_count)
            ),
            None => println!("  S{source}=(A{left_pos},C{right_pos}) -> unmatched"),
        }
    }
}

fn canonical_matching_with_trace(
    edges_by_source: &[Vec<usize>],
    mixed_count: usize,
) -> (Vec<Option<usize>>, Vec<TraceEvent>) {
    let mut target_match = vec![None; mixed_count * mixed_count];
    let mut source_match = vec![None; edges_by_source.len()];
    let mut events = Vec::new();
    for source in 0..edges_by_source.len() {
        let mut seen = vec![false; target_match.len()];
        if !augment_with_trace(
            source,
            edges_by_source,
            &mut seen,
            &mut target_match,
            &mut source_match,
            &mut events,
            0,
        ) {
            events.push(TraceEvent::Move {
                depth: 0,
                source,
                old_target: source_match[source],
                new_target: usize::MAX,
            });
        }
    }
    (source_match, events)
}

fn augment_with_trace(
    source: usize,
    edges_by_source: &[Vec<usize>],
    seen: &mut [bool],
    target_match: &mut [Option<usize>],
    source_match: &mut [Option<usize>],
    events: &mut Vec<TraceEvent>,
    depth: usize,
) -> bool {
    for &target in &edges_by_source[source] {
        if seen[target] {
            continue;
        }
        seen[target] = true;
        events.push(TraceEvent::Try {
            depth,
            source,
            target,
        });
        let can_use = match target_match[target] {
            None => true,
            Some(previous) => {
                events.push(TraceEvent::Conflict {
                    depth,
                    target,
                    previous,
                });
                augment_with_trace(
                    previous,
                    edges_by_source,
                    seen,
                    target_match,
                    source_match,
                    events,
                    depth + 1,
                )
            }
        };
        if can_use {
            let old_target = source_match[source];
            target_match[target] = Some(source);
            source_match[source] = Some(target);
            events.push(TraceEvent::Move {
                depth,
                source,
                old_target,
                new_target: target,
            });
            return true;
        }
    }
    false
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
    out.sort_by_key(|target| (target.left_pos, target.right_pos));
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
        carrier[level] = [0, 1];
        carriers.push(carrier.clone());
        carrier[level] = [1, 0];
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

fn print_tableaux(label: &str, shape: &SkewShape, tableaux: &[TableauRecord], indices: &[usize]) {
    println!("{label} tableaux:");
    for (pos, &idx) in indices.iter().enumerate() {
        println!(
            "{label}{pos}: {}",
            format_tableau(shape, &tableaux[idx].values)
        );
    }
    println!();
}

fn format_carrier(carrier: &[[i32; 2]]) -> String {
    let parts: Vec<_> = carrier
        .iter()
        .enumerate()
        .filter(|(_, delta)| **delta != [0, 0])
        .map(|(level, delta)| format!("D{level}=({},{})", delta[0], delta[1]))
        .collect();
    if parts.is_empty() {
        "empty".to_string()
    } else {
        parts.join(", ")
    }
}

fn target_pair_string(target: usize, mixed_count: usize) -> String {
    if target == usize::MAX {
        return "failed".to_string();
    }
    format!("(M{},M{})", target / mixed_count, target % mixed_count)
}

fn print_z_groups(
    gt_rows: &[GtRows],
    ii_indices: &[usize],
    mixed_indices: &[usize],
    jj_indices: &[usize],
    source_pairs: &[(usize, usize)],
    mixed_count: usize,
) {
    let mut sources_by_z = BTreeMap::<GtRows, Vec<usize>>::new();
    for (source, &(left_pos, right_pos)) in source_pairs.iter().enumerate() {
        let z = add_gt_rows(
            &gt_rows[ii_indices[left_pos]],
            &gt_rows[jj_indices[right_pos]],
        );
        sources_by_z.entry(z).or_default().push(source);
    }

    let mut targets_by_z = BTreeMap::<GtRows, Vec<usize>>::new();
    for (left_pos, &left_idx) in mixed_indices.iter().enumerate() {
        for (right_pos, &right_idx) in mixed_indices.iter().enumerate() {
            let z = add_gt_rows(&gt_rows[left_idx], &gt_rows[right_idx]);
            targets_by_z
                .entry(z)
                .or_default()
                .push(left_pos * mixed_count + right_pos);
        }
    }

    println!("fixed pair-sum groups:");
    for (z, sources) in sources_by_z {
        let target_strings: Vec<_> = targets_by_z
            .get(&z)
            .into_iter()
            .flat_map(|targets| targets.iter())
            .map(|&target| target_pair_string(target, mixed_count))
            .collect();
        let source_strings: Vec<_> = sources
            .iter()
            .map(|&source| {
                let (left_pos, right_pos) = source_pairs[source];
                format!("S{source}=(A{left_pos},C{right_pos})")
            })
            .collect();
        println!(
            "  Z {}: sources [{}], targets [{}]",
            summarize_gt(&z),
            source_strings.join(", "),
            target_strings.join(", ")
        );
    }
}

fn add_gt_rows(left: &[Vec<u32>], right: &[Vec<u32>]) -> GtRows {
    left.iter()
        .zip(right)
        .map(|(left, right)| left.iter().zip(right).map(|(&a, &b)| a + b).collect())
        .collect()
}

fn summarize_gt(rows: &[Vec<u32>]) -> String {
    rows.iter()
        .enumerate()
        .filter(|(_, row)| row.iter().any(|&entry| entry != 0))
        .map(|(level, row)| format!("r{level}={row:?}"))
        .collect::<Vec<_>>()
        .join(";")
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

fn parse_args() -> Result<Args, String> {
    let mut args = Args::default();
    let mut iter = std::env::args().skip(1);
    while let Some(flag) = iter.next() {
        match flag.as_str() {
            "--help" | "-h" => {
                println!("{}", help_text());
                std::process::exit(0);
            }
            "--lambda" => args.lambda = parse_u32_vec(&take_value(&mut iter, &flag)?)?,
            "--mu" => args.mu = parse_u32_vec(&take_value(&mut iter, &flag)?)?,
            "--row-flags" => args.row_flags = Some(parse_u32_vec(&take_value(&mut iter, &flag)?)?),
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
            "--beta" => args.beta = parse_u32_vec(&take_value(&mut iter, &flag)?)?,
            "--show-z-groups" => args.show_z_groups = true,
            other => return Err(format!("unknown argument `{other}`")),
        }
    }
    if args.lower_label == 0 || args.lower_label >= args.alphabet {
        return Err("--lower-label must be between 1 and alphabet-1".to_string());
    }
    if args.beta.len() != args.alphabet {
        return Err("--beta must have one entry per alphabet letter".to_string());
    }
    Ok(args)
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
    "Dump the carrier graph and canonical augmenting matching for one two-row fiber.

USAGE:
  two_row_carrier_fiber_dump [OPTIONS]

OPTIONS:
  --lambda PARTS      Outer partition. Default: 6,4.
  --mu PARTS          Inner partition. Default: 2.
  --row-flags PARTS   Row upper flags. Omit for ordinary skew tableaux.
  --alphabet N        Alphabet size. Default: 5.
  --lower-label I     Active labels are I,I+1. Default: 4.
  --beta CONTENT      Base content. Default: 0,0,1,2,3.
  --show-z-groups     Also print source/target groups with the same pair-sum.
  --help              Print this help."
}
