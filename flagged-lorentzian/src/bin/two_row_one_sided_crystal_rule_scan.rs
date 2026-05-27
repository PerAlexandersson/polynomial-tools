use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use flagged_lorentzian::{
    crystal_e, crystal_f, enumerate_tableaux, families::skew_shapes_of_size, RowFlaggedSkewShape,
    SkewGtPattern, SkewShape,
};

type Content = Vec<u32>;
type GtRows = Vec<Vec<u32>>;

#[derive(Debug, Clone)]
struct Args {
    lambda: Option<Vec<u32>>,
    mu: Vec<u32>,
    row_flags: Option<Vec<u32>>,
    alphabet: usize,
    lower_label: u32,
    max_skew_size: u32,
    max_outer_extra: u32,
    exactly_two_rows: bool,
    rule: Rule,
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
            exactly_two_rows: true,
            rule: Rule::LeftThenRight,
            tableau_limit: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Rule {
    LeftOnly,
    RightOnly,
    LeftThenRight,
    RightThenLeft,
    PreferLowerTarget,
    PreferUpperTarget,
    LeftThenCarrierMin,
    LeftThenCarrierMax,
    LeftThenShortestCarrierMin,
    LeftThenShortestCarrierMax,
}

#[derive(Debug, Clone)]
struct ShapeOutcome {
    fibers_checked: usize,
    negative_pairs: usize,
    failure: Option<RuleFailure>,
    limit_exceeded: bool,
}

#[derive(Debug, Clone)]
struct RuleFailure {
    beta: Content,
    z: GtRows,
    source: (usize, usize),
    message: String,
    left_target: Option<(usize, usize)>,
    right_target: Option<(usize, usize)>,
    chosen_target: Option<(usize, usize)>,
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
    let mut fibers_checked = 0usize;
    let mut negative_pairs = 0usize;

    for skew_size in 0..=args.max_skew_size {
        let max_outer_size = skew_size + args.max_outer_extra;
        for shape in skew_shapes_of_size(skew_size, max_outer_size) {
            if args.exactly_two_rows && shape.row_count() != 2 {
                continue;
            }
            if shape.row_count() > 2 {
                continue;
            }
            let flagged_shape = flagged_shape_for_args(shape, &args);
            let outcome = scan_one_shape(&flagged_shape, &args);
            if outcome.limit_exceeded {
                skipped_by_limit += 1;
                continue;
            }
            shapes_checked += 1;
            fibers_checked += outcome.fibers_checked;
            negative_pairs += outcome.negative_pairs;
            if let Some(failure) = outcome.failure {
                println!("FAIL");
                println!("outer={:?}", flagged_shape.shape().outer().parts());
                println!("inner={:?}", flagged_shape.shape().inner().parts());
                println!("flags={:?}", flagged_shape.row_flags());
                println!(
                    "beta={:?}, {}, source={:?}, left_target={:?}, right_target={:?}, chosen_target={:?}, z={}",
                    failure.beta,
                    failure.message,
                    failure.source,
                    failure.left_target,
                    failure.right_target,
                    failure.chosen_target,
                    summarize_gt(&failure.z)
                );
                println!(
                    "shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, fibers_checked={fibers_checked}, negative_pairs={negative_pairs}, elapsed={:.3}s",
                    started.elapsed().as_secs_f64()
                );
                std::process::exit(1);
            }
        }
        println!(
            "skew_size={skew_size}: shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, fibers_checked={fibers_checked}, negative_pairs={negative_pairs}, elapsed={:.3}s",
            started.elapsed().as_secs_f64()
        );
    }

    println!("PASS");
    println!(
        "shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, fibers_checked={fibers_checked}, negative_pairs={negative_pairs}, elapsed={:.3}s",
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
                "--at-most-two-rows" => args.exactly_two_rows = false,
                "--rule" => args.rule = Rule::parse(&take_value(&mut iter, &flag)?)?,
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

impl Rule {
    fn parse(input: &str) -> Result<Self, String> {
        match input {
            "left-only" => Ok(Self::LeftOnly),
            "right-only" => Ok(Self::RightOnly),
            "left-then-right" => Ok(Self::LeftThenRight),
            "right-then-left" => Ok(Self::RightThenLeft),
            "prefer-lower-target" => Ok(Self::PreferLowerTarget),
            "prefer-upper-target" => Ok(Self::PreferUpperTarget),
            "left-then-carrier-min" => Ok(Self::LeftThenCarrierMin),
            "left-then-carrier-max" => Ok(Self::LeftThenCarrierMax),
            "left-then-shortest-carrier-min" => Ok(Self::LeftThenShortestCarrierMin),
            "left-then-shortest-carrier-max" => Ok(Self::LeftThenShortestCarrierMax),
            other => Err(format!(
                "invalid --rule `{other}`; expected left-only, right-only, left-then-right, right-then-left, prefer-lower-target, prefer-upper-target, left-then-carrier-min, left-then-carrier-max, left-then-shortest-carrier-min, or left-then-shortest-carrier-max"
            )),
        }
    }
}

fn scan_one_shape(flagged_shape: &RowFlaggedSkewShape, args: &Args) -> ShapeOutcome {
    let tableaux = match enumerate_tableaux(flagged_shape, args.tableau_limit) {
        Ok(tableaux) => tableaux,
        Err(_) => {
            return ShapeOutcome {
                fibers_checked: 0,
                negative_pairs: 0,
                failure: None,
                limit_exceeded: true,
            }
        }
    };

    let shape = flagged_shape.shape();
    let gt_rows: Vec<_> = tableaux
        .iter()
        .map(|tableau| {
            SkewGtPattern::from_tableau(shape, &tableau.values, args.alphabet)
                .rows()
                .to_vec()
        })
        .collect();
    let gt_to_pos: BTreeMap<_, _> = gt_rows
        .iter()
        .enumerate()
        .map(|(pos, rows)| (rows.clone(), pos))
        .collect();
    let values_to_pos: BTreeMap<_, _> = tableaux
        .iter()
        .enumerate()
        .map(|(pos, tableau)| (tableau.values.clone(), pos))
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
    let mut fibers_checked = 0usize;
    let mut negative_pairs = 0usize;

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
        negative_pairs += ii_indices.len() * jj_indices.len();
        if let Some(mut failure) = scan_fiber(
            shape,
            ii_indices,
            mixed_indices,
            jj_indices,
            &gt_rows,
            &gt_to_pos,
            &values_to_pos,
            args,
        ) {
            failure.beta = beta;
            return ShapeOutcome {
                fibers_checked,
                negative_pairs,
                failure: Some(failure),
                limit_exceeded: false,
            };
        }
    }

    ShapeOutcome {
        fibers_checked,
        negative_pairs,
        failure: None,
        limit_exceeded: false,
    }
}

#[allow(clippy::too_many_arguments)]
fn scan_fiber(
    shape: &SkewShape,
    ii_indices: &[usize],
    mixed_indices: &[usize],
    jj_indices: &[usize],
    gt_rows: &[GtRows],
    gt_to_pos: &BTreeMap<GtRows, usize>,
    values_to_pos: &BTreeMap<Vec<u32>, usize>,
    args: &Args,
) -> Option<RuleFailure> {
    let mixed_pos_by_index: BTreeMap<_, _> = mixed_indices
        .iter()
        .copied()
        .enumerate()
        .map(|(pos, idx)| (idx, pos))
        .collect();
    let mut used_by_z = BTreeMap::<GtRows, BTreeSet<(usize, usize)>>::new();

    for (left_pos, &left_idx) in ii_indices.iter().enumerate() {
        for (right_pos, &right_idx) in jj_indices.iter().enumerate() {
            let z = add_gt_rows(&gt_rows[left_idx], &gt_rows[right_idx]);
            let left_target = left_crystal_target(
                shape,
                left_idx,
                &z,
                mixed_indices,
                gt_rows,
                gt_to_pos,
                values_to_pos,
                &mixed_pos_by_index,
                args.lower_label,
            );
            let right_target = right_crystal_target(
                shape,
                right_idx,
                &z,
                mixed_indices,
                gt_rows,
                gt_to_pos,
                values_to_pos,
                &mixed_pos_by_index,
                args.lower_label,
            );
            let carrier_targets = carrier_targets_for_source(
                &gt_rows[left_idx],
                &z,
                mixed_indices,
                gt_rows,
                gt_to_pos,
                &mixed_pos_by_index,
                args.lower_label as usize,
            );
            let shortest_carrier_targets = shortest_carrier_targets(&carrier_targets);
            let chosen = choose_target(left_target, right_target, &carrier_targets, args.rule);
            let chosen =
                match args.rule {
                    Rule::LeftThenShortestCarrierMin => left_target
                        .or_else(|| shortest_carrier_targets.first().map(|target| target.0)),
                    Rule::LeftThenShortestCarrierMax => left_target
                        .or_else(|| shortest_carrier_targets.last().map(|target| target.0)),
                    _ => chosen,
                };
            let Some(chosen) = chosen else {
                return Some(RuleFailure {
                    beta: Vec::new(),
                    z,
                    source: (left_pos, right_pos),
                    message: "rule has no target".to_string(),
                    left_target,
                    right_target,
                    chosen_target: None,
                });
            };
            if !used_by_z.entry(z.clone()).or_default().insert(chosen) {
                return Some(RuleFailure {
                    beta: Vec::new(),
                    z,
                    source: (left_pos, right_pos),
                    message: format!("{:?} is not injective on this fixed-Z fiber", args.rule),
                    left_target,
                    right_target,
                    chosen_target: Some(chosen),
                });
            }
        }
    }
    None
}

#[allow(clippy::too_many_arguments)]
fn left_crystal_target(
    shape: &SkewShape,
    left_idx: usize,
    z: &[Vec<u32>],
    mixed_indices: &[usize],
    gt_rows: &[GtRows],
    gt_to_pos: &BTreeMap<GtRows, usize>,
    values_to_pos: &BTreeMap<Vec<u32>, usize>,
    mixed_pos_by_index: &BTreeMap<usize, usize>,
    lower_label: u32,
) -> Option<(usize, usize)> {
    let values = values_to_pos
        .iter()
        .find_map(|(values, &idx)| (idx == left_idx).then_some(values))?;
    let image_values = crystal_f(shape, values, lower_label)?;
    let &left_image_idx = values_to_pos.get(&image_values)?;
    let &left_image_pos = mixed_pos_by_index.get(&left_image_idx)?;
    let complement = subtract_gt_rows(z, &gt_rows[left_image_idx])?;
    let &right_image_idx = gt_to_pos.get(&complement)?;
    let &right_image_pos = mixed_pos_by_index.get(&right_image_idx)?;
    mixed_indices.get(left_image_pos)?;
    mixed_indices.get(right_image_pos)?;
    Some((left_image_pos, right_image_pos))
}

#[allow(clippy::too_many_arguments)]
fn right_crystal_target(
    shape: &SkewShape,
    right_idx: usize,
    z: &[Vec<u32>],
    mixed_indices: &[usize],
    gt_rows: &[GtRows],
    gt_to_pos: &BTreeMap<GtRows, usize>,
    values_to_pos: &BTreeMap<Vec<u32>, usize>,
    mixed_pos_by_index: &BTreeMap<usize, usize>,
    lower_label: u32,
) -> Option<(usize, usize)> {
    let values = values_to_pos
        .iter()
        .find_map(|(values, &idx)| (idx == right_idx).then_some(values))?;
    let image_values = crystal_e(shape, values, lower_label)?;
    let &right_image_idx = values_to_pos.get(&image_values)?;
    let &right_image_pos = mixed_pos_by_index.get(&right_image_idx)?;
    let complement = subtract_gt_rows(z, &gt_rows[right_image_idx])?;
    let &left_image_idx = gt_to_pos.get(&complement)?;
    let &left_image_pos = mixed_pos_by_index.get(&left_image_idx)?;
    mixed_indices.get(left_image_pos)?;
    mixed_indices.get(right_image_pos)?;
    Some((left_image_pos, right_image_pos))
}

fn choose_target(
    left: Option<(usize, usize)>,
    right: Option<(usize, usize)>,
    carrier_targets: &[((usize, usize), usize)],
    rule: Rule,
) -> Option<(usize, usize)> {
    match rule {
        Rule::LeftOnly => left,
        Rule::RightOnly => right,
        Rule::LeftThenRight => left.or(right),
        Rule::RightThenLeft => right.or(left),
        Rule::PreferLowerTarget => match (left, right) {
            (Some(left), Some(right)) => Some(left.min(right)),
            (Some(left), None) => Some(left),
            (None, Some(right)) => Some(right),
            (None, None) => None,
        },
        Rule::PreferUpperTarget => match (left, right) {
            (Some(left), Some(right)) => Some(left.max(right)),
            (Some(left), None) => Some(left),
            (None, Some(right)) => Some(right),
            (None, None) => None,
        },
        Rule::LeftThenCarrierMin => left.or_else(|| carrier_targets.first().map(|target| target.0)),
        Rule::LeftThenCarrierMax => left.or_else(|| carrier_targets.last().map(|target| target.0)),
        Rule::LeftThenShortestCarrierMin => {
            left.or_else(|| carrier_targets.first().map(|target| target.0))
        }
        Rule::LeftThenShortestCarrierMax => {
            left.or_else(|| carrier_targets.last().map(|target| target.0))
        }
    }
}

fn carrier_targets_for_source(
    left: &[Vec<u32>],
    z: &[Vec<u32>],
    mixed_indices: &[usize],
    gt_rows: &[GtRows],
    gt_to_pos: &BTreeMap<GtRows, usize>,
    mixed_pos_by_index: &BTreeMap<usize, usize>,
    lower_label: usize,
) -> Vec<((usize, usize), usize)> {
    let mut out = Vec::new();
    for (left_pos, &left_idx) in mixed_indices.iter().enumerate() {
        let Some(carrier) = carrier_delta(left, &gt_rows[left_idx], lower_label) else {
            continue;
        };
        let Some(complement) = subtract_gt_rows(z, &gt_rows[left_idx]) else {
            continue;
        };
        let Some(&right_idx) = gt_to_pos.get(&complement) else {
            continue;
        };
        let Some(&right_pos) = mixed_pos_by_index.get(&right_idx) else {
            continue;
        };
        let start_level = carrier
            .iter()
            .position(|&delta| delta != [0, 0])
            .unwrap_or(0);
        out.push(((left_pos, right_pos), start_level));
    }
    out.sort_unstable();
    out.dedup();
    out
}

fn shortest_carrier_targets(
    carrier_targets: &[((usize, usize), usize)],
) -> Vec<((usize, usize), usize)> {
    let Some(max_start) = carrier_targets.iter().map(|target| target.1).max() else {
        return Vec::new();
    };
    carrier_targets
        .iter()
        .copied()
        .filter(|target| target.1 == max_start)
        .collect()
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
            "beta={:?}, {}, source={:?}, left_target={:?}, right_target={:?}, chosen_target={:?}, z={}",
            failure.beta,
            failure.message,
            failure.source,
            failure.left_target,
            failure.right_target,
            failure.chosen_target,
            summarize_gt(&failure.z)
        );
    } else {
        println!("PASS");
    }
    println!(
        "fibers_checked={}, negative_pairs={}, elapsed={elapsed_seconds:.3}s",
        outcome.fibers_checked, outcome.negative_pairs
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
    "Check source-local one-sided crystal rules inside fixed-Z two-row fibers.

USAGE:
  two_row_one_sided_crystal_rule_scan [OPTIONS]

OPTIONS:
  --lambda PARTS          Outer partition. If omitted, scan a family.
  --mu PARTS              Inner partition. Default: empty.
  --row-flags PARTS       Row upper flags. Omit for ordinary skew tableaux.
  --alphabet N            Alphabet size. Default: 5.
  --lower-label I         Active labels are I,I+1. Default: 4.
  --max-skew-size N       Maximum skew size for family scans. Default: 7.
  --max-outer-extra N     Family outer sizes go up to skew_size + this. Default: 8.
  --at-most-two-rows      Include one-row shapes in family scans.
  --rule RULE             left-only, right-only, left-then-right,
                          right-then-left, prefer-lower-target,
                          prefer-upper-target, left-then-carrier-min,
                          left-then-carrier-max,
                          left-then-shortest-carrier-min,
                          or left-then-shortest-carrier-max.
                          Default: left-then-right.
  --tableau-limit N       Skip shapes whose tableau enumeration exceeds N.
  --help                  Print this help."
}
