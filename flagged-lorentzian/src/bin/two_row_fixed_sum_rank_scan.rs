use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use flagged_lorentzian::{
    enumerate_tableaux, families::skew_shapes_of_size, RowFlaggedSkewShape, SkewGtPattern,
    SkewShape,
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
    require_carrier: bool,
    selection: SelectionMode,
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
            require_carrier: true,
            selection: SelectionMode::SameRank,
            tableau_limit: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SelectionMode {
    SameRank,
    Nearest,
    Ceiling,
    Floor,
    FirstEdge,
    LastEdge,
    ClosestSource,
    ClosestScaledSource,
}

#[derive(Debug, Clone)]
struct ShapeOutcome {
    fibers_checked: usize,
    negative_pairs: usize,
    failure: Option<RankFailure>,
    limit_exceeded: bool,
}

#[derive(Debug, Clone)]
struct RankFailure {
    beta: Content,
    message: String,
    z: GtRows,
    source_rank: usize,
    source: (usize, usize),
    target: Option<(usize, usize)>,
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
                    "beta={:?}, {}, source_rank={}, source={:?}, target={:?}, z={}",
                    failure.beta,
                    failure.message,
                    failure.source_rank,
                    failure.source,
                    failure.target,
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
                "--allow-noncarrier" => args.require_carrier = false,
                "--selection" => {
                    args.selection = SelectionMode::parse(&take_value(&mut iter, &flag)?)?
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
        if args.lower_label as usize >= args.alphabet {
            return Err("--lower-label must be smaller than --alphabet".to_string());
        }
        Ok(args)
    }
}

impl SelectionMode {
    fn parse(input: &str) -> Result<Self, String> {
        match input {
            "same-rank" => Ok(Self::SameRank),
            "nearest" => Ok(Self::Nearest),
            "ceiling" => Ok(Self::Ceiling),
            "floor" => Ok(Self::Floor),
            "first-edge" => Ok(Self::FirstEdge),
            "last-edge" => Ok(Self::LastEdge),
            "closest-source" => Ok(Self::ClosestSource),
            "closest-scaled-source" => Ok(Self::ClosestScaledSource),
            other => Err(format!(
                "invalid --selection `{other}`; expected same-rank, nearest, ceiling, floor, first-edge, last-edge, closest-source, or closest-scaled-source"
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

    let gt_rows: Vec<_> = tableaux
        .iter()
        .map(|tableau| {
            SkewGtPattern::from_tableau(flagged_shape.shape(), &tableau.values, args.alphabet)
                .rows()
                .to_vec()
        })
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
        let outcome = scan_fiber(
            ii_indices,
            mixed_indices,
            jj_indices,
            &gt_rows,
            args.lower_label as usize,
            args.require_carrier,
            args.selection,
        );
        negative_pairs += outcome.negative_pairs;
        if let Some(mut failure) = outcome.failure {
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

#[derive(Debug, Clone)]
struct FiberOutcome {
    negative_pairs: usize,
    failure: Option<RankFailure>,
}

fn scan_fiber(
    ii_indices: &[usize],
    mixed_indices: &[usize],
    jj_indices: &[usize],
    gt_rows: &[GtRows],
    lower_label: usize,
    require_carrier: bool,
    selection: SelectionMode,
) -> FiberOutcome {
    let mut sources_by_z = BTreeMap::<GtRows, Vec<(usize, usize)>>::new();
    for (left_pos, &left_idx) in ii_indices.iter().enumerate() {
        for (right_pos, &right_idx) in jj_indices.iter().enumerate() {
            let z = add_gt_rows(&gt_rows[left_idx], &gt_rows[right_idx]);
            sources_by_z
                .entry(z)
                .or_default()
                .push((left_pos, right_pos));
        }
    }

    let mut targets_by_z = BTreeMap::<GtRows, Vec<(usize, usize)>>::new();
    for (left_pos, &left_idx) in mixed_indices.iter().enumerate() {
        for (right_pos, &right_idx) in mixed_indices.iter().enumerate() {
            let z = add_gt_rows(&gt_rows[left_idx], &gt_rows[right_idx]);
            targets_by_z
                .entry(z)
                .or_default()
                .push((left_pos, right_pos));
        }
    }

    let negative_pairs = ii_indices.len() * jj_indices.len();
    for (z, sources) in sources_by_z {
        let Some(targets) = targets_by_z.get(&z) else {
            return FiberOutcome {
                negative_pairs,
                failure: Some(RankFailure {
                    beta: Vec::new(),
                    message: "no targets with same pair-sum".to_string(),
                    z,
                    source_rank: 0,
                    source: sources[0],
                    target: None,
                }),
            };
        };
        if targets.len() < sources.len() {
            return FiberOutcome {
                negative_pairs,
                failure: Some(RankFailure {
                    beta: Vec::new(),
                    message: format!(
                        "not enough fixed-sum targets: sources={}, targets={}",
                        sources.len(),
                        targets.len()
                    ),
                    z,
                    source_rank: targets.len(),
                    source: sources[targets.len()],
                    target: None,
                }),
            };
        }
        let target_rank_by_pair: BTreeMap<_, _> = targets
            .iter()
            .copied()
            .enumerate()
            .map(|(rank, pair)| (pair, rank))
            .collect();
        let mut used_targets = BTreeSet::new();
        for (rank, &source) in sources.iter().enumerate() {
            let target = if require_carrier {
                let edge_targets: Vec<_> = targets
                    .iter()
                    .copied()
                    .filter(|target| {
                        is_carrier_edge(
                            &gt_rows[ii_indices[source.0]],
                            &gt_rows[jj_indices[source.1]],
                            &gt_rows[mixed_indices[target.0]],
                            &gt_rows[mixed_indices[target.1]],
                            lower_label,
                        )
                    })
                    .collect();
                let Some(target) = select_target(
                    rank,
                    source,
                    ii_indices.len(),
                    jj_indices.len(),
                    mixed_indices.len(),
                    &edge_targets,
                    &target_rank_by_pair,
                    selection,
                ) else {
                    return FiberOutcome {
                        negative_pairs,
                        failure: Some(RankFailure {
                            beta: Vec::new(),
                            message: "source has no carrier target in fixed-sum fiber".to_string(),
                            z,
                            source_rank: rank,
                            source,
                            target: None,
                        }),
                    };
                };
                target
            } else {
                targets[rank]
            };

            if !used_targets.insert(target) {
                return FiberOutcome {
                    negative_pairs,
                    failure: Some(RankFailure {
                        beta: Vec::new(),
                        message: format!("{selection:?} selection is not injective"),
                        z,
                        source_rank: rank,
                        source,
                        target: Some(target),
                    }),
                };
            }
        }
    }

    FiberOutcome {
        negative_pairs,
        failure: None,
    }
}

fn select_target(
    source_rank: usize,
    source: (usize, usize),
    ii_count: usize,
    jj_count: usize,
    mixed_count: usize,
    edge_targets: &[(usize, usize)],
    target_rank_by_pair: &BTreeMap<(usize, usize), usize>,
    selection: SelectionMode,
) -> Option<(usize, usize)> {
    match selection {
        SelectionMode::SameRank => edge_targets
            .iter()
            .copied()
            .find(|target| target_rank_by_pair[target] == source_rank),
        SelectionMode::Nearest => edge_targets.iter().copied().min_by_key(|target| {
            let target_rank = target_rank_by_pair[target];
            (target_rank.abs_diff(source_rank), target_rank, *target)
        }),
        SelectionMode::Ceiling => edge_targets
            .iter()
            .copied()
            .filter(|target| target_rank_by_pair[target] >= source_rank)
            .min_by_key(|target| (target_rank_by_pair[target], *target))
            .or_else(|| {
                edge_targets
                    .iter()
                    .copied()
                    .max_by_key(|target| (target_rank_by_pair[target], *target))
            }),
        SelectionMode::Floor => edge_targets
            .iter()
            .copied()
            .filter(|target| target_rank_by_pair[target] <= source_rank)
            .max_by_key(|target| (target_rank_by_pair[target], *target))
            .or_else(|| {
                edge_targets
                    .iter()
                    .copied()
                    .min_by_key(|target| (target_rank_by_pair[target], *target))
            }),
        SelectionMode::FirstEdge => edge_targets.first().copied(),
        SelectionMode::LastEdge => edge_targets.last().copied(),
        SelectionMode::ClosestSource => edge_targets.iter().copied().min_by_key(|target| {
            (
                target.0.abs_diff(source.0) + target.1.abs_diff(source.1),
                target.0.abs_diff(source.0),
                target.1.abs_diff(source.1),
                *target,
            )
        }),
        SelectionMode::ClosestScaledSource => {
            let source_left = scaled_position(source.0, ii_count, mixed_count);
            let source_right = scaled_position(source.1, jj_count, mixed_count);
            edge_targets.iter().copied().min_by_key(|target| {
                let target_left = 2 * target.0 as i64;
                let target_right = 2 * target.1 as i64;
                (
                    (target_left - source_left).abs() + (target_right - source_right).abs(),
                    (target_left - source_left).abs(),
                    (target_right - source_right).abs(),
                    *target,
                )
            })
        }
    }
}

fn scaled_position(pos: usize, source_count: usize, mixed_count: usize) -> i64 {
    if source_count <= 1 || mixed_count <= 1 {
        return 0;
    }
    // Twice the rounded linear embedding of a source index into mixed-index space.
    2 * pos as i64 * (mixed_count as i64 - 1) / (source_count as i64 - 1)
}

fn is_carrier_edge(
    left: &[Vec<u32>],
    right: &[Vec<u32>],
    target_left: &[Vec<u32>],
    target_right: &[Vec<u32>],
    lower_label: usize,
) -> bool {
    if add_gt_rows(left, right) != add_gt_rows(target_left, target_right) {
        return false;
    }
    let d: Vec<[i32; 2]> = left
        .iter()
        .zip(target_left)
        .map(|(left_row, target_row)| {
            [
                left_row[0] as i32 - target_row[0] as i32,
                left_row[1] as i32 - target_row[1] as i32,
            ]
        })
        .collect();
    if d.len() != right.len() || d.len() != target_right.len() || lower_label >= d.len() {
        return false;
    }
    for level in 0..d.len() {
        let expected_from_right = [
            target_right[level][0] as i32 - right[level][0] as i32,
            target_right[level][1] as i32 - right[level][1] as i32,
        ];
        if d[level] != expected_from_right {
            return false;
        }
    }
    for row in d.iter().take(lower_label + 1) {
        if row.iter().any(|entry| entry.abs() > 1) {
            return false;
        }
    }
    if d.iter().skip(lower_label + 1).any(|&row| row != [0, 0]) {
        return false;
    }
    if d[lower_label] != [1, 0] && d[lower_label] != [0, 1] {
        return false;
    }
    let start = d.iter().position(|&row| row != [0, 0]);
    let Some(start) = start else {
        return false;
    };
    if start > lower_label {
        return false;
    }
    d.iter()
        .take(lower_label)
        .skip(start)
        .all(|&row| row == [1, -1] || row == [-1, 1])
        && d.iter().take(start).all(|&row| row == [0, 0])
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
    println!("alphabet={}", flagged_shape.alphabet());
    println!("fibers_checked={}", outcome.fibers_checked);
    println!("negative_pairs={}", outcome.negative_pairs);
    println!("elapsed={elapsed_seconds:.3}s");
    if let Some(failure) = &outcome.failure {
        println!("FAIL");
        println!(
            "beta={:?}, {}, source_rank={}, source={:?}, target={:?}, z={}",
            failure.beta,
            failure.message,
            failure.source_rank,
            failure.source,
            failure.target,
            summarize_gt(&failure.z)
        );
    } else {
        println!("PASS");
    }
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
    "Check the fixed-pair-sum rank map for two-row carrier fibers.

USAGE:
  two_row_fixed_sum_rank_scan [OPTIONS]

OPTIONS:
  --lambda PARTS          Outer partition, e.g. 4,2. If omitted, scan a family.
  --mu PARTS              Inner partition, e.g. 2. Default: empty.
  --row-flags PARTS       Row upper flags. Omit for ordinary skew tableaux.
  --alphabet N            Alphabet size. Default: 5.
  --lower-label I         Use active labels I,I+1. Default: 4.
  --max-skew-size N       Maximum skew size for family scans. Default: 7.
  --max-outer-extra N     Family outer sizes go up to skew_size + this. Default: 8.
  --at-most-two-rows      Include one-row shapes.
  --allow-noncarrier      Only check fixed-pair-sum target availability.
  --selection MODE        Carrier target choice: same-rank, nearest, ceiling,
                          floor, first-edge, last-edge, closest-source, or
                          closest-scaled-source. Default: same-rank.
  --tableau-limit N       Skip shapes whose tableau enumeration exceeds N.
  --help                  Print this help."
}
