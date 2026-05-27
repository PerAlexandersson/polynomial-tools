use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use flagged_lorentzian::{
    enumerate_tableaux, families::skew_shapes_of_size, is_gt_array, RowFlaggedSkewShape,
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
    lower_label: usize,
    max_skew_size: u32,
    max_outer_extra: u32,
    max_rows: usize,
    exactly_three_rows: bool,
    check_matching: bool,
    adjacent_roots: bool,
    allow_pauses: bool,
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
            max_skew_size: 5,
            max_outer_extra: 7,
            max_rows: 3,
            exactly_three_rows: false,
            check_matching: false,
            adjacent_roots: false,
            allow_pauses: false,
            tableau_limit: None,
        }
    }
}

#[derive(Debug, Clone)]
struct ShapeOutcome {
    fibers_checked: usize,
    negative_pairs: usize,
    failures: Vec<CarrierFailure>,
    limit_exceeded: bool,
}

#[derive(Debug, Clone)]
struct CarrierFailure {
    beta: Content,
    message: String,
    left_pos: usize,
    right_pos: usize,
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
    let started = Instant::now();

    if let Some(lambda) = args.lambda.clone() {
        let shape = SkewShape::from_parts(lambda, args.mu.clone());
        let flagged_shape = flagged_shape_for_args(shape, &args);
        let outcome = scan_one_shape(&flagged_shape, &args);
        print_single_outcome(&flagged_shape, &outcome, started.elapsed().as_secs_f64());
        std::process::exit(if outcome.failures.is_empty() { 0 } else { 1 });
    }

    let mut shapes_checked = 0usize;
    let mut skipped_by_limit = 0usize;
    let mut fibers_checked = 0usize;
    let mut negative_pairs = 0usize;

    for skew_size in 0..=args.max_skew_size {
        let max_outer_size = skew_size + args.max_outer_extra;
        for shape in skew_shapes_of_size(skew_size, max_outer_size) {
            if shape.row_count() > args.max_rows {
                continue;
            }
            if args.exactly_three_rows && shape.row_count() != 3 {
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
            if let Some(failure) = outcome.failures.first() {
                println!("FAIL");
                println!("outer={:?}", flagged_shape.shape().outer().parts());
                println!("inner={:?}", flagged_shape.shape().inner().parts());
                println!("flags={:?}", flagged_shape.row_flags());
                println!(
                    "beta={:?}, {}, first_pair=({}, {})",
                    failure.beta, failure.message, failure.left_pos, failure.right_pos
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

fn scan_one_shape(flagged_shape: &RowFlaggedSkewShape, args: &Args) -> ShapeOutcome {
    let tableaux = match enumerate_tableaux(flagged_shape, args.tableau_limit) {
        Ok(tableaux) => tableaux,
        Err(_) => {
            return ShapeOutcome {
                fibers_checked: 0,
                negative_pairs: 0,
                failures: Vec::new(),
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

    let lower = args.lower_label - 1;
    let upper = lower + 1;
    let mut seen_beta = BTreeSet::new();
    let mut fibers_checked = 0usize;
    let mut negative_pairs = 0usize;
    let mut failures = Vec::new();

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

        let (Some(ii_indices), Some(jj_indices), Some(mixed_indices)) = (
            by_content.get(&ii),
            by_content.get(&jj),
            by_content.get(&mixed),
        ) else {
            continue;
        };

        let mixed_gt_to_pos: BTreeMap<_, _> = mixed_indices
            .iter()
            .enumerate()
            .map(|(pos, &idx)| (gt_rows[idx].clone(), pos))
            .collect();

        fibers_checked += 1;
        let source_count = ii_indices.len() * jj_indices.len();
        negative_pairs += source_count;
        let failure = if args.check_matching {
            carrier_matching_failure(
                ii_indices,
                jj_indices,
                mixed_indices.len(),
                &gt_rows,
                &mixed_gt_to_pos,
                args,
            )
        } else {
            carrier_existence_failure(ii_indices, jj_indices, &gt_rows, &mixed_gt_to_pos, args)
        };
        if let Some((left_pos, right_pos, message)) = failure {
            failures.push(CarrierFailure {
                beta,
                message,
                left_pos,
                right_pos,
            });
            return ShapeOutcome {
                fibers_checked,
                negative_pairs,
                failures,
                limit_exceeded: false,
            };
        }
    }

    ShapeOutcome {
        fibers_checked,
        negative_pairs,
        failures,
        limit_exceeded: false,
    }
}

fn carrier_existence_failure(
    ii_indices: &[usize],
    jj_indices: &[usize],
    gt_rows: &[GtRows],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    args: &Args,
) -> Option<(usize, usize, String)> {
    for (left_pos, &left_idx) in ii_indices.iter().enumerate() {
        for (right_pos, &right_idx) in jj_indices.iter().enumerate() {
            let targets = carrier_targets(
                &gt_rows[left_idx],
                &gt_rows[right_idx],
                mixed_gt_to_pos,
                1,
                args,
            );
            if targets.is_empty() {
                return Some((left_pos, right_pos, "no root-carrier target".to_string()));
            }
        }
    }
    None
}

fn carrier_matching_failure(
    ii_indices: &[usize],
    jj_indices: &[usize],
    mixed_count: usize,
    gt_rows: &[GtRows],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    args: &Args,
) -> Option<(usize, usize, String)> {
    let source_pairs: Vec<_> = (0..ii_indices.len())
        .flat_map(|left_pos| (0..jj_indices.len()).map(move |right_pos| (left_pos, right_pos)))
        .collect();
    let mut edges_by_source = Vec::with_capacity(source_pairs.len());
    for &(left_pos, right_pos) in &source_pairs {
        let left_idx = ii_indices[left_pos];
        let right_idx = jj_indices[right_pos];
        let edges = carrier_targets(
            &gt_rows[left_idx],
            &gt_rows[right_idx],
            mixed_gt_to_pos,
            mixed_count,
            args,
        );
        if edges.is_empty() {
            return Some((left_pos, right_pos, "no root-carrier target".to_string()));
        }
        edges_by_source.push(edges);
    }

    let mut target_match = vec![None; mixed_count * mixed_count];
    for source in 0..source_pairs.len() {
        let mut seen = vec![false; target_match.len()];
        if !augment_matching(source, &edges_by_source, &mut seen, &mut target_match) {
            let (left_pos, right_pos) = source_pairs[source];
            return Some((
                left_pos,
                right_pos,
                format!("root-carrier graph has no augmenting match at source {source}"),
            ));
        }
    }
    None
}

fn augment_matching(
    source: usize,
    edges_by_source: &[Vec<usize>],
    seen: &mut [bool],
    target_match: &mut [Option<usize>],
) -> bool {
    for &target in &edges_by_source[source] {
        if seen[target] {
            continue;
        }
        seen[target] = true;
        let can_use = match target_match[target] {
            None => true,
            Some(previous) => augment_matching(previous, edges_by_source, seen, target_match),
        };
        if can_use {
            target_match[target] = Some(source);
            return true;
        }
    }
    false
}

fn carrier_targets(
    left: &[Vec<u32>],
    right: &[Vec<u32>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    mixed_count: usize,
    args: &Args,
) -> Vec<usize> {
    if left.is_empty()
        || right.is_empty()
        || left[0].len() != right[0].len()
        || args.lower_label >= left.len()
    {
        return Vec::new();
    }

    let mut out = BTreeSet::new();
    for carrier in carrier_words(
        left[0].len(),
        args.lower_label,
        args.adjacent_roots,
        args.allow_pauses,
    ) {
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
        out.insert(left_pos * mixed_count + right_pos);
    }
    out.into_iter().collect()
}

fn carrier_words(
    width: usize,
    lower_label: usize,
    adjacent_roots: bool,
    allow_pauses: bool,
) -> Vec<Vec<Vec<i32>>> {
    let mut carriers = Vec::new();
    let roots = root_vectors(width, adjacent_roots);
    let terminal_units = unit_vectors(width);
    for start_level in (1..=lower_label).rev() {
        let mut carrier = vec![vec![0; width]; lower_label + 1];
        extend_carrier(
            &mut carriers,
            &mut carrier,
            &roots,
            &terminal_units,
            start_level,
            lower_label,
            allow_pauses,
        );
    }
    carriers
}

fn extend_carrier(
    carriers: &mut Vec<Vec<Vec<i32>>>,
    carrier: &mut Vec<Vec<i32>>,
    roots: &[Vec<i32>],
    terminal_units: &[Vec<i32>],
    level: usize,
    lower_label: usize,
    allow_pauses: bool,
) {
    if level == lower_label {
        for unit in terminal_units {
            carrier[level] = unit.clone();
            carriers.push(carrier.clone());
        }
        carrier[level].fill(0);
        return;
    }

    for root in roots {
        carrier[level] = root.clone();
        extend_carrier(
            carriers,
            carrier,
            roots,
            terminal_units,
            level + 1,
            lower_label,
            allow_pauses,
        );
    }
    if allow_pauses {
        carrier[level].fill(0);
        extend_carrier(
            carriers,
            carrier,
            roots,
            terminal_units,
            level + 1,
            lower_label,
            allow_pauses,
        );
    }
    carrier[level].fill(0);
}

fn root_vectors(width: usize, adjacent_only: bool) -> Vec<Vec<i32>> {
    let mut roots = Vec::new();
    for source in 0..width {
        for target in 0..width {
            if source == target {
                continue;
            }
            if adjacent_only && source.abs_diff(target) != 1 {
                continue;
            }
            let mut root = vec![0; width];
            root[source] = 1;
            root[target] = -1;
            roots.push(root);
        }
    }
    roots
}

fn unit_vectors(width: usize) -> Vec<Vec<i32>> {
    (0..width)
        .map(|idx| {
            let mut unit = vec![0; width];
            unit[idx] = 1;
            unit
        })
        .collect()
}

fn apply_carrier(rows: &[Vec<u32>], carrier: &[Vec<i32>], sign: i32) -> Option<GtRows> {
    let mut out = rows.to_vec();
    for (row, delta) in out.iter_mut().zip(carrier) {
        for col in 0..row.len() {
            let value = row[col] as i32 + sign * delta[col];
            if value < 0 {
                return None;
            }
            row[col] = value as u32;
        }
    }
    is_gt_array(&out).then_some(out)
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
    if let Some(failure) = outcome.failures.first() {
        println!("FAIL");
        println!(
            "beta={:?}, {}, first_pair=({}, {})",
            failure.beta, failure.message, failure.left_pos, failure.right_pos
        );
    } else {
        println!("PASS");
    }
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
            "--lambda" => args.lambda = Some(parse_u32_vec(&take_value(&mut iter, &flag)?)?),
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
            "--max-rows" => {
                args.max_rows = take_value(&mut iter, &flag)?
                    .parse()
                    .map_err(|err| format!("invalid --max-rows: {err}"))?
            }
            "--exactly-three-rows" => args.exactly_three_rows = true,
            "--check-matching" => args.check_matching = true,
            "--adjacent-roots" => args.adjacent_roots = true,
            "--allow-pauses" => args.allow_pauses = true,
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
    if args.lower_label == 0 || args.lower_label >= args.alphabet {
        return Err("--lower-label must be between 1 and alphabet-1".to_string());
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
    "Check root-valued carrier variants for up to three-row shapes.

USAGE:
  three_row_root_carrier_scan [OPTIONS]

OPTIONS:
  --lambda PARTS          Outer partition. If omitted, scan a family.
  --mu PARTS              Inner partition. Default: empty.
  --row-flags PARTS       Row upper flags. Omit for ordinary skew tableaux.
  --alphabet N            Alphabet size. Default: 5.
  --lower-label I         Active labels are I,I+1. Default: 4.
  --max-skew-size N       Maximum skew size for family scans. Default: 5.
  --max-outer-extra N     Family outer sizes go up to skew_size + this. Default: 7.
  --max-rows N            Restrict family scans to at most N rows. Default: 3.
  --exactly-three-rows    Skip shapes that do not have exactly three rows.
  --check-matching        Check that the root-carrier graph has a matching.
  --adjacent-roots        Only allow root states e_a-e_b with adjacent rows.
  --allow-pauses          Allow zero carrier states between start and active level.
  --tableau-limit N       Skip shapes whose tableau enumeration exceeds N.
  --help                  Print this help."
}
