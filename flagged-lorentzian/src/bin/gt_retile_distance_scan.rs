use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use flagged_lorentzian::{
    active_subword_descent_data_for_values, add_patterns, descent_data_for_values,
    elementary_row_exchange_neighbors, enumerate_tableaux, families::skew_shapes_of_size, gt,
    pair_envelope, sharp_flag, DescentData, DescentStatistic, RowFlaggedSkewShape, SkewGtPattern,
    SkewShape, TableauRecord,
};

type Content = Vec<u32>;
type GtSum = Vec<Vec<u32>>;
type ActiveRow = Vec<u32>;
type Envelope = Vec<u32>;
type DescentPair = (DescentData, DescentData);

#[derive(Debug, Clone)]
struct Args {
    lambda: Option<Vec<u32>>,
    mu: Vec<u32>,
    row_flags: Option<Vec<u32>>,
    alphabet: usize,
    lower_label: u32,
    max_skew_size: u32,
    max_outer_extra: u32,
    max_rows: Option<usize>,
    max_exchange_depth: usize,
    stop_at_l1: Option<u32>,
    connected_only: bool,
    descent_mode: DescentMode,
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
            max_skew_size: 6,
            max_outer_extra: 7,
            max_rows: None,
            max_exchange_depth: 3,
            stop_at_l1: None,
            connected_only: false,
            descent_mode: DescentMode::None,
            tableau_limit: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DescentMode {
    None,
    Global,
    Componentwise,
    ActiveGlobal,
    ActiveComponentwise,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct FullKey {
    active_row: ActiveRow,
    envelope: Envelope,
    descent_pair: Option<DescentPair>,
    gt_sum: GtSum,
}

#[derive(Debug, Clone)]
struct PairData {
    active_row: ActiveRow,
    envelope: Envelope,
    descent_pair: Option<DescentPair>,
    gt_sum: GtSum,
}

#[derive(Debug, Clone)]
struct DistanceSummary {
    l1: u32,
    rows_changed: usize,
    diff: Vec<Vec<i32>>,
}

#[derive(Debug, Clone)]
struct ShapeOutcome {
    fibers_checked: usize,
    full_sum_deficits: usize,
    worst_distance: Option<DistanceSummary>,
    worst_exchange_depth: Option<usize>,
    failure: Option<String>,
    limit_exceeded: bool,
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
    let mut fibers_checked = 0usize;
    let mut full_sum_deficits = 0usize;
    let mut skipped_by_limit = 0usize;
    let mut worst_distance: Option<DistanceSummary> = None;
    let mut worst_exchange_depth = None;

    for skew_size in 0..=args.max_skew_size {
        let max_outer_size = skew_size + args.max_outer_extra;
        for shape in skew_shapes_of_size(skew_size, max_outer_size) {
            if let Some(max_rows) = args.max_rows {
                if shape.row_count() > max_rows {
                    continue;
                }
            }
            if args.connected_only && !shape.is_connected() {
                continue;
            }
            let flagged_shape = flagged_shape_for_args(shape, &args);
            let outcome = scan_one_shape(&flagged_shape, &args);
            if outcome.limit_exceeded {
                skipped_by_limit += 1;
                continue;
            }
            if let Some(failure) = outcome.failure {
                println!("FAIL");
                println!("outer={:?}", flagged_shape.shape().outer().parts());
                println!("inner={:?}", flagged_shape.shape().inner().parts());
                println!("flags={:?}", flagged_shape.row_flags());
                println!("{failure}");
                std::process::exit(1);
            }

            shapes_checked += 1;
            fibers_checked += outcome.fibers_checked;
            full_sum_deficits += outcome.full_sum_deficits;
            worst_distance = max_distance(worst_distance, outcome.worst_distance);
            worst_exchange_depth =
                max_optional_usize(worst_exchange_depth, outcome.worst_exchange_depth);
        }
        println!(
            "skew_size={skew_size}: shapes_checked={shapes_checked}, fibers_checked={fibers_checked}, full_sum_deficits={full_sum_deficits}, skipped_by_limit={skipped_by_limit}, worst={}, worst_exchange_depth={}, elapsed={:.3}s",
            format_distance(worst_distance.as_ref()),
            format_optional_usize(worst_exchange_depth),
            started.elapsed().as_secs_f64()
        );
    }

    println!("PASS");
    println!(
        "shapes_checked={shapes_checked}, fibers_checked={fibers_checked}, full_sum_deficits={full_sum_deficits}, skipped_by_limit={skipped_by_limit}, worst={}, worst_exchange_depth={}, elapsed={:.3}s",
        format_distance(worst_distance.as_ref()),
        format_optional_usize(worst_exchange_depth),
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
                "--max-rows" => {
                    args.max_rows = Some(
                        take_value(&mut iter, &flag)?
                            .parse()
                            .map_err(|err| format!("invalid --max-rows: {err}"))?,
                    )
                }
                "--max-exchange-depth" => {
                    args.max_exchange_depth = take_value(&mut iter, &flag)?
                        .parse()
                        .map_err(|err| format!("invalid --max-exchange-depth: {err}"))?
                }
                "--stop-at-l1" => {
                    args.stop_at_l1 = Some(
                        take_value(&mut iter, &flag)?
                            .parse()
                            .map_err(|err| format!("invalid --stop-at-l1: {err}"))?,
                    )
                }
                "--connected-only" => args.connected_only = true,
                "--descent-mode" => {
                    args.descent_mode = DescentMode::parse(&take_value(&mut iter, &flag)?)?
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

impl DescentMode {
    fn parse(input: &str) -> Result<Self, String> {
        match input {
            "none" => Ok(Self::None),
            "global" => Ok(Self::Global),
            "componentwise" => Ok(Self::Componentwise),
            "active-global" => Ok(Self::ActiveGlobal),
            "active-componentwise" => Ok(Self::ActiveComponentwise),
            other => Err(format!(
                "invalid --descent-mode `{other}`; expected `none`, `global`, `componentwise`, `active-global`, or `active-componentwise`"
            )),
        }
    }

    fn reading_orders(self, shape: &SkewShape) -> Option<Vec<Vec<usize>>> {
        match self {
            Self::None => None,
            Self::Global | Self::ActiveGlobal => {
                Some(DescentStatistic::Global.reading_orders(shape))
            }
            Self::Componentwise | Self::ActiveComponentwise => {
                Some(DescentStatistic::Componentwise.reading_orders(shape))
            }
        }
    }

    fn descent_data(
        self,
        values: &[u32],
        reading_orders: &[Vec<usize>],
        lower_label: u32,
    ) -> Option<DescentData> {
        match self {
            Self::None => None,
            Self::Global | Self::Componentwise => {
                Some(descent_data_for_values(values, reading_orders))
            }
            Self::ActiveGlobal | Self::ActiveComponentwise => Some(
                active_subword_descent_data_for_values(values, reading_orders, lower_label),
            ),
        }
    }
}

fn scan_one_shape(flagged_shape: &RowFlaggedSkewShape, args: &Args) -> ShapeOutcome {
    let tableaux = match enumerate_tableaux(flagged_shape, args.tableau_limit) {
        Ok(tableaux) => tableaux,
        Err(_) => {
            return ShapeOutcome {
                fibers_checked: 0,
                full_sum_deficits: 0,
                worst_distance: None,
                worst_exchange_depth: None,
                failure: None,
                limit_exceeded: true,
            }
        }
    };

    let reading_orders = args.descent_mode.reading_orders(flagged_shape.shape());
    let pair_data: Vec<_> = tableaux
        .iter()
        .map(|tableau| {
            single_data(
                flagged_shape.shape(),
                tableau,
                args.alphabet,
                args.lower_label,
                args.descent_mode,
                reading_orders.as_deref(),
            )
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
    let mut full_sum_deficits = 0usize;
    let mut worst_distance = None;
    let mut worst_exchange_depth = None;

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
            &pair_data,
            ii_indices,
            mixed_indices,
            jj_indices,
            args.lower_label,
            args.max_exchange_depth,
            args.stop_at_l1,
        );
        if let Some(failure) = outcome.failure {
            return ShapeOutcome {
                fibers_checked,
                full_sum_deficits,
                worst_distance,
                worst_exchange_depth,
                failure: Some(format!("beta={beta:?}: {failure}")),
                limit_exceeded: false,
            };
        }
        full_sum_deficits += outcome.full_sum_deficits;
        worst_distance = max_distance(worst_distance, outcome.worst_distance);
        worst_exchange_depth =
            max_optional_usize(worst_exchange_depth, outcome.worst_exchange_depth);
    }

    ShapeOutcome {
        fibers_checked,
        full_sum_deficits,
        worst_distance,
        worst_exchange_depth,
        failure: None,
        limit_exceeded: false,
    }
}

#[derive(Debug, Clone)]
struct FiberOutcome {
    full_sum_deficits: usize,
    worst_distance: Option<DistanceSummary>,
    worst_exchange_depth: Option<usize>,
    failure: Option<String>,
}

fn scan_fiber(
    single_data: &[SingleData],
    ii_indices: &[usize],
    mixed_indices: &[usize],
    jj_indices: &[usize],
    lower_label: u32,
    max_exchange_depth: usize,
    stop_at_l1: Option<u32>,
) -> FiberOutcome {
    let negative_pairs = pair_data(single_data, ii_indices, jj_indices, lower_label);
    let positive_pairs = pair_data(single_data, mixed_indices, mixed_indices, lower_label);

    let mut positive_by_full_key = BTreeMap::<FullKey, usize>::new();
    for pair in &positive_pairs {
        *positive_by_full_key.entry(pair.full_key()).or_insert(0) += 1;
    }

    let mut negative_by_full_key = BTreeMap::<FullKey, usize>::new();
    for pair in &negative_pairs {
        *negative_by_full_key.entry(pair.full_key()).or_insert(0) += 1;
    }

    let mut deficits = 0usize;
    let mut worst_distance = None;
    let mut worst_exchange_depth = None;

    for (key, negative_count) in negative_by_full_key {
        let positive_count = positive_by_full_key.get(&key).copied().unwrap_or(0);
        if negative_count <= positive_count {
            continue;
        }
        deficits += negative_count - positive_count;

        let candidate_targets = active_envelope_targets_excluding_full_key(&key, &positive_pairs);
        let Some(distance) = closest_distance(&key, &candidate_targets) else {
            return FiberOutcome {
                full_sum_deficits: deficits,
                worst_distance,
                worst_exchange_depth,
                failure: Some(format!(
                    "no alternate active-row/envelope target for full-key deficit {key:?}"
                )),
            };
        };
        if stop_at_l1.is_some_and(|threshold| distance.l1 >= threshold) {
            worst_distance = max_distance(worst_distance, Some(distance.clone()));
            return FiberOutcome {
                full_sum_deficits: deficits,
                worst_distance,
                worst_exchange_depth,
                failure: Some(format!(
                    "stopped at distance {}; key_active={:?}, key_envelope={:?}, diff={}",
                    distance.l1,
                    key.active_row,
                    key.envelope,
                    summarize_diff(&distance.diff)
                )),
            };
        }
        worst_distance = max_distance(worst_distance, Some(distance));

        let Some(exchange_depth) = exchange_depth_to_any(
            &key.gt_sum,
            &candidate_targets,
            lower_label as usize,
            max_exchange_depth,
        ) else {
            return FiberOutcome {
                full_sum_deficits: deficits,
                worst_distance,
                worst_exchange_depth,
                failure: Some(format!(
                    "no elementary-exchange path within depth {max_exchange_depth} for full-key deficit {key:?}"
                )),
            };
        };
        worst_exchange_depth = max_optional_usize(worst_exchange_depth, Some(exchange_depth));
    }

    FiberOutcome {
        full_sum_deficits: deficits,
        worst_distance,
        worst_exchange_depth,
        failure: None,
    }
}

fn active_envelope_targets_excluding_full_key(
    key: &FullKey,
    positive_pairs: &[PairData],
) -> Vec<GtSum> {
    let mut targets: Vec<_> = positive_pairs
        .iter()
        .filter(|pair| pair.active_row == key.active_row && pair.envelope == key.envelope)
        .filter(|pair| pair.descent_pair == key.descent_pair)
        .filter(|pair| pair.gt_sum != key.gt_sum)
        .map(|pair| pair.gt_sum.clone())
        .collect();
    targets.sort();
    targets.dedup();
    targets
}

fn closest_distance(key: &FullKey, targets: &[GtSum]) -> Option<DistanceSummary> {
    targets
        .iter()
        .map(|target| distance_between(&key.gt_sum, target))
        .min_by_key(|distance| (distance.l1, distance.rows_changed))
}

fn exchange_depth_to_any(
    start: &[Vec<u32>],
    targets: &[GtSum],
    fixed_level: usize,
    max_depth: usize,
) -> Option<usize> {
    let target_set: BTreeSet<_> = targets.iter().cloned().collect();
    if target_set.contains(start) {
        return Some(0);
    }

    let mut seen = BTreeSet::from([start.to_vec()]);
    let mut frontier = vec![start.to_vec()];
    for depth in 1..=max_depth {
        let mut next = Vec::new();
        for array in frontier {
            for neighbor in elementary_row_exchange_neighbors(&array, fixed_level) {
                if !seen.insert(neighbor.clone()) {
                    continue;
                }
                if target_set.contains(&neighbor) {
                    return Some(depth);
                }
                next.push(neighbor);
            }
        }
        frontier = next;
    }
    None
}

#[derive(Debug, Clone)]
struct SingleData {
    gt: SkewGtPattern,
    sharp_flag: Envelope,
    descent: Option<DescentData>,
}

fn single_data(
    shape: &SkewShape,
    tableau: &TableauRecord,
    alphabet: usize,
    lower_label: u32,
    descent_mode: DescentMode,
    reading_orders: Option<&[Vec<usize>]>,
) -> SingleData {
    SingleData {
        gt: SkewGtPattern::from_tableau(shape, &tableau.values, alphabet),
        sharp_flag: sharp_flag(shape, &tableau.values),
        descent: reading_orders.and_then(|reading_orders| {
            descent_mode.descent_data(&tableau.values, reading_orders, lower_label)
        }),
    }
}

impl PairData {
    fn full_key(&self) -> FullKey {
        FullKey {
            active_row: self.active_row.clone(),
            envelope: self.envelope.clone(),
            descent_pair: self.descent_pair.clone(),
            gt_sum: self.gt_sum.clone(),
        }
    }
}

fn pair_data(
    single_data: &[SingleData],
    left_indices: &[usize],
    right_indices: &[usize],
    lower_label: u32,
) -> Vec<PairData> {
    let mut out = Vec::new();
    for &left_idx in left_indices {
        for &right_idx in right_indices {
            let left = &single_data[left_idx];
            let right = &single_data[right_idx];
            let gt_sum = add_patterns(&left.gt, &right.gt);
            out.push(PairData {
                active_row: gt_sum[lower_label as usize].clone(),
                envelope: pair_envelope(&left.sharp_flag, &right.sharp_flag),
                descent_pair: unordered_descent_pair(left.descent.as_ref(), right.descent.as_ref()),
                gt_sum,
            });
        }
    }
    out
}

fn unordered_descent_pair(
    left: Option<&DescentData>,
    right: Option<&DescentData>,
) -> Option<DescentPair> {
    let (left, right) = (left?.clone(), right?.clone());
    if left <= right {
        Some((left, right))
    } else {
        Some((right, left))
    }
}

fn distance_between(left: &[Vec<u32>], right: &[Vec<u32>]) -> DistanceSummary {
    let diff = gt::subtract_pattern_sums(left, right);
    DistanceSummary {
        l1: diff
            .iter()
            .flatten()
            .map(|entry| entry.unsigned_abs())
            .sum(),
        rows_changed: diff
            .iter()
            .filter(|row| row.iter().any(|&entry| entry != 0))
            .count(),
        diff,
    }
}

fn max_distance(
    left: Option<DistanceSummary>,
    right: Option<DistanceSummary>,
) -> Option<DistanceSummary> {
    match (left, right) {
        (None, None) => None,
        (Some(distance), None) | (None, Some(distance)) => Some(distance),
        (Some(left), Some(right)) => {
            if (right.l1, right.rows_changed) > (left.l1, left.rows_changed) {
                Some(right)
            } else {
                Some(left)
            }
        }
    }
}

fn max_optional_usize(left: Option<usize>, right: Option<usize>) -> Option<usize> {
    match (left, right) {
        (None, None) => None,
        (Some(value), None) | (None, Some(value)) => Some(value),
        (Some(left), Some(right)) => Some(left.max(right)),
    }
}

fn format_optional_usize(value: Option<usize>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "none".to_string())
}

fn format_distance(distance: Option<&DistanceSummary>) -> String {
    match distance {
        Some(distance) => format!(
            "l1={}, rows={}, diff={}",
            distance.l1,
            distance.rows_changed,
            summarize_diff(&distance.diff)
        ),
        None => "none".to_string(),
    }
}

fn summarize_diff(diff: &[Vec<i32>]) -> String {
    let parts: Vec<_> = diff
        .iter()
        .enumerate()
        .filter(|(_, row)| row.iter().any(|&entry| entry != 0))
        .map(|(level, row)| format!("r{level}={row:?}"))
        .collect();
    if parts.is_empty() {
        "0".to_string()
    } else {
        parts.join(";")
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
    println!("full_sum_deficits={}", outcome.full_sum_deficits);
    println!("worst={}", format_distance(outcome.worst_distance.as_ref()));
    println!(
        "worst_exchange_depth={}",
        format_optional_usize(outcome.worst_exchange_depth)
    );
    println!("elapsed={elapsed_seconds:.3}s");
    if let Some(failure) = &outcome.failure {
        println!("FAIL");
        println!("{failure}");
    } else {
        println!("PASS");
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
    "Measure how far full-pair-GT-sum deficits are from active-row/envelope targets.

USAGE:
  gt_retile_distance_scan [OPTIONS]

OPTIONS:
  --lambda PARTS          Outer partition, e.g. 4,2,1. If omitted, scan a family.
  --mu PARTS              Inner partition, e.g. 2. Default: empty.
  --row-flags PARTS       Row upper flags. Omit for ordinary skew tableaux.
  --alphabet N            Alphabet size. Default: 5.
  --lower-label I         Use active row for labels I,I+1. Default: 4.
  --max-skew-size N       Maximum skew size for family scans. Default: 6.
  --max-outer-extra N     Family outer sizes go up to skew_size + this. Default: 7.
  --max-rows N            Restrict family scans to shapes with at most N rows.
  --max-exchange-depth N  Maximum elementary GT exchange depth. Default: 3.
  --stop-at-l1 N          Stop and print the first deficit with nearest L1 at least N.
  --connected-only        Restrict family scans to connected skew shapes.
  --descent-mode MODE     Preserve no descents, `global`, `componentwise`,
                          `active-global`, or `active-componentwise`. Default: none.
  --tableau-limit N       Skip shapes whose tableau enumeration exceeds N.
  --help                  Print this help."
}
