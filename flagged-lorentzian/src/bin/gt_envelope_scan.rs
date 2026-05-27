use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use flagged_lorentzian::{
    enumerate_tableaux, families::skew_shapes_of_size, gt, RowFlaggedSkewShape, SkewShape,
    TableauRecord,
};

type Content = Vec<u32>;
type ActiveRow = Vec<u32>;
type Envelope = Vec<u32>;
type GtSum = Vec<Vec<u32>>;

#[derive(Debug, Clone)]
struct Args {
    lambda: Option<Vec<u32>>,
    mu: Vec<u32>,
    row_flags: Option<Vec<u32>>,
    alphabet: usize,
    lower_label: u32,
    mode: Mode,
    max_skew_size: u32,
    max_outer_extra: u32,
    connected_only: bool,
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
            mode: Mode::ExactEnvelope,
            max_skew_size: 7,
            max_outer_extra: 8,
            connected_only: false,
            tableau_limit: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    ActiveRow,
    ExactEnvelope,
    FullGtSum,
    FullGtSumEnvelope,
    NonincreaseEnvelope,
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
    let mut negative_pairs = 0u128;
    let mut positive_pairs = 0u128;
    let mut skipped_by_limit = 0usize;

    for skew_size in 0..=args.max_skew_size {
        let max_outer_size = skew_size + args.max_outer_extra;
        for shape in skew_shapes_of_size(skew_size, max_outer_size) {
            if args.connected_only && !shape.is_connected() {
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
            positive_pairs += outcome.positive_pairs;

            if let Some(failure) = outcome.failure {
                println!("FAIL");
                println!("outer={:?}", flagged_shape.shape().outer().parts());
                println!("inner={:?}", flagged_shape.shape().inner().parts());
                println!("flags={:?}", flagged_shape.row_flags());
                println!("{failure}");
                println!(
                    "shapes_checked={shapes_checked}, fibers_checked={fibers_checked}, negative_pairs={negative_pairs}, positive_pairs={positive_pairs}, skipped_by_limit={skipped_by_limit}, elapsed={:.3}s",
                    started.elapsed().as_secs_f64()
                );
                std::process::exit(1);
            }
        }
        println!(
            "skew_size={skew_size}: shapes_checked={shapes_checked}, fibers_checked={fibers_checked}, negative_pairs={negative_pairs}, positive_pairs={positive_pairs}, skipped_by_limit={skipped_by_limit}, elapsed={:.3}s",
            started.elapsed().as_secs_f64()
        );
    }

    println!("PASS");
    println!(
        "shapes_checked={shapes_checked}, fibers_checked={fibers_checked}, negative_pairs={negative_pairs}, positive_pairs={positive_pairs}, skipped_by_limit={skipped_by_limit}, elapsed={:.3}s",
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
                "--mode" => {
                    args.mode = match take_value(&mut iter, &flag)?.as_str() {
                        "active-row" => Mode::ActiveRow,
                        "exact-envelope" => Mode::ExactEnvelope,
                        "full-gt-sum" => Mode::FullGtSum,
                        "full-gt-sum-envelope" => Mode::FullGtSumEnvelope,
                        "nonincrease-envelope" => Mode::NonincreaseEnvelope,
                        other => {
                            return Err(format!(
                                "invalid --mode `{other}`; expected active-row, exact-envelope, full-gt-sum, full-gt-sum-envelope, or nonincrease-envelope"
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
    failure: Option<EnvelopeFailure>,
    limit_exceeded: bool,
}

#[derive(Debug, Clone)]
struct EnvelopeFailure {
    beta: Content,
    negative_pairs: u128,
    positive_pairs: u128,
    matched_pairs: u128,
    mode: Mode,
    witness: Option<Witness>,
}

#[derive(Debug, Clone)]
struct Witness {
    active_row: ActiveRow,
    gt_sum: GtSum,
    envelope: Envelope,
    negative_count: u128,
    positive_count: u128,
}

impl std::fmt::Display for EnvelopeFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "beta={:?}, mode={:?}, negative_pairs={}, positive_pairs={}, matched_pairs={}",
            self.beta, self.mode, self.negative_pairs, self.positive_pairs, self.matched_pairs
        )?;
        if let Some(witness) = &self.witness {
            write!(
                f,
                ", witness_active_row={:?}, witness_gt_sum={:?}, witness_envelope={:?}, witness_negative={}, witness_positive={}",
                witness.active_row,
                witness.gt_sum,
                witness.envelope,
                witness.negative_count,
                witness.positive_count
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
                failure: None,
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
            args.alphabet,
            args.lower_label,
            args.mode,
        );
        negative_pairs += result.negative_pairs;
        positive_pairs += result.positive_pairs;
        if result.matched_pairs < result.negative_pairs {
            return ShapeScanOutcome {
                tableaux: tableaux.len(),
                fibers_checked,
                negative_pairs,
                positive_pairs,
                failure: Some(result),
                limit_exceeded: false,
            };
        }
    }

    ShapeScanOutcome {
        tableaux: tableaux.len(),
        fibers_checked,
        negative_pairs,
        positive_pairs,
        failure: None,
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
    alphabet: usize,
    lower_label: u32,
    mode: Mode,
) -> EnvelopeFailure {
    let negative_pairs = ii_indices.len() as u128 * jj_indices.len() as u128;
    let positive_pairs = mixed_indices.len() as u128 * mixed_indices.len() as u128;

    match mode {
        Mode::ActiveRow | Mode::ExactEnvelope | Mode::FullGtSum | Mode::FullGtSumEnvelope => {
            let negative_counts = exact_counts(
                shape,
                tableaux,
                ii_indices,
                jj_indices,
                alphabet,
                lower_label,
                mode,
            );
            let positive_counts = exact_counts(
                shape,
                tableaux,
                mixed_indices,
                mixed_indices,
                alphabet,
                lower_label,
                mode,
            );
            let mut matched_pairs = 0u128;
            let mut witness = None;
            for (key, &negative_count) in &negative_counts {
                let positive_count = positive_counts.get(key).copied().unwrap_or(0);
                matched_pairs += negative_count.min(positive_count);
                if negative_count > positive_count && witness.is_none() {
                    witness = Some(Witness {
                        active_row: key.active_row.clone(),
                        gt_sum: key.gt_sum.clone(),
                        envelope: key.envelope.clone(),
                        negative_count,
                        positive_count,
                    });
                }
            }
            EnvelopeFailure {
                beta: beta.to_vec(),
                negative_pairs,
                positive_pairs,
                matched_pairs,
                mode,
                witness,
            }
        }
        Mode::NonincreaseEnvelope => {
            let negative_pairs_data = pair_data(
                shape,
                tableaux,
                ii_indices,
                jj_indices,
                alphabet,
                lower_label,
            );
            let positive_pairs_data = pair_data(
                shape,
                tableaux,
                mixed_indices,
                mixed_indices,
                alphabet,
                lower_label,
            );
            let edges = build_nonincrease_edges(&negative_pairs_data, &positive_pairs_data);
            let matched_pairs = max_bipartite_matching(&edges, positive_pairs_data.len()) as u128;
            EnvelopeFailure {
                beta: beta.to_vec(),
                negative_pairs,
                positive_pairs,
                matched_pairs,
                mode,
                witness: None,
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct ExactKey {
    active_row: ActiveRow,
    gt_sum: GtSum,
    envelope: Envelope,
}

#[derive(Debug, Clone)]
struct PairData {
    active_row: ActiveRow,
    gt_sum: GtSum,
    envelope: Envelope,
}

fn exact_counts(
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    left_indices: &[usize],
    right_indices: &[usize],
    alphabet: usize,
    lower_label: u32,
    mode: Mode,
) -> BTreeMap<ExactKey, u128> {
    let mut counts = BTreeMap::new();
    for &left_idx in left_indices {
        for &right_idx in right_indices {
            let data = pair_datum(
                shape,
                &tableaux[left_idx],
                &tableaux[right_idx],
                alphabet,
                lower_label,
            );
            let gt_sum = match mode {
                Mode::ActiveRow | Mode::ExactEnvelope => Vec::new(),
                Mode::FullGtSum | Mode::FullGtSumEnvelope => data.gt_sum,
                Mode::NonincreaseEnvelope => unreachable!(),
            };
            let envelope = match mode {
                Mode::ActiveRow | Mode::FullGtSum => Vec::new(),
                Mode::ExactEnvelope | Mode::FullGtSumEnvelope => data.envelope,
                Mode::NonincreaseEnvelope => unreachable!(),
            };
            *counts
                .entry(ExactKey {
                    active_row: data.active_row,
                    gt_sum,
                    envelope,
                })
                .or_insert(0) += 1;
        }
    }
    counts
}

fn pair_data(
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    left_indices: &[usize],
    right_indices: &[usize],
    alphabet: usize,
    lower_label: u32,
) -> Vec<PairData> {
    let mut out = Vec::new();
    for &left_idx in left_indices {
        for &right_idx in right_indices {
            out.push(pair_datum(
                shape,
                &tableaux[left_idx],
                &tableaux[right_idx],
                alphabet,
                lower_label,
            ));
        }
    }
    out
}

fn pair_datum(
    shape: &SkewShape,
    left: &TableauRecord,
    right: &TableauRecord,
    alphabet: usize,
    lower_label: u32,
) -> PairData {
    let left_gt = gt::SkewGtPattern::from_tableau(shape, &left.values, alphabet);
    let right_gt = gt::SkewGtPattern::from_tableau(shape, &right.values, alphabet);
    let gt_sum = gt::add_patterns(&left_gt, &right_gt);
    PairData {
        active_row: gt_sum[lower_label as usize].clone(),
        gt_sum,
        envelope: gt::pair_envelope(
            &gt::sharp_flag(shape, &left.values),
            &gt::sharp_flag(shape, &right.values),
        ),
    }
}

fn build_nonincrease_edges(negative: &[PairData], positive: &[PairData]) -> Vec<Vec<usize>> {
    negative
        .iter()
        .map(|neg| {
            positive
                .iter()
                .enumerate()
                .filter_map(|(pos_idx, pos)| {
                    (neg.active_row == pos.active_row && less_equal(&pos.envelope, &neg.envelope))
                        .then_some(pos_idx)
                })
                .collect()
        })
        .collect()
}

fn less_equal(left: &[u32], right: &[u32]) -> bool {
    left.iter().zip(right).all(|(&left, &right)| left <= right)
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

fn flagged_shape_for_args(shape: SkewShape, args: &Args) -> RowFlaggedSkewShape {
    match &args.row_flags {
        Some(row_flags) => RowFlaggedSkewShape::new(shape, row_flags.clone(), args.alphabet),
        None => RowFlaggedSkewShape::ordinary(shape, args.alphabet),
    }
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
    println!("flags={:?}", flagged_shape.row_flags());
    println!("alphabet={}", flagged_shape.alphabet());
    println!("tableaux={}", outcome.tableaux);
    println!("fibers_checked={}", outcome.fibers_checked);
    println!("negative_pairs={}", outcome.negative_pairs);
    println!("positive_pairs={}", outcome.positive_pairs);
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
    "Check GT active-row and flag-envelope matching refinements for 2x2 fibers.

USAGE:
  gt_envelope_scan [OPTIONS]

OPTIONS:
  --lambda PARTS          Outer partition, e.g. 6,2. If omitted, scan a family.
  --mu PARTS              Inner partition, e.g. 2. Default: empty.
  --row-flags PARTS       Row upper flags, e.g. 4,5. Omit for ordinary skew tableaux.
  --alphabet N            Alphabet size. Default: 5.
  --lower-label I         Use top-row invariant for labels I,I+1. Default: 4.
  --mode KIND             active-row, exact-envelope, full-gt-sum,
                          full-gt-sum-envelope, or nonincrease-envelope.
                          Default: exact-envelope.
  --max-skew-size N       Maximum skew size for family scans. Default: 7.
  --max-outer-extra N     Family outer sizes go up to skew_size + this. Default: 8.
  --connected-only        Restrict family scans to connected skew shapes.
  --tableau-limit N       Skip any shape whose tableau enumeration exceeds N.
  --help                  Print this help.

EXAMPLE:
  cargo run -q -p flagged-lorentzian --bin gt_envelope_scan -- \\
    --lambda 6,2 --mu 2 --alphabet 5 --lower-label 4 --mode exact-envelope"
}
