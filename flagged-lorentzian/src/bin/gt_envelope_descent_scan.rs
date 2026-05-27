use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use flagged_lorentzian::descent::{
    active_subword_descent_data_for_values, descent_data_for_values, DescentData, DescentStatistic,
};
use flagged_lorentzian::{
    enumerate_tableaux, families::skew_shapes_of_size, RowFlaggedSkewShape, SkewShape,
    TableauRecord,
};

type Content = Vec<u32>;
type ActiveRow = Vec<u32>;
type Envelope = Vec<u32>;

#[derive(Debug, Clone)]
struct Args {
    lambda: Option<Vec<u32>>,
    mu: Vec<u32>,
    row_flags: Option<Vec<u32>>,
    alphabet: usize,
    lower_label: u32,
    descent_mode: DescentMode,
    exact_envelope: bool,
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
            descent_mode: DescentMode::Componentwise,
            exact_envelope: true,
            max_skew_size: 7,
            max_outer_extra: 8,
            connected_only: false,
            tableau_limit: None,
        }
    }
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
                "--descent" => {
                    args.descent_mode = match take_value(&mut iter, &flag)?.as_str() {
                        "global" => DescentMode::Global,
                        "componentwise" => DescentMode::Componentwise,
                        "active-global" => DescentMode::ActiveGlobal,
                        "active-componentwise" => DescentMode::ActiveComponentwise,
                        other => {
                            return Err(format!(
                                "invalid --descent `{other}`; expected global, componentwise, active-global, or active-componentwise"
                            ))
                        }
                    }
                }
                "--active-row-only" => args.exact_envelope = false,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DescentMode {
    Global,
    Componentwise,
    ActiveGlobal,
    ActiveComponentwise,
}

impl DescentMode {
    fn reading_orders(self, shape: &SkewShape) -> Vec<Vec<usize>> {
        match self {
            DescentMode::Global | DescentMode::ActiveGlobal => {
                DescentStatistic::Global.reading_orders(shape)
            }
            DescentMode::Componentwise | DescentMode::ActiveComponentwise => {
                DescentStatistic::Componentwise.reading_orders(shape)
            }
        }
    }

    fn descent_data(
        self,
        values: &[u32],
        reading_orders: &[Vec<usize>],
        lower_label: u32,
    ) -> DescentData {
        match self {
            DescentMode::Global | DescentMode::Componentwise => {
                descent_data_for_values(values, reading_orders)
            }
            DescentMode::ActiveGlobal | DescentMode::ActiveComponentwise => {
                active_subword_descent_data_for_values(values, reading_orders, lower_label)
            }
        }
    }
}

#[derive(Debug)]
struct ShapeScanOutcome {
    tableaux: usize,
    fibers_checked: usize,
    negative_pairs: u128,
    positive_pairs: u128,
    failure: Option<CombinedFailure>,
    limit_exceeded: bool,
}

#[derive(Debug, Clone)]
struct CombinedFailure {
    beta: Content,
    key: CombinedKey,
    negative_count: u128,
    positive_count: u128,
}

impl std::fmt::Display for CombinedFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "beta={:?}, active_row={:?}, envelope={:?}, descents=({},{}), negative={}, positive={}",
            self.beta,
            self.key.active_row,
            self.key.envelope,
            self.key.first_descent,
            self.key.second_descent,
            self.negative_count,
            self.positive_count
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct CombinedKey {
    active_row: ActiveRow,
    envelope: Envelope,
    first_descent: DescentData,
    second_descent: DescentData,
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

    let reading_orders = args.descent_mode.reading_orders(flagged_shape.shape());
    let descents: Vec<_> = tableaux
        .iter()
        .map(|tableau| {
            args.descent_mode
                .descent_data(&tableau.values, &reading_orders, args.lower_label)
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
        negative_pairs += ii_indices.len() as u128 * jj_indices.len() as u128;
        positive_pairs += mixed_indices.len() as u128 * mixed_indices.len() as u128;

        if let Some(failure) = check_fiber(
            flagged_shape.shape(),
            &tableaux,
            &descents,
            &beta,
            ii_indices,
            mixed_indices,
            jj_indices,
            args.lower_label,
            args.exact_envelope,
        ) {
            return ShapeScanOutcome {
                tableaux: tableaux.len(),
                fibers_checked,
                negative_pairs,
                positive_pairs,
                failure: Some(failure),
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

#[allow(clippy::too_many_arguments)]
fn check_fiber(
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    descents: &[DescentData],
    beta: &[u32],
    ii_indices: &[usize],
    mixed_indices: &[usize],
    jj_indices: &[usize],
    lower_label: u32,
    exact_envelope: bool,
) -> Option<CombinedFailure> {
    let negative_counts = combined_counts(
        shape,
        tableaux,
        descents,
        ii_indices,
        jj_indices,
        lower_label,
        exact_envelope,
    );
    let positive_counts = combined_counts(
        shape,
        tableaux,
        descents,
        mixed_indices,
        mixed_indices,
        lower_label,
        exact_envelope,
    );

    for (key, &negative_count) in &negative_counts {
        let positive_count = positive_counts.get(key).copied().unwrap_or(0);
        if negative_count > positive_count {
            return Some(CombinedFailure {
                beta: beta.to_vec(),
                key: key.clone(),
                negative_count,
                positive_count,
            });
        }
    }

    None
}

fn combined_counts(
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    descents: &[DescentData],
    left_indices: &[usize],
    right_indices: &[usize],
    lower_label: u32,
    exact_envelope: bool,
) -> BTreeMap<CombinedKey, u128> {
    let mut counts = BTreeMap::new();
    for &left_idx in left_indices {
        for &right_idx in right_indices {
            let mut first_descent = descents[left_idx].clone();
            let mut second_descent = descents[right_idx].clone();
            if second_descent < first_descent {
                std::mem::swap(&mut first_descent, &mut second_descent);
            }

            let envelope = if exact_envelope {
                envelope(
                    &sharp_flag(shape, &tableaux[left_idx].values),
                    &sharp_flag(shape, &tableaux[right_idx].values),
                )
            } else {
                Vec::new()
            };

            *counts
                .entry(CombinedKey {
                    active_row: add_vectors(
                        &gt_active_row(shape, &tableaux[left_idx].values, lower_label),
                        &gt_active_row(shape, &tableaux[right_idx].values, lower_label),
                    ),
                    envelope,
                    first_descent,
                    second_descent,
                })
                .or_insert(0) += 1;
        }
    }
    counts
}

fn gt_active_row(shape: &SkewShape, values: &[u32], lower_label: u32) -> ActiveRow {
    let mut row = (0..shape.row_count())
        .map(|idx| shape.inner().part(idx))
        .collect::<Vec<_>>();
    for (cell, &value) in shape.cells().iter().zip(values) {
        if value <= lower_label {
            row[cell.row] += 1;
        }
    }
    row
}

fn sharp_flag(shape: &SkewShape, values: &[u32]) -> Envelope {
    let mut raw = vec![0u32; shape.row_count()];
    for (cell, &value) in shape.cells().iter().zip(values) {
        raw[cell.row] = raw[cell.row].max(value);
    }

    let mut sharp = raw;
    for row in 1..sharp.len() {
        sharp[row] = sharp[row].max(sharp[row - 1]);
    }
    sharp
}

fn envelope(left: &[u32], right: &[u32]) -> Envelope {
    left.iter()
        .zip(right)
        .map(|(&left, &right)| left.max(right))
        .collect()
}

fn add_vectors(left: &[u32], right: &[u32]) -> ActiveRow {
    left.iter()
        .zip(right)
        .map(|(&left, &right)| left + right)
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
    "Check whether descent data remains preservable with the GT active-row/envelope invariant.

USAGE:
  gt_envelope_descent_scan [OPTIONS]

OPTIONS:
  --lambda PARTS          Outer partition, e.g. 6,2. If omitted, scan a family.
  --mu PARTS              Inner partition, e.g. 2. Default: empty.
  --row-flags PARTS       Row upper flags, e.g. 4,5. Omit for ordinary skew tableaux.
  --alphabet N            Alphabet size. Default: 5.
  --lower-label I         Use active GT row for labels I,I+1. Default: 4.
  --descent KIND          global, componentwise, active-global, or active-componentwise.
                          Default: componentwise.
  --active-row-only       Do not include exact pair envelope in the key.
  --max-skew-size N       Maximum skew size for family scans. Default: 7.
  --max-outer-extra N     Family outer sizes go up to skew_size + this. Default: 8.
  --connected-only        Restrict family scans to connected skew shapes.
  --tableau-limit N       Skip any shape whose tableau enumeration exceeds N.
  --help                  Print this help.

EXAMPLE:
  cargo run -q -p flagged-lorentzian --bin gt_envelope_descent_scan -- \\
    --lambda 6,2 --mu 2 --alphabet 5 --lower-label 4 --descent componentwise"
}
