use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use flagged_lorentzian::{
    enumerate_tableaux, families::skew_shapes_of_size, RowFlaggedSkewShape, SkewShape,
};

type Content = Vec<u32>;

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
            max_skew_size: 7,
            max_outer_extra: 8,
            max_rows: None,
            connected_only: false,
            tableau_limit: None,
        }
    }
}

#[derive(Debug, Clone)]
struct ShapeOutcome {
    fibers_checked: usize,
    failure: Option<FiberFailure>,
    limit_exceeded: bool,
}

#[derive(Debug, Clone)]
struct FiberFailure {
    beta: Content,
    ii_count: usize,
    mixed_count: usize,
    jj_count: usize,
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
            shapes_checked += 1;
            fibers_checked += outcome.fibers_checked;
            if let Some(failure) = outcome.failure {
                println!("FAIL");
                println!("outer={:?}", flagged_shape.shape().outer().parts());
                println!("inner={:?}", flagged_shape.shape().inner().parts());
                println!("flags={:?}", flagged_shape.row_flags());
                println!(
                    "beta={:?}, ii={}, mixed={}, jj={}",
                    failure.beta, failure.ii_count, failure.mixed_count, failure.jj_count
                );
                println!(
                    "shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, fibers_checked={fibers_checked}, elapsed={:.3}s",
                    started.elapsed().as_secs_f64()
                );
                std::process::exit(1);
            }
        }
        println!(
            "skew_size={skew_size}: shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, fibers_checked={fibers_checked}, elapsed={:.3}s",
            started.elapsed().as_secs_f64()
        );
    }

    println!("PASS");
    println!(
        "shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, fibers_checked={fibers_checked}, elapsed={:.3}s",
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

fn scan_one_shape(flagged_shape: &RowFlaggedSkewShape, args: &Args) -> ShapeOutcome {
    let tableaux = match enumerate_tableaux(flagged_shape, args.tableau_limit) {
        Ok(tableaux) => tableaux,
        Err(_) => {
            return ShapeOutcome {
                fibers_checked: 0,
                failure: None,
                limit_exceeded: true,
            }
        }
    };

    let mut by_content = BTreeMap::<Content, usize>::new();
    for tableau in tableaux {
        *by_content.entry(tableau.content).or_insert(0) += 1;
    }

    let lower = args.lower_label as usize - 1;
    let upper = lower + 1;
    let mut seen_beta = BTreeSet::new();
    let mut fibers_checked = 0usize;

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
        let (Some(&ii_count), Some(&mixed_count), Some(&jj_count)) = (
            by_content.get(&ii),
            by_content.get(&mixed),
            by_content.get(&jj),
        ) else {
            continue;
        };
        fibers_checked += 1;
        if ii_count > mixed_count || jj_count > mixed_count {
            return ShapeOutcome {
                fibers_checked,
                failure: Some(FiberFailure {
                    beta,
                    ii_count,
                    mixed_count,
                    jj_count,
                }),
                limit_exceeded: false,
            };
        }
    }

    ShapeOutcome {
        fibers_checked,
        failure: None,
        limit_exceeded: false,
    }
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
    println!("elapsed={elapsed_seconds:.3}s");
    if let Some(failure) = &outcome.failure {
        println!("FAIL");
        println!(
            "beta={:?}, ii={}, mixed={}, jj={}",
            failure.beta, failure.ii_count, failure.mixed_count, failure.jj_count
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
    "Check whether both one-tableau fibers inject into the mixed fiber by count.

USAGE:
  one_factor_count_scan [OPTIONS]

OPTIONS:
  --lambda PARTS          Outer partition, e.g. 4,2. If omitted, scan a family.
  --mu PARTS              Inner partition, e.g. 2. Default: empty.
  --row-flags PARTS       Row upper flags. Omit for ordinary skew tableaux.
  --alphabet N            Alphabet size. Default: 5.
  --lower-label I         Use active labels I,I+1. Default: 4.
  --max-skew-size N       Maximum skew size for family scans. Default: 7.
  --max-outer-extra N     Family outer sizes go up to skew_size + this. Default: 8.
  --max-rows N            Restrict family scans to shapes with at most N rows.
  --connected-only        Restrict family scans to connected skew shapes.
  --tableau-limit N       Skip shapes whose tableau enumeration exceeds N.
  --help                  Print this help."
}
