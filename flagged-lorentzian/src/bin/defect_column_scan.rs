use std::time::Instant;

use flagged_lorentzian::{
    check_two_by_two_defect_column_count_fibers, check_two_by_two_defect_column_fibers,
    check_two_by_two_pair_defect_column_count_fibers, check_two_by_two_pair_defect_column_fibers,
    enumerate_tableaux, families::skew_shapes_of_size, RowFlaggedSkewShape, SkewShape,
};

#[derive(Debug, Clone)]
struct Args {
    lambda: Option<Vec<u32>>,
    mu: Vec<u32>,
    row_flags: Option<Vec<u32>>,
    alphabet: usize,
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
            mode: Mode::Set,
            max_skew_size: 7,
            max_outer_extra: 8,
            connected_only: false,
            tableau_limit: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    Set,
    Count,
    PairSet,
    PairCount,
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
    let mut tests_checked = 0usize;
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
            tests_checked += outcome.tests_checked;
            if let Some(failure) = outcome.failure {
                println!("FAIL");
                println!("outer={:?}", flagged_shape.shape().outer().parts());
                println!("inner={:?}", flagged_shape.shape().inner().parts());
                println!("flags={:?}", flagged_shape.row_flags());
                println!("{failure}");
                println!(
                    "shapes_checked={shapes_checked}, tests_checked={tests_checked}, skipped_by_limit={skipped_by_limit}, elapsed={:.3}s",
                    started.elapsed().as_secs_f64()
                );
                std::process::exit(1);
            }
        }
        println!(
            "skew_size={skew_size}: shapes_checked={shapes_checked}, tests_checked={tests_checked}, skipped_by_limit={skipped_by_limit}, elapsed={:.3}s",
            started.elapsed().as_secs_f64()
        );
    }

    println!("PASS");
    println!(
        "shapes_checked={shapes_checked}, tests_checked={tests_checked}, skipped_by_limit={skipped_by_limit}, elapsed={:.3}s",
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
                "--mode" => {
                    args.mode = match take_value(&mut iter, &flag)?.as_str() {
                        "set" => Mode::Set,
                        "count" => Mode::Count,
                        "pair-set" => Mode::PairSet,
                        "pair-count" => Mode::PairCount,
                        other => {
                            return Err(format!(
                        "invalid --mode `{other}`; expected set, count, pair-set, or pair-count"
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
        Ok(args)
    }
}

#[derive(Debug)]
struct ShapeScanOutcome {
    tests_checked: usize,
    total_tableaux: usize,
    failure: Option<String>,
    limit_exceeded: bool,
}

fn scan_one_shape(flagged_shape: &RowFlaggedSkewShape, args: &Args) -> ShapeScanOutcome {
    match enumerate_tableaux(flagged_shape, args.tableau_limit) {
        Ok(tableaux) => {
            let (tests_checked, failure) = match args.mode {
                Mode::Set => {
                    let report = check_two_by_two_defect_column_fibers(
                        flagged_shape.shape(),
                        &tableaux,
                        flagged_shape.alphabet(),
                    );
                    (
                        report.tests_checked,
                        report.failure.map(|failure| failure.to_string()),
                    )
                }
                Mode::Count => {
                    let report = check_two_by_two_defect_column_count_fibers(
                        flagged_shape.shape(),
                        &tableaux,
                        flagged_shape.alphabet(),
                    );
                    (
                        report.tests_checked,
                        report.failure.map(|failure| failure.to_string()),
                    )
                }
                Mode::PairSet => {
                    let report = check_two_by_two_pair_defect_column_fibers(
                        flagged_shape.shape(),
                        &tableaux,
                        flagged_shape.alphabet(),
                    );
                    (
                        report.tests_checked,
                        report.failure.map(|failure| failure.to_string()),
                    )
                }
                Mode::PairCount => {
                    let report = check_two_by_two_pair_defect_column_count_fibers(
                        flagged_shape.shape(),
                        &tableaux,
                        flagged_shape.alphabet(),
                    );
                    (
                        report.tests_checked,
                        report.failure.map(|failure| failure.to_string()),
                    )
                }
            };
            ShapeScanOutcome {
                tests_checked,
                total_tableaux: tableaux.len(),
                failure,
                limit_exceeded: false,
            }
        }
        Err(_) => ShapeScanOutcome {
            tests_checked: 0,
            total_tableaux: 0,
            failure: None,
            limit_exceeded: true,
        },
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
    println!("tableaux={}", outcome.total_tableaux);
    println!("tests_checked={}", outcome.tests_checked);
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
    "Check 2x2 inequalities after fixing columns that contain both defect labels.

USAGE:
  defect_column_scan [OPTIONS]

OPTIONS:
  --lambda PARTS          Outer partition, e.g. 4,3,1. If omitted, scan a family.
  --mu PARTS              Inner partition, e.g. 3,1. Default: empty.
  --row-flags PARTS       Row upper flags, e.g. 3,4,5. Omit for ordinary skew tableaux.
  --alphabet N            Alphabet size. Default: 5.
  --mode KIND             set, count, pair-set, or pair-count. Default: set.
  --max-skew-size N       Maximum skew size for family scans. Default: 7.
  --max-outer-extra N     Family outer sizes go up to skew_size + this. Default: 8.
  --connected-only        Restrict family scans to connected skew shapes.
  --tableau-limit N       Skip any shape whose tableau enumeration exceeds N.
  --help                  Print this help.

EXAMPLE:
  cargo run -q -p flagged-lorentzian --bin defect_column_scan -- \\
    --lambda 4,3,1 --mu 3,1 --alphabet 5"
}
