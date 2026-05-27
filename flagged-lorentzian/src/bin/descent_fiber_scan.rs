use std::time::Instant;

use flagged_lorentzian::{
    check_two_by_two_fiber_inequalities, enumerate_content_statistic_counts,
    families::skew_shapes_of_size, DescentStatistic, EnumerationOptions, RowFlaggedSkewShape,
    SkewShape,
};

#[derive(Debug, Clone)]
struct Args {
    lambda: Option<Vec<u32>>,
    mu: Vec<u32>,
    row_flags: Option<Vec<u32>>,
    alphabet: usize,
    stat: StatisticArg,
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
            stat: StatisticArg::Componentwise,
            max_skew_size: 7,
            max_outer_extra: 8,
            connected_only: false,
            tableau_limit: None,
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum StatisticArg {
    Global,
    Componentwise,
}

impl From<StatisticArg> for DescentStatistic {
    fn from(value: StatisticArg) -> Self {
        match value {
            StatisticArg::Global => DescentStatistic::Global,
            StatisticArg::Componentwise => DescentStatistic::Componentwise,
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
                "--stat" => {
                    args.stat = match take_value(&mut iter, &flag)?.as_str() {
                        "global" => StatisticArg::Global,
                        "componentwise" => StatisticArg::Componentwise,
                        other => {
                            return Err(format!(
                                "invalid --stat `{other}`; expected global or componentwise"
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

fn take_value(iter: &mut impl Iterator<Item = String>, flag: &str) -> Result<String, String> {
    iter.next()
        .ok_or_else(|| format!("missing value after {flag}"))
}

fn help_text() -> &'static str {
    "Check descent-fiber 2x2 inequalities for flagged skew tableaux.

USAGE:
  descent_fiber_scan [OPTIONS]

OPTIONS:
  --lambda PARTS          Outer partition, e.g. 4,3,1. If omitted, scan a family.
  --mu PARTS              Inner partition, e.g. 3,1. Default: empty.
  --row-flags PARTS       Row upper flags, e.g. 3,4,5. Omit for ordinary skew tableaux.
  --alphabet N            Alphabet size. Default: 5.
  --stat KIND             global or componentwise. Default: componentwise.
  --max-skew-size N       Maximum skew size for family scans. Default: 7.
  --max-outer-extra N     Family outer sizes go up to skew_size + this. Default: 8.
  --connected-only        Restrict family scans to connected skew shapes.
  --tableau-limit N       Skip any shape whose tableau enumeration exceeds N.
  --help                  Print this help.

EXAMPLES:
  cargo run -q -p flagged-lorentzian --bin descent_fiber_scan -- \\
    --lambda 4,3,1 --mu 3,1 --alphabet 5 --stat global

  cargo run -q -p flagged-lorentzian --bin descent_fiber_scan -- \\
    --lambda 4,3,1 --mu 3,1 --alphabet 5 --stat componentwise"
}

#[derive(Debug)]
struct ShapeScanOutcome {
    tests_checked: usize,
    total_tableaux: u128,
    failure: Option<flagged_lorentzian::FiberFailure>,
    limit_exceeded: bool,
}

fn scan_one_shape(flagged_shape: &RowFlaggedSkewShape, args: &Args) -> ShapeScanOutcome {
    let options = EnumerationOptions::new(args.stat.into()).with_tableau_limit(args.tableau_limit);
    match enumerate_content_statistic_counts(flagged_shape, options) {
        Ok(counts) => {
            let total_tableaux = flagged_lorentzian::enumeration::total_tableau_count(&counts);
            let report = check_two_by_two_fiber_inequalities(&counts, flagged_shape.alphabet());
            ShapeScanOutcome {
                tests_checked: report.tests_checked,
                total_tableaux,
                failure: report.failure,
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
