//! CLI for polynomial-tools: real-rootedness, interlacing, recurrence search.
//!
//! Reads polynomials from stdin as comma-separated integer coefficients
//! in ascending degree order, one polynomial per line.

use num_bigint::BigInt;
use polynomial_tools::recurrence::BigRational as RecurrenceBigRational;
use polynomial_tools::recurrence::*;
use polynomial_tools::*;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::hint::black_box;
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::Instant;

fn is_help_arg(arg: &str) -> bool {
    matches!(arg, "-h" | "--help" | "help")
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum OutputFormat {
    Text,
    Json,
}

impl OutputFormat {
    fn from_args(args: &[String]) -> Result<Self, String> {
        let mut format = Self::Text;
        for arg in args {
            match arg.as_str() {
                "--json" => format = Self::Json,
                other => return Err(format!("unknown option: {other}")),
            }
        }
        Ok(format)
    }
}

fn json_escape(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    for ch in input.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

fn json_string(input: &str) -> String {
    format!("\"{}\"", json_escape(input))
}

fn json_bigint_vec(values: &[BigInt]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|value| json_string(&value.to_string()))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn json_usize_vec(values: &[usize]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn json_string_vec(values: &[String]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|value| json_string(value))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn json_bool_option(value: Option<bool>) -> String {
    value
        .map(|v| if v { "true" } else { "false" }.to_string())
        .unwrap_or_else(|| "null".to_string())
}

fn print_coefficient_input_help() {
    println!("Input:");
    println!("  Read one polynomial per line from stdin.");
    println!("  Dense coefficient lists are in ascending degree order:");
    println!("    a_0, a_1, ..., a_d  represents  a_0 + a_1 t + ... + a_d t^d");
    println!("  Brackets, whitespace-separated coefficients, and expanded notation");
    println!("  such as 1 + 2t + 3t^2 are also accepted.");
    println!("  Blank lines and lines starting with # are ignored.");
}

fn print_rational_coefficient_input_help() {
    println!("Input:");
    println!("  Read dense coefficient lists from stdin, one polynomial per line.");
    println!("  Coefficients are in ascending degree order and may be arbitrary-size");
    println!("  integers or exact rationals such as 42, -17, or 3/7.");
    println!("  Brackets and whitespace-separated coefficients are accepted.");
    println!("  Expanded polynomial expressions are not accepted by recurrence.");
    println!("  Blank lines and lines starting with # are ignored.");
}

fn print_top_level_help() {
    println!("polytool {}", env!("CARGO_PKG_VERSION"));
    println!("Dense univariate polynomial tools for combinatorial research.");
    println!();
    println!("Usage:");
    println!("  polytool <command> [options]");
    println!("  polytool --help");
    println!("  polytool help <command>");
    println!();
    println!("Commands:");
    println!("  real-rooted       Check real-rootedness of each polynomial");
    println!("  interlacing       Check interlacing of consecutive polynomial pairs");
    println!("  interlacing-profile");
    println!("                    Count consecutive previous interlacings until first fail");
    println!("  properties        Show real-rootedness, unimodality, and related properties");
    println!("  gamma-expansion   Expand palindromic polynomials in the gamma basis");
    println!("  family-check      Check properties and consecutive interlacing together");
    println!("  sequence          Generate standard polynomial sequences");
    println!("  recurrence        Search for a polynomial recurrence");
    println!("  recurrence-generate");
    println!("                    Generate coefficient rows from recurrence JSON");
    println!("  bench             Run built-in timing suites");
    println!("  bkw-scout         Scout BKW equal-modulus loci for a recurrence symbol");
    println!("  resultant         Compute the resultant of two polynomials");
    println!("  discriminant      Compute the discriminant of each polynomial");
    println!("  hstar-to-ehrhart  Convert h*-vectors to Ehrhart polynomials");
    println!("  ehrhart-to-hstar  Convert Ehrhart polynomials to h*-vectors");
    println!("  stapledon         Compute a Stapledon decomposition");
    println!();
    println!("Options:");
    println!("  -h, --help        Print help text");
    println!();
    println!("Run `polytool help <command>` for command-specific help.");
}

fn print_stdin_command_help(command: &str, summary: &str, extra: &[&str]) {
    println!("Usage:");
    println!("  polytool {command}");
    println!("  polytool {command} --help");
    println!();
    println!("{summary}");
    println!();
    print_coefficient_input_help();
    if !extra.is_empty() {
        println!();
        for line in extra {
            println!("{line}");
        }
    }
}

fn print_recurrence_help() {
    println!("Usage:");
    println!("  polytool recurrence [options]");
    println!("  polytool recurrence --help");
    println!();
    println!("Search adaptively for a linear differential-polynomial recurrence");
    println!("among the input polynomial sequence.");
    println!();
    print_rational_coefficient_input_help();
    println!();
    println!("Search options:");
    println!("  --skip-prefix <n>          Ignore the first n input polynomials");
    println!("  --min-rec-len <n>          Minimum number of previous rows to use");
    println!("  --max-rec-len <n>          Maximum number of previous rows to use");
    println!("  --min-var-deg <d>          Minimum t-degree for recurrence coefficients");
    println!("  --max-var-deg <d>          Maximum t-degree for recurrence coefficients");
    println!("  --min-idx-deg <d>          Minimum n-degree for recurrence coefficients");
    println!("  --max-idx-deg <d>          Maximum n-degree for recurrence coefficients");
    println!("  --min-diff-deg <d>         Minimum derivative order");
    println!("  --max-diff-deg <d>         Maximum derivative order");
    println!("  --inhomogeneous            Try inhomogeneous terms");
    println!("  --min-inhomo-var-deg <d>   Minimum t-degree for inhomogeneous terms");
    println!("  --max-inhomo-var-deg <d>   Maximum t-degree for inhomogeneous terms");
    println!("  --min-inhomo-idx-deg <d>   Minimum n-degree for inhomogeneous terms");
    println!("  --max-inhomo-idx-deg <d>   Maximum n-degree for inhomogeneous terms");
    println!("  --denominator              Try LHS denominators");
    println!("  --try-denominator          Alias for --denominator");
    println!("  --max-denom-var-deg <d>    Maximum t-degree for denominators");
    println!("  --max-denom-idx-deg <d>    Maximum n-degree for denominators");
    println!("  --alternating-sign         Also try alternating-sign terms");
    println!("  --min-margin <n>           Require equations >= unknowns + n");
    println!("  --fit-extra-rows <n>       Extra rows beyond the first solvable prefix");
    println!("  --no-verify                Use all input rows for fitting");
    println!("  --modular-prefilter        Probabilistically reject candidates modulo primes");
    println!("  --json                     Emit recurrence JSON with initial conditions");
    println!("  --format json              Alias for --json");
    println!("  --python                   Emit exact Python code for the recurrence");
    println!("  --format python            Alias for --python");
    println!("  --verbose                  Print candidate search details");
    println!("  -h, --help                 Print this help text");
    println!();
    println!("Example:");
    println!("  printf '1\\n1\\n2\\n3\\n5\\n8\\n' | polytool recurrence");
    println!("  polytool recurrence --json < rows.txt > recurrence.json");
}

fn print_recurrence_generate_help() {
    println!("Usage:");
    println!("  polytool recurrence-generate --recurrence <file|-> --rows <n>");
    println!("  polytool recurrence-generate --recurrence <file|-> --additional <n>");
    println!("  polytool recurrence-generate --help");
    println!();
    println!("Generate coefficient rows from recurrence JSON emitted by");
    println!("`polytool recurrence --json`.");
    println!();
    println!("Options:");
    println!("  --recurrence <file|->      Read recurrence JSON from file, or stdin with -");
    println!("  --rows <n>                 Output n total rows, including initial rows");
    println!("  --additional <n>           Output initial rows plus n generated rows");
    println!("  --json                     Emit generated rows as JSON");
    println!("  -h, --help                 Print this help text");
    println!();
    println!("Example:");
    println!("  polytool recurrence-generate --recurrence recurrence.json --rows 50");
}

fn next_option_value<'a>(
    args: &'a [String],
    index: &mut usize,
    option: &str,
) -> Result<&'a str, String> {
    *index += 1;
    args.get(*index)
        .map(String::as_str)
        .ok_or_else(|| format!("{option} expects a value"))
}

fn parse_usize_option(args: &[String], index: &mut usize, option: &str) -> Result<usize, String> {
    let value = next_option_value(args, index, option)?;
    value
        .parse()
        .map_err(|_| format!("{option} expects a nonnegative integer, got '{value}'"))
}

fn print_sequence_help() {
    println!("Usage:");
    println!("  polytool sequence <name> <max-n> [--json]");
    println!("  polytool sequence --help");
    println!();
    println!("Generate standard polynomial sequences.");
    println!();
    println!("Names:");
    println!("  eulerian");
    println!("  narayana");
    println!("  type-b-eulerian  (aliases: type-b, type_b_eulerian)");
    println!("  chebyshev-t      (aliases: chebyshev_t, t)");
    println!("  chebyshev-u      (aliases: chebyshev_u, u)");
    println!("  hermite");
    println!();
    println!("Options:");
    println!("  --json           Emit machine-readable JSON");
    println!("  -h, --help       Print this help text");
    println!();
    println!("Note: generated coefficients use arbitrary-size integers.");
    println!();
    println!("Example:");
    println!("  polytool sequence eulerian 5");
}

fn print_family_check_help() {
    println!("Usage:");
    println!("  polytool family-check [options]");
    println!("  polytool family-check --help");
    println!();
    println!("Check a polynomial family in one pass: properties, consecutive");
    println!("interlacing, and optionally recurrence search.");
    println!();
    print_coefficient_input_help();
    println!();
    println!("Options:");
    println!("  --recurrence                  Search for an adaptive recurrence");
    println!("  --require-real-rooted         Fail if any row is not real-rooted");
    println!("  --require-palindromic         Fail if any row is not palindromic");
    println!("  --require-gamma-positive      Fail if any row is not gamma-positive");
    println!("  --require-unimodal            Fail if any row is not unimodal");
    println!("  --require-log-concave         Fail if any row is not log-concave");
    println!("  --require-ultra-log-concave   Fail if any row is not ultra-log-concave");
    println!("  --require-weak-interlacing    Fail if consecutive rows do not weakly interlace");
    println!("  --json                        Emit machine-readable JSON");
    println!("  -h, --help                    Print this help text");
    println!();
    println!("Note: non-recurrence checks accept arbitrary-size integer coefficients;");
    println!("      family-check --recurrence currently requires coefficients that fit in i64.");
}

fn print_stapledon_help() {
    println!("Usage:");
    println!("  polytool stapledon <n>");
    println!("  polytool stapledon --help");
    println!();
    println!("Compute the Stapledon decomposition with respect to degree bound n.");
    println!();
    print_coefficient_input_help();
    println!();
    println!("Note: accepts arbitrary-size integer coefficients.");
}

fn print_command_help(command: &str) -> bool {
    match command {
        "real-rooted" => print_stdin_command_help(
            "real-rooted",
            "Check whether each input polynomial has only real roots.",
            &[
                "Note:",
                "  Accepts arbitrary-size integer coefficients.",
                "",
                "Example:",
                "  echo '1, 3, 2' | polytool real-rooted",
            ],
        ),
        "interlacing" => print_stdin_command_help(
            "interlacing",
            "Check strict and weak interlacing for consecutive input pairs.",
            &[
                "Options:",
                "  --json    Emit machine-readable JSON",
                "",
                "Note:",
                "  Accepts arbitrary-size integer coefficients.",
                "",
                "Example:",
                "  printf '2,-3,1\\n-1,1\\n' | polytool interlacing",
            ],
        ),
        "interlacing-profile" => print_stdin_command_help(
            "interlacing-profile",
            "For each row, count backward consecutive previous interlacings until first failure.",
            &[
                "Options:",
                "  --json    Emit machine-readable JSON, including checked previous-pair reports",
                "",
                "Note:",
                "  Accepts arbitrary-size integer coefficients.",
                "",
                "Example:",
                "  printf '1\\n1,1\\n1,2,1\\n' | polytool interlacing-profile",
            ],
        ),
        "properties" => print_stdin_command_help(
            "properties",
            concat!(
                "Report real-rootedness, palindromicity, gamma-positivity, ",
                "unimodality, log-concavity, and ultra-log-concavity."
            ),
            &[
                "Options:",
                "  --json    Emit machine-readable JSON",
                "",
                "Note:",
                "  Accepts arbitrary-size integer coefficients.",
                "",
                "Example:",
                "  echo '1, 11, 11, 1' | polytool properties",
            ],
        ),
        "gamma-expansion" | "gamma" => print_stdin_command_help(
            "gamma-expansion",
            concat!(
                "Expand each palindromic input polynomial in the gamma basis ",
                "t^i (1+t)^(d-2i)."
            ),
            &[
                "Options:",
                "  --json    Emit machine-readable JSON",
                "",
                "Note:",
                "  Accepts arbitrary-size integer coefficients.",
                "",
                "Aliases:",
                "  polytool gamma",
                "",
                "Example:",
                "  echo '1, 11, 11, 1' | polytool gamma-expansion",
            ],
        ),
        "family-check" => print_family_check_help(),
        "sequence" => print_sequence_help(),
        "recurrence" => print_recurrence_help(),
        "recurrence-generate" => print_recurrence_generate_help(),
        "bench" => print_bench_help(),
        "bkw-scout" => bkw_scout_usage(),
        "resultant" => print_stdin_command_help(
            "resultant",
            "Compute the resultant of the first two input polynomials.",
            &[
                "Note:",
                "  Accepts arbitrary-size integer coefficients.",
                "",
                "Example:",
                "  printf '1,0,1\\n-1,1\\n' | polytool resultant",
            ],
        ),
        "discriminant" => print_stdin_command_help(
            "discriminant",
            "Compute the discriminant of each input polynomial.",
            &[
                "Note:",
                "  Accepts arbitrary-size integer coefficients.",
                "",
                "Example:",
                "  echo '1, 0, 1' | polytool discriminant",
            ],
        ),
        "hstar-to-ehrhart" => print_stdin_command_help(
            "hstar-to-ehrhart",
            "Convert each h*-vector into Ehrhart polynomial coefficients.",
            &[
                "Options:",
                "  --json    Emit machine-readable JSON",
                "",
                "Note:",
                "  Accepts arbitrary-size integer h*-vector entries.",
                "",
                "Example:",
                "  echo '1, 2, 1' | polytool hstar-to-ehrhart",
            ],
        ),
        "ehrhart-to-hstar" => print_stdin_command_help(
            "ehrhart-to-hstar",
            "Convert each Ehrhart polynomial into an h*-vector.",
            &[
                "Input coefficients may be integers or exact rationals.",
                "",
                "Options:",
                "  --json    Emit machine-readable JSON",
                "",
                "Note:",
                "  Returns arbitrary-size integer h*-vector entries.",
                "",
                "Example:",
                "  echo '1, 2, 2' | polytool ehrhart-to-hstar",
            ],
        ),
        "stapledon" => print_stapledon_help(),
        _ => return false,
    }
    true
}

fn read_polys_bigint() -> Vec<Vec<BigInt>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input).unwrap();
    parse_polynomials_bigint(&input)
        .into_iter()
        .filter_map(|r| match r {
            Ok(p) => Some(p),
            Err(e) => {
                eprintln!("Warning: {}", e);
                None
            }
        })
        .collect()
}

fn read_poly_parse_results_bigint() -> Vec<Result<Vec<BigInt>, String>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input).unwrap();
    parse_polynomials_bigint(&input)
}

fn read_polys_rational() -> Vec<Vec<RecurrenceBigRational>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input).unwrap();
    input
        .lines()
        .filter(|l| {
            let t = l.trim();
            !t.is_empty() && !t.starts_with('#')
        })
        .filter_map(|line| match parse_coeff_list_rational(line) {
            Ok(p) => Some(p),
            Err(e) => {
                eprintln!("Warning: {}", e);
                None
            }
        })
        .collect()
}

fn parse_coeff_list_rational(input: &str) -> Result<Vec<RecurrenceBigRational>, String> {
    let s = input.trim();
    if s.is_empty() {
        return Err("empty input".to_string());
    }
    if s.chars().any(|c| c.is_ascii_alphabetic()) {
        return Err(format!(
            "recurrence input expects coefficient lists, got expanded polynomial '{}'",
            input
        ));
    }

    let s = strip_coeff_list_brackets(s);
    let parts: Vec<&str> = if s.contains(',') {
        s.split(',').collect()
    } else {
        s.split_whitespace().collect()
    };

    parts
        .iter()
        .map(|part| parse_big_rational(part.trim()))
        .collect()
}

fn parse_big_rational(token: &str) -> Result<RecurrenceBigRational, String> {
    if token.is_empty() {
        return Err("empty coefficient".to_string());
    }
    if let Some((num, den)) = token.split_once('/') {
        let num = parse_bigint(num.trim())?;
        let den = parse_bigint(den.trim())?;
        if den == BigInt::from(0) {
            return Err(format!("zero denominator in coefficient '{}'", token));
        }
        Ok(RecurrenceBigRational::new(num, den))
    } else {
        Ok(RecurrenceBigRational::from_integer(parse_bigint(token)?))
    }
}

fn parse_bigint(token: &str) -> Result<BigInt, String> {
    token
        .parse::<BigInt>()
        .map_err(|e| format!("invalid integer '{}': {}", token, e))
}

fn strip_coeff_list_brackets(s: &str) -> &str {
    let s = s.trim();
    if (s.starts_with('[') && s.ends_with(']'))
        || (s.starts_with('{') && s.ends_with('}'))
        || (s.starts_with('(') && s.ends_with(')'))
    {
        &s[1..s.len() - 1]
    } else {
        s
    }
}

fn strip_trailing_zeros_bigint(coeffs: &[BigInt]) -> &[BigInt] {
    let end = coeffs
        .iter()
        .rposition(|c| c != &BigInt::from(0))
        .map_or(0, |i| i + 1);
    &coeffs[..end]
}

fn polynomial_degree(coeffs: &[i64]) -> usize {
    coeffs.iter().rposition(|&c| c != 0).unwrap_or(0)
}

fn polynomial_degree_bigint(coeffs: &[BigInt]) -> usize {
    coeffs
        .iter()
        .rposition(|c| c != &BigInt::from(0))
        .unwrap_or(0)
}

#[derive(Clone, Debug)]
struct PropertyReport {
    index: usize,
    coefficients: Vec<BigInt>,
    polynomial: String,
    degree: usize,
    real_rooted: bool,
    simple_roots: bool,
    palindromic: bool,
    gamma_positive: bool,
    gamma_coefficients: Option<Vec<BigInt>>,
    unimodal: bool,
    log_concave: bool,
    ultra_log_concave: bool,
}

fn property_report_i64(index: usize, coeffs: &[i64]) -> PropertyReport {
    let gamma_coefficients = gamma_coefficients_ignoring_initial_zeros(coeffs)
        .map(|gamma| gamma.into_iter().map(BigInt::from).collect::<Vec<_>>());
    let gamma_positive = gamma_coefficients
        .as_ref()
        .is_some_and(|gamma| gamma.iter().all(|g| g >= &BigInt::from(0)));
    PropertyReport {
        index,
        coefficients: coeffs.iter().map(|&c| BigInt::from(c)).collect(),
        polynomial: format_poly(coeffs),
        degree: polynomial_degree(coeffs),
        real_rooted: is_real_rooted(coeffs),
        simple_roots: has_simple_roots(coeffs),
        palindromic: is_palindromic_ignoring_initial_zeros(coeffs),
        gamma_positive,
        gamma_coefficients,
        unimodal: is_unimodal(coeffs),
        log_concave: is_log_concave(coeffs),
        ultra_log_concave: is_ultra_log_concave(coeffs),
    }
}

fn property_report(index: usize, coeffs: &[BigInt]) -> PropertyReport {
    if let Some(coeffs_i64) = bigint_coeffs_to_i64(coeffs) {
        return property_report_i64(index, &coeffs_i64);
    }

    let gamma_coefficients = gamma_coefficients_ignoring_initial_zeros_bigint_coeffs(coeffs);
    let gamma_positive = gamma_coefficients
        .as_ref()
        .is_some_and(|gamma| gamma.iter().all(|g| g >= &BigInt::from(0)));
    PropertyReport {
        index,
        coefficients: coeffs.to_vec(),
        polynomial: format_poly_bigint_coeffs(coeffs),
        degree: polynomial_degree_bigint(coeffs),
        real_rooted: is_real_rooted_bigint_coeffs(coeffs),
        simple_roots: has_simple_roots_bigint_coeffs(coeffs),
        palindromic: is_palindromic_ignoring_initial_zeros_bigint_coeffs(coeffs),
        gamma_positive,
        gamma_coefficients,
        unimodal: is_unimodal_bigint_coeffs(coeffs),
        log_concave: is_log_concave_bigint_coeffs(coeffs),
        ultra_log_concave: is_ultra_log_concave_bigint_coeffs(coeffs),
    }
}

fn property_labels(report: &PropertyReport) -> Vec<String> {
    let mut props = Vec::new();
    if report.real_rooted {
        props.push("real-rooted".to_string());
    }
    if report.palindromic {
        props.push("palindromic".to_string());
    }
    if report.gamma_positive {
        if let Some(gamma) = &report.gamma_coefficients {
            props.push(format!("gamma-positive {:?}", gamma));
        }
    }
    if report.unimodal {
        props.push("unimodal".to_string());
    }
    if report.log_concave {
        props.push("log-concave".to_string());
    }
    if report.ultra_log_concave {
        props.push("ultra-log-concave".to_string());
    }
    if props.is_empty() {
        props.push("(none)".to_string());
    }
    props
}

fn property_report_json(report: &PropertyReport) -> String {
    let gamma = report
        .gamma_coefficients
        .as_ref()
        .map(|coeffs| json_bigint_vec(coeffs))
        .unwrap_or_else(|| "null".to_string());
    format!(
        "{{\"index\":{},\"polynomial\":{},\"coefficients\":{},\"degree\":{},\
         \"real_rooted\":{},\"simple_roots\":{},\"palindromic\":{},\
         \"gamma_positive\":{},\"gamma_coefficients\":{},\"unimodal\":{},\
         \"log_concave\":{},\"ultra_log_concave\":{}}}",
        report.index,
        json_string(&report.polynomial),
        json_bigint_vec(&report.coefficients),
        report.degree,
        report.real_rooted,
        report.simple_roots,
        report.palindromic,
        report.gamma_positive,
        gamma,
        report.unimodal,
        report.log_concave,
        report.ultra_log_concave
    )
}

fn print_property_reports_json(reports: &[PropertyReport]) {
    println!(
        "{{\"items\":[{}]}}",
        reports
            .iter()
            .map(property_report_json)
            .collect::<Vec<_>>()
            .join(",")
    );
}

fn cmd_real_rooted() {
    for coeffs in read_polys_bigint() {
        let c = strip_trailing_zeros_bigint(&coeffs);
        let rr = is_real_rooted_bigint_coeffs(c);
        println!(
            "{}: {}",
            format_poly_bigint_coeffs(c),
            if rr { "real-rooted" } else { "NOT real-rooted" }
        );
    }
}

#[derive(Clone, Debug)]
struct InterlacingReport {
    pair_index: usize,
    left_index: usize,
    right_index: usize,
    p: String,
    q: String,
    strict: Option<bool>,
    weak: Option<bool>,
    status: String,
}

fn interlacing_status(strict: Option<bool>, weak: Option<bool>) -> String {
    match (strict, weak) {
        (Some(true), _) => "strictly_interlace",
        (_, Some(true)) => "weakly_interlace",
        (Some(false), Some(false)) => "do_not_interlace",
        (Some(false), None) => "not_real_rooted_or_incompatible",
        (None, Some(false)) => "not_real_rooted_or_incompatible",
        (None, None) => "not_real_rooted_or_incompatible",
    }
    .to_string()
}

fn interlacing_report(
    pair_index: usize,
    left_index: usize,
    right_index: usize,
    p: &[BigInt],
    q: &[BigInt],
) -> InterlacingReport {
    let strict = check_interlacing_bigint_coeffs(p, q);
    let weak = check_weak_interlacing_bigint_coeffs(p, q);
    InterlacingReport {
        pair_index,
        left_index,
        right_index,
        p: format_poly_bigint_coeffs(p),
        q: format_poly_bigint_coeffs(q),
        strict,
        weak,
        status: interlacing_status(strict, weak),
    }
}

fn interlacing_report_json(report: &InterlacingReport) -> String {
    format!(
        "{{\"pair_index\":{},\"left_index\":{},\"right_index\":{},\
         \"p\":{},\"q\":{},\"strict\":{},\"weak\":{},\"status\":{}}}",
        report.pair_index,
        report.left_index,
        report.right_index,
        json_string(&report.p),
        json_string(&report.q),
        json_bool_option(report.strict),
        json_bool_option(report.weak),
        json_string(&report.status)
    )
}

fn print_interlacing_reports_json(reports: &[InterlacingReport]) {
    println!(
        "{{\"pairs\":[{}]}}",
        reports
            .iter()
            .map(interlacing_report_json)
            .collect::<Vec<_>>()
            .join(",")
    );
}

#[derive(Clone, Debug)]
struct InterlacingProfileReport {
    index: usize,
    polynomial: String,
    previous_count: usize,
    checked_previous_count: usize,
    interlacing_previous_count: usize,
    strict_previous_count: usize,
    weak_previous_count: usize,
    previous_interlacing_indices: Vec<usize>,
    previous: Vec<InterlacingReport>,
}

fn interlacing_report_has_interlacing(report: &InterlacingReport) -> bool {
    matches!(
        report.status.as_str(),
        "strictly_interlace" | "weakly_interlace"
    )
}

fn interlacing_profile_report(
    index: usize,
    polynomial: &[BigInt],
    previous: Vec<InterlacingReport>,
) -> InterlacingProfileReport {
    let previous_interlacing_indices = previous
        .iter()
        .filter(|report| interlacing_report_has_interlacing(report))
        .map(|report| report.left_index)
        .collect::<Vec<_>>();
    InterlacingProfileReport {
        index,
        polynomial: format_poly_bigint_coeffs(polynomial),
        previous_count: index,
        checked_previous_count: previous.len(),
        interlacing_previous_count: previous_interlacing_indices.len(),
        strict_previous_count: previous
            .iter()
            .filter(|report| report.strict == Some(true))
            .count(),
        weak_previous_count: previous
            .iter()
            .filter(|report| report.weak == Some(true))
            .count(),
        previous_interlacing_indices,
        previous,
    }
}

fn interlacing_profile_report_json(report: &InterlacingProfileReport) -> String {
    format!(
        "{{\"index\":{},\"polynomial\":{},\"previous_count\":{},\
         \"checked_previous_count\":{},\"interlacing_previous_count\":{},\
         \"strict_previous_count\":{},\"weak_previous_count\":{},\
         \"previous_interlacing_indices\":{},\"previous\":[{}]}}",
        report.index,
        json_string(&report.polynomial),
        report.previous_count,
        report.checked_previous_count,
        report.interlacing_previous_count,
        report.strict_previous_count,
        report.weak_previous_count,
        json_usize_vec(&report.previous_interlacing_indices),
        report
            .previous
            .iter()
            .map(interlacing_report_json)
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn print_interlacing_profile_reports_json(reports: &[InterlacingProfileReport]) {
    println!(
        "{{\"items\":[{}]}}",
        reports
            .iter()
            .map(interlacing_profile_report_json)
            .collect::<Vec<_>>()
            .join(",")
    );
}

fn cmd_interlacing(args: &[String]) {
    let format = match OutputFormat::from_args(args) {
        Ok(format) => format,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let polys = read_polys_bigint();
    if polys.len() < 2 {
        eprintln!("Need at least two polynomials for interlacing check.");
        return;
    }
    let reports = polys
        .windows(2)
        .enumerate()
        .map(|(i, pair)| {
            let p = strip_trailing_zeros_bigint(&pair[0]);
            let q = strip_trailing_zeros_bigint(&pair[1]);
            interlacing_report(i, i, i + 1, p, q)
        })
        .collect::<Vec<_>>();
    if format == OutputFormat::Json {
        print_interlacing_reports_json(&reports);
        return;
    }
    for report in reports {
        let status = match report.status.as_str() {
            "strictly_interlace" => "strictly interlace",
            "weakly_interlace" => "weakly interlace (shared roots)",
            "do_not_interlace" => "do NOT interlace",
            _ => "incompatible degrees",
        };
        println!("{} & {}: {}", report.p, report.q, status);
    }
}

fn cmd_interlacing_profile(args: &[String]) {
    let format = match OutputFormat::from_args(args) {
        Ok(format) => format,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let polys = read_polys_bigint()
        .into_iter()
        .map(|coeffs| strip_trailing_zeros_bigint(&coeffs).to_vec())
        .collect::<Vec<_>>();
    if polys.is_empty() {
        eprintln!("Need at least one polynomial for interlacing profile.");
        return;
    }

    let mut pair_index = 0;
    let mut reports = Vec::with_capacity(polys.len());
    for current in 0..polys.len() {
        let mut previous = Vec::new();
        for left in (0..current).rev() {
            let report =
                interlacing_report(pair_index, left, current, &polys[left], &polys[current]);
            pair_index += 1;
            let interlaces = interlacing_report_has_interlacing(&report);
            previous.push(report);
            if !interlaces {
                break;
            }
        }
        reports.push(interlacing_profile_report(
            current,
            &polys[current],
            previous,
        ));
    }

    if format == OutputFormat::Json {
        print_interlacing_profile_reports_json(&reports);
        return;
    }

    for report in reports {
        println!(
            "row {}: {} previous rows interlace before first failure (checked {}/{}); \
             strict={}, weak={}; indices={:?}; {}",
            report.index,
            report.interlacing_previous_count,
            report.checked_previous_count,
            report.previous_count,
            report.strict_previous_count,
            report.weak_previous_count,
            report.previous_interlacing_indices,
            report.polynomial
        );
    }
}

fn cmd_properties(args: &[String]) {
    let format = match OutputFormat::from_args(args) {
        Ok(format) => format,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let reports = read_polys_bigint()
        .into_iter()
        .enumerate()
        .map(|(index, coeffs)| {
            let c = strip_trailing_zeros_bigint(&coeffs);
            property_report(index, c)
        })
        .collect::<Vec<_>>();
    if format == OutputFormat::Json {
        print_property_reports_json(&reports);
        return;
    }
    for report in reports {
        println!(
            "{}: {}",
            report.polynomial,
            property_labels(&report).join(", ")
        );
    }
}

fn cmd_gamma_expansion(args: &[String]) {
    let format = match OutputFormat::from_args(args) {
        Ok(format) => format,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let mut had_error = false;
    let mut json_items = Vec::new();
    for (index, coeffs) in read_polys_bigint().into_iter().enumerate() {
        let c = strip_trailing_zeros_bigint(&coeffs);
        match gamma_coefficients_bigint_coeffs(c) {
            Some(gamma) => {
                let polynomial = format_poly_bigint_coeffs(c);
                let expansion = format_gamma_expansion(&gamma, c.len().saturating_sub(1));
                if format == OutputFormat::Json {
                    json_items.push(format!(
                        "{{\"index\":{},\"ok\":true,\"polynomial\":{},\"coefficients\":{},\
                         \"degree\":{},\"gamma\":{},\"expansion\":{}}}",
                        index,
                        json_string(&polynomial),
                        json_bigint_vec(c),
                        polynomial_degree_bigint(c),
                        json_bigint_vec(&gamma),
                        json_string(&expansion)
                    ));
                } else {
                    println!("{polynomial}: gamma {:?}; expansion: {expansion}", gamma);
                }
            }
            None => {
                let polynomial = format_poly_bigint_coeffs(c);
                let error = "not palindromic; no gamma expansion";
                if format == OutputFormat::Json {
                    json_items.push(format!(
                        "{{\"index\":{},\"ok\":false,\"polynomial\":{},\"coefficients\":{},\
                         \"degree\":{},\"error\":{}}}",
                        index,
                        json_string(&polynomial),
                        json_bigint_vec(c),
                        polynomial_degree_bigint(c),
                        json_string(error)
                    ));
                } else {
                    eprintln!("{polynomial}: {error}");
                }
                had_error = true;
            }
        }
    }
    if format == OutputFormat::Json {
        println!(
            "{{\"ok\":{},\"items\":[{}]}}",
            !had_error,
            json_items.join(",")
        );
    }
    if had_error {
        std::process::exit(1);
    }
}

fn format_gamma_expansion(gamma: &[BigInt], degree: usize) -> String {
    let mut out = String::new();
    let zero = BigInt::from(0);
    let one = BigInt::from(1);
    for (i, coeff) in gamma.iter().enumerate() {
        if coeff == &zero {
            continue;
        }
        let negative = coeff < &zero;
        let abs_coeff = if negative {
            -coeff.clone()
        } else {
            coeff.clone()
        };
        let factor = gamma_basis_factor(i, degree.saturating_sub(2 * i));
        let body = if factor == "1" {
            abs_coeff.to_string()
        } else if abs_coeff == one {
            factor
        } else {
            format!("{abs_coeff} {factor}")
        };

        if out.is_empty() {
            if negative {
                out.push('-');
            }
            out.push_str(&body);
        } else if negative {
            out.push_str(" - ");
            out.push_str(&body);
        } else {
            out.push_str(" + ");
            out.push_str(&body);
        }
    }
    if out.is_empty() {
        "0".to_string()
    } else {
        out
    }
}

fn gamma_basis_factor(t_power: usize, one_plus_t_power: usize) -> String {
    let mut factors = Vec::new();
    match t_power {
        0 => {}
        1 => factors.push("t".to_string()),
        n => factors.push(format!("t^{n}")),
    }
    match one_plus_t_power {
        0 => {}
        1 => factors.push("(1+t)".to_string()),
        n => factors.push(format!("(1+t)^{n}")),
    }
    if factors.is_empty() {
        "1".to_string()
    } else {
        factors.join(" ")
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SequenceKind {
    Eulerian,
    Narayana,
    TypeBEulerian,
    ChebyshevT,
    ChebyshevU,
    Hermite,
}

impl SequenceKind {
    fn parse(input: &str) -> Option<Self> {
        match input {
            "eulerian" => Some(Self::Eulerian),
            "narayana" => Some(Self::Narayana),
            "type-b-eulerian" | "type-b" | "type_b_eulerian" | "type_b" => {
                Some(Self::TypeBEulerian)
            }
            "chebyshev-t" | "chebyshev_t" | "t" => Some(Self::ChebyshevT),
            "chebyshev-u" | "chebyshev_u" | "u" => Some(Self::ChebyshevU),
            "hermite" => Some(Self::Hermite),
            _ => None,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Eulerian => "eulerian",
            Self::Narayana => "narayana",
            Self::TypeBEulerian => "type_b_eulerian",
            Self::ChebyshevT => "chebyshev_t",
            Self::ChebyshevU => "chebyshev_u",
            Self::Hermite => "hermite",
        }
    }

    fn polynomials(self, max_n: usize) -> Vec<Vec<BigInt>> {
        match self {
            Self::Eulerian => polynomial_tools::sequences::eulerian_polynomials_bigint(max_n),
            Self::Narayana => polynomial_tools::sequences::narayana_polynomials_bigint(max_n),
            Self::TypeBEulerian => {
                polynomial_tools::sequences::type_b_eulerian_polynomials_bigint(max_n)
            }
            Self::ChebyshevT => polynomial_tools::sequences::chebyshev_polynomials_t_bigint(max_n),
            Self::ChebyshevU => polynomial_tools::sequences::chebyshev_polynomials_u_bigint(max_n),
            Self::Hermite => polynomial_tools::sequences::hermite_polynomials_bigint(max_n),
        }
    }

    fn label(self, index: usize) -> String {
        match self {
            Self::Eulerian => format!("A_{}(t)", index + 1),
            Self::Narayana => format!("N_{}(t)", index + 1),
            Self::TypeBEulerian => format!("B_{index}(t)"),
            Self::ChebyshevT => format!("T_{index}(t)"),
            Self::ChebyshevU => format!("U_{index}(t)"),
            Self::Hermite => format!("He_{index}(t)"),
        }
    }
}

fn cmd_sequence(args: &[String]) {
    let mut positional = Vec::new();
    let mut format = OutputFormat::Text;
    for arg in args {
        match arg.as_str() {
            "--json" => format = OutputFormat::Json,
            other => positional.push(other),
        }
    }
    if positional.len() != 2 {
        print_sequence_help();
        std::process::exit(1);
    }
    let Some(kind) = SequenceKind::parse(positional[0]) else {
        eprintln!("unknown sequence: {}", positional[0]);
        print_sequence_help();
        std::process::exit(1);
    };
    let max_n = positional[1].parse::<usize>().unwrap_or_else(|_| {
        eprintln!(
            "expected a nonnegative integer max-n, got '{}'",
            positional[1]
        );
        std::process::exit(1);
    });
    let polynomials = kind.polynomials(max_n);
    if format == OutputFormat::Json {
        let items = polynomials
            .iter()
            .enumerate()
            .map(|(index, coeffs)| {
                let c = strip_trailing_zeros_bigint(coeffs);
                format!(
                    "{{\"index\":{},\"label\":{},\"polynomial\":{},\"coefficients\":{},\
                     \"degree\":{}}}",
                    index,
                    json_string(&kind.label(index)),
                    json_string(&format_poly_bigint_coeffs(c)),
                    json_bigint_vec(c),
                    polynomial_degree_bigint(c)
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "{{\"sequence\":{},\"max_n\":{},\"polynomials\":[{}]}}",
            json_string(kind.name()),
            max_n,
            items
        );
        return;
    }
    for (index, coeffs) in polynomials.iter().enumerate() {
        let c = strip_trailing_zeros_bigint(coeffs);
        println!("{} = {}", kind.label(index), format_poly_bigint_coeffs(c));
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RecurrenceOutputFormat {
    Text,
    Json,
    Python,
}

fn cmd_recurrence(args: &[String]) {
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_recurrence_help();
        return;
    }

    let mut search = AdaptiveSearchOptions::default();
    let mut format = RecurrenceOutputFormat::Text;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--json" => {
                format = RecurrenceOutputFormat::Json;
            }
            "--python" => {
                format = RecurrenceOutputFormat::Python;
            }
            "--format" => {
                let value = match next_option_value(args, &mut i, "--format") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
                format = match value {
                    "text" => RecurrenceOutputFormat::Text,
                    "json" => RecurrenceOutputFormat::Json,
                    "python" => RecurrenceOutputFormat::Python,
                    other => {
                        eprintln!("unknown recurrence output format: {other}");
                        return;
                    }
                };
            }
            "--skip-prefix" => {
                search.skip_prefix = match parse_usize_option(args, &mut i, "--skip-prefix") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--min-rec-len" => {
                search.min_rec_len = match parse_usize_option(args, &mut i, "--min-rec-len") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--max-rec-len" => {
                search.max_rec_len = match parse_usize_option(args, &mut i, "--max-rec-len") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--min-var-deg" => {
                search.min_var_deg = match parse_usize_option(args, &mut i, "--min-var-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--max-var-deg" => {
                search.max_var_deg = match parse_usize_option(args, &mut i, "--max-var-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--min-idx-deg" => {
                search.min_idx_deg = match parse_usize_option(args, &mut i, "--min-idx-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--max-idx-deg" => {
                search.max_idx_deg = match parse_usize_option(args, &mut i, "--max-idx-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--min-diff-deg" => {
                search.min_diff_deg = match parse_usize_option(args, &mut i, "--min-diff-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--max-diff-deg" => {
                search.max_diff_deg = match parse_usize_option(args, &mut i, "--max-diff-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--inhomogeneous" => {
                search.try_inhomogeneous = true;
            }
            "--min-inhomo-var-deg" => {
                search.try_inhomogeneous = true;
                search.min_inhomo_var_deg =
                    match parse_usize_option(args, &mut i, "--min-inhomo-var-deg") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            return;
                        }
                    };
            }
            "--max-inhomo-var-deg" => {
                search.try_inhomogeneous = true;
                search.max_inhomo_var_deg =
                    match parse_usize_option(args, &mut i, "--max-inhomo-var-deg") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            return;
                        }
                    };
            }
            "--min-inhomo-idx-deg" => {
                search.try_inhomogeneous = true;
                search.min_inhomo_idx_deg =
                    match parse_usize_option(args, &mut i, "--min-inhomo-idx-deg") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            return;
                        }
                    };
            }
            "--max-inhomo-idx-deg" => {
                search.try_inhomogeneous = true;
                search.max_inhomo_idx_deg =
                    match parse_usize_option(args, &mut i, "--max-inhomo-idx-deg") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            return;
                        }
                    };
            }
            "--denominator" | "--try-denominator" => {
                search.try_denominator = true;
            }
            "--alternating-sign" => {
                search.try_alternating_sign = true;
            }
            "--max-denom-var-deg" => {
                search.try_denominator = true;
                search.max_denom_var_deg =
                    match parse_usize_option(args, &mut i, "--max-denom-var-deg") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            return;
                        }
                    };
            }
            "--max-denom-idx-deg" => {
                search.try_denominator = true;
                search.max_denom_idx_deg =
                    match parse_usize_option(args, &mut i, "--max-denom-idx-deg") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            return;
                        }
                    };
            }
            "--min-margin" => {
                search.min_margin = match parse_usize_option(args, &mut i, "--min-margin") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--fit-extra-rows" => {
                search.fit_extra_rows = match parse_usize_option(args, &mut i, "--fit-extra-rows") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--no-verify" => {
                search.no_verify = true;
            }
            "--modular-prefilter" => {
                search.modular_prefilter = true;
            }
            "--verbose" => {
                search.verbose = true;
            }
            other => {
                eprintln!("unknown recurrence option: {other}");
                return;
            }
        }
        i += 1;
    }

    let polys = read_polys_rational();
    if polys.len() < 3 {
        eprintln!("Need at least 3 polynomials for recurrence search.");
        return;
    }

    eprintln!(
        "Searching for recurrence among {} polynomials...",
        polys.len()
    );
    match find_recurrence_adaptive_rational(&polys, &search) {
        Some(res) => {
            let searched_polys = polys.get(search.skip_prefix..).unwrap_or(&[]);
            let initial_count = res.recurrence.max_offset().min(searched_polys.len());
            let initial_polys = &searched_polys[..initial_count];
            match format {
                RecurrenceOutputFormat::Json => {
                    let recurrence_json = RecurrenceJson::from_recurrence_rational(
                        &res.recurrence,
                        1,
                        initial_polys,
                        Some(RecurrenceJsonSearch {
                            recurrence_text: res.recurrence.to_string(),
                            source_rows: polys.len(),
                            skip_prefix: search.skip_prefix,
                            unknowns: res.num_unknowns,
                            weighted_unknowns: res.weighted_unknowns,
                            equations: res.num_equations,
                            fit_polynomials: res.fit_polynomials,
                            verification_polynomials: res.verification_polynomials,
                            candidates_tried: res.candidates_tried,
                            options: RecurrenceOptionsJson::from(&res.opts),
                        }),
                    );
                    println!(
                        "{}",
                        serde_json::to_string_pretty(&recurrence_json).unwrap()
                    );
                }
                RecurrenceOutputFormat::Python => {
                    println!(
                        "{}",
                        res.recurrence.to_python_definition_rational(initial_polys)
                    );
                }
                RecurrenceOutputFormat::Text => {
                    println!("{}", res.recurrence);
                }
            }
            eprintln!(
                "Found with {} unknowns, weighted score {}, {} equations, {} fit rows, \
                 {} verification rows ({} candidates tried)",
                res.num_unknowns,
                res.weighted_unknowns,
                res.num_equations,
                res.fit_polynomials,
                res.verification_polynomials,
                res.candidates_tried
            );
        }
        None => {
            eprintln!("No recurrence found within the search bounds.");
        }
    }
}

fn format_rational_row(coeffs: &[RecurrenceBigRational]) -> String {
    coeffs
        .iter()
        .map(format_rational_coeff)
        .collect::<Vec<_>>()
        .join(", ")
}

fn read_recurrence_json(path: &str) -> Result<RecurrenceJson, String> {
    let mut input = String::new();
    if path == "-" {
        io::stdin()
            .read_to_string(&mut input)
            .map_err(|e| format!("failed to read recurrence JSON from stdin: {e}"))?;
    } else {
        input = fs::read_to_string(path)
            .map_err(|e| format!("failed to read recurrence JSON `{path}`: {e}"))?;
    }
    serde_json::from_str(&input).map_err(|e| format!("failed to parse recurrence JSON: {e}"))
}

fn cmd_recurrence_generate(args: &[String]) {
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_recurrence_generate_help();
        return;
    }

    let mut recurrence_path: Option<String> = None;
    let mut rows: Option<usize> = None;
    let mut additional: Option<usize> = None;
    let mut format = OutputFormat::Text;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--recurrence" => {
                i += 1;
                recurrence_path = args.get(i).cloned();
            }
            "--rows" => {
                i += 1;
                rows = args.get(i).and_then(|s| s.parse().ok());
            }
            "--additional" => {
                i += 1;
                additional = args.get(i).and_then(|s| s.parse().ok());
            }
            "--json" => {
                format = OutputFormat::Json;
            }
            other => {
                eprintln!("unknown recurrence-generate option: {other}");
                return;
            }
        }
        i += 1;
    }

    let Some(path) = recurrence_path else {
        eprintln!("recurrence-generate needs --recurrence <file|->");
        return;
    };
    if rows.is_some() && additional.is_some() {
        eprintln!("use only one of --rows or --additional");
        return;
    }

    let recurrence_json = match read_recurrence_json(&path) {
        Ok(value) => value,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let (recurrence, first_index, initial_polys) = match recurrence_json.to_recurrence_parts() {
        Ok(parts) => parts,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };

    let total_rows = if let Some(rows) = rows {
        rows
    } else if let Some(additional) = additional {
        initial_polys.len() + additional
    } else {
        eprintln!("recurrence-generate needs --rows <n> or --additional <n>");
        return;
    };

    let generated = match recurrence.generate_rows_rational(&initial_polys, first_index, total_rows)
    {
        Ok(rows) => rows,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };

    if format == OutputFormat::Json {
        let rows_json: Vec<Vec<String>> = generated
            .iter()
            .map(|row| row.iter().map(format_rational_coeff).collect())
            .collect();
        println!(
            "{}",
            serde_json::to_string_pretty(&serde_json::json!({
                "first_index": first_index,
                "polynomials": rows_json,
            }))
            .unwrap()
        );
    } else {
        for row in &generated {
            println!("{}", format_rational_row(row));
        }
    }
}

fn print_bench_help() {
    println!("Usage:");
    println!("  polytool bench recurrence-fixtures [options]");
    println!("  polytool bench interlacing [options]");
    println!("  polytool bench --help");
    println!();
    println!("Run built-in timing suites and print tab-separated results.");
    println!();
    println!("Subcommands:");
    println!("  recurrence-fixtures   Time adaptive recurrence search on fixture rows");
    println!("  compare               Compare two recurrence-fixture JSON benchmark runs");
    println!("  interlacing           Time consecutive interlacing checks on a sequence");
    println!();
    println!("Recurrence fixture options:");
    println!("  --base <dir>              Fixture directory containing manifest.tsv");
    println!("  --only <substring>        Run only fixture slugs containing substring");
    println!("  --repeat <n>              Repeat each fixture n times (default: 1)");
    println!("  --no-modular-prefilter    Disable modular recurrence prefilter");
    println!("  --summary                 Append fixture and category summaries to stdout");
    println!("  --report <path.md>        Write a Markdown benchmark report");
    println!("  --format <tsv|json>       Output format (default: tsv)");
    println!("  --json                    Alias for --format json");
    println!();
    println!("Compare options:");
    println!("  polytool bench compare [options] <old.json> <new.json>");
    println!("  --top <n>                 Number of worst regressions to show (default: 10)");
    println!("  --format <tsv|json>       Output format (default: tsv)");
    println!();
    println!("Interlacing options:");
    println!("  --sequence <name>         Sequence name (default: eulerian)");
    println!("  --max-n <n>               Maximum n/index to generate (default: 20)");
    println!("  --repeat <n>              Repeat each pair n times (default: 3)");
}

#[derive(Debug)]
struct BenchRecurrenceOptions {
    base: PathBuf,
    repeat: usize,
    modular_prefilter: bool,
    only: Option<String>,
    summary: bool,
    report: Option<PathBuf>,
    format: BenchOutputFormat,
}

impl Default for BenchRecurrenceOptions {
    fn default() -> Self {
        Self {
            base: Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("fixtures")
                .join("recurrence-benchmarks"),
            repeat: 1,
            modular_prefilter: true,
            only: None,
            summary: false,
            report: None,
            format: BenchOutputFormat::Tsv,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum BenchOutputFormat {
    Tsv,
    Json,
}

impl BenchOutputFormat {
    fn parse(value: &str) -> Result<Self, String> {
        match value {
            "tsv" | "text" => Ok(Self::Tsv),
            "json" => Ok(Self::Json),
            other => Err(format!("unknown bench output format `{other}`")),
        }
    }
}

#[derive(Debug)]
struct RecurrenceFixtureManifestRow {
    slug: String,
    title: String,
    description: String,
    suggested_args: String,
    rows_file: String,
}

#[derive(Debug)]
struct BenchRecurrenceRun {
    slug: String,
    title: String,
    description: String,
    categories: Vec<String>,
    run: usize,
    found: bool,
    elapsed_ms: f64,
    candidates: usize,
    unknowns: usize,
    weighted_unknowns: usize,
    fit_rows: usize,
    verify_rows: usize,
    recurrence: String,
    diagnostics: Option<AdaptiveSearchDiagnostics>,
}

#[derive(Debug)]
struct BenchRecurrenceFixtureSummary {
    slug: String,
    title: String,
    description: String,
    categories: Vec<String>,
    runs: usize,
    found_runs: usize,
    min_ms: f64,
    median_ms: f64,
    mean_ms: f64,
    max_ms: f64,
    candidates: usize,
    unknowns: usize,
    weighted_unknowns: usize,
    fit_rows: usize,
    verify_rows: usize,
    recurrence: String,
}

#[derive(Debug)]
struct BenchRecurrenceCategorySummary {
    category: String,
    fixtures: usize,
    runs: usize,
    found_runs: usize,
    mean_ms: f64,
    max_ms: f64,
}

fn cmd_bench(args: &[String]) {
    let Some((subcommand, rest)) = args.split_first() else {
        print_bench_help();
        std::process::exit(1);
    };
    match subcommand.as_str() {
        "recurrence-fixtures" | "recurrence" => cmd_bench_recurrence_fixtures(rest),
        "compare" => cmd_bench_compare(rest),
        "interlacing" => cmd_bench_interlacing(rest),
        _ => {
            eprintln!("unknown bench subcommand: {subcommand}");
            print_bench_help();
            std::process::exit(1);
        }
    }
}

fn parse_bench_recurrence_options(args: &[String]) -> Result<BenchRecurrenceOptions, String> {
    let mut options = BenchRecurrenceOptions::default();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--base" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--base expects a directory".to_string())?;
                options.base = PathBuf::from(value);
            }
            "--only" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--only expects a substring".to_string())?;
                options.only = Some(value.clone());
            }
            "--repeat" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--repeat expects a positive integer".to_string())?;
                options.repeat = value
                    .parse::<usize>()
                    .map_err(|_| format!("--repeat expects a positive integer, got '{value}'"))?
                    .max(1);
            }
            "--no-modular-prefilter" => options.modular_prefilter = false,
            "--summary" => options.summary = true,
            "--report" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--report expects a path".to_string())?;
                options.report = Some(PathBuf::from(value));
            }
            "--format" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--format expects tsv or json".to_string())?;
                options.format = BenchOutputFormat::parse(value)?;
            }
            "--json" => options.format = BenchOutputFormat::Json,
            other => return Err(format!("unknown recurrence benchmark option: {other}")),
        }
        i += 1;
    }
    Ok(options)
}

fn parse_recurrence_fixture_manifest(
    base: &Path,
) -> Result<Vec<RecurrenceFixtureManifestRow>, String> {
    let manifest_path = base.join("manifest.tsv");
    let manifest = fs::read_to_string(&manifest_path)
        .map_err(|e| format!("failed to read `{}`: {e}", manifest_path.display()))?;
    manifest
        .lines()
        .skip(1)
        .map(|line| {
            let cols = line.split('\t').collect::<Vec<_>>();
            if cols.len() != 7 {
                return Err(format!("manifest row should have 7 columns: {line}"));
            }
            Ok(RecurrenceFixtureManifestRow {
                slug: cols[0].to_string(),
                title: cols[1].to_string(),
                description: cols[2].to_string(),
                suggested_args: cols[3].to_string(),
                rows_file: cols[5].to_string(),
            })
        })
        .collect()
}

fn read_recurrence_fixture_rows(path: &Path) -> Result<Vec<Vec<RecurrenceBigRational>>, String> {
    let input = fs::read_to_string(path)
        .map_err(|e| format!("failed to read `{}`: {e}", path.display()))?;
    input
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|coeff| {
                    parse_rational_coeff(coeff.trim())
                        .map_err(|e| format!("failed to parse `{}`: {e}", line.trim()))
                })
                .collect()
        })
        .collect()
}

fn parse_fixture_search_options(
    suggested_args: &str,
    modular_prefilter: bool,
) -> Result<AdaptiveSearchOptions, String> {
    let mut search = AdaptiveSearchOptions {
        modular_prefilter,
        ..Default::default()
    };
    let mut args = suggested_args.split_whitespace();
    while let Some(arg) = args.next() {
        match arg {
            "--skip-prefix" => search.skip_prefix = parse_fixture_usize(arg, args.next())?,
            "--min-rec-len" => search.min_rec_len = parse_fixture_usize(arg, args.next())?,
            "--max-rec-len" => search.max_rec_len = parse_fixture_usize(arg, args.next())?,
            "--min-var-deg" => search.min_var_deg = parse_fixture_usize(arg, args.next())?,
            "--max-var-deg" => search.max_var_deg = parse_fixture_usize(arg, args.next())?,
            "--min-idx-deg" => search.min_idx_deg = parse_fixture_usize(arg, args.next())?,
            "--max-idx-deg" => search.max_idx_deg = parse_fixture_usize(arg, args.next())?,
            "--min-diff-deg" => search.min_diff_deg = parse_fixture_usize(arg, args.next())?,
            "--max-diff-deg" => search.max_diff_deg = parse_fixture_usize(arg, args.next())?,
            "--inhomogeneous" => search.try_inhomogeneous = true,
            "--min-inhomo-var-deg" => {
                search.min_inhomo_var_deg = parse_fixture_usize(arg, args.next())?;
            }
            "--max-inhomo-var-deg" => {
                search.max_inhomo_var_deg = parse_fixture_usize(arg, args.next())?;
            }
            "--min-inhomo-idx-deg" => {
                search.min_inhomo_idx_deg = parse_fixture_usize(arg, args.next())?;
            }
            "--max-inhomo-idx-deg" => {
                search.max_inhomo_idx_deg = parse_fixture_usize(arg, args.next())?;
            }
            "--denominator" | "--try-denominator" => search.try_denominator = true,
            "--max-denom-var-deg" => {
                search.try_denominator = true;
                search.max_denom_var_deg = parse_fixture_usize(arg, args.next())?;
            }
            "--max-denom-idx-deg" => {
                search.try_denominator = true;
                search.max_denom_idx_deg = parse_fixture_usize(arg, args.next())?;
            }
            "--alternating-sign" => search.try_alternating_sign = true,
            "--min-margin" => search.min_margin = parse_fixture_usize(arg, args.next())?,
            "--fit-extra-rows" => search.fit_extra_rows = parse_fixture_usize(arg, args.next())?,
            "--no-verify" => search.no_verify = true,
            other => return Err(format!("unsupported fixture search option `{other}`")),
        }
    }
    Ok(search)
}

fn recurrence_fixture_categories(
    fixture: &RecurrenceFixtureManifestRow,
    search: &AdaptiveSearchOptions,
) -> Vec<String> {
    let mut categories = vec![if fixture.slug.contains("_oeis_") {
        "oeis".to_string()
    } else {
        "synthetic".to_string()
    }];
    if search.skip_prefix > 0 {
        categories.push("skip-prefix".to_string());
    }
    if search.max_rec_len >= 4 {
        categories.push("high-order".to_string());
    }
    if search.max_diff_deg > 0 {
        categories.push("derivative".to_string());
    }
    if search.max_diff_deg >= 2 {
        categories.push("higher-derivative".to_string());
    }
    if search.max_idx_deg > 0
        || (search.try_denominator && search.max_denom_idx_deg > 0)
        || (search.try_inhomogeneous && search.max_inhomo_idx_deg > 0)
    {
        categories.push("index-dependent".to_string());
    }
    if search.try_denominator {
        categories.push("denominator".to_string());
    }
    if search.try_inhomogeneous {
        categories.push("inhomogeneous".to_string());
    }
    if search.try_alternating_sign {
        categories.push("alternating-sign".to_string());
    }
    if fixture.slug.contains("closed_form") {
        categories.push("closed-form".to_string());
    }
    categories
}

fn bench_mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn bench_median(values: &mut [f64]) -> f64 {
    values.sort_by(|a, b| a.total_cmp(b));
    let middle = values.len() / 2;
    if values.len().is_multiple_of(2) {
        (values[middle - 1] + values[middle]) / 2.0
    } else {
        values[middle]
    }
}

fn markdown_cell(text: &str) -> String {
    text.replace('\n', " ").replace('|', "\\|")
}

fn recurrence_bench_fixture_summaries(
    records: &[BenchRecurrenceRun],
) -> Vec<BenchRecurrenceFixtureSummary> {
    let mut by_fixture: BTreeMap<&str, Vec<&BenchRecurrenceRun>> = BTreeMap::new();
    for record in records {
        by_fixture.entry(&record.slug).or_default().push(record);
    }

    by_fixture
        .into_iter()
        .map(|(slug, runs)| {
            let found_runs = runs.iter().filter(|record| record.found).count();
            let mut elapsed = runs
                .iter()
                .map(|record| record.elapsed_ms)
                .collect::<Vec<_>>();
            let min_ms = elapsed.iter().copied().fold(f64::INFINITY, f64::min);
            let max_ms = elapsed.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let mean_ms = bench_mean(&elapsed);
            let median_ms = bench_median(&mut elapsed);
            let sample = runs[0];
            BenchRecurrenceFixtureSummary {
                slug: slug.to_string(),
                title: sample.title.clone(),
                description: sample.description.clone(),
                categories: sample.categories.clone(),
                runs: runs.len(),
                found_runs,
                min_ms,
                median_ms,
                mean_ms,
                max_ms,
                candidates: sample.candidates,
                unknowns: sample.unknowns,
                weighted_unknowns: sample.weighted_unknowns,
                fit_rows: sample.fit_rows,
                verify_rows: sample.verify_rows,
                recurrence: sample.recurrence.clone(),
            }
        })
        .collect()
}

fn recurrence_bench_category_summary_values(
    records: &[BenchRecurrenceRun],
) -> Vec<BenchRecurrenceCategorySummary> {
    recurrence_bench_category_summaries(records)
        .into_iter()
        .map(|(category, fixtures, runs)| {
            let found_runs = runs.iter().filter(|record| record.found).count();
            let elapsed = runs
                .iter()
                .map(|record| record.elapsed_ms)
                .collect::<Vec<_>>();
            let max_ms = elapsed.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            BenchRecurrenceCategorySummary {
                category,
                fixtures: fixtures.len(),
                runs: runs.len(),
                found_runs,
                mean_ms: bench_mean(&elapsed),
                max_ms,
            }
        })
        .collect()
}

fn print_recurrence_bench_summary(records: &[BenchRecurrenceRun]) {
    if records.is_empty() {
        println!("# fixture_summary");
        println!("slug\ttitle\tcategories\truns\tfound_runs\tmin_ms\tmedian_ms\tmean_ms\tmax_ms\tcandidates\tunknowns\tweighted\tfit_rows\tverify_rows");
        println!("# category_summary");
        println!("category\tfixtures\truns\tfound_runs\tmean_ms\tmax_ms");
        return;
    }

    println!("# fixture_summary");
    println!("slug\ttitle\tcategories\truns\tfound_runs\tmin_ms\tmedian_ms\tmean_ms\tmax_ms\tcandidates\tunknowns\tweighted\tfit_rows\tverify_rows");
    for summary in recurrence_bench_fixture_summaries(records) {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{}\t{}\t{}\t{}\t{}",
            summary.slug,
            summary.title,
            summary.categories.join(","),
            summary.runs,
            summary.found_runs,
            summary.min_ms,
            summary.median_ms,
            summary.mean_ms,
            summary.max_ms,
            summary.candidates,
            summary.unknowns,
            summary.weighted_unknowns,
            summary.fit_rows,
            summary.verify_rows,
        );
    }

    println!("# category_summary");
    println!("category\tfixtures\truns\tfound_runs\tmean_ms\tmax_ms");
    for summary in recurrence_bench_category_summary_values(records) {
        println!(
            "{}\t{}\t{}\t{}\t{:.3}\t{:.3}",
            summary.category,
            summary.fixtures,
            summary.runs,
            summary.found_runs,
            summary.mean_ms,
            summary.max_ms,
        );
    }
}

fn recurrence_bench_category_summaries(
    records: &[BenchRecurrenceRun],
) -> Vec<(String, BTreeSet<&str>, Vec<&BenchRecurrenceRun>)> {
    let mut categories: BTreeMap<String, (BTreeSet<&str>, Vec<&BenchRecurrenceRun>)> =
        BTreeMap::new();
    for record in records {
        for category in &record.categories {
            let entry = categories.entry(category.clone()).or_default();
            entry.0.insert(&record.slug);
            entry.1.push(record);
        }
    }
    categories
        .into_iter()
        .map(|(category, (fixtures, runs))| (category, fixtures, runs))
        .collect()
}

fn write_recurrence_bench_report(
    options: &BenchRecurrenceOptions,
    records: &[BenchRecurrenceRun],
) -> Result<(), String> {
    let Some(path) = &options.report else {
        return Ok(());
    };

    let mut report = String::new();
    report.push_str("# Recurrence Fixture Benchmark Report\n\n");
    report.push_str(&format!("- Base: `{}`\n", options.base.display()));
    report.push_str(&format!("- Repeat: `{}`\n", options.repeat));
    report.push_str(&format!(
        "- Modular prefilter: `{}`\n",
        if options.modular_prefilter {
            "enabled"
        } else {
            "disabled"
        }
    ));
    if let Some(only) = &options.only {
        report.push_str(&format!("- Filter: `{}`\n", only));
    }
    report.push_str("\n## Category Summary\n\n");
    report.push_str("| category | fixtures | runs | found | mean ms | max ms |\n");
    report.push_str("|---|---:|---:|---:|---:|---:|\n");
    for summary in recurrence_bench_category_summary_values(records) {
        report.push_str(&format!(
            "| {} | {} | {} | {} | {:.3} | {:.3} |\n",
            markdown_cell(&summary.category),
            summary.fixtures,
            summary.runs,
            summary.found_runs,
            summary.mean_ms,
            summary.max_ms,
        ));
    }

    report.push_str("\n## Fixture Summary\n\n");
    report.push_str("| slug | title | description | categories | runs | found | min ms | median ms | mean ms | max ms | unknowns | weighted | recurrence |\n");
    report.push_str("|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|\n");
    for summary in recurrence_bench_fixture_summaries(records) {
        report.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {:.3} | {:.3} | {:.3} | {:.3} | {} | {} | {} |\n",
            markdown_cell(&summary.slug),
            markdown_cell(&summary.title),
            markdown_cell(&summary.description),
            markdown_cell(&summary.categories.join(", ")),
            summary.runs,
            summary.found_runs,
            summary.min_ms,
            summary.median_ms,
            summary.mean_ms,
            summary.max_ms,
            summary.unknowns,
            summary.weighted_unknowns,
            markdown_cell(&summary.recurrence),
        ));
    }

    if records.is_empty() {
        report.push_str("\nNo fixtures matched this benchmark filter.\n");
    }

    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)
            .map_err(|e| format!("failed to create `{}`: {e}", parent.display()))?;
    }
    fs::write(path, report).map_err(|e| format!("failed to write `{}`: {e}", path.display()))
}

fn recurrence_diagnostics_json(diagnostics: &AdaptiveSearchDiagnostics) -> Value {
    json!({
        "generated_candidates": diagnostics.generated_candidates,
        "considered_candidates": diagnostics.considered_candidates,
        "insufficient_fit_rows": diagnostics.insufficient_fit_rows,
        "equation_bound_rejections": diagnostics.equation_bound_rejections,
        "degree_bound_rejections": diagnostics.degree_bound_rejections,
        "modular_prefilter_rejections": diagnostics.modular_prefilter_rejections,
        "exact_solve_attempts": diagnostics.exact_solve_attempts,
        "failed_exact_solves": diagnostics.failed_exact_solves,
        "heldout_verification_failures": diagnostics.heldout_verification_failures,
        "denominator_escalation_entered": diagnostics.denominator_escalation_entered,
    })
}

fn recurrence_bench_run_json(record: &BenchRecurrenceRun) -> Value {
    json!({
        "slug": &record.slug,
        "title": &record.title,
        "description": &record.description,
        "categories": &record.categories,
        "run": record.run,
        "found": record.found,
        "elapsed_ms": record.elapsed_ms,
        "candidates": record.candidates,
        "unknowns": record.unknowns,
        "weighted_unknowns": record.weighted_unknowns,
        "fit_rows": record.fit_rows,
        "verify_rows": record.verify_rows,
        "recurrence": &record.recurrence,
        "diagnostics": record
            .diagnostics
            .as_ref()
            .map(recurrence_diagnostics_json)
            .unwrap_or(Value::Null),
    })
}

fn recurrence_bench_fixture_summary_json(summary: &BenchRecurrenceFixtureSummary) -> Value {
    json!({
        "slug": &summary.slug,
        "title": &summary.title,
        "description": &summary.description,
        "categories": &summary.categories,
        "runs": summary.runs,
        "found_runs": summary.found_runs,
        "min_ms": summary.min_ms,
        "median_ms": summary.median_ms,
        "mean_ms": summary.mean_ms,
        "max_ms": summary.max_ms,
        "candidates": summary.candidates,
        "unknowns": summary.unknowns,
        "weighted_unknowns": summary.weighted_unknowns,
        "fit_rows": summary.fit_rows,
        "verify_rows": summary.verify_rows,
        "recurrence": &summary.recurrence,
    })
}

fn recurrence_bench_category_summary_json(summary: &BenchRecurrenceCategorySummary) -> Value {
    json!({
        "category": &summary.category,
        "fixtures": summary.fixtures,
        "runs": summary.runs,
        "found_runs": summary.found_runs,
        "mean_ms": summary.mean_ms,
        "max_ms": summary.max_ms,
    })
}

fn recurrence_bench_json(
    options: &BenchRecurrenceOptions,
    records: &[BenchRecurrenceRun],
) -> Value {
    let fixture_summaries = recurrence_bench_fixture_summaries(records);
    let category_summaries = recurrence_bench_category_summary_values(records);
    json!({
        "schema": "polynomial-tools.bench.recurrence-fixtures.v1",
        "command": "recurrence-fixtures",
        "options": {
            "base": options.base.display().to_string(),
            "repeat": options.repeat,
            "modular_prefilter": options.modular_prefilter,
            "only": &options.only,
        },
        "runs": records.iter().map(recurrence_bench_run_json).collect::<Vec<_>>(),
        "fixture_summaries": fixture_summaries
            .iter()
            .map(recurrence_bench_fixture_summary_json)
            .collect::<Vec<_>>(),
        "category_summaries": category_summaries
            .iter()
            .map(recurrence_bench_category_summary_json)
            .collect::<Vec<_>>(),
    })
}

#[derive(Debug)]
struct BenchCompareOptions {
    old_path: PathBuf,
    new_path: PathBuf,
    top: usize,
    format: BenchOutputFormat,
}

#[derive(Debug, Clone)]
struct BenchComparePoint {
    mean_ms: f64,
    found_runs: usize,
    runs: usize,
}

#[derive(Debug)]
struct BenchCompareRow {
    kind: String,
    name: String,
    old_mean_ms: Option<f64>,
    new_mean_ms: Option<f64>,
    speedup: Option<f64>,
    delta_ms: Option<f64>,
    old_found_runs: Option<usize>,
    old_runs: Option<usize>,
    new_found_runs: Option<usize>,
    new_runs: Option<usize>,
    status: String,
}

fn parse_bench_compare_options(args: &[String]) -> Result<BenchCompareOptions, String> {
    let mut top = 10;
    let mut format = BenchOutputFormat::Tsv;
    let mut paths = Vec::new();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--top" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--top expects a positive integer".to_string())?;
                top = value
                    .parse::<usize>()
                    .map_err(|_| format!("--top expects a positive integer, got '{value}'"))?
                    .max(1);
            }
            "--format" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--format expects tsv or json".to_string())?;
                format = BenchOutputFormat::parse(value)?;
            }
            "--json" => format = BenchOutputFormat::Json,
            other if other.starts_with('-') => {
                return Err(format!("unknown benchmark compare option: {other}"));
            }
            path => paths.push(PathBuf::from(path)),
        }
        i += 1;
    }

    if paths.len() != 2 {
        return Err("bench compare expects exactly two JSON paths".to_string());
    }
    Ok(BenchCompareOptions {
        old_path: paths.remove(0),
        new_path: paths.remove(0),
        top,
        format,
    })
}

fn read_benchmark_json(path: &Path) -> Result<Value, String> {
    let text = fs::read_to_string(path)
        .map_err(|e| format!("failed to read benchmark JSON `{}`: {e}", path.display()))?;
    let value: Value = serde_json::from_str(&text)
        .map_err(|e| format!("failed to parse benchmark JSON `{}`: {e}", path.display()))?;
    let schema = value
        .get("schema")
        .and_then(Value::as_str)
        .ok_or_else(|| format!("benchmark JSON `{}` has no schema", path.display()))?;
    if schema != "polynomial-tools.bench.recurrence-fixtures.v1" {
        return Err(format!(
            "unsupported benchmark JSON schema `{schema}` in `{}`",
            path.display()
        ));
    }
    Ok(value)
}

fn value_usize(value: &Value, key: &str) -> Result<usize, String> {
    value
        .get(key)
        .and_then(Value::as_u64)
        .and_then(|n| usize::try_from(n).ok())
        .ok_or_else(|| format!("benchmark summary entry is missing integer `{key}`"))
}

fn value_f64(value: &Value, key: &str) -> Result<f64, String> {
    value
        .get(key)
        .and_then(Value::as_f64)
        .ok_or_else(|| format!("benchmark summary entry is missing number `{key}`"))
}

fn value_string(value: &Value, key: &str) -> Result<String, String> {
    value
        .get(key)
        .and_then(Value::as_str)
        .map(ToString::to_string)
        .ok_or_else(|| format!("benchmark summary entry is missing string `{key}`"))
}

fn extract_benchmark_points(
    report: &Value,
    array_key: &str,
    name_key: &str,
) -> Result<BTreeMap<String, BenchComparePoint>, String> {
    let entries = report
        .get(array_key)
        .and_then(Value::as_array)
        .ok_or_else(|| format!("benchmark JSON missing array `{array_key}`"))?;
    let mut map = BTreeMap::new();
    for entry in entries {
        let name = value_string(entry, name_key)?;
        map.insert(
            name,
            BenchComparePoint {
                mean_ms: value_f64(entry, "mean_ms")?,
                found_runs: value_usize(entry, "found_runs")?,
                runs: value_usize(entry, "runs")?,
            },
        );
    }
    Ok(map)
}

fn compare_status(
    old: Option<&BenchComparePoint>,
    new: Option<&BenchComparePoint>,
    speedup: Option<f64>,
) -> String {
    match (old, new) {
        (None, Some(_)) => "new".to_string(),
        (Some(_), None) => "removed".to_string(),
        (None, None) => "missing".to_string(),
        (Some(old), Some(new)) if old.found_runs > 0 && new.found_runs == 0 => "lost".to_string(),
        (Some(old), Some(new)) if old.found_runs == 0 && new.found_runs > 0 => {
            "newly_found".to_string()
        }
        (Some(_), Some(_)) => {
            let speedup = speedup.unwrap_or(1.0);
            if speedup >= 1.10 {
                "faster".to_string()
            } else if speedup <= 0.90 {
                "slower".to_string()
            } else {
                "flat".to_string()
            }
        }
    }
}

fn benchmark_compare_rows(
    kind: &str,
    old: &BTreeMap<String, BenchComparePoint>,
    new: &BTreeMap<String, BenchComparePoint>,
) -> Vec<BenchCompareRow> {
    let names = old
        .keys()
        .chain(new.keys())
        .cloned()
        .collect::<BTreeSet<_>>();
    names
        .into_iter()
        .map(|name| {
            let old_point = old.get(&name);
            let new_point = new.get(&name);
            let speedup = match (old_point, new_point) {
                (Some(old), Some(new)) if new.mean_ms > 0.0 => Some(old.mean_ms / new.mean_ms),
                _ => None,
            };
            let delta_ms = match (old_point, new_point) {
                (Some(old), Some(new)) => Some(new.mean_ms - old.mean_ms),
                _ => None,
            };
            BenchCompareRow {
                kind: kind.to_string(),
                name,
                old_mean_ms: old_point.map(|point| point.mean_ms),
                new_mean_ms: new_point.map(|point| point.mean_ms),
                speedup,
                delta_ms,
                old_found_runs: old_point.map(|point| point.found_runs),
                old_runs: old_point.map(|point| point.runs),
                new_found_runs: new_point.map(|point| point.found_runs),
                new_runs: new_point.map(|point| point.runs),
                status: compare_status(old_point, new_point, speedup),
            }
        })
        .collect()
}

fn format_optional_f64(value: Option<f64>) -> String {
    value.map(|value| format!("{value:.3}")).unwrap_or_default()
}

fn format_optional_usize(value: Option<usize>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}

fn print_compare_rows(section: &str, rows: &[&BenchCompareRow]) {
    println!("# {section}");
    println!(
        "kind\tname\told_mean_ms\tnew_mean_ms\tspeedup\tdelta_ms\told_found\told_runs\tnew_found\tnew_runs\tstatus"
    );
    for row in rows {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.kind,
            row.name,
            format_optional_f64(row.old_mean_ms),
            format_optional_f64(row.new_mean_ms),
            format_optional_f64(row.speedup),
            format_optional_f64(row.delta_ms),
            format_optional_usize(row.old_found_runs),
            format_optional_usize(row.old_runs),
            format_optional_usize(row.new_found_runs),
            format_optional_usize(row.new_runs),
            row.status,
        );
    }
}

fn compare_row_json(row: &BenchCompareRow) -> Value {
    json!({
        "kind": &row.kind,
        "name": &row.name,
        "old_mean_ms": row.old_mean_ms,
        "new_mean_ms": row.new_mean_ms,
        "speedup": row.speedup,
        "delta_ms": row.delta_ms,
        "old_found_runs": row.old_found_runs,
        "old_runs": row.old_runs,
        "new_found_runs": row.new_found_runs,
        "new_runs": row.new_runs,
        "status": &row.status,
    })
}

fn cmd_bench_compare(args: &[String]) {
    let options = match parse_bench_compare_options(args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("{error}");
            print_bench_help();
            std::process::exit(1);
        }
    };
    let old_report = match read_benchmark_json(&options.old_path) {
        Ok(report) => report,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(1);
        }
    };
    let new_report = match read_benchmark_json(&options.new_path) {
        Ok(report) => report,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(1);
        }
    };

    let old_fixtures = match extract_benchmark_points(&old_report, "fixture_summaries", "slug") {
        Ok(points) => points,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(1);
        }
    };
    let new_fixtures = match extract_benchmark_points(&new_report, "fixture_summaries", "slug") {
        Ok(points) => points,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(1);
        }
    };
    let old_categories =
        match extract_benchmark_points(&old_report, "category_summaries", "category") {
            Ok(points) => points,
            Err(error) => {
                eprintln!("{error}");
                std::process::exit(1);
            }
        };
    let new_categories =
        match extract_benchmark_points(&new_report, "category_summaries", "category") {
            Ok(points) => points,
            Err(error) => {
                eprintln!("{error}");
                std::process::exit(1);
            }
        };

    let fixture_rows = benchmark_compare_rows("fixture", &old_fixtures, &new_fixtures);
    let category_rows = benchmark_compare_rows("category", &old_categories, &new_categories);
    let mut worst_regressions = fixture_rows
        .iter()
        .filter(|row| row.old_mean_ms.is_some() && row.new_mean_ms.is_some())
        .collect::<Vec<_>>();
    worst_regressions.sort_by(|a, b| {
        let a_ratio = a.new_mean_ms.unwrap_or(0.0) / a.old_mean_ms.unwrap_or(1.0);
        let b_ratio = b.new_mean_ms.unwrap_or(0.0) / b.old_mean_ms.unwrap_or(1.0);
        b_ratio.total_cmp(&a_ratio)
    });
    worst_regressions.truncate(options.top);

    if options.format == BenchOutputFormat::Json {
        let output = json!({
            "schema": "polynomial-tools.bench.compare.v1",
            "old": options.old_path.display().to_string(),
            "new": options.new_path.display().to_string(),
            "fixture_compare": fixture_rows.iter().map(compare_row_json).collect::<Vec<_>>(),
            "category_compare": category_rows.iter().map(compare_row_json).collect::<Vec<_>>(),
            "worst_regressions": worst_regressions
                .iter()
                .map(|row| compare_row_json(row))
                .collect::<Vec<_>>(),
        });
        println!("{}", serde_json::to_string_pretty(&output).unwrap());
    } else {
        print_compare_rows("fixture_compare", &fixture_rows.iter().collect::<Vec<_>>());
        print_compare_rows(
            "category_compare",
            &category_rows.iter().collect::<Vec<_>>(),
        );
        print_compare_rows("worst_regressions", &worst_regressions);
    }
}

fn parse_fixture_usize(flag: &str, value: Option<&str>) -> Result<usize, String> {
    let value = value.ok_or_else(|| format!("{flag} expects a value"))?;
    value
        .parse()
        .map_err(|_| format!("{flag} expects a nonnegative integer, got '{value}'"))
}

fn cmd_bench_recurrence_fixtures(args: &[String]) {
    let options = match parse_bench_recurrence_options(args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("{error}");
            print_bench_help();
            std::process::exit(1);
        }
    };
    let manifest = match parse_recurrence_fixture_manifest(&options.base) {
        Ok(manifest) => manifest,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(1);
        }
    };
    if options.format == BenchOutputFormat::Tsv {
        println!(
            "slug\trun\tfound\telapsed_ms\tcandidates\tunknowns\tweighted\tfit_rows\tverify_rows\trecurrence"
        );
    }
    let mut records = Vec::new();
    for fixture in manifest {
        if options
            .only
            .as_ref()
            .is_some_and(|needle| !fixture.slug.contains(needle))
        {
            continue;
        }
        let rows_path = options.base.join(&fixture.rows_file);
        let rows = match read_recurrence_fixture_rows(&rows_path) {
            Ok(rows) => rows,
            Err(error) => {
                eprintln!("{error}");
                std::process::exit(1);
            }
        };
        let search = match parse_fixture_search_options(
            &fixture.suggested_args,
            options.modular_prefilter,
        ) {
            Ok(search) => search,
            Err(error) => {
                eprintln!("{}: {error}", fixture.slug);
                std::process::exit(1);
            }
        };
        let categories = recurrence_fixture_categories(&fixture, &search);
        for run in 1..=options.repeat {
            let start = Instant::now();
            let result = find_recurrence_adaptive_rational(black_box(&rows), black_box(&search));
            let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
            if let Some(result) = result {
                let recurrence = result.recurrence.to_string();
                if options.format == BenchOutputFormat::Tsv {
                    println!(
                        "{}\t{}\ttrue\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}",
                        fixture.slug,
                        run,
                        elapsed_ms,
                        result.candidates_tried,
                        result.num_unknowns,
                        result.weighted_unknowns,
                        result.fit_polynomials,
                        result.verification_polynomials,
                        recurrence,
                    );
                }
                records.push(BenchRecurrenceRun {
                    slug: fixture.slug.clone(),
                    title: fixture.title.clone(),
                    description: fixture.description.clone(),
                    categories: categories.clone(),
                    run,
                    found: true,
                    elapsed_ms,
                    candidates: result.candidates_tried,
                    unknowns: result.num_unknowns,
                    weighted_unknowns: result.weighted_unknowns,
                    fit_rows: result.fit_polynomials,
                    verify_rows: result.verification_polynomials,
                    recurrence,
                    diagnostics: Some(result.diagnostics),
                });
            } else {
                if options.format == BenchOutputFormat::Tsv {
                    println!(
                        "{}\t{}\tfalse\t{:.3}\t0\t0\t0\t0\t0\t",
                        fixture.slug, run, elapsed_ms
                    );
                }
                records.push(BenchRecurrenceRun {
                    slug: fixture.slug.clone(),
                    title: fixture.title.clone(),
                    description: fixture.description.clone(),
                    categories: categories.clone(),
                    run,
                    found: false,
                    elapsed_ms,
                    candidates: 0,
                    unknowns: 0,
                    weighted_unknowns: 0,
                    fit_rows: 0,
                    verify_rows: 0,
                    recurrence: String::new(),
                    diagnostics: None,
                });
            }
        }
    }
    if options.format == BenchOutputFormat::Tsv && options.summary {
        print_recurrence_bench_summary(&records);
    }
    if let Err(error) = write_recurrence_bench_report(&options, &records) {
        eprintln!("{error}");
        std::process::exit(1);
    }
    if options.format == BenchOutputFormat::Json {
        println!(
            "{}",
            serde_json::to_string_pretty(&recurrence_bench_json(&options, &records)).unwrap()
        );
    }
}

#[derive(Debug)]
struct BenchInterlacingOptions {
    sequence: SequenceKind,
    max_n: usize,
    repeat: usize,
}

impl Default for BenchInterlacingOptions {
    fn default() -> Self {
        Self {
            sequence: SequenceKind::Eulerian,
            max_n: 20,
            repeat: 3,
        }
    }
}

fn parse_bench_interlacing_options(args: &[String]) -> Result<BenchInterlacingOptions, String> {
    let mut options = BenchInterlacingOptions::default();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--sequence" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--sequence expects a sequence name".to_string())?;
                options.sequence = SequenceKind::parse(value)
                    .ok_or_else(|| format!("unknown sequence: {value}"))?;
            }
            "--max-n" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--max-n expects a nonnegative integer".to_string())?;
                options.max_n = value
                    .parse()
                    .map_err(|_| format!("--max-n expects a nonnegative integer, got '{value}'"))?;
            }
            "--repeat" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--repeat expects a positive integer".to_string())?;
                options.repeat = value
                    .parse::<usize>()
                    .map_err(|_| format!("--repeat expects a positive integer, got '{value}'"))?
                    .max(1);
            }
            other => return Err(format!("unknown interlacing benchmark option: {other}")),
        }
        i += 1;
    }
    Ok(options)
}

fn cmd_bench_interlacing(args: &[String]) {
    let options = match parse_bench_interlacing_options(args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("{error}");
            print_bench_help();
            std::process::exit(1);
        }
    };
    let polynomials = options.sequence.polynomials(options.max_n);
    println!("sequence\tleft_index\tright_index\tdegree\trepeat\tavg_us\tresult");
    for right in 1..polynomials.len() {
        let left = right - 1;
        let p = strip_trailing_zeros_bigint(&polynomials[left]);
        let q = strip_trailing_zeros_bigint(&polynomials[right]);
        let mut result = None;
        let start = Instant::now();
        for _ in 0..options.repeat {
            result = check_interlacing_bigint_coeffs(black_box(p), black_box(q));
        }
        let avg_us = start.elapsed().as_secs_f64() * 1_000_000.0 / options.repeat as f64;
        println!(
            "{}\t{}\t{}\t{}\t{}\t{:.3}\t{}",
            options.sequence.name(),
            left,
            right,
            polynomial_degree_bigint(q),
            options.repeat,
            avg_us,
            json_bool_option(result),
        );
    }
}

fn cmd_resultant() {
    let polys = read_polys_bigint();
    if polys.len() < 2 {
        eprintln!("Need exactly two polynomials for resultant.");
        return;
    }
    let p = strip_trailing_zeros_bigint(&polys[0]);
    let q = strip_trailing_zeros_bigint(&polys[1]);
    let r = resultant_bigint_coeffs(p, q);
    println!(
        "Res({}, {}) = {}",
        format_poly_bigint_coeffs(p),
        format_poly_bigint_coeffs(q),
        r
    );
}

fn cmd_discriminant() {
    for coeffs in read_polys_bigint() {
        let c = strip_trailing_zeros_bigint(&coeffs);
        let d = discriminant_bigint_coeffs(c);
        println!("disc({}) = {}", format_poly_bigint_coeffs(c), d);
    }
}

fn cmd_hstar_to_ehrhart(args: &[String]) {
    let format = match OutputFormat::from_args(args) {
        Ok(format) => format,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let mut json_items = Vec::new();
    for (index, coeffs) in read_polys_bigint().into_iter().enumerate() {
        let c = strip_trailing_zeros_bigint(&coeffs);
        let ehrhart = hstar_to_ehrhart_bigint_coeffs(c);
        let display: Vec<String> = ehrhart.iter().map(|r| format!("{}", r)).collect();
        if format == OutputFormat::Json {
            json_items.push(format!(
                "{{\"index\":{},\"hstar\":{},\"ehrhart_coefficients\":{}}}",
                index,
                json_bigint_vec(c),
                json_string_vec(&display)
            ));
        } else {
            println!(
                "h*={} => L(n) coeffs: [{}]",
                format_poly_bigint_coeffs(c),
                display.join(", ")
            );
        }
    }
    if format == OutputFormat::Json {
        println!("{{\"items\":[{}]}}", json_items.join(","));
    }
}

fn cmd_ehrhart_to_hstar(args: &[String]) {
    let format = match OutputFormat::from_args(args) {
        Ok(format) => format,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let mut json_items = Vec::new();
    for (index, coeffs) in read_polys_rational().into_iter().enumerate() {
        let hstar = ehrhart_to_hstar_bigint(&coeffs);
        let display: Vec<String> = coeffs.iter().map(ToString::to_string).collect();
        if format == OutputFormat::Json {
            json_items.push(format!(
                "{{\"index\":{},\"ehrhart_coefficients\":{},\"hstar\":{}}}",
                index,
                json_string_vec(&display),
                json_bigint_vec(&hstar)
            ));
        } else {
            println!(
                "L(n) coeffs: [{}] => h*={}",
                display.join(", "),
                format_poly_bigint_coeffs(&hstar)
            );
        }
    }
    if format == OutputFormat::Json {
        println!("{{\"items\":[{}]}}", json_items.join(","));
    }
}

#[derive(Clone, Debug)]
struct FamilyCheckOptions {
    format: OutputFormat,
    recurrence: bool,
    require_real_rooted: bool,
    require_palindromic: bool,
    require_gamma_positive: bool,
    require_unimodal: bool,
    require_log_concave: bool,
    require_ultra_log_concave: bool,
    require_weak_interlacing: bool,
}

#[derive(Clone, Debug)]
struct ParseErrorReport {
    index: usize,
    error: String,
}

impl Default for FamilyCheckOptions {
    fn default() -> Self {
        Self {
            format: OutputFormat::Text,
            recurrence: false,
            require_real_rooted: false,
            require_palindromic: false,
            require_gamma_positive: false,
            require_unimodal: false,
            require_log_concave: false,
            require_ultra_log_concave: false,
            require_weak_interlacing: false,
        }
    }
}

#[derive(Clone, Debug)]
struct RecurrenceSummary {
    found: bool,
    recurrence: Option<String>,
    unknowns: Option<usize>,
    weighted_unknowns: Option<usize>,
    equations: Option<usize>,
    fit_polynomials: Option<usize>,
    verification_polynomials: Option<usize>,
    candidates_tried: Option<usize>,
    error: Option<String>,
}

fn parse_family_check_args(args: &[String]) -> Result<FamilyCheckOptions, String> {
    let mut options = FamilyCheckOptions::default();
    for arg in args {
        match arg.as_str() {
            "--json" => options.format = OutputFormat::Json,
            "--recurrence" => options.recurrence = true,
            "--require-real-rooted" => options.require_real_rooted = true,
            "--require-palindromic" => options.require_palindromic = true,
            "--require-gamma-positive" => options.require_gamma_positive = true,
            "--require-unimodal" => options.require_unimodal = true,
            "--require-log-concave" => options.require_log_concave = true,
            "--require-ultra-log-concave" => options.require_ultra_log_concave = true,
            "--require-weak-interlacing" => options.require_weak_interlacing = true,
            other => return Err(format!("unknown family-check option: {other}")),
        }
    }
    Ok(options)
}

fn first_family_failure(
    reports: &[PropertyReport],
    pairs: &[InterlacingReport],
    recurrence: Option<&RecurrenceSummary>,
    options: &FamilyCheckOptions,
) -> Option<String> {
    for report in reports {
        if options.require_real_rooted && !report.real_rooted {
            return Some(format!("polynomial {} is not real-rooted", report.index));
        }
        if options.require_palindromic && !report.palindromic {
            return Some(format!("polynomial {} is not palindromic", report.index));
        }
        if options.require_gamma_positive && !report.gamma_positive {
            return Some(format!("polynomial {} is not gamma-positive", report.index));
        }
        if options.require_unimodal && !report.unimodal {
            return Some(format!("polynomial {} is not unimodal", report.index));
        }
        if options.require_log_concave && !report.log_concave {
            return Some(format!("polynomial {} is not log-concave", report.index));
        }
        if options.require_ultra_log_concave && !report.ultra_log_concave {
            return Some(format!(
                "polynomial {} is not ultra-log-concave",
                report.index
            ));
        }
    }
    if options.require_weak_interlacing {
        for pair in pairs {
            if pair.weak != Some(true) {
                return Some(format!(
                    "polynomials {} and {} do not weakly interlace",
                    pair.left_index, pair.right_index
                ));
            }
        }
    }
    if options.recurrence {
        if let Some(summary) = recurrence {
            if !summary.found {
                return Some(
                    summary
                        .error
                        .clone()
                        .unwrap_or_else(|| "no recurrence found".to_string()),
                );
            }
        }
    }
    None
}

fn find_family_recurrence(polys: &[Vec<i64>]) -> RecurrenceSummary {
    if polys.len() < 3 {
        return RecurrenceSummary {
            found: false,
            recurrence: None,
            unknowns: None,
            weighted_unknowns: None,
            equations: None,
            fit_polynomials: None,
            verification_polynomials: None,
            candidates_tried: None,
            error: Some("need at least 3 polynomials".to_string()),
        };
    }
    let search = AdaptiveSearchOptions::default();
    match find_recurrence_adaptive(polys, &search) {
        Some(result) => RecurrenceSummary {
            found: true,
            recurrence: Some(result.recurrence.to_string()),
            unknowns: Some(result.num_unknowns),
            weighted_unknowns: Some(result.weighted_unknowns),
            equations: Some(result.num_equations),
            fit_polynomials: Some(result.fit_polynomials),
            verification_polynomials: Some(result.verification_polynomials),
            candidates_tried: Some(result.candidates_tried),
            error: None,
        },
        None => RecurrenceSummary {
            found: false,
            recurrence: None,
            unknowns: None,
            weighted_unknowns: None,
            equations: None,
            fit_polynomials: None,
            verification_polynomials: None,
            candidates_tried: None,
            error: Some("no recurrence found within the search bounds".to_string()),
        },
    }
}

fn family_recurrence_unavailable(error: impl Into<String>) -> RecurrenceSummary {
    RecurrenceSummary {
        found: false,
        recurrence: None,
        unknowns: None,
        weighted_unknowns: None,
        equations: None,
        fit_polynomials: None,
        verification_polynomials: None,
        candidates_tried: None,
        error: Some(error.into()),
    }
}

fn parse_error_report_json(report: &ParseErrorReport) -> String {
    format!(
        "{{\"index\":{},\"error\":{}}}",
        report.index,
        json_string(&report.error)
    )
}

fn recurrence_summary_json(summary: &RecurrenceSummary) -> String {
    format!(
        "{{\"found\":{},\"recurrence\":{},\"unknowns\":{},\"weighted_unknowns\":{},\
         \"equations\":{},\"fit_polynomials\":{},\"verification_polynomials\":{},\
         \"candidates_tried\":{},\"error\":{}}}",
        summary.found,
        summary
            .recurrence
            .as_ref()
            .map(|r| json_string(r))
            .unwrap_or_else(|| "null".to_string()),
        summary
            .unknowns
            .map(|n| n.to_string())
            .unwrap_or_else(|| "null".to_string()),
        summary
            .weighted_unknowns
            .map(|n| n.to_string())
            .unwrap_or_else(|| "null".to_string()),
        summary
            .equations
            .map(|n| n.to_string())
            .unwrap_or_else(|| "null".to_string()),
        summary
            .fit_polynomials
            .map(|n| n.to_string())
            .unwrap_or_else(|| "null".to_string()),
        summary
            .verification_polynomials
            .map(|n| n.to_string())
            .unwrap_or_else(|| "null".to_string()),
        summary
            .candidates_tried
            .map(|n| n.to_string())
            .unwrap_or_else(|| "null".to_string()),
        summary
            .error
            .as_ref()
            .map(|e| json_string(e))
            .unwrap_or_else(|| "null".to_string())
    )
}

fn cmd_family_check(args: &[String]) {
    let options = match parse_family_check_args(args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let mut parse_errors = Vec::new();
    let mut polys = Vec::new();
    for (index, item) in read_poly_parse_results_bigint().into_iter().enumerate() {
        match item {
            Ok(coeffs) => {
                polys.push((index, strip_trailing_zeros_bigint(&coeffs).to_vec()));
            }
            Err(error) => parse_errors.push(ParseErrorReport { index, error }),
        }
    }
    let reports = polys
        .iter()
        .map(|(index, coeffs)| property_report(*index, coeffs))
        .collect::<Vec<_>>();
    let i64_polys = polys
        .iter()
        .map(|(_, coeffs)| bigint_coeffs_to_i64(coeffs))
        .collect::<Vec<_>>();
    let pairs = if parse_errors.is_empty() {
        (0..polys.len().saturating_sub(1))
            .map(|pair_index| {
                let left_index = polys[pair_index].0;
                let right_index = polys[pair_index + 1].0;
                interlacing_report(
                    pair_index,
                    left_index,
                    right_index,
                    &polys[pair_index].1,
                    &polys[pair_index + 1].1,
                )
            })
            .collect::<Vec<_>>()
    } else {
        Vec::new()
    };
    let recurrence = options.recurrence.then(|| {
        if !parse_errors.is_empty() {
            return family_recurrence_unavailable(
                "recurrence skipped because one or more polynomials failed to parse",
            );
        }
        let Some(polys_i64) = i64_polys.iter().cloned().collect::<Option<Vec<_>>>() else {
            return family_recurrence_unavailable(
                "recurrence search currently requires coefficients that fit in i64",
            );
        };
        find_family_recurrence(&polys_i64)
    });
    let first_failure = if let Some(error) = parse_errors.first() {
        Some(format!(
            "polynomial {} failed to parse: {}",
            error.index, error.error
        ))
    } else {
        first_family_failure(&reports, &pairs, recurrence.as_ref(), &options)
    };
    let passed = first_failure.is_none();

    if options.format == OutputFormat::Json {
        let parse_error_items = parse_errors
            .iter()
            .map(parse_error_report_json)
            .collect::<Vec<_>>()
            .join(",");
        let items = reports
            .iter()
            .map(property_report_json)
            .collect::<Vec<_>>()
            .join(",");
        let pair_items = pairs
            .iter()
            .map(interlacing_report_json)
            .collect::<Vec<_>>()
            .join(",");
        let recurrence_json = recurrence
            .as_ref()
            .map(recurrence_summary_json)
            .unwrap_or_else(|| "null".to_string());
        let failure_json = first_failure
            .as_ref()
            .map(|failure| json_string(failure))
            .unwrap_or_else(|| "null".to_string());
        println!(
            "{{\"item_count\":{},\"all_required_checks_passed\":{},\
             \"first_failure\":{},\"parse_errors\":[{}],\"items\":[{}],\
             \"consecutive_pairs\":[{}],\"recurrence\":{}}}",
            reports.len(),
            passed,
            failure_json,
            parse_error_items,
            items,
            pair_items,
            recurrence_json
        );
        if !passed {
            std::process::exit(1);
        }
        return;
    }

    println!("Polynomial family check");
    println!("- Polynomials: {}", reports.len());
    println!(
        "- Required checks: {}",
        if passed { "passed" } else { "failed" }
    );
    if let Some(failure) = &first_failure {
        println!("- First failure: {failure}");
    }
    if !parse_errors.is_empty() {
        println!();
        println!("Parse errors:");
        for error in &parse_errors {
            println!("{}: {}", error.index, error.error);
        }
    }
    println!();
    println!("Properties:");
    for report in &reports {
        println!("{}: {}", report.index, property_labels(report).join(", "));
    }
    if !pairs.is_empty() {
        println!();
        println!("Consecutive interlacing:");
        for pair in &pairs {
            println!(
                "{}-{}: weak={}, strict={}, status={}",
                pair.left_index,
                pair.right_index,
                json_bool_option(pair.weak),
                json_bool_option(pair.strict),
                pair.status
            );
        }
    }
    if let Some(summary) = &recurrence {
        println!();
        println!("Recurrence:");
        if let Some(rec) = &summary.recurrence {
            println!("- Found: {rec}");
        } else if let Some(error) = &summary.error {
            println!("- Not found: {error}");
        } else {
            println!("- Not found");
        }
    }
    if !passed {
        std::process::exit(1);
    }
}

fn cmd_stapledon(args: &[String]) {
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_stapledon_help();
        return;
    }

    if args.len() != 1 {
        print_stapledon_help();
        return;
    }

    let n: usize = match args[0].parse() {
        Ok(n) => n,
        Err(_) => {
            eprintln!("Expected a nonnegative integer degree bound.");
            return;
        }
    };

    for coeffs in read_polys_bigint() {
        let c = strip_trailing_zeros_bigint(&coeffs);
        match stapledon_decomposition_bigint_coeffs(c, n) {
            Some((a, b)) => {
                println!(
                    "{} = {} + x ({})",
                    format_poly_bigint_coeffs(c),
                    format_poly_bigint_coeffs(&a),
                    format_poly_bigint_coeffs(&b),
                );
            }
            None => {
                eprintln!(
                    "{} has degree greater than the requested bound {}.",
                    format_poly_bigint_coeffs(c),
                    n,
                );
            }
        }
    }
}

fn bkw_scout_usage() {
    println!("Usage:");
    println!("  polytool bkw-scout [options]");
    println!("  polytool bkw-scout --help");
    println!();
    println!("Input symbol: z-coefficient polynomials in x, ascending z-degree.");
    println!("Example: F(x,z)=1-xz+z^2 is `1; -x; 1`.");
    println!();
    println!("Options:");
    println!("  --symbol <s>             Symbol coefficients, e.g. '1; -x; 1'");
    println!("                           If omitted, read the symbol from stdin.");
    println!("  --box <r0> <r1> <i0> <i1>");
    println!("                           Complex x rectangle (default: -3 3 -3 3)");
    println!("  --grid <n>               Use an n by n grid (default: 61)");
    println!("  --grid-re <n>            Number of real-axis grid samples");
    println!("  --grid-im <n>            Number of imaginary-axis grid samples");
    println!("  --top <n>                Number of candidates to print (default: 20)");
    println!("  --include-real-axis      Include samples with Im(x)=0");
    println!("  --min-imag <eps>         Minimum |Im(x)| unless real axis is included");
    println!("  --refine-steps <n>       Local coordinate-refinement steps (default: 8)");
    println!("  --no-refine              Disable local refinement");
    println!("  --tol <eps>              Durand--Kerner root tolerance");
    println!("  --max-iter <n>           Durand--Kerner max iterations");
    println!("  --format <text|json>     Output format (default: text)");
    println!("  --mathematica            Also print an exact Mathematica Reduce skeleton");
}

fn cmd_bkw_scout(args: &[String]) {
    let mut options = BkwScoutOptions::default();
    let mut symbol_input: Option<String> = None;
    let mut format = "text".to_string();
    let mut print_mathematica = false;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                bkw_scout_usage();
                return;
            }
            "--symbol" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--symbol requires an argument.");
                    return;
                }
                symbol_input = Some(args[i].clone());
            }
            "--box" => {
                if i + 4 >= args.len() {
                    eprintln!("--box requires four numeric arguments.");
                    return;
                }
                options.re_min = parse_f64_arg("--box re_min", &args[i + 1]);
                options.re_max = parse_f64_arg("--box re_max", &args[i + 2]);
                options.im_min = parse_f64_arg("--box im_min", &args[i + 3]);
                options.im_max = parse_f64_arg("--box im_max", &args[i + 4]);
                i += 4;
            }
            "--grid" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--grid requires a positive integer.");
                    return;
                }
                let grid = parse_usize_arg("--grid", &args[i]);
                options.grid_re = grid;
                options.grid_im = grid;
            }
            "--grid-re" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--grid-re requires a positive integer.");
                    return;
                }
                options.grid_re = parse_usize_arg("--grid-re", &args[i]);
            }
            "--grid-im" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--grid-im requires a positive integer.");
                    return;
                }
                options.grid_im = parse_usize_arg("--grid-im", &args[i]);
            }
            "--top" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--top requires a positive integer.");
                    return;
                }
                options.max_results = parse_usize_arg("--top", &args[i]);
            }
            "--include-real-axis" => {
                options.include_real_axis = true;
            }
            "--min-imag" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--min-imag requires a numeric argument.");
                    return;
                }
                options.min_imaginary_abs = parse_f64_arg("--min-imag", &args[i]);
            }
            "--refine-steps" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--refine-steps requires a nonnegative integer.");
                    return;
                }
                options.refine_steps = parse_usize_arg("--refine-steps", &args[i]);
            }
            "--no-refine" => {
                options.refine_steps = 0;
            }
            "--tol" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--tol requires a numeric argument.");
                    return;
                }
                options.root_options.tolerance = parse_f64_arg("--tol", &args[i]);
            }
            "--max-iter" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--max-iter requires a positive integer.");
                    return;
                }
                options.root_options.max_iterations = parse_usize_arg("--max-iter", &args[i]);
            }
            "--format" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--format requires text or json.");
                    return;
                }
                format = args[i].clone();
            }
            "--mathematica" => {
                print_mathematica = true;
            }
            other => {
                eprintln!("Unknown bkw-scout option: {}", other);
                bkw_scout_usage();
                return;
            }
        }
        i += 1;
    }

    if options.re_min > options.re_max || options.im_min > options.im_max {
        eprintln!("Invalid box: lower bounds must be <= upper bounds.");
        return;
    }
    if options.grid_re == 0 || options.grid_im == 0 || options.max_results == 0 {
        eprintln!("Grid sizes and --top must be positive.");
        return;
    }

    let input = match symbol_input {
        Some(s) => s,
        None => {
            let mut input = String::new();
            io::stdin().read_to_string(&mut input).unwrap();
            input
        }
    };
    let symbol = match BkwSymbol::parse_z_coefficient_symbol(&input) {
        Ok(symbol) => symbol,
        Err(e) => {
            eprintln!("{}", e);
            return;
        }
    };

    let candidates = symbol.scout_equal_modulus_locus(&options);
    match format.as_str() {
        "text" => print_bkw_text(&symbol, &options, &candidates),
        "json" => print_bkw_json(&symbol, &options, &candidates),
        other => {
            eprintln!("Unknown output format: {}", other);
            return;
        }
    }

    if print_mathematica {
        println!();
        println!("Mathematica exact equal-modulus skeleton:");
        println!("{}", symbol.mathematica_equal_modulus_query());
    }
}

fn parse_f64_arg(name: &str, value: &str) -> f64 {
    value.parse::<f64>().unwrap_or_else(|_| {
        eprintln!("{} expects a floating-point number, got '{}'.", name, value);
        std::process::exit(2);
    })
}

fn parse_usize_arg(name: &str, value: &str) -> usize {
    value.parse::<usize>().unwrap_or_else(|_| {
        eprintln!("{} expects a nonnegative integer, got '{}'.", name, value);
        std::process::exit(2);
    })
}

fn print_bkw_text(symbol: &BkwSymbol, options: &BkwScoutOptions, candidates: &[BkwScoutCandidate]) {
    println!("BKW equal-modulus scout");
    println!("Symbol: F(x,z) = {}", symbol.format_symbol());
    println!(
        "Box: Re(x) in [{:.6}, {:.6}], Im(x) in [{:.6}, {:.6}], grid {} x {}",
        options.re_min,
        options.re_max,
        options.im_min,
        options.im_max,
        options.grid_re,
        options.grid_im
    );
    println!(
        "Candidates: {} (ranked by dominant-root relative modulus gap)",
        candidates.len()
    );
    if candidates.is_empty() {
        println!("No candidate points found.  Try a larger box or finer grid.");
        return;
    }
    for (rank, candidate) in candidates.iter().enumerate() {
        let dominance = candidate
            .dominance_ratio
            .map(|r| format!("{:.6e}", r))
            .unwrap_or_else(|| "n/a".to_string());
        println!(
            "{}. x = {}   gap = {:.6e}   log_gap = {:.6e}   dominance(second/third) = {}",
            rank + 1,
            candidate.x,
            candidate.relative_modulus_gap,
            candidate.log_modulus_gap,
            dominance
        );
        println!(
            "   degree_at_x = {}, converged = {}, iterations = {}, residual = {:.3e}",
            candidate.z_degree_at_x,
            candidate.converged,
            candidate.iterations,
            candidate.root_residual
        );
        if let Some((left, right)) = candidate.tied_roots() {
            println!(
                "   tied roots: z{} = {} (|z|={:.12}), z{} = {} (|z|={:.12})",
                left.index, left.root, left.modulus, right.index, right.root, right.modulus
            );
        }
        let root_line = candidate
            .roots_by_modulus
            .iter()
            .map(|root| format!("z{}:{} |.|={:.6}", root.index, root.root, root.modulus))
            .collect::<Vec<_>>()
            .join("; ");
        println!("   roots by modulus: {}", root_line);
    }
    println!();
    println!(
        "{}",
        concat!(
            "Scout warning: this ranks numerical equal-modulus candidates only; ",
            "it does not prove BKW dominance, amplitude nonvanishing, ",
            "or eventual non-real-rootedness."
        )
    );
}

fn print_bkw_json(symbol: &BkwSymbol, options: &BkwScoutOptions, candidates: &[BkwScoutCandidate]) {
    println!("{{");
    println!("  \"symbol\": {:?},", symbol.format_symbol());
    println!(
        "  \"box\": {{\"re_min\": {}, \"re_max\": {}, \"im_min\": {}, \"im_max\": {}}},",
        options.re_min, options.re_max, options.im_min, options.im_max
    );
    println!(
        "  \"grid\": {{\"re\": {}, \"im\": {}}},",
        options.grid_re, options.grid_im
    );
    println!("  \"candidates\": [");
    for (i, candidate) in candidates.iter().enumerate() {
        let comma = if i + 1 == candidates.len() { "" } else { "," };
        println!("    {{");
        println!(
            "      \"x\": {{\"re\": {}, \"im\": {}}},",
            candidate.x.re, candidate.x.im
        );
        println!("      \"z_degree_at_x\": {},", candidate.z_degree_at_x);
        println!(
            "      \"relative_modulus_gap\": {},",
            candidate.relative_modulus_gap
        );
        println!("      \"log_modulus_gap\": {},", candidate.log_modulus_gap);
        match candidate.dominance_ratio {
            Some(ratio) if ratio.is_finite() => println!("      \"dominance_ratio\": {},", ratio),
            Some(_) => println!("      \"dominance_ratio\": \"Infinity\","),
            None => println!("      \"dominance_ratio\": null,"),
        }
        println!("      \"root_residual\": {},", candidate.root_residual);
        println!("      \"converged\": {},", candidate.converged);
        println!("      \"iterations\": {},", candidate.iterations);
        println!("      \"roots_by_modulus\": [");
        for (j, root) in candidate.roots_by_modulus.iter().enumerate() {
            let root_comma = if j + 1 == candidate.roots_by_modulus.len() {
                ""
            } else {
                ","
            };
            println!(
                "        {{\"index\": {}, \"re\": {}, \"im\": {}, \"modulus\": {}}}{}",
                root.index, root.root.re, root.root.im, root.modulus, root_comma
            );
        }
        println!("      ]");
        println!("    }}{}", comma);
    }
    println!("  ],");
    println!(
        "{}",
        concat!(
            "  \"warning\": \"numerical scout only; dominance and BKW ",
            "amplitude conditions are not certified\""
        )
    );
    println!("}}");
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        print_top_level_help();
        std::process::exit(1);
    }

    let cmd = &args[1];
    let rest = &args[2..];
    if is_help_arg(cmd) {
        if cmd == "help" && !rest.is_empty() {
            if is_help_arg(&rest[0]) {
                print_top_level_help();
            } else if !print_command_help(&rest[0]) {
                eprintln!("Unknown command: {}", rest[0]);
                std::process::exit(1);
            }
        } else {
            print_top_level_help();
        }
        return;
    }
    if rest.iter().any(|arg| is_help_arg(arg)) {
        if !print_command_help(cmd) {
            eprintln!("Unknown command: {}", cmd);
            std::process::exit(1);
        }
        return;
    }

    match cmd.as_str() {
        "real-rooted" => cmd_real_rooted(),
        "interlacing" => cmd_interlacing(rest),
        "interlacing-profile" => cmd_interlacing_profile(rest),
        "properties" => cmd_properties(rest),
        "gamma-expansion" | "gamma" => cmd_gamma_expansion(rest),
        "family-check" => cmd_family_check(rest),
        "sequence" => cmd_sequence(rest),
        "recurrence" => cmd_recurrence(rest),
        "recurrence-generate" => cmd_recurrence_generate(rest),
        "bench" => cmd_bench(rest),
        "bkw-scout" => cmd_bkw_scout(rest),
        "resultant" => cmd_resultant(),
        "discriminant" => cmd_discriminant(),
        "hstar-to-ehrhart" => cmd_hstar_to_ehrhart(rest),
        "ehrhart-to-hstar" => cmd_ehrhart_to_hstar(rest),
        "stapledon" => cmd_stapledon(rest),
        _ => {
            eprintln!("Unknown command: {}", cmd);
            std::process::exit(1);
        }
    }
}
