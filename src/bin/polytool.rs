//! CLI for polynomial-tools: real-rootedness, interlacing, recurrence search.
//!
//! Reads polynomials from stdin as comma-separated integer coefficients
//! in ascending degree order, one polynomial per line.

use num_bigint::BigInt;
use polynomial_tools::recurrence::BigRational as RecurrenceBigRational;
use polynomial_tools::recurrence::*;
use polynomial_tools::*;
use std::fs;
use std::io::{self, Read};

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

fn json_i64_vec(values: &[i64]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(ToString::to_string)
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
    println!("  Read dense coefficient lists from stdin, one polynomial per line.");
    println!("  Coefficients are in ascending degree order:");
    println!("    a_0, a_1, ..., a_d  represents  a_0 + a_1 t + ... + a_d t^d");
    println!("  Brackets and whitespace-separated coefficients are also accepted.");
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
}

fn print_stapledon_help() {
    println!("Usage:");
    println!("  polytool stapledon <n>");
    println!("  polytool stapledon --help");
    println!();
    println!("Compute the Stapledon decomposition with respect to degree bound n.");
    println!();
    print_coefficient_input_help();
}

fn print_command_help(command: &str) -> bool {
    match command {
        "real-rooted" => print_stdin_command_help(
            "real-rooted",
            "Check whether each input polynomial has only real roots.",
            &["Example:", "  echo '1, 3, 2' | polytool real-rooted"],
        ),
        "interlacing" => print_stdin_command_help(
            "interlacing",
            "Check strict and weak interlacing for consecutive input pairs.",
            &[
                "Options:",
                "  --json    Emit machine-readable JSON",
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
        "bkw-scout" => bkw_scout_usage(),
        "resultant" => print_stdin_command_help(
            "resultant",
            "Compute the resultant of the first two input polynomials.",
            &[
                "Example:",
                "  printf '1,0,1\\n-1,1\\n' | polytool resultant",
            ],
        ),
        "discriminant" => print_stdin_command_help(
            "discriminant",
            "Compute the discriminant of each input polynomial.",
            &["Example:", "  echo '1, 0, 1' | polytool discriminant"],
        ),
        "hstar-to-ehrhart" => print_stdin_command_help(
            "hstar-to-ehrhart",
            "Convert each h*-vector into Ehrhart polynomial coefficients.",
            &[
                "Options:",
                "  --json    Emit machine-readable JSON",
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
                "Example:",
                "  echo '1, 2, 2' | polytool ehrhart-to-hstar",
            ],
        ),
        "stapledon" => print_stapledon_help(),
        _ => return false,
    }
    true
}

fn read_polys() -> Vec<Vec<i64>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input).unwrap();
    parse_polynomials(&input)
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

fn strip_trailing_zeros(coeffs: &[i64]) -> &[i64] {
    let end = coeffs.iter().rposition(|&c| c != 0).map_or(0, |i| i + 1);
    &coeffs[..end]
}

fn polynomial_degree(coeffs: &[i64]) -> usize {
    coeffs.iter().rposition(|&c| c != 0).unwrap_or(0)
}

#[derive(Clone, Debug)]
struct PropertyReport {
    index: usize,
    coefficients: Vec<i64>,
    polynomial: String,
    degree: usize,
    real_rooted: bool,
    simple_roots: bool,
    palindromic: bool,
    gamma_positive: bool,
    gamma_coefficients: Option<Vec<i64>>,
    unimodal: bool,
    log_concave: bool,
    ultra_log_concave: bool,
}

fn property_report(index: usize, coeffs: &[i64]) -> PropertyReport {
    PropertyReport {
        index,
        coefficients: coeffs.to_vec(),
        polynomial: format_poly(coeffs),
        degree: polynomial_degree(coeffs),
        real_rooted: is_real_rooted(coeffs),
        simple_roots: has_simple_roots(coeffs),
        palindromic: is_palindromic_ignoring_initial_zeros(coeffs),
        gamma_positive: is_gamma_positive_ignoring_initial_zeros(coeffs),
        gamma_coefficients: gamma_coefficients_ignoring_initial_zeros(coeffs),
        unimodal: is_unimodal(coeffs),
        log_concave: is_log_concave(coeffs),
        ultra_log_concave: is_ultra_log_concave(coeffs),
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
        .map(|coeffs| json_i64_vec(coeffs))
        .unwrap_or_else(|| "null".to_string());
    format!(
        "{{\"index\":{},\"polynomial\":{},\"coefficients\":{},\"degree\":{},\
         \"real_rooted\":{},\"simple_roots\":{},\"palindromic\":{},\
         \"gamma_positive\":{},\"gamma_coefficients\":{},\"unimodal\":{},\
         \"log_concave\":{},\"ultra_log_concave\":{}}}",
        report.index,
        json_string(&report.polynomial),
        json_i64_vec(&report.coefficients),
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
    for coeffs in read_polys() {
        let c = strip_trailing_zeros(&coeffs);
        let rr = is_real_rooted(c);
        println!(
            "{}: {}",
            format_poly(c),
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
    p: &[i64],
    q: &[i64],
) -> InterlacingReport {
    let strict = check_interlacing(p, q);
    let weak = check_weak_interlacing(p, q);
    InterlacingReport {
        pair_index,
        left_index,
        right_index,
        p: format_poly(p),
        q: format_poly(q),
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
    polynomial: &[i64],
    previous: Vec<InterlacingReport>,
) -> InterlacingProfileReport {
    let previous_interlacing_indices = previous
        .iter()
        .filter(|report| interlacing_report_has_interlacing(report))
        .map(|report| report.left_index)
        .collect::<Vec<_>>();
    InterlacingProfileReport {
        index,
        polynomial: format_poly(polynomial),
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
    let polys = read_polys();
    if polys.len() < 2 {
        eprintln!("Need at least two polynomials for interlacing check.");
        return;
    }
    let reports = polys
        .windows(2)
        .enumerate()
        .map(|(i, pair)| {
            let p = strip_trailing_zeros(&pair[0]);
            let q = strip_trailing_zeros(&pair[1]);
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
    let polys = read_polys()
        .into_iter()
        .map(|coeffs| strip_trailing_zeros(&coeffs).to_vec())
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
    let reports = read_polys()
        .into_iter()
        .enumerate()
        .map(|(index, coeffs)| {
            let c = strip_trailing_zeros(&coeffs);
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
    for (index, coeffs) in read_polys().into_iter().enumerate() {
        let c = strip_trailing_zeros(&coeffs);
        match gamma_coefficients(c) {
            Some(gamma) => {
                let polynomial = format_poly(c);
                let expansion = format_gamma_expansion(&gamma, c.len().saturating_sub(1));
                if format == OutputFormat::Json {
                    json_items.push(format!(
                        "{{\"index\":{},\"ok\":true,\"polynomial\":{},\"coefficients\":{},\
                         \"degree\":{},\"gamma\":{},\"expansion\":{}}}",
                        index,
                        json_string(&polynomial),
                        json_i64_vec(c),
                        polynomial_degree(c),
                        json_i64_vec(&gamma),
                        json_string(&expansion)
                    ));
                } else {
                    println!("{polynomial}: gamma {:?}; expansion: {expansion}", gamma);
                }
            }
            None => {
                let polynomial = format_poly(c);
                let error = "not palindromic; no gamma expansion";
                if format == OutputFormat::Json {
                    json_items.push(format!(
                        "{{\"index\":{},\"ok\":false,\"polynomial\":{},\"coefficients\":{},\
                         \"degree\":{},\"error\":{}}}",
                        index,
                        json_string(&polynomial),
                        json_i64_vec(c),
                        polynomial_degree(c),
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

fn format_gamma_expansion(gamma: &[i64], degree: usize) -> String {
    let mut out = String::new();
    for (i, &coeff) in gamma.iter().enumerate() {
        if coeff == 0 {
            continue;
        }
        let abs_coeff = coeff.abs();
        let factor = gamma_basis_factor(i, degree.saturating_sub(2 * i));
        let body = if factor == "1" {
            abs_coeff.to_string()
        } else if abs_coeff == 1 {
            factor
        } else {
            format!("{abs_coeff} {factor}")
        };

        if out.is_empty() {
            if coeff < 0 {
                out.push('-');
            }
            out.push_str(&body);
        } else if coeff < 0 {
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

    fn polynomials(self, max_n: usize) -> Vec<Vec<i64>> {
        match self {
            Self::Eulerian => polynomial_tools::sequences::eulerian_polynomials(max_n),
            Self::Narayana => polynomial_tools::sequences::narayana_polynomials(max_n),
            Self::TypeBEulerian => polynomial_tools::sequences::type_b_eulerian_polynomials(max_n),
            Self::ChebyshevT => polynomial_tools::sequences::chebyshev_polynomials_t(max_n),
            Self::ChebyshevU => polynomial_tools::sequences::chebyshev_polynomials_u(max_n),
            Self::Hermite => polynomial_tools::sequences::hermite_polynomials(max_n),
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
                let c = strip_trailing_zeros(coeffs);
                format!(
                    "{{\"index\":{},\"label\":{},\"polynomial\":{},\"coefficients\":{},\
                     \"degree\":{}}}",
                    index,
                    json_string(&kind.label(index)),
                    json_string(&format_poly(c)),
                    json_i64_vec(c),
                    polynomial_degree(c)
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
        let c = strip_trailing_zeros(coeffs);
        println!("{} = {}", kind.label(index), format_poly(c));
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
                i += 1;
                format = match args[i].as_str() {
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
                i += 1;
                search.skip_prefix = args[i].parse().unwrap();
            }
            "--min-rec-len" => {
                i += 1;
                search.min_rec_len = args[i].parse().unwrap();
            }
            "--max-rec-len" => {
                i += 1;
                search.max_rec_len = args[i].parse().unwrap();
            }
            "--min-var-deg" => {
                i += 1;
                search.min_var_deg = args[i].parse().unwrap();
            }
            "--max-var-deg" => {
                i += 1;
                search.max_var_deg = args[i].parse().unwrap();
            }
            "--min-idx-deg" => {
                i += 1;
                search.min_idx_deg = args[i].parse().unwrap();
            }
            "--max-idx-deg" => {
                i += 1;
                search.max_idx_deg = args[i].parse().unwrap();
            }
            "--min-diff-deg" => {
                i += 1;
                search.min_diff_deg = args[i].parse().unwrap();
            }
            "--max-diff-deg" => {
                i += 1;
                search.max_diff_deg = args[i].parse().unwrap();
            }
            "--inhomogeneous" => {
                search.try_inhomogeneous = true;
            }
            "--min-inhomo-var-deg" => {
                i += 1;
                search.try_inhomogeneous = true;
                search.min_inhomo_var_deg = args[i].parse().unwrap();
            }
            "--max-inhomo-var-deg" => {
                i += 1;
                search.try_inhomogeneous = true;
                search.max_inhomo_var_deg = args[i].parse().unwrap();
            }
            "--min-inhomo-idx-deg" => {
                i += 1;
                search.try_inhomogeneous = true;
                search.min_inhomo_idx_deg = args[i].parse().unwrap();
            }
            "--max-inhomo-idx-deg" => {
                i += 1;
                search.try_inhomogeneous = true;
                search.max_inhomo_idx_deg = args[i].parse().unwrap();
            }
            "--denominator" | "--try-denominator" => {
                search.try_denominator = true;
            }
            "--alternating-sign" => {
                search.try_alternating_sign = true;
            }
            "--max-denom-var-deg" => {
                i += 1;
                search.try_denominator = true;
                search.max_denom_var_deg = args[i].parse().unwrap();
            }
            "--max-denom-idx-deg" => {
                i += 1;
                search.try_denominator = true;
                search.max_denom_idx_deg = args[i].parse().unwrap();
            }
            "--min-margin" => {
                i += 1;
                search.min_margin = args[i].parse().unwrap();
            }
            "--fit-extra-rows" => {
                i += 1;
                search.fit_extra_rows = args[i].parse().unwrap();
            }
            "--no-verify" => {
                search.no_verify = true;
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

fn cmd_resultant() {
    let polys = read_polys();
    if polys.len() < 2 {
        eprintln!("Need exactly two polynomials for resultant.");
        return;
    }
    let r = resultant(&polys[0], &polys[1]);
    println!(
        "Res({}, {}) = {}",
        format_poly(&polys[0]),
        format_poly(&polys[1]),
        r
    );
}

fn cmd_discriminant() {
    for coeffs in read_polys() {
        let c = strip_trailing_zeros(&coeffs);
        let d = discriminant(c);
        println!("disc({}) = {}", format_poly(c), d);
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
    for (index, coeffs) in read_polys().into_iter().enumerate() {
        let ehrhart = hstar_to_ehrhart(&coeffs);
        let display: Vec<String> = ehrhart.iter().map(|r| format!("{}", r)).collect();
        if format == OutputFormat::Json {
            json_items.push(format!(
                "{{\"index\":{},\"hstar\":{},\"ehrhart_coefficients\":{}}}",
                index,
                json_i64_vec(&coeffs),
                json_string_vec(&display)
            ));
        } else {
            println!(
                "h*={} => L(n) coeffs: [{}]",
                format_poly(&coeffs),
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
        let hstar = ehrhart_to_hstar(&coeffs);
        let display: Vec<String> = coeffs.iter().map(ToString::to_string).collect();
        if format == OutputFormat::Json {
            json_items.push(format!(
                "{{\"index\":{},\"ehrhart_coefficients\":{},\"hstar\":{}}}",
                index,
                json_string_vec(&display),
                json_i64_vec(&hstar)
            ));
        } else {
            println!(
                "L(n) coeffs: [{}] => h*={}",
                display.join(", "),
                format_poly(&hstar)
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
    let polys = read_polys()
        .into_iter()
        .map(|coeffs| strip_trailing_zeros(&coeffs).to_vec())
        .collect::<Vec<_>>();
    let reports = polys
        .iter()
        .enumerate()
        .map(|(index, coeffs)| property_report(index, coeffs))
        .collect::<Vec<_>>();
    let pairs = polys
        .windows(2)
        .enumerate()
        .map(|(index, pair)| interlacing_report(index, index, index + 1, &pair[0], &pair[1]))
        .collect::<Vec<_>>();
    let recurrence = options.recurrence.then(|| find_family_recurrence(&polys));
    let first_failure = first_family_failure(&reports, &pairs, recurrence.as_ref(), &options);
    let passed = first_failure.is_none();

    if options.format == OutputFormat::Json {
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
             \"first_failure\":{},\"items\":[{}],\"consecutive_pairs\":[{}],\
             \"recurrence\":{}}}",
            reports.len(),
            passed,
            failure_json,
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

    for coeffs in read_polys() {
        let c = strip_trailing_zeros(&coeffs);
        match stapledon_decomposition(c, n) {
            Some((a, b)) => {
                println!(
                    "{} = {} + x ({})",
                    format_poly(c),
                    format_poly(&a),
                    format_poly(&b),
                );
            }
            None => {
                eprintln!(
                    "{} has degree greater than the requested bound {}.",
                    format_poly(c),
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
