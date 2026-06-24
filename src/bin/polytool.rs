//! CLI for polynomial-tools: real-rootedness, interlacing, recurrence search.
//!
//! Reads polynomials from stdin as comma-separated integer coefficients
//! in ascending degree order, one polynomial per line.

use num_bigint::BigInt;
use polynomial_tools::recurrence::BigRational as RecurrenceBigRational;
use polynomial_tools::recurrence::*;
use polynomial_tools::*;
use std::io::{self, Read};

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

fn cmd_interlacing() {
    let polys = read_polys();
    if polys.len() < 2 {
        eprintln!("Need at least two polynomials for interlacing check.");
        return;
    }
    for pair in polys.windows(2) {
        let p = strip_trailing_zeros(&pair[0]);
        let q = strip_trailing_zeros(&pair[1]);
        let result = check_interlacing(p, q);
        let weak = check_weak_interlacing(p, q);
        let status = match (result, weak) {
            (Some(true), _) => "strictly interlace",
            (_, Some(true)) => "weakly interlace (shared roots)",
            (Some(false), _) => "do NOT interlace",
            _ => "incompatible degrees",
        };
        println!("{} & {}: {}", format_poly(p), format_poly(q), status);
    }
}

fn cmd_properties() {
    for coeffs in read_polys() {
        let c = strip_trailing_zeros(&coeffs);
        let mut props = Vec::new();

        if is_real_rooted(c) {
            props.push("real-rooted".to_string());
        }
        if is_palindromic_ignoring_initial_zeros(c) {
            props.push("palindromic".to_string());
        }
        if is_gamma_positive_ignoring_initial_zeros(c) {
            if let Some(gamma) = gamma_coefficients_ignoring_initial_zeros(c) {
                props.push(format!("gamma-positive {:?}", gamma));
            }
        }
        if is_log_concave(c) {
            props.push("log-concave".to_string());
        }
        if is_ultra_log_concave(c) {
            props.push("ultra-log-concave".to_string());
        }

        if props.is_empty() {
            props.push("(none)".to_string());
        }
        println!("{}: {}", format_poly(c), props.join(", "));
    }
}

fn cmd_recurrence(args: &[String]) {
    let mut search = AdaptiveSearchOptions::default();

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--skip-prefix" => {
                i += 1;
                search.skip_prefix = args[i].parse().unwrap();
            }
            "--full-depth" => {
                search.require_all_offsets = true;
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
            "--verbose" => {
                search.verbose = true;
            }
            _ => {}
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
            println!("{}", res.recurrence);
            eprintln!(
                "Found with {} unknowns, {} equations ({} candidates tried)",
                res.num_unknowns, res.num_equations, res.candidates_tried
            );
        }
        None => {
            eprintln!("No recurrence found within the search bounds.");
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

fn cmd_hstar_to_ehrhart() {
    for coeffs in read_polys() {
        let ehrhart = hstar_to_ehrhart(&coeffs);
        let display: Vec<String> = ehrhart.iter().map(|r| format!("{}", r)).collect();
        println!(
            "h*={} => L(n) coeffs: [{}]",
            format_poly(&coeffs),
            display.join(", ")
        );
    }
}

fn cmd_stapledon(args: &[String]) {
    if args.len() != 1 {
        eprintln!("Usage: polytool stapledon <n>");
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

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: polytool <command> [options]");
        eprintln!();
        eprintln!("Commands:");
        eprintln!("  real-rooted     Check real-rootedness of each polynomial");
        eprintln!("  interlacing     Check interlacing of consecutive polynomial pairs");
        eprintln!(
            "  properties      Show all properties (real-rooted, palindromic, gamma, log-concave)"
        );
        eprintln!("  recurrence      Search for a polynomial recurrence");
        eprintln!("  resultant       Compute resultant of two polynomials");
        eprintln!("  discriminant    Compute discriminant of each polynomial");
        eprintln!("  hstar-to-ehrhart  Convert h*-vector to Ehrhart polynomial");
        eprintln!(
            "  stapledon       Compute the Stapledon decomposition with respect to a bound n"
        );
        eprintln!();
        eprintln!("Input: polynomials as comma-separated integer coefficients (ascending degree),");
        eprintln!("       one per line on stdin. Lines starting with # are ignored.");
        std::process::exit(1);
    }

    let cmd = &args[1];
    let rest = &args[2..];
    match cmd.as_str() {
        "real-rooted" => cmd_real_rooted(),
        "interlacing" => cmd_interlacing(),
        "properties" => cmd_properties(),
        "recurrence" => cmd_recurrence(rest),
        "resultant" => cmd_resultant(),
        "discriminant" => cmd_discriminant(),
        "hstar-to-ehrhart" => cmd_hstar_to_ehrhart(),
        "stapledon" => cmd_stapledon(rest),
        _ => {
            eprintln!("Unknown command: {}", cmd);
            std::process::exit(1);
        }
    }
}
