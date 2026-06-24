//! CLI for polynomial-tools: real-rootedness, interlacing, recurrence search.
//!
//! Reads polynomials from stdin as comma-separated integer coefficients
//! in ascending degree order, one polynomial per line.

use num_bigint::BigInt;
use polynomial_tools::recurrence::BigRational as RecurrenceBigRational;
use polynomial_tools::recurrence::*;
use polynomial_tools::*;
use std::io::{self, Read};

fn is_help_arg(arg: &str) -> bool {
    matches!(arg, "-h" | "--help" | "help")
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
    println!("  properties        Show real-rootedness, gamma, and related properties");
    println!("  recurrence        Search for a polynomial recurrence");
    println!("  bkw-scout         Scout BKW equal-modulus loci for a recurrence symbol");
    println!("  resultant         Compute the resultant of two polynomials");
    println!("  discriminant      Compute the discriminant of each polynomial");
    println!("  hstar-to-ehrhart  Convert h*-vectors to Ehrhart polynomials");
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
    println!("  --full-depth               Require all recurrence offsets to be used");
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
    println!("  --verbose                  Print candidate search details");
    println!("  -h, --help                 Print this help text");
    println!();
    println!("Example:");
    println!("  printf '1\\n1\\n2\\n3\\n5\\n8\\n' | polytool recurrence");
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
                "Example:",
                "  printf '2,-3,1\\n-1,1\\n' | polytool interlacing",
            ],
        ),
        "properties" => print_stdin_command_help(
            "properties",
            concat!(
                "Report real-rootedness, palindromicity, gamma-positivity, ",
                "log-concavity, and ultra-log-concavity."
            ),
            &["Example:", "  echo '1, 11, 11, 1' | polytool properties"],
        ),
        "recurrence" => print_recurrence_help(),
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
            &["Example:", "  echo '1, 2, 1' | polytool hstar-to-ehrhart"],
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
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_recurrence_help();
        return;
    }

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
        "interlacing" => cmd_interlacing(),
        "properties" => cmd_properties(),
        "recurrence" => cmd_recurrence(rest),
        "bkw-scout" => cmd_bkw_scout(rest),
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
