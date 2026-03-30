/// Explore backtrack permutation polynomials from the paper.
///
/// For each combination of (pattern-avoidance set, statistic), computes the polynomial
///   Σ_{π ∈ S_n} t^{stat(backtrack_π(S))}
/// and searches for recurrences.
use combpoly::permutation::{
    all_permutations, backtrack_stat_poly, backtrack_stat_poly_subset, contains_pattern,
    is_alternating, is_derangement, is_grassmann, is_involution, is_vexillary, parse_sequence,
};
use polynomial_tools::real_rootedness as polynomial;
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};
use combpoly::statistics::Stat;
use std::env;
use std::io::{self, Write};

fn main() {
    let max_n: u8 = env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    let patterns = ["123", "132", "213", "231", "312", "321"];
    let stats = [
        ("des", Stat::Des),
        ("exc", Stat::Exc),
        ("peak", Stat::Peak),
        ("lrmax", Stat::Lrmax),
        ("inv", Stat::Inv),
        ("maj", Stat::Maj),
        ("fix", Stat::Fix),
        ("cyc", Stat::Cyc),
    ];

    // Run Catalan triangle example: sum over Av_132(n) of t^{exc(backtrack_π(Av_321))}
    let run_catalan = env::args().any(|a| a == "--catalan" || a == "-c");
    // Also run Grassmann permutations with exc
    let run_grassmann = env::args().any(|a| a == "--grassmann" || a == "-g");
    let run_recurrence = env::args().any(|a| a == "--recurrence" || a == "-r");
    let only_pattern: Option<String> = env::args()
        .position(|a| a == "--pattern" || a == "-p")
        .and_then(|i| env::args().nth(i + 1));
    let only_stat: Option<String> = env::args()
        .position(|a| a == "--stat" || a == "-s")
        .and_then(|i| env::args().nth(i + 1));

    for &(stat_name, stat) in &stats {
        if let Some(ref s) = only_stat {
            if s != stat_name {
                continue;
            }
        }

        for pat_str in &patterns {
            if let Some(ref p) = only_pattern {
                if p != pat_str {
                    continue;
                }
            }

            let pat = parse_sequence(pat_str);
            println!("\n=== Av_{} + {} (n=1..{}) ===", pat_str, stat_name, max_n);

            let mut polys: Vec<Vec<i64>> = Vec::new();

            for n in 1..=max_n {
                // Build S = Av_pat(n)
                let all = all_permutations(n);
                let s_set: Vec<Vec<u8>> = all
                    .into_iter()
                    .filter(|p| !contains_pattern(p, &pat))
                    .collect();

                let coeffs = backtrack_stat_poly(n, &s_set, stat);
                eprint!(
                    "\r  n={}: |Av|={}, poly = {}",
                    n,
                    s_set.len(),
                    polynomial::format_poly(&coeffs)
                );
                io::stderr().flush().ok();

                // Print as coefficient row
                print!("n={}: ", n);
                for (i, &c) in coeffs.iter().enumerate() {
                    if i > 0 {
                        print!(" ");
                    }
                    print!("{}", c);
                }
                println!();

                polys.push(coeffs);
            }
            eprintln!();

            // Print the triangle
            println!("--- Triangle ---");
            for (i, p) in polys.iter().enumerate() {
                print!("n={}: ", i + 1);
                for (j, &c) in p.iter().enumerate() {
                    if j > 0 {
                        print!("\t");
                    }
                    print!("{}", c);
                }
                println!();
            }

            // Check real-rootedness
            let mut first_non_rr = None;
            for (i, p) in polys.iter().enumerate() {
                if !polynomial::is_real_rooted(p) {
                    first_non_rr = Some(i + 1);
                    break;
                }
            }
            match first_non_rr {
                Some(n) => println!("Real-rooted: FAILS at n={}", n),
                None => println!("Real-rooted: yes for all n <= {}", max_n),
            }

            if run_recurrence && polys.len() >= 4 {
                search_recurrence(&polys, &format!("Av_{} + {}", pat_str, stat_name));
            }
        }
    }

    // Catalan triangle: Σ_{π ∈ Av_132(n)} t^{exc(backtrack_π(Av_321))}
    if run_catalan {
        use combpoly::permutation::backtrack_stat_poly_subset;
        let pat132 = parse_sequence("132");
        let pat321 = parse_sequence("321");

        println!("\n=== Catalan triangle: Av_132 π-set + Av_321 S-set + exc (n=1..{}) ===", max_n);
        let mut polys: Vec<Vec<i64>> = Vec::new();

        for n in 1..=max_n {
            let all = all_permutations(n);
            let pi_set: Vec<Vec<u8>> = all
                .iter()
                .filter(|p| !contains_pattern(p, &pat132))
                .cloned()
                .collect();
            let s_set: Vec<Vec<u8>> = all
                .into_iter()
                .filter(|p| !contains_pattern(p, &pat321))
                .collect();

            let coeffs = backtrack_stat_poly_subset(&pi_set, &s_set, Stat::Exc);
            print!("n={}: ", n);
            for (i, &c) in coeffs.iter().enumerate() {
                if i > 0 {
                    print!("\t");
                }
                print!("{}", c);
            }
            println!();
            polys.push(coeffs);
        }

        // Also try other combinations of (π-set pattern, S-set pattern, stat)
        let combos = [
            ("132", "321", "des"),
            ("213", "321", "exc"),
            ("321", "132", "exc"),
            ("132", "213", "des"),
        ];
        for &(pi_pat, s_pat, stat_name) in &combos {
            let stat = match stat_name {
                "des" => Stat::Des,
                "exc" => Stat::Exc,
                _ => continue,
            };
            let pi_p = parse_sequence(pi_pat);
            let s_p = parse_sequence(s_pat);

            println!(
                "\n=== Av_{} π-set + Av_{} S-set + {} (n=1..{}) ===",
                pi_pat, s_pat, stat_name, max_n
            );
            for n in 1..=max_n {
                let all = all_permutations(n);
                let pi_set: Vec<Vec<u8>> = all
                    .iter()
                    .filter(|p| !contains_pattern(p, &pi_p))
                    .cloned()
                    .collect();
                let s_set: Vec<Vec<u8>> = all
                    .into_iter()
                    .filter(|p| !contains_pattern(p, &s_p))
                    .collect();

                let coeffs = backtrack_stat_poly_subset(&pi_set, &s_set, stat);
                print!("n={}: ", n);
                for (i, &c) in coeffs.iter().enumerate() {
                    if i > 0 {
                        print!("\t");
                    }
                    print!("{}", c);
                }
                println!();
            }
        }
    }

    // Named S-set families: --sset involutions|derangements|alternating|vexillary
    let sset_name: Option<String> = env::args()
        .position(|a| a == "--sset")
        .and_then(|i| env::args().nth(i + 1));

    if let Some(ref sset) = sset_name {
        let filter_fn: Box<dyn Fn(&[u8]) -> bool> = match sset.as_str() {
            "involutions" => Box::new(|p: &[u8]| is_involution(p)),
            "derangements" => Box::new(|p: &[u8]| is_derangement(p)),
            "alternating" => Box::new(|p: &[u8]| is_alternating(p)),
            "vexillary" => Box::new(|p: &[u8]| is_vexillary(p)),
            "grassmann" => Box::new(|p: &[u8]| is_grassmann(p)),
            _ => {
                eprintln!("Unknown sset: {}. Use: involutions, derangements, alternating, vexillary, grassmann", sset);
                std::process::exit(1);
            }
        };

        for &(stat_name, stat) in &stats {
            if let Some(ref s) = only_stat {
                if s != stat_name {
                    continue;
                }
            }

            println!("\n=== {} + {} (n=1..{}) ===", sset, stat_name, max_n);
            let mut polys: Vec<Vec<i64>> = Vec::new();

            for n in 1..=max_n {
                let all = all_permutations(n);
                let s_set: Vec<Vec<u8>> = all.into_iter().filter(|p| filter_fn(p)).collect();
                if s_set.is_empty() {
                    println!("n={}: (empty S-set)", n);
                    continue;
                }

                let coeffs = backtrack_stat_poly(n, &s_set, stat);
                eprint!(
                    "\r  n={}: |{}|={}, poly = {}",
                    n,
                    sset,
                    s_set.len(),
                    polynomial::format_poly(&coeffs)
                );
                io::stderr().flush().ok();

                print!("n={}: ", n);
                for (i, &c) in coeffs.iter().enumerate() {
                    if i > 0 {
                        print!("\t");
                    }
                    print!("{}", c);
                }
                println!();

                polys.push(coeffs);
            }
            eprintln!();

            // Print the triangle
            println!("--- Triangle ---");
            for (i, p) in polys.iter().enumerate() {
                print!("n={}: ", i + 1);
                for (j, &c) in p.iter().enumerate() {
                    if j > 0 {
                        print!("\t");
                    }
                    print!("{}", c);
                }
                println!();
            }

            // Check real-rootedness
            let mut first_non_rr = None;
            for (i, p) in polys.iter().enumerate() {
                if !polynomial::is_real_rooted(p) {
                    first_non_rr = Some(i + 1);
                    break;
                }
            }
            match first_non_rr {
                Some(n) => println!("Real-rooted: FAILS at n={}", n),
                None => println!("Real-rooted: yes for all n <= {}", max_n),
            }

            if run_recurrence && polys.len() >= 4 {
                search_recurrence(&polys, &format!("{} + {}", sset, stat_name));
            }
        }
    }

    // Restricted π-set: --piset <pattern> --sset-pat <pattern> --combo-stat <stat>
    let piset_pat: Option<String> = env::args()
        .position(|a| a == "--piset")
        .and_then(|i| env::args().nth(i + 1));
    let sset_pat: Option<String> = env::args()
        .position(|a| a == "--sset-pat")
        .and_then(|i| env::args().nth(i + 1));
    let combo_stat: Option<String> = env::args()
        .position(|a| a == "--combo-stat")
        .and_then(|i| env::args().nth(i + 1));

    if let (Some(pi_pat_str), Some(s_pat_str)) = (piset_pat, sset_pat) {
        let pi_p = parse_sequence(&pi_pat_str);
        let s_p = parse_sequence(&s_pat_str);
        let stat_name = combo_stat.as_deref().unwrap_or("des");
        let stat = match stat_name {
            "des" => Stat::Des,
            "exc" => Stat::Exc,
            "peak" => Stat::Peak,
            "lrmax" => Stat::Lrmax,
            "inv" => Stat::Inv,
            "maj" => Stat::Maj,
            "fix" => Stat::Fix,
            "cyc" => Stat::Cyc,
            _ => {
                eprintln!("Unknown stat: {}", stat_name);
                std::process::exit(1);
            }
        };

        println!(
            "\n=== Av_{} π-set + Av_{} S-set + {} (n=1..{}) ===",
            pi_pat_str, s_pat_str, stat_name, max_n
        );
        let mut polys: Vec<Vec<i64>> = Vec::new();
        for n in 1..=max_n {
            let all = all_permutations(n);
            let pi_set: Vec<Vec<u8>> = all
                .iter()
                .filter(|p| !contains_pattern(p, &pi_p))
                .cloned()
                .collect();
            let s_set: Vec<Vec<u8>> = all
                .into_iter()
                .filter(|p| !contains_pattern(p, &s_p))
                .collect();

            let coeffs = backtrack_stat_poly_subset(&pi_set, &s_set, stat);
            print!("n={}: ", n);
            for (i, &c) in coeffs.iter().enumerate() {
                if i > 0 {
                    print!("\t");
                }
                print!("{}", c);
            }
            println!();
            polys.push(coeffs);
        }

        // Check real-rootedness
        let mut first_non_rr = None;
        for (i, p) in polys.iter().enumerate() {
            if !polynomial::is_real_rooted(p) {
                first_non_rr = Some(i + 1);
                break;
            }
        }
        match first_non_rr {
            Some(n) => println!("Real-rooted: FAILS at n={}", n),
            None => println!("Real-rooted: yes for all n <= {}", max_n),
        }

        if run_recurrence && polys.len() >= 4 {
            search_recurrence(
                &polys,
                &format!("Av_{} pi + Av_{} S + {}", pi_pat_str, s_pat_str, stat_name),
            );
        }
    }

    // Grassmann permutations with excedances
    if run_grassmann {
        println!("\n=== Grassmann + exc (n=1..{}) ===", max_n);
        let mut polys: Vec<Vec<i64>> = Vec::new();

        for n in 1..=max_n {
            let all = all_permutations(n);
            let s_set: Vec<Vec<u8>> = all.into_iter().filter(|p| is_grassmann(p)).collect();

            let coeffs = backtrack_stat_poly(n, &s_set, Stat::Exc);
            eprint!(
                "\r  n={}: |Grassmann|={}, poly = {}",
                n,
                s_set.len(),
                polynomial::format_poly(&coeffs)
            );
            io::stderr().flush().ok();

            print!("n={}: ", n);
            for (i, &c) in coeffs.iter().enumerate() {
                if i > 0 {
                    print!(" ");
                }
                print!("{}", c);
            }
            println!();

            polys.push(coeffs);
        }
        eprintln!();

        println!("--- Triangle ---");
        for (i, p) in polys.iter().enumerate() {
            print!("n={}: ", i + 1);
            for (j, &c) in p.iter().enumerate() {
                if j > 0 {
                    print!("\t");
                }
                print!("{}", c);
            }
            println!();
        }

        if run_recurrence && polys.len() >= 4 {
            search_recurrence(&polys, "Grassmann + exc");
        }
    }
}

fn search_recurrence(polys: &[Vec<i64>], label: &str) {
    let search = AdaptiveSearchOptions {
        max_rec_len: 4,
        max_var_deg: 3,
        max_idx_deg: 3,
        max_diff_deg: 2,
        verbose: false,
        ..Default::default()
    };

    eprintln!("Searching recurrence for {} ...", label);
    match find_recurrence_adaptive(polys, &search) {
        Some(result) => {
            println!(
                "RECURRENCE (rec_len={}, var_deg={}, idx_deg={}, diff_deg={}): {}",
                result.opts.rec_len,
                result.opts.var_deg,
                result.opts.idx_deg,
                result.opts.diff_deg,
                result.recurrence
            );
        }
        None => println!("No recurrence found."),
    }
}
