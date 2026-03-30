/// Catalan backtrack: restrict search orders to 231-avoiding permutations.
///
/// Computes P_n(t) = Σ_{π ∈ Av_231(n)} t^{stat(backtrack_π(S_n))}
/// for various statistics, checks real-rootedness, searches for recurrences.
use combpoly::permutation::{all_permutations, backtrack_image, contains_pattern};
use polynomial_tools::real_rootedness as polynomial;
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};
use combpoly::statistics::{self, Stat};
use std::env;
use std::io::{self, Write};

fn main() {
    let max_n: u8 = env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    let pattern_231: Vec<u8> = vec![2, 3, 1];

    let stats: Vec<(&str, Stat)> = vec![
        ("des", Stat::Des),
        ("exc", Stat::Exc),
        ("peak", Stat::Peak),
        ("valley", Stat::Valley),
        ("lrmax", Stat::Lrmax),
        ("fix", Stat::Fix),
        ("cyc", Stat::Cyc),
        ("inv", Stat::Inv),
        ("maj", Stat::Maj),
    ];

    println!("Catalan backtrack: π ∈ Av_231(n), S = S_n");
    println!("Computing through n = {}\n", max_n);

    for &(stat_name, stat) in &stats {
        println!("{}", "=".repeat(60));
        println!("  stat = {}", stat_name);
        println!("{}", "=".repeat(60));

        let mut polys: Vec<Vec<i64>> = Vec::new();

        for n in 1..=max_n {
            let all = all_permutations(n);
            // Search orders: 231-avoiding permutations
            let pi_set: Vec<Vec<u8>> = all
                .iter()
                .filter(|p| !contains_pattern(p, &pattern_231))
                .cloned()
                .collect();

            // Constraint set: full S_n
            let s_set = &all;

            let mut coeffs: Vec<i64> = Vec::new();
            for pi in &pi_set {
                let sigma = backtrack_image(pi, s_set);
                let val = statistics::compute(&sigma, stat);
                while coeffs.len() <= val {
                    coeffs.push(0);
                }
                coeffs[val] += 1;
            }
            if coeffs.is_empty() {
                coeffs.push(0);
            }

            let rr = polynomial::is_real_rooted(&coeffs);
            let sum: i64 = coeffs.iter().sum();

            print!("  n={:>2}: [", n);
            for (j, &c) in coeffs.iter().enumerate() {
                if j > 0 {
                    print!(", ");
                }
                print!("{}", c);
            }
            print!("]  sum={}  RR={}", sum, if rr { "yes" } else { "NO" });
            println!();

            polys.push(coeffs);
            eprint!("\r  n={} done   ", n);
            io::stderr().flush().ok();
        }
        eprintln!();

        // Recurrence search: linear (no derivatives) first, then with derivatives
        println!("\n  Recurrence search:");
        for (label, opts) in &[
            // Linear recurrences (no derivatives)
            (
                "2-term linear",
                AdaptiveSearchOptions {
                    max_rec_len: 1,
                    max_var_deg: 3,
                    max_idx_deg: 3,
                    max_diff_deg: 0,
                    verbose: false,
                    ..Default::default()
                },
            ),
            (
                "3-term linear",
                AdaptiveSearchOptions {
                    max_rec_len: 2,
                    max_var_deg: 3,
                    max_idx_deg: 3,
                    max_diff_deg: 0,
                    verbose: false,
                    ..Default::default()
                },
            ),
            (
                "4-term linear",
                AdaptiveSearchOptions {
                    max_rec_len: 3,
                    max_var_deg: 3,
                    max_idx_deg: 3,
                    max_diff_deg: 0,
                    verbose: false,
                    ..Default::default()
                },
            ),
            // With derivatives
            (
                "2-term + deriv",
                AdaptiveSearchOptions {
                    max_rec_len: 1,
                    max_var_deg: 3,
                    max_idx_deg: 3,
                    max_diff_deg: 1,
                    verbose: false,
                    ..Default::default()
                },
            ),
            (
                "3-term + deriv",
                AdaptiveSearchOptions {
                    max_rec_len: 2,
                    max_var_deg: 3,
                    max_idx_deg: 3,
                    max_diff_deg: 1,
                    verbose: false,
                    ..Default::default()
                },
            ),
        ] {
            eprint!("\r  Searching {} ...         ", label);
            io::stderr().flush().ok();
            match find_recurrence_adaptive(&polys, opts) {
                Some(r) => println!("    {}: {}", label, r.recurrence),
                None => println!("    {}: (none found)", label),
            }
        }
        eprintln!();
        println!();
    }

    // Also try other Catalan families as search orders
    let other_patterns: Vec<(&str, Vec<u8>)> = vec![
        ("132", vec![1, 3, 2]),
        ("213", vec![2, 1, 3]),
        ("312", vec![3, 1, 2]),
        ("321", vec![3, 2, 1]),
    ];

    let core_stats: Vec<(&str, Stat)> = vec![
        ("des", Stat::Des),
        ("exc", Stat::Exc),
        ("peak", Stat::Peak),
    ];

    println!("\n{}", "#".repeat(60));
    println!("  Comparison: π ∈ Av_τ(n), S = S_n, for other τ");
    println!("{}\n", "#".repeat(60));

    for (pat_name, pat) in &other_patterns {
        for &(stat_name, stat) in &core_stats {
            println!("--- π ∈ Av_{}(n), stat = {} ---", pat_name, stat_name);

            let mut polys: Vec<Vec<i64>> = Vec::new();

            for n in 1..=max_n {
                let all = all_permutations(n);
                let pi_set: Vec<Vec<u8>> = all
                    .iter()
                    .filter(|p| !contains_pattern(p, pat))
                    .cloned()
                    .collect();
                let s_set = &all;

                let mut coeffs: Vec<i64> = Vec::new();
                for pi in &pi_set {
                    let sigma = backtrack_image(pi, s_set);
                    let val = statistics::compute(&sigma, stat);
                    while coeffs.len() <= val {
                        coeffs.push(0);
                    }
                    coeffs[val] += 1;
                }
                if coeffs.is_empty() {
                    coeffs.push(0);
                }

                let rr = polynomial::is_real_rooted(&coeffs);
                let sum: i64 = coeffs.iter().sum();

                print!("  n={:>2}: [", n);
                for (j, &c) in coeffs.iter().enumerate() {
                    if j > 0 {
                        print!(", ");
                    }
                    print!("{}", c);
                }
                print!("]  sum={}  RR={}", sum, if rr { "yes" } else { "NO" });
                println!();

                polys.push(coeffs);
            }

            // Recurrence search: linear first, then with derivatives
            for (label, opts) in &[
                (
                    "2-term linear",
                    AdaptiveSearchOptions {
                        max_rec_len: 1,
                        max_var_deg: 3,
                        max_idx_deg: 3,
                        max_diff_deg: 0,
                        verbose: false,
                        ..Default::default()
                    },
                ),
                (
                    "3-term linear",
                    AdaptiveSearchOptions {
                        max_rec_len: 2,
                        max_var_deg: 3,
                        max_idx_deg: 3,
                        max_diff_deg: 0,
                        verbose: false,
                        ..Default::default()
                    },
                ),
                (
                    "3-term + deriv",
                    AdaptiveSearchOptions {
                        max_rec_len: 2,
                        max_var_deg: 3,
                        max_idx_deg: 3,
                        max_diff_deg: 1,
                        verbose: false,
                        ..Default::default()
                    },
                ),
            ] {
                match find_recurrence_adaptive(&polys, opts) {
                    Some(r) => println!("    {}: {}", label, r.recurrence),
                    None => println!("    {}: (none found)", label),
                }
            }
            println!();
        }
    }
}
