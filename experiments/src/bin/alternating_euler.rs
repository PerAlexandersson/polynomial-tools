/// Alternating Eulerian polynomials H_n(t).
///
/// H_n(t) = sum_{σ ∈ A_n} t^{swaps(σ)}
///
/// where swaps(σ) = |{i ∈ [n-1] : σ^{-1}(i) < σ^{-1}(i+1) - 1}|
///
/// This statistic counts entries i in σ such that i appears somewhere
/// left of i+1, and swapping i and i+1 gives another alternating perm.
///
/// Conjectured: H_n(t) is real-rooted and H_n ≪ H_{n+1}.
use combpoly::permutation::alternating_permutations;
use combpoly::statistics::{compute, Stat};
use polynomial_tools::real_rootedness::{check_interlacing_sturm, format_poly, is_real_rooted};
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};
use std::time::Instant;

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(14);

    let mut polys: Vec<Vec<i64>> = Vec::new();

    println!("=== Alternating Eulerian polynomials H_n(t) ===\n");

    for n in 1..=max_n {
        let t = Instant::now();
        let alt = alternating_permutations(n);
        let count = alt.len();

        // Build generating polynomial
        let max_swap = alt
            .iter()
            .map(|s| compute(s, Stat::Swaps))
            .max()
            .unwrap_or(0);
        let mut coeffs = vec![0i64; max_swap + 1];
        for s in &alt {
            coeffs[compute(s, Stat::Swaps)] += 1;
        }

        let elapsed = t.elapsed();

        let rr = if coeffs.len() <= 2 {
            true
        } else {
            is_real_rooted(&coeffs)
        };

        println!(
            "H_{:<2}(t) = {:<50} |A_{}|={:<8} {:>8.2?} {}",
            n,
            format_poly(&coeffs),
            n,
            count,
            elapsed,
            if rr {
                "✓ real-rooted"
            } else {
                "✗ NOT real-rooted"
            }
        );

        polys.push(coeffs);
    }

    // Check interlacing
    println!("\n=== Interlacing check ===");
    for i in 0..(polys.len() - 1) {
        if polys[i].len() <= 1 || polys[i + 1].len() <= 1 {
            println!("H_{} ≪ H_{}: trivial", i + 1, i + 2);
            continue;
        }
        let (f, g) = if polys[i].len() <= polys[i + 1].len() {
            (&polys[i], &polys[i + 1])
        } else {
            (&polys[i + 1], &polys[i])
        };
        match check_interlacing_sturm(f, g) {
            Some(true) => println!("H_{} ≪ H_{}: ✓", i + 1, i + 2),
            Some(false) => println!("H_{} ≪ H_{}: ✗ FAILS", i + 1, i + 2),
            None => println!("H_{} ≪ H_{}: ? (could not determine)", i + 1, i + 2),
        }
    }

    // Search for recurrence
    println!("\n=== Recurrence search ===");
    let opts = AdaptiveSearchOptions::default();
    match find_recurrence_adaptive(&polys, &opts) {
        Some(res) => println!("{}", res.recurrence),
        None => println!("(no recurrence found)"),
    }

    // Print coefficient table for comparison with paper
    println!("\n=== Coefficient table ===");
    for (i, p) in polys.iter().enumerate() {
        let coeffs_str: Vec<String> = p.iter().map(|c| c.to_string()).collect();
        println!("n={:>2}: [{}]", i + 1, coeffs_str.join(", "));
    }
}
