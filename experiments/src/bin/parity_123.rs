/// Count 123-avoiding permutations of [n] that preserve parity:
/// w(i) ≡ i (mod 2) for all i (even positions map to even values,
/// odd positions map to odd values).
///
/// Strategy: independently generate 123-avoiding perms on odd values
/// and even values, then combine and keep only those that are 123-avoiding
/// on the full permutation.
///
/// Computes generating polynomials for des, exc, peak, inv, maj, fix, cyc, lrmin.
/// Checks for recurrences via polynomial-tools.
use combpoly::permutation::avoiding_permutations;
use combpoly::polynomial_builder::build_generating_polynomial;
use combpoly::statistics::Stat;
use polynomial_tools::real_rootedness::{format_poly, is_real_rooted};
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};
use std::time::Instant;

/// Generate all parity-preserving 123-avoiding permutations of [n].
///
/// 1. Generate 123-avoiding perms of {1,3,5,...} (odd values, placed at odd positions)
/// 2. Generate 123-avoiding perms of {2,4,6,...} (even values, placed at even positions)
/// 3. Interleave each pair and keep only globally 123-avoiding results.
fn parity_123_avoiding(n: u8) -> Vec<Vec<u8>> {
    if n == 0 {
        return vec![vec![]];
    }

    let n_odd = (n + 1) / 2; // positions 1,3,5,... (count)
    let n_even = n / 2; // positions 2,4,6,... (count)

    // 123-avoiding perms of [1..n_odd], which we'll map to odd values
    let odd_perms = avoiding_permutations(n_odd, &[vec![1, 2, 3]]);
    // 123-avoiding perms of [1..n_even], which we'll map to even values
    let even_perms = if n_even > 0 {
        avoiding_permutations(n_even, &[vec![1, 2, 3]])
    } else {
        vec![vec![]]
    };

    let mut results = Vec::new();

    for op in &odd_perms {
        // Map to actual odd values: op[i] -> 2*op[i] - 1
        for ep in &even_perms {
            // Map to actual even values: ep[i] -> 2*ep[i]
            // Interleave: position 1 (odd) gets odd_perm[0], position 2 (even) gets even_perm[0], etc.
            let mut w = vec![0u8; n as usize];
            for i in 0..n_odd as usize {
                w[2 * i] = 2 * op[i] - 1;
            }
            for i in 0..n_even as usize {
                w[2 * i + 1] = 2 * ep[i];
            }

            // Check if the full permutation is 123-avoiding
            if is_123_avoiding(&w) {
                results.push(w);
            }
        }
    }

    results
}

/// Check if a permutation (1-indexed values in 0-indexed array) avoids 123.
fn is_123_avoiding(w: &[u8]) -> bool {
    let n = w.len();
    // Track minimum of suffix w[j+1..] that is > w[i] for some i < j
    // Equivalently: no i < j < k with w[i] < w[j] < w[k]
    // Use a stack-based O(n) check: w is 123-avoiding iff sortable by one stack
    // Simpler: track max of "second elements" seen so far
    // Actually, use the standard O(n log n) or O(n^2) check.
    // For our sizes, O(n^2) is fine.
    for i in 0..n {
        for j in (i + 1)..n {
            if w[j] > w[i] {
                // Check if any w[k] > w[j] for k > j
                for k in (j + 1)..n {
                    if w[k] > w[j] {
                        return false;
                    }
                }
            }
        }
    }
    true
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(20);

    let stats = [
        Stat::Des,
        Stat::Exc,
        Stat::Peak,
        Stat::Inv,
        Stat::Maj,
        Stat::Fix,
        Stat::Cyc,
        Stat::Lrmin,
        Stat::Lrmax,
        Stat::Valley,
    ];

    let mut all_polys: Vec<Vec<Vec<i64>>> = stats.iter().map(|_| Vec::new()).collect();
    let mut counts: Vec<usize> = Vec::new();

    println!("=== Parity-preserving 123-avoiding permutations ===\n");
    println!("{:>3}  {:>10}  {:>10}", "n", "count", "time");
    println!("{}", "-".repeat(28));

    for n in 1..=max_n {
        let t = Instant::now();
        let filtered = parity_123_avoiding(n);
        let elapsed = t.elapsed();
        let count = filtered.len();
        counts.push(count);

        for (i, &stat) in stats.iter().enumerate() {
            let poly = build_generating_polynomial(&filtered, stat);
            all_polys[i].push(poly);
        }

        println!("{:>3}  {:>10}  {:>10.2?}", n, count, elapsed);
    }

    // Print generating polynomials and check properties
    for (i, &stat) in stats.iter().enumerate() {
        println!("\n=== stat: {} ===", stat);
        for (j, poly) in all_polys[i].iter().enumerate() {
            let n = j + 1;
            let rr = if poly.len() <= 2 {
                true
            } else {
                is_real_rooted(poly)
            };
            println!(
                "  n={:>2}: {}{}",
                n,
                format_poly(poly),
                if rr { "" } else { "  [NOT real-rooted]" }
            );
        }

        // Search for recurrence
        println!("\n  Recurrence search:");
        let search_opts = AdaptiveSearchOptions::default();
        match find_recurrence_adaptive(&all_polys[i], &search_opts) {
            Some(res) => println!("  {}", res.recurrence),
            None => println!("  (no recurrence found)"),
        }
    }

    // Print count sequence for OEIS lookup
    println!("\n=== Count sequence ===");
    let seq: Vec<String> = counts.iter().map(|c| c.to_string()).collect();
    println!("{}", seq.join(", "));
}
