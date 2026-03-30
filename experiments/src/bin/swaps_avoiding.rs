/// Generalized swaps statistic on pattern-avoiding permutations.
///
/// For a set T ⊆ S_n and π ∈ T:
///   swaps_T(π) = |{i ∈ [n-1] : swapping values i and i+1 in π gives another element of T}|
///
/// Compute H_T(t) = Σ_{π ∈ T} t^{swaps_T(π)} for the six Catalan families
/// Av(τ) with τ of length 3.
use combpoly::permutation::avoiding_permutations;
use polynomial_tools::real_rootedness::{
    check_interlacing_sturm, format_poly, is_real_rooted,
};
use std::collections::HashSet;

/// Compute swaps_T(sigma): number of values i such that swapping
/// i and i+1 in sigma produces another permutation in T.
fn swaps_in_set(sigma: &[u8], set: &HashSet<Vec<u8>>) -> usize {
    let n = sigma.len();
    if n <= 1 { return 0; }
    let mut count = 0;
    for i in 1..n as u8 {
        // Try swapping values i and i+1
        let mut swapped = sigma.to_vec();
        let pos_i = sigma.iter().position(|&v| v == i).unwrap();
        let pos_ip1 = sigma.iter().position(|&v| v == i + 1).unwrap();
        swapped[pos_i] = i + 1;
        swapped[pos_ip1] = i;
        if set.contains(&swapped) {
            count += 1;
        }
    }
    count
}

fn build_poly(perms: &[Vec<u8>], set: &HashSet<Vec<u8>>) -> Vec<i64> {
    if perms.is_empty() { return vec![0]; }
    let max_s = perms.iter().map(|s| swaps_in_set(s, set)).max().unwrap_or(0);
    let mut coeffs = vec![0i64; max_s + 1];
    for s in perms {
        coeffs[swaps_in_set(s, set)] += 1;
    }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 { coeffs.pop(); }
    coeffs
}

fn is_palindromic(p: &[i64]) -> bool {
    let n = p.len();
    (0..n / 2).all(|i| p[i] == p[n - 1 - i])
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(10);

    let patterns: Vec<(&str, Vec<u8>)> = vec![
        ("123", vec![1, 2, 3]),
        ("132", vec![1, 3, 2]),
        ("213", vec![2, 1, 3]),
        ("231", vec![2, 3, 1]),
        ("312", vec![3, 1, 2]),
        ("321", vec![3, 2, 1]),
    ];

    for (name, pattern) in &patterns {
        println!("==============================");
        println!("=== Av({}) ===", name);
        println!("==============================\n");

        let mut polys: Vec<Vec<i64>> = Vec::new();

        for n in 1..=max_n {
            let perms = avoiding_permutations(n, &[pattern.clone()]);
            let set: HashSet<Vec<u8>> = perms.iter().cloned().collect();
            let poly = build_poly(&perms, &set);
            let rr = if poly.len() <= 2 { true } else { is_real_rooted(&poly) };
            let pal = is_palindromic(&poly);
            let total: i64 = poly.iter().sum();

            println!(
                "  n={:>2}: |Av|={:<6} {:<55} {} {}",
                n, total,
                format_poly(&poly),
                if rr { "✓rr" } else { "✗rr" },
                if pal { "pal" } else { "   " },
            );

            polys.push(poly);
        }

        // Check interlacing
        println!("\n  Interlacing:");
        for i in 0..polys.len() - 1 {
            let f = &polys[i];
            let g = &polys[i + 1];
            if f.len() <= 1 || g.len() <= 1 {
                continue;
            }
            let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
            match check_interlacing_sturm(small, large) {
                Some(true) => println!("    n={} ≪ n={}: ✓", i + 1, i + 2),
                Some(false) => println!("    n={} ≪ n={}: ✗", i + 1, i + 2),
                None => println!("    n={} ≪ n={}: ?", i + 1, i + 2),
            }
        }
        println!();
    }
}
