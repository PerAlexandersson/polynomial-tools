/// k-alternating Eulerian polynomials.
///
/// For each k, compute H^{(k)}_n(t) = sum_{σ ∈ A^{(k)}_n} t^{swaps(σ)}
/// where A^{(k)}_n = k-alternating permutations (descents at positions divisible by k).
///
/// Check real-rootedness, palindromicity, and interlacing.
use combpoly::permutation::k_alternating_permutations;
use combpoly::statistics::{compute, Stat};
use polynomial_tools::real_rootedness::{
    check_interlacing_sturm, format_poly, is_real_rooted,
};
use std::time::Instant;

fn build_poly(perms: &[Vec<u8>]) -> Vec<i64> {
    if perms.is_empty() { return vec![0]; }
    let max_s = perms.iter().map(|s| compute(s, Stat::Swaps)).max().unwrap_or(0);
    let mut coeffs = vec![0i64; max_s + 1];
    for s in perms { coeffs[compute(s, Stat::Swaps)] += 1; }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 { coeffs.pop(); }
    coeffs
}

fn is_palindromic(p: &[i64]) -> bool {
    let n = p.len();
    (0..n / 2).all(|i| p[i] == p[n - 1 - i])
}

fn main() {
    let max_k: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(5);

    let max_n: u8 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(12);

    for k in 2..=max_k {
        println!("==============================");
        println!("=== k={} alternating ===", k);
        println!("==============================\n");

        let mut polys: Vec<Vec<i64>> = Vec::new();

        for n in 1..=max_n {
            let t = Instant::now();
            let perms = k_alternating_permutations(n, k);
            let elapsed = t.elapsed();
            let count = perms.len();

            if count == 0 {
                println!("n={:>2}: (empty)", n);
                polys.push(vec![0]);
                continue;
            }

            let poly = build_poly(&perms);
            let rr = if poly.len() <= 2 { true } else { is_real_rooted(&poly) };
            let pal = is_palindromic(&poly);

            println!(
                "n={:>2}: count={:<8} {:<60} {} {} {:>8.2?}",
                n, count,
                format_poly(&poly),
                if rr { "✓rr" } else { "✗rr" },
                if pal { "pal" } else { "   " },
                elapsed,
            );

            polys.push(poly);
        }

        // Check interlacing of consecutive H^(k)_n
        println!("\nInterlacing H^({})_n ≪ H^({})_{{n+1}}:", k, k);
        for i in 0..polys.len() - 1 {
            let f = &polys[i];
            let g = &polys[i + 1];
            if f.len() <= 1 || g.len() <= 1 {
                println!("  n={} → n={}: trivial", i + 1, i + 2);
                continue;
            }
            let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
            match check_interlacing_sturm(small, large) {
                Some(true) => println!("  n={} ≪ n={}: ✓", i + 1, i + 2),
                Some(false) => println!("  n={} ≪ n={}: ✗", i + 1, i + 2),
                None => println!("  n={} ≪ n={}: ?", i + 1, i + 2),
            }
        }

        // Print count sequence
        let counts: Vec<String> = (1..=max_n)
            .map(|n| {
                let perms = k_alternating_permutations(n, k);
                perms.len().to_string()
            })
            .collect();
        println!("\nCount sequence (k={}): {}\n", k, counts.join(", "));
    }
}
