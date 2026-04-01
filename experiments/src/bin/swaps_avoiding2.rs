/// Long swaps on pattern-avoiding permutations.
///
/// long_swaps(σ) = |{i ∈ [n-1] : σ⁻¹(i) < σ⁻¹(i+1) - 1}|
/// i.e., value i appears left of value i+1 and they are NOT adjacent.
///
/// This is the original positional definition from the alternating paper,
/// which doesn't require checking membership in a set.
use combpoly::permutation::avoiding_permutations;
use combpoly::statistics::{compute, Stat};
use polynomial_tools::real_rootedness::{
    check_interlacing_sturm, check_weak_interlacing, format_poly, is_real_rooted,
};

fn build_poly(perms: &[Vec<u8>]) -> Vec<i64> {
    if perms.is_empty() {
        return vec![0];
    }
    let max_s = perms
        .iter()
        .map(|s| compute(s, Stat::Swaps))
        .max()
        .unwrap_or(0);
    let mut coeffs = vec![0i64; max_s + 1];
    for s in perms {
        coeffs[compute(s, Stat::Swaps)] += 1;
    }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
        coeffs.pop();
    }
    coeffs
}

fn is_palindromic(p: &[i64]) -> bool {
    let n = p.len();
    (0..n / 2).all(|i| p[i] == p[n - 1 - i])
}

fn check_il(f: &[i64], g: &[i64]) -> &'static str {
    if f.len() <= 1 || g.len() <= 1 {
        return "triv";
    }
    let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
    match check_interlacing_sturm(small, large) {
        Some(true) => "✓",
        Some(false) => match check_weak_interlacing(small, large) {
            Some(true) => "✓w",
            _ => "✗",
        },
        None => "?",
    }
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(11);

    let patterns: Vec<(&str, Vec<u8>)> = vec![
        ("123", vec![1, 2, 3]),
        ("132", vec![1, 3, 2]),
        ("213", vec![2, 1, 3]),
        ("231", vec![2, 3, 1]),
        ("312", vec![3, 1, 2]),
        ("321", vec![3, 2, 1]),
    ];

    for (name, pattern) in &patterns {
        println!("=== Av({}) with long_swaps ===", name);

        let mut polys: Vec<Vec<i64>> = Vec::new();

        for n in 1..=max_n {
            let perms = avoiding_permutations(n, &[pattern.clone()]);
            let poly = build_poly(&perms);
            let rr = if poly.len() <= 2 {
                true
            } else {
                is_real_rooted(&poly)
            };
            let pal = is_palindromic(&poly);
            let total: i64 = poly.iter().sum();

            println!(
                "  n={:>2}: {:<6} {:<55} {} {}",
                n,
                total,
                format_poly(&poly),
                if rr { "✓rr" } else { "✗rr" },
                if pal { "pal" } else { "   " },
            );
            polys.push(poly);
        }

        println!("  Interlacing:");
        for i in 0..polys.len() - 1 {
            let r = check_il(&polys[i], &polys[i + 1]);
            if r != "triv" {
                println!("    n={} ≪ n={}: {}", i + 1, i + 2, r);
            }
        }
        println!();
    }
}
