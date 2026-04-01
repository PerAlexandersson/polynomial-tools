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

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(10);
    let pattern = vec![1u8, 2, 3, 4];
    let mut polys = Vec::new();
    println!("=== Av(1234) with long_swaps ===\n");
    for n in 1..=max_n {
        let perms = avoiding_permutations(n, &[pattern.clone()]);
        let poly = build_poly(&perms);
        let rr = if poly.len() <= 2 {
            true
        } else {
            is_real_rooted(&poly)
        };
        let total: i64 = poly.iter().sum();
        println!(
            "n={:>2}: {:<6} {:<65} {}",
            n,
            total,
            format_poly(&poly),
            if rr { "✓rr" } else { "✗rr" }
        );
        polys.push(poly);
    }
    println!("\nInterlacing:");
    for i in 0..polys.len() - 1 {
        let (f, g) = (&polys[i], &polys[i + 1]);
        if f.len() <= 1 || g.len() <= 1 {
            continue;
        }
        let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
        let r = match check_interlacing_sturm(small, large) {
            Some(true) => "✓",
            Some(false) => match check_weak_interlacing(small, large) {
                Some(true) => "✓w",
                _ => "✗",
            },
            None => "?",
        };
        println!("  n={} ≪ n={}: {}", i + 1, i + 2, r);
    }
}
