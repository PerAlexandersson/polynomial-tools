/// Refine H_n(t) by last entry: H_{n,k}(t) = sum_{σ ∈ A_n, σ(n)=k} t^{swaps(σ)}.
///
/// Look for matrix recurrences and interlacing structure.
use combpoly::permutation::alternating_permutations;
use combpoly::statistics::{compute, Stat};
use polynomial_tools::real_rootedness::{check_interlacing_sturm, format_poly, is_real_rooted};
use std::time::Instant;

/// Build generating polynomial from a subset of permutations.
fn build_poly(perms: &[&Vec<u8>]) -> Vec<i64> {
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
    // Trim trailing zeros
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
        coeffs.pop();
    }
    coeffs
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(12);

    for n in 1..=max_n {
        let t = Instant::now();
        let alt = alternating_permutations(n);
        let elapsed = t.elapsed();

        println!(
            "=== n={} ({} alternating perms, {:?}) ===",
            n,
            alt.len(),
            elapsed
        );

        // Refine by last entry
        let mut by_last: Vec<Vec<&Vec<u8>>> = vec![vec![]; n as usize + 1];
        for s in &alt {
            let last = *s.last().unwrap() as usize;
            by_last[last].push(s);
        }

        // Also refine by first entry
        let mut by_first: Vec<Vec<&Vec<u8>>> = vec![vec![]; n as usize + 1];
        for s in &alt {
            let first = s[0] as usize;
            by_first[first].push(s);
        }

        println!("  Refinement by LAST entry:");
        let mut last_polys: Vec<(u8, Vec<i64>)> = Vec::new();
        for k in 1..=n {
            if by_last[k as usize].is_empty() {
                continue;
            }
            let poly = build_poly(&by_last[k as usize]);
            let count: i64 = poly.iter().sum();
            let rr = if poly.len() <= 2 {
                true
            } else {
                is_real_rooted(&poly)
            };
            println!(
                "    k={:>2}: count={:<6} {} {}",
                k,
                count,
                format_poly(&poly),
                if rr { "✓" } else { "✗" }
            );
            last_polys.push((k, poly));
        }

        // Check interlacing among the refined polys
        if last_polys.len() >= 2 {
            println!("  Interlacing (last entry, consecutive k):");
            for i in 0..last_polys.len() - 1 {
                let (k1, ref f) = last_polys[i];
                let (k2, ref g) = last_polys[i + 1];
                if f.len() <= 1 || g.len() <= 1 {
                    continue;
                }
                let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
                match check_interlacing_sturm(small, large) {
                    Some(true) => println!("    H_{{n,{}}} ≪ H_{{n,{}}}: ✓", k1, k2),
                    Some(false) => println!("    H_{{n,{}}} ≪ H_{{n,{}}}: ✗", k1, k2),
                    None => println!("    H_{{n,{}}} ≪ H_{{n,{}}}: ?", k1, k2),
                }
            }
        }

        println!("  Refinement by FIRST entry:");
        let mut first_polys: Vec<(u8, Vec<i64>)> = Vec::new();
        for k in 1..=n {
            if by_first[k as usize].is_empty() {
                continue;
            }
            let poly = build_poly(&by_first[k as usize]);
            let count: i64 = poly.iter().sum();
            let rr = if poly.len() <= 2 {
                true
            } else {
                is_real_rooted(&poly)
            };
            println!(
                "    k={:>2}: count={:<6} {} {}",
                k,
                count,
                format_poly(&poly),
                if rr { "✓" } else { "✗" }
            );
            first_polys.push((k, poly));
        }

        // Also try refining by (first, last) pair for small n
        if n <= 8 {
            println!("  Refinement by (first, last):");
            for f in 1..=n {
                for l in 1..=n {
                    let subset: Vec<&Vec<u8>> = alt
                        .iter()
                        .filter(|s| s[0] == f && *s.last().unwrap() == l)
                        .collect();
                    if subset.is_empty() {
                        continue;
                    }
                    let poly = build_poly(&subset);
                    let count: i64 = poly.iter().sum();
                    println!(
                        "    ({},{}): count={:<4} {}",
                        f,
                        l,
                        count,
                        format_poly(&poly)
                    );
                }
            }
        }

        println!();
    }
}
