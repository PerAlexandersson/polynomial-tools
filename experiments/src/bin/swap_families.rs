/// Long swaps on derangements, involutions, and permutations by cycle count.
use combpoly::permutation::all_permutations;
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

fn is_derangement(sigma: &[u8]) -> bool {
    sigma.iter().enumerate().all(|(i, &v)| v != (i as u8 + 1))
}

fn is_involution(sigma: &[u8]) -> bool {
    // σ² = id: σ(σ(i)) = i for all i
    sigma
        .iter()
        .enumerate()
        .all(|(i, &v)| sigma[(v - 1) as usize] == (i as u8 + 1))
}

fn num_cycles(sigma: &[u8]) -> usize {
    let n = sigma.len();
    let mut visited = vec![false; n];
    let mut cycles = 0;
    for i in 0..n {
        if visited[i] {
            continue;
        }
        cycles += 1;
        let mut j = i;
        while !visited[j] {
            visited[j] = true;
            j = (sigma[j] - 1) as usize;
        }
    }
    cycles
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

fn is_palindromic(p: &[i64]) -> bool {
    let n = p.len();
    (0..n / 2).all(|i| p[i] == p[n - 1 - i])
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(10);

    // === Derangements ===
    println!("=== Derangements ===\n");
    let mut der_polys = Vec::new();
    for n in 1..=max_n {
        let all = all_permutations(n);
        let ders: Vec<Vec<u8>> = all.into_iter().filter(|s| is_derangement(s)).collect();
        let poly = build_poly(&ders);
        let rr = if poly.len() <= 2 {
            true
        } else {
            is_real_rooted(&poly)
        };
        let pal = is_palindromic(&poly);
        println!(
            "n={:>2}: {:<7} {:<60} {} {}",
            n,
            ders.len(),
            format_poly(&poly),
            if rr { "✓rr" } else { "✗rr" },
            if pal { "pal" } else { "   " },
        );
        der_polys.push(poly);
    }
    println!("\nInterlacing:");
    for i in 0..der_polys.len() - 1 {
        let r = check_il(&der_polys[i], &der_polys[i + 1]);
        if r != "triv" {
            println!("  n={} ≪ n={}: {}", i + 1, i + 2, r);
        }
    }

    // === Involutions ===
    println!("\n\n=== Involutions ===\n");
    let mut inv_polys = Vec::new();
    for n in 1..=max_n {
        let all = all_permutations(n);
        let invs: Vec<Vec<u8>> = all.into_iter().filter(|s| is_involution(s)).collect();
        let poly = build_poly(&invs);
        let rr = if poly.len() <= 2 {
            true
        } else {
            is_real_rooted(&poly)
        };
        let pal = is_palindromic(&poly);
        println!(
            "n={:>2}: {:<7} {:<60} {} {}",
            n,
            invs.len(),
            format_poly(&poly),
            if rr { "✓rr" } else { "✗rr" },
            if pal { "pal" } else { "   " },
        );
        inv_polys.push(poly);
    }
    println!("\nInterlacing:");
    for i in 0..inv_polys.len() - 1 {
        let r = check_il(&inv_polys[i], &inv_polys[i + 1]);
        if r != "triv" {
            println!("  n={} ≪ n={}: {}", i + 1, i + 2, r);
        }
    }

    // === Permutations by cycle count ===
    println!("\n\n=== By cycle count (n=8) ===\n");
    let n = 8u8.min(max_n);
    let all = all_permutations(n);
    let max_cyc = all.iter().map(|s| num_cycles(s)).max().unwrap_or(0);
    for c in 1..=max_cyc {
        let subset: Vec<Vec<u8>> = all.iter().filter(|s| num_cycles(s) == c).cloned().collect();
        let poly = build_poly(&subset);
        let rr = if poly.len() <= 2 {
            true
        } else {
            is_real_rooted(&poly)
        };
        println!(
            "  cyc={}: {:<6} {:<55} {}",
            c,
            subset.len(),
            format_poly(&poly),
            if rr { "✓rr" } else { "✗rr" },
        );
    }

    // Cycle count for multiple n values, fixed c
    println!("\n=== Fixed cycle count c, varying n ===\n");
    for c in 1..=3usize {
        println!("--- c={} ---", c);
        let mut cyc_polys = Vec::new();
        for n in 1..=max_n {
            let all = all_permutations(n);
            let subset: Vec<Vec<u8>> = all.into_iter().filter(|s| num_cycles(s) == c).collect();
            if subset.is_empty() {
                cyc_polys.push(vec![0]);
                continue;
            }
            let poly = build_poly(&subset);
            let rr = if poly.len() <= 2 {
                true
            } else {
                is_real_rooted(&poly)
            };
            println!(
                "  n={:>2}: {:<6} {:<55} {}",
                n,
                subset.len(),
                format_poly(&poly),
                if rr { "✓rr" } else { "✗rr" },
            );
            cyc_polys.push(poly);
        }
        println!("  Interlacing:");
        for i in 0..cyc_polys.len() - 1 {
            if cyc_polys[i] == vec![0] || cyc_polys[i + 1] == vec![0] {
                continue;
            }
            let r = check_il(&cyc_polys[i], &cyc_polys[i + 1]);
            if r != "triv" {
                println!("    n={} ≪ n={}: {}", i + 1, i + 2, r);
            }
        }
        println!();
    }
}
