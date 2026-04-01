/// Long swaps on all length-4 pattern-avoiding permutations.
/// Group by symmetry (complement + reverse) and identify distinct families.
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

fn complement(p: &[u8]) -> Vec<u8> {
    let n = p.len() as u8;
    p.iter().map(|&v| n + 1 - v).collect()
}

fn reverse(p: &[u8]) -> Vec<u8> {
    p.iter().rev().copied().collect()
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(10);

    // All 24 permutations of [4]
    let all4: Vec<Vec<u8>> = {
        let mut v = Vec::new();
        for a in 1..=4u8 {
            for b in 1..=4u8 {
                if b == a {
                    continue;
                }
                for c in 1..=4u8 {
                    if c == a || c == b {
                        continue;
                    }
                    for d in 1..=4u8 {
                        if d == a || d == b || d == c {
                            continue;
                        }
                        v.push(vec![a, b, c, d]);
                    }
                }
            }
        }
        v
    };

    // Group patterns by complement+reverse equivalence
    let mut groups: Vec<Vec<Vec<u8>>> = Vec::new();
    let mut assigned: Vec<bool> = vec![false; 24];

    for i in 0..24 {
        if assigned[i] {
            continue;
        }
        let p = &all4[i];
        let comp = complement(p);
        let rev = reverse(p);
        let comp_rev = complement(&rev);

        let mut group = vec![p.clone()];
        for q in [&comp, &rev, &comp_rev] {
            if !group.contains(q) {
                group.push(q.clone());
            }
        }
        for j in 0..24 {
            if group.contains(&all4[j]) {
                assigned[j] = true;
            }
        }
        groups.push(group);
    }

    println!(
        "=== {} symmetry classes of length-4 patterns ===\n",
        groups.len()
    );

    // For each group, compute polynomials and check properties
    // Collect results to compare polynomial families
    let mut all_results: Vec<(String, Vec<Vec<i64>>)> = Vec::new();

    for group in &groups {
        let repr = &group[0];
        let name: String = group
            .iter()
            .map(|p| p.iter().map(|v| v.to_string()).collect::<String>())
            .collect::<Vec<_>>()
            .join(" = ");

        let mut polys = Vec::new();
        println!("--- {} ---", name);

        for n in 1..=max_n {
            let perms = avoiding_permutations(n, &[repr.clone()]);
            let poly = build_poly(&perms);
            let rr = if poly.len() <= 2 {
                true
            } else {
                is_real_rooted(&poly)
            };
            let total: i64 = poly.iter().sum();
            println!(
                "  n={:>2}: {:<7} {:<65} {}",
                n,
                total,
                format_poly(&poly),
                if rr { "✓" } else { "✗" }
            );
            polys.push(poly);
        }

        // Interlacing
        let mut all_il = true;
        for i in 0..polys.len() - 1 {
            let (f, g) = (&polys[i], &polys[i + 1]);
            if f.len() <= 1 || g.len() <= 1 {
                continue;
            }
            let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
            let ok = match check_interlacing_sturm(small, large) {
                Some(true) => true,
                Some(false) => match check_weak_interlacing(small, large) {
                    Some(true) => true,
                    _ => false,
                },
                None => {
                    all_il = false;
                    false
                }
            };
            if !ok {
                all_il = false;
            }
        }
        println!(
            "  Interlacing: {}\n",
            if all_il { "ALL ✓" } else { "SOME ✗" }
        );

        all_results.push((name, polys));
    }

    // Check which groups give the same polynomials
    println!("=== Polynomial equivalence classes ===\n");
    let mut equiv: Vec<Vec<usize>> = Vec::new();
    let mut matched: Vec<bool> = vec![false; all_results.len()];

    for i in 0..all_results.len() {
        if matched[i] {
            continue;
        }
        let mut class = vec![i];
        for j in (i + 1)..all_results.len() {
            if matched[j] {
                continue;
            }
            if all_results[i].1 == all_results[j].1 {
                class.push(j);
                matched[j] = true;
            }
        }
        matched[i] = true;
        equiv.push(class);
    }

    for (idx, class) in equiv.iter().enumerate() {
        let names: Vec<&str> = class.iter().map(|&i| all_results[i].0.as_str()).collect();
        println!("Class {}: {}", idx + 1, names.join("  |  "));
    }
}
