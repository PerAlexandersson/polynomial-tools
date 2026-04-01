/// Compute the average search depth to reach an element of S under backtracking.
///
/// For search order π, the backtrack enumeration lists all n! permutations
/// as σ_1∘π^{-1}, σ_2∘π^{-1}, ... where σ_j are in lex order.
/// The "depth" is the rank k of the first element in S.
///
/// We compute:
///   - Total depth = Σ_π depth(π, S)
///   - Average depth = Total / n!
///   - Distribution: how many π achieve each depth value
///
/// Also tries: generating polynomial Σ_π t^{depth(π,S)} (truncated).
use combpoly::permutation::{
    all_permutations, avoiding_permutations, backtrack_image, contains_pattern, is_derangement,
    is_involution, parse_sequence,
};
use std::collections::BTreeMap;
use std::env;

/// Lex rank of a permutation (0-based). Uses Lehmer code / factorial number system.
fn lex_rank(perm: &[u8]) -> u64 {
    let n = perm.len();
    let mut rank: u64 = 0;
    let mut available: Vec<bool> = vec![true; n + 1]; // 1-indexed

    for i in 0..n {
        // Count how many available values are less than perm[i]
        let val = perm[i] as usize;
        let smaller = (1..val).filter(|&v| available[v]).count() as u64;
        // Multiply by (n-1-i)!
        let factorial = (1..=(n - 1 - i) as u64).product::<u64>();
        rank += smaller * factorial;
        available[val] = false;
    }
    rank
}

/// Compose two permutations: (a ∘ b)(i) = a(b(i)). Both 1-indexed stored in 0-indexed arrays.
fn compose(a: &[u8], b: &[u8]) -> Vec<u8> {
    b.iter().map(|&bi| a[(bi - 1) as usize]).collect()
}

fn main() {
    let max_n: u8 = env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(7);

    let patterns = ["123", "132", "213", "231", "312", "321"];

    println!("=== Backtrack search depth ===\n");

    // Pattern-avoiding sets
    println!("--- Pattern-avoiding constraint sets ---\n");
    println!(
        "{:<12} {:>4} {:>10} {:>14} {:>10} {:>8}",
        "S-set", "n", "|S|", "total_depth", "avg_depth", "max_dep"
    );
    println!("{}", "-".repeat(65));

    for pat_str in &patterns {
        let pat = parse_sequence(pat_str);
        for n in 1..=max_n {
            let all_pi = all_permutations(n);
            let s_set = avoiding_permutations(n, &[pat.clone()]);

            let mut total_depth: u64 = 0;
            let mut max_depth: u64 = 0;

            for pi in &all_pi {
                let sigma = backtrack_image(pi, &s_set);
                let sigma_pi = compose(&sigma, pi);
                let depth = lex_rank(&sigma_pi) + 1; // 1-based
                total_depth += depth;
                max_depth = max_depth.max(depth);
            }

            let n_fact = (1..=n as u64).product::<u64>();
            let avg = total_depth as f64 / n_fact as f64;

            println!(
                "Av_{:<8} {:>4} {:>10} {:>14} {:>10.4} {:>8}",
                pat_str,
                n,
                s_set.len(),
                total_depth,
                avg,
                max_depth
            );
        }
        println!();
    }

    // Look at the total_depth and avg_depth sequences more carefully
    println!("\n--- Sequences: total depth and average depth ---\n");

    for pat_str in &patterns {
        let pat = parse_sequence(pat_str);
        let mut totals: Vec<u64> = Vec::new();
        let mut avgs: Vec<String> = Vec::new();

        for n in 1..=max_n {
            let all_pi = all_permutations(n);
            let s_set = avoiding_permutations(n, &[pat.clone()]);

            let mut total_depth: u64 = 0;
            for pi in &all_pi {
                let sigma = backtrack_image(pi, &s_set);
                let sigma_pi = compose(&sigma, pi);
                total_depth += lex_rank(&sigma_pi) + 1;
            }

            let n_fact = (1..=n as u64).product::<u64>();
            totals.push(total_depth);

            // Express average as fraction
            let g = gcd(total_depth, n_fact);
            avgs.push(format!("{}/{}", total_depth / g, n_fact / g));
        }

        println!("Av_{}: total = {:?}", pat_str, totals);
        println!("Av_{}: avg   = {:?}", pat_str, avgs);
        println!();
    }

    // Named families
    println!("\n--- Named constraint sets ---\n");

    let named: Vec<(&str, Box<dyn Fn(&[u8]) -> bool>)> = vec![
        ("derangements", Box::new(|p: &[u8]| is_derangement(p))),
        ("involutions", Box::new(|p: &[u8]| is_involution(p))),
    ];

    for (name, filter_fn) in &named {
        let mut totals: Vec<u64> = Vec::new();
        let mut avgs: Vec<String> = Vec::new();

        for n in 1..=max_n {
            let all_pi = all_permutations(n);
            let s_set: Vec<Vec<u8>> = all_pi.iter().filter(|p| filter_fn(p)).cloned().collect();
            if s_set.is_empty() {
                totals.push(0);
                avgs.push("N/A".to_string());
                continue;
            }

            let mut total_depth: u64 = 0;
            for pi in &all_pi {
                let sigma = backtrack_image(pi, &s_set);
                let sigma_pi = compose(&sigma, pi);
                total_depth += lex_rank(&sigma_pi) + 1;
            }

            let n_fact = (1..=n as u64).product::<u64>();
            totals.push(total_depth);
            let g = gcd(total_depth, n_fact);
            avgs.push(format!("{}/{}", total_depth / g, n_fact / g));
        }

        println!("{}: total = {:?}", name, totals);
        println!("{}: avg   = {:?}", name, avgs);
        println!();
    }

    // Depth distribution for select cases
    println!("\n--- Depth distribution for select cases ---\n");

    for pat_str in &["213", "132"] {
        let pat = parse_sequence(pat_str);
        for n in [3, 4, 5, 6].iter().copied().filter(|&n| n <= max_n) {
            let all_pi = all_permutations(n);
            let s_set = avoiding_permutations(n, &[pat.clone()]);

            let mut dist: BTreeMap<u64, u64> = BTreeMap::new();
            for pi in &all_pi {
                let sigma = backtrack_image(pi, &s_set);
                let sigma_pi = compose(&sigma, pi);
                let depth = lex_rank(&sigma_pi) + 1;
                *dist.entry(depth).or_default() += 1;
            }

            println!("Av_{}, n={}: depth distribution (depth: count)", pat_str, n);
            for (d, c) in &dist {
                print!("  {}:{}", d, c);
            }
            println!("\n");
        }
    }
}

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}
