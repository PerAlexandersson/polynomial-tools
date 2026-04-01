/// Cross-pattern search depth: average depth to find an element of Av_τ₂
/// when the search order π is restricted to Av_τ₁ (or all of S_n).
///
/// Computes a 7×6 matrix:
///   rows = π-set (all S_n, Av_123, ..., Av_321)
///   cols = S-set (Av_123, ..., Av_321)
///   entry = average depth
///
/// Also checks: total depth over ALL π (not restricted), and restricted to Av_τ₁.
use combpoly::permutation::{
    all_permutations, avoiding_permutations, backtrack_image, parse_sequence,
};
use std::env;

fn lex_rank(perm: &[u8]) -> u64 {
    let n = perm.len();
    let mut rank: u64 = 0;
    let mut available: Vec<bool> = vec![true; n + 1];
    for i in 0..n {
        let val = perm[i] as usize;
        let smaller = (1..val).filter(|&v| available[v]).count() as u64;
        let factorial = (1..=(n - 1 - i) as u64).product::<u64>();
        rank += smaller * factorial;
        available[val] = false;
    }
    rank
}

fn compose(a: &[u8], b: &[u8]) -> Vec<u8> {
    b.iter().map(|&bi| a[(bi - 1) as usize]).collect()
}

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn avg_str(total: u64, count: u64) -> String {
    if count == 0 {
        return "N/A".to_string();
    }
    let g = gcd(total, count);
    let num = total / g;
    let den = count / g;
    if den == 1 {
        format!("{}", num)
    } else {
        format!("{}/{}", num, den)
    }
}

fn main() {
    let max_n: u8 = env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(7);

    let patterns = ["123", "132", "213", "231", "312", "321"];

    // Precompute avoidance sets for each n
    for n in 3..=max_n {
        println!("=== n = {} ===\n", n);

        let all_pi = all_permutations(n);
        let n_fact = all_pi.len() as u64;

        // Precompute S-sets
        let s_sets: Vec<Vec<Vec<u8>>> = patterns
            .iter()
            .map(|p| avoiding_permutations(n, &[parse_sequence(p)]))
            .collect();

        // Precompute pi-sets (including "all")
        let pi_sets: Vec<(&str, Vec<Vec<u8>>)> = {
            let mut v: Vec<(&str, Vec<Vec<u8>>)> = vec![("S_n", all_pi.clone())];
            for (i, p) in patterns.iter().enumerate() {
                v.push((p, s_sets[i].clone())); // Av_τ = same sets
            }
            v
        };

        // For each (pi-set, S-set) pair, compute total and average depth
        // Header
        print!("{:>10}", "π \\ S");
        for p in &patterns {
            print!("{:>14}", format!("Av_{}", p));
        }
        println!();
        print!("{:>10}", "");
        for _ in &patterns {
            print!("{:>14}", "----------");
        }
        println!();

        for (pi_name, pi_set) in &pi_sets {
            print!(
                "{:>10}",
                if *pi_name == "S_n" {
                    "S_n".to_string()
                } else {
                    format!("Av_{}", pi_name)
                }
            );

            for (s_idx, _) in patterns.iter().enumerate() {
                let s_set = &s_sets[s_idx];
                let mut total: u64 = 0;
                for pi in pi_set {
                    let sigma = backtrack_image(pi, s_set);
                    let sigma_pi = compose(&sigma, pi);
                    total += lex_rank(&sigma_pi) + 1;
                }
                let count = pi_set.len() as u64;
                print!("{:>14}", avg_str(total, count));
            }
            println!();
        }
        println!();

        // Also show the total depths (not averaged) for the S_n row
        print!("{:>10}", "total");
        for s_idx in 0..patterns.len() {
            let s_set = &s_sets[s_idx];
            let mut total: u64 = 0;
            for pi in &all_pi {
                let sigma = backtrack_image(pi, s_set);
                let sigma_pi = compose(&sigma, pi);
                total += lex_rank(&sigma_pi) + 1;
            }
            print!("{:>14}", total);
        }
        println!("\n");
    }

    // Summary: collect sequences for interesting cross-combinations
    println!("=== Sequences: average depth for cross-combinations ===\n");

    let interesting = [
        ("S_n", "213"),
        ("S_n", "123"),
        ("123", "312"),
        ("312", "123"),
        ("213", "132"),
        ("132", "213"),
        ("123", "213"),
        ("213", "123"),
        ("132", "132"),
        ("213", "213"),
        ("123", "123"),
        ("321", "321"),
    ];

    for &(pi_pat, s_pat) in &interesting {
        let pi_p = if pi_pat == "S_n" {
            None
        } else {
            Some(parse_sequence(pi_pat))
        };
        let s_p = parse_sequence(s_pat);

        let mut avgs: Vec<String> = Vec::new();
        let mut totals: Vec<u64> = Vec::new();

        for n in 1..=max_n {
            let all_pi = all_permutations(n);
            let s_set = avoiding_permutations(n, &[s_p.clone()]);

            let pi_set: Vec<Vec<u8>> = match &pi_p {
                None => all_pi,
                Some(p) => avoiding_permutations(n, &[p.clone()]),
            };

            let mut total: u64 = 0;
            for pi in &pi_set {
                let sigma = backtrack_image(pi, &s_set);
                let sigma_pi = compose(&sigma, pi);
                total += lex_rank(&sigma_pi) + 1;
            }

            let count = pi_set.len() as u64;
            totals.push(total);
            avgs.push(avg_str(total, count));
        }

        let pi_label = if pi_pat == "S_n" {
            "S_n".to_string()
        } else {
            format!("Av_{}", pi_pat)
        };
        println!(
            "π ∈ {:>6}, S = Av_{}: totals = {:?}",
            pi_label, s_pat, totals
        );
        println!("{:>24}  avgs   = {:?}", "", avgs);
        println!();
    }
}
