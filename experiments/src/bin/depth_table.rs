/// Print a clean table of ALL total search depth sequences (= n! * avg).
/// For all 6×6 cross-pattern combinations (length 3) plus S_n rows,
/// plus select cross-length combinations.
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

fn total_depth(pi_set: &[Vec<u8>], s_set: &[Vec<u8>]) -> u64 {
    let mut total: u64 = 0;
    for pi in pi_set {
        let sigma = backtrack_image(pi, s_set);
        let sigma_pi = compose(&sigma, pi);
        total += lex_rank(&sigma_pi) + 1;
    }
    total
}

fn main() {
    let max_n: u8 = env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(8);

    let pats3 = ["123", "132", "213", "231", "312", "321"];

    // Precompute for each n
    println!("Total search depth D_n(Π, S) = Σ_{{π∈Π}} depth(π, S)\n");
    println!("Each row: sequence for n = 1, 2, ..., {}\n", max_n);

    // Header
    println!(
        "{:<28} {}",
        "Combination",
        (1..=max_n)
            .map(|n| format!("{:>12}", n))
            .collect::<Vec<_>>()
            .join("")
    );
    println!("{}", "-".repeat(28 + 12 * max_n as usize));

    // S_n rows (length-3 S-sets)
    for s_pat in &pats3 {
        let s_p = parse_sequence(s_pat);
        let mut seq = Vec::new();
        for n in 1..=max_n {
            let all_pi = all_permutations(n);
            let s_set = avoiding_permutations(n, &[s_p.clone()]);
            seq.push(total_depth(&all_pi, &s_set));
        }
        println!(
            "{:<28} {}",
            format!("S_n → Av_{}", s_pat),
            seq.iter()
                .map(|v| format!("{:>12}", v))
                .collect::<Vec<_>>()
                .join("")
        );
    }

    println!("{}", "-".repeat(28 + 12 * max_n as usize));

    // 6×6 cross-pattern (length 3)
    for pi_pat in &pats3 {
        let pi_p = parse_sequence(pi_pat);
        for s_pat in &pats3 {
            let s_p = parse_sequence(s_pat);
            let mut seq = Vec::new();
            for n in 1..=max_n {
                let pi_set = avoiding_permutations(n, &[pi_p.clone()]);
                let s_set = avoiding_permutations(n, &[s_p.clone()]);
                seq.push(total_depth(&pi_set, &s_set));
            }
            println!(
                "{:<28} {}",
                format!("Av_{} → Av_{}", pi_pat, s_pat),
                seq.iter()
                    .map(|v| format!("{:>12}", v))
                    .collect::<Vec<_>>()
                    .join("")
            );
        }
    }

    println!("{}", "-".repeat(28 + 12 * max_n as usize));

    // Cross-length: S_n → Av_τ₄
    let pats4_select = ["1234", "1243", "1324", "1342"];
    for s_pat in &pats4_select {
        let s_p = parse_sequence(s_pat);
        let mut seq = Vec::new();
        for n in 1..=max_n {
            let all_pi = all_permutations(n);
            let s_set = avoiding_permutations(n, &[s_p.clone()]);
            seq.push(total_depth(&all_pi, &s_set));
        }
        println!(
            "{:<28} {}",
            format!("S_n → Av_{}", s_pat),
            seq.iter()
                .map(|v| format!("{:>12}", v))
                .collect::<Vec<_>>()
                .join("")
        );
    }

    println!("{}", "-".repeat(28 + 12 * max_n as usize));

    // Cross-length: Av_τ₄ → Av_τ₃
    let cross = [
        ("1234", "123"),
        ("1234", "132"),
        ("1234", "213"),
        ("1243", "123"),
        ("1243", "132"),
        ("1243", "213"),
        ("2134", "123"),
        ("2134", "132"),
        ("2134", "213"),
        ("1324", "123"),
        ("1324", "132"),
        ("1324", "213"),
    ];
    let max_n4 = max_n.min(7);
    for &(pi_pat, s_pat) in &cross {
        let pi_p = parse_sequence(pi_pat);
        let s_p = parse_sequence(s_pat);
        let mut seq = Vec::new();
        for n in 1..=max_n4 {
            let pi_set = avoiding_permutations(n, &[pi_p.clone()]);
            let s_set = avoiding_permutations(n, &[s_p.clone()]);
            seq.push(total_depth(&pi_set, &s_set));
        }
        println!(
            "{:<28} {}",
            format!("Av_{} → Av_{}", pi_pat, s_pat),
            seq.iter()
                .map(|v| format!("{:>12}", v))
                .collect::<Vec<_>>()
                .join("")
        );
    }
}
