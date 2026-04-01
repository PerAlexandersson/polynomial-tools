/// Cross-length search depth: average depth when pattern lengths differ.
///
/// E.g., π ∈ Av_1234 and S = Av_123, or π ∈ Av_123 and S = Av_1234.
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

fn avg_frac(total: u64, count: u64) -> String {
    if count == 0 {
        return "N/A".into();
    }
    let g = gcd(total, count);
    let (n, d) = (total / g, count / g);
    if d == 1 {
        format!("{}", n)
    } else {
        format!("{}/{}", n, d)
    }
}

fn compute_depths(pi_set: &[Vec<u8>], s_set: &[Vec<u8>]) -> (u64, u64, u64) {
    let mut total: u64 = 0;
    let mut max_d: u64 = 0;
    for pi in pi_set {
        let sigma = backtrack_image(pi, s_set);
        let sigma_pi = compose(&sigma, pi);
        let d = lex_rank(&sigma_pi) + 1;
        total += d;
        max_d = max_d.max(d);
    }
    (total, pi_set.len() as u64, max_d)
}

fn main() {
    let max_n: u8 = env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(7);

    let pats3 = ["123", "132", "213"];
    let pats4 = ["1234", "1243", "1324", "1342", "2134", "2143"];

    // ============================================================
    // π ∈ Av_τ₃, S = Av_τ₄ (length-3 search order, length-4 target)
    // ============================================================
    println!("=== π ∈ Av_τ₃(n), S = Av_τ₄(n) ===\n");

    for pi_pat in &pats3 {
        let pi_p = parse_sequence(pi_pat);
        for s_pat in &pats4 {
            let s_p = parse_sequence(s_pat);
            print!("π∈Av_{}, S=Av_{}: ", pi_pat, s_pat);
            for n in 1..=max_n {
                let pi_set = avoiding_permutations(n, &[pi_p.clone()]);
                let s_set = avoiding_permutations(n, &[s_p.clone()]);
                let (total, count, _) = compute_depths(&pi_set, &s_set);
                print!("{} ", avg_frac(total, count));
            }
            println!();
        }
        println!();
    }

    // ============================================================
    // π ∈ Av_τ₄, S = Av_τ₃ (length-4 search order, length-3 target)
    // ============================================================
    println!("\n=== π ∈ Av_τ₄(n), S = Av_τ₃(n) ===\n");

    for pi_pat in &pats4 {
        let pi_p = parse_sequence(pi_pat);
        for s_pat in &pats3 {
            let s_p = parse_sequence(s_pat);
            print!("π∈Av_{}, S=Av_{}: ", pi_pat, s_pat);
            let mut totals = Vec::new();
            for n in 1..=max_n {
                let pi_set = avoiding_permutations(n, &[pi_p.clone()]);
                let s_set = avoiding_permutations(n, &[s_p.clone()]);
                let (total, count, _) = compute_depths(&pi_set, &s_set);
                totals.push(total);
                print!("{} ", avg_frac(total, count));
            }
            println!("  totals={:?}", totals);
        }
        println!();
    }

    // ============================================================
    // π ∈ S_n, S = Av_τ₄ (all search orders, length-4 target)
    // ============================================================
    println!("\n=== π ∈ S_n, S = Av_τ₄(n) ===\n");

    for s_pat in &pats4 {
        let s_p = parse_sequence(s_pat);
        print!("S=Av_{}: ", s_pat);
        let mut totals = Vec::new();
        for n in 1..=max_n {
            let pi_set = all_permutations(n);
            let s_set = avoiding_permutations(n, &[s_p.clone()]);
            let (total, count, max_d) = compute_depths(&pi_set, &s_set);
            totals.push(total);
            print!("{}(max={}) ", avg_frac(total, count), max_d);
        }
        println!("\n  totals={:?}", totals);
    }

    // ============================================================
    // Same length-4 pattern for both π and S
    // ============================================================
    println!("\n=== π ∈ Av_τ₄, S = Av_τ₄ (same pattern) ===\n");

    for pat in &pats4 {
        let p = parse_sequence(pat);
        print!("Av_{} → Av_{}: ", pat, pat);
        for n in 1..=max_n {
            let s_set = avoiding_permutations(n, &[p.clone()]);
            let (total, count, _) = compute_depths(&s_set, &s_set);
            print!("{} ", avg_frac(total, count));
        }
        println!();
    }

    // ============================================================
    // Cross length-4 pairs (a few interesting ones)
    // ============================================================
    println!("\n=== π ∈ Av_τ₄, S = Av_σ₄ (cross, select pairs) ===\n");

    let cross4 = [
        ("1234", "1243"),
        ("1234", "2134"),
        ("1243", "1234"),
        ("1324", "1234"),
        ("1234", "1324"),
    ];
    for &(pi_pat, s_pat) in &cross4 {
        let pi_p = parse_sequence(pi_pat);
        let s_p = parse_sequence(s_pat);
        print!("π∈Av_{}, S=Av_{}: ", pi_pat, s_pat);
        let mut totals = Vec::new();
        for n in 1..=max_n {
            let pi_set = avoiding_permutations(n, &[pi_p.clone()]);
            let s_set = avoiding_permutations(n, &[s_p.clone()]);
            let (total, count, _) = compute_depths(&pi_set, &s_set);
            totals.push(total);
            print!("{} ", avg_frac(total, count));
        }
        println!("  totals={:?}", totals);
    }
}
