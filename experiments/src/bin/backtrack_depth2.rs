/// Deeper analysis of backtrack search depth.
///
/// 1. EGF coefficients (= average depth)
/// 2. Why are 123=132=231=321? Check pointwise vs only totals.
/// 3. Ratio analysis: avg(n)/avg(n-1), limit behavior.
/// 4. Try to push n=8.
use combpoly::permutation::{
    all_permutations, avoiding_permutations, backtrack_image, parse_sequence,
};
use std::collections::BTreeMap;
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

fn inverse(perm: &[u8]) -> Vec<u8> {
    let n = perm.len();
    let mut inv = vec![0u8; n];
    for (i, &v) in perm.iter().enumerate() {
        inv[(v - 1) as usize] = (i + 1) as u8;
    }
    inv
}

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 { a } else { gcd(b, a % b) }
}

fn main() {
    let max_n: u8 = env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    let patterns = ["123", "132", "213", "231", "312", "321"];

    // ============================================================
    // 1. Check whether the four-pattern equality is pointwise or just totals
    // ============================================================
    println!("=== Pointwise check: is depth(π, Av_123) = depth(π, Av_132) for all π? ===\n");

    for n in 3..=max_n.min(6) {
        let all_pi = all_permutations(n);
        let mut depths: BTreeMap<&str, Vec<u64>> = BTreeMap::new();

        for pat_str in &patterns {
            let pat = parse_sequence(pat_str);
            let s_set = avoiding_permutations(n, &[pat.clone()]);
            let mut dvec = Vec::new();
            for pi in &all_pi {
                let sigma = backtrack_image(pi, &s_set);
                let sigma_pi = compose(&sigma, pi);
                dvec.push(lex_rank(&sigma_pi) + 1);
            }
            depths.insert(pat_str, dvec);
        }

        // Check which patterns have identical depth vectors
        let mut groups: Vec<Vec<&str>> = Vec::new();
        let mut used: Vec<bool> = vec![false; patterns.len()];
        for i in 0..patterns.len() {
            if used[i] { continue; }
            let mut group = vec![patterns[i]];
            used[i] = true;
            for j in i + 1..patterns.len() {
                if !used[j] && depths[patterns[i]] == depths[patterns[j]] {
                    group.push(patterns[j]);
                    used[j] = true;
                }
            }
            groups.push(group);
        }

        println!("n={}: Identical depth vectors:", n);
        for g in &groups {
            println!("  {}", g.join(" = "));
        }

        // For the equal pairs, check if there's a permutation of π that maps one to the other
        // Check: is depth(π, Av_123) = depth(π, Av_132) for all π? (same π)
        let d123 = &depths["123"];
        let d132 = &depths["132"];
        let same_pi = d123 == d132;
        println!("  depth(π, Av_123) = depth(π, Av_132) same π? {}", same_pi);

        // Check: depth(π, Av_123) = depth(π⁻¹, Av_132)?
        let mut inv_match = true;
        for (idx, pi) in all_pi.iter().enumerate() {
            let pi_inv = inverse(pi);
            // Find index of pi_inv in all_pi
            let inv_idx = all_pi.iter().position(|p| p == &pi_inv).unwrap();
            if d123[idx] != d132[inv_idx] {
                inv_match = false;
                break;
            }
        }
        println!("  depth(π, Av_123) = depth(π⁻¹, Av_132)? {}", inv_match);

        // Also check: depth(π, Av_123) = depth(π^r, Av_132)?
        let n_val = n;
        let mut rev_match = true;
        for (idx, pi) in all_pi.iter().enumerate() {
            let pi_rev: Vec<u8> = pi.iter().rev().copied().collect();
            let rev_idx = all_pi.iter().position(|p| p == &pi_rev).unwrap();
            if d123[idx] != d132[rev_idx] {
                rev_match = false;
                break;
            }
        }
        println!("  depth(π, Av_123) = depth(π^rev, Av_132)? {}", rev_match);

        // Check comp
        let mut comp_match = true;
        for (idx, pi) in all_pi.iter().enumerate() {
            let pi_comp: Vec<u8> = pi.iter().map(|&x| n_val + 1 - x).collect();
            let comp_idx = all_pi.iter().position(|p| p == &pi_comp).unwrap();
            if d123[idx] != d132[comp_idx] {
                comp_match = false;
                break;
            }
        }
        println!("  depth(π, Av_123) = depth(π^comp, Av_132)? {}", comp_match);

        // Check: depth(π, Av_123) = depth(π^rev, Av_321)?
        let d321 = &depths["321"];
        let mut rev_321 = true;
        for (idx, _pi) in all_pi.iter().enumerate() {
            let pi_rev: Vec<u8> = all_pi[idx].iter().rev().copied().collect();
            let rev_idx = all_pi.iter().position(|p| p == &pi_rev).unwrap();
            if d123[idx] != d321[rev_idx] {
                rev_321 = false;
                break;
            }
        }
        println!("  depth(π, Av_123) = depth(π^rev, Av_321)? {}", rev_321);

        // Check: depth(π, Av_123) = depth(π^comp, Av_321)?
        let mut comp_321 = true;
        for (idx, pi) in all_pi.iter().enumerate() {
            let pi_comp: Vec<u8> = pi.iter().map(|&x| n_val + 1 - x).collect();
            let comp_idx = all_pi.iter().position(|p| p == &pi_comp).unwrap();
            if d123[idx] != d321[comp_idx] {
                comp_321 = false;
                break;
            }
        }
        println!("  depth(π, Av_123) = depth(π^comp, Av_321)? {}", comp_321);

        // Check depth(π, Av_132) = depth(π^rev, Av_231)?
        let d231 = &depths["231"];
        let mut rev_231 = true;
        for (idx, pi) in all_pi.iter().enumerate() {
            let pi_rev: Vec<u8> = pi.iter().rev().copied().collect();
            let rev_idx = all_pi.iter().position(|p| p == &pi_rev).unwrap();
            if d132[idx] != d231[rev_idx] {
                rev_231 = false;
                break;
            }
        }
        println!("  depth(π, Av_132) = depth(π^rev, Av_231)? {}", rev_231);

        println!();
    }

    // ============================================================
    // 2. EGF coefficients and ratio analysis
    // ============================================================
    println!("\n=== EGF coefficients and ratio analysis (n=1..{}) ===\n", max_n);

    for pat_str in &["213", "123"] {
        let pat = parse_sequence(pat_str);
        let mut avgs: Vec<f64> = Vec::new();
        let mut avg_fracs: Vec<(u64, u64)> = Vec::new();

        for n in 1..=max_n {
            let all_pi = all_permutations(n);
            let s_set = avoiding_permutations(n, &[pat.clone()]);
            let n_fact = (1..=n as u64).product::<u64>();

            let mut total: u64 = 0;
            for pi in &all_pi {
                let sigma = backtrack_image(pi, &s_set);
                let sigma_pi = compose(&sigma, pi);
                total += lex_rank(&sigma_pi) + 1;
            }

            let avg = total as f64 / n_fact as f64;
            let g = gcd(total, n_fact);
            avg_fracs.push((total / g, n_fact / g));
            avgs.push(avg);
        }

        println!("Av_{}: EGF coefficients = average depths:", pat_str);
        for (i, (num, den)) in avg_fracs.iter().enumerate() {
            println!("  n={}: {}/{} ≈ {:.6}", i + 1, num, den, avgs[i]);
        }

        println!("\nRatios avg(n)/avg(n-1):");
        for i in 1..avgs.len() {
            println!(
                "  n={}: {:.6}   (compare n-1={})",
                i + 1,
                avgs[i] / avgs[i - 1],
                i
            );
        }

        println!("\nRatio avg(n) / ((n-1)! / C_{{n-1}}):");
        let catalan: Vec<u64> = (0..=max_n)
            .map(|n| {
                if n == 0 { 1 }
                else {
                    let n = n as u64;
                    (1..=n).fold(1u64, |acc, k| acc * (n + k) / k)  / (n + 1)
                }
            })
            .collect();
        for i in 1..avgs.len() {
            let n = i + 1;
            let n_minus_1_fact = (1..=i as u64).product::<u64>() as f64;
            let c_n_minus_1 = catalan[i] as f64;
            let ratio = avgs[i] / (n_minus_1_fact / c_n_minus_1);
            println!("  n={}: {:.6}   ((n-1)!/C_{{n-1}} = {:.2})", n, ratio, n_minus_1_fact / c_n_minus_1);
        }

        // avg(n) * C_n / n!  (= expected fraction if depth were uniform)
        println!("\navg(n) * C_n / n!:");
        for i in 0..avgs.len() {
            let n = (i + 1) as u64;
            let n_fact = (1..=n).product::<u64>() as f64;
            let c_n = catalan[i + 1] as f64;
            println!("  n={}: {:.6}", n, avgs[i] * c_n / n_fact);
        }

        println!();
    }

    // ============================================================
    // 3. Distribution comparison: sorted depth multisets
    // ============================================================
    println!("\n=== Depth distributions (sorted) for n=5 ===\n");

    let n = 5u8.min(max_n);
    let all_pi = all_permutations(n);
    for pat_str in &patterns {
        let pat = parse_sequence(pat_str);
        let s_set = avoiding_permutations(n, &[pat.clone()]);
        let mut dist: BTreeMap<u64, u64> = BTreeMap::new();
        for pi in &all_pi {
            let sigma = backtrack_image(pi, &s_set);
            let sigma_pi = compose(&sigma, pi);
            let depth = lex_rank(&sigma_pi) + 1;
            *dist.entry(depth).or_default() += 1;
        }
        let vals: Vec<String> = dist.iter().map(|(d, c)| format!("{}:{}", d, c)).collect();
        println!("Av_{}: {}", pat_str, vals.join(" "));
    }
}
