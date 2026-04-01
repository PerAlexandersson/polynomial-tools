//! Produce the three tables requested:
//! 1. FPF Cayley + exc, n=1..12
//! 2. All Cayley by max value + exc, n=9
//! 3. FPF Cayley by last entry + exc, n=9

use combpoly::cayley::{all_cayley_permutations, for_each_fpf_cayley};
use combpoly::statistics;
use combpoly::statistics::Stat;
use polynomial_tools::{format_poly, is_real_rooted};

fn accumulate_poly(coeffs: &mut Vec<i64>, val: usize) {
    if val >= coeffs.len() {
        coeffs.resize(val + 1, 0);
    }
    coeffs[val] += 1;
}

fn main() {
    let stat = Stat::Exc;

    // === Table 1: FPF Cayley + exc, n=1..12 ===
    println!("=== Table 1: FPF Cayley + exc, n=1..12 ===\n");
    for n in 1..=12u8 {
        let mut coeffs: Vec<i64> = Vec::new();
        let mut count: u64 = 0;
        for_each_fpf_cayley(n, |w| {
            accumulate_poly(&mut coeffs, statistics::compute(w, stat));
            count += 1;
        });
        if coeffs.is_empty() {
            coeffs.push(0);
        }
        let rr = if coeffs.iter().all(|&c| c == 0) {
            true
        } else {
            is_real_rooted(&coeffs)
        };
        println!(
            "n={:2}: {:>12} fpf, rr={}, poly = {}",
            n,
            count,
            rr,
            format_poly(&coeffs)
        );
    }

    // === Table 2: All Cayley, exc by max value, n=9 ===
    println!("\n=== Table 2: All Cayley, exc by max value, n=9 ===\n");
    {
        let n = 9u8;
        let all = all_cayley_permutations(n);
        let m = n as usize;
        let mut by_max: Vec<Vec<i64>> = vec![vec![]; m + 1];
        for w in &all {
            let max_val = *w.iter().max().unwrap() as usize;
            accumulate_poly(&mut by_max[max_val], statistics::compute(w, stat));
        }
        for mv in 1..=m {
            if by_max[mv].is_empty() || by_max[mv].iter().all(|&c| c == 0) {
                continue;
            }
            let count: i64 = by_max[mv].iter().sum();
            let rr = is_real_rooted(&by_max[mv]);
            println!(
                "m={}: {:>8} perms, rr={}, poly = {}",
                mv,
                count,
                rr,
                format_poly(&by_max[mv])
            );
        }
    }

    // === Table 3: FPF Cayley, exc by last entry, n=9 ===
    println!("\n=== Table 3: FPF Cayley, exc by last entry, n=9 ===\n");
    {
        let n = 9u8;
        let m = n as usize;
        let mut by_last: Vec<Vec<i64>> = vec![vec![]; m + 1];
        for_each_fpf_cayley(n, |w| {
            if !w.is_empty() {
                let last = *w.last().unwrap() as usize;
                accumulate_poly(&mut by_last[last], statistics::compute(w, stat));
            }
        });
        for k in 1..=m {
            if by_last[k].is_empty() || by_last[k].iter().all(|&c| c == 0) {
                continue;
            }
            let count: i64 = by_last[k].iter().sum();
            let rr = is_real_rooted(&by_last[k]);
            println!(
                "last={}: {:>8} perms, rr={}, poly = {}",
                k,
                count,
                rr,
                format_poly(&by_last[k])
            );
        }
    }
}
