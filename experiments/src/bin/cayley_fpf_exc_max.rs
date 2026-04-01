//! FPF Cayley perms, exc stat, refined by max value.
//! Check palindromicity and real-rootedness.

use combpoly::cayley::for_each_fpf_cayley;
use combpoly::statistics;
use combpoly::statistics::Stat;
use polynomial_tools::{format_poly, is_real_rooted};

fn accumulate_poly(coeffs: &mut Vec<i64>, val: usize) {
    if val >= coeffs.len() {
        coeffs.resize(val + 1, 0);
    }
    coeffs[val] += 1;
}

fn is_palindromic(p: &[i64]) -> bool {
    let n = p.len();
    (0..n / 2).all(|i| p[i] == p[n - 1 - i])
}

fn main() {
    let max_n: u8 = 10;

    println!("=== FPF Cayley, exc, by max value ===\n");
    for n in 2..=max_n {
        let m = n as usize;
        let mut by_max: Vec<Vec<i64>> = vec![vec![]; m + 1];

        for_each_fpf_cayley(n, |w| {
            let max_val = *w.iter().max().unwrap() as usize;
            accumulate_poly(&mut by_max[max_val], statistics::compute(w, Stat::Exc));
        });

        println!("n = {}:", n);
        for mv in 1..=m {
            if by_max[mv].is_empty() || by_max[mv].iter().all(|&c| c == 0) {
                continue;
            }
            let count: i64 = by_max[mv].iter().sum();
            let rr = is_real_rooted(&by_max[mv]);
            let pal = is_palindromic(&by_max[mv]);
            println!(
                "  max={}: {:>8} perms, rr={:5}, pal={:5}, poly = {}",
                mv,
                count,
                rr,
                pal,
                format_poly(&by_max[mv])
            );
        }
    }

    // Also check: all Cayley, exc by (max value, last entry) — is this doubly refined still rr?
    println!("\n=== All Cayley, exc, by (max, last) — real-rootedness ===\n");
    use combpoly::cayley::all_cayley_permutations;
    for n in 1..=8u8 {
        let all = all_cayley_permutations(n);
        let m = n as usize;
        let mut by_ml: Vec<Vec<Vec<i64>>> = vec![vec![vec![]; m + 1]; m + 1]; // [max][last]

        for w in &all {
            let max_val = *w.iter().max().unwrap() as usize;
            let last = *w.last().unwrap() as usize;
            accumulate_poly(&mut by_ml[max_val][last], statistics::compute(w, Stat::Exc));
        }

        let mut all_rr = true;
        let mut failures = Vec::new();
        for mv in 1..=m {
            for k in 1..=m {
                let p = &by_ml[mv][k];
                if p.is_empty() || p.iter().all(|&c| c == 0) {
                    continue;
                }
                if !is_real_rooted(p) {
                    all_rr = false;
                    let count: i64 = p.iter().sum();
                    failures.push(format!(
                        "  max={}, last={}: {} perms, poly = {}",
                        mv,
                        k,
                        count,
                        format_poly(p)
                    ));
                }
            }
        }
        if all_rr {
            println!("n={}: all (max, last) pairs real-rooted", n);
        } else {
            println!("n={}: FAILURES:", n);
            for f in &failures {
                println!("{}", f);
            }
        }
    }
}
