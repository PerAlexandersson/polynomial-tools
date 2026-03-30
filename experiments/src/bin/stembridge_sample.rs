/// Sample random 4321-avoiding permutations and check real-rootedness
/// of the descent polynomial over the weak order lower ideal.
///
/// Generates random permutations by rejection sampling from S_n,
/// or by random walks on the 4321-avoiding class.
use combpoly::order;
use combpoly::permutation::contains_pattern;
use polynomial_tools::real_rootedness as polynomial;
use combpoly::statistics::Stat;
use rayon::prelude::*;
use std::env;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

/// Generate a random 4321-avoiding permutation of [1..n]
/// using the "random insertion" method for the dedicated generator structure.
///
/// For 321-avoiding: insert values 1..n one at a time at random valid positions.
/// For 4321-avoiding: build incrementally, at each step choosing a random
/// value from those that don't create a 4321 pattern.
fn random_av4321(n: u8, rng: &mut u64) -> Vec<u8> {
    loop {
        // Build left-to-right, choosing random valid values
        let mut perm: Vec<u8> = Vec::with_capacity(n as usize);
        let mut available: Vec<u8> = (1..=n).collect();
        let mut ok = true;

        for _ in 0..n {
            // Find valid values (those that don't create 4321)
            let pat = [4u8, 3, 2, 1];
            let valid: Vec<usize> = (0..available.len())
                .filter(|&i| {
                    perm.push(available[i]);
                    let good = !has_pattern_ending_here(&perm, &pat);
                    perm.pop();
                    good
                })
                .collect();

            if valid.is_empty() {
                ok = false;
                break;
            }

            let idx = valid[xorshift(rng) as usize % valid.len()];
            perm.push(available[idx]);
            available.remove(idx);
        }

        if ok {
            debug_assert!(!contains_pattern(&perm, &[4, 3, 2, 1]));
            return perm;
        }
    }
}

/// Check if `perm` contains `pattern` using the last element.
fn has_pattern_ending_here(perm: &[u8], pattern: &[u8]) -> bool {
    let k = pattern.len();
    let n = perm.len();
    if k > n || k <= 1 {
        return k <= n && k == 1;
    }
    search_ending(perm, pattern, 0, &mut Vec::with_capacity(k))
}

fn search_ending(perm: &[u8], pattern: &[u8], start: usize, chosen: &mut Vec<u8>) -> bool {
    let k = pattern.len();
    let n = perm.len();
    let last_val = perm[n - 1];

    if chosen.len() == k - 1 {
        for j in 0..chosen.len() {
            let pat_cmp = pattern[j].cmp(&pattern[k - 1]);
            let perm_cmp = chosen[j].cmp(&last_val);
            if pat_cmp != perm_cmp {
                return false;
            }
        }
        return true;
    }

    let remaining = (k - 1) - chosen.len();
    if n < 2 {
        return false;
    }
    let end = (n - 1).saturating_sub(remaining);
    for i in start..=end {
        let pos = chosen.len();
        let mut valid = true;
        for j in 0..pos {
            let pat_cmp = pattern[j].cmp(&pattern[pos]);
            let perm_cmp = chosen[j].cmp(&perm[i]);
            if pat_cmp != perm_cmp {
                valid = false;
                break;
            }
        }
        if valid {
            chosen.push(perm[i]);
            if search_ending(perm, pattern, i + 1, chosen) {
                chosen.pop();
                return true;
            }
            chosen.pop();
        }
    }
    false
}

fn xorshift(state: &mut u64) -> u64 {
    let mut x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *state = x;
    x
}

fn main() {
    let n: u8 = env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(17);

    let num_samples: usize = env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(1000);

    let stat = Stat::Des;
    let found = AtomicBool::new(false);
    let progress = AtomicUsize::new(0);

    eprintln!(
        "Sampling {} random Av_4321({}) permutations...",
        num_samples, n
    );
    let t0 = std::time::Instant::now();

    // Generate samples in parallel (each thread gets its own RNG)
    let samples: Vec<Vec<u8>> = (0..num_samples)
        .into_par_iter()
        .map(|i| {
            let mut rng = 123456789u64 ^ (i as u64 * 6364136223846793005 + 1);
            random_av4321(n, &mut rng)
        })
        .collect();

    eprintln!("Generated {} samples in {:?}", samples.len(), t0.elapsed());
    eprintln!("Scanning for real-rootedness failures...");
    let t0 = std::time::Instant::now();

    let result: Option<(usize, Vec<u8>, Vec<i64>)> = samples
        .par_iter()
        .enumerate()
        .find_map_any(|(i, w)| {
            if found.load(Ordering::Relaxed) {
                return None;
            }

            let done = progress.fetch_add(1, Ordering::Relaxed);
            if done % 100 == 0 && done > 0 {
                eprintln!("  {}/{} elapsed {:?}", done, num_samples, t0.elapsed());
            }

            let coeffs = order::weak_ideal_poly(w, stat);

            if !polynomial::is_real_rooted(&coeffs) {
                found.store(true, Ordering::Relaxed);
                Some((i, w.clone(), coeffs))
            } else {
                None
            }
        });

    eprintln!();
    match result {
        Some((_i, w, coeffs)) => {
            let ideal_size: i64 = coeffs.iter().sum();
            println!("COUNTEREXAMPLE FOUND!");
            println!(
                "w = [{}]",
                w.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            );
            println!(
                "Weak ideal size: {}, deg: {}",
                ideal_size,
                coeffs.len() - 1
            );
            println!("Coefficients: {:?}", coeffs);
            println!("Polynomial: {}", polynomial::format_poly(&coeffs));
            println!("Log-concave: {}", polynomial::is_log_concave(&coeffs));
            println!("Real-rooted: {}", polynomial::is_real_rooted(&coeffs));
        }
        None => {
            eprintln!(
                "All {} samples pass real-rootedness. Time: {:?}",
                num_samples,
                t0.elapsed()
            );
        }
    }
    eprintln!("Total time: {:?}", t0.elapsed());
}
