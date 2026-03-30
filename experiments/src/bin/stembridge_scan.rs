/// Scan 4321-avoiding permutations for failure of real-rootedness
/// of the descent polynomial over the weak order lower ideal.
///
/// Uses rayon for parallel scanning. Exits at the first counterexample found.
use combpoly::order;
use combpoly::permutation::avoiding_permutations;
use polynomial_tools::real_rootedness as polynomial;
use combpoly::statistics::Stat;
use rayon::prelude::*;
use std::env;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

fn main() {
    let n: u8 = env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(10);

    let stat = Stat::Des;
    let found = AtomicBool::new(false);
    let progress = AtomicUsize::new(0);

    eprintln!("Generating Av_4321({})...", n);
    let t0 = std::time::Instant::now();
    let perms = avoiding_permutations(n, &[vec![4, 3, 2, 1]]);
    eprintln!(
        "Generated {} permutations in {:?}",
        perms.len(),
        t0.elapsed()
    );

    let total = perms.len();
    eprintln!(
        "Scanning for failure of real-rootedness (weak ideal + des, parallel)..."
    );
    let t0 = std::time::Instant::now();

    let result: Option<(usize, Vec<u8>, Vec<i64>)> = perms
        .par_iter()
        .enumerate()
        .find_map_any(|(i, w)| {
            if found.load(Ordering::Relaxed) {
                return None;
            }

            let done = progress.fetch_add(1, Ordering::Relaxed);
            if done % 50000 == 0 && done > 0 {
                eprintln!(
                    "  {}/{} ({:.1}%) elapsed {:?}",
                    done,
                    total,
                    100.0 * done as f64 / total as f64,
                    t0.elapsed()
                );
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
        Some((i, w, coeffs)) => {
            let ideal_size: i64 = coeffs.iter().sum();
            println!("COUNTEREXAMPLE FOUND at permutation #{}", i);
            println!(
                "w = [{}]",
                w.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            );
            println!(
                "Weak ideal size: {}, polynomial degree: {}",
                ideal_size,
                coeffs.len() - 1
            );
            println!("Coefficients: {:?}", coeffs);
            println!("Polynomial: {}", polynomial::format_poly(&coeffs));
            println!("Log-concave: {}", polynomial::is_log_concave(&coeffs));
            println!("Real-rooted: {}", polynomial::is_real_rooted(&coeffs));
            println!("Time: {:?}", t0.elapsed());
        }
        None => {
            eprintln!(
                "All {} permutations pass real-rootedness. Time: {:?}",
                total,
                t0.elapsed()
            );
        }
    }
}
