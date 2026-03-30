use combpoly::order;
use combpoly::permutation::avoiding_permutations;
use polynomial_tools::real_rootedness as polynomial;
use combpoly::statistics::Stat;
use rayon::prelude::*;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

fn main() {
    let stat = Stat::Des;
    let found = AtomicBool::new(false);
    let progress = AtomicUsize::new(0);

    for n in 14..=17u8 {
        found.store(false, Ordering::Relaxed);
        progress.store(0, Ordering::Relaxed);

        let t0 = std::time::Instant::now();
        let perms = avoiding_permutations(n, &[vec![3, 2, 1]]);
        let total = perms.len();
        eprint!("n={:>2}: {:>8} perms, scanning RR... ", n, total);

        let t1 = std::time::Instant::now();
        let result: Option<(usize, Vec<u8>, Vec<i64>)> = perms
            .par_iter()
            .enumerate()
            .find_map_any(|(i, w)| {
                if found.load(Ordering::Relaxed) { return None; }
                let done = progress.fetch_add(1, Ordering::Relaxed);
                if done % 100000 == 0 && done > 0 {
                    eprint!("{:.0}% ", 100.0 * done as f64 / total as f64);
                }
                let coeffs = order::weak_ideal_poly(w, stat);
                if !polynomial::is_real_rooted(&coeffs) {
                    found.store(true, Ordering::Relaxed);
                    Some((i, w.clone(), coeffs))
                } else {
                    None
                }
            });

        match result {
            Some((_i, w, coeffs)) => {
                let ideal_size: i64 = coeffs.iter().sum();
                eprintln!("FAIL! ({:?})", t1.elapsed());
                println!("  w = {:?}", w);
                println!("  ideal={}, deg={}, coeffs={:?}", ideal_size, coeffs.len()-1, coeffs);
                println!("  log-concave: {}, ULC: {}",
                    polynomial::is_log_concave(&coeffs),
                    polynomial::is_ultra_log_concave(&coeffs));
            }
            None => {
                eprintln!("all pass ({:?})", t1.elapsed());
            }
        }
    }
}
