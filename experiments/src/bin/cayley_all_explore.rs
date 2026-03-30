//! Explore ALL Cayley permutations (not just FPF) with des/exc,
//! refined by last entry and by maximum value.

use combpoly::cayley::all_cayley_permutations;
use polynomial_tools::{check_interlacing, format_poly, is_real_rooted};
use combpoly::statistics;
use combpoly::statistics::Stat;

fn accumulate_poly(coeffs: &mut Vec<i64>, val: usize) {
    if val >= coeffs.len() {
        coeffs.resize(val + 1, 0);
    }
    coeffs[val] += 1;
}

fn main() {
    let max_n: u8 = 9; // Fubini(9) = 7M, manageable

    // === 1. All Cayley perms, overall exc and des ===
    println!("=== All Cayley perms, exc ===\n");
    let mut exc_polys: Vec<Vec<i64>> = Vec::new();
    let mut des_polys: Vec<Vec<i64>> = Vec::new();

    for n in 0..=max_n {
        let all = all_cayley_permutations(n);
        let mut exc_coeffs: Vec<i64> = Vec::new();
        let mut des_coeffs: Vec<i64> = Vec::new();
        for w in &all {
            accumulate_poly(&mut exc_coeffs, statistics::compute(w, Stat::Exc));
            accumulate_poly(&mut des_coeffs, statistics::compute(w, Stat::Des));
        }
        if exc_coeffs.is_empty() { exc_coeffs.push(0); }
        if des_coeffs.is_empty() { des_coeffs.push(0); }

        let exc_rr = is_real_rooted(&exc_coeffs);
        let des_rr = is_real_rooted(&des_coeffs);
        let exc_il: String = if exc_polys.is_empty() { "".into() } else {
            match check_interlacing(exc_polys.last().unwrap(), &exc_coeffs) {
                Some(true) => ", il=Y".into(),
                Some(false) => ", il=N".into(),
                None => ", il=?".into(),
            }
        };
        let des_il: String = if des_polys.is_empty() { "".into() } else {
            match check_interlacing(des_polys.last().unwrap(), &des_coeffs) {
                Some(true) => ", il=Y".into(),
                Some(false) => ", il=N".into(),
                None => ", il=?".into(),
            }
        };
        println!("n={}: exc rr={}{}, poly = {}", n, exc_rr, exc_il, format_poly(&exc_coeffs));
        println!("      des rr={}{}, poly = {}", des_rr, des_il, format_poly(&des_coeffs));

        exc_polys.push(exc_coeffs);
        des_polys.push(des_coeffs);
    }

    // === 2. All Cayley perms refined by max value m ===
    println!("\n=== All Cayley perms, des, by max value ===\n");
    for n in 1..=max_n {
        let all = all_cayley_permutations(n);
        let m = n as usize;
        let mut by_max: Vec<Vec<i64>> = vec![vec![]; m + 1];

        for w in &all {
            let max_val = *w.iter().max().unwrap() as usize;
            accumulate_poly(&mut by_max[max_val], statistics::compute(w, Stat::Des));
        }

        println!("n = {}:", n);
        for mv in 1..=m {
            if by_max[mv].is_empty() || by_max[mv].iter().all(|&c| c == 0) { continue; }
            let count: i64 = by_max[mv].iter().sum();
            let rr = is_real_rooted(&by_max[mv]);
            println!("  max={}: {:>8} perms, rr={:5}, poly = {}", mv, count, rr, format_poly(&by_max[mv]));
        }
    }

    // === 3. All Cayley perms refined by max value, exc stat ===
    println!("\n=== All Cayley perms, exc, by max value ===\n");
    for n in 1..=max_n {
        let all = all_cayley_permutations(n);
        let m = n as usize;
        let mut by_max: Vec<Vec<i64>> = vec![vec![]; m + 1];

        for w in &all {
            let max_val = *w.iter().max().unwrap() as usize;
            accumulate_poly(&mut by_max[max_val], statistics::compute(w, Stat::Exc));
        }

        println!("n = {}:", n);
        for mv in 1..=m {
            if by_max[mv].is_empty() || by_max[mv].iter().all(|&c| c == 0) { continue; }
            let count: i64 = by_max[mv].iter().sum();
            let rr = is_real_rooted(&by_max[mv]);
            println!("  max={}: {:>8} perms, rr={:5}, poly = {}", mv, count, rr, format_poly(&by_max[mv]));
        }
    }

    // === 4. All Cayley, des by last entry ===
    println!("\n=== All Cayley perms, des, by last entry ===\n");
    for n in 1..=max_n {
        let all = all_cayley_permutations(n);
        let m = n as usize;
        let mut by_last: Vec<Vec<i64>> = vec![vec![]; m + 1];

        for w in &all {
            let last = *w.last().unwrap() as usize;
            accumulate_poly(&mut by_last[last], statistics::compute(w, Stat::Des));
        }

        println!("n = {}:", n);
        let mut all_rr = true;
        for k in 1..=m {
            if by_last[k].is_empty() || by_last[k].iter().all(|&c| c == 0) { continue; }
            let count: i64 = by_last[k].iter().sum();
            let rr = is_real_rooted(&by_last[k]);
            if !rr { all_rr = false; }
            println!("  last={}: {:>8} perms, rr={:5}, poly = {}", k, count, rr, format_poly(&by_last[k]));
        }
        if all_rr {
            println!("  (all real-rooted)");
        }
    }

    // === 5. FPF Cayley, des by max value ===
    println!("\n=== FPF Cayley perms, des, by max value ===\n");
    for n in 2..=max_n {
        let all = all_cayley_permutations(n);
        let m = n as usize;
        let mut by_max: Vec<Vec<i64>> = vec![vec![]; m + 1];

        for w in &all {
            // Check FPF
            if w.iter().enumerate().any(|(i, &v)| v as usize == i + 1) { continue; }
            let max_val = *w.iter().max().unwrap() as usize;
            accumulate_poly(&mut by_max[max_val], statistics::compute(w, Stat::Des));
        }

        println!("n = {}:", n);
        let mut all_rr = true;
        for mv in 1..=m {
            if by_max[mv].is_empty() || by_max[mv].iter().all(|&c| c == 0) { continue; }
            let count: i64 = by_max[mv].iter().sum();
            let rr = is_real_rooted(&by_max[mv]);
            if !rr { all_rr = false; }
            println!("  max={}: {:>8} perms, rr={:5}, poly = {}", mv, count, rr, format_poly(&by_max[mv]));
        }
        if all_rr {
            println!("  (all real-rooted)");
        }
    }
}
