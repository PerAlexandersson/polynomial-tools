use combpoly::cayley::for_each_fpf_cayley;
use polynomial_tools::{check_interlacing, format_poly, is_real_rooted};
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};
use combpoly::statistics;
use combpoly::statistics::Stat;

fn main() {
    let max_n: u8 = 11;
    let stat = Stat::Exc;
    let mut polys: Vec<Vec<i64>> = Vec::new();

    for n in 0..=max_n {
        let mut coeffs: Vec<i64> = Vec::new();
        let mut count: u64 = 0;
        for_each_fpf_cayley(n, |w| {
            let v = statistics::compute(w, stat);
            if v >= coeffs.len() {
                coeffs.resize(v + 1, 0);
            }
            coeffs[v] += 1;
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
        let interlace = if polys.is_empty() {
            "".to_string()
        } else {
            match check_interlacing(polys.last().unwrap(), &coeffs) {
                Some(true) => ", interlaces prev".to_string(),
                Some(false) => ", DOES NOT interlace prev".to_string(),
                None => ", interlace check failed".to_string(),
            }
        };
        println!(
            "n={:2}: {:>10} fpf, rr = {}{}, poly = {}",
            n,
            count,
            rr,
            interlace,
            format_poly(&coeffs)
        );
        polys.push(coeffs);
    }

    // Also check des
    println!("\n=== Des stat ===");
    let mut des_polys: Vec<Vec<i64>> = Vec::new();
    for n in 0..=max_n {
        let mut coeffs: Vec<i64> = Vec::new();
        let mut count: u64 = 0;
        for_each_fpf_cayley(n, |w| {
            let v = statistics::compute(w, Stat::Des);
            if v >= coeffs.len() {
                coeffs.resize(v + 1, 0);
            }
            coeffs[v] += 1;
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
        let interlace = if des_polys.is_empty() {
            "".to_string()
        } else {
            match check_interlacing(des_polys.last().unwrap(), &coeffs) {
                Some(true) => ", interlaces prev".to_string(),
                Some(false) => ", DOES NOT interlace prev".to_string(),
                None => ", interlace check failed".to_string(),
            }
        };
        println!(
            "n={:2}: {:>10} fpf, rr = {}{}, poly = {}",
            n,
            count,
            rr,
            interlace,
            format_poly(&coeffs)
        );
        des_polys.push(coeffs);
    }

    println!("\n--- Adaptive recurrence search (exc) ---");
    let search = AdaptiveSearchOptions {
        max_rec_len: 3,
        max_var_deg: 2,
        max_idx_deg: 2,
        max_diff_deg: 1,
        try_denominator: true,
        max_denom_var_deg: 1,
        max_denom_idx_deg: 2,
        min_margin: 2,
        try_inhomogeneous: false,
        verbose: false,
    };
    match find_recurrence_adaptive(&polys, &search) {
        Some(result) => {
            println!(
                "Found recurrence (rec_len={}, var_deg={}, idx_deg={}, diff_deg={}):",
                result.opts.rec_len,
                result.opts.var_deg,
                result.opts.idx_deg,
                result.opts.diff_deg
            );
            println!("  {}", result.recurrence);
        }
        None => println!("No recurrence found."),
    }

    println!("\n--- Adaptive recurrence search (des) ---");
    match find_recurrence_adaptive(&des_polys, &search) {
        Some(result) => {
            println!(
                "Found recurrence (rec_len={}, var_deg={}, idx_deg={}, diff_deg={}):",
                result.opts.rec_len,
                result.opts.var_deg,
                result.opts.idx_deg,
                result.opts.diff_deg
            );
            println!("  {}", result.recurrence);
        }
        None => println!("No recurrence found."),
    }
}
