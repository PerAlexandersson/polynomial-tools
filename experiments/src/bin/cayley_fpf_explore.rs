//! Explore FPF Cayley permutations refined by last entry.
//!
//! For each n and each possible last value k, compute the exc polynomial
//! over FPF Cayley perms of size n with last entry = k.
//! Check real-rootedness, interlacing, and look for recurrences.

use combpoly::cayley::for_each_fpf_cayley;
use polynomial_tools::{check_interlacing, format_poly, is_real_rooted};
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};
use combpoly::statistics;
use combpoly::statistics::Stat;

fn accumulate_poly(coeffs: &mut Vec<i64>, val: usize) {
    if val >= coeffs.len() {
        coeffs.resize(val + 1, 0);
    }
    coeffs[val] += 1;
}

fn main() {
    let max_n: u8 = 10; // keep manageable
    let stat = Stat::Exc;

    println!("=== FPF Cayley permutations by last entry, stat = exc ===\n");

    // For each n, collect polynomials by last entry value
    // polys_by_last[n][k] = polynomial for FPF Cayley perms of size n with last entry k
    let mut polys_by_last: Vec<Vec<Vec<i64>>> = Vec::new();
    // Also collect the overall polynomial for reference
    let mut overall_polys: Vec<Vec<i64>> = Vec::new();

    for n in 0..=max_n {
        let m = n as usize;
        // last entry can be 1..=n (but not n for FPF since that's a fixed point when last pos = n)
        let max_last = if m == 0 { 0 } else { m };
        let mut by_last: Vec<Vec<i64>> = vec![vec![]; max_last + 1]; // index by last value (0 unused)
        let mut overall: Vec<i64> = Vec::new();

        for_each_fpf_cayley(n, |w| {
            let v = statistics::compute(w, stat);
            accumulate_poly(&mut overall, v);
            if !w.is_empty() {
                let last = *w.last().unwrap() as usize;
                accumulate_poly(&mut by_last[last], v);
            }
        });

        if overall.is_empty() {
            overall.push(0);
        }

        println!("n = {}:  total = {}", n, overall.iter().sum::<i64>());

        for k in 1..=max_last {
            if by_last[k].is_empty() {
                continue;
            }
            if by_last[k].iter().all(|&c| c == 0) {
                continue;
            }
            let count: i64 = by_last[k].iter().sum();
            let rr = is_real_rooted(&by_last[k]);
            let is_fpf_last = k != m; // last entry != last position (1-indexed)
            println!(
                "  last={}: {:>8} perms, rr = {:5}, fpf_pos = {:5}, poly = {}",
                k,
                count,
                rr,
                is_fpf_last,
                format_poly(&by_last[k])
            );
        }

        polys_by_last.push(by_last);
        overall_polys.push(overall);
        println!();
    }

    // Check interlacing across n for each fixed last entry value
    println!("=== Interlacing across n for fixed last entry ===\n");
    let max_last_check = max_n as usize;
    for k in 1..=max_last_check {
        println!("last = {}:", k);
        let mut prev: Option<&Vec<i64>> = None;
        for n in 0..=max_n as usize {
            if n < k || k >= polys_by_last[n].len() {
                continue;
            }
            let poly = &polys_by_last[n][k];
            if poly.is_empty() || poly.iter().all(|&c| c == 0) {
                continue;
            }
            let interlace_str = match prev {
                Some(p) => match check_interlacing(p, poly) {
                    Some(true) => "YES".to_string(),
                    Some(false) => "NO".to_string(),
                    None => "FAIL".to_string(),
                },
                None => "-".to_string(),
            };
            println!(
                "  n={:2}: interlace_prev={:4}, poly = {}",
                n,
                interlace_str,
                format_poly(poly)
            );
            prev = Some(poly);
        }
        println!();
    }

    // Check interlacing across last values for fixed n
    println!("=== Interlacing across last entry for fixed n ===\n");
    for n in 2..=max_n as usize {
        println!("n = {}:", n);
        let mut prev: Option<&Vec<i64>> = None;
        let mut prev_k = 0;
        for k in 1..=n {
            let poly = &polys_by_last[n][k];
            if poly.is_empty() || poly.iter().all(|&c| c == 0) {
                continue;
            }
            let interlace_str = match prev {
                Some(p) => {
                    let label = format!("vs last={}", prev_k);
                    match check_interlacing(p, poly) {
                        Some(true) => format!("{}: YES", label),
                        Some(false) => format!("{}: NO", label),
                        None => format!("{}: FAIL", label),
                    }
                }
                None => "-".to_string(),
            };
            println!(
                "  last={}: {:4}, poly = {}",
                k,
                interlace_str,
                format_poly(poly)
            );
            prev = Some(poly);
            prev_k = k;
        }
        println!();
    }

    // Try recurrences on the "last = k" columns
    println!("=== Recurrence search for fixed last entry ===\n");
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

    for k in 1..=max_last_check {
        // Collect polys for last=k across n values
        let mut seq: Vec<Vec<i64>> = Vec::new();
        for n in 0..=max_n as usize {
            if k >= polys_by_last[n].len() {
                seq.push(vec![0]);
            } else {
                let p = &polys_by_last[n][k];
                if p.is_empty() {
                    seq.push(vec![0]);
                } else {
                    seq.push(p.clone());
                }
            }
        }
        print!("last = {}: ", k);
        match find_recurrence_adaptive(&seq, &search) {
            Some(result) => {
                println!(
                    "FOUND (rec_len={}, var_deg={}, idx_deg={}, diff_deg={})",
                    result.opts.rec_len,
                    result.opts.var_deg,
                    result.opts.idx_deg,
                    result.opts.diff_deg
                );
                println!("  {}", result.recurrence);
            }
            None => println!("none found"),
        }
    }

    // Also try: refine by first entry
    println!("\n=== FPF Cayley by FIRST entry, stat = exc ===\n");
    let mut polys_by_first: Vec<Vec<Vec<i64>>> = Vec::new();
    for n in 0..=max_n {
        let m = n as usize;
        let max_first = if m == 0 { 0 } else { m };
        let mut by_first: Vec<Vec<i64>> = vec![vec![]; max_first + 1];

        for_each_fpf_cayley(n, |w| {
            let v = statistics::compute(w, stat);
            if !w.is_empty() {
                let first = w[0] as usize;
                accumulate_poly(&mut by_first[first], v);
            }
        });

        println!("n = {}:", n);
        for k in 1..=max_first {
            if by_first[k].is_empty() || by_first[k].iter().all(|&c| c == 0) {
                continue;
            }
            let count: i64 = by_first[k].iter().sum();
            let rr = is_real_rooted(&by_first[k]);
            let is_fpf_first = k != 1; // first entry != 1 (position 1)
            println!(
                "  first={}: {:>8} perms, rr = {:5}, fpf_pos = {:5}, poly = {}",
                k,
                count,
                rr,
                is_fpf_first,
                format_poly(&by_first[k])
            );
        }
        polys_by_first.push(by_first);
        println!();
    }

    // Check: does sum over last entry k != n give back the overall poly? (sanity check)
    println!("=== Sanity: sum over last entries vs overall ===");
    for n in 2..=max_n as usize {
        let mut sum_poly: Vec<i64> = Vec::new();
        for k in 1..polys_by_last[n].len() {
            for (i, &c) in polys_by_last[n][k].iter().enumerate() {
                if i >= sum_poly.len() {
                    sum_poly.resize(i + 1, 0);
                }
                sum_poly[i] += c;
            }
        }
        let ok = sum_poly == overall_polys[n];
        println!("n={}: sum matches overall = {}", n, ok);
    }

    // Try des stat too, briefly
    println!("\n=== FPF Cayley by last entry, stat = des ===\n");
    for n in 0..=max_n {
        let m = n as usize;
        let max_last = if m == 0 { 0 } else { m };
        let mut by_last: Vec<Vec<i64>> = vec![vec![]; max_last + 1];

        for_each_fpf_cayley(n, |w| {
            let v = statistics::compute(w, Stat::Des);
            if !w.is_empty() {
                let last = *w.last().unwrap() as usize;
                accumulate_poly(&mut by_last[last], v);
            }
        });

        if m < 2 { continue; }

        println!("n = {}:", n);
        for k in 1..=max_last {
            if by_last[k].is_empty() || by_last[k].iter().all(|&c| c == 0) {
                continue;
            }
            let count: i64 = by_last[k].iter().sum();
            let rr = is_real_rooted(&by_last[k]);
            println!(
                "  last={}: {:>8} perms, rr = {:5}, poly = {}",
                k,
                count,
                rr,
                format_poly(&by_last[k])
            );
        }
        println!();
    }
}
