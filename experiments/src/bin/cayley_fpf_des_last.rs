//! Deeper exploration of FPF Cayley + des refined by last entry.
//! All (n, last) pairs appear real-rooted — check interlacing patterns.

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
    let max_n: u8 = 11;

    // Collect all polys: polys[n][k] = des poly for FPF Cayley of size n with last entry k
    let mut polys: Vec<Vec<Vec<i64>>> = Vec::new();
    let mut overall: Vec<Vec<i64>> = Vec::new();

    for n in 0..=max_n {
        let m = n as usize;
        let max_last = m;
        let mut by_last: Vec<Vec<i64>> = vec![vec![]; max_last + 1];
        let mut total: Vec<i64> = Vec::new();

        for_each_fpf_cayley(n, |w| {
            let v = statistics::compute(w, Stat::Des);
            accumulate_poly(&mut total, v);
            if !w.is_empty() {
                let last = *w.last().unwrap() as usize;
                accumulate_poly(&mut by_last[last], v);
            }
        });

        if total.is_empty() {
            total.push(0);
        }
        polys.push(by_last);
        overall.push(total);
    }

    // 1. Confirm all real-rooted
    println!("=== Real-rootedness check (des, by last entry) ===\n");
    let mut all_rr = true;
    for n in 2..=max_n as usize {
        for k in 1..polys[n].len() {
            if polys[n][k].is_empty() || polys[n][k].iter().all(|&c| c == 0) {
                continue;
            }
            let rr = is_real_rooted(&polys[n][k]);
            if !rr {
                println!("NOT real-rooted: n={}, last={}, poly = {}", n, k, format_poly(&polys[n][k]));
                all_rr = false;
            }
        }
    }
    if all_rr {
        println!("ALL (n, last) pairs are real-rooted for des!\n");
    }

    // 2. Interlacing across n for fixed last entry k
    println!("=== Interlacing across n, fixed last entry (des) ===\n");
    for k in 1..=max_n as usize {
        let mut entries: Vec<(usize, &Vec<i64>)> = Vec::new();
        for n in k..=max_n as usize {
            if k >= polys[n].len() { continue; }
            let p = &polys[n][k];
            if p.is_empty() || p.iter().all(|&c| c == 0) { continue; }
            entries.push((n, p));
        }
        if entries.len() < 2 { continue; }

        print!("last = {}:", k);
        let mut all_interlace = true;
        for i in 1..entries.len() {
            let (n_prev, p_prev) = entries[i - 1];
            let (n_cur, p_cur) = entries[i];
            match check_interlacing(p_prev, p_cur) {
                Some(true) => {}
                Some(false) => {
                    print!(" n={}->{}:NO", n_prev, n_cur);
                    all_interlace = false;
                }
                None => {
                    print!(" n={}->{}:FAIL", n_prev, n_cur);
                    all_interlace = false;
                }
            }
        }
        if all_interlace {
            print!(" all interlace");
        }
        println!();
    }

    // 3. Interlacing across last entry k for fixed n
    println!("\n=== Interlacing across last entry, fixed n (des) ===\n");
    for n in 2..=max_n as usize {
        let mut entries: Vec<(usize, &Vec<i64>)> = Vec::new();
        for k in 1..polys[n].len() {
            let p = &polys[n][k];
            if p.is_empty() || p.iter().all(|&c| c == 0) { continue; }
            entries.push((k, p));
        }
        if entries.len() < 2 { continue; }

        print!("n = {:2}:", n);
        let mut all_interlace = true;
        for i in 1..entries.len() {
            let (k_prev, p_prev) = entries[i - 1];
            let (k_cur, p_cur) = entries[i];
            match check_interlacing(p_prev, p_cur) {
                Some(true) => {}
                Some(false) => {
                    print!(" k={}->{}:NO", k_prev, k_cur);
                    all_interlace = false;
                }
                None => {
                    print!(" k={}->{}:FAIL", k_prev, k_cur);
                    all_interlace = false;
                }
            }
        }
        if all_interlace {
            print!(" all interlace");
        }
        println!();
    }

    // 4. Check: does P_{n,k}(des) interlace with P_{n,k+1}(des) AND P_{n+1,k}(des)?
    // i.e., is this a doubly interlacing array?
    println!("\n=== Double interlacing check ===\n");
    let mut double_ok = true;
    for n in 2..max_n as usize {
        for k in 1..polys[n].len() {
            let p = &polys[n][k];
            if p.is_empty() || p.iter().all(|&c| c == 0) { continue; }

            // Check vs (n, k+1)
            if k + 1 < polys[n].len() {
                let q = &polys[n][k + 1];
                if !q.is_empty() && !q.iter().all(|&c| c == 0) {
                    match check_interlacing(p, q) {
                        Some(false) | None => {
                            println!("  ({},{}) vs ({},{}) = NO/FAIL", n, k, n, k + 1);
                            double_ok = false;
                        }
                        _ => {}
                    }
                }
            }

            // Check vs (n+1, k)
            if n + 1 <= max_n as usize && k < polys[n + 1].len() {
                let q = &polys[n + 1][k];
                if !q.is_empty() && !q.iter().all(|&c| c == 0) {
                    match check_interlacing(p, q) {
                        Some(false) | None => {
                            println!("  ({},{}) vs ({},{}) = NO/FAIL", n, k, n + 1, k);
                            double_ok = false;
                        }
                        _ => {}
                    }
                }
            }
        }
    }
    if double_ok {
        println!("  Doubly interlacing! (both across n and across k)");
    }

    // 5. Look at the polynomials for last = n-1 (second-to-last value)
    println!("\n=== Polynomials for last = n-1 (des) ===\n");
    let mut last_nm1: Vec<Vec<i64>> = Vec::new();
    for n in 2..=max_n as usize {
        let k = n - 1;
        let p = &polys[n][k];
        let count: i64 = p.iter().sum();
        println!("n={:2}, last={}: {:>10} perms, poly = {}", n, k, count, format_poly(p));
        last_nm1.push(p.clone());
    }

    // 6. Recurrence search on each column
    println!("\n=== Recurrence search by last entry (des) ===\n");
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
    for k in 1..=6.min(max_n as usize) {
        let mut seq: Vec<Vec<i64>> = Vec::new();
        for n in 0..=max_n as usize {
            if k >= polys[n].len() || polys[n][k].is_empty() {
                seq.push(vec![0]);
            } else {
                seq.push(polys[n][k].clone());
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

    // 7. Print the triangle of row sums (counts by last entry)
    println!("\n=== Count triangle: FPF Cayley count by (n, last) ===\n");
    for n in 2..=max_n as usize {
        print!("n={:2}: ", n);
        for k in 1..polys[n].len() {
            let count: i64 = polys[n][k].iter().sum();
            print!("{:>10} ", count);
        }
        println!();
    }

    // 8. Check if overall des poly = sum is real-rooted and interlaces
    println!("\n=== Overall des poly interlacing ===\n");
    for i in 1..overall.len() {
        if overall[i].is_empty() || overall[i].iter().all(|&c| c == 0) { continue; }
        if overall[i - 1].is_empty() || overall[i - 1].iter().all(|&c| c == 0) { continue; }
        let il = check_interlacing(&overall[i - 1], &overall[i]);
        println!("n={} -> n={}: interlace = {:?}", i - 1, i, il);
    }
}
