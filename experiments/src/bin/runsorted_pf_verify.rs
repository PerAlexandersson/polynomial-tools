//! Compute n=1..N for all 6 variants, then verify recurrences with exact shape.
//! Uses fused backtracking that generates only run-sorted PFs (no memory blowup).

use combpoly::parking::{for_each_runsorted_pf, RunBreak, RunSort};
use polynomial_tools::format_poly;
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
    let variants: &[(&str, RunBreak, RunSort)] = &[
        ("strict-asc, strict", RunBreak::StrictAsc, RunSort::StrictMin),
        ("strict-asc, weak", RunBreak::StrictAsc, RunSort::WeakMin),
        ("strict-asc, lex", RunBreak::StrictAsc, RunSort::Lex),
        ("non-decr, strict", RunBreak::NonDecr, RunSort::StrictMin),
        ("non-decr, weak", RunBreak::NonDecr, RunSort::WeakMin),
        ("non-decr, lex", RunBreak::NonDecr, RunSort::Lex),
    ];

    let max_n: u8 = 10;

    // Exact search params: rec_len=3, var_deg=1, idx_deg=1, diff_deg=1
    let search = AdaptiveSearchOptions {
        max_rec_len: 3,
        max_var_deg: 1,
        max_idx_deg: 1,
        max_diff_deg: 1,
        try_denominator: true,
        max_denom_var_deg: 0,
        max_denom_idx_deg: 1,
        min_margin: 3,
        try_inhomogeneous: false,
        verbose: false,
    };

    for (vi, (label, rb, rs)) in variants.iter().enumerate() {
        println!("\n=== Variant {}: {} ===\n", vi + 1, label);

        let mut polys: Vec<Vec<i64>> = Vec::new();
        for n in 1..=max_n {
            let mut coeffs: Vec<i64> = Vec::new();
            let mut count = 0u64;
            for_each_runsorted_pf(n, *rb, *rs, &mut |w| {
                accumulate_poly(&mut coeffs, statistics::compute(w, Stat::Des));
                count += 1;
            });
            if coeffs.is_empty() {
                coeffs.push(0);
            }
            println!(
                "n={:2}: {:>12} objects, poly = {}",
                n, count, format_poly(&coeffs)
            );
            polys.push(coeffs);
        }

        println!("\n--- Recurrence search (exact shape) ---");
        match find_recurrence_adaptive(&polys, &search) {
            Some(result) => {
                println!(
                    "FOUND (rec_len={}, var_deg={}, idx_deg={}, diff_deg={}, margin={}):",
                    result.opts.rec_len,
                    result.opts.var_deg,
                    result.opts.idx_deg,
                    result.opts.diff_deg,
                    result.num_equations as i64 - result.num_unknowns as i64
                );
                println!("  {}", result.recurrence);
            }
            None => println!("NO RECURRENCE with this shape — previous fit was SPURIOUS."),
        }
    }
}
