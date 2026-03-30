//! Recurrence search for run-sorted PF des polynomials, one variant at a time.
//! Pass variant number 1-6 as command line argument.

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
    let args: Vec<String> = std::env::args().collect();
    let variant: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(1);

    let (label, rb, rs): (&str, RunBreak, RunSort) = match variant {
        1 => ("strict-asc runs, strict mins", RunBreak::StrictAsc, RunSort::StrictMin),
        2 => ("strict-asc runs, weak mins", RunBreak::StrictAsc, RunSort::WeakMin),
        3 => ("strict-asc runs, lex", RunBreak::StrictAsc, RunSort::Lex),
        4 => ("non-decr runs, strict mins", RunBreak::NonDecr, RunSort::StrictMin),
        5 => ("non-decr runs, weak mins", RunBreak::NonDecr, RunSort::WeakMin),
        6 => ("non-decr runs, lex", RunBreak::NonDecr, RunSort::Lex),
        _ => panic!("variant must be 1-6"),
    };

    println!("=== Variant {}: {} ===\n", variant, label);

    let max_n: u8 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(10);
    let stat = match args.get(3).map(|s| s.as_str()) {
        Some("exc") => Stat::Exc,
        Some("peak") => Stat::Peak,
        Some("inv") => Stat::Inv,
        Some("maj") => Stat::Maj,
        _ => Stat::Des,
    };
    println!("stat = {:?}, max_n = {}\n", stat, max_n);

    let mut polys: Vec<Vec<i64>> = Vec::new();

    for n in 1..=max_n {
        let mut coeffs: Vec<i64> = Vec::new();
        let mut count = 0u64;
        for_each_runsorted_pf(n, rb, rs, &mut |w| {
            accumulate_poly(&mut coeffs, statistics::compute(w, stat));
            count += 1;
        });
        if coeffs.is_empty() {
            coeffs.push(0);
        }
        println!("n={}: {:>10} objects, poly = {}", n, count, format_poly(&coeffs));
        polys.push(coeffs);
    }

    println!("\n--- Recurrence search ---");
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
        verbose: true,
    };
    match find_recurrence_adaptive(&polys, &search) {
        Some(result) => {
            println!(
                "FOUND (rec_len={}, var_deg={}, idx_deg={}, diff_deg={}):",
                result.opts.rec_len, result.opts.var_deg, result.opts.idx_deg, result.opts.diff_deg
            );
            println!("  {}", result.recurrence);
        }
        None => println!("No recurrence found."),
    }
}
