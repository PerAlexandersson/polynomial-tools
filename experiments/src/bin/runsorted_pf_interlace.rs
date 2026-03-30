//! Check interlacing of consecutive des polynomials for all 6 run-sorted PF variants.

use combpoly::parking::{for_each_runsorted_pf, RunBreak, RunSort};
use polynomial_tools::{check_interlacing, format_poly};
use combpoly::statistics;
use combpoly::statistics::Stat;

fn accumulate_poly(coeffs: &mut Vec<i64>, val: usize) {
    if val >= coeffs.len() {
        coeffs.resize(val + 1, 0);
    }
    coeffs[val] += 1;
}

fn main() {
    let max_n: u8 = 10;
    let stat = Stat::Des;

    let variants: &[(&str, RunBreak, RunSort)] = &[
        ("strict-asc runs, strict mins", RunBreak::StrictAsc, RunSort::StrictMin),
        ("strict-asc runs, weak mins", RunBreak::StrictAsc, RunSort::WeakMin),
        ("strict-asc runs, lex", RunBreak::StrictAsc, RunSort::Lex),
        ("non-decr runs, strict mins", RunBreak::NonDecr, RunSort::StrictMin),
        ("non-decr runs, weak mins", RunBreak::NonDecr, RunSort::WeakMin),
        ("non-decr runs, lex", RunBreak::NonDecr, RunSort::Lex),
    ];

    for (label, rb, rs) in variants {
        println!("=== {} ===\n", label);
        let mut prev_poly: Option<Vec<i64>> = None;

        for n in 1..=max_n {
            let mut coeffs: Vec<i64> = Vec::new();
            for_each_runsorted_pf(n, *rb, *rs, &mut |w| {
                accumulate_poly(&mut coeffs, statistics::compute(w, stat));
            });
            if coeffs.is_empty() {
                coeffs.push(0);
            }

            let il_str = match &prev_poly {
                None => "-".to_string(),
                Some(p) => {
                    if p.iter().all(|&c| c == 0) {
                        "-".to_string()
                    } else {
                        match check_interlacing(p, &coeffs) {
                            Some(true) => "YES".to_string(),
                            Some(false) => "NO".to_string(),
                            None => "FAIL".to_string(),
                        }
                    }
                }
            };

            println!("n={}: il={:4}, poly = {}", n, il_str, format_poly(&coeffs));
            prev_poly = Some(coeffs);
        }
        println!();
    }
}
