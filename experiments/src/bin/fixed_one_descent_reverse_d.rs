//! Explore the fixed-one-descent family for fixed m and varying d, using
//! reversed polynomials so the degree grows with d.
//!
//! For
//!   L_m^(d)(t) = sum_r binom(d,r) binom(m,r) t^r - t,
//! define
//!   R_d^(m)(x) = x^d L_m^(d)(1/x).
//!
//! Likewise for the unperturbed family
//!   B_m^(d)(t) = L_m^(d)(t) + t,
//! define
//!   S_d^(m)(x) = x^d B_m^(d)(1/x).
//!
//! We test:
//! - real-rootedness,
//! - consecutive interlacing in d,
//! - adaptive recurrence search in d.

use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly, is_real_rooted};
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};

fn trim(mut p: Vec<i64>) -> Vec<i64> {
    while p.len() > 1 && p.last() == Some(&0) {
        p.pop();
    }
    if p.is_empty() {
        vec![0]
    } else {
        p
    }
}

fn poly_degree(p: &[i64]) -> usize {
    let mut d = p.len();
    while d > 1 && p[d - 1] == 0 {
        d -= 1;
    }
    d - 1
}

fn poly_mul_t(p: &[i64]) -> Vec<i64> {
    let mut out = vec![0; p.len() + 1];
    for (i, &c) in p.iter().enumerate() {
        out[i + 1] = c;
    }
    trim(out)
}

fn interlaces_weak(f: &[i64], g: &[i64]) -> bool {
    let f = trim(f.to_vec());
    let g = trim(g.to_vec());
    if f == [0] {
        return g == [0] || check_weak_interlacing(&[], &g) == Some(true);
    }
    if g == [0] {
        return true;
    }

    let df = poly_degree(&f);
    let dg = poly_degree(&g);
    if dg == df + 1 {
        check_weak_interlacing(&f, &g) == Some(true)
    } else if dg == df {
        check_weak_interlacing(&g, &poly_mul_t(&f)) == Some(true)
    } else {
        false
    }
}

fn binomial_i64(n: usize, k: usize) -> i64 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut num: i128 = 1;
    let mut den: i128 = 1;
    for i in 0..k {
        num *= (n - i) as i128;
        den *= (i + 1) as i128;
    }
    i64::try_from(num / den).expect("binomial overflow")
}

fn one_descent_unperturbed_polynomial(d: usize, m: usize) -> Vec<i64> {
    let degree = d.min(m);
    let mut coeffs = vec![0; degree + 1];
    for (r, coeff) in coeffs.iter_mut().enumerate().take(degree + 1) {
        *coeff = binomial_i64(d, r) * binomial_i64(m, r);
    }
    trim(coeffs)
}

fn one_descent_polynomial(d: usize, m: usize) -> Vec<i64> {
    let mut coeffs = one_descent_unperturbed_polynomial(d, m);
    if coeffs.len() < 2 {
        coeffs.resize(2, 0);
    }
    coeffs[1] -= 1;
    trim(coeffs)
}

fn reverse_to_degree(poly: &[i64], target_degree: usize) -> Vec<i64> {
    let mut out = vec![0; target_degree + 1];
    for (r, &coeff) in poly.iter().enumerate() {
        out[target_degree - r] = coeff;
    }
    trim(out)
}

fn summarize_sequence(label: &str, polys: &[Vec<i64>], d_start: usize) {
    let real_rooted_count = polys
        .iter()
        .filter(|p| p.len() <= 2 || is_real_rooted(p))
        .count();

    let mut forward_ok = 0usize;
    let mut backward_ok = 0usize;
    let mut first_forward_failure = None;
    let mut first_backward_failure = None;

    for (idx, pair) in polys.windows(2).enumerate() {
        if interlaces_weak(&pair[0], &pair[1]) {
            forward_ok += 1;
        } else if first_forward_failure.is_none() {
            first_forward_failure = Some((d_start + idx, pair[0].clone(), pair[1].clone()));
        }

        if interlaces_weak(&pair[1], &pair[0]) {
            backward_ok += 1;
        } else if first_backward_failure.is_none() {
            first_backward_failure = Some((d_start + idx, pair[0].clone(), pair[1].clone()));
        }
    }

    println!(
        "  {} real-rooted: {}/{}",
        label,
        real_rooted_count,
        polys.len()
    );
    if polys.len() >= 2 {
        println!(
            "  {} consecutive interlacing: forward {}/{}, backward {}/{}",
            label,
            forward_ok,
            polys.len() - 1,
            backward_ok,
            polys.len() - 1
        );
    }

    if let Some((d, left, right)) = first_forward_failure {
        println!(
            "  {} first forward failure at d={} -> d+1={}: {} / {}",
            label,
            d,
            d + 1,
            format_poly(&left),
            format_poly(&right)
        );
    }
    if let Some((d, left, right)) = first_backward_failure {
        println!(
            "  {} first backward failure at d={} -> d+1={}: {} / {}",
            label,
            d,
            d + 1,
            format_poly(&right),
            format_poly(&left)
        );
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let max_m: usize = args.first().and_then(|s| s.parse().ok()).unwrap_or(8);
    let max_d: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(20);
    let run_recurrence = args.iter().any(|arg| arg == "--recurrence");

    let search = AdaptiveSearchOptions {
        max_rec_len: 5,
        max_var_deg: 2,
        max_idx_deg: 2,
        max_diff_deg: 1,
        try_inhomogeneous: true,
        try_denominator: true,
        max_denom_var_deg: 1,
        max_denom_idx_deg: 2,
        min_margin: 2,
        verbose: false,
    };

    println!("=== Fixed one-descent reversed-in-d experiment ===");
    println!("m in [1, {}], d in [1, {}]", max_m, max_d);
    println!();

    for m in 1..=max_m {
        let l_reversed: Vec<Vec<i64>> = (1..=max_d)
            .map(|d| reverse_to_degree(&one_descent_polynomial(d, m), d))
            .collect();
        let b_reversed: Vec<Vec<i64>> = (1..=max_d)
            .map(|d| reverse_to_degree(&one_descent_unperturbed_polynomial(d, m), d))
            .collect();

        println!("m={}", m);
        println!(
            "  sample reversed L: d=1 {}, d={} {}, d={} {}",
            format_poly(&l_reversed[0]),
            (max_d / 2).max(1),
            format_poly(&l_reversed[(max_d / 2).saturating_sub(1)]),
            max_d,
            format_poly(&l_reversed[max_d - 1]),
        );

        summarize_sequence("R(L)", &l_reversed, 1);
        summarize_sequence("R(B)", &b_reversed, 1);

        if run_recurrence {
            match find_recurrence_adaptive(&l_reversed, &search) {
                Some(result) => {
                    println!(
                        "  R(L) recurrence: {}  [rec_len={}, var_deg={}, idx_deg={}, diff_deg={}, tried={}]",
                        result.recurrence,
                        result.opts.rec_len,
                        result.opts.var_deg,
                        result.opts.idx_deg,
                        result.opts.diff_deg,
                        result.candidates_tried
                    );
                }
                None => println!("  R(L) recurrence: none found in search bounds"),
            }

            match find_recurrence_adaptive(&b_reversed, &search) {
                Some(result) => {
                    println!(
                        "  R(B) recurrence: {}  [rec_len={}, var_deg={}, idx_deg={}, diff_deg={}, tried={}]",
                        result.recurrence,
                        result.opts.rec_len,
                        result.opts.var_deg,
                        result.opts.idx_deg,
                        result.opts.diff_deg,
                        result.candidates_tried
                    );
                }
                None => println!("  R(B) recurrence: none found in search bounds"),
            }
        }

        println!();
    }
}
