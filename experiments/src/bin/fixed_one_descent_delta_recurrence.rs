//! Explore recurrences for the normalized one-descent step family.
//!
//! For fixed d, define
//!
//!   Delta_d^(m)(x) = sum_{r=1}^m binom(m,r) binom(d,r-1) x^(m-r).
//!
//! This has degree m-1, so varying m gives a genuine degree-growing sequence.
//! We test:
//! - real-rootedness and forward interlacing in m,
//! - the exact derivative recurrence
//!     d/dx Delta_d^(m+1) = (m+1) Delta_d^(m),
//! - the exact mixed recurrence
//!     Delta_d^(m+1) = x Delta_d^(m) + U_d^(m),
//!   where U_d^(m)(x) = sum_r binom(m,r) binom(d,r) x^(m-r),
//! - adaptive recurrence search.

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

fn poly_mul_x(p: &[i64]) -> Vec<i64> {
    let mut out = vec![0; p.len() + 1];
    for (i, &c) in p.iter().enumerate() {
        out[i + 1] = c;
    }
    trim(out)
}

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut out = vec![0; len];
    for i in 0..len {
        out[i] = a.get(i).copied().unwrap_or(0) + b.get(i).copied().unwrap_or(0);
    }
    trim(out)
}

fn derivative(p: &[i64]) -> Vec<i64> {
    if p.len() <= 1 {
        return vec![0];
    }
    let mut out = vec![0; p.len() - 1];
    for i in 1..p.len() {
        out[i - 1] = (i as i64) * p[i];
    }
    trim(out)
}

fn scale_poly(p: &[i64], s: i64) -> Vec<i64> {
    trim(p.iter().map(|&c| c * s).collect())
}

fn poly_add_scaled(a: &[i64], sa: i64, b: &[i64], sb: i64) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut out = vec![0; len];
    for i in 0..len {
        out[i] = sa * a.get(i).copied().unwrap_or(0) + sb * b.get(i).copied().unwrap_or(0);
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
        check_weak_interlacing(&g, &poly_mul_x(&f)) == Some(true)
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

fn delta_poly(m: usize, d: usize) -> Vec<i64> {
    let mut coeffs = vec![0; m];
    for r in 1..=m {
        coeffs[m - r] = binomial_i64(m, r) * binomial_i64(d, r - 1);
    }
    trim(coeffs)
}

fn u_poly(m: usize, d: usize) -> Vec<i64> {
    let mut coeffs = vec![0; m + 1];
    for r in 0..=m {
        coeffs[m - r] = binomial_i64(m, r) * binomial_i64(d, r);
    }
    trim(coeffs)
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let max_d: usize = args.first().and_then(|s| s.parse().ok()).unwrap_or(8);
    let max_m: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(16);
    let run_recurrence = args.iter().any(|arg| arg == "--recurrence");

    let search = AdaptiveSearchOptions {
        max_rec_len: 4,
        max_var_deg: 2,
        max_idx_deg: 1,
        max_diff_deg: 1,
        try_inhomogeneous: true,
        try_denominator: true,
        max_denom_var_deg: 1,
        max_denom_idx_deg: 1,
        min_margin: 2,
        ..Default::default()
    };

    println!("=== One-descent step recurrence experiment ===");
    println!("d in [1, {}], m in [1, {}]", max_d, max_m);
    println!();

    for d in 1..=max_d {
        let polys: Vec<Vec<i64>> = (1..=max_m).map(|m| delta_poly(m, d)).collect();
        let ups: Vec<Vec<i64>> = (1..=max_m).map(|m| u_poly(m, d)).collect();
        let real_rooted_count = polys
            .iter()
            .filter(|p| p.len() <= 2 || is_real_rooted(p))
            .count();
        let mut forward_ok = 0usize;
        let mut first_failure = None;
        for (idx, pair) in polys.windows(2).enumerate() {
            if interlaces_weak(&pair[0], &pair[1]) {
                forward_ok += 1;
            } else if first_failure.is_none() {
                first_failure = Some((idx + 1, pair[0].clone(), pair[1].clone()));
            }
        }

        let mut exact_derivative = 0usize;
        let mut exact_mixed = 0usize;
        let mut exact_order2 = 0usize;
        for m in 1..max_m {
            if derivative(&polys[m]) == scale_poly(&polys[m - 1], (m + 1) as i64) {
                exact_derivative += 1;
            }
            if polys[m] == poly_add(&poly_mul_x(&polys[m - 1]), &ups[m - 1]) {
                exact_mixed += 1;
            }
            if m >= 2 {
                let n = (m + 1) as i64;
                let first = scale_poly(&polys[m], n - 1);
                let term1 = poly_add_scaled(
                    &polys[m - 1],
                    (d as i64) + 2 - n,
                    &poly_mul_x(&polys[m - 1]),
                    2 * (n - 1),
                );
                let x_prev2 = poly_mul_x(&polys[m - 2]);
                let x2_prev2 = poly_mul_x(&x_prev2);
                let term2 = poly_add_scaled(&x_prev2, n - 1, &x2_prev2, -(n - 1));
                let rhs = poly_add(&term1, &term2);
                if trim(first) == trim(rhs) {
                    exact_order2 += 1;
                }
            }
        }

        println!("d={}", d);
        println!(
            "  sample step polys: m=1 {}, m={} {}, m={} {}",
            format_poly(&polys[0]),
            (max_m / 2).max(1),
            format_poly(&polys[(max_m / 2).saturating_sub(1)]),
            max_m,
            format_poly(&polys[max_m - 1]),
        );
        println!("  real-rooted: {}/{}", real_rooted_count, polys.len());
        println!(
            "  forward interlacing in m: {}/{}",
            forward_ok,
            polys.len() - 1
        );
        println!(
            "  exact d/dx Delta_(m+1) = (m+1) Delta_m: {}/{}",
            exact_derivative,
            max_m - 1
        );
        println!(
            "  exact Delta_(m+1) = x Delta_m + U_m: {}/{}",
            exact_mixed,
            max_m - 1
        );
        println!(
            "  exact (n-1)P_n = (d+2-n+2(n-1)x)P_(n-1) + (n-1)x(1-x)P_(n-2): {}/{}",
            exact_order2,
            max_m.saturating_sub(2)
        );
        if let Some((m, left, right)) = first_failure {
            println!(
                "  first forward failure at m={} -> m+1={}: {} / {}",
                m,
                m + 1,
                format_poly(&left),
                format_poly(&right)
            );
        }

        if run_recurrence {
            match find_recurrence_adaptive(&polys, &search) {
                Some(result) => println!(
                    "  adaptive recurrence: {}  [rec_len={}, var_deg={}, idx_deg={}, diff_deg={}, tried={}]",
                    result.recurrence,
                    result.opts.rec_len,
                    result.opts.var_deg,
                    result.opts.idx_deg,
                    result.opts.diff_deg,
                    result.candidates_tried
                ),
                None => println!("  adaptive recurrence: none found in search bounds"),
            }
        }

        println!();
    }
}
