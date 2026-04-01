//! Explore the fixed-one-descent family
//!
//!   L_m^(d)(t) = sum_r binom(d,r) binom(m,r) t^r - t,
//!
//! where m = n-d and d is fixed.
//!
//! We compare this against the unperturbed family
//!
//!   B_m^(d)(t) = L_m^(d)(t) + t
//!              = sum_r binom(d,r) binom(m,r) t^r,
//!
//! and test:
//! - real-rootedness,
//! - consecutive interlacing in m,
//! - adaptive polynomial recurrences in m.

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
    let value = num / den;
    i64::try_from(value).expect("binomial coefficient overflowed i64")
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

fn poly_sub(a: &[i64], b: &[i64]) -> Vec<i64> {
    let mut out = vec![0; a.len().max(b.len())];
    for (i, &c) in a.iter().enumerate() {
        out[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        out[i] -= c;
    }
    trim(out)
}

fn poly_scale(a: &[i64], scalar: i64) -> Vec<i64> {
    trim(a.iter().map(|&c| c * scalar).collect())
}

fn verify_binomial_constant_recurrence(polys: &[Vec<i64>], order: usize) -> bool {
    if polys.len() <= order {
        return false;
    }
    for n in order..polys.len() {
        let mut lhs = polys[n].clone();
        for j in 1..=order {
            let sign = if j % 2 == 1 { 1 } else { -1 };
            let coeff = sign * binomial_i64(order, j);
            lhs = poly_sub(&lhs, &poly_scale(&polys[n - j], coeff));
        }
        if lhs != [0] {
            return false;
        }
    }
    true
}

fn summarize_sequence(label: &str, polys: &[Vec<i64>], m_start: usize) {
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
            first_forward_failure = Some((m_start + idx, pair[0].clone(), pair[1].clone()));
        }

        if interlaces_weak(&pair[1], &pair[0]) {
            backward_ok += 1;
        } else if first_backward_failure.is_none() {
            first_backward_failure = Some((m_start + idx, pair[0].clone(), pair[1].clone()));
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

    if let Some((m, left, right)) = first_forward_failure {
        println!(
            "  {} first forward failure at m={} -> m+1={}: {}  /  {}",
            label,
            m,
            m + 1,
            format_poly(&left),
            format_poly(&right)
        );
    }
    if let Some((m, left, right)) = first_backward_failure {
        println!(
            "  {} first backward failure at m={} -> m+1={}: {}  /  {}",
            label,
            m,
            m + 1,
            format_poly(&right),
            format_poly(&left)
        );
    }
}

fn main() {
    let max_d: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);
    let max_m: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(20);

    println!("=== Fixed one-descent sequence experiment ===");
    println!("d in [1, {}], m in [1, {}]", max_d, max_m);
    println!();

    for d in 1..=max_d {
        let l_polys: Vec<Vec<i64>> = (1..=max_m).map(|m| one_descent_polynomial(d, m)).collect();
        let b_polys: Vec<Vec<i64>> = (1..=max_m)
            .map(|m| one_descent_unperturbed_polynomial(d, m))
            .collect();

        println!("d={}", d);
        println!(
            "  sample L: m=1 {}, m={} {}, m={} {}",
            format_poly(&l_polys[0]),
            (max_m / 2).max(1),
            format_poly(&l_polys[(max_m / 2).saturating_sub(1)]),
            max_m,
            format_poly(&l_polys[max_m - 1]),
        );

        summarize_sequence("L", &l_polys, 1);
        summarize_sequence("B", &b_polys, 1);

        let expected_order = d + 1;
        println!(
            "  exact binomial recurrence of order {}: L={}, B={}",
            expected_order,
            verify_binomial_constant_recurrence(&l_polys, expected_order),
            verify_binomial_constant_recurrence(&b_polys, expected_order),
        );

        let simple_search = AdaptiveSearchOptions {
            max_rec_len: expected_order,
            max_var_deg: 0,
            max_idx_deg: 0,
            max_diff_deg: 0,
            try_inhomogeneous: false,
            try_denominator: false,
            min_margin: 1,
            verbose: false,
            ..Default::default()
        };

        match find_recurrence_adaptive(&l_polys, &simple_search) {
            Some(result) => {
                println!(
                    "  L simple recurrence: {}  [rec_len={}, tried={}]",
                    result.recurrence, result.opts.rec_len, result.candidates_tried
                );
            }
            None => println!("  L simple recurrence: none found"),
        }

        match find_recurrence_adaptive(&b_polys, &simple_search) {
            Some(result) => {
                println!(
                    "  B simple recurrence: {}  [rec_len={}, tried={}]",
                    result.recurrence, result.opts.rec_len, result.candidates_tried
                );
            }
            None => println!("  B simple recurrence: none found"),
        }

        println!();
    }
}
