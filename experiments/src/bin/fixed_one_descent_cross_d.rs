//! Explore the cross-d Pascal step for the fixed one-descent family.
//!
//! The key identities are
//!
//!   B_{m+1}^{(d)} = B_m^{(d)} + d * Integral(B_m^{(d-1)}),
//!   L_{m+1}^{(d)} = L_m^{(d)} + d * Integral(B_m^{(d-1)}),
//!
//! Here
//!   B_m^{(d)}(t) = sum_r binom(d,r) binom(m,r) t^r,
//!   L_m^{(d)}(t) = B_m^{(d)}(t) - t.
//!
//! This binary checks:
//! - the exact identities,
//! - real-rootedness of all pieces,
//! - interlacing/compatibility between the summands,
//! - whether the integral cross-d step plausibly explains the forward interlacing
//!   L_m^{(d)} << L_{m+1}^{(d)}.

use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly, is_real_rooted};

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

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let mut out = vec![0; a.len().max(b.len())];
    for (i, &c) in a.iter().enumerate() {
        out[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        out[i] += c;
    }
    trim(out)
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

fn compatible(f: &[i64], g: &[i64]) -> bool {
    if f == [0] || g == [0] {
        return true;
    }
    let sum = poly_add(f, g);
    sum.len() <= 2 || is_real_rooted(&sum)
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

fn b_polynomial(d: usize, m: usize) -> Vec<i64> {
    let degree = d.min(m);
    let mut coeffs = vec![0; degree + 1];
    for (r, coeff) in coeffs.iter_mut().enumerate().take(degree + 1) {
        *coeff = binomial_i64(d, r) * binomial_i64(m, r);
    }
    trim(coeffs)
}

fn l_polynomial(d: usize, m: usize) -> Vec<i64> {
    let mut coeffs = b_polynomial(d, m);
    if coeffs.len() < 2 {
        coeffs.resize(2, 0);
    }
    coeffs[1] -= 1;
    trim(coeffs)
}

fn cross_d_step_polynomial(d: usize, m: usize) -> Vec<i64> {
    let degree = d.min(m + 1);
    let mut coeffs = vec![0; degree + 1];
    for (r, coeff) in coeffs.iter_mut().enumerate().take(degree + 1).skip(1) {
        *coeff = binomial_i64(d, r) * binomial_i64(m, r - 1);
    }
    trim(coeffs)
}

fn main() {
    let max_d: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(10);
    let max_m: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(20);

    let mut identity_ok = 0usize;
    let mut total_cases = 0usize;

    let mut l_to_next_counts = (0usize, 0usize);
    let mut l_vs_step_counts = (0usize, 0usize);
    let mut step_vs_lnext_counts = (0usize, 0usize);
    let mut b_vs_step_counts = (0usize, 0usize);
    let mut compat_l_step_counts = (0usize, 0usize);
    let mut compat_b_step_counts = (0usize, 0usize);
    let mut b_prev_vs_b_counts = (0usize, 0usize);
    let mut step_real_rooted_counts = (0usize, 0usize);

    let mut first_failures: Vec<String> = Vec::new();

    for d in 1..=max_d {
        for m in 1..=max_m {
            let l_dm = l_polynomial(d, m);
            let l_dmp1 = l_polynomial(d, m + 1);
            let b_dm = b_polynomial(d, m);
            let b_dmp1 = b_polynomial(d, m + 1);
            let b_dminus1_m = if d >= 1 {
                b_polynomial(d - 1, m)
            } else {
                vec![0]
            };
            let step_term = cross_d_step_polynomial(d, m);

            total_cases += 1;

            let b_identity = poly_sub(&b_dmp1, &poly_add(&b_dm, &step_term)) == [0];
            let l_identity = poly_sub(&l_dmp1, &poly_add(&l_dm, &step_term)) == [0];
            if b_identity && l_identity {
                identity_ok += 1;
            } else if first_failures.is_empty() {
                first_failures.push(format!(
                    "identity fail at d={}, m={}: B_next={}, B+step={}; L_next={}, L+step={}",
                    d,
                    m,
                    format_poly(&b_dmp1),
                    format_poly(&poly_add(&b_dm, &step_term)),
                    format_poly(&l_dmp1),
                    format_poly(&poly_add(&l_dm, &step_term)),
                ));
            }

            l_to_next_counts.0 += 1;
            if interlaces_weak(&l_dm, &l_dmp1) {
                l_to_next_counts.1 += 1;
            } else if first_failures.len() < 6 {
                first_failures.push(format!(
                    "L_m << L_(m+1) fail at d={}, m={}: {} / {}",
                    d,
                    m,
                    format_poly(&l_dm),
                    format_poly(&l_dmp1)
                ));
            }

            l_vs_step_counts.0 += 1;
            if interlaces_weak(&l_dm, &step_term) {
                l_vs_step_counts.1 += 1;
            } else if first_failures.len() < 6 {
                first_failures.push(format!(
                    "L_m << step_term fail at d={}, m={}: {} / {}",
                    d,
                    m,
                    format_poly(&l_dm),
                    format_poly(&step_term)
                ));
            }

            step_vs_lnext_counts.0 += 1;
            if interlaces_weak(&step_term, &l_dmp1) {
                step_vs_lnext_counts.1 += 1;
            } else if first_failures.len() < 6 {
                first_failures.push(format!(
                    "step_term << L_(m+1) fail at d={}, m={}: {} / {}",
                    d,
                    m,
                    format_poly(&step_term),
                    format_poly(&l_dmp1)
                ));
            }

            b_vs_step_counts.0 += 1;
            if interlaces_weak(&b_dm, &step_term) {
                b_vs_step_counts.1 += 1;
            }

            compat_l_step_counts.0 += 1;
            if compatible(&l_dm, &step_term) {
                compat_l_step_counts.1 += 1;
            }

            compat_b_step_counts.0 += 1;
            if compatible(&b_dm, &step_term) {
                compat_b_step_counts.1 += 1;
            }

            b_prev_vs_b_counts.0 += 1;
            if d >= 1 && interlaces_weak(&b_dminus1_m, &b_dm) {
                b_prev_vs_b_counts.1 += 1;
            }

            step_real_rooted_counts.0 += 1;
            if step_term.len() <= 2 || is_real_rooted(&step_term) {
                step_real_rooted_counts.1 += 1;
            }
        }
    }

    println!("=== Fixed one-descent cross-d experiment ===");
    println!("d in [1, {}], m in [1, {}]", max_d, max_m);
    println!();
    println!(
        "exact identities with step_term = d * Integral(B_(d-1,m)): {}/{}",
        identity_ok, total_cases
    );
    println!(
        "L_m << L_(m+1): {}/{}",
        l_to_next_counts.1, l_to_next_counts.0
    );
    println!(
        "L_m << step_term: {}/{}",
        l_vs_step_counts.1, l_vs_step_counts.0
    );
    println!(
        "step_term << L_(m+1): {}/{}",
        step_vs_lnext_counts.1, step_vs_lnext_counts.0
    );
    println!(
        "B_m << step_term: {}/{}",
        b_vs_step_counts.1, b_vs_step_counts.0
    );
    println!(
        "B_(d-1,m) << B_(d,m): {}/{}",
        b_prev_vs_b_counts.1, b_prev_vs_b_counts.0
    );
    println!(
        "compatibility of L_m and step_term: {}/{}",
        compat_l_step_counts.1, compat_l_step_counts.0
    );
    println!(
        "compatibility of B_m and step_term: {}/{}",
        compat_b_step_counts.1, compat_b_step_counts.0
    );
    println!(
        "real-rootedness of step_term: {}/{}",
        step_real_rooted_counts.1, step_real_rooted_counts.0
    );
    println!();

    if !first_failures.is_empty() {
        println!("Sample failures:");
        for line in first_failures {
            println!("  {}", line);
        }
    }
}
