//! Explore the normalized reversed one-descent family
//!
//!   Q_d^(m)(x) = x^m L_m^(d)(1/x)
//!              = sum_r binom(m,r) binom(d,r) x^(m-r) - x^(m-1),
//!
//! and the unperturbed companion
//!
//!   U_d^(m)(x) = x^m B_m^(d)(1/x)
//!              = sum_r binom(m,r) binom(d,r) x^(m-r).
//!
//! This removes the automatic zero roots present in x^d L_m^(d)(1/x) when d>m.
//! We test:
//! - real-rootedness,
//! - forward interlacing in d,
//! - directed interlacing for the step Delta_d = Q_(d+1) - Q_d,
//! - derivative relations across m.

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

fn poly_sub(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut out = vec![0; len];
    for i in 0..len {
        let ai = a.get(i).copied().unwrap_or(0);
        let bi = b.get(i).copied().unwrap_or(0);
        out[i] = ai - bi;
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

fn add_monomial(p: &[i64], degree: usize, coeff: i64) -> Vec<i64> {
    let mut out = p.to_vec();
    if out.len() <= degree {
        out.resize(degree + 1, 0);
    }
    out[degree] += coeff;
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
    } else if df == dg + 1 {
        check_weak_interlacing(&g, &f) == Some(true)
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

fn normalized_unperturbed(m: usize, d: usize) -> Vec<i64> {
    let mut coeffs = vec![0; m + 1];
    for r in 0..=m {
        coeffs[m - r] = binomial_i64(m, r) * binomial_i64(d, r);
    }
    trim(coeffs)
}

fn normalized_perturbed(m: usize, d: usize) -> Vec<i64> {
    let mut coeffs = normalized_unperturbed(m, d);
    if m >= 1 {
        coeffs[m - 1] -= 1;
    }
    trim(coeffs)
}

fn normalized_step(m: usize, d: usize) -> Vec<i64> {
    let mut coeffs = vec![0; m];
    for r in 1..=m {
        coeffs[m - r] = binomial_i64(m, r) * binomial_i64(d, r - 1);
    }
    trim(coeffs)
}

fn summarize_forward(label: &str, polys: &[Vec<i64>], d_start: usize) {
    let real_rooted_count = polys
        .iter()
        .filter(|p| p.len() <= 2 || is_real_rooted(p))
        .count();

    let mut forward_ok = 0usize;
    let mut first_forward_failure = None;
    for (idx, pair) in polys.windows(2).enumerate() {
        if interlaces_weak(&pair[0], &pair[1]) {
            forward_ok += 1;
        } else if first_forward_failure.is_none() {
            first_forward_failure = Some((d_start + idx, pair[0].clone(), pair[1].clone()));
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
            "  {} consecutive interlacing: forward {}/{}",
            label,
            forward_ok,
            polys.len() - 1,
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
}

fn summarize_step(
    label: &str,
    current: &[Vec<i64>],
    next: &[Vec<i64>],
    steps: &[Vec<i64>],
    d_start: usize,
) {
    let mut step_real_rooted = 0usize;
    let mut step_into_current = 0usize;
    let mut step_into_next = 0usize;
    let mut current_compatible = 0usize;
    let mut first_current_failure = None;
    let mut first_next_failure = None;

    for idx in 0..steps.len() {
        if steps[idx].len() <= 2 || is_real_rooted(&steps[idx]) {
            step_real_rooted += 1;
        }
        if interlaces_weak(&steps[idx], &current[idx]) {
            step_into_current += 1;
        } else if first_current_failure.is_none() {
            first_current_failure = Some((d_start + idx, steps[idx].clone(), current[idx].clone()));
        }
        if interlaces_weak(&steps[idx], &next[idx]) {
            step_into_next += 1;
        } else if first_next_failure.is_none() {
            first_next_failure = Some((d_start + idx, steps[idx].clone(), next[idx].clone()));
        }
        if is_real_rooted(&trim(current[idx].clone()))
            && is_real_rooted(&trim(steps[idx].clone()))
            && is_real_rooted(&trim({
                let len = current[idx].len().max(steps[idx].len());
                let mut sum = vec![0; len];
                for i in 0..len {
                    sum[i] = current[idx].get(i).copied().unwrap_or(0)
                        + steps[idx].get(i).copied().unwrap_or(0);
                }
                sum
            }))
        {
            current_compatible += 1;
        }
    }

    println!(
        "  {} step real-rooted: {}/{}",
        label,
        step_real_rooted,
        steps.len()
    );
    println!(
        "  {} step interlaces current/next: {}/{}, {}/{}",
        label,
        step_into_current,
        steps.len(),
        step_into_next,
        steps.len()
    );
    println!(
        "  {} current + step real-rooted: {}/{}",
        label,
        current_compatible,
        steps.len()
    );
    if let Some((d, step, cur)) = first_current_failure {
        println!(
            "  {} first step -> current failure at d={}: {} / {}",
            label,
            d,
            format_poly(&step),
            format_poly(&cur)
        );
    }
    if let Some((d, step, nxt)) = first_next_failure {
        println!(
            "  {} first step -> next failure at d={}: {} / {}",
            label,
            d,
            format_poly(&step),
            format_poly(&nxt)
        );
    }
}

fn summarize_derivative_relations(max_m: usize, max_d: usize) {
    let mut u_exact = 0usize;
    let mut q_exact = 0usize;
    let mut delta_exact = 0usize;
    let mut u_interlacing = 0usize;
    let mut q_interlacing = 0usize;
    let mut delta_interlacing = 0usize;
    let mut u_prev_into_delta = 0usize;
    let mut delta_into_u_prev = 0usize;
    let mut q_prev_into_delta = 0usize;
    let mut delta_into_q_prev = 0usize;
    let total = max_d * max_m.saturating_sub(1);

    for m in 2..=max_m {
        for d in 1..=max_d {
            let u = normalized_unperturbed(m, d);
            let u_prev = normalized_unperturbed(m - 1, d);
            let u_der = derivative(&u);
            let u_rhs = u_prev.iter().map(|&c| c * (m as i64)).collect::<Vec<_>>();
            if trim(u_der.clone()) == trim(u_rhs) {
                u_exact += 1;
            }
            if interlaces_weak(&u_prev, &u) {
                u_interlacing += 1;
            }

            let q = normalized_perturbed(m, d);
            let q_prev = normalized_perturbed(m - 1, d);
            let q_der = derivative(&q);
            let q_rhs = add_monomial(
                &q_prev.iter().map(|&c| c * (m as i64)).collect::<Vec<_>>(),
                m - 2,
                1,
            );
            if trim(q_der.clone()) == trim(q_rhs) {
                q_exact += 1;
            }
            if interlaces_weak(&q_prev, &q) {
                q_interlacing += 1;
            }

            let delta = normalized_step(m, d);
            let delta_prev = normalized_step(m - 1, d);
            let delta_der = derivative(&delta);
            let delta_rhs = delta_prev
                .iter()
                .map(|&c| c * (m as i64))
                .collect::<Vec<_>>();
            if trim(delta_der) == trim(delta_rhs) {
                delta_exact += 1;
            }
            if interlaces_weak(&delta_prev, &delta) {
                delta_interlacing += 1;
            }
            if interlaces_weak(&u_prev, &delta) {
                u_prev_into_delta += 1;
            }
            if interlaces_weak(&delta, &u_prev) {
                delta_into_u_prev += 1;
            }
            if interlaces_weak(&q_prev, &delta) {
                q_prev_into_delta += 1;
            }
            if interlaces_weak(&delta, &q_prev) {
                delta_into_q_prev += 1;
            }
        }
    }

    println!("Derivative relations across m");
    println!("  U'_d^(m) = m U_d^(m-1): {}/{}", u_exact, total);
    println!("  Q'_d^(m) = m Q_d^(m-1) + x^(m-2): {}/{}", q_exact, total);
    println!(
        "  Delta'_d^(m) = m Delta_d^(m-1): {}/{}",
        delta_exact, total
    );
    println!("  U_d^(m-1) << U_d^(m): {}/{}", u_interlacing, total);
    println!("  Q_d^(m-1) << Q_d^(m): {}/{}", q_interlacing, total);
    println!(
        "  Delta_d^(m-1) << Delta_d^(m): {}/{}",
        delta_interlacing, total
    );
    println!(
        "  U_d^(m-1) << Delta_d^(m): {}/{}",
        u_prev_into_delta, total
    );
    println!(
        "  Delta_d^(m) << U_d^(m-1): {}/{}",
        delta_into_u_prev, total
    );
    println!(
        "  Q_d^(m-1) << Delta_d^(m): {}/{}",
        q_prev_into_delta, total
    );
    println!(
        "  Delta_d^(m) << Q_d^(m-1): {}/{}",
        delta_into_q_prev, total
    );
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let max_m: usize = args.first().and_then(|s| s.parse().ok()).unwrap_or(8);
    let max_d: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(20);

    println!("=== Normalized reversed one-descent experiment ===");
    println!("m in [1, {}], d in [1, {}]", max_m, max_d);
    println!();

    for m in 1..=max_m {
        let q: Vec<Vec<i64>> = (1..=max_d).map(|d| normalized_perturbed(m, d)).collect();
        let u: Vec<Vec<i64>> = (1..=max_d).map(|d| normalized_unperturbed(m, d)).collect();
        let q_steps: Vec<Vec<i64>> = q
            .windows(2)
            .map(|pair| poly_sub(&pair[1], &pair[0]))
            .collect();
        let u_steps: Vec<Vec<i64>> = u
            .windows(2)
            .map(|pair| poly_sub(&pair[1], &pair[0]))
            .collect();

        println!("m={}", m);
        println!(
            "  sample Q: d=1 {}, d={} {}, d={} {}",
            format_poly(&q[0]),
            (max_d / 2).max(1),
            format_poly(&q[(max_d / 2).saturating_sub(1)]),
            max_d,
            format_poly(&q[max_d - 1]),
        );
        summarize_forward("Q", &q, 1);
        summarize_forward("U", &u, 1);
        summarize_forward("Delta(Q)", &q_steps, 1);
        summarize_forward("Delta(U)", &u_steps, 1);
        summarize_step("Q", &q[..q.len() - 1], &q[1..], &q_steps, 1);
        summarize_step("U", &u[..u.len() - 1], &u[1..], &u_steps, 1);
        println!();
    }

    summarize_derivative_relations(max_m, max_d);
}
