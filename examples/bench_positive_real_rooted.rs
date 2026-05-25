//! Benchmark exact real-rootedness checks for positive-coefficient polynomials.
//!
//! Run with:
//!
//! ```text
//! cargo run --release -p polynomial-tools --example bench_positive_real_rooted
//! ```

use num_bigint::BigInt;
use polynomial_tools::{
    is_real_rooted_bezout_bigint_coeffs, is_real_rooted_bigint_coeffs,
    is_real_rooted_one_signed_bigint_coeffs, is_real_rooted_prs_bigint_coeffs,
    primitive_sturm_max_coefficient_bits,
};
use std::hint::black_box;
use std::time::{Duration, Instant};

fn trim(mut p: Vec<BigInt>) -> Vec<BigInt> {
    while p.last().is_some_and(|c| c == &BigInt::from(0)) {
        p.pop();
    }
    p
}

fn add(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let n = a.len().max(b.len());
    let mut r = vec![BigInt::from(0); n];
    for (i, c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, c) in b.iter().enumerate() {
        r[i] += c;
    }
    trim(r)
}

fn mul(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![BigInt::from(0); a.len() + b.len() - 1];
    for (i, ca) in a.iter().enumerate() {
        for (j, cb) in b.iter().enumerate() {
            r[i + j] += ca * cb;
        }
    }
    trim(r)
}

fn product_linear(deg: usize) -> Vec<BigInt> {
    let mut p = vec![BigInt::from(1)];
    for a in 1..=deg {
        p = mul(&p, &[BigInt::from(a), BigInt::from(1)]);
    }
    p
}

fn eulerian_polynomial(n: usize) -> Vec<BigInt> {
    let mut p = vec![BigInt::from(1)];
    for m in 2..=n {
        let mut dp = vec![BigInt::from(0); p.len().saturating_sub(1)];
        for k in 1..p.len() {
            dp[k - 1] = &p[k] * BigInt::from(k);
        }

        let mut term1 = vec![BigInt::from(0); p.len() + 1];
        for k in 0..p.len() {
            term1[k] += &p[k];
            term1[k + 1] += &p[k] * BigInt::from(m - 1);
        }

        let mut term2 = vec![BigInt::from(0); dp.len() + 2];
        for k in 0..dp.len() {
            term2[k + 1] += &dp[k];
            term2[k + 2] -= &dp[k];
        }
        p = add(&term1, &term2);
    }
    p
}

fn time_it(name: &str, f: impl FnOnce() -> bool) -> (bool, Duration) {
    let t0 = Instant::now();
    let result = black_box(f());
    let elapsed = t0.elapsed();
    println!(
        "  {:<18} {:<5} {:>10.3} ms",
        name,
        result,
        elapsed.as_secs_f64() * 1000.0
    );
    (result, elapsed)
}

fn time_silent(f: impl FnOnce() -> bool) -> (bool, Duration) {
    let t0 = Instant::now();
    let result = black_box(f());
    (result, t0.elapsed())
}

fn bench_case(name: &str, p: &[BigInt], run_bezout: bool) {
    let degree = p.len().saturating_sub(1);
    let coeff_bits = p.iter().map(|c| c.bits()).max().unwrap_or(0);
    let prs_bits = primitive_sturm_max_coefficient_bits(p);
    println!();
    println!("{name}: degree={degree}, max_coeff_bits={coeff_bits}, max_prs_bits={prs_bits}");
    let (fast, _) = time_it("default-fast", || is_real_rooted_bigint_coeffs(p));
    let (one_signed, _) = time_it("one-signed", || {
        is_real_rooted_one_signed_bigint_coeffs(p).expect("case should be one-signed")
    });
    let (prs, _) = time_it("primitive-prs", || is_real_rooted_prs_bigint_coeffs(p));
    assert_eq!(fast, one_signed);
    assert_eq!(fast, prs);

    if run_bezout {
        let (bezout, _) = time_it("bezout-bigint", || is_real_rooted_bezout_bigint_coeffs(p));
        assert_eq!(fast, bezout);
    } else {
        println!("  {:<18} skipped", "bezout-bigint");
    }
}

fn cutoff_sweep(name: &str, degrees: &[usize], make_poly: impl Fn(usize) -> Vec<BigInt>) {
    println!();
    println!("{name} cutoff sweep");
    println!(
        "{:>6} {:>8} {:>14} {:>14} {:>14} {:>10}",
        "deg", "bits", "default ms", "prs ms", "bezout ms", "winner"
    );
    println!("{}", "-".repeat(76));

    for &degree in degrees {
        let p = make_poly(degree);
        let bits = p.iter().map(|c| c.bits()).max().unwrap_or(0);
        let (default_rr, default_t) = time_silent(|| is_real_rooted_bigint_coeffs(&p));
        let (prs_rr, prs_t) = time_silent(|| is_real_rooted_prs_bigint_coeffs(&p));
        assert_eq!(default_rr, prs_rr);

        let run_bezout = degree <= 45;
        let (bezout_rr, bezout_t) = if run_bezout {
            let r = time_silent(|| is_real_rooted_bezout_bigint_coeffs(&p));
            assert_eq!(default_rr, r.0);
            (Some(r.0), Some(r.1))
        } else {
            (None, None)
        };
        black_box(bezout_rr);

        let default_ms = default_t.as_secs_f64() * 1000.0;
        let prs_ms = prs_t.as_secs_f64() * 1000.0;
        let bezout_ms = bezout_t.map(|t| t.as_secs_f64() * 1000.0);
        let winner = match bezout_ms {
            Some(b) if b < default_ms && b < prs_ms => "bezout",
            _ if prs_ms < default_ms => "prs",
            _ => "default",
        };

        let bezout_display = bezout_ms
            .map(|ms| format!("{ms:14.3}"))
            .unwrap_or_else(|| format!("{:>14}", "skipped"));
        println!(
            "{degree:>6} {bits:>8} {default_ms:>14.3} {prs_ms:>14.3} {bezout_display} {winner:>10}",
        );
    }
}

fn main() {
    println!("Exact real-rootedness benchmark for positive coefficients");
    bench_case("prod_{a=1}^{30} (x+a)", &product_linear(30), true);
    bench_case("prod_{a=1}^{80} (x+a)", &product_linear(80), false);
    bench_case("Eulerian A_35", &eulerian_polynomial(35), true);
    bench_case("Eulerian A_80", &eulerian_polynomial(80), false);
    bench_case(
        "small non-real-rooted",
        &[
            BigInt::from(1),
            BigInt::from(43),
            BigInt::from(196),
            BigInt::from(168),
            BigInt::from(23),
            BigInt::from(1),
        ],
        true,
    );

    cutoff_sweep(
        "prod_{a=1}^{d} (x+a)",
        &[5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80],
        product_linear,
    );
    cutoff_sweep("Eulerian A_{d+1}", &[5, 10, 15, 20, 25, 30, 35, 40], |d| {
        eulerian_polynomial(d + 1)
    });
}
