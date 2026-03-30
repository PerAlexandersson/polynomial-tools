//! Benchmark: Bézout matrix vs Sturm chain for real-rootedness and interlacing.
//!
//! Run with: cargo run --release --example bench_bezout_vs_sturm

use polynomial_tools::real_rootedness::*;
use std::time::Instant;

/// Compute Eulerian polynomials A_1(t), A_2(t), ..., A_n(t) using the recurrence
/// A_n(t) = (1 + (n-1)t) A_{n-1}(t) + t(1-t) A'_{n-1}(t).
fn eulerian_polynomials(max_n: usize) -> Vec<Vec<i64>> {
    let mut polys = Vec::new();
    polys.push(vec![1i64]); // A_1 = 1

    for n in 2..=max_n {
        let prev = &polys[n - 2];
        let d = prev.len();

        // Derivative of prev
        let mut dp = vec![0i64; d.saturating_sub(1)];
        for k in 1..d {
            dp[k - 1] = prev[k] * k as i64;
        }

        // (1 + (n-1)t) * prev
        let mut term1 = vec![0i64; d + 1];
        for k in 0..d {
            term1[k] += prev[k];
            term1[k + 1] += prev[k] * (n as i64 - 1);
        }

        // t(1-t) * dp = t*dp - t^2*dp
        let dp_len = dp.len();
        let mut term2 = vec![0i64; dp_len + 2];
        for k in 0..dp_len {
            term2[k + 1] += dp[k];     // t * dp
            term2[k + 2] -= dp[k];     // -t^2 * dp
        }

        // Sum
        let len = term1.len().max(term2.len());
        let mut result = vec![0i64; len];
        for k in 0..term1.len() {
            result[k] += term1[k];
        }
        for k in 0..term2.len() {
            result[k] += term2[k];
        }
        // Trim trailing zeros
        while result.last() == Some(&0) {
            result.pop();
        }
        polys.push(result);
    }
    polys
}

fn main() {
    let max_n = 30;
    let polys = eulerian_polynomials(max_n);

    println!("=== Real-rootedness: Sturm vs Bézout ===");
    println!("{:>4} {:>8} {:>12} {:>12}", "n", "deg", "Sturm (μs)", "Bézout (μs)");
    println!("{}", "-".repeat(42));

    for (i, p) in polys.iter().enumerate() {
        let n = i + 1;
        if n < 3 {
            continue;
        }

        let iters = if n < 10 { 100 } else if n < 20 { 10 } else { 3 };

        let t0 = Instant::now();
        for _ in 0..iters {
            let _ = is_real_rooted_sturm(p);
        }
        let sturm_us = t0.elapsed().as_micros() as f64 / iters as f64;

        let t0 = Instant::now();
        for _ in 0..iters {
            let _ = is_real_rooted(p);
        }
        let bezout_us = t0.elapsed().as_micros() as f64 / iters as f64;

        let deg = p.len() - 1;
        println!(
            "{:>4} {:>8} {:>12.1} {:>12.1}  {}",
            n,
            deg,
            sturm_us,
            bezout_us,
            if bezout_us < sturm_us { "← Bézout" } else { "← Sturm" }
        );
    }

    println!();
    println!("=== Interlacing: Sturm vs Bézout (consecutive Eulerian polynomials) ===");
    println!("{:>4} {:>8} {:>12} {:>12}", "n", "deg", "Sturm (μs)", "Bézout (μs)");
    println!("{}", "-".repeat(42));

    for i in 1..polys.len() {
        let n = i + 1;
        if n < 4 {
            continue;
        }
        let f = &polys[i];   // A_{n}(t), degree n-1
        let g = &polys[i - 1]; // A_{n-1}(t), degree n-2

        let iters = if n < 10 { 100 } else if n < 20 { 10 } else { 3 };

        let t0 = Instant::now();
        for _ in 0..iters {
            let _ = check_interlacing_sturm(f, g);
        }
        let sturm_us = t0.elapsed().as_micros() as f64 / iters as f64;

        let t0 = Instant::now();
        for _ in 0..iters {
            let _ = check_interlacing(f, g);
        }
        let bezout_us = t0.elapsed().as_micros() as f64 / iters as f64;

        let deg = f.len() - 1;
        let result = check_interlacing(f, g);
        println!(
            "{:>4} {:>8} {:>12.1} {:>12.1}  {} (interlaces: {:?})",
            n,
            deg,
            sturm_us,
            bezout_us,
            if bezout_us < sturm_us { "← Bézout" } else { "← Sturm" },
            result
        );
    }
}
