//! Test Neville TNN check against brute-force, and benchmark.
//! Also test with larger matrices and higher minor sizes.

use polynomial_tools::{check_total_positivity, check_tnn_neville, check_tnn_neville_bigint};
use num_bigint::BigInt;
use std::time::Instant;

type Poly = Vec<i128>;

fn poly_zero() -> Poly { vec![0] }
fn poly_is_zero(p: &Poly) -> bool { p.iter().all(|&c| c == 0) }
fn poly_add(a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let mut r = vec![0i128; n];
    for (i, &c) in a.iter().enumerate() { r[i] += c; }
    for (i, &c) in b.iter().enumerate() { r[i] += c; }
    while r.len() > 1 && *r.last().unwrap() == 0 { r.pop(); }
    r
}
fn poly_mul(a: &Poly, b: &Poly) -> Poly {
    if poly_is_zero(a) || poly_is_zero(b) { return poly_zero(); }
    let n = a.len() + b.len() - 1;
    let mut r = vec![0i128; n];
    for (i, &ca) in a.iter().enumerate() {
        if ca == 0 { continue; }
        for (j, &cb) in b.iter().enumerate() { r[i + j] += ca * cb; }
    }
    while r.len() > 1 && *r.last().unwrap() == 0 { r.pop(); }
    r
}

fn compute_pn(mu: &[usize], d: usize, max_n: usize) -> Vec<Poly> {
    let ell = mu.len() - 1;
    let base = d + 1;
    let ns = base.pow(ell as u32);
    let decode = |mut idx: usize| -> Vec<usize> {
        let mut s = vec![0usize; ell];
        for i in (0..ell).rev() { s[i] = idx % base; idx /= base; }
        s
    };
    let encode = |state: &[usize]| -> usize {
        let mut idx = 0;
        for &s in state { idx = idx * base + s; }
        idx
    };
    let mut t_mat = vec![vec![poly_zero(); ns]; ns];
    for j in 0..ns {
        let old = decode(j);
        for v in 0..=d {
            let mut ok = true;
            if v > 0 {
                for m in 1..=ell {
                    let prev = old[m - 1];
                    if prev > 0 && v < prev + mu[m] { ok = false; break; }
                }
            }
            if ok {
                let mut new_s = vec![0usize; ell];
                new_s[0] = v;
                for k in 1..ell { new_s[k] = old[k - 1]; }
                let i = encode(&new_s);
                let w: Poly = if v > 0 { vec![0, 1] } else { vec![1] };
                t_mat[i][j] = poly_add(&t_mat[i][j], &w);
            }
        }
    }
    let mut results = Vec::with_capacity(max_n + 1);
    let mut vec_s: Vec<Poly> = vec![poly_zero(); ns];
    vec_s[0] = vec![1];
    for n in 0..=max_n {
        results.push(vec_s[0].clone());
        if n < max_n {
            let mut new_vec: Vec<Poly> = vec![poly_zero(); ns];
            for i in 0..ns {
                for j in 0..ns {
                    if poly_is_zero(&t_mat[i][j]) || poly_is_zero(&vec_s[j]) { continue; }
                    let term = poly_mul(&t_mat[i][j], &vec_s[j]);
                    new_vec[i] = poly_add(&new_vec[i], &term);
                }
            }
            vec_s = new_vec;
        }
    }
    results
}

fn coeff_matrix_i64(pn: &[Poly]) -> Vec<Vec<i64>> {
    let max_j = pn.iter().map(|p| p.len()).max().unwrap_or(1);
    pn.iter().map(|p| {
        let mut row = vec![0i64; max_j];
        for (j, &c) in p.iter().enumerate() { row[j] = c as i64; }
        row
    }).collect()
}

fn main() {
    for h in 2..=5 {
        println!("=== h = {} ===", h);
        let mu = vec![h, h, 1, 1];
        let max_n = 40;
        let pn = compute_pn(&mu, h, max_n);
        let mat = coeff_matrix_i64(&pn);

        // Neville check (full TNN, all minors)
        let t0 = Instant::now();
        let neville_result = check_tnn_neville(&mat);
        let neville_time = t0.elapsed();
        match &neville_result {
            Ok(()) => println!("  Neville TNN: PASS ({:.2?})", neville_time),
            Err(msg) => println!("  Neville TNN: FAIL ({:.2?}) — {}", neville_time, msg),
        }

        // Brute force for comparison (only up to 4x4 minors, smaller matrix)
        let small_mat: Vec<Vec<i64>> = mat[..20.min(mat.len())].to_vec();
        let t1 = Instant::now();
        let brute_result = check_total_positivity(&small_mat, 4, false);
        let brute_time = t1.elapsed();
        match &brute_result {
            Ok(()) => println!("  Brute 4x4 (n<=19): PASS ({:.2?})", brute_time),
            Err(msg) => println!("  Brute 4x4 (n<=19): FAIL ({:.2?}) — {}", brute_time, msg),
        }

        // Now try larger matrix with Neville (BigInt to avoid overflow)
        let max_n_large = 100;
        let pn_large = compute_pn(&mu, h, max_n_large);
        let mat_big: Vec<Vec<BigInt>> = pn_large.iter().map(|p| {
            let max_j = pn_large.iter().map(|q| q.len()).max().unwrap_or(1);
            let mut row = vec![BigInt::from(0); max_j];
            for (j, &c) in p.iter().enumerate() { row[j] = BigInt::from(c); }
            row
        }).collect();
        let t2 = Instant::now();
        let neville_large = check_tnn_neville_bigint(&mat_big);
        let neville_large_time = t2.elapsed();
        match &neville_large {
            Ok(()) => println!("  Neville BigInt (n<={}): PASS ({:.2?})", max_n_large, neville_large_time),
            Err(msg) => println!("  Neville BigInt (n<={}): FAIL ({:.2?}) — {}", max_n_large, neville_large_time, msg),
        }

        println!();
    }

    // Also test M^(1) which should FAIL
    println!("=== Negative control: M^(1) for h=2 (should fail) ===");
    fn binom(n: i64, k: i64) -> i64 {
        if k < 0 || k > n || n < 0 { return 0; }
        let k = k.min(n - k) as usize;
        let n = n as usize;
        let mut r = 1i64;
        for i in 0..k { r = r * (n - i) as i64 / (i + 1) as i64; }
        r
    }
    let max_n = 25;
    let max_j = 8;
    let mut m1 = Vec::new();
    for n in 0..=max_n {
        let mut row = vec![0i64; (max_j + 1) as usize];
        for j in 0..=max_j {
            let mut val = 0i64;
            for r in 0..=50i64 {
                let s = n - 2*j - 2*r;
                if s < 0 { break; }
                val += binom(2*r, j) * binom(j, s);
            }
            row[j as usize] = val;
        }
        m1.push(row);
    }
    match check_tnn_neville(&m1) {
        Ok(()) => println!("  M^(1): PASS (unexpected!)"),
        Err(msg) => println!("  M^(1): FAIL (expected) — {}", msg),
    }
}
