//! Systematic scan: check if the fault-free coefficient matrix [x^n t^j] B_{mu,d}(x,t)
//! is TNN and has real-rooted rows, for all partitions mu up to given bounds.
//!
//! Usage: cargo run --release --bin tiling_B_tnn_scan [max_n] [max_parts] [max_val]

use num_bigint::BigInt;
use polynomial_tools::{check_tnn_neville_bigint, is_real_rooted};
use std::time::Instant;

type Poly = Vec<i128>;

fn poly_zero() -> Poly {
    vec![0]
}
fn poly_is_zero(p: &Poly) -> bool {
    p.iter().all(|&c| c == 0)
}
fn poly_add(a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let mut r = vec![0i128; n];
    for (i, &c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        r[i] += c;
    }
    while r.len() > 1 && *r.last().unwrap() == 0 {
        r.pop();
    }
    r
}
fn poly_sub(a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let mut r = vec![0i128; n];
    for (i, &c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        r[i] -= c;
    }
    while r.len() > 1 && *r.last().unwrap() == 0 {
        r.pop();
    }
    r
}
fn poly_mul(a: &Poly, b: &Poly) -> Poly {
    if poly_is_zero(a) || poly_is_zero(b) {
        return poly_zero();
    }
    let n = a.len() + b.len() - 1;
    let mut r = vec![0i128; n];
    for (i, &ca) in a.iter().enumerate() {
        if ca == 0 {
            continue;
        }
        for (j, &cb) in b.iter().enumerate() {
            r[i + j] += ca * cb;
        }
    }
    while r.len() > 1 && *r.last().unwrap() == 0 {
        r.pop();
    }
    r
}

fn compute_pn(mu: &[usize], d: usize, max_n: usize) -> Vec<Poly> {
    let ell = mu.len() - 1;
    let base = d + 1;
    let ns = base.pow(ell as u32);
    let decode = |mut idx: usize| -> Vec<usize> {
        let mut s = vec![0usize; ell];
        for i in (0..ell).rev() {
            s[i] = idx % base;
            idx /= base;
        }
        s
    };
    let encode = |state: &[usize]| -> usize {
        let mut idx = 0;
        for &s in state {
            idx = idx * base + s;
        }
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
                    if prev > 0 && v < prev + mu[m] {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                let mut new_s = vec![0usize; ell];
                new_s[0] = v;
                for k in 1..ell {
                    new_s[k] = old[k - 1];
                }
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
                    if poly_is_zero(&t_mat[i][j]) || poly_is_zero(&vec_s[j]) {
                        continue;
                    }
                    let term = poly_mul(&t_mat[i][j], &vec_s[j]);
                    new_vec[i] = poly_add(&new_vec[i], &term);
                }
            }
            vec_s = new_vec;
        }
    }
    results
}

fn extract_b(pn: &[Poly]) -> Vec<Poly> {
    let max_n = pn.len() - 1;
    let mut bn: Vec<Poly> = vec![poly_zero(); max_n + 1];
    for n in 1..=max_n {
        let mut sum = poly_zero();
        for k in 1..n {
            if poly_is_zero(&bn[n - k]) {
                continue;
            }
            let term = poly_mul(&pn[k], &bn[n - k]);
            sum = poly_add(&sum, &term);
        }
        bn[n] = poly_sub(&pn[n], &sum);
    }
    bn
}

fn to_bigint_matrix(bn: &[Poly]) -> Vec<Vec<BigInt>> {
    let max_j = bn.iter().map(|p| p.len()).max().unwrap_or(1);
    bn.iter()
        .map(|p| {
            let mut row = vec![BigInt::from(0); max_j];
            for (j, &c) in p.iter().enumerate() {
                row[j] = BigInt::from(c);
            }
            row
        })
        .collect()
}

/// Check if every row of bn is real-rooted (as a polynomial in t).
/// Returns (all_rr, first_failing_n)
fn check_rows_rr(bn: &[Poly]) -> (bool, Option<(usize, Vec<i64>)>) {
    for (n, p) in bn.iter().enumerate() {
        if poly_is_zero(p) {
            continue;
        }
        let coeffs64: Vec<i64> = p.iter().map(|&c| c as i64).collect();
        if !is_real_rooted(&coeffs64) {
            return (false, Some((n, coeffs64)));
        }
    }
    (true, None)
}

/// Generate all partitions (weakly decreasing) with exactly `k` parts,
/// each part between 1 and `max_val`.
fn partitions(k: usize, max_val: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    fn helper(
        k: usize,
        max_val: usize,
        min_val: usize,
        current: &mut Vec<usize>,
        result: &mut Vec<Vec<usize>>,
    ) {
        if k == 0 {
            result.push(current.clone());
            return;
        }
        // We build in decreasing order, so next part <= current last (or max_val if empty)
        let upper = if current.is_empty() {
            max_val
        } else {
            *current.last().unwrap()
        };
        for v in (min_val..=upper).rev() {
            current.push(v);
            helper(k - 1, max_val, min_val, current, result);
            current.pop();
        }
    }
    let mut current = Vec::new();
    helper(k, max_val, 1, &mut current, &mut result);
    result
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let max_n: usize = args.get(0).and_then(|s| s.parse().ok()).unwrap_or(40);
    let max_ell_plus_1: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(5);
    let max_val: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(6);
    let max_d: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(7);

    println!(
        "Scanning: max_n={}, max_parts={}, max_val={}, max_d={}",
        max_n, max_ell_plus_1, max_val, max_d
    );
    println!(
        "{:<22} {:>2} {:>6} {:>5} {:>5}  {}",
        "mu", "d", "mu1>=d", "TNN", "RR", "notes"
    );
    println!("{}", "-".repeat(80));

    let mut total = 0;
    let mut tnn_pass = 0;
    let mut rr_pass = 0;
    let mut tnn_fail_list = Vec::new();
    let mut rr_fail_list = Vec::new();

    for num_parts in 2..=max_ell_plus_1 {
        for mu in partitions(num_parts, max_val) {
            // mu is weakly decreasing, at least 2 parts
            let ell = mu.len() - 1;
            for d in 1..=max_d.min(mu[0]) {
                // Skip if state space too large
                let ns = (d + 1).pow(ell as u32);
                if ns > 50_000 {
                    continue;
                }

                total += 1;
                let t0 = Instant::now();
                let pn = compute_pn(&mu, d, max_n);
                let bn = extract_b(&pn);

                // Check max_j
                let mj = bn.iter().map(|p| p.len()).max().unwrap_or(1);
                if mj <= 1 {
                    // Trivially TNN and RR
                    tnn_pass += 1;
                    rr_pass += 1;
                    continue;
                }

                let mat = to_bigint_matrix(&bn);
                let tnn_ok = check_tnn_neville_bigint(&mat).is_ok();
                let (rr_ok, rr_fail_info) = check_rows_rr(&bn);
                let elapsed = t0.elapsed();

                if tnn_ok {
                    tnn_pass += 1;
                }
                if rr_ok {
                    rr_pass += 1;
                }

                let cond = if mu[1] >= d { "yes" } else { "no" };

                if !tnn_ok || !rr_ok {
                    let mut notes = String::new();
                    if !tnn_ok {
                        notes.push_str("TNN FAIL ");
                    }
                    if !rr_ok {
                        if let Some((n, _coeffs)) = &rr_fail_info {
                            notes.push_str(&format!("RR FAIL at n={}", n));
                        }
                    }
                    println!(
                        "{:<22} {:>2} {:>6} {:>5} {:>5}  {} ({:.2?})",
                        format!("{:?}", mu),
                        d,
                        cond,
                        if tnn_ok { "PASS" } else { "FAIL" },
                        if rr_ok { "PASS" } else { "FAIL" },
                        notes,
                        elapsed
                    );
                    if !tnn_ok {
                        tnn_fail_list.push((mu.clone(), d));
                    }
                    if !rr_ok {
                        rr_fail_list.push((mu.clone(), d, rr_fail_info));
                    }
                }
            }
        }
    }

    println!("\n{}", "=".repeat(80));
    println!("Total cases tested: {}", total);
    println!("TNN: {} pass, {} fail", tnn_pass, tnn_fail_list.len());
    println!("RR:  {} pass, {} fail", rr_pass, rr_fail_list.len());

    if !tnn_fail_list.is_empty() {
        println!("\nTNN failures:");
        for (mu, d) in &tnn_fail_list {
            println!("  mu={:?}, d={}", mu, d);
        }
    }
    if !rr_fail_list.is_empty() {
        println!(
            "\nRR failures (first {} shown):",
            rr_fail_list.len().min(20)
        );
        for (mu, d, info) in rr_fail_list.iter().take(20) {
            if let Some((n, coeffs)) = info {
                println!("  mu={:?}, d={}: B_{} = {:?}", mu, d, n, coeffs);
            }
        }
    }
}
