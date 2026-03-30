//! Exact interlacing test for the peak polynomial recursion.
//!
//! Tests Q ≼ tR where Q = mP - R and P+ = Q + tR.
//! Uses exact arithmetic: Sturm chains for root counting, polynomial GCD
//! to remove common factors before checking interlacing.

use std::collections::BTreeSet;


/// Check if polynomial is real-rooted with all roots in (-inf, 0].
fn is_rr_nonpos(coeffs: &[i64]) -> bool {
    // Check real-rootedness first
    if !polynomial_tools::is_real_rooted(coeffs) {
        return false;
    }
    // Check all roots are non-positive by evaluating at small positive value
    // If all roots ≤ 0 and leading coefficient has consistent sign, f(ε) has sign = lc * ε^d
    // More robust: use the Sturm-based root finder and check all roots ≤ 0
    match polynomial_tools::real_roots(coeffs) {
        None => false,
        Some(roots) => {
            use num_traits::Zero;
            roots.iter().all(|r| *r <= num_rational::Ratio::zero())
        }
    }
}

/// Exact interlacing check: f ≼ g.
/// Removes common factors via polynomial GCD first,
/// then checks root interlacing using Sturm's method.
fn exact_interlaces(fc: &[i64], gc: &[i64]) -> Result<bool, String> {
    match polynomial_tools::check_weak_interlacing(fc, gc) {
        Some(result) => Ok(result),
        None => Err("degree difference not 1 after GCD removal".to_string()),
    }
}

/// Compute GCD of two polynomials over BigRational.

// ── Board / Bruhat / Peaks ──────────────────────────────────────────

fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> {
    let n = perm.len();
    let mut visited: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut queue: BTreeSet<Vec<u8>> = BTreeSet::new();
    queue.insert(perm.to_vec());
    while let Some(current) = queue.pop_last() {
        for i in 0..n {
            for j in i + 1..n {
                if current[i] > current[j] {
                    let mut child = current.clone();
                    child.swap(i, j);
                    if !visited.contains(&child) {
                        queue.insert(child);
                    }
                }
            }
        }
        visited.insert(current);
    }
    visited.into_iter().collect()
}

fn board_to_312_perm(board: &[u8]) -> Vec<u8> {
    let n = board.len();
    let mut perm = vec![0u8; n];
    let mut used = vec![false; n + 1];
    for i in 0..n {
        let max_col = (board[i] as usize).min(n);
        for c in (1..=max_col).rev() {
            if !used[c] {
                perm[i] = c as u8;
                used[c] = true;
                break;
            }
        }
    }
    perm
}

fn peaks(w: &[u8]) -> usize {
    if w.len() < 3 {
        return 0;
    }
    (1..w.len() - 1)
        .filter(|&i| w[i - 1] < w[i] && w[i] > w[i + 1])
        .count()
}

fn poly_trim(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
    v
}
fn poly_is_zero(p: &[i64]) -> bool {
    p.iter().all(|&c| c == 0)
}
fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] += v;
    }
    poly_trim(&r)
}
fn poly_sub(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] -= v;
    }
    poly_trim(&r)
}
fn poly_scale(p: &[i64], c: i64) -> Vec<i64> {
    poly_trim(&p.iter().map(|&x| x * c).collect::<Vec<_>>())
}
fn poly_mul_t(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; p.len() + 1];
    for (i, &v) in p.iter().enumerate() {
        r[i + 1] = v;
    }
    poly_trim(&r)
}
fn poly_fmt(p: &[i64]) -> String {
    let p = poly_trim(p);
    if poly_is_zero(&p) {
        return "0".into();
    }
    let mut terms = vec![];
    for (i, &c) in p.iter().enumerate() {
        if c == 0 {
            continue;
        }
        match (c, i) {
            (c, 0) => terms.push(format!("{}", c)),
            (1, 1) => terms.push("t".into()),
            (c, 1) => terms.push(format!("{}t", c)),
            (1, e) => terms.push(format!("t^{}", e)),
            (c, e) => terms.push(format!("{}t^{}", c, e)),
        }
    }
    terms.join(" + ")
}

fn generate_boards(n: usize) -> Vec<Vec<u8>> {
    let mut results = vec![];
    let mut current = vec![];
    gen_rec(n, n, 0, &mut current, &mut results);
    results
}
fn gen_rec(
    n: usize,
    max_col: usize,
    depth: usize,
    current: &mut Vec<u8>,
    results: &mut Vec<Vec<u8>>,
) {
    if depth == n {
        results.push(current.clone());
        return;
    }
    let min_val = (depth + 1).max(if depth > 0 { current[depth - 1] as usize } else { 1 });
    for v in min_val..=max_col {
        current.push(v as u8);
        gen_rec(n, max_col, depth + 1, current, results);
        current.pop();
    }
}

fn main() {
    let max_n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(7);
    println!(
        "EXACT interlacing test: (mP-R) ≼ tR for all boards with n <= {}\n",
        max_n
    );

    let mut total = 0;
    let mut q_rr_fails = 0;
    let mut interl_fails = 0;
    let mut pplus_rr_fails = 0;

    for n in 2..=max_n {
        let boards = generate_boards(n);
        for board in &boards {
            total += 1;
            let m = (board[0] as usize).min(n);
            let perm = board_to_312_perm(board);
            let ideal = bruhat_lower_ideal(&perm);

            let mut d_polys = vec![vec![0i64]; m + 1];
            let mut p_poly = vec![0i64];
            for pi in &ideal {
                let pk = peaks(pi);
                while p_poly.len() <= pk {
                    p_poly.push(0);
                }
                p_poly[pk] += 1;
                if pi.len() >= 2 && pi[0] > pi[1] {
                    let k = pi[0] as usize;
                    if k <= m {
                        while d_polys[k].len() <= pk {
                            d_polys[k].push(0);
                        }
                        d_polys[k][pk] += 1;
                    }
                }
            }
            p_poly = poly_trim(&p_poly);

            let mut r_poly = vec![0i64];
            for j in 1..=m {
                r_poly = poly_add(&r_poly, &poly_scale(&d_polys[j], j as i64));
            }

            let q = poly_sub(&poly_scale(&p_poly, m as i64), &r_poly);
            let tr = poly_mul_t(&r_poly);

            if !poly_is_zero(&q) && !is_rr_nonpos(&q) {
                q_rr_fails += 1;
                println!("FAIL: Q not RR for {:?}: Q={}", board, poly_fmt(&q));
            }

            if !poly_is_zero(&q) && !poly_is_zero(&tr) {
                match exact_interlaces(&q, &tr) {
                    Ok(true) => {}
                    Ok(false) => {
                        interl_fails += 1;
                        println!("FAIL: Q ≼ tR for {:?}", board);
                        println!("  Q  = {}", poly_fmt(&q));
                        println!("  tR = {}", poly_fmt(&tr));
                    }
                    Err(e) => {
                        interl_fails += 1;
                        println!("ERR: Q ≼ tR for {:?}: {}", board, e);
                        println!("  Q  = {}", poly_fmt(&q));
                        println!("  tR = {}", poly_fmt(&tr));
                    }
                }
            }

            let pplus = poly_add(&q, &tr);
            if !poly_is_zero(&pplus) && !is_rr_nonpos(&pplus) {
                pplus_rr_fails += 1;
                println!("FAIL: P+ not RR for {:?}: {}", board, poly_fmt(&pplus));
            }
        }
        println!(
            "n={}: total={}, Q_RR={}, Q≼tR={}, P+_RR={}",
            n, total, q_rr_fails, interl_fails, pplus_rr_fails
        );
    }

    println!("\n=== SUMMARY ===");
    println!(
        "Boards: {}, Q not RR: {}, Q≼tR: {}, P+ not RR: {}",
        total, q_rr_fails, interl_fails, pplus_rr_fails
    );
    if q_rr_fails == 0 && interl_fails == 0 && pplus_rr_fails == 0 {
        println!("ALL CONDITIONS HOLD.");
    }
}
