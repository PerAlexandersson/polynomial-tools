//! Find the minimal recurrence for tiling polynomials over Q(t).
//! i.e., find a_0(t)*P_n + a_1(t)*P_{n-1} + ... + a_k(t)*P_{n-k} = 0
//! with a_j(t) polynomials in t, minimal k.
//!
//! Also check interlacing depth.
//!
//! Usage: cargo run --release --bin tiling_minrec -- 3,2,2,1 4

use polynomial_tools::real_roots;
use std::env;

type Poly = Vec<i128>;

fn poly_zero() -> Poly { vec![0] }
fn poly_one() -> Poly { vec![1] }
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

fn poly_to_i64(p: &Poly) -> Option<Vec<i64>> {
    p.iter().map(|&c| i64::try_from(c).ok()).collect()
}

fn deg(p: &Poly) -> usize {
    let mut d = p.len() - 1;
    while d > 0 && p[d] == 0 { d -= 1; }
    d
}

fn format_tpoly(p: &[i128]) -> String {
    let mut terms = vec![];
    for (i, &c) in p.iter().enumerate() {
        if c == 0 { continue; }
        let s = match i {
            0 => format!("{}", c),
            1 if c == 1 => "t".into(),
            1 if c == -1 => "-t".into(),
            1 => format!("{}t", c),
            _ if c == 1 => format!("t^{}", i),
            _ if c == -1 => format!("-t^{}", i),
            _ => format!("{}t^{}", c, i),
        };
        terms.push(s);
    }
    if terms.is_empty() { "0".into() } else { terms.join(" + ") }
}

fn encode(state: &[usize], d: usize) -> usize {
    let base = d + 1;
    let mut idx = 0;
    for &s in state { idx = idx * base + s; }
    idx
}

fn decode(mut idx: usize, ell: usize, d: usize) -> Vec<usize> {
    let base = d + 1;
    let mut s = vec![0usize; ell];
    for i in (0..ell).rev() { s[i] = idx % base; idx /= base; }
    s
}

fn compute_polys(mu: &[usize], d: usize, max_n: usize) -> Vec<Poly> {
    let ell = mu.len() - 1;
    let base = d + 1;
    let ns = base.pow(ell as u32);
    let mut t_mat = vec![vec![poly_zero(); ns]; ns];
    for j in 0..ns {
        let old = decode(j, ell, d);
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
                let i = encode(&new_s, d);
                let w: Poly = if v > 0 { vec![0, 1] } else { vec![1] };
                t_mat[i][j] = poly_add(&t_mat[i][j], &w);
            }
        }
    }
    let mut vec_s: Vec<Poly> = vec![poly_zero(); ns];
    vec_s[0] = poly_one();
    let mut polys = vec![];
    for _n in 0..=max_n {
        polys.push(vec_s[0].clone());
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
    polys
}

/// Search for recurrence: a_0(t)*P_n + ... + a_k(t)*P_{n-k} = 0
/// where a_j(t) are polynomials of degree <= max_deg_a.
/// Returns the coefficients a_0, ..., a_k if found.
fn find_recurrence(polys: &[Poly], order: usize, max_deg_a: usize) -> Option<Vec<Poly>> {
    // Unknowns: a_{j,m} for j=0..order, m=0..max_deg_a
    // Total unknowns: (order+1)*(max_deg_a+1)
    let n_unknowns = (order + 1) * (max_deg_a + 1);

    // Equations: for each n >= order, for each t^p coefficient:
    // sum_{j=0}^{order} sum_{m=0}^{max_deg_a} a_{j,m} * [t^{p-m}] P_{n-j} = 0
    let mut rows: Vec<Vec<i128>> = vec![];
    let n_start = order;
    let n_end = polys.len().min(n_start + n_unknowns + 10);

    let max_t_deg = polys.iter().take(n_end).map(|p| deg(p)).max().unwrap_or(0);

    for n in n_start..n_end {
        for p in 0..=(max_t_deg + max_deg_a) {
            let mut row = vec![0i128; n_unknowns];
            for j in 0..=order {
                let pnj = &polys[n - j];
                for m in 0..=max_deg_a {
                    if p >= m {
                        let t_idx = p - m;
                        let coeff = if t_idx < pnj.len() { pnj[t_idx] } else { 0 };
                        let var_idx = j * (max_deg_a + 1) + m;
                        row[var_idx] = coeff;
                    }
                }
            }
            // Skip all-zero rows
            if row.iter().all(|&c| c == 0) { continue; }
            rows.push(row);
        }
    }

    // Gaussian elimination over Q to find the null space
    // We need a nontrivial solution with a_0 != 0 (monic in some sense)
    let n_rows = rows.len();
    let n_cols = n_unknowns;

    // Use fraction-free Gaussian elimination (integers only)
    // Actually, let's do row reduction with rational arithmetic
    // For simplicity, use f64 and check if the rank is < n_unknowns

    // Convert to f64
    let mut mat: Vec<Vec<f64>> = rows.iter()
        .map(|r| r.iter().map(|&c| c as f64).collect())
        .collect();

    // Row echelon form
    let mut pivot_col = 0;
    let mut pivot_row = 0;
    let mut pivots = vec![];

    while pivot_row < n_rows && pivot_col < n_cols {
        // Find pivot
        let mut max_val = 0.0f64;
        let mut max_row = pivot_row;
        for i in pivot_row..n_rows {
            if mat[i][pivot_col].abs() > max_val {
                max_val = mat[i][pivot_col].abs();
                max_row = i;
            }
        }
        if max_val < 1e-10 {
            pivot_col += 1;
            continue;
        }
        mat.swap(pivot_row, max_row);
        let scale = mat[pivot_row][pivot_col];
        for c in 0..n_cols { mat[pivot_row][c] /= scale; }
        for i in 0..n_rows {
            if i == pivot_row { continue; }
            let factor = mat[i][pivot_col];
            if factor.abs() < 1e-15 { continue; }
            for c in 0..n_cols {
                mat[i][c] -= factor * mat[pivot_row][c];
            }
        }
        pivots.push(pivot_col);
        pivot_row += 1;
        pivot_col += 1;
    }

    let rank = pivots.len();
    let null_dim = n_cols - rank;

    if null_dim == 0 {
        return None; // No recurrence of this order
    }

    // Find a null vector: choose the first free variable
    let free_cols: Vec<usize> = (0..n_cols).filter(|c| !pivots.contains(c)).collect();
    let free_col = free_cols[0];

    // Back-substitute: set free_col = 1, others = 0 among free vars
    let mut sol = vec![0.0f64; n_cols];
    sol[free_col] = 1.0;
    for (&pc, pr) in pivots.iter().zip(0..rank).rev() {
        let mut val = 0.0;
        for c in (pc+1)..n_cols {
            val += mat[pr][c] * sol[c];
        }
        sol[pc] = -val;
    }

    // Convert back to rational: find LCM of denominators
    // Round to nearest rational with small denominator
    // For now, just round and verify

    // Scale to make integer-ish
    let max_abs = sol.iter().map(|&x| x.abs()).fold(0.0f64, f64::max);
    if max_abs < 1e-10 { return None; }

    // Try to find integer ratios
    // Find the entry with smallest nonzero absolute value
    let min_nonzero = sol.iter().map(|&x| x.abs()).filter(|&x| x > 1e-10).fold(f64::MAX, f64::min);
    let scaled: Vec<f64> = sol.iter().map(|&x| x / min_nonzero).collect();

    // Round to integers
    let int_sol: Vec<i128> = scaled.iter().map(|&x| x.round() as i128).collect();

    // Extract a_j(t)
    let mut coeffs = vec![];
    for j in 0..=order {
        let mut aj = vec![0i128; max_deg_a + 1];
        for m in 0..=max_deg_a {
            aj[m] = int_sol[j * (max_deg_a + 1) + m];
        }
        while aj.len() > 1 && *aj.last().unwrap() == 0 { aj.pop(); }
        coeffs.push(aj);
    }

    // Verify: check that sum a_j(t) * P_{n-j} = 0 for all n
    let mut verified = true;
    for n in order..polys.len() {
        let mut sum = poly_zero();
        for j in 0..=order {
            let term = poly_mul(&coeffs[j], &polys[n - j]);
            sum = poly_add(&sum, &term);
        }
        if !poly_is_zero(&sum) {
            verified = false;
            break;
        }
    }

    if verified {
        Some(coeffs)
    } else {
        None
    }
}

fn check_interleaves(f: &[i64], g: &[i64]) -> Option<bool> {
    let df = { let mut d = f.len()-1; while d>0 && f[d]==0 { d-=1; } d };
    let dg = { let mut d = g.len()-1; while d>0 && g[d]==0 { d-=1; } d };
    if df == 0 { return Some(true); }
    if dg == 0 { return Some(true); }
    if df + 1 != dg && df != dg { return Some(false); }
    let rf = real_roots(f)?;
    let rg = real_roots(g)?;
    if rf.len() != df || rg.len() != dg { return Some(false); }
    let mut rf = rf; let mut rg = rg;
    rf.sort(); rg.sort();
    if df + 1 == dg {
        for i in 0..df {
            if rf[i] < rg[i] || rf[i] > rg[i+1] { return Some(false); }
        }
        Some(true)
    } else {
        for i in 0..df {
            if rf[i] > rg[i] { return Some(false); }
            if i+1 < df && rg[i] > rf[i+1] { return Some(false); }
        }
        Some(true)
    }
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mu: Vec<usize> = args[0].split(',').map(|s| s.trim().parse().unwrap()).collect();
    let d: usize = args[1].parse().unwrap();
    let max_n: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(40);

    println!("mu = {:?}, d = {}", mu, d);
    let polys = compute_polys(&mu, d, max_n);

    // Part 1: Search for minimal recurrence
    println!("\n=== Minimal recurrence search ===");
    for order in 1..=13 {
        for max_deg in 0..=4 {
            if let Some(coeffs) = find_recurrence(&polys, order, max_deg) {
                println!("Found recurrence of order {} with max t-degree {}:", order, max_deg);
                for (j, c) in coeffs.iter().enumerate() {
                    if !poly_is_zero(c) {
                        println!("  a_{}(t) = {}", j, format_tpoly(c));
                    }
                }
                println!("  (verified against all {} polynomials)", polys.len());

                // Print in readable form
                print!("  ");
                let mut first = true;
                for (j, c) in coeffs.iter().enumerate() {
                    if poly_is_zero(c) { continue; }
                    if !first { print!(" + "); }
                    print!("({})*P_{{n-{}}}", format_tpoly(c), j);
                    first = false;
                }
                println!(" = 0");

                println!("\nMinimal order: {}", order);
                return;
            }
        }
    }
    println!("No short recurrence found up to order 13");
    print_interlacing_depth(&polys);
}

fn print_interlacing_depth(polys: &[Poly]) {
    println!("\n=== Interlacing depth ===");
    for k in 1..=12 {
        let mut first_fail: Option<(usize, usize, usize)> = None;
        let mut last_ok = 0;
        let mut n_checked = 0;

        for n in 0..polys.len() - k {
            let f = match poly_to_i64(&polys[n]) { Some(v) => v, None => break };
            let g = match poly_to_i64(&polys[n + k]) { Some(v) => v, None => break };
            let df = deg(&polys[n]);
            let dg = deg(&polys[n + k]);
            if df == 0 || dg == 0 { last_ok = n; n_checked += 1; continue; }
            if dg != df && dg != df + 1 { last_ok = n; n_checked += 1; continue; }

            match check_interleaves(&f, &g) {
                Some(true) => { last_ok = n; n_checked += 1; }
                Some(false) => {
                    first_fail = Some((n, df, dg));
                    break;
                }
                None => break,
            }
        }
        match first_fail {
            None => println!("  k={:>2}: P_n << P_{{n+{}}} holds for all {} compatible pairs up to n={}",
                           k, k, n_checked, last_ok),
            Some((n, df, dg)) => println!("  k={:>2}: FAILS at n={} (deg {}→{})",
                           k, n, df, dg),
        }
    }
}
