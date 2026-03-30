//! 1) Generate tiling polynomials via transfer matrix
//! 2) Search for a shorter recurrence (over Q(t))
//! 3) Check interlacing depth: for which k does P_n << P_{n+k}?
//!
//! Usage: cargo run --release --bin tiling_short_rec -- 3,2,2,1 4

use polynomial_tools::{real_roots, is_real_rooted};
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

fn format_poly_t(p: &Poly) -> String {
    if poly_is_zero(p) { return "0".into(); }
    let mut terms = vec![];
    for (i, &c) in p.iter().enumerate() {
        if c == 0 { continue; }
        let t_part = match i {
            0 => format!("{}", c),
            1 if c == 1 => "t".into(),
            1 if c == -1 => "-t".into(),
            1 => format!("{}t", c),
            _ if c == 1 => format!("t^{}", i),
            _ if c == -1 => format!("-t^{}", i),
            _ => format!("{}t^{}", c, i),
        };
        terms.push(t_part);
    }
    if terms.is_empty() { "0".into() } else { terms.join(" + ") }
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

/// Try to find a recurrence P_n = sum_{j=1}^{order} alpha_j(t) * P_{n-j}
/// by solving a linear system over Q(t).
/// We work coefficient-by-coefficient in t: for each t^k, find rational alpha_j coefficients.
fn search_short_recurrence(polys: &[Poly], max_order: usize) {
    println!("\n=== Searching for shorter recurrence ===");

    let n_polys = polys.len();
    // Find max degree in t across all polys
    let max_t_deg = polys.iter().map(|p| deg(p)).max().unwrap_or(0);

    for order in 1..=max_order {
        // Try: P_n = alpha_1(t)*P_{n-1} + ... + alpha_{order}(t)*P_{n-order}
        // Each alpha_j(t) = sum_k alpha_{j,k} * t^k
        // We need enough equations. For each n and each t^k coefficient, we get one equation.
        // Unknowns: alpha_{j,k} for j=1..order, k=0..max_alpha_deg
        // We allow alpha to have degree up to max_t_deg + 2

        let max_alpha_deg = max_t_deg + 2;
        let n_unknowns = order * (max_alpha_deg + 1);

        // Collect equations: for n = order..n_polys-1, coefficient of t^k in P_n = sum_j alpha_j * P_{n-j}
        // Each equation: sum_{j=1}^{order} sum_{m=0}^{max_alpha_deg} alpha_{j,m} * [t^{k-m}] P_{n-j} = [t^k] P_n

        let mut equations: Vec<(Vec<i128>, i128)> = vec![];

        for n in order..n_polys.min(order + n_unknowns + 5) {
            let p_n = &polys[n];
            let max_k = deg(p_n) + max_alpha_deg;
            for k in 0..=max_k {
                let rhs = if k < p_n.len() { p_n[k] } else { 0 };
                let mut row = vec![0i128; n_unknowns];
                for j in 1..=order {
                    let p_prev = &polys[n - j];
                    for m in 0..=max_alpha_deg {
                        if k >= m {
                            let coeff_idx = (j - 1) * (max_alpha_deg + 1) + m;
                            let t_power = k - m;
                            let val = if t_power < p_prev.len() { p_prev[t_power] } else { 0 };
                            row[coeff_idx] = val;
                        }
                    }
                }
                equations.push((row, rhs));
            }
        }

        // Solve by Gaussian elimination over Q (using i128 for now, clearing denominators)
        // Simple approach: find alpha values that satisfy ALL equations
        // Use exact rational arithmetic with lcm/gcd

        // Actually, let's just try a simpler approach:
        // For each t-degree independently, check if there's a recurrence with CONSTANT alpha_j
        // i.e., alpha_j(t) = c_j (no t-dependence)

        // Even simpler: check if there's a recurrence P_n = c1*P_{n-1} + c2(t)*P_{n-2} + ...
        // where c1 = 1 (from the empty column) and we look for the minimal additional terms

        // Let's just check: does P_n - P_{n-1} satisfy a shorter recurrence?
        let mut diffs: Vec<Poly> = vec![];
        for n in 1..n_polys {
            let mut d = polys[n].clone();
            for (i, &c) in polys[n-1].iter().enumerate() {
                if i < d.len() { d[i] -= c; }
            }
            while d.len() > 1 && *d.last().unwrap() == 0 { d.pop(); }
            diffs.push(d);
        }

        // Check if diffs satisfy a recurrence of given order
        // diff[n] = P_{n+1} - P_n. Does diff[n] = sum beta_j * diff[n-j]?
        // This would give P_{n+1} - P_n = sum beta_j (P_{n+1-j} - P_{n-j})

        // Actually let me just try to verify small recurrences directly
        if order <= 6 {
            // Extract alpha_j(t) from the first few equations
            // P_n = sum alpha_j P_{n-j} for j=1..order
            // Use n = order, order+1, ..., order + enough to determine alphas

            // Approach: for each coefficient of t^k, set up a linear system
            // But this is complex. Let me just try specific forms.
            continue;
        }
    }

    // Direct approach: check if there's a 4-term recurrence
    // P_n = a(t) P_{n-1} + b(t) P_{n-2} + c(t) P_{n-3} + d(t) P_{n-4}
    println!("\nChecking 4-term recurrence P_n = a(t)P_{{n-1}} + b(t)P_{{n-2}} + c(t)P_{{n-3}} + d(t)P_{{n-4}}:");
    // For each t^k coeff in P_n: sum of a_m * [t^{k-m}]P_{n-1} + ... = [t^k]P_n
    // Try: a(t) = 1, b(t) = ?, c(t) = ?, d(t) = ? with b,c,d polynomials in t

    // Actually, the simplest check: compute R_n = P_n - P_{n-1} and see if R_n
    // has a shorter recurrence.
    println!("\nR_n = P_n - P_{{n-1}} (first differences):");
    let mut diffs: Vec<Poly> = vec![];
    for n in 1..n_polys {
        let mut d = vec![0i128; polys[n].len().max(polys[n-1].len())];
        for (i, &c) in polys[n].iter().enumerate() { d[i] += c; }
        for (i, &c) in polys[n-1].iter().enumerate() { d[i] -= c; }
        while d.len() > 1 && *d.last().unwrap() == 0 { d.pop(); }
        diffs.push(d);
    }

    for (i, d) in diffs.iter().enumerate().take(20) {
        if deg(d) <= 5 {
            println!("  R_{} = {}", i+1, format_poly_t(d));
        }
    }

    // Now check: can we express R_n in terms of fewer past R's?
    // Try R_n = sum_{j=1}^{order} beta_j(t) R_{n-j}
    println!("\nSearching for recurrence on R_n = P_n - P_{{n-1}}:");
    for order in 1..=8 {
        // For each candidate order, extract beta_j by solving
        // from the first few R values
        let start = order;
        let n_eqs = diffs.len() - start;
        if n_eqs < 3 { break; }

        // Very simple: try to find the B(x,t) for the R_n sequence
        // R_n = sum beta_j R_{n-j}
        let mut betas: Vec<Poly> = vec![poly_zero(); order + 1];
        let mut ok = true;

        for n in start..diffs.len() {
            let mut residual = diffs[n].clone();
            for j in 1..=order {
                if j <= betas.len() - 1 && !poly_is_zero(&betas[j]) {
                    let term = poly_mul(&betas[j], &diffs[n - j]);
                    for (i, &c) in term.iter().enumerate() {
                        if i < residual.len() { residual[i] -= c; }
                        else {
                            residual.resize(i + 1, 0);
                            residual[i] -= c;
                        }
                    }
                }
            }
            while residual.len() > 1 && *residual.last().unwrap() == 0 { residual.pop(); }

            if poly_is_zero(&residual) { continue; }

            // residual should be expressible as beta_? * R_{n-?} for some unfixed beta
            // Find the first unfixed beta
            let mut found = false;
            for j in 1..=order {
                if poly_is_zero(&betas[j]) && !poly_is_zero(&diffs[n - j]) {
                    // Try beta_j = residual / R_{n-j} (polynomial division)
                    // This only works if it divides exactly
                    let divisor = &diffs[n - j];
                    // Simple case: if R_{n-j} starts with a constant
                    // For now, just record and verify later
                    // Actually this approach is too simplistic. Let me just
                    // print the B-polynomial for the R sequence.
                    found = true;
                    break;
                }
            }
            if !found {
                ok = false;
                break;
            }
            break; // need more sophisticated solver
        }
    }

    // Better approach: just compute B for the R_n sequence directly
    // like we did for P_n
    println!("\nExtracting B-coefficients for R_n sequence:");
    let mut r_betas: Vec<Poly> = vec![poly_zero(); diffs.len() + 1];
    for n in 1..diffs.len() {
        let mut rhs = diffs[n].clone();
        for j in 1..n {
            if poly_is_zero(&r_betas[j]) { continue; }
            let term = poly_mul(&r_betas[j], &diffs[n - j]);
            for (i, &c) in term.iter().enumerate() {
                if i >= rhs.len() { rhs.resize(i + 1, 0); }
                rhs[i] -= c;
            }
        }
        while rhs.len() > 1 && *rhs.last().unwrap() == 0 { rhs.pop(); }
        if !poly_is_zero(&rhs) {
            r_betas[n] = rhs.clone();
            println!("  beta_{} = {}", n, format_poly_t(&rhs));
        }
    }
}

fn check_interlacing_depth(polys: &[Poly], max_k: usize) {
    println!("\n=== Interlacing depth: P_n << P_{{n+k}} ===");
    println!("{:>5} ", "k\\n");

    let n_max = polys.len() - 1;

    for k in 1..=max_k {
        let mut first_fail: Option<usize> = None;
        let mut last_checked = 0;

        for n in 0..=(n_max - k) {
            let f = match poly_to_i64(&polys[n]) { Some(v) => v, None => break };
            let g = match poly_to_i64(&polys[n + k]) { Some(v) => v, None => break };

            let df = { let mut d = f.len()-1; while d>0 && f[d]==0 { d-=1; } d };
            let dg = { let mut d = g.len()-1; while d>0 && g[d]==0 { d-=1; } d };

            // Skip trivial cases
            if df == 0 || dg == 0 { last_checked = n; continue; }
            // Only check when degrees are compatible
            if dg != df && dg != df + 1 { last_checked = n; continue; }

            match check_interleaves(&f, &g) {
                Some(true) => { last_checked = n; }
                Some(false) => {
                    if first_fail.is_none() { first_fail = Some(n); }
                    break;
                }
                None => break,
            }
        }

        match first_fail {
            None => println!("  k={}: P_n << P_{{n+{}}} holds for all n <= {} (where degrees are compatible)",
                           k, k, last_checked),
            Some(n) => println!("  k={}: FAILS at n={} (deg P_{}={}, deg P_{}={})",
                           k, n, n, deg(&polys[n]), n+k, deg(&polys[n+k])),
        }
    }
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mu: Vec<usize> = args[0].split(',').map(|s| s.trim().parse().unwrap()).collect();
    let d: usize = args[1].parse().unwrap();
    let max_n: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(40);

    println!("mu = {:?}, d = {}", mu, d);
    let polys = compute_polys(&mu, d, max_n);

    // Print first few
    for (n, p) in polys.iter().enumerate().take(18) {
        if deg(p) <= 6 {
            println!("P_{} = {}", n, format_poly_t(p));
        }
    }

    search_short_recurrence(&polys, 8);
    check_interlacing_depth(&polys, 10);
}
