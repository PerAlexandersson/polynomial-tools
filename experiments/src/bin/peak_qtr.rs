//! Test the condition (mP - R) ≼ tR for all Ferrers boards.
//! Where P = total peak poly, R = Σ j·D_j, and P+ = mP + (t-1)R = (mP-R) + tR.
//! If (mP-R) ≼ tR, then P+ is real-rooted by the interlacing cone property.

use std::collections::BTreeSet;

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
                    if !visited.contains(&child) { queue.insert(child); }
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
            if !used[c] { perm[i] = c as u8; used[c] = true; break; }
        }
    }
    perm
}

fn peaks(w: &[u8]) -> usize {
    if w.len() < 3 { return 0; }
    (1..w.len() - 1).filter(|&i| w[i - 1] < w[i] && w[i] > w[i + 1]).count()
}

fn poly_trim(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.len() > 1 && *v.last().unwrap() == 0 { v.pop(); }
    v
}
fn poly_is_zero(p: &[i64]) -> bool { p.iter().all(|&c| c == 0) }
fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &v) in a.iter().enumerate() { r[i] += v; }
    for (i, &v) in b.iter().enumerate() { r[i] += v; }
    poly_trim(&r)
}
fn poly_sub(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &v) in a.iter().enumerate() { r[i] += v; }
    for (i, &v) in b.iter().enumerate() { r[i] -= v; }
    poly_trim(&r)
}
fn poly_scale(p: &[i64], c: i64) -> Vec<i64> {
    poly_trim(&p.iter().map(|&x| x * c).collect::<Vec<_>>())
}
fn poly_mul_t(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; p.len() + 1];
    for (i, &v) in p.iter().enumerate() { r[i + 1] = v; }
    poly_trim(&r)
}
fn poly_fmt(p: &[i64]) -> String {
    let p = poly_trim(p);
    if poly_is_zero(&p) { return "0".into(); }
    let mut terms = vec![];
    for (i, &c) in p.iter().enumerate() {
        if c == 0 { continue; }
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

fn find_real_roots(coeffs: &[i64]) -> Option<Vec<f64>> {
    let roots = polynomial_tools::real_roots(coeffs)?;
    let mut float_roots: Vec<f64> = roots.iter().map(|r| {
        let n = r.numer().to_str_radix(10).parse::<f64>().unwrap();
        let d = r.denom().to_str_radix(10).parse::<f64>().unwrap();
        n / d
    }).collect();
    float_roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Some(float_roots)
}

/// Check f ≼ g (f interlaces g): roots of f are weakly LEFT of roots of g.
fn interlaces(f: &[i64], g: &[i64]) -> bool {
    let f = poly_trim(f);
    let g = poly_trim(g);
    if poly_is_zero(&f) { return true; }
    if poly_is_zero(&g) { return false; }
    let fr = match find_real_roots(&f) { Some(r) => r, None => return false };
    let gr = match find_real_roots(&g) { Some(r) => r, None => return false };
    let df = fr.len();
    let dg = gr.len();
    if dg != df && dg != df + 1 { return false; }
    let eps = 1e-3; // generous tolerance for numerical interlacing
    if df == 0 { return true; }
    if dg == df + 1 {
        for i in 0..df {
            if gr[i] > fr[i] + eps || fr[i] > gr[i + 1] + eps { return false; }
        }
        true
    } else {
        // Same degree: f roots ≤ g roots componentwise
        for i in 0..df {
            if fr[i] > gr[i] + eps { return false; }
        }
        true
    }
}

fn generate_boards(n: usize) -> Vec<Vec<u8>> {
    let mut results = vec![];
    let mut current = vec![];
    gen_rec(n, n, 0, &mut current, &mut results);
    results
}
fn gen_rec(n: usize, max_col: usize, depth: usize, current: &mut Vec<u8>, results: &mut Vec<Vec<u8>>) {
    if depth == n { results.push(current.clone()); return; }
    let min_val = (depth + 1).max(if depth > 0 { current[depth - 1] as usize } else { 1 });
    for v in min_val..=max_col {
        current.push(v as u8);
        gen_rec(n, max_col, depth + 1, current, results);
        current.pop();
    }
}

fn main() {
    let max_n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(7);
    println!("Testing (mP-R) ≼ tR on all Ferrers boards with n <= {}\n", max_n);

    let mut total = 0;
    let mut q_rr_fails = 0;     // Q = mP-R not RR
    let mut tr_rr_fails = 0;    // tR not RR
    let mut interl_fails = 0;   // Q ≼ tR fails
    let mut pplus_rr_fails = 0; // P+ not RR

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
                while p_poly.len() <= pk { p_poly.push(0); }
                p_poly[pk] += 1;
                if pi.len() >= 2 && pi[0] > pi[1] {
                    let k = pi[0] as usize;
                    if k <= m {
                        while d_polys[k].len() <= pk { d_polys[k].push(0); }
                        d_polys[k][pk] += 1;
                    }
                }
            }
            p_poly = poly_trim(&p_poly);

            let mut r_poly = vec![0i64];
            for j in 1..=m {
                r_poly = poly_add(&r_poly, &poly_scale(&d_polys[j], j as i64));
            }

            // Q = mP - R
            let q = poly_sub(&poly_scale(&p_poly, m as i64), &r_poly);
            // tR
            let tr = poly_mul_t(&r_poly);

            // Check Q real-rooted
            if !poly_is_zero(&q) && find_real_roots(&q).is_none() {
                q_rr_fails += 1;
                if q_rr_fails <= 3 {
                    println!("FAIL: Q=mP-R not RR for {:?}: Q={}", board, poly_fmt(&q));
                }
            }

            // Check tR real-rooted
            if !poly_is_zero(&tr) && find_real_roots(&tr).is_none() {
                tr_rr_fails += 1;
            }

            // Check Q ≼ tR (Q roots are more negative than tR roots)
            if !poly_is_zero(&q) && !poly_is_zero(&tr) {
                if !interlaces(&q, &tr) {
                    interl_fails += 1;
                    if interl_fails <= 5 {
                        let qr = find_real_roots(&q).map(|r| format!("{:.4?}", r)).unwrap_or("NOT-RR".into());
                        let trr = find_real_roots(&tr).map(|r| format!("{:.4?}", r)).unwrap_or("NOT-RR".into());
                        println!("FAIL: Q ≼ tR for {:?}: Q={} roots={}, tR={} roots={}",
                            board, poly_fmt(&q), qr, poly_fmt(&tr), trr);
                    }
                }
            }

            // Check P+ = Q + tR is RR
            let pplus = poly_add(&q, &tr);
            if !poly_is_zero(&pplus) && find_real_roots(&pplus).is_none() {
                pplus_rr_fails += 1;
                if pplus_rr_fails <= 3 {
                    println!("FAIL: P+ not RR for {:?}: P+={}", board, poly_fmt(&pplus));
                }
            }
        }
        println!("n={}: total={}, Q_RR={}, tR_RR={}, Q≼tR={}, P+_RR={}",
            n, total, q_rr_fails, tr_rr_fails, interl_fails, pplus_rr_fails);
    }

    println!("\n=== SUMMARY ===");
    println!("Total: {}, Q=mP-R not RR: {}, tR not RR: {}, Q≼tR fails: {}, P+ not RR: {}",
        total, q_rr_fails, tr_rr_fails, interl_fails, pplus_rr_fails);
}
