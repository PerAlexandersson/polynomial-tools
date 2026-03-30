//! Test the recursion P^{λ+} = m·P^λ + (t-1)·R^λ
//! where R^λ = Σ_j j·D_j^λ, and study P/R interlacing.

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
fn poly_scale(p: &[i64], c: i64) -> Vec<i64> {
    poly_trim(&p.iter().map(|&x| x * c).collect::<Vec<_>>())
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
    let eps = 1e-6;
    if df == 0 { return true; }
    if dg == df + 1 {
        for i in 0..df {
            if gr[i] > fr[i] + eps || fr[i] > gr[i + 1] + eps { return false; }
        }
        true
    } else {
        for i in 0..df {
            if gr[i] > fr[i] + eps { return false; }
            if i + 1 < df && fr[i] > gr[i + 1] + eps { return false; }
        }
        true
    }
}

fn roots_str(p: &[i64]) -> String {
    match find_real_roots(p) {
        Some(r) if r.is_empty() => "[]".into(),
        Some(r) => format!("{:?}", r.iter().map(|x| format!("{:.4}", x)).collect::<Vec<_>>()),
        None => "NOT-RR".into(),
    }
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
    println!("Testing P/R recursion on all Ferrers boards with n <= {}\n", max_n);

    let mut total = 0;
    let mut r_rr_fails = 0;
    let mut r_interl_p_fails = 0; // R ≼ P
    let mut p_interl_r_fails = 0; // P ≼ R
    let mut recursion_rr_fails = 0;

    for n in 2..=max_n {
        let boards = generate_boards(n);
        for board in &boards {
            total += 1;
            let m = (board[0] as usize).min(n);
            let perm = board_to_312_perm(board);
            let ideal = bruhat_lower_ideal(&perm);

            // Compute D_k for this board
            let mut d_polys = vec![vec![0i64]; m + 1]; // 1-indexed
            for pi in &ideal {
                if pi.len() < 2 { continue; }
                let k = pi[0] as usize;
                if k > m { continue; }
                let pk = peaks(pi);
                let is_descent = pi[0] > pi[1];
                if is_descent {
                    while d_polys[k].len() <= pk { d_polys[k].push(0); }
                    d_polys[k][pk] += 1;
                }
            }

            // Compute R = Σ j*D_j and P = total peak poly
            let mut r_poly = vec![0i64];
            let mut p_poly = vec![0i64];
            for pi in &ideal {
                let pk = peaks(pi);
                while p_poly.len() <= pk { p_poly.push(0); }
                p_poly[pk] += 1;
            }
            for j in 1..=m {
                r_poly = poly_add(&r_poly, &poly_scale(&d_polys[j], j as i64));
            }

            p_poly = poly_trim(&p_poly);
            r_poly = poly_trim(&r_poly);

            // Check R is real-rooted
            if !poly_is_zero(&r_poly) && find_real_roots(&r_poly).is_none() {
                r_rr_fails += 1;
                if r_rr_fails <= 5 {
                    println!("FAIL: R not RR for {:?}: R={}", board, poly_fmt(&r_poly));
                }
            }

            // Check interlacing: R ≼ P?
            if !poly_is_zero(&r_poly) && !poly_is_zero(&p_poly) {
                if !interlaces(&r_poly, &p_poly) {
                    r_interl_p_fails += 1;
                    if r_interl_p_fails <= 5 {
                        println!("FAIL: R ≼ P for {:?}: R={} roots={}, P={} roots={}",
                            board, poly_fmt(&r_poly), roots_str(&r_poly),
                            poly_fmt(&p_poly), roots_str(&p_poly));
                    }
                }
                // Check P ≼ R?
                if !interlaces(&p_poly, &r_poly) {
                    p_interl_r_fails += 1;
                }
            }

            // Verify the recursion: for each valid m' ≤ m, check P+ is RR
            // (just test m' = m)
            // P+ = m*P + (t-1)*R
            let mp = poly_scale(&p_poly, m as i64);
            // (t-1)*R: shift R by t, then subtract R
            let mut tr = vec![0i64; r_poly.len() + 1];
            for (i, &v) in r_poly.iter().enumerate() { tr[i + 1] = v; }
            let neg_r: Vec<i64> = r_poly.iter().map(|&x| -x).collect();
            let t_minus_1_r = poly_add(&tr, &neg_r);
            let p_plus = poly_trim(&poly_add(&mp, &t_minus_1_r));

            if !poly_is_zero(&p_plus) && find_real_roots(&p_plus).is_none() {
                recursion_rr_fails += 1;
                if recursion_rr_fails <= 3 {
                    println!("FAIL: P+ not RR for {:?}: P+={}", board, poly_fmt(&p_plus));
                }
            }
        }
        println!("n={}: boards={}, R_RR_fails={}, R≼P_fails={}, P≼R_fails={}, P+_RR_fails={}",
            n, total, r_rr_fails, r_interl_p_fails, p_interl_r_fails, recursion_rr_fails);
    }

    println!("\n=== SUMMARY ===");
    println!("Total boards: {}", total);
    println!("R not real-rooted: {}", r_rr_fails);
    println!("R ≼ P fails: {}", r_interl_p_fails);
    println!("P ≼ R fails: {}", p_interl_r_fails);
    println!("P+ = mP+(t-1)R not RR: {}", recursion_rr_fails);
}
