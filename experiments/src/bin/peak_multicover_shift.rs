//! Test: for multi-cover chains nu_k -> nu,
//! does sum(C_i) << H_nu hold?
//! This would let the shift lemma handle the composite.

use polynomial_tools::real_rootedness::check_weak_interlacing;

fn all_perms(n: u8) -> Vec<Vec<u8>> {
    if n <= 1 { return vec![(1..=n).collect()]; }
    let mut r = Vec::new();
    for p in all_perms(n-1) {
        for i in 0..=p.len() {
            let mut q = p.clone(); q.insert(i, n); r.push(q);
        }
    }
    r
}

fn ferrers_perms(board: &[usize]) -> Vec<Vec<u8>> {
    let n = board.len();
    all_perms(n as u8).into_iter()
        .filter(|p| (0..n).all(|i| (p[i] as usize) <= board[i]))
        .collect()
}

fn hit_poly(board: &[usize], mu: &[usize]) -> Vec<i64> {
    let n = board.len();
    let perms = ferrers_perms(board);
    let mut coeffs = vec![0i64; n + 1];
    for p in &perms {
        let hits = (0..n).filter(|&i| p[i] as usize > mu[i]).count();
        coeffs[hits] += 1;
    }
    while coeffs.last() == Some(&0) { coeffs.pop(); }
    coeffs
}

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let n = a.len().max(b.len());
    let mut r = vec![0i64; n];
    for i in 0..a.len() { r[i] += a[i]; }
    for i in 0..b.len() { r[i] += b[i]; }
    r
}

fn poly_tmul(a: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; a.len() + 1];
    for i in 0..a.len() { r[i+1] = a[i]; }
    r
}

fn trim(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.last() == Some(&0) { v.pop(); }
    v
}

fn interlaces_weak(f: &[i64], g: &[i64]) -> bool {
    let f = trim(f); let g = trim(g);
    if f.is_empty() { return true; }
    if g.is_empty() { return false; }
    let df = f.len() - 1; let dg = g.len() - 1;
    if dg == df + 1 {
        check_weak_interlacing(&f, &g) == Some(true)
    } else if dg == df {
        let tf = poly_tmul(&f);
        check_weak_interlacing(&g, &tf) == Some(true)
    } else { false }
}

/// Get cofactor: delete row j, col k from board/mu
fn cofactor_board_mu(
    board: &[usize], mu: &[usize], j: usize, k: usize
) -> (Vec<usize>, Vec<usize>) {
    let n = board.len();
    let mut new_board = Vec::new();
    let mut new_mu = Vec::new();
    for i in 0..n {
        if i == j { continue; }
        // Delete column k: board[i] -> board[i] - [k <= board[i]]
        let b = board[i] - if k <= board[i] { 1 } else { 0 };
        let m = mu[i] - if k <= mu[i] { 1 } else { 0 };
        new_board.push(b);
        new_mu.push(m);
    }
    (new_board, new_mu)
}

fn boards_312(n: usize) -> Vec<Vec<usize>> {
    fn gen(n: usize, b: &mut Vec<usize>, r: &mut Vec<Vec<usize>>) {
        if b.len() == n { r.push(b.clone()); return; }
        let i = b.len();
        let prev = b.last().copied().unwrap_or(i+1).max(i+1);
        for v in prev..=n { b.push(v); gen(n, b, r); b.pop(); }
    }
    let mut r = Vec::new(); let mut b = Vec::new();
    gen(n, &mut b, &mut r); r
}

fn main() {
    println!("=== Multi-cover shift lemma test ===\n");
    let mut total = 0;
    let mut pass = 0;

    for n in 2..=5 {
        let boards = boards_312(n);
        for board in &boards {
            let m = *board.last().unwrap();
            // Generate all valid sub-partitions
            fn gen_mu(board: &[usize], idx: usize,
                      prev: usize, mu: &mut Vec<usize>,
                      result: &mut Vec<Vec<usize>>) {
                if idx == board.len() {
                    result.push(mu.clone()); return;
                }
                let max_val = prev.min(board[idx]);
                for v in 0..=max_val {
                    mu.push(v);
                    gen_mu(board, idx+1, v, mu, result);
                    mu.pop();
                }
            }
            let mut all_mu = Vec::new();
            gen_mu(board, 0, m, &mut Vec::new(), &mut all_mu);

            for mu in &all_mu {
                // Find all multi-cover pairs (nu_k, nu)
                // with |nu - nu_k| >= 2
                for nu in &all_mu {
                    // Check nu >= mu componentwise
                    // and diff >= 2
                    let diff: usize = (0..n)
                        .map(|i| nu[i] - mu[i])
                        .sum();
                    if diff < 2 { continue; }
                    if (0..n).any(|i| nu[i] < mu[i]) { continue; }
                    if (0..n).any(|i| nu[i] > board[i]) { continue; }

                    // Compute H_nu and the telescoping sum of cofactors
                    let h_nu = hit_poly(board, nu);
                    let h_mu = hit_poly(board, mu);
                    if h_nu.is_empty() || h_mu.is_empty() { continue; }

                    // H_mu - H_nu = (t-1) * sum_cofactors
                    // So sum_cofactors = (H_mu - H_nu) / (t-1)
                    let mut diff_poly = vec![0i64;
                        h_mu.len().max(h_nu.len())];
                    for i in 0..h_mu.len() { diff_poly[i] += h_mu[i]; }
                    for i in 0..h_nu.len() { diff_poly[i] -= h_nu[i]; }
                    let diff_poly = trim(&diff_poly);

                    // Divide by (t-1): if p = (t-1)*q,
                    // then q[0] = -p[0] and
                    // q[i] = q[i-1] - p[i] ... wait, synthetic division
                    // p(t) = (t-1) * q(t)
                    // q[d-1] = p[d]
                    // q[i] = p[i+1] + q[i+1] for i < d-1
                    if diff_poly.is_empty() { continue; }
                    let d = diff_poly.len() - 1;
                    let mut q = vec![0i64; d];
                    if d == 0 { continue; }
                    q[d-1] = diff_poly[d];
                    for i in (0..d-1).rev() {
                        q[i] = diff_poly[i+1] + q[i+1];
                    }
                    // Verify: diff_poly[0] should equal -q[0]
                    if diff_poly[0] != -q[0] {
                        continue; // not divisible by (t-1)
                    }
                    let sum_cof = trim(&q);
                    if sum_cof.is_empty() { continue; }

                    total += 1;
                    if interlaces_weak(&sum_cof, &h_nu) {
                        pass += 1;
                    } else {
                        if total - pass <= 5 {
                            println!("  FAIL: board={:?} mu={:?} nu={:?}",
                                board, mu, nu);
                            println!("    sum_cof={:?} H_nu={:?}",
                                sum_cof, h_nu);
                        }
                    }
                }
            }
        }
        println!("n={}: {}/{}", n, pass, total);
    }
    println!("\nTotal: {}/{}", pass, total);
}
