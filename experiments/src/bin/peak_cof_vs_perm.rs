//! Test: does EVERY cofactor of M_mu interlace perm(M_mu)?
//! I.e., for every (j,k) with M[j,k] != 0:
//!   cof_{j,k}(M_mu) << perm(M_mu)
//! This is the DIRECT condition, avoiding the row expansion.

use polynomial_tools::real_rootedness::check_weak_interlacing;

fn all_perms(n: u8) -> Vec<Vec<u8>> {
    if n <= 1 {
        return vec![(1..=n).collect()];
    }
    let mut r = Vec::new();
    for p in all_perms(n - 1) {
        for i in 0..=p.len() {
            let mut q = p.clone();
            q.insert(i, n);
            r.push(q);
        }
    }
    r
}

fn ferrers_perms(board: &[usize]) -> Vec<Vec<u8>> {
    let n = board.len();
    all_perms(n as u8)
        .into_iter()
        .filter(|p| (0..n).all(|i| (p[i] as usize) <= board[i]))
        .collect()
}

fn hit_poly(board: &[usize], mu: &[usize]) -> Vec<i64> {
    let n = board.len();
    let perms = ferrers_perms(board);
    let mut coeffs = vec![0i64; n + 2];
    for p in &perms {
        let hits = (0..n).filter(|&i| p[i] as usize > mu[i]).count();
        coeffs[hits] += 1;
    }
    while coeffs.last() == Some(&0) {
        coeffs.pop();
    }
    coeffs
}

fn poly_tmul(a: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; a.len() + 1];
    for i in 0..a.len() {
        r[i + 1] = a[i];
    }
    r
}

fn trim(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.last() == Some(&0) {
        v.pop();
    }
    v
}

fn interlaces_weak(f: &[i64], g: &[i64]) -> bool {
    let f = trim(f);
    let g = trim(g);
    if f.is_empty() {
        return true;
    }
    if g.is_empty() {
        return false;
    }
    let df = f.len() - 1;
    let dg = g.len() - 1;
    if dg == df + 1 {
        check_weak_interlacing(&f, &g) == Some(true)
    } else if dg == df {
        let tf = poly_tmul(&f);
        check_weak_interlacing(&g, &tf) == Some(true)
    } else {
        false
    }
}

fn boards_312(n: usize) -> Vec<Vec<usize>> {
    fn gen(n: usize, b: &mut Vec<usize>, r: &mut Vec<Vec<usize>>) {
        if b.len() == n {
            r.push(b.clone());
            return;
        }
        let i = b.len();
        let prev = b.last().copied().unwrap_or(i + 1).max(i + 1);
        for v in prev..=n {
            b.push(v);
            gen(n, b, r);
            b.pop();
        }
    }
    let mut r = Vec::new();
    let mut b = Vec::new();
    gen(n, &mut b, &mut r);
    r
}

fn main() {
    println!("=== Cofactor << Permanent (all positions) ===\n");
    let mut total = 0;
    let mut pass = 0;

    for n in 2..=6 {
        let boards = boards_312(n);
        let (mut nt, mut np) = (0, 0);

        for board in &boards {
            let m = *board.last().unwrap();
            fn gen_mu(
                board: &[usize],
                idx: usize,
                prev: usize,
                mu: &mut Vec<usize>,
                result: &mut Vec<Vec<usize>>,
            ) {
                if idx == board.len() {
                    result.push(mu.clone());
                    return;
                }
                let max_val = prev.min(board[idx]);
                for v in 0..=max_val {
                    mu.push(v);
                    gen_mu(board, idx + 1, v, mu, result);
                    mu.pop();
                }
            }
            let mut all_mu = Vec::new();
            gen_mu(board, 0, m, &mut Vec::new(), &mut all_mu);

            for mu in &all_mu {
                let h = hit_poly(board, mu);
                if h.is_empty() || h.len() <= 1 {
                    continue;
                }

                // Test every position (j, k) with M[j,k] != 0
                for j in 0..n {
                    for k in 1..=board[j] {
                        // Cofactor: delete row j, col k
                        let board_del: Vec<usize> = (0..n)
                            .filter(|&i| i != j)
                            .map(|i| board[i] - if k <= board[i] { 1 } else { 0 })
                            .collect();
                        let mu_del: Vec<usize> = (0..n)
                            .filter(|&i| i != j)
                            .map(|i| mu[i] - if k <= mu[i] { 1 } else { 0 })
                            .collect();

                        let cof = hit_poly(&board_del, &mu_del);
                        if cof.is_empty() {
                            continue;
                        }

                        nt += 1;
                        if interlaces_weak(&cof, &h) {
                            np += 1;
                        } else {
                            if nt - np <= 3 {
                                println!(
                                    "  FAIL: board={:?} mu={:?} j={} k={}",
                                    board,
                                    mu,
                                    j + 1,
                                    k
                                );
                                println!("    cof={:?} H={:?}", cof, h);
                            }
                        }
                    }
                }
            }
        }
        total += nt;
        pass += np;
        println!("n={}: {}/{}", n, np, nt);
    }
    println!("\nTotal: {}/{}", pass, total);
}
