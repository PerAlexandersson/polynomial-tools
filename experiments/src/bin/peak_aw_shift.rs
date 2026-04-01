//! Test: AW at λ+ via the shift lemma.
//! W_l^+ = A_k^+ + (t-1)h_{k,l} where:
//!   k > l: h = S_l + Σ_{j=l}^{k-1} D_j
//!   k < l: h = S_k + Σ_{j=k}^{l-1} U_j
//!   k = l: h = S_l
//! Need: h_{k,l} ≪ A_k^+ and A_k^+(0) ≥ h_{k,l}(0).

use polynomial_tools::real_rootedness::{check_weak_interlacing, is_real_rooted};

fn peaks(w: &[u8]) -> usize {
    if w.len() < 3 {
        return 0;
    }
    (1..w.len() - 1)
        .filter(|&i| w[i - 1] < w[i] && w[i] > w[i + 1])
        .count()
}
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
fn compute_du(board: &[usize]) -> (Vec<Vec<i64>>, Vec<Vec<i64>>) {
    let perms = ferrers_perms(board);
    let n = board.len();
    let m = *board.last().unwrap();
    let mut d = vec![vec![]; m + 1];
    let mut u = vec![vec![]; m + 1];
    for p in &perms {
        if n < 2 {
            continue;
        }
        let k = p[0] as usize;
        let pk = peaks(p);
        let poly = if p[0] > p[1] { &mut d[k] } else { &mut u[k] };
        while poly.len() <= pk {
            poly.push(0);
        }
        poly[pk] += 1;
    }
    (d, u)
}
fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let n = a.len().max(b.len());
    let mut r = vec![0i64; n];
    for i in 0..a.len() {
        r[i] += a[i];
    }
    for i in 0..b.len() {
        r[i] += b[i];
    }
    r
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
fn deg(p: &[i64]) -> usize {
    let t = trim(p);
    if t.is_empty() {
        0
    } else {
        t.len() - 1
    }
}
fn interlaces_weak(f: &[i64], g: &[i64]) -> bool {
    let f = trim(f);
    let g = trim(g);
    if f.is_empty() {
        return is_real_rooted(&g);
    }
    if g.is_empty() {
        return false;
    }
    let (df, dg) = (deg(&f), deg(&g));
    if dg == df + 1 {
        check_weak_interlacing(&f, &g) == Some(true)
    } else if dg == df {
        let tf = poly_tmul(&f);
        check_weak_interlacing(&g, &tf) == Some(true)
    } else {
        false
    }
}
fn eval0(p: &[i64]) -> i64 {
    if p.is_empty() {
        0
    } else {
        p[0]
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
    println!("=== AW via shift lemma: h_{{k,l}} ≪ A_k^+ ===\n");
    let (mut total, mut pass) = (0usize, 0usize);
    let (mut total_bd, mut pass_bd) = (0usize, 0usize); // boundary condition

    for n in 2..=7 {
        let boards = boards_312(n);
        let (mut nt, mut np, mut ntb, mut npb) = (0, 0, 0, 0);
        for board in &boards {
            let m = *board.last().unwrap();
            let (dp, up) = compute_du(board);
            let mut a_poly = vec![vec![]; m + 1];
            let mut w_poly = vec![vec![]; m + 1];
            for j in 1..=m {
                a_poly[j] = poly_add(&dp[j], &up[j]);
                w_poly[j] = poly_add(&poly_tmul(&dp[j]), &up[j]);
            }
            let m_prime = m + 1;
            let mut s = vec![vec![]; m_prime + 1];
            let mut t_sum = vec![vec![]; m_prime + 1];
            let mut a_plus = vec![vec![]; m_prime + 1];
            for k in 1..=m_prime {
                let mut sk = vec![];
                let mut tk = vec![];
                for j in 1..=m {
                    if j < k {
                        sk = poly_add(&sk, &a_poly[j]);
                    }
                    if j >= k {
                        tk = poly_add(&tk, &w_poly[j]);
                    }
                }
                s[k] = sk.clone();
                t_sum[k] = tk.clone();
                a_plus[k] = poly_add(&sk, &tk);
            }

            for k in 1..=m_prime {
                for l in 1..=m_prime {
                    let apk = trim(&a_plus[k]);
                    if apk.is_empty() {
                        continue;
                    }

                    // Compute h_{k,l}
                    let h = if k > l {
                        // h = S_l + Σ_{j=l}^{k-1} D_j
                        let mut hh = s[l].clone();
                        for j in l..k {
                            if j <= m {
                                hh = poly_add(&hh, &dp[j]);
                            }
                        }
                        hh
                    } else if k < l {
                        // h = S_k + Σ_{j=k}^{l-1} U_j
                        let mut hh = s[k].clone();
                        for j in k..l {
                            if j <= m {
                                hh = poly_add(&hh, &up[j]);
                            }
                        }
                        hh
                    } else {
                        // h = S_l
                        s[l].clone()
                    };
                    let h = trim(&h);
                    if h.is_empty() {
                        continue;
                    }

                    // Test h ≪ A_k^+
                    nt += 1;
                    if interlaces_weak(&h, &apk) {
                        np += 1;
                    } else {
                        if nt - np <= 3 {
                            println!("  h≪A FAIL: board={:?}, k={}, l={}", board, k, l);
                        }
                    }

                    // Test boundary: A_k^+(0) ≥ h(0)
                    ntb += 1;
                    if eval0(&apk) >= eval0(&h) {
                        npb += 1;
                    }
                }
            }
        }
        total += nt;
        pass += np;
        total_bd += ntb;
        pass_bd += npb;
        println!("n={}: h≪A+={}/{}, boundary={}/{}", n, np, nt, npb, ntb);
    }
    println!("\nh≪A+:     {}/{}", pass, total);
    println!("boundary: {}/{}", pass_bd, total_bd);
}
