//! Test AW propagation: does AW at λ imply AW at λ+?
//! A_k^+ = S_k + T_k, W_l^+ = tS_l + T_l.
//! Need: A_k^+ ≪ W_l^+ for all k, l.
//!
//! We know from the shift lemma:
//!   W_l^+ = A_l^+ + (t-1)S_l, so if S_l ≪ A_l^+, then A_l^+ ≪ W_l^+ (diagonal).
//!
//! For k ≠ l: need A_k^+ ≪ W_l^+ = tS_l + T_l.
//! By cone: A_k^+ ≪ tS_l (equiv: S_l ≪ A_k^+ by Wagner) AND A_k^+ ≪ T_l.
//!
//! Test which sub-conditions hold.

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
    println!("=== AW propagation sub-conditions ===\n");
    // A_k^+ ≪ W_l^+ = tS_l + T_l.
    // Decompose: A_k^+ ≪ T_l AND S_l ≪ A_k^+ (Wagner for tS_l part)

    let mut aw_plus = (0usize, 0usize); // A_k^+ ≪ W_l^+ (the full condition)
    let mut ak_tl = (0usize, 0usize); // A_k^+ ≪ T_l
    let mut sl_ak = (0usize, 0usize); // S_l ≪ A_k^+ (equiv to A_k^+ ≪ tS_l via Wagner)

    for n in 2..=7 {
        let boards = boards_312(n);
        for board in &boards {
            let m = *board.last().unwrap();
            let (dp, up) = compute_du(board);
            let mut a = vec![vec![]; m + 1];
            let mut w = vec![vec![]; m + 1];
            for j in 1..=m {
                a[j] = poly_add(&dp[j], &up[j]);
                w[j] = poly_add(&poly_tmul(&dp[j]), &up[j]);
            }

            let m_prime = m + 1;
            let mut s = vec![vec![]; m_prime + 1];
            let mut t_sum = vec![vec![]; m_prime + 1];
            let mut a_plus = vec![vec![]; m_prime + 1];
            let mut w_plus = vec![vec![]; m_prime + 1];
            for k in 1..=m_prime {
                let mut sk = vec![];
                let mut tk = vec![];
                for j in 1..=m {
                    if j < k {
                        sk = poly_add(&sk, &a[j]);
                    }
                    if j >= k {
                        tk = poly_add(&tk, &w[j]);
                    }
                }
                s[k] = sk.clone();
                t_sum[k] = tk.clone();
                a_plus[k] = poly_add(&sk, &tk);
                w_plus[k] = poly_add(&poly_tmul(&sk), &tk);
            }

            for k in 1..=m_prime {
                for l in 1..=m_prime {
                    let apk = trim(&a_plus[k]);
                    let wpl = trim(&w_plus[l]);
                    let tl = trim(&t_sum[l]);
                    let sl = trim(&s[l]);
                    if apk.is_empty() || wpl.is_empty() {
                        continue;
                    }

                    // Full AW at λ+
                    aw_plus.0 += 1;
                    if interlaces_weak(&apk, &wpl) {
                        aw_plus.1 += 1;
                    }

                    // A_k^+ ≪ T_l
                    if !tl.is_empty() {
                        ak_tl.0 += 1;
                        if interlaces_weak(&apk, &tl) {
                            ak_tl.1 += 1;
                        }
                    }

                    // S_l ≪ A_k^+
                    if !sl.is_empty() {
                        sl_ak.0 += 1;
                        if interlaces_weak(&sl, &apk) {
                            sl_ak.1 += 1;
                        }
                    }
                }
            }
        }
        println!(
            "n={}: AW+={}/{}, A+≪T={}/{}, S≪A+={}/{}",
            n, aw_plus.1, aw_plus.0, ak_tl.1, ak_tl.0, sl_ak.1, sl_ak.0
        );
    }
    println!("\nAW at λ+:   {}/{}", aw_plus.1, aw_plus.0);
    println!("A_k^+≪T_l:  {}/{}", ak_tl.1, ak_tl.0);
    println!("S_l≪A_k^+:  {}/{}", sl_ak.1, sl_ak.0);
}
