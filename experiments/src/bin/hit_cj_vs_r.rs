//! Test Idea 3: C_j ≪ R directly, where R = H_{μ'} - C_j.
//! Also test: R ≪ H_{μ'} (since H_{μ'} = C_j + R, this is trivial from cone IF C_j ≪ R).
//! And test the BIVARIATE stability: is R + s*C_j real-rooted for s = 2, 3, 10?
//! (If it's real-rooted for large s, that's evidence for C_j ≪ R.)
use polynomial_tools::real_rootedness::{check_weak_interlacing, is_real_rooted};
use std::collections::BTreeSet;
fn pt(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
    v
}
fn pz(p: &[i64]) -> bool {
    p.iter().all(|&c| c == 0)
}
fn pa(a: &[i64], b: &[i64]) -> Vec<i64> {
    let l = a.len().max(b.len());
    let mut r = vec![0i64; l];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] += v;
    }
    pt(&r)
}
fn psub(a: &[i64], b: &[i64]) -> Vec<i64> {
    let l = a.len().max(b.len());
    let mut r = vec![0i64; l];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] -= v;
    }
    pt(&r)
}
fn pscale(p: &[i64], c: i64) -> Vec<i64> {
    pt(&p.iter().map(|&x| x * c).collect::<Vec<_>>())
}
fn interlaces(f: &[i64], g: &[i64]) -> bool {
    let f = pt(f);
    let g = pt(g);
    if pz(&f) {
        return true;
    }
    if pz(&g) {
        return false;
    }
    check_weak_interlacing(&f, &g).unwrap_or(false)
}
fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> {
    let n = perm.len();
    let mut vis: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut q: BTreeSet<Vec<u8>> = BTreeSet::new();
    q.insert(perm.to_vec());
    while let Some(cur) = q.pop_last() {
        for i in 0..n {
            for j in i + 1..n {
                if cur[i] > cur[j] {
                    let mut c = cur.clone();
                    c.swap(i, j);
                    if !vis.contains(&c) {
                        q.insert(c);
                    }
                }
            }
        }
        vis.insert(cur);
    }
    vis.into_iter().collect()
}
fn board_to_perm(b: &[u8]) -> Vec<u8> {
    let n = b.len();
    let mut p = vec![0u8; n];
    let mut u = vec![false; n + 1];
    for i in 0..n {
        for c in (1..=(b[i] as usize).min(n)).rev() {
            if !u[c] {
                p[i] = c as u8;
                u[c] = true;
                break;
            }
        }
    }
    p
}
fn is_312_avoiding(perm: &[u8]) -> bool {
    let n = perm.len();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                if perm[k] < perm[i] && perm[i] < perm[j] {
                    return false;
                }
            }
        }
    }
    true
}
fn gen_boards(n: usize) -> Vec<Vec<u8>> {
    let mut r = vec![];
    let mut c = vec![];
    gb(n, n, 0, &mut c, &mut r);
    r
}
fn gb(n: usize, mx: usize, d: usize, c: &mut Vec<u8>, r: &mut Vec<Vec<u8>>) {
    if d == n {
        r.push(c.clone());
        return;
    }
    for v in (d + 1).max(if d > 0 { c[d - 1] as usize } else { 1 })..=mx {
        c.push(v as u8);
        gb(n, mx, d + 1, c, r);
        c.pop();
    }
}
fn sub_partitions(lambda: &[u8]) -> Vec<Vec<u8>> {
    let n = lambda.len();
    let mut result = Vec::new();
    let mut mu = vec![0u8; n];
    fn gen(lam: &[u8], mu: &mut Vec<u8>, pos: usize, mx: u8, res: &mut Vec<Vec<u8>>) {
        if pos == lam.len() {
            res.push(mu.clone());
            return;
        }
        let u = lam[pos].min(mx);
        for v in 0..=u {
            mu[pos] = v;
            gen(lam, mu, pos + 1, v, res);
        }
    }
    gen(lambda, &mut mu, 0, lambda[0], &mut result);
    result
}
fn hit_poly(ideal: &[Vec<u8>], mu: &[u8]) -> Vec<i64> {
    let n = mu.len();
    let mut p = vec![0i64];
    for sigma in ideal {
        let hits = (0..n)
            .filter(|&i| sigma[i] as usize > mu[i] as usize)
            .count();
        while p.len() <= hits {
            p.push(0);
        }
        p[hits] += 1;
    }
    pt(&p)
}
fn cofactor_poly(ideal: &[Vec<u8>], mu: &[u8], row_j: usize, col: usize) -> Vec<i64> {
    let n = mu.len();
    let mut p = vec![0i64];
    for sigma in ideal {
        if sigma[row_j] as usize != col {
            continue;
        }
        let mut hits = 0;
        for i in 0..n {
            if i == row_j {
                continue;
            }
            if sigma[i] as usize > mu[i] as usize {
                hits += 1;
            }
        }
        while p.len() <= hits {
            p.push(0);
        }
        p[hits] += 1;
    }
    pt(&p)
}
fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);
    let mut cj_r = [0u64; 2]; // C_j ≪ R
    let mut r_rr = [0u64; 2]; // R real-rooted
    let mut rs2 = [0u64; 2]; // R + 2*C_j real-rooted
    let mut rs5 = [0u64; 2]; // R + 5*C_j real-rooted
    let mut rs100 = [0u64; 2]; // R + 100*C_j real-rooted
    for n in 2..=max_n {
        for board in &gen_boards(n) {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            let ideal = bruhat_lower_ideal(&perm);
            for mu in &sub_partitions(board) {
                for j in 0..n {
                    if mu[j] >= board[j] {
                        continue;
                    }
                    let mut mu_prime = mu.clone();
                    mu_prime[j] += 1;
                    if j > 0 && mu_prime[j] > mu_prime[j - 1] {
                        continue;
                    }
                    let h_mp = hit_poly(&ideal, &mu_prime);
                    let col_j = mu[j] as usize + 1;
                    let c_j = cofactor_poly(&ideal, &mu_prime, j, col_j);
                    if pz(&c_j) {
                        continue;
                    }
                    let r = psub(&h_mp, &c_j);
                    if pz(&r) {
                        continue;
                    }
                    // R real-rooted?
                    r_rr[0] += 1;
                    if is_real_rooted(&r) {
                        r_rr[1] += 1;
                    }
                    // C_j ≪ R?
                    cj_r[0] += 1;
                    if interlaces(&c_j, &r) {
                        cj_r[1] += 1;
                    }
                    // R + s*C_j real-rooted for various s?
                    for (s, ctr) in [(2i64, &mut rs2), (5, &mut rs5), (100, &mut rs100)] {
                        let f = pa(&r, &pscale(&c_j, s));
                        ctr[0] += 1;
                        if is_real_rooted(&f) {
                            ctr[1] += 1;
                        }
                    }
                }
            }
        }
    }
    println!("=== Idea 3: C_j ≪ R directly (n ≤ {}) ===", max_n);
    let show = |name: &str, c: [u64; 2]| {
        if c[0] == 0 {
            println!("  {}: (no data)", name);
        } else if c[1] == c[0] {
            println!("  {}: {}/{} ALL PASS <<<", name, c[0], c[0]);
        } else {
            println!("  {}: {}/{} pass ({} FAIL)", name, c[1], c[0], c[0] - c[1]);
        }
    };
    show("R real-rooted", r_rr);
    show("C_j ≪ R", cj_r);
    show("R + 2·C_j real-rooted", rs2);
    show("R + 5·C_j real-rooted", rs5);
    show("R + 100·C_j real-rooted", rs100);
}
