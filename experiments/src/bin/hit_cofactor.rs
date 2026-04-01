//! Verify the cofactor comparison: for each cover μ' = μ+e_j,
//! expand H_{μ'} along row j and check:
//! (1) C_j ≪ C_{j,k} for k ≤ μ_j (cover interlacing on smaller board)
//! (2) C_{j,k} ≪ C_j for k > μ_j+1 (reverse cover on smaller board)
//! (3) C_j ≪ t·C_{j,k} for k > μ_j+1 (Wagner from (2))
//! (4) deg gap between C_j and H_{μ'} is ≤ 1
use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly};
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
fn pmt(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; p.len() + 1];
    for (i, &v) in p.iter().enumerate() {
        r[i + 1] = v;
    }
    pt(&r)
}
fn pdeg(p: &[i64]) -> usize {
    let v = pt(p);
    if pz(&v) {
        0
    } else {
        v.len() - 1
    }
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
fn hit_poly_on_ideal(ideal: &[Vec<u8>], mu: &[u8]) -> Vec<i64> {
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
/// Compute cofactor: permanent with row j and column col removed.
/// mu_reduced: the sub-partition for the remaining rows (adjusted for column removal).
fn cofactor_poly(ideal: &[Vec<u8>], mu: &[u8], row_j: usize, col: usize) -> Vec<i64> {
    let n = mu.len();
    // Filter to perms where sigma[row_j] = col, then count hits on remaining positions
    let mut p = vec![0i64];
    for sigma in ideal {
        if sigma[row_j] as usize != col {
            continue;
        }
        // Count hits on positions i ≠ row_j, with adjusted mu
        let mut hits = 0;
        for i in 0..n {
            if i == row_j {
                continue;
            }
            // After removing column col: the effective value is sigma[i] adjusted
            // The hit condition is sigma[i] > mu[i], but we also need to account
            // for the removed column in the sub-partition adjustment.
            // Actually, the cofactor is just: restrict to sigma[row_j]=col,
            // then count hits normally for all OTHER positions.
            // The sub-partition for position i stays mu[i] in the ORIGINAL matrix.
            // Wait no -- the cofactor is the permanent of the matrix with row j
            // and column col deleted. The hit at position i (i≠j) is whether
            // sigma[i] (in the original column numbering) exceeds mu[i].
            // Since we're just restricting sigma[row_j]=col and counting hits
            // on other positions, the hit condition stays sigma[i] > mu[i].
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
    let mut deg_gap = [0u64; 2]; // deg(H_{μ'}) - deg(C_j) ≤ 1?
    let mut cj_lt_cjk = [0u64; 2]; // C_j ≪ C_{j,k} for k ≤ μ_j
    let mut cjk_lt_cj = [0u64; 2]; // C_{j,k} ≪ C_j for k > μ_j+1
    let mut cj_lt_tcjk = [0u64; 2]; // C_j ≪ t·C_{j,k} for k > μ_j+1
    let mut cj_lt_hmp = [0u64; 2]; // C_j ≪ H_{μ'} (the main condition)
    for n in 2..=max_n {
        for board in &gen_boards(n) {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            let m = *board.iter().max().unwrap() as usize;
            let ideal = bruhat_lower_ideal(&perm);
            let subs = sub_partitions(board);
            for mu in &subs {
                for j in 0..n {
                    // row j (0-indexed)
                    if mu[j] >= board[j] {
                        continue;
                    } // can't increment
                    let mut mu_prime = mu.clone();
                    mu_prime[j] += 1;
                    // Check mu_prime is still weakly decreasing
                    if j > 0 && mu_prime[j] > mu_prime[j - 1] {
                        continue;
                    }
                    let h_mp = hit_poly_on_ideal(&ideal, &mu_prime);
                    let c_j = cofactor_poly(&ideal, &mu_prime, j, mu[j] as usize + 1);
                    if pz(&c_j) {
                        continue;
                    }
                    // Degree gap
                    let dh = pdeg(&h_mp);
                    let dc = pdeg(&c_j);
                    deg_gap[0] += 1;
                    if dh <= dc + 1 {
                        deg_gap[1] += 1;
                    }
                    // C_j ≪ H_{μ'}
                    cj_lt_hmp[0] += 1;
                    if interlaces(&c_j, &h_mp) {
                        cj_lt_hmp[1] += 1;
                    }
                    // Cofactor comparisons
                    for k in 1..=m {
                        if k == mu[j] as usize + 1 {
                            continue;
                        } // skip C_j itself
                        if k > board[j] as usize {
                            continue;
                        }
                        let c_jk = cofactor_poly(&ideal, &mu_prime, j, k);
                        if pz(&c_jk) {
                            continue;
                        }
                        if k <= mu[j] as usize {
                            // k ≤ μ_j: should have C_j ≪ C_{j,k}
                            cj_lt_cjk[0] += 1;
                            if interlaces(&c_j, &c_jk) {
                                cj_lt_cjk[1] += 1;
                            }
                        } else {
                            // k > μ_j+1: should have C_{j,k} ≪ C_j and C_j ≪ t·C_{j,k}
                            cjk_lt_cj[0] += 1;
                            if interlaces(&c_jk, &c_j) {
                                cjk_lt_cj[1] += 1;
                            }
                            let t_cjk = pmt(&c_jk);
                            cj_lt_tcjk[0] += 1;
                            if interlaces(&c_j, &t_cjk) {
                                cj_lt_tcjk[1] += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    println!("=== Cofactor proof structure (n ≤ {}) ===", max_n);
    let show = |name: &str, c: [u64; 2]| {
        if c[0] == 0 {
            println!("  {}: (no data)", name);
        } else if c[1] == c[0] {
            println!("  {}: {}/{} ALL PASS <<<", name, c[0], c[0]);
        } else {
            println!("  {}: {}/{} pass ({} FAIL)", name, c[1], c[0], c[0] - c[1]);
        }
    };
    show("deg(H_{μ'}) ≤ deg(C_j) + 1", deg_gap);
    show("C_j ≪ H_{μ'} (main condition)", cj_lt_hmp);
    show("C_j ≪ C_{j,k} for k ≤ μ_j", cj_lt_cjk);
    show("C_{j,k} ≪ C_j for k > μ_j+1", cjk_lt_cj);
    show("C_j ≪ t·C_{j,k} for k > μ_j+1 (Wagner)", cj_lt_tcjk);
}
