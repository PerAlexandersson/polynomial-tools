//! Check: does h = (A_{k+1}^+ - A_k^+)/(t-1) ever have negative coefficients?
//! Print examples.
use polynomial_tools::real_rootedness::format_poly;
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

fn main() {
    let mut neg_count = 0u64;
    let mut total = 0u64;
    let mut printed = 0;
    for n in 2..=6usize {
        for board in &gen_boards(n) {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            let m = board[0] as usize;
            let ideal = bruhat_lower_ideal(&perm);
            let subs = sub_partitions(board);
            for mu in &subs {
                // Compute A_k for each first entry k
                let mut a_k: Vec<Vec<i64>> = vec![vec![0i64]; m + 1];
                for sigma in &ideal {
                    let k = sigma[0] as usize;
                    let hits = (0..n)
                        .filter(|&i| sigma[i] as usize > mu[i] as usize)
                        .count();
                    while a_k[k].len() <= hits {
                        a_k[k].push(0);
                    }
                    a_k[k][hits] += 1;
                }
                // Compute consecutive differences A_{k+1} - A_k and divide by (t-1)
                for k in 1..m {
                    let diff = psub(&a_k[k + 1], &a_k[k]);
                    if pz(&diff) {
                        continue;
                    }
                    // Check (t-1) divides
                    let dp = pt(&diff);
                    let sum: i64 = dp.iter().sum();
                    if sum != 0 {
                        continue;
                    } // doesn't divide
                      // Divide by (t-1)
                    let mut h = vec![0i64; dp.len()];
                    h[0] = -dp[0];
                    for i in 1..dp.len() {
                        h[i] = h[i - 1] - dp[i];
                    }
                    let h = pt(&h);
                    if pz(&h) {
                        continue;
                    }
                    total += 1;
                    let has_neg = h.iter().any(|&c| c < 0);
                    if has_neg {
                        neg_count += 1;
                        if printed < 5 {
                            println!("NEG: board={:?} mu={:?} k={}", board, mu, k);
                            println!("  A_{} = {}", k, format_poly(&pt(&a_k[k])));
                            println!("  A_{} = {}", k + 1, format_poly(&pt(&a_k[k + 1])));
                            println!("  diff = {}", format_poly(&dp));
                            println!("  h = {}", format_poly(&h));
                            printed += 1;
                        }
                    }
                }
            }
        }
    }
    println!("\n=== h = (A_{{k+1}} - A_k)/(t-1) sign check ===");
    println!("Total h polynomials: {}", total);
    println!(
        "With negative coefficients: {} ({:.1}%)",
        neg_count,
        100.0 * neg_count as f64 / total as f64
    );
    println!(
        "All non-negative: {} ({:.1}%)",
        total - neg_count,
        100.0 * (total - neg_count) as f64 / total as f64
    );
}
