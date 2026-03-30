//! Check: A_k - A_{k+1} = (t-1)·h with h ≥ 0 (correct sign convention!)
use polynomial_tools::real_rootedness::{format_poly, check_weak_interlacing};
use std::collections::BTreeSet;
fn pt(p: &[i64]) -> Vec<i64> { let mut v = p.to_vec(); while v.len() > 1 && *v.last().unwrap() == 0 { v.pop(); } v }
fn pz(p: &[i64]) -> bool { p.iter().all(|&c| c == 0) }
fn psub(a: &[i64], b: &[i64]) -> Vec<i64> { let l = a.len().max(b.len()); let mut r = vec![0i64; l]; for (i, &v) in a.iter().enumerate() { r[i] += v; } for (i, &v) in b.iter().enumerate() { r[i] -= v; } pt(&r) }
fn pmt(p: &[i64]) -> Vec<i64> { let mut r = vec![0i64; p.len() + 1]; for (i, &v) in p.iter().enumerate() { r[i + 1] = v; } pt(&r) }
fn pdeg(p: &[i64]) -> Option<usize> { let v = pt(p); if pz(&v) { None } else { Some(v.len() - 1) } }
fn interlaces(f: &[i64], g: &[i64]) -> bool { let f = pt(f); let g = pt(g); if pz(&f) { return true; } if pz(&g) { return false; }
    match check_weak_interlacing(&f, &g) { Some(true) => true, Some(false) => false,
        None => { match (pdeg(&f), pdeg(&g)) { (Some(df), Some(dg)) if df == dg => { let tf = pmt(&f); check_weak_interlacing(&g, &tf).unwrap_or(false) }, _ => false, } } } }
fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> { let n = perm.len(); let mut vis: BTreeSet<Vec<u8>> = BTreeSet::new(); let mut q: BTreeSet<Vec<u8>> = BTreeSet::new(); q.insert(perm.to_vec()); while let Some(cur) = q.pop_last() { for i in 0..n { for j in i+1..n { if cur[i] > cur[j] { let mut c = cur.clone(); c.swap(i, j); if !vis.contains(&c) { q.insert(c); } } } } vis.insert(cur); } vis.into_iter().collect() }
fn board_to_perm(b: &[u8]) -> Vec<u8> { let n = b.len(); let mut p = vec![0u8; n]; let mut u = vec![false; n+1]; for i in 0..n { for c in (1..=(b[i] as usize).min(n)).rev() { if !u[c] { p[i] = c as u8; u[c] = true; break; } } } p }
fn is_312_avoiding(perm: &[u8]) -> bool { let n = perm.len(); for i in 0..n { for j in i+1..n { for k in j+1..n { if perm[k] < perm[i] && perm[i] < perm[j] { return false; } } } } true }
fn gen_boards(n: usize) -> Vec<Vec<u8>> { let mut r = vec![]; let mut c = vec![]; gb(n, n, 0, &mut c, &mut r); r }
fn gb(n: usize, mx: usize, d: usize, c: &mut Vec<u8>, r: &mut Vec<Vec<u8>>) { if d == n { r.push(c.clone()); return; } for v in (d+1).max(if d > 0 { c[d-1] as usize } else { 1 })..=mx { c.push(v as u8); gb(n, mx, d+1, c, r); c.pop(); } }
fn sub_partitions(lambda: &[u8]) -> Vec<Vec<u8>> { let n = lambda.len(); let mut result = Vec::new(); let mut mu = vec![0u8; n];
    fn gen(lam: &[u8], mu: &mut Vec<u8>, pos: usize, mx: u8, res: &mut Vec<Vec<u8>>) { if pos == lam.len() { res.push(mu.clone()); return; } let u = lam[pos].min(mx); for v in 0..=u { mu[pos] = v; gen(lam, mu, pos+1, v, res); } }
    gen(lambda, &mut mu, 0, lambda[0], &mut result); result }

fn main() {
    let mut nonneg = [0u64; 2]; // h has non-negative coefficients
    let mut h_interl_ak1 = [0u64; 2]; // h ≪ A_{k+1} (shift lemma condition)
    let mut boundary_ok = [0u64; 2]; // A_{k+1}(0) ≥ h(0)
    let mut total = 0u64;
    for n in 2..=6usize {
        for board in &gen_boards(n) {
            let perm = board_to_perm(board); if !is_312_avoiding(&perm) { continue; }
            let m = board[0] as usize;
            let ideal = bruhat_lower_ideal(&perm);
            let subs = sub_partitions(board);
            for mu in &subs {
                let mut a_k: Vec<Vec<i64>> = vec![vec![0i64]; m+1];
                for sigma in &ideal {
                    let k = sigma[0] as usize;
                    let hits = (0..n).filter(|&i| sigma[i] as usize > mu[i] as usize).count();
                    while a_k[k].len() <= hits { a_k[k].push(0); }
                    a_k[k][hits] += 1;
                }
                for k in 1..m {
                    // A_k - A_{k+1} = (t-1)·h' (correct sign: LARGER index MINUS smaller)
                    let diff = psub(&a_k[k], &a_k[k+1]);
                    if pz(&diff) { continue; }
                    let dp = pt(&diff);
                    let sum: i64 = dp.iter().sum();
                    if sum != 0 { continue; }
                    let mut h = vec![0i64; dp.len()];
                    h[0] = -dp[0];
                    for i in 1..dp.len() { h[i] = h[i-1] - dp[i]; }
                    let h = pt(&h);
                    if pz(&h) { continue; }
                    total += 1;
                    let all_nonneg = h.iter().all(|&c| c >= 0);
                    nonneg[0] += 1;
                    if all_nonneg { nonneg[1] += 1; }
                    // h ≪ A_{k+1}?
                    if !pz(&a_k[k+1]) {
                        h_interl_ak1[0] += 1;
                        if interlaces(&h, &a_k[k+1]) { h_interl_ak1[1] += 1; }
                    }
                    // A_{k+1}(0) ≥ h(0)?
                    let ak1_0 = if a_k[k+1].is_empty() { 0 } else { a_k[k+1][0] };
                    let h_0 = if h.is_empty() { 0 } else { h[0] };
                    boundary_ok[0] += 1;
                    if ak1_0 >= h_0 { boundary_ok[1] += 1; }
                }
            }
        }
    }
    println!("=== A_k - A_{{k+1}} = (t-1)·h (correct sign) ===");
    println!("Total: {}", total);
    let show = |name: &str, c: [u64;2]| { if c[1]==c[0] { println!("  {}: {}/{} ALL PASS <<<", name, c[0], c[0]); } else { println!("  {}: {}/{} pass ({} FAIL)", name, c[1], c[0], c[0]-c[1]); } };
    show("h has non-negative coefficients", nonneg);
    show("h ≪ A_{k+1} (shift lemma condition)", h_interl_ak1);
    show("A_{k+1}(0) ≥ h(0) (boundary)", boundary_ok);
}
