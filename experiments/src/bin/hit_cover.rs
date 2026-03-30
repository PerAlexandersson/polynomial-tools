//! Test: (1) multi-cover interlacing, (2) H_μ - H_{μ+e_j} structure
use polynomial_tools::real_rootedness::{check_weak_interlacing, is_real_rooted, format_poly};
use std::collections::BTreeSet;
fn pt(p: &[i64]) -> Vec<i64> { let mut v = p.to_vec(); while v.len() > 1 && *v.last().unwrap() == 0 { v.pop(); } v }
fn pz(p: &[i64]) -> bool { p.iter().all(|&c| c == 0) }
fn pa(a: &[i64], b: &[i64]) -> Vec<i64> { let l = a.len().max(b.len()); let mut r = vec![0i64; l]; for (i, &v) in a.iter().enumerate() { r[i] += v; } for (i, &v) in b.iter().enumerate() { r[i] += v; } pt(&r) }
fn pmt(p: &[i64]) -> Vec<i64> { let mut r = vec![0i64; p.len() + 1]; for (i, &v) in p.iter().enumerate() { r[i + 1] = v; } pt(&r) }
fn psub(a: &[i64], b: &[i64]) -> Vec<i64> { let l = a.len().max(b.len()); let mut r = vec![0i64; l]; for (i, &v) in a.iter().enumerate() { r[i] += v; } for (i, &v) in b.iter().enumerate() { r[i] -= v; } pt(&r) }
fn pdeg(p: &[i64]) -> Option<usize> { let v = pt(p); if pz(&v) { None } else { Some(v.len() - 1) } }
fn interlaces(f: &[i64], g: &[i64]) -> bool { let f = pt(f); let g = pt(g); if pz(&f) { return true; } if pz(&g) { return false; }
    check_weak_interlacing(&f, &g).unwrap_or(false) }
fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> { let n = perm.len(); let mut vis: BTreeSet<Vec<u8>> = BTreeSet::new(); let mut q: BTreeSet<Vec<u8>> = BTreeSet::new(); q.insert(perm.to_vec()); while let Some(cur) = q.pop_last() { for i in 0..n { for j in i+1..n { if cur[i] > cur[j] { let mut c = cur.clone(); c.swap(i, j); if !vis.contains(&c) { q.insert(c); } } } } vis.insert(cur); } vis.into_iter().collect() }
fn board_to_perm(b: &[u8]) -> Vec<u8> { let n = b.len(); let mut p = vec![0u8; n]; let mut u = vec![false; n+1]; for i in 0..n { for c in (1..=(b[i] as usize).min(n)).rev() { if !u[c] { p[i] = c as u8; u[c] = true; break; } } } p }
fn is_312_avoiding(perm: &[u8]) -> bool { let n = perm.len(); for i in 0..n { for j in i+1..n { for k in j+1..n { if perm[k] < perm[i] && perm[i] < perm[j] { return false; } } } } true }
fn gen_boards(n: usize) -> Vec<Vec<u8>> { let mut r = vec![]; let mut c = vec![]; gb(n, n, 0, &mut c, &mut r); r }
fn gb(n: usize, mx: usize, d: usize, c: &mut Vec<u8>, r: &mut Vec<Vec<u8>>) { if d == n { r.push(c.clone()); return; } for v in (d+1).max(if d > 0 { c[d-1] as usize } else { 1 })..=mx { c.push(v as u8); gb(n, mx, d+1, c, r); c.pop(); } }
fn sub_partitions(lambda: &[u8]) -> Vec<Vec<u8>> { let n = lambda.len(); let mut result = Vec::new(); let mut mu = vec![0u8; n];
    fn gen(lam: &[u8], mu: &mut Vec<u8>, pos: usize, mx: u8, res: &mut Vec<Vec<u8>>) { if pos == lam.len() { res.push(mu.clone()); return; } let u = lam[pos].min(mx); for v in 0..=u { mu[pos] = v; gen(lam, mu, pos+1, v, res); } }
    gen(lambda, &mut mu, 0, lambda[0], &mut result); result }
fn hit_poly(ideal: &[Vec<u8>], mu: &[u8]) -> Vec<i64> { let n = mu.len(); let mut p = vec![0i64];
    for sigma in ideal { let hits = (0..n).filter(|&i| sigma[i] as usize > mu[i] as usize).count();
        while p.len() <= hits { p.push(0); } p[hits] += 1; } pt(&p) }

fn main() {
    let max_n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(6);
    // Test 1: H_μ - H_{μ+e_j} = (t-1) · C_j for some C_j with nonneg coefficients?
    let mut diff_t1 = [0u64; 2]; // Does (H_μ - H_{μ'})/(t-1) have nonneg coeffs?
    let mut diff_rr = [0u64; 2]; // Is H_μ - H_{μ'} real-rooted?
    let mut cofactor_rr = [0u64; 2]; // Is C_j = (H_μ - H_{μ'})/(t-1) real-rooted?
    // Test 2: Multi-cover (differ in exactly 2 parts by 1 each)
    let mut multi2 = [0u64; 2]; // H_μ ≪ H_{μ'} for 2-cover
    let mut multi2_rev = [0u64; 2];
    // Test 3: Shift lemma structure: H_{μ'} = H_μ + (t-1)·C, C ≪ H_μ?
    // Actually H_μ = H_{μ'} + (t-1)·C (since H_μ has MORE t's)
    let mut shift_lemma = [0u64; 2]; // C ≪ H_{μ'}?
    
    for n in 2..=max_n {
        for board in &gen_boards(n) {
            let perm = board_to_perm(board); if !is_312_avoiding(&perm) { continue; }
            let ideal = bruhat_lower_ideal(&perm);
            let subs = sub_partitions(board);
            let polys: Vec<Vec<i64>> = subs.iter().map(|mu| hit_poly(&ideal, mu)).collect();
            
            for i in 0..subs.len() { for j in 0..subs.len() { if i == j { continue; }
                let diff: Vec<(usize, u8)> = (0..n).filter(|&k| subs[i][k] != subs[j][k])
                    .map(|k| (k, subs[j][k] - subs[i][k])).collect();
                // Cover: exactly one part increases by 1
                if diff.len() == 1 && diff[0].1 == 1 {
                    // H_μ has MORE hits than H_{μ'} (since μ < μ'), so H_μ has higher degree
                    // H_μ - H_{μ'}: check if = (t-1)·C
                    let d = psub(&polys[i], &polys[j]); // H_μ - H_{μ'}
                    if !pz(&d) {
                        diff_rr[0] += 1;
                        if is_real_rooted(&pt(&d)) { diff_rr[1] += 1; }
                        // Check (t-1) divisibility: d = (t-1)·c means d(1) = 0
                        let d_at_1: i64 = pt(&d).iter().sum();
                        diff_t1[0] += 1;
                        if d_at_1 == 0 {
                            // Divide by (t-1): polynomial long division
                            let dp = pt(&d);
                            let mut c = vec![0i64; dp.len()];
                            // d = (t-1)·c: c[k] = -(d[0] + c[0] + c[1] + ... + c[k-1]) / (-1)
                            // Actually: (t-1)·c = t·c - c. Coeff of t^k: c[k-1] - c[k] = d[k].
                            // So c[k] = c[k-1] - d[k], with c[-1] = 0.
                            // c[0] = -d[0]. c[1] = c[0]-d[1]. etc.
                            c[0] = -dp[0];
                            for k in 1..dp.len() { c[k] = c[k-1] - dp[k]; }
                            let c = pt(&c);
                            cofactor_rr[0] += 1;
                            if !pz(&c) && (is_real_rooted(&c) || c.len() <= 2) { cofactor_rr[1] += 1; }
                            // Shift lemma: C ≪ H_{μ'}?
                            if !pz(&c) && !pz(&polys[j]) {
                                shift_lemma[0] += 1;
                                if interlaces(&c, &polys[j]) { shift_lemma[1] += 1; }
                            }
                            diff_t1[1] += 1; // (t-1) divides
                        }
                    }
                }
                // 2-cover: exactly 2 parts increase by 1
                if diff.len() == 2 && diff.iter().all(|&(_, d)| d == 1) {
                    if !pz(&polys[i]) && !pz(&polys[j]) {
                        multi2[0] += 1;
                        if interlaces(&polys[i], &polys[j]) { multi2[1] += 1; }
                        multi2_rev[0] += 1;
                        if interlaces(&polys[j], &polys[i]) { multi2_rev[1] += 1; }
                    }
                }
            }}
        }
    }
    println!("=== Cover structure (n ≤ {}) ===", max_n);
    let show = |name: &str, c: [u64;2]| { if c[0]==0 { println!("  {}: (no data)", name); } else if c[1]==0 { println!("  {}: 0/{} (ALL FAIL)", name, c[0]); } else if c[1]==c[0] { println!("  {}: {}/{} ALL PASS <<<", name, c[0], c[0]); } else { println!("  {}: {}/{} pass ({} FAIL)", name, c[1], c[0], c[0]-c[1]); } };
    show("H_μ - H_{μ'} divisible by (t-1)", diff_t1);
    show("H_μ - H_{μ'} real-rooted", diff_rr);
    show("Cofactor C = (H_μ - H_{μ'})/(t-1) real-rooted", cofactor_rr);
    show("Shift lemma: C ≪ H_{μ'}", shift_lemma);
    show("2-cover: H_μ ≪ H_{μ'}", multi2);
    show("2-cover: H_{μ'} ≪ H_μ", multi2_rev);
}
