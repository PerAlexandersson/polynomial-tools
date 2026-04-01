//! Test: does H_μ^λ ≪ H_{μ'}^λ when μ ⊆ μ' (μ_i ≤ μ'_i for all i)?
//! "Monotonicity" in the sub-partition.
//! Also test: does H_{μ'}^λ ≪ t·H_μ^λ (Wagner-like)?
use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly, is_real_rooted};
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
fn pdeg(p: &[i64]) -> Option<usize> {
    let v = pt(p);
    if pz(&v) {
        None
    } else {
        Some(v.len() - 1)
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

fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);
    let mut mono = [0u64; 2]; // H_μ ≪ H_{μ'} when μ ⊆ μ' (μ_i ≤ μ'_i), μ ≠ μ'
    let mut mono_rev = [0u64; 2]; // H_{μ'} ≪ H_μ (reverse direction)
    let mut mono_cover = [0u64; 2]; // μ' covers μ (differ in exactly one part by 1)
    let mut mono_cover_rev = [0u64; 2];
    // Also: H_{μ'} ≪ t·H_μ (Wagner-like for covers)
    let mut wagner_cover = [0u64; 2];

    for n in 2..=max_n {
        for board in &gen_boards(n) {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            let m = board[0] as usize;
            if n <= 1 {
                continue;
            }
            let ideal = bruhat_lower_ideal(&perm);
            let subs = sub_partitions(board);
            let polys: Vec<Vec<i64>> = subs.iter().map(|mu| hit_poly(&ideal, mu)).collect();
            // Test all pairs (μ, μ') with μ ⊆ μ'
            for i in 0..subs.len() {
                for j in 0..subs.len() {
                    if i == j {
                        continue;
                    }
                    // Check μ_i ⊆ μ_j (subs[i]_k ≤ subs[j]_k for all k)
                    let contained = (0..n).all(|k| subs[i][k] <= subs[j][k]);
                    if !contained {
                        continue;
                    }
                    let f = &polys[i];
                    let g = &polys[j];
                    if pz(f) || pz(g) {
                        continue;
                    }
                    // H_μ ≪ H_{μ'} (smaller sub-partition on the left)
                    mono[0] += 1;
                    if !interlaces(f, g) {
                        mono[1] += 1;
                    }
                    // H_{μ'} ≪ H_μ (reverse)
                    mono_rev[0] += 1;
                    if !interlaces(g, f) {
                        mono_rev[1] += 1;
                    }
                    // Check if cover (differ in exactly one part by 1)
                    let diff: Vec<usize> = (0..n).filter(|&k| subs[i][k] != subs[j][k]).collect();
                    if diff.len() == 1 && subs[j][diff[0]] == subs[i][diff[0]] + 1 {
                        mono_cover[0] += 1;
                        if !interlaces(f, g) {
                            mono_cover[1] += 1;
                        }
                        mono_cover_rev[0] += 1;
                        if !interlaces(g, f) {
                            mono_cover_rev[1] += 1;
                        }
                        // Wagner: H_{μ'} ≪ t·H_μ?
                        let tf = pmt(f);
                        wagner_cover[0] += 1;
                        if !interlaces(g, &tf) {
                            wagner_cover[1] += 1;
                        }
                    }
                }
            }
        }
    }
    println!("=== Monotonicity in sub-partition (n ≤ {}) ===", max_n);
    let show = |name: &str, c: [u64; 2]| {
        if c[0] == 0 {
            println!("  {}: (no data)", name);
        } else if c[1] == 0 {
            println!("  {}: {}/{} ALL PASS <<<", name, c[0], c[0]);
        } else {
            println!("  {}: {}/{} pass ({} FAIL)", name, c[0] - c[1], c[0], c[1]);
        }
    };
    show("H_μ ≪ H_{μ'} (μ ⊆ μ', all pairs)", mono);
    show("H_{μ'} ≪ H_μ (reverse, all pairs)", mono_rev);
    show("H_μ ≪ H_{μ'} (covers only)", mono_cover);
    show("H_{μ'} ≪ H_μ (covers reverse)", mono_cover_rev);
    show("H_{μ'} ≪ t·H_μ (Wagner, covers)", wagner_cover);
}
