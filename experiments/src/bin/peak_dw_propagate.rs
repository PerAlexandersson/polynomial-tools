//! Test: does D_j ≪ W_{j'} for all j,j' propagate?
//! At λ+: D_k^+ = S_k, W_l^+ = tS_l + T_l.
//! Need: S_k ≪ tS_l + T_l for all k, l.
//! This is exactly AW restricted to the D-part: D_k^+ ≪ W_l^+.

use polynomial_tools::real_rootedness::{check_weak_interlacing, is_real_rooted};

fn peaks(w: &[u8]) -> usize {
    if w.len() < 3 { return 0; }
    (1..w.len()-1).filter(|&i| w[i-1] < w[i] && w[i] > w[i+1]).count()
}
fn all_perms(n: u8) -> Vec<Vec<u8>> {
    if n <= 1 { return vec![(1..=n).collect()]; }
    let mut r = Vec::new();
    for p in all_perms(n-1) { for i in 0..=p.len() { let mut q=p.clone(); q.insert(i,n); r.push(q); } }
    r
}
fn ferrers_perms(board: &[usize]) -> Vec<Vec<u8>> {
    let n = board.len();
    all_perms(n as u8).into_iter().filter(|p| (0..n).all(|i| (p[i] as usize) <= board[i])).collect()
}
fn compute_du(board: &[usize]) -> (Vec<Vec<i64>>, Vec<Vec<i64>>) {
    let perms = ferrers_perms(board);
    let n = board.len(); let m = *board.last().unwrap();
    let mut d = vec![vec![]; m+1]; let mut u = vec![vec![]; m+1];
    for p in &perms {
        if n < 2 { continue; }
        let k = p[0] as usize; let pk = peaks(p);
        let poly = if p[0] > p[1] { &mut d[k] } else { &mut u[k] };
        while poly.len() <= pk { poly.push(0); } poly[pk] += 1;
    }
    (d, u)
}
fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let n = a.len().max(b.len()); let mut r = vec![0i64; n];
    for i in 0..a.len() { r[i] += a[i]; } for i in 0..b.len() { r[i] += b[i]; } r
}
fn poly_tmul(a: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; a.len() + 1];
    for i in 0..a.len() { r[i+1] = a[i]; } r
}
fn trim(p: &[i64]) -> Vec<i64> { let mut v=p.to_vec(); while v.last()==Some(&0) { v.pop(); } v }
fn deg(p: &[i64]) -> usize { let t = trim(p); if t.is_empty() { 0 } else { t.len()-1 } }
fn interlaces_weak(f: &[i64], g: &[i64]) -> bool {
    let f = trim(f); let g = trim(g);
    if f.is_empty() { return is_real_rooted(&g); }
    if g.is_empty() { return false; }
    let (df, dg) = (deg(&f), deg(&g));
    if dg == df + 1 { check_weak_interlacing(&f, &g) == Some(true) }
    else if dg == df { let tf = poly_tmul(&f); check_weak_interlacing(&g, &tf) == Some(true) }
    else { false }
}
fn boards_312(n: usize) -> Vec<Vec<usize>> {
    fn gen(n: usize, b: &mut Vec<usize>, r: &mut Vec<Vec<usize>>) {
        if b.len() == n { r.push(b.clone()); return; }
        let i = b.len(); let prev = b.last().copied().unwrap_or(i+1).max(i+1);
        for v in prev..=n { b.push(v); gen(n, b, r); b.pop(); }
    }
    let mut r = Vec::new(); let mut b = Vec::new(); gen(n, &mut b, &mut r); r
}

fn main() {
    println!("=== D≪W propagation: S_k ≪ W_l^+ for all k,l ===\n");
    let (mut total, mut pass) = (0usize, 0usize);
    
    for n in 2..=7 {
        let boards = boards_312(n);
        let (mut nt, mut np) = (0, 0);
        for board in &boards {
            let m = *board.last().unwrap();
            let (dp, up) = compute_du(board);
            let mut a = vec![vec![]; m+1];
            let mut w = vec![vec![]; m+1];
            for j in 1..=m {
                a[j] = poly_add(&dp[j], &up[j]);
                w[j] = poly_add(&poly_tmul(&dp[j]), &up[j]);
            }
            let m_prime = m + 1;
            let mut s = vec![vec![]; m_prime+1];
            let mut w_plus = vec![vec![]; m_prime+1];
            for k in 1..=m_prime {
                let mut sk = vec![]; let mut tk = vec![];
                for j in 1..=m {
                    if j < k { sk = poly_add(&sk, &a[j]); }
                    if j >= k { tk = poly_add(&tk, &w[j]); }
                }
                s[k] = sk.clone();
                w_plus[k] = poly_add(&poly_tmul(&sk), &tk);
            }
            for k in 1..=m_prime { for l in 1..=m_prime {
                let sk = trim(&s[k]); let wl = trim(&w_plus[l]);
                if sk.is_empty() || wl.is_empty() { continue; }
                nt += 1;
                if interlaces_weak(&sk, &wl) { np += 1; }
            }}
        }
        total += nt; pass += np;
        println!("n={}: S≪W+ {}/{}", n, np, nt);
    }
    println!("\nTotal: {}/{}", pass, total);
    
    // Also: can we PROVE S_k ≪ W_l^+ from AW + DU at λ?
    // S_k = Σ_{j<k} A_j. Each A_j ≪ W_{j'} (AW) for all j'.
    // W_l^+ = tS_l + T_l. By cone: A_j ≪ T_l (from AW, cone). ✓
    // And A_j ≪ tS_l? Wagner: S_l ≪ A_j. S_l = Σ_{j'<l} A_{j'}.
    // S_l ≪ A_j needs each A_{j'} ≪ A_j for j' < l. NOT available from AW.
    // BUT: S_k is a partial sum (not individual D_j), so maybe cone helps.
    // S_k ≪ W_l^+ = tS_l + T_l. By cone: S_k ≪ tS_l AND S_k ≪ T_l.
    // S_k ≪ T_l: DU at λ+ (Step 1). ✓
    // S_k ≪ tS_l: Wagner gives S_l ≪ S_k.
    //   For l ≤ k: S_l ⊂ S_k, S_l ≪ S_k by cone (add more terms). 
    //     Each extra A_j (for l ≤ j < k): S_l ≪ A_j? Same issue.
    //   For l > k: S_k ⊂ S_l, so S_l has MORE terms. S_l ≪ S_k?
    //     S_l = S_k + extra, so BIGGER ≪ SMALLER. This is reversed DD.
    println!("\nNote: S_k ≪ tS_l requires S_l ≪ S_k (Wagner).");
    println!("This needs reversed DD at λ+ (S≪S(rev)), which holds computationally.");
}
