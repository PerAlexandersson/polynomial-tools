//! Test: does R ≪ Σ_W imply R ≪ P, given Σ_W = P + (t-1)Σ_D?
//! More generally: if R ≪ f+(t-1)h and h ≪ f, does R ≪ f follow?
//! ("Reverse shift lemma")

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
fn poly_scale(a: &[i64], c: i64) -> Vec<i64> { a.iter().map(|&x| x*c).collect() }
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
    println!("=== Test: R ≪ Σ_W (provable) vs R ≪ P (target) ===\n");
    
    // Also test: does Σ_D ≪ P hold? (needed for shift lemma P ≪ Σ_W)
    let (mut t1, mut p1) = (0, 0); // Σ_D ≪ P
    let (mut t2, mut p2) = (0, 0); // R ≪ Σ_W
    let (mut t3, mut p3) = (0, 0); // R ≪ P

    for n in 2..=7 {
        let boards = boards_312(n);
        for board in &boards {
            let m = *board.last().unwrap();
            let (dp, up) = compute_du(board);

            let mut sigma_d = vec![]; let mut sigma_u = vec![];
            let mut p_lam = vec![];
            for j in 1..=m {
                sigma_d = poly_add(&sigma_d, &dp[j]);
                sigma_u = poly_add(&sigma_u, &up[j]);
                p_lam = poly_add(&p_lam, &dp[j]);
                p_lam = poly_add(&p_lam, &up[j]);
            }
            // Σ_W = tΣ_D + Σ_U
            let sigma_w = poly_add(&poly_tmul(&sigma_d), &sigma_u);

            // Σ_D ≪ P
            let sd = trim(&sigma_d); let pl = trim(&p_lam);
            if !sd.is_empty() && !pl.is_empty() {
                t1 += 1;
                if interlaces_weak(&sd, &pl) { p1 += 1; }
            }

            for m_prime in [m, m+1] {
                let mut r_mp = vec![];
                for j in 1..=m {
                    let c = std::cmp::min(j, m_prime) as i64;
                    r_mp = poly_add(&r_mp, &poly_scale(&dp[j], c));
                }
                let r_mp = trim(&r_mp);
                if r_mp.is_empty() { continue; }

                // R ≪ Σ_W
                let sw = trim(&sigma_w);
                t2 += 1;
                if interlaces_weak(&r_mp, &sw) { p2 += 1; }

                // R ≪ P
                t3 += 1;
                if interlaces_weak(&r_mp, &pl) { p3 += 1; }
            }
        }
        println!("n={}: Σ_D≪P {}/{}, R≪Σ_W {}/{}, R≪P {}/{}", 
                 n, p1, t1, p2, t2, p3, t3);
    }
    println!("\nΣ_D≪P: {}/{}", p1, t1);
    println!("R≪Σ_W: {}/{}", p2, t2);
    println!("R≪P:   {}/{}", p3, t3);
}
