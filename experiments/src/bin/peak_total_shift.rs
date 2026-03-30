//! Test R_{m'} ≪ P^λ for the total peak recursion P^{λ+} = m'P^λ + (t-1)R_{m'}.
//! Uses exact Bézout interlacing from polynomial-tools.

use polynomial_tools::real_rootedness::{check_weak_interlacing, is_real_rooted};

fn peaks(w: &[u8]) -> usize {
    if w.len() < 3 { return 0; }
    (1..w.len()-1).filter(|&i| w[i-1] < w[i] && w[i] > w[i+1]).count()
}

fn all_perms(n: u8) -> Vec<Vec<u8>> {
    if n <= 1 { return vec![(1..=n).collect()]; }
    let mut r = Vec::new();
    for p in all_perms(n-1) {
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
    all_perms(n as u8).into_iter()
        .filter(|p| (0..n).all(|i| (p[i] as usize) <= board[i]))
        .collect()
}

fn compute_du(board: &[usize]) -> (Vec<Vec<i64>>, Vec<Vec<i64>>) {
    let perms = ferrers_perms(board);
    let n = board.len();
    let m = *board.last().unwrap();
    let mut d = vec![vec![]; m+1];
    let mut u = vec![vec![]; m+1];
    for p in &perms {
        if n < 2 { continue; }
        let k = p[0] as usize;
        let pk = peaks(p);
        let poly = if p[0] > p[1] { &mut d[k] } else { &mut u[k] };
        while poly.len() <= pk { poly.push(0); }
        poly[pk] += 1;
    }
    (d, u)
}

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let n = a.len().max(b.len());
    let mut r = vec![0i64; n];
    for i in 0..a.len() { r[i] += a[i]; }
    for i in 0..b.len() { r[i] += b[i]; }
    r
}

fn poly_scale(a: &[i64], c: i64) -> Vec<i64> {
    a.iter().map(|&x| x * c).collect()
}

/// Multiply by t (shift coefficients right, insert 0 at index 0)
fn poly_tmul(a: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; a.len() + 1];
    for i in 0..a.len() { r[i+1] = a[i]; }
    r
}

fn trim(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.last() == Some(&0) { v.pop(); }
    v
}

fn deg(p: &[i64]) -> usize {
    let t = trim(p);
    if t.is_empty() { 0 } else { t.len() - 1 }
}

/// Check f ≪ g (weakly) for polynomials with all nonpositive roots.
/// Same-degree case: f ≪ g ⟺ g ≪ tf (Wagner), reducing to degree diff 1.
fn interlaces_weak(f: &[i64], g: &[i64]) -> bool {
    let f = trim(f);
    let g = trim(g);
    if f.is_empty() { return is_real_rooted(&g); }
    if g.is_empty() { return false; }

    let df = deg(&f);
    let dg = deg(&g);

    if dg == df + 1 {
        // Direct: deg g = deg f + 1
        check_weak_interlacing(&f, &g) == Some(true)
    } else if dg == df {
        // Same degree: f ≪ g ⟺ g ≪ tf (Wagner's lemma)
        let tf = poly_tmul(&f);
        check_weak_interlacing(&g, &tf) == Some(true)
    } else {
        false
    }
}

fn boards_312(n: usize) -> Vec<Vec<usize>> {
    fn gen(n: usize, b: &mut Vec<usize>, r: &mut Vec<Vec<usize>>) {
        if b.len() == n { r.push(b.clone()); return; }
        let i = b.len();
        let prev = b.last().copied().unwrap_or(i+1).max(i+1);
        for v in prev..=n { b.push(v); gen(n, b, r); b.pop(); }
    }
    let mut r = Vec::new();
    let mut b = Vec::new();
    gen(n, &mut b, &mut r);
    r
}

fn main() {
    println!("=== R_m' ≪ P^λ (exact Bézout) ===\n");
    let mut total = 0;
    let mut pass = 0;

    for n in 2..=8 {
        let boards = boards_312(n);
        let (mut nt, mut np) = (0, 0);

        for board in &boards {
            let m = *board.last().unwrap();
            let (dp, up) = compute_du(board);

            // P^λ = Σ_j (D_j + U_j)
            let mut p_lam = vec![];
            for j in 1..=m {
                p_lam = poly_add(&p_lam, &dp[j]);
                p_lam = poly_add(&p_lam, &up[j]);
            }
            let p_lam = trim(&p_lam);
            if p_lam.is_empty() { continue; }

            for m_prime in [m, m + 1] {
                // R_{m'} = Σ_j min(j, m') · D_j
                let mut r_mp = vec![];
                for j in 1..=m {
                    let c = std::cmp::min(j, m_prime) as i64;
                    r_mp = poly_add(&r_mp, &poly_scale(&dp[j], c));
                }
                let r_mp = trim(&r_mp);
                if r_mp.is_empty() { continue; }

                nt += 1;
                if interlaces_weak(&r_mp, &p_lam) {
                    np += 1;
                } else {
                    println!("  FAIL: board={:?}, m'={}, R={:?}, P={:?}",
                             board, m_prime, &r_mp, &p_lam);
                }
            }
        }
        total += nt;
        pass += np;
        println!("n={}: R≪P {}/{} ({} boards)", n, np, nt, boards.len());
    }
    println!("\nTotal: {}/{}", pass, total);
}
