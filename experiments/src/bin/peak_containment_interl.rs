//! Test the stronger containment interlacing:
//! For mu' ⊇ mu with |deg H_mu - deg H_mu'| <= 1:
//!   H_{mu'} ≪ H_mu
//! on 312-avoiding Ferrers boards.

use polynomial_tools::real_rootedness::check_weak_interlacing;

fn all_perms(n: u8) -> Vec<Vec<u8>> {
    if n <= 1 { return vec![(1..=n).collect()]; }
    let mut r = Vec::new();
    for p in all_perms(n-1) {
        for i in 0..=p.len() {
            let mut q = p.clone(); q.insert(i, n); r.push(q);
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

fn hit_poly(board: &[usize], mu: &[usize]) -> Vec<i64> {
    let n = board.len();
    let perms = ferrers_perms(board);
    let mut coeffs = vec![0i64; n + 1];
    for p in &perms {
        let hits = (0..n).filter(|&i| p[i] as usize > mu[i]).count();
        coeffs[hits] += 1;
    }
    while coeffs.last() == Some(&0) { coeffs.pop(); }
    coeffs
}

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

fn deg(p: &[i64]) -> i64 {
    let t = trim(p);
    if t.is_empty() { -1 } else { (t.len() - 1) as i64 }
}

fn interlaces_weak(f: &[i64], g: &[i64]) -> bool {
    let f = trim(f); let g = trim(g);
    if f.is_empty() { return true; }
    if g.is_empty() { return false; }
    let df = f.len() - 1; let dg = g.len() - 1;
    if dg == df + 1 {
        check_weak_interlacing(&f, &g) == Some(true)
    } else if dg == df {
        let tf = poly_tmul(&f);
        check_weak_interlacing(&g, &tf) == Some(true)
    } else { false }
}

fn boards_312(n: usize) -> Vec<Vec<usize>> {
    fn gen(n: usize, b: &mut Vec<usize>, r: &mut Vec<Vec<usize>>) {
        if b.len() == n { r.push(b.clone()); return; }
        let i = b.len();
        let prev = b.last().copied().unwrap_or(i+1).max(i+1);
        for v in prev..=n { b.push(v); gen(n, b, r); b.pop(); }
    }
    let mut r = Vec::new(); let mut b = Vec::new();
    gen(n, &mut b, &mut r); r
}

fn main() {
    println!("=== Containment interlacing (deg gap <= 1) ===\n");
    let mut total = 0;
    let mut pass = 0;
    let mut fail_examples = 0;

    for n in 2..=5 {
        let boards = boards_312(n);
        let (mut nt, mut np) = (0, 0);

        for board in &boards {
            let m = *board.last().unwrap();

            // Generate all valid sub-partitions
            fn gen_mu(board: &[usize], idx: usize,
                      prev: usize, mu: &mut Vec<usize>,
                      result: &mut Vec<Vec<usize>>) {
                if idx == board.len() {
                    result.push(mu.clone()); return;
                }
                let max_val = prev.min(board[idx]);
                for v in 0..=max_val {
                    mu.push(v);
                    gen_mu(board, idx+1, v, mu, result);
                    mu.pop();
                }
            }
            let mut all_mu = Vec::new();
            gen_mu(board, 0, m, &mut Vec::new(), &mut all_mu);

            // Test all pairs mu' ⊇ mu
            for mu in &all_mu {
                let h_mu = hit_poly(board, mu);
                if h_mu.is_empty() { continue; }

                for mu_prime in &all_mu {
                    // Check mu' >= mu componentwise
                    if (0..n).any(|i| mu_prime[i] < mu[i]) {
                        continue;
                    }
                    if mu_prime == mu { continue; }

                    let h_mp = hit_poly(board, mu_prime);
                    if h_mp.is_empty() { continue; }

                    let d_mu = deg(&h_mu);
                    let d_mp = deg(&h_mp);

                    // Only test if degree gap <= 1
                    if (d_mu - d_mp).abs() > 1 { continue; }

                    nt += 1;
                    // H_{mu'} ≪ H_mu (larger sub-partition
                    // has fewer hits, roots further left)
                    if interlaces_weak(&h_mp, &h_mu) {
                        np += 1;
                    } else {
                        fail_examples += 1;
                        if fail_examples <= 5 {
                            println!("  FAIL: board={:?}", board);
                            println!("    mu ={:?} deg={}",
                                mu, d_mu);
                            println!("    mu'={:?} deg={}",
                                mu_prime, d_mp);
                            println!("    H_mu ={:?}", h_mu);
                            println!("    H_mu'={:?}", h_mp);
                        }
                    }
                }
            }
        }
        total += nt; pass += np;
        println!("n={}: {}/{}", n, np, nt);
    }
    println!("\nTotal: {}/{}", pass, total);
}
