//! Test: for the cofactor expansion along row j at a cover step,
//! is |deg(C_{j,col}) - deg(C_{j,k})| <= 1 for ALL k?
//! This is the degree-gap lemma needed for the proof.
//! Use brute force on Bruhat ideals.
//! Test for n <= 8.
use polynomial_tools::real_rootedness::check_weak_interlacing;
use std::collections::BTreeSet;

fn pt(p: &[i64]) -> Vec<i64> { let mut v = p.to_vec(); while v.len() > 1 && *v.last().unwrap() == 0 { v.pop(); } v }
fn pz(p: &[i64]) -> bool { p.iter().all(|&c| c == 0) }
fn pdeg(p: &[i64]) -> usize { let v = pt(p); if pz(&v) { 0 } else { v.len() - 1 } }
fn interlaces(f: &[i64], g: &[i64]) -> bool {
    let f = pt(f); let g = pt(g);
    if pz(&f) { return true; } if pz(&g) { return false; }
    check_weak_interlacing(&f, &g).unwrap_or(false)
}

fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> {
    let n = perm.len();
    let mut vis: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut q: BTreeSet<Vec<u8>> = BTreeSet::new();
    q.insert(perm.to_vec());
    while let Some(cur) = q.pop_last() {
        for i in 0..n { for j in i+1..n {
            if cur[i] > cur[j] { let mut c = cur.clone(); c.swap(i, j);
                if !vis.contains(&c) { q.insert(c); } }
        }}
        vis.insert(cur);
    }
    vis.into_iter().collect()
}
fn board_to_perm(b: &[u8]) -> Vec<u8> {
    let n = b.len(); let mut p = vec![0u8; n]; let mut u = vec![false; n+1];
    for i in 0..n { for c in (1..=(b[i] as usize).min(n)).rev() {
        if !u[c] { p[i] = c as u8; u[c] = true; break; }
    }} p
}
fn is_312_avoiding(perm: &[u8]) -> bool {
    let n = perm.len();
    for i in 0..n { for j in i+1..n { for k in j+1..n {
        if perm[k] < perm[i] && perm[i] < perm[j] { return false; }
    }}} true
}
fn gen_boards(n: usize) -> Vec<Vec<u8>> {
    let mut r = vec![]; let mut c = vec![];
    gb(n, n, 0, &mut c, &mut r); r
}
fn gb(n: usize, mx: usize, d: usize, c: &mut Vec<u8>, r: &mut Vec<Vec<u8>>) {
    if d == n { r.push(c.clone()); return; }
    for v in (d+1).max(if d > 0 { c[d-1] as usize } else { 1 })..=mx {
        c.push(v as u8); gb(n, mx, d+1, c, r); c.pop();
    }
}
fn sub_partitions(lambda: &[u8]) -> Vec<Vec<u8>> {
    let n = lambda.len(); let mut result = Vec::new(); let mut mu = vec![0u8; n];
    fn gen(lam: &[u8], mu: &mut Vec<u8>, pos: usize, mx: u8, res: &mut Vec<Vec<u8>>) {
        if pos == lam.len() { res.push(mu.clone()); return; }
        let u = lam[pos].min(mx);
        for v in 0..=u { mu[pos] = v; gen(lam, mu, pos+1, v, res); }
    }
    gen(lambda, &mut mu, 0, lambda[0], &mut result); result
}

/// Compute cofactor polynomial: restrict to sigma[row_j]=col,
/// count hits on remaining positions against mu_prime.
fn cofactor_poly(ideal: &[Vec<u8>], mu_prime: &[u8], row_j: usize, col: usize) -> Vec<i64> {
    let n = mu_prime.len();
    let mut p = vec![0i64];
    for sigma in ideal {
        if sigma[row_j] as usize != col { continue; }
        let mut hits = 0;
        for i in 0..n {
            if i == row_j { continue; }
            if sigma[i] as usize > mu_prime[i] as usize { hits += 1; }
        }
        while p.len() <= hits { p.push(0); }
        p[hits] += 1;
    }
    pt(&p)
}

fn main() {
    let max_n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(8);

    let mut deg_gap_tests = 0u64;
    let mut deg_gap_pass = 0u64; // |deg C_j - deg C_{j,k}| <= 1
    let mut deg_gap_fail = 0u64;
    let mut max_gap_seen = 0i64;
    let mut interl_tests = 0u64;
    let mut interl_pass = 0u64;

    for n in 2..=max_n {
        let boards = gen_boards(n);
        let mut board_count = 0;
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) { continue; }
            let m = *board.iter().max().unwrap() as usize;
            let ideal = bruhat_lower_ideal(&perm);
            board_count += 1;

            // Only test a subset of sub-partitions for n >= 7 to keep runtime manageable
            let subs = sub_partitions(board);
            let test_subs: Vec<&Vec<u8>> = if n <= 6 {
                subs.iter().collect()
            } else {
                // For n=7,8: test every sub-partition still (but fewer boards)
                subs.iter().collect()
            };

            for mu in &test_subs {
                for j in 0..n {
                    if mu[j] >= board[j] { continue; }
                    let mut mu_prime = (*mu).clone();
                    mu_prime[j] += 1;
                    if j > 0 && mu_prime[j] > mu_prime[j-1] { continue; }
                    let col_j = mu[j] as usize + 1; // the cover column
                    let c_j = cofactor_poly(&ideal, &mu_prime, j, col_j);
                    if pz(&c_j) { continue; }
                    let deg_cj = pdeg(&c_j);

                    for k in 1..=m {
                        if k == col_j { continue; }
                        if k > board[j] as usize { continue; }
                        let c_jk = cofactor_poly(&ideal, &mu_prime, j, k);
                        if pz(&c_jk) { continue; }
                        let deg_cjk = pdeg(&c_jk);
                        let gap = (deg_cj as i64 - deg_cjk as i64).abs();
                        deg_gap_tests += 1;
                        if gap <= 1 {
                            deg_gap_pass += 1;
                        } else {
                            deg_gap_fail += 1;
                            if gap > max_gap_seen {
                                max_gap_seen = gap;
                                eprintln!("GAP={} board={:?} mu={:?} j={} col_j={} k={} deg_cj={} deg_cjk={}",
                                    gap, board, mu, j, col_j, k, deg_cj, deg_cjk);
                            }
                        }
                        // Also test actual interlacing
                        interl_tests += 1;
                        let ok = if k <= mu[j] as usize {
                            interlaces(&c_j, &c_jk)
                        } else {
                            interlaces(&c_jk, &c_j)
                        };
                        if ok { interl_pass += 1; }
                    }
                }
            }
        }
        eprintln!("n={}: {} boards done, {} deg-gap tests so far", n, board_count, deg_gap_tests);
    }

    println!("=== Degree gap lemma (n <= {}) ===", max_n);
    println!("  Degree gap |deg C_j - deg C_{{j,k}}| <= 1: {}/{} ({} FAIL)",
        deg_gap_pass, deg_gap_tests, deg_gap_fail);
    println!("  Max gap seen: {}", max_gap_seen);
    println!("  Cofactor interlacing: {}/{} ({} FAIL)",
        interl_pass, interl_tests, interl_tests - interl_pass);
}
