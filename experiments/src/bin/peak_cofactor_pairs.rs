//! Test C_j ≪ C_{j,k} for ONLY the specific pairs
//! arising in the proof of Theorem (thm:hit_rr).
//! These are: board λ', sub-partitions ν ⊇ ν_k,
//! where ν and ν_k come from the cofactor construction.

use polynomial_tools::real_rootedness::check_weak_interlacing;

fn all_perms(n: u8) -> Vec<Vec<u8>> {
    if n <= 1 {
        return vec![(1..=n).collect()];
    }
    let mut r = Vec::new();
    for p in all_perms(n - 1) {
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
    all_perms(n as u8)
        .into_iter()
        .filter(|p| (0..n).all(|i| (p[i] as usize) <= board[i]))
        .collect()
}

fn hit_poly(board: &[usize], mu: &[usize]) -> Vec<i64> {
    let n = board.len();
    let perms = ferrers_perms(board);
    let mut coeffs = vec![0i64; n + 2];
    for p in &perms {
        let hits = (0..n).filter(|&i| p[i] as usize > mu[i]).count();
        coeffs[hits] += 1;
    }
    while coeffs.last() == Some(&0) {
        coeffs.pop();
    }
    coeffs
}

fn poly_tmul(a: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; a.len() + 1];
    for i in 0..a.len() {
        r[i + 1] = a[i];
    }
    r
}

fn trim(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.last() == Some(&0) {
        v.pop();
    }
    v
}

fn deg(p: &[i64]) -> i64 {
    let t = trim(p);
    if t.is_empty() {
        -1
    } else {
        (t.len() - 1) as i64
    }
}

fn interlaces_weak(f: &[i64], g: &[i64]) -> bool {
    let f = trim(f);
    let g = trim(g);
    if f.is_empty() {
        return true;
    }
    if g.is_empty() {
        return false;
    }
    let df = f.len() - 1;
    let dg = g.len() - 1;
    if dg == df + 1 {
        check_weak_interlacing(&f, &g) == Some(true)
    } else if dg == df {
        let tf = poly_tmul(&f);
        check_weak_interlacing(&g, &tf) == Some(true)
    } else {
        false
    }
}

fn boards_312(n: usize) -> Vec<Vec<usize>> {
    fn gen(n: usize, b: &mut Vec<usize>, r: &mut Vec<Vec<usize>>) {
        if b.len() == n {
            r.push(b.clone());
            return;
        }
        let i = b.len();
        let prev = b.last().copied().unwrap_or(i + 1).max(i + 1);
        for v in prev..=n {
            b.push(v);
            gen(n, b, r);
            b.pop();
        }
    }
    let mut r = Vec::new();
    let mut b = Vec::new();
    gen(n, &mut b, &mut r);
    r
}

fn main() {
    println!("=== Cofactor-pair interlacing (proof-specific) ===\n");
    let mut total = 0;
    let mut pass = 0;

    for n in 2..=6 {
        let boards = boards_312(n);
        let (mut nt, mut np) = (0, 0);

        for board in &boards {
            let m = *board.last().unwrap();
            // Generate all valid sub-partitions
            fn gen_mu(
                board: &[usize],
                idx: usize,
                prev: usize,
                mu: &mut Vec<usize>,
                result: &mut Vec<Vec<usize>>,
            ) {
                if idx == board.len() {
                    result.push(mu.clone());
                    return;
                }
                let max_val = prev.min(board[idx]);
                for v in 0..=max_val {
                    mu.push(v);
                    gen_mu(board, idx + 1, v, mu, result);
                    mu.pop();
                }
            }
            let mut all_mu = Vec::new();
            gen_mu(board, 0, m, &mut Vec::new(), &mut all_mu);

            for mu in &all_mu {
                // For each valid cover j
                for j in 0..n {
                    if mu[j] >= board[j] {
                        continue;
                    }
                    if j > 0 && mu[j - 1] <= mu[j] {
                        continue;
                    }
                    let mut mu_prime = mu.clone();
                    mu_prime[j] += 1;

                    // Board λ' = board with row j deleted
                    let board_del: Vec<usize> = (0..n)
                        .filter(|&i| i != j)
                        .map(|i| board[i] - 1) // delete col mu[j]+1
                        .collect();
                    // Actually: delete row j, delete col (mu[j]+1)
                    // Board' rows: board[i] for i != j,
                    // each reduced by 1 if mu[j]+1 <= board[i]
                    let col_del = mu[j] + 1;
                    let board_prime: Vec<usize> = (0..n)
                        .filter(|&i| i != j)
                        .map(|i| board[i] - if col_del <= board[i] { 1 } else { 0 })
                        .collect();

                    // ν (for C_j = C_{j, mu[j]+1}):
                    let nu: Vec<usize> = (0..n)
                        .filter(|&i| i != j)
                        .map(|i| mu_prime[i] - if col_del <= mu_prime[i] { 1 } else { 0 })
                        .collect();

                    // For each k <= mu[j] (part (i)):
                    // ν_k has parts: mu'_i - [k <= mu'_i]
                    for k in 1..=mu[j] {
                        let nu_k: Vec<usize> = (0..n)
                            .filter(|&ii| ii != j)
                            .map(|ii| mu_prime[ii] - if k <= mu_prime[ii] { 1 } else { 0 })
                            .collect();

                        // Check nu >= nu_k componentwise
                        if (0..nu.len()).any(|i| nu[i] < nu_k[i]) {
                            continue;
                        }
                        let diff: usize = nu.iter().zip(nu_k.iter()).map(|(a, b)| a - b).sum();
                        if diff <= 1 {
                            continue;
                        } // single cover, OK

                        // Multi-cover: test H_nu ≪ H_{nu_k}
                        let h_nu = hit_poly(&board_prime, &nu);
                        let h_nk = hit_poly(&board_prime, &nu_k);
                        if h_nu.is_empty() || h_nk.is_empty() {
                            continue;
                        }

                        nt += 1;
                        if interlaces_weak(&h_nu, &h_nk) {
                            np += 1;
                        } else {
                            if nt - np <= 3 {
                                println!("  FAIL: board={:?} j={} k={}", board, j + 1, k);
                                println!("    nu={:?} nu_k={:?} diff={}", nu, nu_k, diff);
                                println!("    H_nu={:?} (deg {})", h_nu, deg(&h_nu));
                                println!("    H_nk={:?} (deg {})", h_nk, deg(&h_nk));
                            }
                        }
                    }

                    // For each k > mu[j]+1 (part (ii)):
                    for k in (mu[j] + 2)..=board[j] {
                        let nu_k: Vec<usize> = (0..n)
                            .filter(|&ii| ii != j)
                            .map(|ii| mu_prime[ii] - if k <= mu_prime[ii] { 1 } else { 0 })
                            .collect();

                        if (0..nu.len()).any(|i| nu_k[i] < nu[i]) {
                            continue;
                        }
                        let diff: usize = nu_k.iter().zip(nu.iter()).map(|(a, b)| a - b).sum();
                        if diff <= 1 {
                            continue;
                        }

                        let h_nu = hit_poly(&board_prime, &nu);
                        let h_nk = hit_poly(&board_prime, &nu_k);
                        if h_nu.is_empty() || h_nk.is_empty() {
                            continue;
                        }

                        nt += 1;
                        // Part (ii): C_{j,k} ≪ C_j, i.e. H_{nu_k} ≪ H_nu
                        if interlaces_weak(&h_nk, &h_nu) {
                            np += 1;
                        } else {
                            if nt - np <= 3 {
                                println!("  FAIL(ii): board={:?} j={} k={}", board, j + 1, k);
                                println!("    nu={:?} nu_k={:?} diff={}", nu, nu_k, diff);
                                println!("    H_nu={:?}", h_nu);
                                println!("    H_nk={:?}", h_nk);
                            }
                        }
                    }
                }
            }
        }
        total += nt;
        pass += np;
        println!("n={}: {}/{}", n, np, nt);
    }
    println!("\nTotal: {}/{}", pass, total);
}
