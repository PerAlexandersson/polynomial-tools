//! Test: at roots of C_j, does R'(r_i) alternate in sign?
//! Also test: is |R_1(r_i)| > |r_i * R_2(r_i)| (dominance)?

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

fn eval_poly(p: &[i64], t: f64) -> f64 {
    p.iter()
        .enumerate()
        .map(|(i, &c)| c as f64 * t.powi(i as i32))
        .sum()
}

fn find_roots(p: &[i64]) -> Vec<f64> {
    let d = p.len() - 1;
    if d == 0 {
        return vec![];
    }
    if d == 1 {
        return vec![-(p[0] as f64) / (p[1] as f64)];
    }
    // Companion matrix + QR
    let lc = p[d] as f64;
    let mut comp = vec![vec![0.0f64; d]; d];
    for i in 0..d - 1 {
        comp[i + 1][i] = 1.0;
    }
    for i in 0..d {
        comp[i][d - 1] = -(p[i] as f64) / lc;
    }
    let mut eigs = qr_eigs(&mut comp);
    eigs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    eigs
}

fn qr_eigs(h: &mut Vec<Vec<f64>>) -> Vec<f64> {
    let n = h.len();
    let mut eigs = Vec::new();
    let mut sz = n;
    for _ in 0..10000 {
        if sz <= 1 {
            if sz == 1 {
                eigs.push(h[0][0]);
            }
            break;
        }
        if h[sz - 1][sz - 2].abs()
            < 1e-12 * (h[sz - 1][sz - 1].abs() + h[sz - 2][sz - 2].abs() + 1e-30)
        {
            eigs.push(h[sz - 1][sz - 1]);
            sz -= 1;
            continue;
        }
        let a = h[sz - 2][sz - 2];
        let b = h[sz - 2][sz - 1];
        let c = h[sz - 1][sz - 2];
        let d = h[sz - 1][sz - 1];
        let disc = ((a - d) * (a - d) + 4.0 * b * c).max(0.0).sqrt();
        let (e1, e2) = ((a + d + disc) / 2.0, (a + d - disc) / 2.0);
        let mu = if (e1 - d).abs() < (e2 - d).abs() {
            e1
        } else {
            e2
        };
        for i in 0..sz {
            h[i][i] -= mu;
        }
        let mut cs = Vec::new();
        for i in 0..sz - 1 {
            let (a, b) = (h[i][i], h[i + 1][i]);
            let r = (a * a + b * b).sqrt();
            if r < 1e-30 {
                cs.push((1.0, 0.0));
                continue;
            }
            let (c, s) = (a / r, b / r);
            cs.push((c, s));
            for j in 0..sz {
                let (t1, t2) = (h[i][j], h[i + 1][j]);
                h[i][j] = c * t1 + s * t2;
                h[i + 1][j] = -s * t1 + c * t2;
            }
        }
        for (i, &(c, s)) in cs.iter().enumerate() {
            for j in 0..sz {
                let (t1, t2) = (h[j][i], h[j][i + 1]);
                h[j][i] = c * t1 + s * t2;
                h[j][i + 1] = -s * t1 + c * t2;
            }
        }
        for i in 0..sz {
            h[i][i] += mu;
        }
    }
    while eigs.len() < n {
        eigs.push(h[eigs.len()][eigs.len()]);
    }
    eigs
}

fn trim(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.last() == Some(&0) {
        v.pop();
    }
    v
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

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let n = a.len().max(b.len());
    let mut r = vec![0i64; n];
    for i in 0..a.len() {
        r[i] += a[i];
    }
    for i in 0..b.len() {
        r[i] += b[i];
    }
    r
}

fn main() {
    println!("=== Sign alternation & dominance test ===\n");
    let mut total_covers = 0;
    let mut sign_ok = 0;
    let mut dominance_ok = 0;

    for n in 2..=5 {
        let boards = boards_312(n);
        for board in &boards {
            let m = *board.last().unwrap();
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
                for j in 0..n {
                    if mu[j] >= board[j] {
                        continue;
                    }
                    if j > 0 && mu[j - 1] <= mu[j] {
                        continue;
                    }

                    let mut mu_prime = mu.clone();
                    mu_prime[j] += 1;

                    // C_j: cofactor at (j, mu[j]+1)
                    let col_del = mu[j] + 1;
                    let board_del: Vec<usize> = (0..n)
                        .filter(|&i| i != j)
                        .map(|i| board[i] - if col_del <= board[i] { 1 } else { 0 })
                        .collect();
                    let mu_del: Vec<usize> = (0..n)
                        .filter(|&i| i != j)
                        .map(|i| mu_prime[i] - if col_del <= mu_prime[i] { 1 } else { 0 })
                        .collect();
                    let cj = hit_poly(&board_del, &mu_del);
                    let cj = trim(&cj);
                    if cj.len() <= 1 {
                        continue;
                    } // degree 0, trivial

                    // R1 = sum_{k <= mu[j]} C_{j,k} (non-hit cofactors)
                    // R2 = sum_{k > mu[j]+1} C_{j,k} (hit cofactors)
                    let mut r1 = vec![0i64; 1];
                    let mut r2 = vec![0i64; 1];
                    for k in 1..=board[j] {
                        if k == col_del {
                            continue;
                        }
                        let bd: Vec<usize> = (0..n)
                            .filter(|&i| i != j)
                            .map(|i| board[i] - if k <= board[i] { 1 } else { 0 })
                            .collect();
                        let md: Vec<usize> = (0..n)
                            .filter(|&i| i != j)
                            .map(|i| mu_prime[i] - if k <= mu_prime[i] { 1 } else { 0 })
                            .collect();
                        let ck = hit_poly(&bd, &md);
                        if k <= mu[j] {
                            r1 = poly_add(&r1, &ck);
                        } else {
                            r2 = poly_add(&r2, &ck);
                        }
                    }

                    // Find roots of C_j
                    let roots = find_roots(&cj);
                    let d = cj.len() - 1;

                    total_covers += 1;
                    let mut this_sign_ok = true;
                    let mut this_dom_ok = true;

                    for (idx, &ri) in roots.iter().enumerate() {
                        let r1_val = eval_poly(&r1, ri);
                        let r2_val = eval_poly(&r2, ri);
                        let rprime_val = r1_val + ri * r2_val;

                        // Expected sign: (-1)^(d - idx - 1)
                        // for f << g same degree
                        let expected_positive = (d - idx - 1) % 2 == 0;
                        let actual_positive = rprime_val > 0.0;

                        if (expected_positive && rprime_val < -1e-8)
                            || (!expected_positive && rprime_val > 1e-8)
                        {
                            this_sign_ok = false;
                        }

                        // Dominance: |R1(ri)| > |ri * R2(ri)|
                        if r1_val.abs() < (ri * r2_val).abs() - 1e-8 {
                            this_dom_ok = false;
                        }
                    }

                    if this_sign_ok {
                        sign_ok += 1;
                    }
                    if this_dom_ok {
                        dominance_ok += 1;
                    }
                }
            }
        }
        println!(
            "n={}: sign_alt {}/{}, dominance {}/{}",
            n, sign_ok, total_covers, dominance_ok, total_covers
        );
    }
    println!(
        "\nTotal: sign_alt {}/{}, dominance {}/{}",
        sign_ok, total_covers, dominance_ok, total_covers
    );
}
