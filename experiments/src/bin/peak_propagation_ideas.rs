//! Explore ideas for propagating DU + reversed DD.
//!
//! The gap: need S_k ≪ S_{k'} for k > k' at λ+, where S_k = Σ_{j<k} A_j.
//! Idea 1: Maybe S_k ≪ S_{k'} follows from R ≪ P directly?
//!         Since R = Σ_s Σ_{D,≥s}, can we extract S_k ≪ S_{k'} from this?
//!
//! Idea 2: The shift lemma on partial sums.
//!         S_k = S_{k-1} + A_{k-1} = S_{k-1} + W_{k-1} - (t-1)D_{k-1}.
//!         So S_{k-1} = S_k - W_{k-1} + (t-1)D_{k-1} = S_k + (t-1)D_{k-1} - W_{k-1}.
//!         Not a clean (t-1) structure.
//!
//! Idea 3: Work with W-partial sums instead.
//!         T_k = Σ_{j≥k} W_j. T_{k-1} = T_k + W_{k-1}.
//!         T_{k-1} - T_k = W_{k-1} = A_{k-1} + (t-1)D_{k-1}.
//!         So T_k = T_{k-1} - A_{k-1} - (t-1)D_{k-1}.
//!         Still not clean.
//!
//! Idea 4: Maybe the AW condition A_j ≪ W_l at λ can be proved
//!         from R ≪ P + some additional lemma?
//!
//! Idea 5: Can we prove S_k ≪ T_k (DU at λ+) and S_k ≪ S_{k'} (rev DD at λ+)
//!         DIRECTLY from R ≪ P at λ, without going through individual conditions?
//!
//! Idea 6: Use the recursion P+ = m'P + (t-1)R to express S_k and T_l at λ+
//!         in terms of λ-level quantities, then apply R ≪ P.
//!
//! Idea 7: What conditions on A_j and W_j TOGETHER would self-propagate?
//!         Test: A_j ≪ W_l for j ≤ l ("forward AW"). Does this propagate?

use polynomial_tools::real_rootedness::{check_weak_interlacing, is_real_rooted};

fn peaks(w: &[u8]) -> usize {
    if w.len() < 3 {
        return 0;
    }
    (1..w.len() - 1)
        .filter(|&i| w[i - 1] < w[i] && w[i] > w[i + 1])
        .count()
}
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
fn compute_du(board: &[usize]) -> (Vec<Vec<i64>>, Vec<Vec<i64>>) {
    let perms = ferrers_perms(board);
    let n = board.len();
    let m = *board.last().unwrap();
    let mut d = vec![vec![]; m + 1];
    let mut u = vec![vec![]; m + 1];
    for p in &perms {
        if n < 2 {
            continue;
        }
        let k = p[0] as usize;
        let pk = peaks(p);
        let poly = if p[0] > p[1] { &mut d[k] } else { &mut u[k] };
        while poly.len() <= pk {
            poly.push(0);
        }
        poly[pk] += 1;
    }
    (d, u)
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
fn deg(p: &[i64]) -> usize {
    let t = trim(p);
    if t.is_empty() {
        0
    } else {
        t.len() - 1
    }
}

fn interlaces_weak(f: &[i64], g: &[i64]) -> bool {
    let f = trim(f);
    let g = trim(g);
    if f.is_empty() {
        return is_real_rooted(&g);
    }
    if g.is_empty() {
        return false;
    }
    let (df, dg) = (deg(&f), deg(&g));
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
    println!("=== Propagation ideas ===\n");

    // Test various conditions that might self-propagate
    let mut counts: std::collections::HashMap<&str, (usize, usize)> =
        std::collections::HashMap::new();
    let tests = [
        "D≪W(all)",
        "A≪W(j≤l)",
        "A≪W(all)",
        "U≪W(j≤l)",
        "U≪W(all)",
        "A≪tD(all)",
        "Σ_D≪U(all)",
        "R≪P",
        // At λ+ level:
        "S≪T(all)",
        "S≪S(rev)",
    ];
    for t in &tests {
        counts.insert(t, (0, 0));
    }

    for n in 2..=7 {
        for board in &boards_312(n) {
            let m = *board.last().unwrap();
            let (dp, up) = compute_du(board);

            // A_j = D_j + U_j, W_j = tD_j + U_j
            let mut a = vec![vec![]; m + 1];
            let mut w = vec![vec![]; m + 1];
            for j in 1..=m {
                a[j] = poly_add(&dp[j], &up[j]);
                w[j] = poly_add(&poly_tmul(&dp[j]), &up[j]);
            }

            // D≪W(all): D_j ≪ W_{j'} for all j, j'
            for j in 1..=m {
                for jp in 1..=m {
                    let dj = trim(&dp[j]);
                    let wjp = trim(&w[jp]);
                    if dj.is_empty() || wjp.is_empty() {
                        continue;
                    }
                    let e = counts.get_mut("D≪W(all)").unwrap();
                    e.0 += 1;
                    if interlaces_weak(&dj, &wjp) {
                        e.1 += 1;
                    }
                }
            }

            // A≪W(j≤l)
            for j in 1..=m {
                for l in j..=m {
                    let aj = trim(&a[j]);
                    let wl = trim(&w[l]);
                    if aj.is_empty() || wl.is_empty() {
                        continue;
                    }
                    let e = counts.get_mut("A≪W(j≤l)").unwrap();
                    e.0 += 1;
                    if interlaces_weak(&aj, &wl) {
                        e.1 += 1;
                    }
                }
            }

            // A≪W(all)
            for j in 1..=m {
                for l in 1..=m {
                    let aj = trim(&a[j]);
                    let wl = trim(&w[l]);
                    if aj.is_empty() || wl.is_empty() {
                        continue;
                    }
                    let e = counts.get_mut("A≪W(all)").unwrap();
                    e.0 += 1;
                    if interlaces_weak(&aj, &wl) {
                        e.1 += 1;
                    }
                }
            }

            // U≪W(j≤l)
            for j in 1..=m {
                for l in j..=m {
                    let uj = trim(&up[j]);
                    let wl = trim(&w[l]);
                    if uj.is_empty() || wl.is_empty() {
                        continue;
                    }
                    let e = counts.get_mut("U≪W(j≤l)").unwrap();
                    e.0 += 1;
                    if interlaces_weak(&uj, &wl) {
                        e.1 += 1;
                    }
                }
            }

            // U≪W(all)
            for j in 1..=m {
                for l in 1..=m {
                    let uj = trim(&up[j]);
                    let wl = trim(&w[l]);
                    if uj.is_empty() || wl.is_empty() {
                        continue;
                    }
                    let e = counts.get_mut("U≪W(all)").unwrap();
                    e.0 += 1;
                    if interlaces_weak(&uj, &wl) {
                        e.1 += 1;
                    }
                }
            }

            // A≪tD(all): A_j ≪ tD_{j'} for all j,j'
            for j in 1..=m {
                for jp in 1..=m {
                    let aj = trim(&a[j]);
                    let tdj = trim(&poly_tmul(&dp[jp]));
                    if aj.is_empty() || tdj.is_empty() {
                        continue;
                    }
                    let e = counts.get_mut("A≪tD(all)").unwrap();
                    e.0 += 1;
                    if interlaces_weak(&aj, &tdj) {
                        e.1 += 1;
                    }
                }
            }

            // At λ+ level: compute S_k and T_l
            let m_prime = m + 1;
            let mut s = vec![vec![]; m_prime + 1]; // S_k = Σ_{j<k} A_j
            let mut t_sum = vec![vec![]; m_prime + 1]; // T_k = Σ_{j≥k} W_j
            for k in 1..=m_prime {
                let mut sk = vec![];
                let mut tk = vec![];
                for j in 1..=m {
                    if j < k {
                        sk = poly_add(&sk, &a[j]);
                    }
                    if j >= k {
                        tk = poly_add(&tk, &w[j]);
                    }
                }
                s[k] = sk;
                t_sum[k] = tk;
            }

            // S≪T(all): S_k ≪ T_l for all k, l
            for k in 1..=m_prime {
                for l in 1..=m_prime {
                    let sk = trim(&s[k]);
                    let tl = trim(&t_sum[l]);
                    if sk.is_empty() || tl.is_empty() {
                        continue;
                    }
                    let e = counts.get_mut("S≪T(all)").unwrap();
                    e.0 += 1;
                    if interlaces_weak(&sk, &tl) {
                        e.1 += 1;
                    }
                }
            }

            // S≪S(rev): S_k ≪ S_{k'} for k > k'
            for k in 2..=m_prime {
                for kp in 1..k {
                    let sk = trim(&s[k]);
                    let skp = trim(&s[kp]);
                    if sk.is_empty() || skp.is_empty() {
                        continue;
                    }
                    let e = counts.get_mut("S≪S(rev)").unwrap();
                    e.0 += 1;
                    if interlaces_weak(&sk, &skp) {
                        e.1 += 1;
                    }
                }
            }
        }
        print!("n={}: ", n);
        let mut items: Vec<_> = counts.iter().collect();
        items.sort_by_key(|(k, _)| k.to_string());
        for (name, (total, pass)) in &items {
            if *total > 0 {
                let status = if pass == total { "✓" } else { "✗" };
                print!("{}={}/{}{} ", name, pass, total, status);
            }
        }
        println!();
    }

    println!("\n=== Summary ===");
    let mut items: Vec<_> = counts.iter().collect();
    items.sort_by_key(|(k, _)| k.to_string());
    for (name, (total, pass)) in &items {
        let fail = total - pass;
        let status = if fail == 0 {
            "ALL PASS".to_string()
        } else {
            format!("{} FAIL", fail)
        };
        println!("  {}: {}/{} ({})", name, pass, total, status);
    }
}
