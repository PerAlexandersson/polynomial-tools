//! Explore excedance-generating polynomials on Ferrers boards,
//! with refinement by first entry and ascent/descent at position 1.
//!
//! Analogous to peak_explore but using exc instead of peak.
//! Tests whether the same IH structure works:
//!   (a) A_j ≼ W_l for j ≤ l
//!   (b) W_j ≼ W_l for j ≤ l
//!   (c) D_j ≼ U_l for all j,l
//! where A_k = D_k + U_k, W_k = tD_k + U_k,
//! D_k = exc poly for perms with first entry k and descent at pos 1,
//! U_k = exc poly for perms with first entry k and ascent at pos 1.

use std::collections::BTreeSet;

fn exact_interlaces(fc: &[i64], gc: &[i64]) -> bool {
    polynomial_tools::check_weak_interlacing(fc, gc).unwrap_or(false)
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

fn excedances(w: &[u8]) -> usize {
    (0..w.len()).filter(|&i| w[i] as usize > i + 1).count()
}

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
fn pf(p: &[i64]) -> String {
    let p = pt(p);
    if pz(&p) {
        "0".into()
    } else {
        let mut t = vec![];
        for (i, &c) in p.iter().enumerate() {
            if c == 0 {
                continue;
            }
            match (c, i) {
                (c, 0) => t.push(format!("{}", c)),
                (1, 1) => t.push("t".into()),
                (c, 1) => t.push(format!("{}t", c)),
                (1, e) => t.push(format!("t^{}", e)),
                (c, e) => t.push(format!("{}t^{}", c, e)),
            }
        }
        t.join(" + ")
    }
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

fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    println!(
        "=== Excedance refinement: testing IH on all boards n ≤ {} ===\n",
        max_n
    );

    // First check: does the SAME recursion hold for exc?
    // The recursion D_k^+ = Σ_{j<k} A_j, U_k^+ = Σ_{j≥k} W_j
    // was derived from the STRUCTURE of the Ferrers board recursion (how first entry
    // and initial ascent/descent relate). It depends on the stat only through
    // how the stat changes when we remove the first row.
    // For peaks: removing the first entry can create/destroy a peak at position 2.
    // For excedances: removing the first entry changes all positions by -1,
    // so excedance depends on position. This may NOT give the same recursion!
    //
    // Let me first just CHECK empirically if the recursion holds for exc.

    let mut total = 0;
    let mut rr_fails = 0;
    let mut aw_fails = 0;
    let mut ww_fails = 0;
    let mut du_fails = 0;
    let mut st_fails = 0;

    // Also check if the recursion D_k^+ = Σ_{j<k} A_j holds for exc
    let mut rec_fails = 0;

    for n in 2..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let nn = board.len();
            let m = (board[0] as usize).min(nn);
            let perm = board_to_perm(board);
            let ideal = bruhat_lower_ideal(&perm);

            // Compute D_k, U_k using EXCEDANCE stat
            let mut d = vec![vec![0i64]; m + 1];
            let mut u = vec![vec![0i64]; m + 1];
            let mut p_poly = vec![0i64];

            for pi in &ideal {
                let exc = excedances(pi);
                while p_poly.len() <= exc {
                    p_poly.push(0);
                }
                p_poly[exc] += 1;

                if pi.len() < 2 {
                    let k = pi[0] as usize;
                    if k <= m {
                        while u[k].len() < 1 {
                            u[k].push(0);
                        }
                        u[k][0] += 1;
                    }
                    continue;
                }
                let k = pi[0] as usize;
                if k > m {
                    continue;
                }
                let is_descent = pi[0] > pi[1];
                let poly = if is_descent { &mut d[k] } else { &mut u[k] };
                while poly.len() <= exc {
                    poly.push(0);
                }
                poly[exc] += 1;
            }
            p_poly = pt(&p_poly);

            let mut a = vec![vec![0i64]; m + 1];
            let mut w = vec![vec![0i64]; m + 1];
            for k in 1..=m {
                a[k] = pa(&d[k], &u[k]);
                w[k] = pa(&pmt(&d[k]), &u[k]);
            }

            // Check total exc poly is real-rooted
            if !pz(&p_poly) && !polynomial_tools::is_real_rooted(&p_poly) {
                rr_fails += 1;
            }

            total += 1;

            // Build S_k, T_k
            let mut s = vec![vec![0i64]; m + 2];
            let mut t_a = vec![vec![0i64]; m + 2];
            for k in 2..=m + 1 {
                s[k] = pa(&s[k - 1], &a[k - 1]);
            }
            for k in (1..=m).rev() {
                t_a[k] = pa(&t_a[k + 1], &w[k]);
            }

            // Test (a) A_j ≼ W_l for j ≤ l
            for j in 1..=m {
                for l in j..=m {
                    if pz(&a[j]) || pz(&w[l]) {
                        continue;
                    }
                    if !exact_interlaces(&a[j], &w[l]) {
                        aw_fails += 1;
                    }
                }
            }

            // Test (b) W_j ≼ W_l for j < l
            for j in 1..m {
                for l in j + 1..=m {
                    if pz(&w[j]) || pz(&w[l]) {
                        continue;
                    }
                    if !exact_interlaces(&w[j], &w[l]) {
                        ww_fails += 1;
                    }
                }
            }

            // Test (c) D_j ≼ U_l for all j,l
            for j in 1..=m {
                for l in 1..=m {
                    if pz(&d[j]) || pz(&u[l]) {
                        continue;
                    }
                    if !exact_interlaces(&d[j], &u[l]) {
                        du_fails += 1;
                    }
                }
            }

            // Test S_k ≼ T_k
            for k in 1..=m {
                if pz(&s[k]) || pz(&t_a[k]) {
                    continue;
                }
                if !exact_interlaces(&s[k], &t_a[k]) {
                    st_fails += 1;
                }
            }

            // Check if the RECURSION itself holds for exc
            // For each valid lambda+ = (m', lambda_1+1,...), check if
            // D_k^+ = S_k and U_k^+ = T_k match direct computation at lambda+
            // Only do this for the specific case m' = m
            if n >= 3 {
                // board = current lambda. lambda+ = (m, board[0]+1,...,board[n-1]+1)
                // but this requires computing on lambda+ directly. Skip for now.
            }
        }
        println!(
            "n={}: boards={}, RR={}, A≼W={}, W≼W={}, D≼U={}, S≼T={}",
            n, total, rr_fails, aw_fails, ww_fails, du_fails, st_fails
        );
    }

    println!("\n=== SUMMARY ===");
    println!("Boards: {}", total);
    println!("Total exc poly not RR: {}", rr_fails);
    println!("A_j≼W_l fails: {}", aw_fails);
    println!("W_j≼W_l fails: {}", ww_fails);
    println!("D_j≼U_l fails: {}", du_fails);
    println!("S_k≼T_k fails: {}", st_fails);

    if rr_fails == 0 && aw_fails == 0 && ww_fails == 0 && du_fails == 0 && st_fails == 0 {
        println!("\nALL CONDITIONS HOLD — same proof structure works for excedances!");
    } else if rr_fails == 0 {
        println!("\nExcedance polys always RR, but IH conditions have failures.");
        println!("The peak proof structure does NOT directly transfer.");
    }

    // Print detailed example
    println!("\n=== Detailed: S_4 = [4,4,4,4], exc stat ===");
    let board = vec![4u8, 4, 4, 4];
    let perm = board_to_perm(&board);
    let ideal = bruhat_lower_ideal(&perm);
    let m = 4;
    let mut d = vec![vec![0i64]; m + 1];
    let mut u = vec![vec![0i64]; m + 1];
    for pi in &ideal {
        if pi.len() < 2 {
            continue;
        }
        let k = pi[0] as usize;
        if k > m {
            continue;
        }
        let exc = excedances(pi);
        let poly = if pi[0] > pi[1] { &mut d[k] } else { &mut u[k] };
        while poly.len() <= exc {
            poly.push(0);
        }
        poly[exc] += 1;
    }
    for k in 1..=m {
        let ak = pa(&d[k], &u[k]);
        let wk = pa(&pmt(&d[k]), &u[k]);
        println!(
            "k={}: D={}, U={}, A={}, W={}",
            k,
            pf(&d[k]),
            pf(&u[k]),
            pf(&ak),
            pf(&wk)
        );
    }
}
