//! Test the (t-1) Lemma conditions for the P_{k,l} refinement.
//!
//! The (t-1) Lemma says: if h << f and f(0) >= h(0), then f << f + (t-1)h.
//!
//! For the peak proof, this reduces interlacing P_k << P_j (fixed last entry l)
//! to showing h = Σ_{i=k}^{j-1} D_{i,l} << P_k.
//!
//! This binary checks:
//! 1. Reversed DD for fixed last entry: D_{i,l} << D_{i',l} for i > i'
//! 2. DU for fixed last entry: D_{i,l} << U_{j,l} for all i,j
//! 3. h << P_k for each relevant (k,j,l) triple
//! 4. The (t-1) Lemma conclusion: P_k << P_j
//! 5. Whether all conditions needed for the proof hold
//!
//! Usage: cargo run --release --bin peak_t1_lemma -- 7

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
fn pf(p: &[i64]) -> String {
    let p = pt(p);
    if pz(&p) {
        return "0".into();
    }
    format_poly(&p)
}

/// Check f << g (handles same-degree via Wagner: f << g iff g << tf).
fn interlaces(f: &[i64], g: &[i64]) -> bool {
    let f = pt(f);
    let g = pt(g);
    if pz(&f) {
        return true;
    }
    if pz(&g) {
        return false;
    }
    match check_weak_interlacing(&f, &g) {
        Some(true) => true,
        Some(false) => false,
        None => {
            let df = f.len() - 1;
            let dg = g.len() - 1;
            if df == dg {
                // Same degree: f << g iff g << tf (Wagner)
                let tf = pmt(&f);
                check_weak_interlacing(&g, &tf) == Some(true)
            } else {
                false
            }
        }
    }
}

// ── Combinatorial infrastructure ──

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

fn peaks(w: &[u8]) -> usize {
    if w.len() < 3 {
        return 0;
    }
    (1..w.len() - 1)
        .filter(|&i| w[i - 1] < w[i] && w[i] > w[i + 1])
        .count()
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
    println!("================================================================");
    println!("  (t-1) LEMMA VERIFICATION for P_{{k,l}} peak polynomials");
    println!("  312-avoiding Ferrers boards, n <= {}", max_n);
    println!("================================================================\n");

    // Counters: [tested, failed]
    let mut dd_rev_fixed_l = [0u64; 2]; // D_{i,l} << D_{i',l} for i > i', fixed l
    let mut du_fixed_l = [0u64; 2]; // D_{i,l} << U_{j,l} for all i,j, fixed l
    let mut uu_fixed_l = [0u64; 2]; // U_{j,l} << U_{j',l} for j < j', fixed l
    let mut h_ll_pk = [0u64; 2]; // h = Σ D_{i,l} << P_k for each (k,j,l)
    let mut pk_ll_pj = [0u64; 2]; // P_k << P_j (the conclusion)
    let mut di_ll_pk = [0u64; 2]; // D_{i,l} << P_k for each i >= k
    let mut f0_ge_h0 = [0u64; 2]; // f(0) >= h(0) check
    let mut pkl_rr = [0u64; 2]; // P_{k,l} real-rooted

    let mut valid_boards = 0u64;

    for n in 1..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            valid_boards += 1;
            let m = (board[0] as usize).min(n);
            let ideal = bruhat_lower_ideal(&perm);

            // Compute D_{k,l} and U_{k,l} for all k,l in 1..=m
            // d_kl[k][l] and u_kl[k][l]
            let mut d_kl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            let mut u_kl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];

            for pi in &ideal {
                let k = pi[0] as usize;
                if k > m {
                    continue;
                }
                let last = *pi.last().unwrap() as usize;
                if last > m {
                    continue;
                }
                let pk = peaks(pi);

                let is_desc = pi.len() >= 2 && pi[0] > pi[1];
                let poly = if is_desc {
                    &mut d_kl[k][last]
                } else {
                    &mut u_kl[k][last]
                };
                while poly.len() <= pk {
                    poly.push(0);
                }
                poly[pk] += 1;
            }

            // Compute A_{k,l} = D_{k,l} + U_{k,l} and W_{k,l} = t*D_{k,l} + U_{k,l}
            let mut a_kl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            let mut w_kl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            for k in 1..=m {
                for l in 1..=m {
                    a_kl[k][l] = pa(&d_kl[k][l], &u_kl[k][l]);
                    w_kl[k][l] = pa(&pmt(&d_kl[k][l]), &u_kl[k][l]);
                }
            }

            // Check P_{k,l} real-rootedness
            for k in 1..=m {
                for l in 1..=m {
                    if !pz(&a_kl[k][l]) {
                        pkl_rr[0] += 1;
                        if !is_real_rooted(&a_kl[k][l]) {
                            pkl_rr[1] += 1;
                            println!(
                                "FAIL P_{{k,l}} real-rooted: board={:?} k={} l={} P={}",
                                board,
                                k,
                                l,
                                pf(&a_kl[k][l])
                            );
                        }
                    }
                }
            }

            // Test 1: Reversed DD for fixed last entry
            for l in 1..=m {
                for i in 1..=m {
                    for ip in 1..i {
                        // i > ip, check D_{i,l} << D_{ip,l}
                        if !pz(&d_kl[i][l]) && !pz(&d_kl[ip][l]) {
                            dd_rev_fixed_l[0] += 1;
                            if !interlaces(&d_kl[i][l], &d_kl[ip][l]) {
                                dd_rev_fixed_l[1] += 1;
                                if dd_rev_fixed_l[1] <= 5 {
                                    println!(
                                        "FAIL DD_rev(i={},ip={},l={}): {} << {} board={:?}",
                                        i,
                                        ip,
                                        l,
                                        pf(&d_kl[i][l]),
                                        pf(&d_kl[ip][l]),
                                        board
                                    );
                                }
                            }
                        }
                    }
                }
            }

            // Test 2: DU for fixed last entry
            for l in 1..=m {
                for i in 1..=m {
                    for j in 1..=m {
                        if !pz(&d_kl[i][l]) && !pz(&u_kl[j][l]) {
                            du_fixed_l[0] += 1;
                            if !interlaces(&d_kl[i][l], &u_kl[j][l]) {
                                du_fixed_l[1] += 1;
                                if du_fixed_l[1] <= 5 {
                                    println!(
                                        "FAIL DU(i={},j={},l={}): {} << {} board={:?}",
                                        i,
                                        j,
                                        l,
                                        pf(&d_kl[i][l]),
                                        pf(&u_kl[j][l]),
                                        board
                                    );
                                }
                            }
                        }
                    }
                }
            }

            // Test 2b: UU for fixed last entry
            for l in 1..=m {
                for j in 1..=m {
                    for jp in j + 1..=m {
                        if !pz(&u_kl[j][l]) && !pz(&u_kl[jp][l]) {
                            uu_fixed_l[0] += 1;
                            if !interlaces(&u_kl[j][l], &u_kl[jp][l]) {
                                uu_fixed_l[1] += 1;
                                if uu_fixed_l[1] <= 5 {
                                    println!(
                                        "FAIL UU(j={},jp={},l={}): {} << {} board={:?}",
                                        j,
                                        jp,
                                        l,
                                        pf(&u_kl[j][l]),
                                        pf(&u_kl[jp][l]),
                                        board
                                    );
                                }
                            }
                        }
                    }
                }
            }

            // Compute P_k for fixed l using the Q-recursion:
            // P_k^{(l)} = Σ_{j<k} D_{j,l} + t·Σ_{j>=k} D_{j,l} + Σ_j U_{j,l}
            // (This is P_{k,l} at the NEXT level λ+, but we can also just check
            //  within-column interlacing at the CURRENT level using A_{k,l} = P_{k,l})

            // Test 3: D_{i,l} << P_{k,l} = A_{k,l} for i >= k (at current level)
            for l in 1..=m {
                for k in 1..=m {
                    if pz(&a_kl[k][l]) {
                        continue;
                    }
                    for i in k..=m {
                        if i == l {
                            continue;
                        } // D_{l,l} = 0
                        if !pz(&d_kl[i][l]) {
                            di_ll_pk[0] += 1;
                            if !interlaces(&d_kl[i][l], &a_kl[k][l]) {
                                di_ll_pk[1] += 1;
                                if di_ll_pk[1] <= 5 {
                                    println!("FAIL D_{{i,l}}<<P_{{k,l}}: i={} k={} l={} D={} P={} board={:?}",
                                        i, k, l, pf(&d_kl[i][l]), pf(&a_kl[k][l]), board);
                                }
                            }
                        }
                    }
                }
            }

            // Test 4: h = Σ_{i=k}^{j-1} D_{i,l} << P_{k,l} = A_{k,l}
            // and Test 5: P_k << P_j (the conclusion)
            for l in 1..=m {
                for k in 1..=m {
                    if pz(&a_kl[k][l]) {
                        continue;
                    }
                    for j in k + 1..=m {
                        if j == l {
                            continue;
                        } // A_{l,l} = 0
                        if pz(&a_kl[j][l]) {
                            continue;
                        }

                        // h = Σ_{i=k}^{j-1} D_{i,l}
                        let mut h = vec![0i64];
                        for i in k..j {
                            if i == l {
                                continue;
                            } // D_{l,l} = 0
                            h = pa(&h, &d_kl[i][l]);
                        }

                        if !pz(&h) {
                            // Check h << P_k
                            h_ll_pk[0] += 1;
                            if !interlaces(&h, &a_kl[k][l]) {
                                h_ll_pk[1] += 1;
                                if h_ll_pk[1] <= 5 {
                                    println!(
                                        "FAIL h<<P_k: k={} j={} l={} h={} P_k={} board={:?}",
                                        k,
                                        j,
                                        l,
                                        pf(&h),
                                        pf(&a_kl[k][l]),
                                        board
                                    );
                                }
                            }

                            // Check f(0) >= h(0)
                            let f0 = a_kl[k][l][0];
                            let h0 = h[0];
                            f0_ge_h0[0] += 1;
                            if f0 < h0 {
                                f0_ge_h0[1] += 1;
                                if f0_ge_h0[1] <= 5 {
                                    println!("FAIL f(0)>=h(0): k={} j={} l={} f(0)={} h(0)={} board={:?}",
                                        k, j, l, f0, h0, board);
                                }
                            }
                        }

                        // Check P_k << P_j (the goal)
                        pk_ll_pj[0] += 1;
                        if !interlaces(&a_kl[k][l], &a_kl[j][l]) {
                            pk_ll_pj[1] += 1;
                            if pk_ll_pj[1] <= 5 {
                                println!(
                                    "FAIL P_k<<P_j: k={} j={} l={} P_k={} P_j={} board={:?}",
                                    k,
                                    j,
                                    l,
                                    pf(&a_kl[k][l]),
                                    pf(&a_kl[j][l]),
                                    board
                                );
                            }
                        }
                    }
                }
            }
        }
        println!(
            "n={}: {} valid boards (cumulative: {})",
            n,
            boards
                .iter()
                .filter(|b| is_312_avoiding(&board_to_perm(b)))
                .count(),
            valid_boards
        );
    }

    println!("\n================================================================");
    println!("  RESULTS (n <= {})", max_n);
    println!("================================================================\n");

    let show = |name: &str, c: [u64; 2]| {
        if c[0] == 0 {
            println!("  {}: (no tests)", name);
        } else {
            let status = if c[1] == 0 { "ALL PASS" } else { "FAILURES" };
            println!("  {}: {}/{} pass {}", name, c[0] - c[1], c[0], status);
        }
    };

    println!("=== Conditions at level lambda (fixed last entry l) ===");
    show("P_{k,l} real-rooted", pkl_rr);
    show(
        "Reversed DD: D_{i,l} << D_{i',l} for i > i'",
        dd_rev_fixed_l,
    );
    show("DU: D_{i,l} << U_{j,l} for all i,j", du_fixed_l);
    show("UU: U_{j,l} << U_{j',l} for j < j'", uu_fixed_l);

    println!("\n=== (t-1) Lemma conditions ===");
    show("D_{i,l} << P_{k,l} for i >= k", di_ll_pk);
    show("h = Σ D_{i,l} << P_{k,l}", h_ll_pk);
    show("f(0) >= h(0)", f0_ge_h0);
    show("CONCLUSION: P_k << P_j (reversed, fixed l)", pk_ll_pj);

    println!("\n=== Proof status ===");
    if dd_rev_fixed_l[1] == 0
        && du_fixed_l[1] == 0
        && h_ll_pk[1] == 0
        && f0_ge_h0[1] == 0
        && pk_ll_pj[1] == 0
    {
        println!("  ALL CONDITIONS HOLD! The (t-1) Lemma proof strategy is viable.");
        println!("  Proof chain: reversed DD + DU + cone => h << P_k => (t-1) Lemma => P_k << P_j");
    } else {
        println!("  Some conditions fail. Details above.");
    }
}
