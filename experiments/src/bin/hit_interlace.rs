//! Test interlacing conditions for restricted hit polynomials H_μ^λ(t).
//!
//! For each 312-avoiding Ferrers board λ with n ≤ max_n and each sub-partition μ ⊆ λ,
//! compute D_k, U_k, A_k (split by first entry k and descent/ascent at position 1)
//! where each permutation σ contributes t^{hits(σ)} with hits = #{i : σ_i > μ_i}.
//!
//! Test: DU, reversed DD, reversed AA, forward AW, D≪W, and more.
//!
//! Usage: cargo run --release --bin hit_interlace -- 6

use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly, is_real_rooted};
use std::collections::BTreeSet;

// ── Polynomial helpers (i64 coefficient vectors, ascending degree) ──

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
fn pdeg(p: &[i64]) -> Option<usize> {
    let v = pt(p);
    if pz(&v) {
        None
    } else {
        Some(v.len() - 1)
    }
}

/// Check f ≪ g (weak interlacing).
fn interlaces(f: &[i64], g: &[i64]) -> bool {
    let f = pt(f);
    let g = pt(g);
    if pz(&f) {
        return true;
    }
    if pz(&g) {
        return false;
    }
    check_weak_interlacing(&f, &g).unwrap_or(false)
}

fn show(name: &str, c: [u64; 2]) {
    if c[0] == 0 {
        println!("  {}: (no data)", name);
    } else if c[1] == 0 {
        println!("  {}: {}/{} ALL PASS", name, c[0], c[0]);
    } else {
        println!("  {}: {}/{} pass ({} FAIL)", name, c[0] - c[1], c[0], c[1]);
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

/// Generate all sub-partitions μ with μ_i ≤ λ_i and μ weakly decreasing.
fn sub_partitions(lambda: &[u8]) -> Vec<Vec<u8>> {
    let mut result = Vec::new();
    let mut mu = vec![0u8; lambda.len()];
    fn gen(lambda: &[u8], mu: &mut Vec<u8>, pos: usize, max_val: u8, result: &mut Vec<Vec<u8>>) {
        if pos == lambda.len() {
            result.push(mu.clone());
            return;
        }
        let upper = lambda[pos].min(max_val);
        for v in 0..=upper {
            mu[pos] = v;
            gen(lambda, mu, pos + 1, v, result);
        }
    }
    gen(lambda, &mut mu, 0, lambda[0], &mut result);
    result
}

struct FailInfo {
    board: Vec<u8>,
    mu: Vec<u8>,
    k1: usize,
    k2: usize,
    f: Vec<i64>,
    g: Vec<i64>,
}

fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);
    println!("================================================================");
    println!("  Hit polynomial interlacing conditions");
    println!("  312-avoiding Ferrers boards, n <= {}", max_n);
    println!("================================================================\n");

    // Condition counters: [total, fails]
    let mut c_rr = [0u64; 2]; // real-rootedness of D_k, U_k, A_k
    let mut c_du = [0u64; 2]; // D_k ≪ U_j for all k, j
    let mut c_dd_rev = [0u64; 2]; // D_k ≪ D_{k'} for k > k' (reversed DD)
    let mut c_aa_rev = [0u64; 2]; // A_k ≪ A_{k'} for k > k' (reversed AA)
    let mut c_aw_fwd = [0u64; 2]; // A_j ≪ W_j' for j ≤ j' (forward AW)
    let mut c_aw_all = [0u64; 2]; // A_j ≪ W_j' for all j, j'
    let mut c_dw = [0u64; 2]; // D_k ≪ W_j for all k, j
    let mut c_uu_fwd = [0u64; 2]; // U_j ≪ U_{j'} for j ≤ j' (forward UU)
    let mut c_ww_fwd = [0u64; 2]; // W_j ≪ W_{j'} for j ≤ j' (forward WW)

    // Additional conditions: partial sums
    let mut c_dplus_rr = [0u64; 2]; // D_k^+ = Σ_{j<k} A_j real-rooted
    let mut c_uplus_rr = [0u64; 2]; // U_k^+ = Σ_{j≥k} W_j real-rooted

    // Failure details
    let mut du_fails: Vec<FailInfo> = Vec::new();
    let mut dd_fails: Vec<FailInfo> = Vec::new();
    let mut aa_fails: Vec<FailInfo> = Vec::new();
    let mut aw_all_fails: Vec<FailInfo> = Vec::new();
    let mut dw_fails: Vec<FailInfo> = Vec::new();
    let mut uu_fails: Vec<FailInfo> = Vec::new();

    let mut total_pairs = 0u64;

    for n in 2..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }

            let m = board[0] as usize; // max column value
            let ideal = bruhat_lower_ideal(&perm);

            // For each sub-partition μ ⊆ λ
            let subs = sub_partitions(board);
            for mu in &subs {
                total_pairs += 1;

                // Compute D_k and U_k for each first entry k.
                // D_k: σ_1 = k, σ_1 > σ_2, contributes t^{hits}
                // U_k: σ_1 = k, σ_1 < σ_2, contributes t^{hits}
                // hits = #{i : σ_i > μ_i}
                let mut d_k: Vec<Vec<i64>> = vec![vec![0i64]; m + 1]; // indexed 1..m
                let mut u_k: Vec<Vec<i64>> = vec![vec![0i64]; m + 1];

                for sigma in &ideal {
                    let k = sigma[0] as usize;
                    if k == 0 || k > m {
                        continue;
                    }

                    // Count hits
                    let hits: usize = (0..n)
                        .filter(|&i| {
                            let mu_i = if i < mu.len() { mu[i] as usize } else { 0 };
                            (sigma[i] as usize) > mu_i
                        })
                        .count();

                    let poly = if n >= 2 && sigma[0] > sigma[1] {
                        &mut d_k[k]
                    } else {
                        &mut u_k[k]
                    };
                    while poly.len() <= hits {
                        poly.push(0);
                    }
                    poly[hits] += 1;
                }

                // Trim
                for k in 1..=m {
                    d_k[k] = pt(&d_k[k]);
                    u_k[k] = pt(&u_k[k]);
                }

                // Compute A_k = D_k + U_k
                let mut a_k: Vec<Vec<i64>> = vec![vec![0i64]; m + 1];
                for k in 1..=m {
                    a_k[k] = pa(&d_k[k], &u_k[k]);
                }

                // Compute W_k = t·D_k + U_k (the natural analog of the peak case)
                let mut w_k: Vec<Vec<i64>> = vec![vec![0i64]; m + 1];
                for k in 1..=m {
                    w_k[k] = pa(&pmt(&d_k[k]), &u_k[k]);
                }

                // Collect nonzero indices
                let active: Vec<usize> = (1..=m).filter(|&k| !pz(&a_k[k])).collect();

                // Test real-rootedness of each D_k, U_k, A_k
                for &k in &active {
                    if !pz(&d_k[k]) {
                        c_rr[0] += 1;
                        if !is_real_rooted(&d_k[k]) && d_k[k].len() > 2 {
                            c_rr[1] += 1;
                        }
                    }
                    if !pz(&u_k[k]) {
                        c_rr[0] += 1;
                        if !is_real_rooted(&u_k[k]) && u_k[k].len() > 2 {
                            c_rr[1] += 1;
                        }
                    }
                    if !pz(&a_k[k]) {
                        c_rr[0] += 1;
                        if !is_real_rooted(&a_k[k]) && a_k[k].len() > 2 {
                            c_rr[1] += 1;
                        }
                    }
                }

                // DU: D_k ≪ U_j for all k, j
                for &k in &active {
                    if pz(&d_k[k]) {
                        continue;
                    }
                    for &j in &active {
                        if pz(&u_k[j]) {
                            continue;
                        }
                        c_du[0] += 1;
                        if !interlaces(&d_k[k], &u_k[j]) {
                            c_du[1] += 1;
                            if du_fails.len() < 5 {
                                du_fails.push(FailInfo {
                                    board: board.clone(),
                                    mu: mu.clone(),
                                    k1: k,
                                    k2: j,
                                    f: d_k[k].clone(),
                                    g: u_k[j].clone(),
                                });
                            }
                        }
                    }
                }

                // Reversed DD: D_k ≪ D_{k'} for k > k'
                for &k in &active {
                    if pz(&d_k[k]) {
                        continue;
                    }
                    for &kp in &active {
                        if kp >= k {
                            continue;
                        }
                        if pz(&d_k[kp]) {
                            continue;
                        }
                        c_dd_rev[0] += 1;
                        if !interlaces(&d_k[k], &d_k[kp]) {
                            c_dd_rev[1] += 1;
                            if dd_fails.len() < 5 {
                                dd_fails.push(FailInfo {
                                    board: board.clone(),
                                    mu: mu.clone(),
                                    k1: k,
                                    k2: kp,
                                    f: d_k[k].clone(),
                                    g: d_k[kp].clone(),
                                });
                            }
                        }
                    }
                }

                // Reversed AA: A_k ≪ A_{k'} for k > k'
                for &k in &active {
                    for &kp in &active {
                        if kp >= k {
                            continue;
                        }
                        c_aa_rev[0] += 1;
                        if !interlaces(&a_k[k], &a_k[kp]) {
                            c_aa_rev[1] += 1;
                            if aa_fails.len() < 5 {
                                aa_fails.push(FailInfo {
                                    board: board.clone(),
                                    mu: mu.clone(),
                                    k1: k,
                                    k2: kp,
                                    f: a_k[k].clone(),
                                    g: a_k[kp].clone(),
                                });
                            }
                        }
                    }
                }

                // AW forward: A_j ≪ W_{j'} for j ≤ j'
                for &j in &active {
                    for &jp in &active {
                        if j > jp {
                            continue;
                        }
                        if pz(&w_k[jp]) {
                            continue;
                        }
                        c_aw_fwd[0] += 1;
                        if !interlaces(&a_k[j], &w_k[jp]) {
                            c_aw_fwd[1] += 1;
                        }
                    }
                }

                // AW all: A_j ≪ W_{j'} for all j, j'
                for &j in &active {
                    for &jp in &active {
                        if pz(&w_k[jp]) {
                            continue;
                        }
                        c_aw_all[0] += 1;
                        if !interlaces(&a_k[j], &w_k[jp]) {
                            c_aw_all[1] += 1;
                            if aw_all_fails.len() < 5 {
                                aw_all_fails.push(FailInfo {
                                    board: board.clone(),
                                    mu: mu.clone(),
                                    k1: j,
                                    k2: jp,
                                    f: a_k[j].clone(),
                                    g: w_k[jp].clone(),
                                });
                            }
                        }
                    }
                }

                // DW: D_k ≪ W_j for all k, j
                for &k in &active {
                    if pz(&d_k[k]) {
                        continue;
                    }
                    for &j in &active {
                        if pz(&w_k[j]) {
                            continue;
                        }
                        c_dw[0] += 1;
                        if !interlaces(&d_k[k], &w_k[j]) {
                            c_dw[1] += 1;
                            if dw_fails.len() < 5 {
                                dw_fails.push(FailInfo {
                                    board: board.clone(),
                                    mu: mu.clone(),
                                    k1: k,
                                    k2: j,
                                    f: d_k[k].clone(),
                                    g: w_k[j].clone(),
                                });
                            }
                        }
                    }
                }

                // Forward UU: U_j ≪ U_{j'} for j ≤ j'
                for &j in &active {
                    if pz(&u_k[j]) {
                        continue;
                    }
                    for &jp in &active {
                        if jp < j {
                            continue;
                        }
                        if pz(&u_k[jp]) {
                            continue;
                        }
                        if j == jp {
                            continue;
                        }
                        c_uu_fwd[0] += 1;
                        if !interlaces(&u_k[j], &u_k[jp]) {
                            c_uu_fwd[1] += 1;
                            if uu_fails.len() < 5 {
                                uu_fails.push(FailInfo {
                                    board: board.clone(),
                                    mu: mu.clone(),
                                    k1: j,
                                    k2: jp,
                                    f: u_k[j].clone(),
                                    g: u_k[jp].clone(),
                                });
                            }
                        }
                    }
                }

                // Forward WW: W_j ≪ W_{j'} for j ≤ j'
                for &j in &active {
                    if pz(&w_k[j]) {
                        continue;
                    }
                    for &jp in &active {
                        if jp < j {
                            continue;
                        }
                        if pz(&w_k[jp]) {
                            continue;
                        }
                        if j == jp {
                            continue;
                        }
                        c_ww_fwd[0] += 1;
                        if !interlaces(&w_k[j], &w_k[jp]) {
                            c_ww_fwd[1] += 1;
                        }
                    }
                }

                // Partial sums: D_k^+ = Σ_{j<k} A_j, U_k^+ = Σ_{j≥k} W_j
                for &k in &active {
                    // D_k^+
                    let mut dplus = vec![0i64];
                    for &j in &active {
                        if j < k {
                            dplus = pa(&dplus, &a_k[j]);
                        }
                    }
                    if !pz(&dplus) {
                        c_dplus_rr[0] += 1;
                        if !is_real_rooted(&dplus) && dplus.len() > 2 {
                            c_dplus_rr[1] += 1;
                        }
                    }

                    // U_k^+
                    let mut uplus = vec![0i64];
                    for &j in &active {
                        if j >= k {
                            uplus = pa(&uplus, &w_k[j]);
                        }
                    }
                    if !pz(&uplus) {
                        c_uplus_rr[0] += 1;
                        if !is_real_rooted(&uplus) && uplus.len() > 2 {
                            c_uplus_rr[1] += 1;
                        }
                    }
                }
            }
        }
    }

    println!("Total (λ, μ) pairs tested: {}\n", total_pairs);

    println!("=== Real-rootedness ===");
    show("D_k, U_k, A_k (individual)", c_rr);
    show("D_k^+ = Σ_{j<k} A_j", c_dplus_rr);
    show("U_k^+ = Σ_{j≥k} W_j", c_uplus_rr);

    println!("\n=== Within-group interlacing ===");
    show("DU:  D_k ≪ U_j  (all k,j)", c_du);
    show("revDD: D_k ≪ D_{k'} (k > k')", c_dd_rev);
    show("revAA: A_k ≪ A_{k'} (k > k')", c_aa_rev);
    show("fwdAW: A_j ≪ W_{j'} (j ≤ j')", c_aw_fwd);
    show("allAW: A_j ≪ W_{j'} (all j,j')", c_aw_all);
    show("DW:  D_k ≪ W_j  (all k,j)", c_dw);
    show("fwdUU: U_j ≪ U_{j'} (j ≤ j')", c_uu_fwd);
    show("fwdWW: W_j ≪ W_{j'} (j ≤ j')", c_ww_fwd);

    // Print failure details
    if !du_fails.is_empty() {
        println!("\n  DU failures (first {}):", du_fails.len());
        for f in &du_fails {
            println!(
                "    λ={:?} μ={:?}: D_{} ≪ U_{} FAILS",
                f.board, f.mu, f.k1, f.k2
            );
            println!("      D = {}", format_poly(&f.f));
            println!("      U = {}", format_poly(&f.g));
        }
    }
    if !dd_fails.is_empty() {
        println!("\n  revDD failures (first {}):", dd_fails.len());
        for f in &dd_fails {
            println!(
                "    λ={:?} μ={:?}: D_{} ≪ D_{} FAILS",
                f.board, f.mu, f.k1, f.k2
            );
            println!("      f = {}", format_poly(&f.f));
            println!("      g = {}", format_poly(&f.g));
        }
    }
    if !aa_fails.is_empty() {
        println!("\n  revAA failures (first {}):", aa_fails.len());
        for f in &aa_fails {
            println!(
                "    λ={:?} μ={:?}: A_{} ≪ A_{} FAILS",
                f.board, f.mu, f.k1, f.k2
            );
            println!("      f = {}", format_poly(&f.f));
            println!("      g = {}", format_poly(&f.g));
        }
    }
    if !aw_all_fails.is_empty() {
        println!("\n  allAW failures (first {}):", aw_all_fails.len());
        for f in &aw_all_fails {
            println!(
                "    λ={:?} μ={:?}: A_{} ≪ W_{} FAILS",
                f.board, f.mu, f.k1, f.k2
            );
            println!("      A = {}", format_poly(&f.f));
            println!("      W = {}", format_poly(&f.g));
        }
    }
    if !dw_fails.is_empty() {
        println!("\n  DW failures (first {}):", dw_fails.len());
        for f in &dw_fails {
            println!(
                "    λ={:?} μ={:?}: D_{} ≪ W_{} FAILS",
                f.board, f.mu, f.k1, f.k2
            );
            println!("      D = {}", format_poly(&f.f));
            println!("      W = {}", format_poly(&f.g));
        }
    }
    if !uu_fails.is_empty() {
        println!("\n  fwdUU failures (first {}):", uu_fails.len());
        for f in &uu_fails {
            println!(
                "    λ={:?} μ={:?}: U_{} ≪ U_{} FAILS",
                f.board, f.mu, f.k1, f.k2
            );
            println!("      f = {}", format_poly(&f.f));
            println!("      g = {}", format_poly(&f.g));
        }
    }

    println!("\nDone.");
}
