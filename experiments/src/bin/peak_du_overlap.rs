//! Test proof strategies for D_{k'}^+ ≪ U_{j'}^+ in the overlapping case k' > j'.
//!
//! At the λ+ level, for a fixed last-column l':
//!   D_{k'}^+ = S_{k'}^{(col)}   (prefix sum of A's at col)
//!   U_{j'}^+ = T_{j'}^{(col)}   (suffix sum of W's at col)
//!   A_{k'}^+ = S_{k'} + T_{k'}
//!   W_{k'}^+ = t·S_{k'} + T_{k'}
//!
//! where col depends on the group (col = l' for k' > l', col = l'-1 for k' < l').
//!
//! The "overlapping" case is k' > j' in the SAME group (same col), i.e.,
//! both k' and j' map to the same underlying column.
//!
//! We test strategies A-D for proving D_{k'}^+ ≪ U_{j'}^+ when k' > j'.
//!
//! Usage: cargo run --release --bin peak_du_overlap -- 7

use polynomial_tools::real_rootedness::{check_weak_interlacing, is_real_rooted, format_poly};
use std::collections::BTreeSet;

// ── Polynomial helpers (i64 coefficient vectors, ascending degree) ──

fn pt(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.len() > 1 && *v.last().unwrap() == 0 { v.pop(); }
    v
}
fn pz(p: &[i64]) -> bool { p.iter().all(|&c| c == 0) }
fn pa(a: &[i64], b: &[i64]) -> Vec<i64> {
    let l = a.len().max(b.len());
    let mut r = vec![0i64; l];
    for (i, &v) in a.iter().enumerate() { r[i] += v; }
    for (i, &v) in b.iter().enumerate() { r[i] += v; }
    pt(&r)
}
fn pmt(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; p.len() + 1];
    for (i, &v) in p.iter().enumerate() { r[i + 1] = v; }
    pt(&r)
}
fn pdeg(p: &[i64]) -> Option<usize> {
    let v = pt(p);
    if pz(&v) { None } else { Some(v.len() - 1) }
}

/// Check f << g (weak interlacing).
fn interlaces(f: &[i64], g: &[i64]) -> bool {
    let f = pt(f);
    let g = pt(g);
    if pz(&f) { return true; }
    if pz(&g) { return false; }
    match check_weak_interlacing(&f, &g) {
        Some(true) => true,
        Some(false) => false,
        None => {
            let df = pdeg(&f);
            let dg = pdeg(&g);
            match (df, dg) {
                (Some(df), Some(dg)) if df == dg => {
                    let tf = pmt(&f);
                    match check_weak_interlacing(&g, &tf) {
                        Some(b) => b,
                        None => false,
                    }
                },
                _ => false,
            }
        }
    }
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
                    if !vis.contains(&c) { q.insert(c); }
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
            if !u[c] { p[i] = c as u8; u[c] = true; break; }
        }
    }
    p
}

fn is_312_avoiding(perm: &[u8]) -> bool {
    let n = perm.len();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                if perm[k] < perm[i] && perm[i] < perm[j] { return false; }
            }
        }
    }
    true
}

fn peaks(w: &[u8]) -> usize {
    if w.len() < 3 { return 0; }
    (1..w.len() - 1).filter(|&i| w[i - 1] < w[i] && w[i] > w[i + 1]).count()
}

fn gen_boards(n: usize) -> Vec<Vec<u8>> {
    let mut r = vec![];
    let mut c = vec![];
    gb(n, n, 0, &mut c, &mut r);
    r
}

fn gb(n: usize, mx: usize, d: usize, c: &mut Vec<u8>, r: &mut Vec<Vec<u8>>) {
    if d == n { r.push(c.clone()); return; }
    for v in (d + 1).max(if d > 0 { c[d - 1] as usize } else { 1 })..=mx {
        c.push(v as u8);
        gb(n, mx, d + 1, c, r);
        c.pop();
    }
}

struct FailDetail {
    board: Vec<u8>,
    col: usize,
    kp: usize,
    jp: usize,
    i_val: Option<usize>,  // for strategy C/D, the index i
    f_poly: Vec<i64>,
    g_poly: Vec<i64>,
}

fn main() {
    let max_n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(7);
    println!("================================================================");
    println!("  D_{{k'}}^+ ≪ U_{{j'}}^+ overlap (k' > j') proof strategies");
    println!("  312-avoiding Ferrers boards, n <= {}", max_n);
    println!("================================================================\n");

    // Strategy A: Induction on k' - j'
    // Base: D_{j'+1}^+ ≪ U_{j'}^+ (gap = 1, which is D_{k'}^+ ≪ U_{k'-1}^+)
    // Step: D_{k'-1}^+ ≪ U_{j'}^+ (IH) AND A_{k'-1} ≪ U_{j'}^+ => D_{k'}^+ ≪ U_{j'}^+
    //   because D_{k'}^+ = S_{k'} = S_{k'-1} + A_{k'-1} and we want S_{k'} ≪ T_{j'}
    //   by left-cone: if S_{k'-1} ≪ T_{j'} and A_{k'-1} ≪ T_{j'}, then S_{k'} ≪ T_{j'}
    // Test: A_{k'-1,col} ≪ U_{j'}^+ = T_{j'}^{(col)} for k' > j'
    let mut sa_base = [0u64; 2];       // base case: D_{j'+1}^+ ≪ U_{j'}^+ (gap=1)
    let mut sa_ak_interl_ujp = [0u64; 2]; // A_{k'-1,col} ≪ T_{j'}  (the key new condition)
    let mut sa_ak_interl_ujp_fails: Vec<FailDetail> = Vec::new();

    // Strategy B: Induction on j' decreasing
    // Base: D_{k'}^+ ≪ U_{k'-1}^+ (gap = 1, same as strat A base)
    // Step: D_{k'}^+ ≪ U_{j'+1}^+ (IH) AND D_{k'}^+ ≪ W_{j'} => D_{k'}^+ ≪ U_{j'}^+
    //   because U_{j'}^+ = T_{j'} = W_{j'} + T_{j'+1} = W_{j',col} + U_{j'+1}^+
    //   by right-cone: if D_{k'}^+ ≪ W_{j',col} and D_{k'}^+ ≪ U_{j'+1}^+, then D_{k'}^+ ≪ U_{j'}^+
    // Test: D_{k'}^+ = S_{k'} ≪ W_{j',col} for k' > j'
    let mut sb_dk_interl_wjp = [0u64; 2];
    let mut sb_dk_interl_wjp_fails: Vec<FailDetail> = Vec::new();

    // Strategy C: Mixed decomposition
    // U_{j'}^+ = T_{j'} = T_{k'} + sum_{j'<=i<k'} W_{i,col}
    // = U_{k'}^+ + sum_{j'<=i<k'} W_{i,col}
    // We already know D_{k'}^+ ≪ U_{k'}^+ (same-index DU).
    // If D_{k'}^+ ≪ W_{i,col} for each j' <= i < k', then by right-cone D_{k'}^+ ≪ U_{j'}^+
    // Test: S_{k'} ≪ W_{i,col} for i < k' (this is the sub-condition)
    let mut sc_same_du = [0u64; 2];      // D_{k'}^+ ≪ U_{k'}^+ (same-index, should pass)
    let mut sc_dk_interl_wi = [0u64; 2];  // S_{k'} ≪ W_{i,col} for i < k'
    let mut sc_dk_interl_wi_fails: Vec<FailDetail> = Vec::new();

    // Strategy D: Direct decomposition via P_{k'}^+
    // D_{k'}^+ + U_{j'}^+ = (D_{k'}^+ + U_{k'}^+) + sum_{j'<=i<k'} W_{i,col}
    //                      = P_{k'}^+ + sum_{j'<=i<k'} W_{i,col}
    // P_{k'}^+ = S_{k'} + T_{k'} = A_{k'}^+, which is real-rooted.
    // If sum_{j'<=i<k'} W_{i,col} ≪ P_{k'}^+, then by right-cone...
    // Actually this gives P_{k'}^+ ≪ D_{k'}^+ + U_{j'}^+, not what we want directly.
    // But if P_{k'}^+ is real-rooted and sum W_i ≪ P_{k'}^+, that means the total
    // D_{k'}^+ + U_{j'}^+ is real-rooted (by cone). That's a weaker statement.
    // Test: real-rootedness of P_{k'}^+ and sum W_i ≪ P_{k'}^+
    let mut sd_pk_rr = [0u64; 2];         // P_{k'}^+ real-rooted
    let mut sd_sumw_interl_pk = [0u64; 2]; // sum_{j'<=i<k'} W_{i,col} ≪ P_{k'}^+
    let mut sd_sumw_interl_pk_fails: Vec<FailDetail> = Vec::new();
    // Also check the reverse: P_{k'}^+ ≪ sum W_i (which would give sum W_i ≪ total by left-cone)
    let mut sd_pk_interl_sumw = [0u64; 2]; // P_{k'}^+ ≪ sum_{j'<=i<k'} W_{i,col}

    // Direct DU overlap verification (ground truth)
    let mut du_overlap = [0u64; 2];
    let mut du_overlap_fails: Vec<FailDetail> = Vec::new();

    let mut valid_boards = 0u64;

    for n in 1..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) { continue; }
            valid_boards += 1;

            let m = board[0] as usize;
            if n <= 1 { continue; }

            let ideal = bruhat_lower_ideal(&perm);

            // ── Compute D_{j,col} and U_{j,col} at lambda level ──
            let mut d_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            let mut u_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];

            for pi in &ideal {
                let j = pi[0] as usize;
                if j > m { continue; }
                let l = *pi.last().unwrap() as usize;
                if l > m { continue; }
                let pk = peaks(pi);
                let poly = if pi.len() >= 2 && pi[0] > pi[1] {
                    &mut d_jl[j][l]
                } else {
                    &mut u_jl[j][l]
                };
                while poly.len() <= pk { poly.push(0); }
                poly[pk] += 1;
            }

            // Compute A_{j,col} and W_{j,col}
            let mut a_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            let mut w_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            for j in 1..=m {
                for col in 1..=m {
                    a_jl[j][col] = pa(&d_jl[j][col], &u_jl[j][col]);
                    w_jl[j][col] = pa(&pmt(&d_jl[j][col]), &u_jl[j][col]);
                }
            }

            // ── Compute prefix/suffix sums for each column ──
            let mp = m + 1;
            // S_{k}^{(col)} = sum_{j=1}^{k-1} A_{j,col}  (D_{k'}^+ when col matches)
            // T_{k}^{(col)} = sum_{j=k}^{m} W_{j,col}    (U_{k'}^+ when col matches)
            let mut s_col: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; mp + 2]; m + 1];
            let mut t_col: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; mp + 2]; m + 1];

            for col in 1..=m {
                for k in 1..=m {
                    while s_col[col].len() <= k + 1 { s_col[col].push(vec![0i64]); }
                    s_col[col][k + 1] = pa(&s_col[col][k], &a_jl[k][col]);
                }
                while t_col[col].len() <= mp + 1 { t_col[col].push(vec![0i64]); }
                for k in (1..=m).rev() {
                    t_col[col][k] = pa(&t_col[col][k + 1], &w_jl[k][col]);
                }
            }

            // ── Now test strategies for each column and each pair k' > j' ──
            // At lambda+ level with last-column l', for k' > l':
            //   D_{k'}^+ = S_{k'}^{(l')}, U_{j'}^+ = T_{j'}^{(l')} if j' is in same group
            // Same group means both k', j' > l' (col = l') or both k', j' < l' (col = l'-1).
            //
            // For the overlap test, we need k' > j' in the same group.
            // Group "above l'": k' > l', col = l'. Valid k' range: l'+1..=m+1, mapping to S/T indices k'_idx in 1..=m
            //   Actually k' ranges over {1,...,m+1}\{l'}. For k' > l': col = l', and the index
            //   into S/T is k' - (if k' > l' then 1 else 0)...
            //
            // Let me be more careful. At lambda+ with n+1 positions:
            //   Positions 1..=m+1, last column l' in 1..=m+1.
            //   For k' != l': col = l' if k' > l', col = l'-1 if k' < l'.
            //   S_{k'}^{(col)} = sum_{j < k', j != l'} A_{...}
            //
            // Actually, from the existing code:
            //   For k' > l': col = l'
            //     D_{k'}^+ = S_{k'}^{(col)} where S uses indices up to m at col
            //     But S_{k'} means sum of A_{j,col} for j = 1..k'-1 (among indices that map to col)
            //
            // The existing code uses:
            //   let sk = if kp < s_col[col].len() { s_col[col][kp].clone() } ...
            //   let tk = if kp < t_col[col].len() { t_col[col][kp].clone() } ...
            // where kp is the original k' index and col is determined by the group.
            //
            // Actually looking more carefully at the existing code, it appears that S_{k'}^{(col)}
            // is just s_col[col][kp] when kp > lp (col = lp), and s_col[col][kp] when kp < lp (col = lp-1).
            // The index kp directly indexes into the prefix sum.
            //
            // For k' > l', col = l':
            //   S_{k'} = s_col[col][k'], which = sum_{j=1}^{k'-1} A_{j,col}
            //   T_{k'} = t_col[col][k'], which = sum_{j=k'}^{m} W_{j,col}
            //   But k' can go up to m+1, and t_col[col][m+1] = 0 (empty tail).
            //
            // For two indices k' > j' > l' in the same group (col = l'):
            //   D_{k'}^+ = S_{k'}^{(col)} = s_col[col][k']
            //   U_{j'}^+ = T_{j'}^{(col)} = t_col[col][j']
            //
            // Similarly for k' > j', both < l' (col = l'-1):
            //   D_{k'}^+ = s_col[col][k']
            //   U_{j'}^+ = t_col[col][j']

            // For simplicity: just iterate over all columns col from 1..=m,
            // and for each col, iterate over k' > j' in 1..=m (the range of valid
            // indices into s_col/t_col). These represent the "within-group overlap"
            // cases regardless of which l' produced them.
            for col in 1..=m {
                for kp in 2..=m {
                    let d_kp = &s_col[col][kp]; // D_{k'}^+
                    if pz(d_kp) { continue; }

                    for jp in 1..kp {
                        let u_jp = &t_col[col][jp]; // U_{j'}^+
                        if pz(u_jp) { continue; }

                        // Ground truth: D_{k'}^+ ≪ U_{j'}^+
                        du_overlap[0] += 1;
                        if !interlaces(d_kp, u_jp) {
                            du_overlap[1] += 1;
                            if du_overlap_fails.len() < 5 {
                                du_overlap_fails.push(FailDetail {
                                    board: board.clone(), col, kp, jp, i_val: None,
                                    f_poly: pt(d_kp), g_poly: pt(u_jp),
                                });
                            }
                        }

                        // ── Strategy A: A_{k'-1,col} ≪ T_{j'} ──
                        // D_{k'}^+ = S_{k'} = S_{k'-1} + A_{k'-1,col}
                        // If S_{k'-1} ≪ T_{j'} (IH for gap k'-1 - j' < k' - j')
                        // and A_{k'-1,col} ≪ T_{j'}, then by left-cone S_{k'} ≪ T_{j'}.
                        // Base case: gap=1, k'=j'+1: D_{j'+1}^+ ≪ U_{j'}^+
                        if kp == jp + 1 {
                            sa_base[0] += 1;
                            if !interlaces(d_kp, u_jp) {
                                sa_base[1] += 1;
                            }
                        }
                        // The sub-condition: A_{k'-1,col} ≪ T_{j'} for all k' > j'
                        // k'-1 >= j' >= 1, so k'-1 is in 1..=m
                        let km1 = kp - 1;
                        if km1 >= 1 && km1 <= m {
                            let a_km1 = &a_jl[km1][col];
                            if !pz(a_km1) {
                                sa_ak_interl_ujp[0] += 1;
                                if !interlaces(a_km1, u_jp) {
                                    sa_ak_interl_ujp[1] += 1;
                                    if sa_ak_interl_ujp_fails.len() < 5 {
                                        sa_ak_interl_ujp_fails.push(FailDetail {
                                            board: board.clone(), col, kp, jp, i_val: Some(km1),
                                            f_poly: pt(a_km1), g_poly: pt(u_jp),
                                        });
                                    }
                                }
                            }
                        }

                        // ── Strategy B: S_{k'} ≪ W_{j',col} ──
                        // U_{j'}^+ = T_{j'} = W_{j',col} + T_{j'+1}
                        // If S_{k'} ≪ T_{j'+1} (IH) and S_{k'} ≪ W_{j',col}, by right-cone done.
                        // Base: gap=1 (j'=k'-1): S_{k'} ≪ T_{k'-1} = W_{k'-1,col} + T_{k'}.
                        //   We know S_{k'} ≪ T_{k'} (same-index). So need S_{k'} ≪ W_{k'-1,col}.
                        if jp >= 1 && jp <= m {
                            let w_jp = &w_jl[jp][col];
                            if !pz(w_jp) {
                                sb_dk_interl_wjp[0] += 1;
                                if !interlaces(d_kp, w_jp) {
                                    sb_dk_interl_wjp[1] += 1;
                                    if sb_dk_interl_wjp_fails.len() < 5 {
                                        sb_dk_interl_wjp_fails.push(FailDetail {
                                            board: board.clone(), col, kp, jp, i_val: None,
                                            f_poly: pt(d_kp), g_poly: pt(w_jp),
                                        });
                                    }
                                }
                            }
                        }

                        // ── Strategy C: S_{k'} ≪ W_{i,col} for each j' <= i < k' ──
                        // Same-index check
                        let u_kp = &t_col[col][kp];
                        if !pz(u_kp) {
                            // Only count once per (col, kp) but that's tricky;
                            // count per (col, kp, jp) pair to keep it simple
                            if jp == 1 { // count same-index once per col,kp
                                sc_same_du[0] += 1;
                                if !interlaces(d_kp, u_kp) {
                                    sc_same_du[1] += 1;
                                }
                            }
                        }
                        // The sub-conditions: S_{k'} ≪ W_{i,col} for j' <= i < k'
                        for i in jp..kp {
                            if i >= 1 && i <= m {
                                let w_i = &w_jl[i][col];
                                if !pz(w_i) && !pz(d_kp) {
                                    sc_dk_interl_wi[0] += 1;
                                    if !interlaces(d_kp, w_i) {
                                        sc_dk_interl_wi[1] += 1;
                                        if sc_dk_interl_wi_fails.len() < 5 {
                                            sc_dk_interl_wi_fails.push(FailDetail {
                                                board: board.clone(), col, kp, jp: jp, i_val: Some(i),
                                                f_poly: pt(d_kp), g_poly: pt(w_i),
                                            });
                                        }
                                    }
                                }
                            }
                        }

                        // ── Strategy D: sum_{j'<=i<k'} W_{i,col} ≪ P_{k'}^+ ──
                        let p_kp = pa(d_kp, u_kp); // P_{k'}^+ = A_{k'}^+ = S_{k'} + T_{k'}
                        if !pz(&p_kp) {
                            sd_pk_rr[0] += 1;
                            if !is_real_rooted(&p_kp) {
                                sd_pk_rr[1] += 1;
                            }
                        }
                        // Compute sum_{j'<=i<k'} W_{i,col}
                        let mut sum_w = vec![0i64];
                        for i in jp..kp {
                            if i >= 1 && i <= m {
                                sum_w = pa(&sum_w, &w_jl[i][col]);
                            }
                        }
                        if !pz(&sum_w) && !pz(&p_kp) {
                            sd_sumw_interl_pk[0] += 1;
                            if !interlaces(&sum_w, &p_kp) {
                                sd_sumw_interl_pk[1] += 1;
                                if sd_sumw_interl_pk_fails.len() < 5 {
                                    sd_sumw_interl_pk_fails.push(FailDetail {
                                        board: board.clone(), col, kp, jp, i_val: None,
                                        f_poly: pt(&sum_w), g_poly: pt(&p_kp),
                                    });
                                }
                            }
                            sd_pk_interl_sumw[0] += 1;
                            if !interlaces(&p_kp, &sum_w) {
                                sd_pk_interl_sumw[1] += 1;
                            }
                        }
                    }
                }
            }
        }
        let nb = boards.iter().filter(|b| is_312_avoiding(&board_to_perm(b))).count();
        println!("n={}: {} boards (cumulative: {})", n, nb, valid_boards);
    }

    println!("\n================================================================");
    println!("  GROUND TRUTH: D_{{k'}}^+ ≪ U_{{j'}}^+ for k' > j' (overlap)");
    println!("================================================================");
    show("D_{k'}^+ ≪ U_{j'}^+", du_overlap);
    if !du_overlap_fails.is_empty() {
        println!("\n  --- First {} DU overlap failures ---", du_overlap_fails.len());
        for (i, f) in du_overlap_fails.iter().enumerate() {
            println!("  #{}: board={:?}, col={}, k'={}, j'={}",
                     i+1, f.board, f.col, f.kp, f.jp);
            println!("    D_{{k'}}^+ = {}", format_poly(&f.f_poly));
            println!("    U_{{j'}}^+ = {}", format_poly(&f.g_poly));
        }
    }

    println!("\n================================================================");
    println!("  STRATEGY A: Induction on k' - j' (left-cone on prefix sums)");
    println!("  Need: A_{{k'-1,col}} ≪ T_{{j'}} for k' > j'");
    println!("================================================================");
    show("Base (gap=1): D_{j'+1}^+ ≪ U_{j'}^+", sa_base);
    show("Key: A_{k'-1,col} ≪ U_{j'}^+", sa_ak_interl_ujp);
    if !sa_ak_interl_ujp_fails.is_empty() {
        println!("\n  --- First {} Strategy A failures ---", sa_ak_interl_ujp_fails.len());
        for (i, f) in sa_ak_interl_ujp_fails.iter().enumerate() {
            println!("  #{}: board={:?}, col={}, k'={}, j'={}, k'-1={}",
                     i+1, f.board, f.col, f.kp, f.jp, f.i_val.unwrap());
            println!("    A_{{k'-1,col}} = {}", format_poly(&f.f_poly));
            println!("    U_{{j'}}^+ = {}", format_poly(&f.g_poly));
        }
    }

    println!("\n================================================================");
    println!("  STRATEGY B: Induction on j' decreasing (right-cone on suffix sums)");
    println!("  Need: S_{{k'}} ≪ W_{{j',col}} for k' > j'");
    println!("================================================================");
    show("Key: D_{k'}^+ ≪ W_{j',col}", sb_dk_interl_wjp);
    if !sb_dk_interl_wjp_fails.is_empty() {
        println!("\n  --- First {} Strategy B failures ---", sb_dk_interl_wjp_fails.len());
        for (i, f) in sb_dk_interl_wjp_fails.iter().enumerate() {
            println!("  #{}: board={:?}, col={}, k'={}, j'={}",
                     i+1, f.board, f.col, f.kp, f.jp);
            println!("    D_{{k'}}^+ = S_{{k'}} = {}", format_poly(&f.f_poly));
            println!("    W_{{j',col}} = {}", format_poly(&f.g_poly));
        }
    }

    println!("\n================================================================");
    println!("  STRATEGY C: Mixed decomposition (right-cone via individual W_i)");
    println!("  Need: S_{{k'}} ≪ W_{{i,col}} for j' <= i < k'");
    println!("================================================================");
    show("Same-index: D_{k'}^+ ≪ U_{k'}^+ (known)", sc_same_du);
    show("Key: S_{k'} ≪ W_{i,col} for i < k'", sc_dk_interl_wi);
    if !sc_dk_interl_wi_fails.is_empty() {
        println!("\n  --- First {} Strategy C failures ---", sc_dk_interl_wi_fails.len());
        for (i, f) in sc_dk_interl_wi_fails.iter().enumerate() {
            println!("  #{}: board={:?}, col={}, k'={}, j'={}, i={}",
                     i+1, f.board, f.col, f.kp, f.jp, f.i_val.unwrap());
            println!("    S_{{k'}} = {}", format_poly(&f.f_poly));
            println!("    W_{{i,col}} = {}", format_poly(&f.g_poly));
        }
    }

    println!("\n================================================================");
    println!("  STRATEGY D: Direct decomposition via P_{{k'}}^+");
    println!("  D+U = P_{{k'}}^+ + sum W. Test sum W ≪ P and P ≪ sum W.");
    println!("================================================================");
    show("P_{k'}^+ real-rooted", sd_pk_rr);
    show("sum W_{i} ≪ P_{k'}^+", sd_sumw_interl_pk);
    show("P_{k'}^+ ≪ sum W_{i}", sd_pk_interl_sumw);
    if !sd_sumw_interl_pk_fails.is_empty() {
        println!("\n  --- First {} Strategy D (sum W ≪ P) failures ---", sd_sumw_interl_pk_fails.len());
        for (i, f) in sd_sumw_interl_pk_fails.iter().enumerate() {
            println!("  #{}: board={:?}, col={}, k'={}, j'={}",
                     i+1, f.board, f.col, f.kp, f.jp);
            println!("    sum W = {}", format_poly(&f.f_poly));
            println!("    P_{{k'}}^+ = {}", format_poly(&f.g_poly));
        }
    }

    println!("\n================================================================");
    println!("  SUMMARY");
    println!("================================================================");
    let strategies = [
        ("Ground truth (DU overlap)", du_overlap),
        ("Strategy A (A_{k'-1} ≪ U_{j'}^+)", sa_ak_interl_ujp),
        ("Strategy B (D_{k'}^+ ≪ W_{j',col})", sb_dk_interl_wjp),
        ("Strategy C (S_{k'} ≪ W_{i,col}, i<k')", sc_dk_interl_wi),
        ("Strategy D (sum W ≪ P_{k'}^+)", sd_sumw_interl_pk),
    ];
    for (name, c) in &strategies {
        let status = if c[0] == 0 { "NO DATA" }
                     else if c[1] == 0 { "ALL PASS  <-- potential proof path!" }
                     else { "HAS FAILURES" };
        println!("  {}: {}", name, status);
    }

    println!("\n  Total boards tested: {}", valid_boards);
    println!("================================================================");
}
