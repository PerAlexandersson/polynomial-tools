//! Verification of Conjecture 13 proof via the full inductive hypothesis.
//!
//! For each 312-avoiding Ferrers board λ with n <= max_n, verifies the
//! 6 within-column conditions (a')-(f') for each column l:
//!
//! (a') A_{j,l} << W_{j',l} for all j, j'
//! (b') W_{j,l} << W_{j',l} for j <= j'
//! (c') D_{j,l} << U_{j',l} for all j, j'   [Conjecture 13 (ii)]
//! (d') U_{j,l} << U_{j',l} for j <= j'
//! (e') D_{k,l} << D_{k',l} for k > k'       [Conjecture 13 (iii)]
//! (f') A_{k,l} << A_{k',l} for k > k'       [reversed AA]
//!
//! Also verifies the INDUCTIVE STEP: if conditions hold at level λ,
//! they hold at level λ^+ (within each group, i.e., for a fixed input column).
//!
//! Additionally checks cross-column conditions needed for cross-group:
//! (g) D_{j,l} << U_{j',l-1} for all j, j'    [cross-column DU]
//! (h) D_{j,l-1} << D_{j',l} for all j, j'    [cross-column DD - to check direction]
//!
//! Usage: cargo run --release --bin peak_conj13_proof -- 7

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

fn main() {
    let max_n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(7);
    println!("================================================================");
    println!("  Conjecture 13 proof verification");
    println!("  312-avoiding Ferrers boards, n <= {}", max_n);
    println!("================================================================\n");

    // Within-column counters
    let mut aw = [0u64; 2];    // (a') A_{j,l} << W_{j',l}
    let mut ww = [0u64; 2];    // (b') W_{j,l} << W_{j',l} for j <= j'
    let mut du = [0u64; 2];    // (c') D_{j,l} << U_{j',l}
    let mut uu = [0u64; 2];    // (d') U_{j,l} << U_{j',l} for j <= j'
    let mut dd_rev = [0u64; 2]; // (e') D_{k,l} << D_{k',l} for k > k'
    let mut aa_rev = [0u64; 2]; // (f') A_{k,l} << A_{k',l} for k > k'

    // New conditions
    let mut aw_rev = [0u64; 2];   // weaker AW: A_{j,l} << W_{j',l} for j >= j' only
    let mut dw = [0u64; 2];       // DW: D_{j,l} << W_{j',l} for all j, j'

    // Failure detail collectors (first 5 each)
    struct FailDetail {
        board: Vec<u8>,
        l: usize,
        j: usize,
        jp: usize,
        f_poly: Vec<i64>,
        g_poly: Vec<i64>,
        f_deg: Option<usize>,
        g_deg: Option<usize>,
    }
    let mut aw_fails: Vec<FailDetail> = Vec::new();
    let mut ww_fails: Vec<FailDetail> = Vec::new();
    let mut uu_fails: Vec<FailDetail> = Vec::new();

    // Cross-column counters
    let mut x_du = [0u64; 2];    // D_{j,l} << U_{j',l-1}
    let mut x_du_rev = [0u64; 2]; // D_{j,l-1} << U_{j',l}
    let mut x_aw = [0u64; 2];    // A_{j,l} << W_{j',l-1}
    let mut x_aw_rev = [0u64; 2]; // A_{j,l-1} << W_{j',l}
    let mut x_dd_fwd = [0u64; 2]; // D_{j,l-1} << D_{j',l} (forward across columns)
    let mut x_dd_rev2 = [0u64; 2]; // D_{j,l} << D_{j',l-1}

    // Inductive step counters (within-group at lambda+ level)
    let mut ind_du = [0u64; 2];    // (c') at lambda+
    let mut ind_dd = [0u64; 2];    // (e') at lambda+
    let mut ind_aw = [0u64; 2];    // (a') at lambda+
    let mut ind_ww = [0u64; 2];    // (b') at lambda+
    let mut ind_uu = [0u64; 2];    // (d') at lambda+
    let mut ind_aa = [0u64; 2];    // (f') at lambda+

    // Cross-group at lambda+ level
    let mut xg_du = [0u64; 2];    // cross-group DU at lambda+
    let mut xg_dd = [0u64; 2];    // cross-group DD at lambda+

    let mut valid_boards = 0u64;

    for n in 1..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) { continue; }
            valid_boards += 1;

            let m = board[0] as usize;

            // Skip trivial cases
            if n <= 1 { continue; }

            let ideal = bruhat_lower_ideal(&perm);

            // ── Compute D_{j,l} and U_{j,l} ──
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

            // Compute A_{j,l} = D + U and W_{j,l} = tD + U
            let mut a_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            let mut w_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            for j in 1..=m {
                for l in 1..=m {
                    a_jl[j][l] = pa(&d_jl[j][l], &u_jl[j][l]);
                    w_jl[j][l] = pa(&pmt(&d_jl[j][l]), &u_jl[j][l]);
                }
            }

            // ── Check within-column conditions for each l ──
            for l in 1..=m {
                for j in 1..=m {
                    for jp in 1..=m {
                        // (a') A_{j,l} << W_{j',l}
                        if !pz(&a_jl[j][l]) && !pz(&w_jl[jp][l]) {
                            aw[0] += 1;
                            if !interlaces(&a_jl[j][l], &w_jl[jp][l]) {
                                aw[1] += 1;
                                if aw_fails.len() < 5 {
                                    aw_fails.push(FailDetail {
                                        board: board.clone(), l, j, jp,
                                        f_poly: pt(&a_jl[j][l]),
                                        g_poly: pt(&w_jl[jp][l]),
                                        f_deg: pdeg(&a_jl[j][l]),
                                        g_deg: pdeg(&w_jl[jp][l]),
                                    });
                                }
                            }
                        }
                        // Weaker AW: only j >= j'
                        if j >= jp && !pz(&a_jl[j][l]) && !pz(&w_jl[jp][l]) {
                            aw_rev[0] += 1;
                            if !interlaces(&a_jl[j][l], &w_jl[jp][l]) {
                                aw_rev[1] += 1;
                            }
                        }
                        // DW: D_{j,l} << W_{j',l}
                        if !pz(&d_jl[j][l]) && !pz(&w_jl[jp][l]) {
                            dw[0] += 1;
                            if !interlaces(&d_jl[j][l], &w_jl[jp][l]) {
                                dw[1] += 1;
                            }
                        }
                        // (c') D_{j,l} << U_{j',l}
                        if !pz(&d_jl[j][l]) && !pz(&u_jl[jp][l]) {
                            du[0] += 1;
                            if !interlaces(&d_jl[j][l], &u_jl[jp][l]) {
                                du[1] += 1;
                            }
                        }
                    }
                    for jp in j..=m {
                        // (b') W_{j,l} << W_{j',l} for j <= j'
                        if !pz(&w_jl[j][l]) && !pz(&w_jl[jp][l]) && j < jp {
                            ww[0] += 1;
                            if !interlaces(&w_jl[j][l], &w_jl[jp][l]) {
                                ww[1] += 1;
                                if ww_fails.len() < 5 {
                                    ww_fails.push(FailDetail {
                                        board: board.clone(), l, j, jp,
                                        f_poly: pt(&w_jl[j][l]),
                                        g_poly: pt(&w_jl[jp][l]),
                                        f_deg: pdeg(&w_jl[j][l]),
                                        g_deg: pdeg(&w_jl[jp][l]),
                                    });
                                }
                            }
                        }
                        // (d') U_{j,l} << U_{j',l} for j <= j'
                        if !pz(&u_jl[j][l]) && !pz(&u_jl[jp][l]) && j < jp {
                            uu[0] += 1;
                            if !interlaces(&u_jl[j][l], &u_jl[jp][l]) {
                                uu[1] += 1;
                                if uu_fails.len() < 5 {
                                    uu_fails.push(FailDetail {
                                        board: board.clone(), l, j, jp,
                                        f_poly: pt(&u_jl[j][l]),
                                        g_poly: pt(&u_jl[jp][l]),
                                        f_deg: pdeg(&u_jl[j][l]),
                                        g_deg: pdeg(&u_jl[jp][l]),
                                    });
                                }
                            }
                        }
                    }
                    // (e') D_{k,l} << D_{k',l} for k > k' (reversed DD)
                    for kp in 1..j {
                        if !pz(&d_jl[j][l]) && !pz(&d_jl[kp][l]) {
                            dd_rev[0] += 1;
                            if !interlaces(&d_jl[j][l], &d_jl[kp][l]) {
                                dd_rev[1] += 1;
                            }
                        }
                    }
                    // (f') A_{k,l} << A_{k',l} for k > k' (reversed AA)
                    for kp in 1..j {
                        if !pz(&a_jl[j][l]) && !pz(&a_jl[kp][l]) {
                            aa_rev[0] += 1;
                            if !interlaces(&a_jl[j][l], &a_jl[kp][l]) {
                                aa_rev[1] += 1;
                            }
                        }
                    }
                }
            }

            // ── Check cross-column conditions ──
            for l in 2..=m {
                for j in 1..=m {
                    for jp in 1..=m {
                        // D_{j,l} << U_{j',l-1}
                        if !pz(&d_jl[j][l]) && !pz(&u_jl[jp][l-1]) {
                            x_du[0] += 1;
                            if !interlaces(&d_jl[j][l], &u_jl[jp][l-1]) {
                                x_du[1] += 1;
                            }
                        }
                        // D_{j,l-1} << U_{j',l}
                        if !pz(&d_jl[j][l-1]) && !pz(&u_jl[jp][l]) {
                            x_du_rev[0] += 1;
                            if !interlaces(&d_jl[j][l-1], &u_jl[jp][l]) {
                                x_du_rev[1] += 1;
                            }
                        }
                        // A_{j,l} << W_{j',l-1}
                        if !pz(&a_jl[j][l]) && !pz(&w_jl[jp][l-1]) {
                            x_aw[0] += 1;
                            if !interlaces(&a_jl[j][l], &w_jl[jp][l-1]) {
                                x_aw[1] += 1;
                            }
                        }
                        // A_{j,l-1} << W_{j',l}
                        if !pz(&a_jl[j][l-1]) && !pz(&w_jl[jp][l]) {
                            x_aw_rev[0] += 1;
                            if !interlaces(&a_jl[j][l-1], &w_jl[jp][l]) {
                                x_aw_rev[1] += 1;
                            }
                        }
                        // D_{j,l-1} << D_{j',l}
                        if !pz(&d_jl[j][l-1]) && !pz(&d_jl[jp][l]) {
                            x_dd_fwd[0] += 1;
                            if !interlaces(&d_jl[j][l-1], &d_jl[jp][l]) {
                                x_dd_fwd[1] += 1;
                            }
                        }
                        // D_{j,l} << D_{j',l-1}
                        if !pz(&d_jl[j][l]) && !pz(&d_jl[jp][l-1]) {
                            x_dd_rev2[0] += 1;
                            if !interlaces(&d_jl[j][l], &d_jl[jp][l-1]) {
                                x_dd_rev2[1] += 1;
                            }
                        }
                    }
                }
            }

            // ── Inductive step: check conditions at lambda+ level ──
            let mp = m + 1; // m' for lambda+

            // Precompute prefix sums (S) and suffix sums (T) for each column
            // S_{k}^{(col)} = sum_{j=1}^{k-1} A_{j,col}
            // T_{k}^{(col)} = sum_{j=k}^{m} W_{j,col}
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

            // D_{k',l'}^+ and U_{k',l'}^+ for each (k', l') with k' != l'
            // D = S_k'^{(col)}, U = T_k'^{(col)}
            // For k' > l': col = l'
            // For k' < l': col = l'-1

            for lp in 1..=mp {
                // Collect all (k', D^+, U^+, A^+, W^+, col) for this l'
                let mut entries: Vec<(usize, Vec<i64>, Vec<i64>, Vec<i64>, Vec<i64>, usize)> = Vec::new();

                for kp in 1..=mp {
                    if kp == lp { continue; }
                    let col = if kp > lp { lp } else { lp - 1 };
                    if col < 1 || col > m { continue; }
                    let sk = if kp < s_col[col].len() { s_col[col][kp].clone() } else { s_col[col][m + 1].clone() };
                    let tk = if kp < t_col[col].len() { t_col[col][kp].clone() } else { vec![0i64] };
                    let ak = pa(&sk, &tk);        // A^+ = S + T
                    let wk = pa(&pmt(&sk), &tk);  // W^+ = tS + T
                    entries.push((kp, sk, tk, ak, wk, col));
                }

                // Check conditions within same group
                for i in 0..entries.len() {
                    let (kp_i, ref d_i, ref u_i, ref a_i, ref w_i, col_i) = entries[i];
                    for j in 0..entries.len() {
                        let (kp_j, ref d_j, ref u_j, ref a_j, ref w_j, col_j) = entries[j];

                        let same_group = col_i == col_j;

                        // (c') D_{k'} << U_{j'} - both within-group and cross-group
                        if !pz(d_i) && !pz(u_j) {
                            if same_group {
                                ind_du[0] += 1;
                                if !interlaces(d_i, u_j) { ind_du[1] += 1; }
                            } else {
                                xg_du[0] += 1;
                                if !interlaces(d_i, u_j) { xg_du[1] += 1; }
                            }
                        }

                        // (e') D_{k'} << D_{k''} for k' > k''
                        if kp_i > kp_j && !pz(d_i) && !pz(d_j) {
                            if same_group {
                                ind_dd[0] += 1;
                                if !interlaces(d_i, d_j) { ind_dd[1] += 1; }
                            } else {
                                xg_dd[0] += 1;
                                if !interlaces(d_i, d_j) { xg_dd[1] += 1; }
                            }
                        }

                        if same_group {
                            // (a') A << W
                            if !pz(a_i) && !pz(w_j) {
                                ind_aw[0] += 1;
                                if !interlaces(a_i, w_j) { ind_aw[1] += 1; }
                            }
                            // (b') W_i << W_j for i <= j
                            if kp_i < kp_j && !pz(w_i) && !pz(w_j) {
                                ind_ww[0] += 1;
                                if !interlaces(w_i, w_j) { ind_ww[1] += 1; }
                            }
                            // (d') U_i << U_j for i <= j
                            if kp_i < kp_j && !pz(u_i) && !pz(u_j) {
                                ind_uu[0] += 1;
                                if !interlaces(u_i, u_j) { ind_uu[1] += 1; }
                            }
                            // (f') A_k << A_k' for k > k'
                            if kp_i > kp_j && !pz(a_i) && !pz(a_j) {
                                ind_aa[0] += 1;
                                if !interlaces(a_i, a_j) { ind_aa[1] += 1; }
                            }
                        }
                    }
                }
            }
        }
        println!("n={}: {} boards (cumulative: {})", n, boards.iter().filter(|b| is_312_avoiding(&board_to_perm(b))).count(), valid_boards);
    }

    println!("\n================================================================");
    println!("  WITHIN-COLUMN CONDITIONS at lambda level");
    println!("================================================================");
    show("(a') A_{j,l} << W_{j',l}  [AW all j,j']", aw);
    show("(a'r) A_{j,l} << W_{j',l} [AW j>=j' only]", aw_rev);
    show("(b') W_{j,l} << W_{j',l}  [WW j<=j']", ww);
    show("(c') D_{j,l} << U_{j',l}  [DU = Conj13(ii)]", du);
    show("(d') U_{j,l} << U_{j',l}  [UU j<=j']", uu);
    show("(e') D_{k,l} << D_{k',l}  [rev DD = Conj13(iii)]", dd_rev);
    show("(f') A_{k,l} << A_{k',l}  [rev AA]", aa_rev);
    show("     D_{j,l} << W_{j',l}  [DW all j,j']", dw);

    // Print failure details
    if !aw_fails.is_empty() {
        println!("\n  --- First {} AW failures ---", aw_fails.len());
        for (i, f) in aw_fails.iter().enumerate() {
            println!("  AW fail #{}: board={:?}, l={}, j={}, j'={}",
                     i+1, f.board, f.l, f.j, f.jp);
            println!("    A_{{j,l}} = {} (deg {:?})", format_poly(&f.f_poly), f.f_deg);
            println!("    W_{{j',l}} = {} (deg {:?})", format_poly(&f.g_poly), f.g_deg);
        }
    }
    if !ww_fails.is_empty() {
        println!("\n  --- First {} WW failures ---", ww_fails.len());
        for (i, f) in ww_fails.iter().enumerate() {
            println!("  WW fail #{}: board={:?}, l={}, j={}, j'={}",
                     i+1, f.board, f.l, f.j, f.jp);
            println!("    W_{{j,l}} = {} (deg {:?})", format_poly(&f.f_poly), f.f_deg);
            println!("    W_{{j',l}} = {} (deg {:?})", format_poly(&f.g_poly), f.g_deg);
        }
    }
    if !uu_fails.is_empty() {
        println!("\n  --- First {} UU failures ---", uu_fails.len());
        for (i, f) in uu_fails.iter().enumerate() {
            println!("  UU fail #{}: board={:?}, l={}, j={}, j'={}",
                     i+1, f.board, f.l, f.j, f.jp);
            println!("    U_{{j,l}} = {} (deg {:?})", format_poly(&f.f_poly), f.f_deg);
            println!("    U_{{j',l}} = {} (deg {:?})", format_poly(&f.g_poly), f.g_deg);
        }
    }

    println!("\n================================================================");
    println!("  CROSS-COLUMN CONDITIONS at lambda level");
    println!("================================================================");
    show("D_{j,l} << U_{j',l-1}    [cross DU →]", x_du);
    show("D_{j,l-1} << U_{j',l}    [cross DU ←]", x_du_rev);
    show("A_{j,l} << W_{j',l-1}    [cross AW →]", x_aw);
    show("A_{j,l-1} << W_{j',l}    [cross AW ←]", x_aw_rev);
    show("D_{j,l-1} << D_{j',l}    [cross DD fwd]", x_dd_fwd);
    show("D_{j,l} << D_{j',l-1}    [cross DD rev]", x_dd_rev2);

    println!("\n================================================================");
    println!("  INDUCTIVE STEP: within-group at lambda+ level");
    println!("================================================================");
    show("(c') D << U [DU within-group]", ind_du);
    show("(e') D << D [rev DD within-group]", ind_dd);
    show("(a') A << W [AW within-group]", ind_aw);
    show("(b') W << W [WW within-group]", ind_ww);
    show("(d') U << U [UU within-group]", ind_uu);
    show("(f') A << A [rev AA within-group]", ind_aa);

    println!("\n================================================================");
    println!("  CROSS-GROUP at lambda+ level");
    println!("================================================================");
    show("D << U [DU cross-group]", xg_du);
    show("D << D [rev DD cross-group]", xg_dd);

    println!("\n================================================================");
    println!("  Total boards tested: {}", valid_boards);
    println!("================================================================");
}
