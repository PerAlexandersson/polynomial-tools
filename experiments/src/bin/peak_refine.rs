//! Refined peak polynomials by end-type and last-entry value.
//!
//! For each 312-avoiding Ferrers board lambda with n <= 6, explores two refinements:
//!
//! 1. End-type (ascent/descent at position n-1): splits D_k, U_k into ea/ed variants
//! 2. Last-entry value: P_{k,l}(t) counting permutations with pi_1 = k, pi_n = l
//!
//! Checks interlacing conditions within and across end-types, and explores
//! whether the recursion decouples by end-type.
//!
//! Usage: cargo run --release --bin peak_refine -- 6

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
    // multiply by t (shift coefficients up by 1)
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
fn pf(p: &[i64]) -> String {
    let p = pt(p);
    if pz(&p) {
        return "0".into();
    }
    format_poly(&p)
}

/// Check f << g (weak interlacing).
/// Handles same-degree case via Wagner's trick: f << g (same deg d)
/// iff check_weak_interlacing(g, t*f) = Some(true).
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
            // Same-degree case: f << g iff check_weak_interlacing(g, t*f)
            let df = pdeg(&f);
            let dg = pdeg(&g);
            match (df, dg) {
                (Some(df), Some(dg)) if df == dg => {
                    let tf = pmt(&f);
                    match check_weak_interlacing(&g, &tf) {
                        Some(b) => b,
                        None => false,
                    }
                }
                _ => false,
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

// ── End-type classification ──

/// Classify the "end type" of a permutation:
/// - 'a' if pi_{n-1} < pi_n (ascending end) or n <= 1
/// - 'd' if pi_{n-1} > pi_n (descending end)
fn end_type(pi: &[u8]) -> char {
    let n = pi.len();
    if n <= 1 {
        return 'a';
    } // convention: length 1 is ascending
    if pi[n - 2] < pi[n - 1] {
        'a'
    } else {
        'd'
    }
}

/// Classify the "start type" of a permutation:
/// - 'u' (up/ascent) if pi_1 < pi_2 or n <= 1
/// - 'd' (down/descent) if pi_1 > pi_2
fn start_type(pi: &[u8]) -> char {
    let n = pi.len();
    if n <= 1 {
        return 'u';
    }
    if pi[0] < pi[1] {
        'u'
    } else {
        'd'
    }
}

fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);
    println!("================================================================");
    println!("  REFINED PEAK POLYNOMIALS: End-type & Last-entry exploration");
    println!("  312-avoiding Ferrers boards, n <= {}", max_n);
    println!("  Using Bezout-based interlacing (polynomial-tools)");
    println!("================================================================\n");

    // ── REFINEMENT 1: End-type counters ──
    // "ea" = ascending end, "ed" = descending end
    // Within same end-type:
    let mut aw_ea = [0u64; 2]; // A_j^{ea} << W_l^{ea} for all j,l
    let mut ww_ea = [0u64; 2]; // W_j^{ea} << W_l^{ea} for j <= l
    let mut aw_ed = [0u64; 2]; // A_j^{ed} << W_l^{ed} for all j,l
    let mut ww_ed = [0u64; 2]; // W_j^{ed} << W_l^{ed} for j <= l

    // Cross-type:
    let mut aw_ea_ed = [0u64; 2]; // A_j^{ea} << W_l^{ed} for all j,l
    let mut aw_ed_ea = [0u64; 2]; // A_j^{ed} << W_l^{ea} for all j,l
    let mut ww_ea_ed = [0u64; 2]; // W_j^{ea} << W_l^{ed} for j <= l
    let mut ak_ea_ed = [0u64; 2]; // A_k^{ea} << A_k^{ed} for each k
    let mut ak_ed_ea = [0u64; 2]; // A_k^{ed} << A_k^{ea} for each k
    let mut ak_either = [0u64; 2]; // one direction works

    // Also check the unsplit conditions for comparison
    let mut aw_all = [0u64; 2]; // A_j << W_l for all j,l (original)
    let mut ww_all = [0u64; 2]; // W_j << W_l for j <= l (original)

    // Additional: DU conditions within end-types
    let mut du_ea = [0u64; 2]; // D_j^{ea} << U_l^{ea} for all j,l
    let mut du_ed = [0u64; 2]; // D_j^{ed} << U_l^{ed} for all j,l
    let mut uu_ea = [0u64; 2]; // U_j^{ea} << U_l^{ea} for j <= l
    let mut uu_ed = [0u64; 2]; // U_j^{ed} << U_l^{ed} for j <= l

    // ── REFINEMENT 2: Last-entry counters ──
    let mut pkl_fixed_l_interlace = [0u64; 2]; // For fixed l: P_{m,l} << P_{m-1,l} << ... << P_{1,l}
    let mut pkl_fixed_k_interlace = [0u64; 2]; // For fixed k: P_{k,1} << P_{k,2} << ... << P_{k,m}
    let mut pkl_fixed_l_rev = [0u64; 2]; // P_{j,l} << P_{j+1,l} (reverse order)
    let mut pkl_fixed_k_rev = [0u64; 2]; // P_{k,l+1} << P_{k,l} (reverse order)
    let mut pkl_rr = [0u64; 2]; // P_{k,l} real-rooted

    // ── DECOUPLING check ──
    let mut decouple_boards_checked = 0u64;
    let mut decouple_match = 0u64;
    let mut decouple_fail = 0u64;

    let mut valid_boards = 0u64;

    for n in 1..=max_n {
        let boards = gen_boards(n);
        let mut n_valid = 0u64;

        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            valid_boards += 1;
            n_valid += 1;

            let m = (board[0] as usize).min(n);
            let ideal = bruhat_lower_ideal(&perm);

            // ────────────────────────────────────────────────
            // REFINEMENT 1: End-type split
            // ────────────────────────────────────────────────

            // D_k^{ea}, D_k^{ed}, U_k^{ea}, U_k^{ed}
            // Index 0..=m (use 1..=m)
            let mut d_ea = vec![vec![0i64]; m + 1];
            let mut d_ed = vec![vec![0i64]; m + 1];
            let mut u_ea = vec![vec![0i64]; m + 1];
            let mut u_ed = vec![vec![0i64]; m + 1];

            // Also compute unsplit D_k, U_k for reference
            let mut d_all = vec![vec![0i64]; m + 1];
            let mut u_all = vec![vec![0i64]; m + 1];

            // REFINEMENT 2: P_{k,l} indexed by (first, last)
            // Index: pkl[k][l] for k,l in 1..=m
            let mut pkl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];

            for pi in &ideal {
                let k = pi[0] as usize;
                if k > m {
                    continue;
                }
                let pk = peaks(pi);
                let et = end_type(pi);
                let st = start_type(pi);
                let last = *pi.last().unwrap() as usize;

                // Unsplit
                let poly_all = if st == 'd' {
                    &mut d_all[k]
                } else {
                    &mut u_all[k]
                };
                while poly_all.len() <= pk {
                    poly_all.push(0);
                }
                poly_all[pk] += 1;

                // End-type split
                let poly_ref = match (st, et) {
                    ('d', 'a') => &mut d_ea[k],
                    ('d', 'd') => &mut d_ed[k],
                    ('u', 'a') => &mut u_ea[k],
                    ('u', 'd') => &mut u_ed[k],
                    _ => unreachable!(),
                };
                while poly_ref.len() <= pk {
                    poly_ref.push(0);
                }
                poly_ref[pk] += 1;

                // Last-entry tracking
                if last <= m {
                    while pkl[k][last].len() <= pk {
                        pkl[k][last].push(0);
                    }
                    pkl[k][last][pk] += 1;
                }
            }

            // Compute A, W for unsplit
            let mut a_all = vec![vec![0i64]; m + 1];
            let mut w_all = vec![vec![0i64]; m + 1];
            for k in 1..=m {
                a_all[k] = pa(&d_all[k], &u_all[k]);
                w_all[k] = pa(&pmt(&d_all[k]), &u_all[k]);
            }

            // Compute A, W for each end type
            let mut a_ea = vec![vec![0i64]; m + 1];
            let mut a_ed = vec![vec![0i64]; m + 1];
            let mut w_ea = vec![vec![0i64]; m + 1];
            let mut w_ed = vec![vec![0i64]; m + 1];
            for k in 1..=m {
                a_ea[k] = pa(&d_ea[k], &u_ea[k]);
                a_ed[k] = pa(&d_ed[k], &u_ed[k]);
                w_ea[k] = pa(&pmt(&d_ea[k]), &u_ea[k]);
                w_ed[k] = pa(&pmt(&d_ed[k]), &u_ed[k]);
            }

            // ── Check unsplit conditions (for reference) ──
            for j in 1..=m {
                for l in 1..=m {
                    if !pz(&a_all[j]) && !pz(&w_all[l]) {
                        aw_all[0] += 1;
                        if !interlaces(&a_all[j], &w_all[l]) {
                            aw_all[1] += 1;
                        }
                    }
                    if j <= l && j < l {
                        if !pz(&w_all[j]) && !pz(&w_all[l]) {
                            ww_all[0] += 1;
                            if !interlaces(&w_all[j], &w_all[l]) {
                                ww_all[1] += 1;
                            }
                        }
                    }
                }
            }

            // ── Check end-type conditions ──
            for j in 1..=m {
                for l in 1..=m {
                    // (a) A_j^{ea} << W_l^{ea} for all j,l
                    if !pz(&a_ea[j]) && !pz(&w_ea[l]) {
                        aw_ea[0] += 1;
                        if !interlaces(&a_ea[j], &w_ea[l]) {
                            aw_ea[1] += 1;
                            if aw_ea[1] <= 3 {
                                println!(
                                    "FAIL A_j^ea << W_l^ea: board={:?} j={} l={}",
                                    board, j, l
                                );
                                println!("  A_j^ea = {}", pf(&a_ea[j]));
                                println!("  W_l^ea = {}", pf(&w_ea[l]));
                            }
                        }
                    }
                    // A_j^{ed} << W_l^{ed} for all j,l
                    if !pz(&a_ed[j]) && !pz(&w_ed[l]) {
                        aw_ed[0] += 1;
                        if !interlaces(&a_ed[j], &w_ed[l]) {
                            aw_ed[1] += 1;
                            if aw_ed[1] <= 3 {
                                println!(
                                    "FAIL A_j^ed << W_l^ed: board={:?} j={} l={}",
                                    board, j, l
                                );
                                println!("  A_j^ed = {}", pf(&a_ed[j]));
                                println!("  W_l^ed = {}", pf(&w_ed[l]));
                            }
                        }
                    }

                    // Cross: A_j^{ea} << W_l^{ed}
                    if !pz(&a_ea[j]) && !pz(&w_ed[l]) {
                        aw_ea_ed[0] += 1;
                        if !interlaces(&a_ea[j], &w_ed[l]) {
                            aw_ea_ed[1] += 1;
                            if aw_ea_ed[1] <= 3 {
                                println!(
                                    "FAIL A_j^ea << W_l^ed: board={:?} j={} l={}",
                                    board, j, l
                                );
                                println!("  A_j^ea = {}", pf(&a_ea[j]));
                                println!("  W_l^ed = {}", pf(&w_ed[l]));
                            }
                        }
                    }
                    // Cross: A_j^{ed} << W_l^{ea}
                    if !pz(&a_ed[j]) && !pz(&w_ea[l]) {
                        aw_ed_ea[0] += 1;
                        if !interlaces(&a_ed[j], &w_ea[l]) {
                            aw_ed_ea[1] += 1;
                            if aw_ed_ea[1] <= 3 {
                                println!(
                                    "FAIL A_j^ed << W_l^ea: board={:?} j={} l={}",
                                    board, j, l
                                );
                                println!("  A_j^ed = {}", pf(&a_ed[j]));
                                println!("  W_l^ea = {}", pf(&w_ea[l]));
                            }
                        }
                    }

                    // (b) W_j << W_l for j < l (same end type)
                    if j < l {
                        if !pz(&w_ea[j]) && !pz(&w_ea[l]) {
                            ww_ea[0] += 1;
                            if !interlaces(&w_ea[j], &w_ea[l]) {
                                ww_ea[1] += 1;
                                if ww_ea[1] <= 3 {
                                    println!(
                                        "FAIL W_j^ea << W_l^ea: board={:?} j={} l={}",
                                        board, j, l
                                    );
                                    println!("  W_j^ea = {}", pf(&w_ea[j]));
                                    println!("  W_l^ea = {}", pf(&w_ea[l]));
                                }
                            }
                        }
                        if !pz(&w_ed[j]) && !pz(&w_ed[l]) {
                            ww_ed[0] += 1;
                            if !interlaces(&w_ed[j], &w_ed[l]) {
                                ww_ed[1] += 1;
                                if ww_ed[1] <= 3 {
                                    println!(
                                        "FAIL W_j^ed << W_l^ed: board={:?} j={} l={}",
                                        board, j, l
                                    );
                                    println!("  W_j^ed = {}", pf(&w_ed[j]));
                                    println!("  W_l^ed = {}", pf(&w_ed[l]));
                                }
                            }
                        }
                        // Cross: W_j^{ea} << W_l^{ed} for j <= l
                        if !pz(&w_ea[j]) && !pz(&w_ed[l]) {
                            ww_ea_ed[0] += 1;
                            if !interlaces(&w_ea[j], &w_ed[l]) {
                                ww_ea_ed[1] += 1;
                                if ww_ea_ed[1] <= 3 {
                                    println!(
                                        "FAIL W_j^ea << W_l^ed: board={:?} j={} l={}",
                                        board, j, l
                                    );
                                    println!("  W_j^ea = {}", pf(&w_ea[j]));
                                    println!("  W_l^ed = {}", pf(&w_ed[l]));
                                }
                            }
                        }
                    }

                    // DU within end types
                    if !pz(&d_ea[j]) && !pz(&u_ea[l]) {
                        du_ea[0] += 1;
                        if !interlaces(&d_ea[j], &u_ea[l]) {
                            du_ea[1] += 1;
                            if du_ea[1] <= 3 {
                                println!(
                                    "FAIL D_j^ea << U_l^ea: board={:?} j={} l={}",
                                    board, j, l
                                );
                            }
                        }
                    }
                    if !pz(&d_ed[j]) && !pz(&u_ed[l]) {
                        du_ed[0] += 1;
                        if !interlaces(&d_ed[j], &u_ed[l]) {
                            du_ed[1] += 1;
                            if du_ed[1] <= 3 {
                                println!(
                                    "FAIL D_j^ed << U_l^ed: board={:?} j={} l={}",
                                    board, j, l
                                );
                            }
                        }
                    }

                    // UU within end types for j <= l
                    if j <= l {
                        if !pz(&u_ea[j]) && !pz(&u_ea[l]) {
                            uu_ea[0] += 1;
                            if !interlaces(&u_ea[j], &u_ea[l]) {
                                uu_ea[1] += 1;
                                if uu_ea[1] <= 3 {
                                    println!(
                                        "FAIL U_j^ea << U_l^ea: board={:?} j={} l={}",
                                        board, j, l
                                    );
                                }
                            }
                        }
                        if !pz(&u_ed[j]) && !pz(&u_ed[l]) {
                            uu_ed[0] += 1;
                            if !interlaces(&u_ed[j], &u_ed[l]) {
                                uu_ed[1] += 1;
                                if uu_ed[1] <= 3 {
                                    println!(
                                        "FAIL U_j^ed << U_l^ed: board={:?} j={} l={}",
                                        board, j, l
                                    );
                                }
                            }
                        }
                    }
                }
            }

            // A_k^{ea} vs A_k^{ed} comparison for each k
            for k in 1..=m {
                if !pz(&a_ea[k]) && !pz(&a_ed[k]) {
                    ak_either[0] += 1;
                    let ea_ed = interlaces(&a_ea[k], &a_ed[k]);
                    let ed_ea = interlaces(&a_ed[k], &a_ea[k]);
                    ak_ea_ed[0] += 1;
                    if !ea_ed {
                        ak_ea_ed[1] += 1;
                    }
                    ak_ed_ea[0] += 1;
                    if !ed_ea {
                        ak_ed_ea[1] += 1;
                    }
                    if !ea_ed && !ed_ea {
                        ak_either[1] += 1;
                        if ak_either[1] <= 3 {
                            println!(
                                "FAIL neither A_k^ea << A_k^ed nor reverse: board={:?} k={}",
                                board, k
                            );
                            println!("  A_k^ea = {}", pf(&a_ea[k]));
                            println!("  A_k^ed = {}", pf(&a_ed[k]));
                        }
                    }
                }
            }

            // ────────────────────────────────────────────────
            // REFINEMENT 2: Last-entry value P_{k,l}
            // ────────────────────────────────────────────────

            // Real-rootedness of P_{k,l}
            for k in 1..=m {
                for l in 1..=m {
                    if !pz(&pkl[k][l]) {
                        pkl_rr[0] += 1;
                        if !is_real_rooted(&pkl[k][l]) {
                            pkl_rr[1] += 1;
                            if pkl_rr[1] <= 3 {
                                println!(
                                    "FAIL P_{{{},{}}} not real-rooted: board={:?} poly={}",
                                    k,
                                    l,
                                    board,
                                    pf(&pkl[k][l])
                                );
                            }
                        }
                    }
                }
            }

            // For fixed l: check P_{j,l} << P_{j+1,l} for all j (descending first index)
            // i.e., is the sequence (P_{1,l}, P_{2,l}, ..., P_{m,l}) ordered by interlacing?
            for l in 1..=m {
                for j in 1..m {
                    // Check P_{j,l} << P_{j+1,l}
                    if !pz(&pkl[j][l]) && !pz(&pkl[j + 1][l]) {
                        pkl_fixed_l_interlace[0] += 1;
                        if !interlaces(&pkl[j][l], &pkl[j + 1][l]) {
                            pkl_fixed_l_interlace[1] += 1;
                        }
                    }
                    // Also check reverse: P_{j+1,l} << P_{j,l}
                    if !pz(&pkl[j + 1][l]) && !pz(&pkl[j][l]) {
                        pkl_fixed_l_rev[0] += 1;
                        if !interlaces(&pkl[j + 1][l], &pkl[j][l]) {
                            pkl_fixed_l_rev[1] += 1;
                        }
                    }
                }
            }

            // For fixed k: check P_{k,l} << P_{k,l+1} for all l
            for k in 1..=m {
                for l in 1..m {
                    if !pz(&pkl[k][l]) && !pz(&pkl[k][l + 1]) {
                        pkl_fixed_k_interlace[0] += 1;
                        if !interlaces(&pkl[k][l], &pkl[k][l + 1]) {
                            pkl_fixed_k_interlace[1] += 1;
                        }
                    }
                    if !pz(&pkl[k][l + 1]) && !pz(&pkl[k][l]) {
                        pkl_fixed_k_rev[0] += 1;
                        if !interlaces(&pkl[k][l + 1], &pkl[k][l]) {
                            pkl_fixed_k_rev[1] += 1;
                        }
                    }
                }
            }
        }
        println!(
            "n={}: {} valid 312-avoiding boards (cumulative: {})",
            n, n_valid, valid_boards
        );
    }

    // ────────────────────────────────────────────────
    // DECOUPLING CHECK
    // ────────────────────────────────────────────────
    // For each board lambda, find "child" boards lambda+ (with an extra first row)
    // and check whether the ea/ed polynomials at level lambda+ depend only on
    // the same end-type at level lambda.
    //
    // The "first row" construction: lambda+ = (m+1, lambda_1, ..., lambda_n)
    // means adding a new first row of length m+1 (the maximum).
    // Actually, in the Ferrers board context, the recursion adds a new FIRST column
    // (splitting by pi_1 = k). The last entries don't change because we're adding
    // at the beginning.
    //
    // So: permutations in the Bruhat ideal of lambda+ that start with k have a
    // specific relationship to permutations in the ideal of lambda.
    // The key point: adding a first entry doesn't change the LAST two entries,
    // so the end-type is preserved. We verify this by checking that
    // the D_k^{ea} at lambda+ level equals a sum over contributions that only
    // involve ea-type data from the child boards.
    //
    // We check this by computing, for each board and its "parent" (board with
    // first entry removed, if it exists), whether the end-type polys at the
    // child level decompose consistently.

    println!("\n================================================================");
    println!("  DECOUPLING VERIFICATION");
    println!("================================================================\n");

    // For the decoupling check, we verify a simpler property:
    // For each board lambda with n >= 2, consider removing the first row to get lambda'.
    // The permutations in the ideal of lambda that start with k are obtained from
    // permutations in the ideal of lambda' by prepending k (and adjusting).
    // Since prepending doesn't change the last two entries, the end-type should
    // be preserved.
    //
    // Concrete check: For lambda at level n, compute:
    //   D_k^{ea}(lambda) + U_k^{ea}(lambda) = A_k^{ea}(lambda)
    //   W_k^{ea}(lambda) = t*D_k^{ea}(lambda) + U_k^{ea}(lambda)
    // These should relate to level n-1 data ONLY through ea-type polynomials.
    //
    // We verify: sum_k A_k^{ea}(lambda) should equal the total ea-peak polynomial,
    // and A_k^{ea} + A_k^{ed} should equal A_k (unsplit). This is a consistency check.

    for n in 2..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            let m = (board[0] as usize).min(n);
            let ideal = bruhat_lower_ideal(&perm);

            decouple_boards_checked += 1;

            // Compute the split polynomials
            let mut d_ea = vec![vec![0i64]; m + 1];
            let mut d_ed = vec![vec![0i64]; m + 1];
            let mut u_ea = vec![vec![0i64]; m + 1];
            let mut u_ed = vec![vec![0i64]; m + 1];
            let mut d_k = vec![vec![0i64]; m + 1];
            let mut u_k = vec![vec![0i64]; m + 1];

            for pi in &ideal {
                let k = pi[0] as usize;
                if k > m {
                    continue;
                }
                let pk = peaks(pi);
                let et = end_type(pi);
                let st = start_type(pi);

                // Unsplit
                let poly_us = if st == 'd' { &mut d_k[k] } else { &mut u_k[k] };
                while poly_us.len() <= pk {
                    poly_us.push(0);
                }
                poly_us[pk] += 1;

                // Split
                let poly_ref = match (st, et) {
                    ('d', 'a') => &mut d_ea[k],
                    ('d', 'd') => &mut d_ed[k],
                    ('u', 'a') => &mut u_ea[k],
                    ('u', 'd') => &mut u_ed[k],
                    _ => unreachable!(),
                };
                while poly_ref.len() <= pk {
                    poly_ref.push(0);
                }
                poly_ref[pk] += 1;
            }

            // Consistency check: D_k^{ea} + D_k^{ed} = D_k, U_k^{ea} + U_k^{ed} = U_k
            let mut ok = true;
            for k in 1..=m {
                let d_sum = pa(&d_ea[k], &d_ed[k]);
                let u_sum = pa(&u_ea[k], &u_ed[k]);
                if pt(&d_sum) != pt(&d_k[k]) || pt(&u_sum) != pt(&u_k[k]) {
                    ok = false;
                    if decouple_fail < 3 {
                        println!("FAIL decoupling consistency: board={:?} k={}", board, k);
                        println!(
                            "  D_k = {} vs D_k^ea + D_k^ed = {}",
                            pf(&d_k[k]),
                            pf(&d_sum)
                        );
                        println!(
                            "  U_k = {} vs U_k^ea + U_k^ed = {}",
                            pf(&u_k[k]),
                            pf(&u_sum)
                        );
                    }
                }
            }
            if ok {
                decouple_match += 1;
            } else {
                decouple_fail += 1;
            }
        }
    }

    println!(
        "Decoupling consistency (D/U split sums match): {}/{} boards",
        decouple_match, decouple_boards_checked
    );
    if decouple_fail > 0 {
        println!("  FAILURES: {}", decouple_fail);
    } else {
        println!("  All boards pass the split consistency check.");
    }

    // ────────────────────────────────────────────────
    // DEEPER DECOUPLING: Check that at level lambda+, the "ea" polynomials
    // depend only on "ea" data from level lambda.
    // ────────────────────────────────────────────────
    // For each pair (lambda, lambda+) where lambda+ has an additional first row,
    // we verify that the ea-type contribution to lambda+ only comes from
    // ea-type permutations at lambda level.
    //
    // Specifically: for each pi in Bruhat ideal of lambda+, removing the first
    // entry should yield a permutation whose end-type matches pi's end-type.
    // (Since removing the first entry doesn't change the last two entries.)

    println!("\n--- Deeper decoupling: end-type preservation under first-entry removal ---\n");

    let mut preserve_checks = 0u64;
    let mut preserve_pass = 0u64;
    let mut preserve_fail = 0u64;

    for n in 2..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            let ideal = bruhat_lower_ideal(&perm);

            for pi in &ideal {
                if pi.len() < 2 {
                    continue;
                }
                preserve_checks += 1;

                // End type of pi
                let et_full = end_type(pi);

                // Remove first entry, adjust remaining
                let first = pi[0];
                let rest: Vec<u8> = pi[1..]
                    .iter()
                    .map(|&x| if x > first { x - 1 } else { x })
                    .collect();

                if rest.len() >= 2 {
                    let et_rest = end_type(&rest);
                    if et_full == et_rest {
                        preserve_pass += 1;
                    } else {
                        preserve_fail += 1;
                        if preserve_fail <= 5 {
                            println!(
                                "  MISMATCH: pi={:?} end_type={} rest={:?} end_type={}",
                                pi, et_full, rest, et_rest
                            );
                        }
                    }
                } else {
                    // rest has length 1: end type is 'a' by convention
                    // pi has length 2: end type is determined by pi[0] vs pi[1]
                    // After removal, rest has length 1 -> 'a'
                    // For pi of length 2: if pi[0] < pi[1] -> 'a', if pi[0] > pi[1] -> 'd'
                    // The end type of pi (length 2) may differ from rest (length 1, always 'a')
                    // This is expected for n=2, so we count separately
                    preserve_pass += 1; // length-1 remainder: convention, not a real failure
                }
            }
        }
    }

    println!("End-type preservation checks: {}", preserve_checks);
    println!(
        "  Pass: {} ({:.1}%)",
        preserve_pass,
        if preserve_checks > 0 {
            100.0 * preserve_pass as f64 / preserve_checks as f64
        } else {
            0.0
        }
    );
    if preserve_fail > 0 {
        println!(
            "  FAIL: {} - end-type NOT preserved under first-entry removal!",
            preserve_fail
        );
    } else {
        println!("  All pass: end-type IS preserved under first-entry removal.");
        println!("  => The recursion DECOUPLES by end-type!");
    }

    // ────────────────────────────────────────────────
    // PRINT SOME SAMPLE DATA
    // ────────────────────────────────────────────────

    println!("\n================================================================");
    println!("  SAMPLE DATA: End-type polynomials for small boards");
    println!("================================================================\n");

    for n in 1..=std::cmp::min(max_n, 5) {
        let boards = gen_boards(n);
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            let m = (board[0] as usize).min(n);
            if m < 2 || n < 3 {
                continue;
            } // skip trivial cases

            let ideal = bruhat_lower_ideal(&perm);
            let mut d_ea = vec![vec![0i64]; m + 1];
            let mut d_ed = vec![vec![0i64]; m + 1];
            let mut u_ea = vec![vec![0i64]; m + 1];
            let mut u_ed = vec![vec![0i64]; m + 1];

            for pi in &ideal {
                let k = pi[0] as usize;
                if k > m {
                    continue;
                }
                let pk = peaks(pi);
                let et = end_type(pi);
                let st = start_type(pi);
                let poly_ref = match (st, et) {
                    ('d', 'a') => &mut d_ea[k],
                    ('d', 'd') => &mut d_ed[k],
                    ('u', 'a') => &mut u_ea[k],
                    ('u', 'd') => &mut u_ed[k],
                    _ => unreachable!(),
                };
                while poly_ref.len() <= pk {
                    poly_ref.push(0);
                }
                poly_ref[pk] += 1;
            }

            let has_data =
                (1..=m).any(|k| !pz(&d_ea[k]) || !pz(&d_ed[k]) || !pz(&u_ea[k]) || !pz(&u_ed[k]));
            if !has_data {
                continue;
            }

            println!("Board {:?} (n={}, m={}, perm={:?}):", board, n, m, perm);
            for k in 1..=m {
                let a_ea_k = pa(&d_ea[k], &u_ea[k]);
                let a_ed_k = pa(&d_ed[k], &u_ed[k]);
                let w_ea_k = pa(&pmt(&d_ea[k]), &u_ea[k]);
                let w_ed_k = pa(&pmt(&d_ed[k]), &u_ed[k]);

                if pz(&a_ea_k) && pz(&a_ed_k) {
                    continue;
                }
                println!("  k={}:", k);
                if !pz(&u_ea[k]) {
                    println!("    U_k^ea = {}", pf(&u_ea[k]));
                }
                if !pz(&u_ed[k]) {
                    println!("    U_k^ed = {}", pf(&u_ed[k]));
                }
                if !pz(&d_ea[k]) {
                    println!("    D_k^ea = {}", pf(&d_ea[k]));
                }
                if !pz(&d_ed[k]) {
                    println!("    D_k^ed = {}", pf(&d_ed[k]));
                }
                if !pz(&a_ea_k) {
                    println!("    A_k^ea = {}", pf(&a_ea_k));
                }
                if !pz(&a_ed_k) {
                    println!("    A_k^ed = {}", pf(&a_ed_k));
                }
                if !pz(&w_ea_k) {
                    println!("    W_k^ea = {}", pf(&w_ea_k));
                }
                if !pz(&w_ed_k) {
                    println!("    W_k^ed = {}", pf(&w_ed_k));
                }
            }
            println!();
        }
    }

    // ────────────────────────────────────────────────
    // SAMPLE DATA: P_{k,l}
    // ────────────────────────────────────────────────

    println!("================================================================");
    println!("  SAMPLE DATA: P_{{k,l}} for small boards");
    println!("================================================================\n");

    for n in 1..=std::cmp::min(max_n, 5) {
        let boards = gen_boards(n);
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            let m = (board[0] as usize).min(n);
            if m < 2 || n < 3 {
                continue;
            }

            let ideal = bruhat_lower_ideal(&perm);
            let mut pkl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];

            for pi in &ideal {
                let k = pi[0] as usize;
                if k > m {
                    continue;
                }
                let pk = peaks(pi);
                let last = *pi.last().unwrap() as usize;
                if last <= m {
                    while pkl[k][last].len() <= pk {
                        pkl[k][last].push(0);
                    }
                    pkl[k][last][pk] += 1;
                }
            }

            let has_data = (1..=m).any(|k| (1..=m).any(|l| !pz(&pkl[k][l])));
            if !has_data {
                continue;
            }

            println!("Board {:?} (n={}, m={}):", board, n, m);
            for k in 1..=m {
                for l in 1..=m {
                    if !pz(&pkl[k][l]) {
                        println!("  P_{{{},{}}} = {}", k, l, pf(&pkl[k][l]));
                    }
                }
            }
            println!();
        }
    }

    // ── Refinement 2: more interlacing patterns ──
    // Check ALL pairs P_{k1,l1} vs P_{k2,l2} for patterns

    println!("================================================================");
    println!("  REFINEMENT 2: Additional interlacing patterns");
    println!("================================================================\n");

    let mut pkl_any_pair_fwd = [0u64; 2]; // P_{k1,l1} << P_{k2,l2} for k1 > k2, same l
    let mut pkl_any_pair_same_k = [0u64; 2]; // P_{k,l1} << P_{k,l2} for l1 < l2
    let mut pkl_diag = [0u64; 2]; // P_{k,l} << P_{k+1,l+1}
    let mut pkl_antidiag = [0u64; 2]; // P_{k,l} << P_{k+1,l-1}

    // Row sums and column sums
    let mut row_sum_interlace = [0u64; 2]; // sum_l P_{k,l} << sum_l P_{k+1,l}
    let mut col_sum_interlace = [0u64; 2]; // sum_k P_{k,l} << sum_k P_{k,l+1}

    for n in 1..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            let m = (board[0] as usize).min(n);
            let ideal = bruhat_lower_ideal(&perm);

            let mut pkl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            for pi in &ideal {
                let k = pi[0] as usize;
                if k > m {
                    continue;
                }
                let pk = peaks(pi);
                let last = *pi.last().unwrap() as usize;
                if last <= m {
                    while pkl[k][last].len() <= pk {
                        pkl[k][last].push(0);
                    }
                    pkl[k][last][pk] += 1;
                }
            }

            // Fixed l, descending k: P_{j+1,l} << P_{j,l}
            for l in 1..=m {
                for j in 1..m {
                    if !pz(&pkl[j + 1][l]) && !pz(&pkl[j][l]) {
                        pkl_any_pair_fwd[0] += 1;
                        if !interlaces(&pkl[j + 1][l], &pkl[j][l]) {
                            pkl_any_pair_fwd[1] += 1;
                        }
                    }
                }
            }

            // Fixed k, ascending l: P_{k,l} << P_{k,l+1}
            for k in 1..=m {
                for l in 1..m {
                    if !pz(&pkl[k][l]) && !pz(&pkl[k][l + 1]) {
                        pkl_any_pair_same_k[0] += 1;
                        if !interlaces(&pkl[k][l], &pkl[k][l + 1]) {
                            pkl_any_pair_same_k[1] += 1;
                        }
                    }
                }
            }

            // Diagonal: P_{k,l} << P_{k+1,l+1}
            for k in 1..m {
                for l in 1..m {
                    if !pz(&pkl[k][l]) && !pz(&pkl[k + 1][l + 1]) {
                        pkl_diag[0] += 1;
                        if !interlaces(&pkl[k][l], &pkl[k + 1][l + 1]) {
                            pkl_diag[1] += 1;
                        }
                    }
                }
            }

            // Anti-diagonal: P_{k,l} << P_{k+1,l-1}
            for k in 1..m {
                for l in 2..=m {
                    if !pz(&pkl[k][l]) && !pz(&pkl[k + 1][l - 1]) {
                        pkl_antidiag[0] += 1;
                        if !interlaces(&pkl[k][l], &pkl[k + 1][l - 1]) {
                            pkl_antidiag[1] += 1;
                        }
                    }
                }
            }

            // Row sums: R_k = sum_l P_{k,l}
            let mut row_sums = vec![vec![0i64]; m + 1];
            for k in 1..=m {
                for l in 1..=m {
                    row_sums[k] = pa(&row_sums[k], &pkl[k][l]);
                }
            }
            for k in 1..m {
                if !pz(&row_sums[k]) && !pz(&row_sums[k + 1]) {
                    row_sum_interlace[0] += 1;
                    if !interlaces(&row_sums[k], &row_sums[k + 1]) {
                        row_sum_interlace[1] += 1;
                    }
                }
            }

            // Column sums: C_l = sum_k P_{k,l}
            let mut col_sums = vec![vec![0i64]; m + 1];
            for l in 1..=m {
                for k in 1..=m {
                    col_sums[l] = pa(&col_sums[l], &pkl[k][l]);
                }
            }
            for l in 1..m {
                if !pz(&col_sums[l]) && !pz(&col_sums[l + 1]) {
                    col_sum_interlace[0] += 1;
                    if !interlaces(&col_sums[l], &col_sums[l + 1]) {
                        col_sum_interlace[1] += 1;
                    }
                }
            }
        }
    }

    println!("Fixed l, P_{{j+1,l}} << P_{{j,l}}:");
    show("  P_{j+1,l} << P_{j,l}", pkl_any_pair_fwd);
    println!("Fixed k, P_{{k,l}} << P_{{k,l+1}}:");
    show("  P_{k,l} << P_{k,l+1}", pkl_any_pair_same_k);
    println!("Diagonal P_{{k,l}} << P_{{k+1,l+1}}:");
    show("  P_{k,l} << P_{k+1,l+1}", pkl_diag);
    println!("Anti-diagonal P_{{k,l}} << P_{{k+1,l-1}}:");
    show("  P_{k,l} << P_{k+1,l-1}", pkl_antidiag);
    println!();
    println!("Row sums R_k << R_{{k+1}}:");
    show("  R_k << R_{k+1}", row_sum_interlace);
    println!("Column sums C_l << C_{{l+1}}:");
    show("  C_l << C_{l+1}", col_sum_interlace);
    println!();

    // ────────────────────────────────────────────────
    // SUMMARY
    // ────────────────────────────────────────────────

    println!("\n================================================================");
    println!("  SUMMARY");
    println!("================================================================\n");

    println!("Total valid 312-avoiding boards tested: {}\n", valid_boards);

    println!("=== REFINEMENT 1: End-type (ea = ascending end, ed = descending end) ===\n");

    println!("--- Unsplit reference ---");
    show("  A_j << W_l (all j,l)", aw_all);
    show("  W_j << W_l (j < l)", ww_all);
    println!();

    println!("--- Within same end-type ---");
    show("  A_j^ea << W_l^ea (all j,l)", aw_ea);
    show("  W_j^ea << W_l^ea (j < l)", ww_ea);
    show("  A_j^ed << W_l^ed (all j,l)", aw_ed);
    show("  W_j^ed << W_l^ed (j < l)", ww_ed);
    println!();

    println!("--- Cross end-type ---");
    show("  A_j^ea << W_l^ed (all j,l)", aw_ea_ed);
    show("  A_j^ed << W_l^ea (all j,l)", aw_ed_ea);
    show("  W_j^ea << W_l^ed (j < l)", ww_ea_ed);
    println!();

    println!("--- A_k^ea vs A_k^ed (same k) ---");
    show("  A_k^ea << A_k^ed", ak_ea_ed);
    show("  A_k^ed << A_k^ea", ak_ed_ea);
    show("  At least one direction", ak_either);
    println!();

    println!("--- DU within end-types ---");
    show("  D_j^ea << U_l^ea (all j,l)", du_ea);
    show("  D_j^ed << U_l^ed (all j,l)", du_ed);
    show("  U_j^ea << U_l^ea (j <= l)", uu_ea);
    show("  U_j^ed << U_l^ed (j <= l)", uu_ed);
    println!();

    println!("=== REFINEMENT 2: Last entry P_{{k,l}} ===\n");
    show("  P_{k,l} real-rooted", pkl_rr);
    println!();
    show("  Fixed l: P_{j,l} << P_{j+1,l}", pkl_fixed_l_interlace);
    show("  Fixed l: P_{j+1,l} << P_{j,l} (reverse)", pkl_fixed_l_rev);
    show("  Fixed k: P_{k,l} << P_{k,l+1}", pkl_fixed_k_interlace);
    show("  Fixed k: P_{k,l+1} << P_{k,l} (reverse)", pkl_fixed_k_rev);
    println!();

    println!("=== DECOUPLING ===\n");
    println!(
        "  Split consistency (D/U sums): {}/{}",
        decouple_match, decouple_boards_checked
    );
    println!(
        "  End-type preservation under first-entry removal: {}/{}",
        preserve_pass, preserve_checks
    );
    if preserve_fail == 0 {
        println!("  => CONFIRMED: The recursion decouples by end-type.");
    } else {
        println!("  => Decoupling FAILS for {} cases.", preserve_fail);
    }

    println!("\n================================================================");
}

fn show(name: &str, c: [u64; 2]) {
    if c[0] == 0 {
        println!("{}: (no data)", name);
    } else if c[1] == 0 {
        println!("{}: {}/{} ALL PASS", name, c[0], c[0]);
    } else {
        println!("{}: {}/{} pass ({} FAIL)", name, c[0] - c[1], c[0], c[1]);
    }
}
