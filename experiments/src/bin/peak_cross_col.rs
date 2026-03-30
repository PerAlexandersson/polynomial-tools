//! Comprehensive cross-column interlacing tests for the cross-group proof
//! of Conjecture 13.
//!
//! For each 312-avoiding Ferrers board with n <= max_n, computes D_{j,l},
//! U_{j,l}, A_{j,l}, W_{j,l} for all j, l and tests:
//!
//! 1. Ground truth: D_{k'}^+(col l') << U_{j'}^+(col l'-1) cross-group at lambda+
//! 2. Cross-column DU: D_{j,l} << U_{j',l-1}
//! 3. Cross-column D<<W: D_{j,l} << W_{j',l-1}
//! 4. Cross-column A<<W: A_{j,l'} << W_{j',l'-1} (all and forward-only)
//! 5. Shift-lemma across columns: h_cross << f_cross at s=1
//! 6. Key decomposition: h = Σ_{i≥j} D_{i,l'-1} << f_s at s=1
//!    where f_s = Σ_{i<k} A_{i,l'} + s·Σ_{i≥j} A_{i,l'-1}
//! 7. Component tests: h << each piece separately
//!
//! Usage: cargo run --release --bin peak_cross_col -- 8

use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly};
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
        println!("  {}: {}/{} ALL PASS  *** POTENTIAL PROOF PATH ***", name, c[0], c[0]);
    } else {
        println!("  {}: {}/{} pass ({} FAIL)", name, c[0] - c[1], c[0], c[1]);
    }
}

struct FailDetail {
    board: Vec<u8>,
    desc: String,
    f_poly: Vec<i64>,
    g_poly: Vec<i64>,
}

fn record_fail(fails: &mut Vec<FailDetail>, board: &[u8], desc: String, f: &[i64], g: &[i64]) {
    if fails.len() < 5 {
        fails.push(FailDetail {
            board: board.to_vec(),
            desc,
            f_poly: pt(f),
            g_poly: pt(g),
        });
    }
}

fn print_fails(label: &str, fails: &[FailDetail]) {
    if fails.is_empty() { return; }
    println!("\n  --- First {} {} failures ---", fails.len(), label);
    for (i, f) in fails.iter().enumerate() {
        println!("    #{}: board={:?}, {}", i+1, f.board, f.desc);
        println!("       f = {} (deg {:?})", format_poly(&f.f_poly), pdeg(&f.f_poly));
        println!("       g = {} (deg {:?})", format_poly(&f.g_poly), pdeg(&f.g_poly));
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
    let max_n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(8);
    println!("================================================================");
    println!("  Cross-column interlacing for cross-group proof (Conjecture 13)");
    println!("  312-avoiding Ferrers boards, n <= {}", max_n);
    println!("================================================================\n");

    // ── Test 1: Ground truth cross-group DU at lambda+ ──
    let mut gt_xg_du = [0u64; 2];
    let mut gt_xg_du_fails: Vec<FailDetail> = Vec::new();

    // ── Test 2: Cross-column DU: D_{j,l} << U_{j',l-1} ──
    let mut x_du = [0u64; 2];
    let mut x_du_fails: Vec<FailDetail> = Vec::new();

    // ── Test 3: Cross-column D<<W: D_{j,l} << W_{j',l-1} ──
    let mut x_dw = [0u64; 2];
    let mut x_dw_fails: Vec<FailDetail> = Vec::new();

    // ── Test 4: Cross-column A<<W: A_{j,l'} << W_{j',l'-1} ──
    let mut x_aw_all = [0u64; 2];
    let mut x_aw_all_fails: Vec<FailDetail> = Vec::new();
    let mut x_aw_fwd = [0u64; 2]; // j <= j' only
    let mut x_aw_fwd_fails: Vec<FailDetail> = Vec::new();

    // ── Test 5: Shift-lemma across columns ──
    // h_cross = Σ_{i≥j'} D_{i,l'-1}, f_cross = Σ_{i<k'} A_{i,l'} + Σ_{i≥j'} A_{i,l'-1}
    let mut t5_h_f = [0u64; 2];
    let mut t5_h_f_fails: Vec<FailDetail> = Vec::new();

    // ── Test 6: Key decomposition h << f_s at s=1 ──
    // h = Σ_{i≥j} D_{i,l'-1}, f_1 = Σ_{i<k} A_{i,l'} + Σ_{i≥j} A_{i,l'-1}
    let mut t6_h_f1 = [0u64; 2];
    let mut t6_h_f1_fails: Vec<FailDetail> = Vec::new();

    // ── Test 7: Component tests ──
    // 7a: h << Σ_{i<k} A_{i,l'} (h from col l'-1, A's from col l')
    let mut t7a = [0u64; 2];
    let mut t7a_fails: Vec<FailDetail> = Vec::new();
    // 7b: h << Σ_{i≥j} A_{i,l'-1} (same column)
    let mut t7b = [0u64; 2];
    let mut t7b_fails: Vec<FailDetail> = Vec::new();
    // 7c: h << Σ_{i≥j} U_{i,l'-1} (same column, Lemma lem:h_interl_U)
    let mut t7c = [0u64; 2];
    let mut t7c_fails: Vec<FailDetail> = Vec::new();

    let mut valid_boards = 0u64;

    for n in 1..=max_n {
        let boards = gen_boards(n);
        let mut board_count = 0u64;
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) { continue; }
            valid_boards += 1;
            board_count += 1;

            let m = board[0] as usize;
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

            // ── Precompute prefix/suffix sums for each column ──
            // S_{k,col} = Σ_{i=1}^{k-1} A_{i,col} (prefix sum of A, exclusive)
            // T_{k,col} = Σ_{i=k}^{m} W_{i,col}   (suffix sum of W)
            // PD_{k,col} = Σ_{i=k}^{m} D_{i,col}  (suffix sum of D)
            // PA_{k,col} = Σ_{i=k}^{m} A_{i,col}  (suffix sum of A)
            // PU_{k,col} = Σ_{i=k}^{m} U_{i,col}  (suffix sum of U)
            let sz = m + 2;
            let mut s_col: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; sz]; m + 1];  // S[col][k]
            let mut t_col: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; sz]; m + 1];  // T[col][k]
            let mut pd_col: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; sz]; m + 1]; // PD[col][k]
            let mut pa_col: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; sz]; m + 1]; // PA[col][k] (suffix of A)
            let mut pu_col: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; sz]; m + 1]; // PU[col][k] (suffix of U)

            for col in 1..=m {
                // Prefix sum of A: S[col][1] = 0, S[col][k] = Σ_{i<k} A_{i,col}
                s_col[col][1] = vec![0i64];
                for k in 1..=m {
                    s_col[col][k + 1] = pa(&s_col[col][k], &a_jl[k][col]);
                }

                // Suffix sums (from m down to 1)
                t_col[col][m + 1] = vec![0i64];
                pd_col[col][m + 1] = vec![0i64];
                pa_col[col][m + 1] = vec![0i64];
                pu_col[col][m + 1] = vec![0i64];
                for k in (1..=m).rev() {
                    t_col[col][k] = pa(&t_col[col][k + 1], &w_jl[k][col]);
                    pd_col[col][k] = pa(&pd_col[col][k + 1], &d_jl[k][col]);
                    pa_col[col][k] = pa(&pa_col[col][k + 1], &a_jl[k][col]);
                    pu_col[col][k] = pa(&pu_col[col][k + 1], &u_jl[k][col]);
                }
            }

            // Helper: get S (prefix sum of A, Σ_{i<k} A_{i,col})
            let get_s = |col: usize, k: usize| -> Vec<i64> {
                if col < 1 || col > m { return vec![0i64]; }
                if k >= sz { s_col[col][m + 1].clone() } else { s_col[col][k].clone() }
            };
            // Helper: get T (suffix sum of W, Σ_{i>=k} W_{i,col})
            let get_t = |col: usize, k: usize| -> Vec<i64> {
                if col < 1 || col > m { return vec![0i64]; }
                if k >= sz { vec![0i64] } else { t_col[col][k].clone() }
            };
            // Helper: get PD (suffix sum of D, Σ_{i>=k} D_{i,col})
            let get_pd = |col: usize, k: usize| -> Vec<i64> {
                if col < 1 || col > m { return vec![0i64]; }
                if k >= sz { vec![0i64] } else { pd_col[col][k].clone() }
            };
            // Helper: get PA (suffix sum of A, Σ_{i>=k} A_{i,col})
            let get_pa = |col: usize, k: usize| -> Vec<i64> {
                if col < 1 || col > m { return vec![0i64]; }
                if k >= sz { vec![0i64] } else { pa_col[col][k].clone() }
            };
            // Helper: get PU (suffix sum of U, Σ_{i>=k} U_{i,col})
            let get_pu = |col: usize, k: usize| -> Vec<i64> {
                if col < 1 || col > m { return vec![0i64]; }
                if k >= sz { vec![0i64] } else { pu_col[col][k].clone() }
            };

            // ── Test 2: Cross-column DU at lambda level ──
            for l in 2..=m {
                for j in 1..=m {
                    for jp in 1..=m {
                        if !pz(&d_jl[j][l]) && !pz(&u_jl[jp][l-1]) {
                            x_du[0] += 1;
                            if !interlaces(&d_jl[j][l], &u_jl[jp][l-1]) {
                                x_du[1] += 1;
                                record_fail(&mut x_du_fails, board,
                                    format!("l={}, j={}, j'={}", l, j, jp),
                                    &d_jl[j][l], &u_jl[jp][l-1]);
                            }
                        }
                    }
                }
            }

            // ── Test 3: Cross-column D<<W at lambda level ──
            for l in 2..=m {
                for j in 1..=m {
                    for jp in 1..=m {
                        if !pz(&d_jl[j][l]) && !pz(&w_jl[jp][l-1]) {
                            x_dw[0] += 1;
                            if !interlaces(&d_jl[j][l], &w_jl[jp][l-1]) {
                                x_dw[1] += 1;
                                record_fail(&mut x_dw_fails, board,
                                    format!("l={}, j={}, j'={}", l, j, jp),
                                    &d_jl[j][l], &w_jl[jp][l-1]);
                            }
                        }
                    }
                }
            }

            // ── Test 4: Cross-column A<<W ──
            for l in 2..=m {
                for j in 1..=m {
                    for jp in 1..=m {
                        if !pz(&a_jl[j][l]) && !pz(&w_jl[jp][l-1]) {
                            x_aw_all[0] += 1;
                            if !interlaces(&a_jl[j][l], &w_jl[jp][l-1]) {
                                x_aw_all[1] += 1;
                                record_fail(&mut x_aw_all_fails, board,
                                    format!("l={}, j={}, j'={}", l, j, jp),
                                    &a_jl[j][l], &w_jl[jp][l-1]);
                            }
                            // Forward only: j <= j'
                            if j <= jp {
                                x_aw_fwd[0] += 1;
                                if !interlaces(&a_jl[j][l], &w_jl[jp][l-1]) {
                                    x_aw_fwd[1] += 1;
                                    record_fail(&mut x_aw_fwd_fails, board,
                                        format!("l={}, j={}, j'={} (fwd)", l, j, jp),
                                        &a_jl[j][l], &w_jl[jp][l-1]);
                                }
                            }
                        }
                    }
                }
            }

            // ── Tests 1, 5, 6, 7: at lambda+ level, cross-group pairs ──
            let mp = m + 1;

            for lp in 2..=mp {
                // Upper group: k' > l', uses col = l'
                // Lower group: k' < l', uses col = l'-1

                // Enumerate all cross-group pairs: k' from upper, j' from lower
                for kp in (lp + 1)..=mp {
                    // k' > l' => upper group, col_upper = l'
                    let col_upper = lp;
                    if col_upper < 1 || col_upper > m { continue; }

                    // D_{k'}^+ = Σ_{i<k'} A_{i,col_upper} = S[col_upper][k']
                    let d_plus_upper = get_s(col_upper, kp);

                    for jp in 1..lp {
                        // j' < l' => lower group, col_lower = l'-1
                        let col_lower = lp - 1;
                        if col_lower < 1 || col_lower > m { continue; }

                        // U_{j'}^+ = Σ_{i>=j'} W_{i,col_lower} = T[col_lower][j']
                        let u_plus_lower = get_t(col_lower, jp);

                        // ── Test 1: Ground truth D^+ << U^+ ──
                        if !pz(&d_plus_upper) && !pz(&u_plus_lower) {
                            gt_xg_du[0] += 1;
                            if !interlaces(&d_plus_upper, &u_plus_lower) {
                                gt_xg_du[1] += 1;
                                record_fail(&mut gt_xg_du_fails, board,
                                    format!("l'={}, k'={} (upper,col{}), j'={} (lower,col{})",
                                            lp, kp, col_upper, jp, col_lower),
                                    &d_plus_upper, &u_plus_lower);
                            }
                        }

                        // ── Test 5: h_cross << f_cross ──
                        // h_cross = Σ_{i>=j'} D_{i,col_lower}
                        let h_cross = get_pd(col_lower, jp);
                        // f_cross = Σ_{i<k'} A_{i,col_upper} + Σ_{i>=j'} A_{i,col_lower}
                        let f_cross = pa(&get_s(col_upper, kp), &get_pa(col_lower, jp));
                        if !pz(&h_cross) && !pz(&f_cross) {
                            t5_h_f[0] += 1;
                            if !interlaces(&h_cross, &f_cross) {
                                t5_h_f[1] += 1;
                                record_fail(&mut t5_h_f_fails, board,
                                    format!("l'={}, k'={}, j'={}", lp, kp, jp),
                                    &h_cross, &f_cross);
                            }
                        }

                        // ── Test 6: h << f_1 (same as Test 5 by definition at s=1) ──
                        // h = Σ_{i>=j'} D_{i,col_lower}
                        // f_1 = Σ_{i<k'} A_{i,col_upper} + 1*Σ_{i>=j'} A_{i,col_lower}
                        // This is the same as Test 5 when s=1, confirming.
                        let h6 = get_pd(col_lower, jp);
                        let f6 = pa(&get_s(col_upper, kp), &get_pa(col_lower, jp));
                        if !pz(&h6) && !pz(&f6) {
                            t6_h_f1[0] += 1;
                            if !interlaces(&h6, &f6) {
                                t6_h_f1[1] += 1;
                                record_fail(&mut t6_h_f1_fails, board,
                                    format!("l'={}, k'={}, j'={}", lp, kp, jp),
                                    &h6, &f6);
                            }
                        }

                        // ── Test 7a: h << Σ_{i<k'} A_{i,col_upper} ──
                        let h7 = get_pd(col_lower, jp);
                        let g7a = get_s(col_upper, kp);
                        if !pz(&h7) && !pz(&g7a) {
                            t7a[0] += 1;
                            if !interlaces(&h7, &g7a) {
                                t7a[1] += 1;
                                record_fail(&mut t7a_fails, board,
                                    format!("l'={}, k'={}, j'={}", lp, kp, jp),
                                    &h7, &g7a);
                            }
                        }

                        // ── Test 7b: h << Σ_{i>=j'} A_{i,col_lower} ──
                        let g7b = get_pa(col_lower, jp);
                        if !pz(&h7) && !pz(&g7b) {
                            t7b[0] += 1;
                            if !interlaces(&h7, &g7b) {
                                t7b[1] += 1;
                                record_fail(&mut t7b_fails, board,
                                    format!("l'={}, k'={}, j'={}", lp, kp, jp),
                                    &h7, &g7b);
                            }
                        }

                        // ── Test 7c: h << Σ_{i>=j'} U_{i,col_lower} ──
                        let g7c = get_pu(col_lower, jp);
                        if !pz(&h7) && !pz(&g7c) {
                            t7c[0] += 1;
                            if !interlaces(&h7, &g7c) {
                                t7c[1] += 1;
                                record_fail(&mut t7c_fails, board,
                                    format!("l'={}, k'={}, j'={}", lp, kp, jp),
                                    &h7, &g7c);
                            }
                        }
                    }
                }
            }
        }
        println!("n={}: {} valid boards (cumulative: {})", n, board_count, valid_boards);
    }

    println!("\n================================================================");
    println!("  TEST 1: GROUND TRUTH - Cross-group DU at lambda+");
    println!("  D_{{k'}}^+(col l') << U_{{j'}}^+(col l'-1)  [k'>l', j'<l']");
    println!("================================================================");
    show("Cross-group D^+ << U^+", gt_xg_du);
    print_fails("Ground-truth cross-group DU", &gt_xg_du_fails);

    println!("\n================================================================");
    println!("  TEST 2: Cross-column DU at lambda level");
    println!("  D_{{j,l}} << U_{{j',l-1}}");
    println!("================================================================");
    show("D_{j,l} << U_{j',l-1}", x_du);
    print_fails("Cross-column DU", &x_du_fails);

    println!("\n================================================================");
    println!("  TEST 3: Cross-column D<<W at lambda level");
    println!("  D_{{j,l}} << W_{{j',l-1}}");
    println!("================================================================");
    show("D_{j,l} << W_{j',l-1}", x_dw);
    print_fails("Cross-column DW", &x_dw_fails);

    println!("\n================================================================");
    println!("  TEST 4: Cross-column A<<W at lambda level");
    println!("  A_{{j,l'}} << W_{{j',l'-1}}");
    println!("================================================================");
    show("A_{j,l} << W_{j',l-1}  [all j,j']", x_aw_all);
    print_fails("Cross-column AW (all)", &x_aw_all_fails);
    show("A_{j,l} << W_{j',l-1}  [j<=j' only]", x_aw_fwd);
    print_fails("Cross-column AW (fwd)", &x_aw_fwd_fails);

    println!("\n================================================================");
    println!("  TEST 5: Shift-lemma across columns");
    println!("  h_cross = sum_{{i>=j'}} D_{{i,l'-1}}");
    println!("  f_cross = sum_{{i<k'}} A_{{i,l'}} + sum_{{i>=j'}} A_{{i,l'-1}}");
    println!("  Test: h_cross << f_cross");
    println!("================================================================");
    show("h_cross << f_cross", t5_h_f);
    print_fails("Shift-lemma cross", &t5_h_f_fails);

    println!("\n================================================================");
    println!("  TEST 6: Key decomposition h << f_s at s=1");
    println!("  h = sum_{{i>=j}} D_{{i,l'-1}}");
    println!("  f_1 = sum_{{i<k}} A_{{i,l'}} + sum_{{i>=j}} A_{{i,l'-1}}");
    println!("  [Same as Test 5 - confirmation]");
    println!("================================================================");
    show("h << f_1", t6_h_f1);
    print_fails("Key decomp h<<f_1", &t6_h_f1_fails);

    println!("\n================================================================");
    println!("  TEST 7: Component tests (h from col l'-1)");
    println!("  h = sum_{{i>=j'}} D_{{i,l'-1}}");
    println!("================================================================");
    show("7a: h << sum_{{i<k'}} A_{{i,l'}}    [cross-col]", t7a);
    print_fails("7a: h << A-prefix(cross)", &t7a_fails);
    show("7b: h << sum_{{i>=j'}} A_{{i,l'-1}} [same col]", t7b);
    print_fails("7b: h << A-suffix(same)", &t7b_fails);
    show("7c: h << sum_{{i>=j'}} U_{{i,l'-1}} [same col, Lemma h<<U]", t7c);
    print_fails("7c: h << U-suffix(same)", &t7c_fails);

    println!("\n================================================================");
    println!("  SUMMARY");
    println!("================================================================");
    println!("  Total boards tested: {}", valid_boards);

    // Collect all results for summary
    let tests: Vec<(&str, [u64; 2])> = vec![
        ("T1: Ground-truth cross-group D+<<U+", gt_xg_du),
        ("T2: Cross-col D_{j,l}<<U_{j',l-1}", x_du),
        ("T3: Cross-col D_{j,l}<<W_{j',l-1}", x_dw),
        ("T4a: Cross-col A_{j,l}<<W_{j',l-1} all", x_aw_all),
        ("T4b: Cross-col A_{j,l}<<W_{j',l-1} j<=j'", x_aw_fwd),
        ("T5: h_cross << f_cross", t5_h_f),
        ("T6: h << f_1 (=T5, confirm)", t6_h_f1),
        ("T7a: h << A-prefix (cross-col)", t7a),
        ("T7b: h << A-suffix (same col)", t7b),
        ("T7c: h << U-suffix (same col)", t7c),
    ];

    println!("\n  Condition                              | Tested  | Pass    | Fail   | Status");
    println!("  -------|------|------|------|------");
    for (name, c) in &tests {
        if c[0] == 0 {
            println!("  {:42} | {:>7} | {:>7} | {:>6} | NO DATA", name, 0, 0, 0);
        } else {
            let status = if c[1] == 0 { "ALL PASS ***" } else { "HAS FAILURES" };
            println!("  {:42} | {:>7} | {:>7} | {:>6} | {}", name, c[0], c[0] - c[1], c[1], status);
        }
    }

    println!();
    println!("  Conditions with ALL PASS are potential proof paths.");
    println!("================================================================");
}
