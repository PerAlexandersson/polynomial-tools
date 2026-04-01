//! Bridge at k'=l' approach for P_{k,l} refinement of peak polynomials.
//!
//! For each 312-avoiding Ferrers board lambda with n <= max_n, computes:
//! 1. D_{j,l} and U_{j,l}: peak polys refined by first entry j AND last entry l,
//!    split by whether pi_1 > pi_2 (D) or pi_1 < pi_2 (U).
//! 2. The lambda+ level polynomials P_{k',l'}^+ for k' != l'.
//! 3. Upper group, lower group, and bridge interlacing conditions.
//! 4. Cross-column interlacing conditions.
//!
//! Usage: cargo run --release --bin peak_bridge -- 6

use combpoly::order::bruhat_lower_ideal;
use experiments::peak_utils::{
    board_to_perm, gen_boards, is_312_avoiding, pa, peak_count, pmt, pt, pz,
};
use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly, is_real_rooted};

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

fn show(name: &str, c: [u64; 2]) {
    if c[0] == 0 {
        println!("{}: (no data)", name);
    } else if c[1] == 0 {
        println!("{}: {}/{} ALL PASS", name, c[0], c[0]);
    } else {
        println!("{}: {}/{} pass ({} FAIL)", name, c[0] - c[1], c[0], c[1]);
    }
}

/// Classify the "start type" of a permutation:
/// - 'd' (down/descent) if pi_1 > pi_2
/// - 'u' (up/ascent) if pi_1 < pi_2 or n <= 1
fn start_type(pi: &[u8]) -> char {
    let n = pi.len();
    if n <= 1 {
        return 'u';
    }
    if pi[0] > pi[1] {
        'd'
    } else {
        'u'
    }
}

fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);
    println!("================================================================");
    println!("  BRIDGE AT k'=l': P_{{k,l}} refinement of peak polynomials");
    println!("  312-avoiding Ferrers boards, n <= {}", max_n);
    println!("================================================================\n");

    // ── Global counters ──

    // Part 1: D_{j,l} and U_{j,l} real-rootedness
    let mut djl_rr = [0u64; 2];
    let mut ujl_rr = [0u64; 2];
    let mut ajl_rr = [0u64; 2];
    let mut wjl_rr = [0u64; 2];

    // Part 3a: Upper group interlacing (k' > l', fixed l')
    let mut upper_group_pair = [0u64; 2];
    let mut upper_group_full = [0u64; 2];

    // Part 3b: Lower group interlacing (k' < l', fixed l')
    let mut lower_group_pair = [0u64; 2];
    let mut lower_group_full = [0u64; 2];

    // Part 3c: Bridge condition P_{l'+1,l'}^+ << P_{l'-1,l'}^+
    let mut bridge_cond = [0u64; 2];

    // Part 4: Bridge decomposition checks
    let mut bridge_s_upper_t_lower = [0u64; 2]; // S_{l'+1}^{(l')} << T_{l'-1}^{(l'-1)}
    let mut bridge_t_upper_t_lower = [0u64; 2]; // T_{l'+1}^{(l')} << T_{l'-1}^{(l'-1)}
    let mut bridge_s_upper_s_lower = [0u64; 2]; // S_{l'+1}^{(l')} << S_{l'-1}^{(l'-1)}
    let mut bridge_t_upper_s_lower = [0u64; 2]; // T_{l'+1}^{(l')} << S_{l'-1}^{(l'-1)}

    // Part 5: Full column interlacing (upper + lower combined)
    let mut full_col_pair = [0u64; 2];
    let mut full_col_seq = [0u64; 2];

    // Part 6: Cross-column conditions
    let mut cross_a_a = [0u64; 2]; // A_{j,l} << A_{j,l-1}
    let mut cross_w_w = [0u64; 2]; // W_{j,l} << W_{j,l-1}
    let mut cross_a_w = [0u64; 2]; // A_{j,l} << W_{j',l-1}
    let mut cross_a_a_rev = [0u64; 2]; // A_{j,l-1} << A_{j,l} (reverse)
    let mut cross_w_w_rev = [0u64; 2]; // W_{j,l-1} << W_{j,l} (reverse)

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
            // PART 1: Compute D_{j,l} and U_{j,l}
            // ────────────────────────────────────────────────
            // D_{j,l}(t) = sum over pi in ideal with pi_1=j, pi_n=l, pi_1>pi_2: t^{peaks(pi)}
            // U_{j,l}(t) = sum over pi in ideal with pi_1=j, pi_n=l, pi_1<pi_2: t^{peaks(pi)}
            // Index: d_jl[j][l], u_jl[j][l] for j,l in 1..=m

            let mut d_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            let mut u_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];

            for pi in &ideal {
                let j = pi[0] as usize;
                if j > m {
                    continue;
                }
                let l = *pi.last().unwrap() as usize;
                if l > m {
                    continue;
                }
                let pk = peak_count(pi);
                let st = start_type(pi);

                let poly = if st == 'd' {
                    &mut d_jl[j][l]
                } else {
                    &mut u_jl[j][l]
                };
                while poly.len() <= pk {
                    poly.push(0);
                }
                poly[pk] += 1;
            }

            // Compute A_{j,l} = D_{j,l} + U_{j,l}
            // Compute W_{j,l} = t*D_{j,l} + U_{j,l}
            let mut a_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            let mut w_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];

            for j in 1..=m {
                for l in 1..=m {
                    a_jl[j][l] = pa(&d_jl[j][l], &u_jl[j][l]);
                    w_jl[j][l] = pa(&pmt(&d_jl[j][l]), &u_jl[j][l]);
                }
            }

            // Real-rootedness checks
            for j in 1..=m {
                for l in 1..=m {
                    if !pz(&d_jl[j][l]) {
                        djl_rr[0] += 1;
                        if !is_real_rooted(&d_jl[j][l]) {
                            djl_rr[1] += 1;
                        }
                    }
                    if !pz(&u_jl[j][l]) {
                        ujl_rr[0] += 1;
                        if !is_real_rooted(&u_jl[j][l]) {
                            ujl_rr[1] += 1;
                        }
                    }
                    if !pz(&a_jl[j][l]) {
                        ajl_rr[0] += 1;
                        if !is_real_rooted(&a_jl[j][l]) {
                            ajl_rr[1] += 1;
                        }
                    }
                    if !pz(&w_jl[j][l]) {
                        wjl_rr[0] += 1;
                        if !is_real_rooted(&w_jl[j][l]) {
                            wjl_rr[1] += 1;
                        }
                    }
                }
            }

            // ────────────────────────────────────────────────
            // PART 2: Compute P_{k',l'}^+ at the lambda+ level
            // ────────────────────────────────────────────────
            // For k' != l':
            //   If k' > l': use input column l' (polynomials A_{j,l'}, W_{j,l'})
            //     P_{k',l'}^+ = sum_{j<k'} A_{j,l'} + sum_{j>=k'} W_{j,l'}
            //   If k' < l': use input column l'-1 (polynomials A_{j,l'-1}, W_{j,l'-1})
            //     P_{k',l'}^+ = sum_{j<k'} A_{j,l'-1} + sum_{j>=k'} W_{j,l'-1}
            //
            // The lambda+ level has m' = m+1 rows, so k' ranges 1..=m+1, l' ranges 1..=m+1
            // But the input columns are indexed 1..=m.

            let mp = m + 1; // m' = m+1

            // We'll also compute the S and T decomposition of each P_{k',l'}^+
            // S_{k'}^{(col)} = sum_{j<k'} A_{j,col}
            // T_{k'}^{(col)} = sum_{j>=k'} W_{j,col}
            // P_{k',l'}^+ = S_{k'}^{(col)} + T_{k'}^{(col)}

            // Precompute prefix sums of A_{j,col} and suffix sums of W_{j,col} for each col
            // S_{k}^{(col)} = sum_{j=1}^{k-1} A_{j,col}  (so S_1 = 0, S_2 = A_1, ...)
            // T_{k}^{(col)} = sum_{j=k}^{m} W_{j,col}  (so T_{m+1} = 0, T_m = W_m, ...)

            // s_col[col][k] = S_k^{(col)} for k in 0..=m+1
            // t_col[col][k] = T_k^{(col)} for k in 0..=m+1
            let mut s_col: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; mp + 1]; m + 1];
            let mut t_col: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; mp + 1]; m + 1];

            for col in 1..=m {
                // Build prefix sums: S_1^{(col)} = 0, S_{k+1}^{(col)} = S_k + A_{k,col}
                for k in 1..=m {
                    s_col[col].push(vec![0i64]); // extend if needed
                    while s_col[col].len() <= k + 1 {
                        s_col[col].push(vec![0i64]);
                    }
                    s_col[col][k + 1] = pa(&s_col[col][k], &a_jl[k][col]);
                }
                // Build suffix sums: T_{m+1}^{(col)} = 0, T_k^{(col)} = T_{k+1} + W_{k,col}
                while t_col[col].len() <= mp {
                    t_col[col].push(vec![0i64]);
                }
                for k in (1..=m).rev() {
                    t_col[col][k] = pa(&t_col[col][k + 1], &w_jl[k][col]);
                }
            }

            // Now compute P_{k',l'}^+ for all k' in 1..=mp, l' in 1..=mp with k' != l'
            // p_plus[k'][l'] = P_{k',l'}^+
            let mut p_plus: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; mp + 1]; mp + 1];
            // Also store the S and T parts separately
            let mut s_part: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; mp + 1]; mp + 1];
            let mut t_part: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; mp + 1]; mp + 1];

            for kp in 1..=mp {
                for lp in 1..=mp {
                    if kp == lp {
                        continue;
                    }
                    if kp > lp {
                        // Use input column l' (if l' <= m)
                        let col = lp;
                        if col >= 1 && col <= m {
                            // S_{k'}^{(col)} needs k' <= m+1
                            let s = if kp <= m + 1 && kp < s_col[col].len() {
                                s_col[col][kp].clone()
                            } else {
                                // k' > m+1: S = sum of all A_{j,col}
                                s_col[col][m + 1].clone()
                            };
                            let t = if kp <= m && kp < t_col[col].len() {
                                t_col[col][kp].clone()
                            } else {
                                vec![0i64]
                            };
                            s_part[kp][lp] = s.clone();
                            t_part[kp][lp] = t.clone();
                            p_plus[kp][lp] = pa(&s, &t);
                        }
                    } else {
                        // kp < lp: use input column l'-1
                        let col = lp - 1;
                        if col >= 1 && col <= m {
                            let s = if kp <= m + 1 && kp < s_col[col].len() {
                                s_col[col][kp].clone()
                            } else {
                                s_col[col][m + 1].clone()
                            };
                            let t = if kp <= m && kp < t_col[col].len() {
                                t_col[col][kp].clone()
                            } else {
                                vec![0i64]
                            };
                            s_part[kp][lp] = s.clone();
                            t_part[kp][lp] = t.clone();
                            p_plus[kp][lp] = pa(&s, &t);
                        }
                    }
                }
            }

            // ────────────────────────────────────────────────
            // PART 3: Interlacing checks for fixed l'
            // ────────────────────────────────────────────────
            for lp in 1..=mp {
                // (a) UPPER GROUP: k' > l'
                // Sequence: P_{mp,l'}^+, P_{mp-1,l'}^+, ..., P_{l'+1,l'}^+
                // Check consecutive interlacing (descending k')
                let upper: Vec<usize> = (lp + 1..=mp).rev().collect();
                for i in 0..upper.len() {
                    if i + 1 < upper.len() {
                        let k1 = upper[i]; // larger k'
                        let k2 = upper[i + 1]; // smaller k'
                        if !pz(&p_plus[k1][lp]) && !pz(&p_plus[k2][lp]) {
                            upper_group_pair[0] += 1;
                            if !interlaces(&p_plus[k1][lp], &p_plus[k2][lp]) {
                                upper_group_pair[1] += 1;
                                if upper_group_pair[1] <= 3 {
                                    println!("FAIL upper pair: board={:?} l'={} P_{{{},{}}}^+ << P_{{{},{}}}^+",
                                             board, lp, k1, lp, k2, lp);
                                    println!("  P_{{{},{}}}^+ = {}", k1, lp, pf(&p_plus[k1][lp]));
                                    println!("  P_{{{},{}}}^+ = {}", k2, lp, pf(&p_plus[k2][lp]));
                                }
                            }
                        }
                    }
                }
                // Full upper group interlacing
                if upper.len() >= 2 {
                    let non_zero: Vec<usize> = upper
                        .iter()
                        .copied()
                        .filter(|&k| !pz(&p_plus[k][lp]))
                        .collect();
                    if non_zero.len() >= 2 {
                        upper_group_full[0] += 1;
                        let mut ok = true;
                        for i in 0..non_zero.len() - 1 {
                            if !interlaces(&p_plus[non_zero[i]][lp], &p_plus[non_zero[i + 1]][lp]) {
                                ok = false;
                                break;
                            }
                        }
                        if !ok {
                            upper_group_full[1] += 1;
                        }
                    }
                }

                // (b) LOWER GROUP: k' < l'
                // Sequence: P_{l'-1,l'}^+, P_{l'-2,l'}^+, ..., P_{1,l'}^+
                let lower: Vec<usize> = (1..lp).rev().collect();
                for i in 0..lower.len() {
                    if i + 1 < lower.len() {
                        let k1 = lower[i]; // larger k' (still < l')
                        let k2 = lower[i + 1]; // smaller k'
                        if !pz(&p_plus[k1][lp]) && !pz(&p_plus[k2][lp]) {
                            lower_group_pair[0] += 1;
                            if !interlaces(&p_plus[k1][lp], &p_plus[k2][lp]) {
                                lower_group_pair[1] += 1;
                                if lower_group_pair[1] <= 3 {
                                    println!("FAIL lower pair: board={:?} l'={} P_{{{},{}}}^+ << P_{{{},{}}}^+",
                                             board, lp, k1, lp, k2, lp);
                                    println!("  P_{{{},{}}}^+ = {}", k1, lp, pf(&p_plus[k1][lp]));
                                    println!("  P_{{{},{}}}^+ = {}", k2, lp, pf(&p_plus[k2][lp]));
                                }
                            }
                        }
                    }
                }
                // Full lower group interlacing
                if lower.len() >= 2 {
                    let non_zero: Vec<usize> = lower
                        .iter()
                        .copied()
                        .filter(|&k| !pz(&p_plus[k][lp]))
                        .collect();
                    if non_zero.len() >= 2 {
                        lower_group_full[0] += 1;
                        let mut ok = true;
                        for i in 0..non_zero.len() - 1 {
                            if !interlaces(&p_plus[non_zero[i]][lp], &p_plus[non_zero[i + 1]][lp]) {
                                ok = false;
                                break;
                            }
                        }
                        if !ok {
                            lower_group_full[1] += 1;
                        }
                    }
                }

                // (c) BRIDGE: P_{l'+1,l'}^+ << P_{l'-1,l'}^+
                if lp >= 2 && lp + 1 <= mp {
                    let kp_upper = lp + 1; // bottom of upper group
                    let kp_lower = lp - 1; // top of lower group
                    if !pz(&p_plus[kp_upper][lp]) && !pz(&p_plus[kp_lower][lp]) {
                        bridge_cond[0] += 1;
                        if !interlaces(&p_plus[kp_upper][lp], &p_plus[kp_lower][lp]) {
                            bridge_cond[1] += 1;
                            if bridge_cond[1] <= 5 {
                                println!(
                                    "FAIL BRIDGE: board={:?} l'={} P_{{{},{}}}^+ << P_{{{},{}}}^+",
                                    board, lp, kp_upper, lp, kp_lower, lp
                                );
                                println!(
                                    "  P_{{{},{}}}^+ = {}",
                                    kp_upper,
                                    lp,
                                    pf(&p_plus[kp_upper][lp])
                                );
                                println!(
                                    "  P_{{{},{}}}^+ = {}",
                                    kp_lower,
                                    lp,
                                    pf(&p_plus[kp_lower][lp])
                                );
                                println!(
                                    "  S_{{{},{}}}^+ = {}",
                                    kp_upper,
                                    lp,
                                    pf(&s_part[kp_upper][lp])
                                );
                                println!(
                                    "  T_{{{},{}}}^+ = {}",
                                    kp_upper,
                                    lp,
                                    pf(&t_part[kp_upper][lp])
                                );
                                println!(
                                    "  S_{{{},{}}}^+ = {}",
                                    kp_lower,
                                    lp,
                                    pf(&s_part[kp_lower][lp])
                                );
                                println!(
                                    "  T_{{{},{}}}^+ = {}",
                                    kp_lower,
                                    lp,
                                    pf(&t_part[kp_lower][lp])
                                );
                            }
                        }
                    }
                }

                // ────────────────────────────────────────────────
                // PART 4: Bridge decomposition (DU-type cross-column checks)
                // ────────────────────────────────────────────────
                if lp >= 2 && lp + 1 <= mp {
                    let kp_upper = lp + 1;
                    let kp_lower = lp - 1;
                    let s_u = &s_part[kp_upper][lp];
                    let t_u = &t_part[kp_upper][lp];
                    let s_l = &s_part[kp_lower][lp];
                    let t_l = &t_part[kp_lower][lp];

                    // S_{l'+1}^{(l')} << T_{l'-1}^{(l'-1)}
                    if !pz(s_u) && !pz(t_l) {
                        bridge_s_upper_t_lower[0] += 1;
                        if !interlaces(s_u, t_l) {
                            bridge_s_upper_t_lower[1] += 1;
                            if bridge_s_upper_t_lower[1] <= 3 {
                                println!("FAIL S_upper << T_lower: board={:?} l'={}", board, lp);
                                println!("  S_upper = {}", pf(s_u));
                                println!("  T_lower = {}", pf(t_l));
                            }
                        }
                    }
                    // T_{l'+1}^{(l')} << T_{l'-1}^{(l'-1)}
                    if !pz(t_u) && !pz(t_l) {
                        bridge_t_upper_t_lower[0] += 1;
                        if !interlaces(t_u, t_l) {
                            bridge_t_upper_t_lower[1] += 1;
                            if bridge_t_upper_t_lower[1] <= 3 {
                                println!("FAIL T_upper << T_lower: board={:?} l'={}", board, lp);
                                println!("  T_upper = {}", pf(t_u));
                                println!("  T_lower = {}", pf(t_l));
                            }
                        }
                    }
                    // S_{l'+1}^{(l')} << S_{l'-1}^{(l'-1)}
                    if !pz(s_u) && !pz(s_l) {
                        bridge_s_upper_s_lower[0] += 1;
                        if !interlaces(s_u, s_l) {
                            bridge_s_upper_s_lower[1] += 1;
                            if bridge_s_upper_s_lower[1] <= 3 {
                                println!("FAIL S_upper << S_lower: board={:?} l'={}", board, lp);
                                println!("  S_upper = {}", pf(s_u));
                                println!("  S_lower = {}", pf(s_l));
                            }
                        }
                    }
                    // T_{l'+1}^{(l')} << S_{l'-1}^{(l'-1)}
                    if !pz(t_u) && !pz(s_l) {
                        bridge_t_upper_s_lower[0] += 1;
                        if !interlaces(t_u, s_l) {
                            bridge_t_upper_s_lower[1] += 1;
                            if bridge_t_upper_s_lower[1] <= 3 {
                                println!("FAIL T_upper << S_lower: board={:?} l'={}", board, lp);
                                println!("  T_upper = {}", pf(t_u));
                                println!("  S_lower = {}", pf(s_l));
                            }
                        }
                    }
                }

                // ────────────────────────────────────────────────
                // PART 5: Full column interlacing
                // ────────────────────────────────────────────────
                // Full sequence: P_{mp,l'}^+, ..., P_{l'+1,l'}^+, P_{l'-1,l'}^+, ..., P_{1,l'}^+
                let mut full_seq: Vec<usize> = Vec::new();
                for kp in (lp + 1..=mp).rev() {
                    full_seq.push(kp);
                }
                for kp in (1..lp).rev() {
                    full_seq.push(kp);
                }

                let non_zero_full: Vec<usize> = full_seq
                    .iter()
                    .copied()
                    .filter(|&k| !pz(&p_plus[k][lp]))
                    .collect();

                // Check consecutive pairs in full sequence
                for i in 0..non_zero_full.len() {
                    if i + 1 < non_zero_full.len() {
                        let k1 = non_zero_full[i];
                        let k2 = non_zero_full[i + 1];
                        full_col_pair[0] += 1;
                        if !interlaces(&p_plus[k1][lp], &p_plus[k2][lp]) {
                            full_col_pair[1] += 1;
                            if full_col_pair[1] <= 5 {
                                println!("FAIL full col pair: board={:?} l'={} P_{{{},{}}}^+ << P_{{{},{}}}^+",
                                         board, lp, k1, lp, k2, lp);
                                println!("  P_{{{},{}}}^+ = {}", k1, lp, pf(&p_plus[k1][lp]));
                                println!("  P_{{{},{}}}^+ = {}", k2, lp, pf(&p_plus[k2][lp]));
                            }
                        }
                    }
                }

                // Full sequence interlacing check
                if non_zero_full.len() >= 2 {
                    full_col_seq[0] += 1;
                    let mut ok = true;
                    for i in 0..non_zero_full.len() - 1 {
                        if !interlaces(
                            &p_plus[non_zero_full[i]][lp],
                            &p_plus[non_zero_full[i + 1]][lp],
                        ) {
                            ok = false;
                            break;
                        }
                    }
                    if !ok {
                        full_col_seq[1] += 1;
                    }
                }
            }

            // ────────────────────────────────────────────────
            // PART 6: Cross-column conditions at lambda level
            // ────────────────────────────────────────────────
            for j in 1..=m {
                for l in 2..=m {
                    // A_{j,l} << A_{j,l-1}
                    if !pz(&a_jl[j][l]) && !pz(&a_jl[j][l - 1]) {
                        cross_a_a[0] += 1;
                        if !interlaces(&a_jl[j][l], &a_jl[j][l - 1]) {
                            cross_a_a[1] += 1;
                            if cross_a_a[1] <= 3 {
                                println!(
                                    "FAIL cross A_{{{},{}}} << A_{{{},{}}}: board={:?}",
                                    j,
                                    l,
                                    j,
                                    l - 1,
                                    board
                                );
                                println!("  A_{{{},{}}} = {}", j, l, pf(&a_jl[j][l]));
                                println!("  A_{{{},{}}} = {}", j, l - 1, pf(&a_jl[j][l - 1]));
                            }
                        }
                    }
                    // A_{j,l-1} << A_{j,l} (reverse direction)
                    if !pz(&a_jl[j][l - 1]) && !pz(&a_jl[j][l]) {
                        cross_a_a_rev[0] += 1;
                        if !interlaces(&a_jl[j][l - 1], &a_jl[j][l]) {
                            cross_a_a_rev[1] += 1;
                        }
                    }
                    // W_{j,l} << W_{j,l-1}
                    if !pz(&w_jl[j][l]) && !pz(&w_jl[j][l - 1]) {
                        cross_w_w[0] += 1;
                        if !interlaces(&w_jl[j][l], &w_jl[j][l - 1]) {
                            cross_w_w[1] += 1;
                            if cross_w_w[1] <= 3 {
                                println!(
                                    "FAIL cross W_{{{},{}}} << W_{{{},{}}}: board={:?}",
                                    j,
                                    l,
                                    j,
                                    l - 1,
                                    board
                                );
                                println!("  W_{{{},{}}} = {}", j, l, pf(&w_jl[j][l]));
                                println!("  W_{{{},{}}} = {}", j, l - 1, pf(&w_jl[j][l - 1]));
                            }
                        }
                    }
                    // W_{j,l-1} << W_{j,l} (reverse)
                    if !pz(&w_jl[j][l - 1]) && !pz(&w_jl[j][l]) {
                        cross_w_w_rev[0] += 1;
                        if !interlaces(&w_jl[j][l - 1], &w_jl[j][l]) {
                            cross_w_w_rev[1] += 1;
                        }
                    }
                }
            }

            // A_{j,l} << W_{j',l-1} for all j, j', l
            for j in 1..=m {
                for jp in 1..=m {
                    for l in 2..=m {
                        if !pz(&a_jl[j][l]) && !pz(&w_jl[jp][l - 1]) {
                            cross_a_w[0] += 1;
                            if !interlaces(&a_jl[j][l], &w_jl[jp][l - 1]) {
                                cross_a_w[1] += 1;
                                if cross_a_w[1] <= 3 {
                                    println!(
                                        "FAIL cross A_{{{},{}}} << W_{{{},{}}}: board={:?}",
                                        j,
                                        l,
                                        jp,
                                        l - 1,
                                        board
                                    );
                                    println!("  A_{{{},{}}} = {}", j, l, pf(&a_jl[j][l]));
                                    println!("  W_{{{},{}}} = {}", jp, l - 1, pf(&w_jl[jp][l - 1]));
                                }
                            }
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
    // SAMPLE DATA: Show D_{j,l}, U_{j,l} for small boards
    // ────────────────────────────────────────────────

    println!("\n================================================================");
    println!("  SAMPLE DATA: D_{{j,l}}, U_{{j,l}}, A_{{j,l}}, W_{{j,l}}");
    println!("================================================================\n");

    for n in 1..=std::cmp::min(max_n, 5) {
        let boards = gen_boards(n);
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            let m = (board[0] as usize).min(n);
            if n < 3 {
                continue;
            } // skip trivial cases

            let ideal = bruhat_lower_ideal(&perm);
            let mut d_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            let mut u_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];

            for pi in &ideal {
                let j = pi[0] as usize;
                if j > m {
                    continue;
                }
                let l = *pi.last().unwrap() as usize;
                if l > m {
                    continue;
                }
                let pk = peak_count(pi);
                let st = start_type(pi);
                let poly = if st == 'd' {
                    &mut d_jl[j][l]
                } else {
                    &mut u_jl[j][l]
                };
                while poly.len() <= pk {
                    poly.push(0);
                }
                poly[pk] += 1;
            }

            let has_data = (1..=m).any(|j| (1..=m).any(|l| !pz(&d_jl[j][l]) || !pz(&u_jl[j][l])));
            if !has_data {
                continue;
            }

            println!("Board {:?} (n={}, m={}, perm={:?}):", board, n, m, perm);
            for j in 1..=m {
                for l in 1..=m {
                    let a = pa(&d_jl[j][l], &u_jl[j][l]);
                    let w = pa(&pmt(&d_jl[j][l]), &u_jl[j][l]);
                    if pz(&a) {
                        continue;
                    }
                    print!("  (j={},l={}): ", j, l);
                    if !pz(&d_jl[j][l]) {
                        print!("D={} ", pf(&d_jl[j][l]));
                    }
                    if !pz(&u_jl[j][l]) {
                        print!("U={} ", pf(&u_jl[j][l]));
                    }
                    print!("A={} ", pf(&a));
                    print!("W={}", pf(&w));
                    println!();
                }
            }
            println!();
        }
    }

    // ────────────────────────────────────────────────
    // SAMPLE DATA: P_{k',l'}^+ for small boards
    // ────────────────────────────────────────────────

    println!("================================================================");
    println!("  SAMPLE DATA: P_{{k',l'}}^+ and S/T decomposition");
    println!("================================================================\n");

    for n in 1..=std::cmp::min(max_n, 4) {
        let boards = gen_boards(n);
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            let m = (board[0] as usize).min(n);
            if n < 2 {
                continue;
            }

            let ideal = bruhat_lower_ideal(&perm);
            let mut d_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            let mut u_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];

            for pi in &ideal {
                let j = pi[0] as usize;
                if j > m {
                    continue;
                }
                let l = *pi.last().unwrap() as usize;
                if l > m {
                    continue;
                }
                let pk = peak_count(pi);
                let st = start_type(pi);
                let poly = if st == 'd' {
                    &mut d_jl[j][l]
                } else {
                    &mut u_jl[j][l]
                };
                while poly.len() <= pk {
                    poly.push(0);
                }
                poly[pk] += 1;
            }

            let mut a_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            let mut w_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            for j in 1..=m {
                for l in 1..=m {
                    a_jl[j][l] = pa(&d_jl[j][l], &u_jl[j][l]);
                    w_jl[j][l] = pa(&pmt(&d_jl[j][l]), &u_jl[j][l]);
                }
            }

            let mp = m + 1;
            let mut s_col: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; mp + 1]; m + 1];
            let mut t_col: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; mp + 1]; m + 1];

            for col in 1..=m {
                for k in 1..=m {
                    while s_col[col].len() <= k + 1 {
                        s_col[col].push(vec![0i64]);
                    }
                    s_col[col][k + 1] = pa(&s_col[col][k], &a_jl[k][col]);
                }
                while t_col[col].len() <= mp {
                    t_col[col].push(vec![0i64]);
                }
                for k in (1..=m).rev() {
                    t_col[col][k] = pa(&t_col[col][k + 1], &w_jl[k][col]);
                }
            }

            println!("Board {:?} (n={}, m={}, perm={:?}):", board, n, m, perm);
            println!("  m' = {} (lambda+ level)", mp);

            for lp in 1..=mp {
                let mut any = false;
                for kp in 1..=mp {
                    if kp == lp {
                        continue;
                    }
                    let (col, label) = if kp > lp {
                        (lp, "upper")
                    } else {
                        (lp - 1, "lower")
                    };
                    if col < 1 || col > m {
                        continue;
                    }

                    let s = if kp < s_col[col].len() {
                        s_col[col][kp].clone()
                    } else {
                        if m + 1 < s_col[col].len() {
                            s_col[col][m + 1].clone()
                        } else {
                            vec![0i64]
                        }
                    };
                    let t = if kp <= m && kp < t_col[col].len() {
                        t_col[col][kp].clone()
                    } else {
                        vec![0i64]
                    };
                    let p = pa(&s, &t);
                    if pz(&p) {
                        continue;
                    }

                    if !any {
                        println!("  l' = {}:", lp);
                        any = true;
                    }
                    println!(
                        "    P_{{{},{}}}^+ = {} [{}; col={}; S={}, T={}]",
                        kp,
                        lp,
                        pf(&p),
                        label,
                        col,
                        pf(&s),
                        pf(&t)
                    );
                }
                if any {
                    println!();
                }
            }
        }
    }

    // ────────────────────────────────────────────────
    // SUMMARY
    // ────────────────────────────────────────────────

    println!("\n================================================================");
    println!("  SUMMARY");
    println!("================================================================\n");

    println!("Total valid 312-avoiding boards tested: {}\n", valid_boards);

    println!("=== PART 1: Real-rootedness of refined polynomials ===\n");
    show("  D_{j,l} real-rooted", djl_rr);
    show("  U_{j,l} real-rooted", ujl_rr);
    show("  A_{j,l} real-rooted", ajl_rr);
    show("  W_{j,l} real-rooted", wjl_rr);
    println!();

    println!("=== PART 3: Group interlacing at lambda+ level ===\n");
    println!("--- (a) UPPER GROUP (k' > l', fixed l') ---");
    show(
        "  Consecutive pairs P_{{k',l'}}^+ << P_{{k'-1,l'}}^+",
        upper_group_pair,
    );
    show("  Full upper group interlacing sequence", upper_group_full);
    println!();
    println!("--- (b) LOWER GROUP (k' < l', fixed l') ---");
    show(
        "  Consecutive pairs P_{{k',l'}}^+ << P_{{k'-1,l'}}^+",
        lower_group_pair,
    );
    show("  Full lower group interlacing sequence", lower_group_full);
    println!();
    println!("--- (c) BRIDGE CONDITION ---");
    show("  P_{{l'+1,l'}}^+ << P_{{l'-1,l'}}^+", bridge_cond);
    println!();

    println!("=== PART 4: Bridge decomposition (cross-column S/T checks) ===\n");
    show(
        "  S_upper << T_lower (S_{{l'+1}}^{{(l')}} << T_{{l'-1}}^{{(l'-1)}})",
        bridge_s_upper_t_lower,
    );
    show(
        "  T_upper << T_lower (T_{{l'+1}}^{{(l')}} << T_{{l'-1}}^{{(l'-1)}})",
        bridge_t_upper_t_lower,
    );
    show(
        "  S_upper << S_lower (S_{{l'+1}}^{{(l')}} << S_{{l'-1}}^{{(l'-1)}})",
        bridge_s_upper_s_lower,
    );
    show(
        "  T_upper << S_lower (T_{{l'+1}}^{{(l')}} << S_{{l'-1}}^{{(l'-1)}})",
        bridge_t_upper_s_lower,
    );
    println!();

    println!("=== PART 5: Full column interlacing (upper + bridge + lower) ===\n");
    show("  Consecutive pairs in full sequence", full_col_pair);
    show("  Full column interlacing sequence", full_col_seq);
    println!();

    println!("=== PART 6: Cross-column conditions at lambda level ===\n");
    show("  A_{{j,l}} << A_{{j,l-1}} (same j)", cross_a_a);
    show("  A_{{j,l-1}} << A_{{j,l}} (reverse)", cross_a_a_rev);
    show("  W_{{j,l}} << W_{{j,l-1}} (same j)", cross_w_w);
    show("  W_{{j,l-1}} << W_{{j,l}} (reverse)", cross_w_w_rev);
    show("  A_{{j,l}} << W_{{j',l-1}} (all j,j')", cross_a_w);
    println!();

    println!("================================================================");
}
