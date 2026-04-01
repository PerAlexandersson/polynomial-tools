//! Deep exploration of the j > l gap in the AW proof for peak polynomials on Ferrers boards.
//!
//! Uses the fast Bezout-based interlacing from polynomial-tools.
//! Tests ALL 312-avoiding Ferrers boards with n <= max_n.
//!
//! Usage: cargo run --release --bin peak_jgl -- 7

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
fn ps(a: &[i64], b: &[i64]) -> Vec<i64> {
    // a - b
    let l = a.len().max(b.len());
    let mut r = vec![0i64; l];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] -= v;
    }
    pt(&r)
}
fn pmt(p: &[i64]) -> Vec<i64> {
    // multiply by t
    let mut r = vec![0i64; p.len() + 1];
    for (i, &v) in p.iter().enumerate() {
        r[i + 1] = v;
    }
    pt(&r)
}
fn pscale(p: &[i64], c: i64) -> Vec<i64> {
    pt(&p.iter().map(|&x| x * c).collect::<Vec<_>>())
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

/// Check f << g (weak interlacing via Bezout + GCD).
/// Returns true if interlacing holds, false otherwise.
/// Treats zero f as trivially interlacing anything.
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
            // Degree difference check: for same-degree, check_weak_interlacing
            // may return None. In that case we also need to handle it.
            // check_weak_interlacing returns None when reduced degree diff != 1.
            // For same-degree interlacing we need a different approach.
            // Let's check both orientations: if deg(f) == deg(g), try f << g
            // by checking if t*f << g (padding).
            // Actually: the Bezout check only works for deg diff = 1.
            // For same degree, we check: does f interlace t*g? No, that changes the relation.
            //
            // For same-degree f, g: f << g means roots of f and g alternate with
            // f_1 <= g_1 <= f_2 <= g_2 <= ...
            // This is equivalent to: f interlaces (f + epsilon * g) for small epsilon,
            // but let's use a direct approach: f << g iff f << g+epsilon*t^{d+1}
            // Actually the standard trick: f << g (same deg d) iff
            //   f << (g + c * t^{d+1}) for the right sign of c.
            // But simpler: same-deg interlacing f << g iff
            //   the Bezout matrix of (g, f) (with g higher-degree by prepending 0) works...
            //
            // The simplest correct approach for same-degree:
            // f << g (same deg d, paper convention) iff g + epsilon * t^{d+1} is
            // interlaced by f for small epsilon > 0.
            // Equivalently: f << (g with a tiny extra leading term).
            //
            // For exact check: f << g (same deg) iff both:
            //   1) f << g (as deg d vs d+1, i.e., check_weak_interlacing(f, [g, 0, ..., 0, c]))
            //   Nah, this is getting complicated. Let's just use the Sturm-based check.
            //
            // Actually re-reading check_weak_interlacing: it does GCD division first.
            // If after GCD division the degrees are equal, it returns None.
            // This means the original polynomials had same degree with a shared factor
            // that doesn't reduce the degree diff to 1.
            //
            // For our use case, same-degree interlacing is actually not needed in the
            // standard way. The paper convention f << g allows deg(f) = deg(g) or
            // deg(f) + 1 = deg(g). The Bezout method only handles deg diff = 1.
            //
            // For same degree: f << g iff the interleaving condition holds.
            // We can check: is g + eps * t^{d+1} interlaced by f?
            // In exact arithmetic: replace eps with a large enough integer.
            // f << g (same deg d) iff there exists M such that f strictly interlaces g + M*t^{d+1}.
            // This is equivalent to: Bezout(g + M*t^{d+1}, f) is positive definite for large M.
            // But we can just try M=1, M=10, etc. Actually for the correct sign:
            // if lc(g) > 0, then g + M*t^{d+1} has degree d+1 with positive leading coeff.
            // Then f (deg d) interlaces g+M*t^{d+1} (deg d+1).
            // For large M, the extra root goes to -infinity (if M>0, lc>0).
            // So f << g (same deg) iff f strictly interlaces g + t^{d+1} * sign(lc(g)) * BIG.
            // But actually we just need: the roots of f alternate with roots of g,
            // with f_1 <= g_1 <= f_2 <= ...
            // This is the same as: g interlaces f (in the deg-diff-1 sense) when we view
            // f as having one more root at -infinity.
            // Wait no. Same degree: f << g means f_1 <= g_1 <= f_2 <= g_2 <= ...
            // This is equivalent to: t*f << g? No. Let me think more carefully.
            //
            // Actually for same degree d: f << g iff the polynomial g - epsilon * f
            // has only negative real roots for small epsilon. Hmm no.
            //
            // Simplest correct approach: pad g with a leading term.
            // If f has degree d and g has degree d (same), then f << g means
            // the roots interlace with f-root first. This is equivalent to
            // checking that f interlaces g' where g' has the roots of g plus
            // one root at -infinity. In practice: f (deg d) interlaces
            // g_padded = epsilon * t^{d+1} + g (deg d+1), for epsilon > 0 with
            // the sign chosen so lc > 0.
            //
            // For exact check: just set epsilon = 1 if lc(g) > 0.
            // Actually that changes the root positions. The correct way:
            // f (deg d) interlaces g_padded (deg d+1) means:
            //   g_padded_root_1 <= f_1 <= g_padded_root_2 <= ... <= f_d <= g_padded_root_{d+1}
            // As epsilon -> 0+, one root of g_padded goes to -infinity and the rest approach roots of g.
            // So: f_1 <= g_1 <= f_2 <= g_2 <= ... <= f_d <= g_d
            // which is exactly same-degree interlacing f << g!
            //
            // So: f << g (same deg d) iff check_weak_interlacing(f, g_padded) = Some(true)
            // where g_padded has degree d+1 with a small positive leading coefficient.
            // But "small" in exact arithmetic? We can use any positive value.
            // The question is whether the Bezout matrix changes character.
            // For epsilon > 0 small enough, B(g+eps*t^{d+1}, f) is positive definite
            // iff the limit (eps->0) condition holds. But for ANY eps > 0, it should work
            // as long as g + eps*t^{d+1} is still real-rooted.
            //
            // Safest: use eps = 1 and check. If g + t^{d+1} is not real-rooted, the
            // Bezout matrix won't be positive definite anyway, and the result is correct.
            // Actually no: g + t^{d+1} might not be real-rooted even if g is.
            //
            // Better: use the largest coefficient of g as scale.
            // Actually let's just use a large multiplier: g_padded = N*g + t^{d+1} * sign(lc(g))
            // For large N, roots of g_padded ≈ roots of N*g plus one near ∓infinity.

            let df = pdeg(&f);
            let dg = pdeg(&g);
            match (df, dg) {
                (Some(df), Some(dg)) if df == dg => {
                    // Same degree: pad g to make it degree d+1
                    let d = df;
                    let lc_g = g[d]; // leading coeff of g (after pt)
                    let sign = if lc_g > 0 { 1i64 } else { -1i64 };
                    // g_padded = BIG * g + sign * t^{d+1}
                    // Use BIG = max(|coeff|) * (d+1) + 1 to ensure root separation
                    let big: i64 =
                        g.iter().map(|&c| c.abs()).max().unwrap_or(1).max(1) * ((d as i64) + 1) + 1;
                    let mut g_padded: Vec<i64> = g.iter().map(|&c| c * big).collect();
                    while g_padded.len() <= d + 1 {
                        g_padded.push(0);
                    }
                    g_padded[d + 1] += sign;
                    let f_scaled: Vec<i64> = f.iter().map(|&c| c * big).collect();
                    // Now deg(g_padded) = d+1, deg(f_scaled) = d
                    match check_weak_interlacing(&f_scaled, &g_padded) {
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

// ── Board data structure ──

struct BoardData {
    board: Vec<u8>,
    n: usize,
    m: usize,
    d: Vec<Vec<i64>>,      // D_k polynomials, index 1..=m
    u: Vec<Vec<i64>>,      // U_k polynomials
    a: Vec<Vec<i64>>,      // A_k = D_k + U_k
    w: Vec<Vec<i64>>,      // W_k = t*D_k + U_k
    s: Vec<Vec<i64>>,      // S_k = sum_{j<k} A_j (partial sum, index 1..=m+1)
    t: Vec<Vec<i64>>,      // T_k = sum_{j>=k} W_j (tail sum, index 1..=m+1)
    a_plus: Vec<Vec<i64>>, // A_k^+ = S_k + T_k
    w_plus: Vec<Vec<i64>>, // W_k^+ = t*S_k + T_k
}

impl BoardData {
    fn compute(board: &[u8]) -> Self {
        let n = board.len();
        let perm = board_to_perm(board);
        let m = (board[0] as usize).min(n);
        let ideal = bruhat_lower_ideal(&perm);

        let mut d = vec![vec![0i64]; m + 1];
        let mut u = vec![vec![0i64]; m + 1];

        for pi in &ideal {
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
            let pk = peaks(pi);
            let poly = if pi[0] > pi[1] { &mut d[k] } else { &mut u[k] };
            while poly.len() <= pk {
                poly.push(0);
            }
            poly[pk] += 1;
        }

        let mut a = vec![vec![0i64]; m + 1];
        let mut w = vec![vec![0i64]; m + 1];
        for k in 1..=m {
            a[k] = pa(&d[k], &u[k]);
            w[k] = pa(&pmt(&d[k]), &u[k]);
        }

        // S_k = sum_{j < k} A_j  (so S_1 = 0, S_2 = A_1, S_3 = A_1 + A_2, ...)
        let mut s = vec![vec![0i64]; m + 2];
        for k in 2..=m + 1 {
            s[k] = pa(&s[k - 1], &a[k - 1]);
        }

        // T_k = sum_{j >= k} W_j  (so T_{m+1} = 0, T_m = W_m, T_{m-1} = W_{m-1} + W_m, ...)
        let mut t = vec![vec![0i64]; m + 2];
        for k in (1..=m).rev() {
            t[k] = pa(&t[k + 1], &w[k]);
        }

        // A_k^+ = S_k + T_k
        let mut a_plus = vec![vec![0i64]; m + 2];
        for k in 1..=m {
            a_plus[k] = pa(&s[k], &t[k]);
        }

        // W_k^+ = t*S_k + T_k
        let mut w_plus = vec![vec![0i64]; m + 2];
        for k in 1..=m {
            w_plus[k] = pa(&pmt(&s[k]), &t[k]);
        }

        BoardData {
            board: board.to_vec(),
            n,
            m,
            d,
            u,
            a,
            w,
            s,
            t,
            a_plus,
            w_plus,
        }
    }
}

fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);
    println!("================================================================");
    println!("  j > l GAP EXPLORATION for peak polynomials on Ferrers boards");
    println!("  312-avoiding boards, n <= {}", max_n);
    println!("  Using Bezout-based interlacing (polynomial-tools)");
    println!("================================================================\n");

    // Counters for Section 1: Degree analysis
    let mut total_pairs = 0u64;
    let mut deg_sl_lt_tl = 0u64; // deg(S_l) < deg(T_l)
    let mut deg_sl_ge_tl = 0u64; // deg(S_l) >= deg(T_l)

    // Counters for Section 2: Transitivity
    let mut chain1_aa = [0u64; 2]; // A_j^+ << A_l^+
    let mut chain1_aw = [0u64; 2]; // A_l^+ << W_l^+
    let mut chain1_direct = [0u64; 2]; // A_j^+ << W_l^+ (direct)
    let mut chain1_trans_eq = [0u64; 2]; // transitivity when deg(A_l^+) == deg(W_l^+)
    let mut chain1_trans_p1 = [0u64; 2]; // when deg(A_l^+) + 1 == deg(W_l^+)

    let mut chain2_at = [0u64; 2]; // A_j^+ << T_l
    let mut chain2_tw = [0u64; 2]; // T_l << W_l^+  (should hold from paper)
    let mut chain2_trans_eq = [0u64; 2]; // transitivity when deg(T_l) == deg(W_l^+)
    let mut chain2_trans_p1 = [0u64; 2]; // when deg(T_l) + 1 == deg(W_l^+)

    // Counters for Section 3: Cone decompositions
    let mut sj_wl = [0u64; 2]; // S_j << W_l^+
    let mut tj_wl = [0u64; 2]; // T_j << W_l^+
    let mut sj_tsl = [0u64; 2]; // S_j << tS_l (equivalently S_l << S_j by Wagner)
    let mut sj_tl = [0u64; 2]; // S_j << T_l
    let mut tj_tsl = [0u64; 2]; // T_j << tS_l (equiv S_l << T_j)
    let mut tj_tl = [0u64; 2]; // T_j << T_l (reverse of UU T_l << T_j)

    // Counters for Section 4: (t-1) structure
    let mut h_rr = [0u64; 2]; // h real-rooted?
    let mut h_direct = [0u64; 2]; // A_j^+ << A_j^+ + (t-1)*h

    // Counters for Section 5: Alternative common interlacers
    let mut ci_al_plus = [0u64; 2]; // A_l^+ as common interlacer (A_j^+ << A_l^+ << W_l^+)
    let mut ci_tl = [0u64; 2]; // T_l as common interlacer (A_j^+ << T_l << W_l^+)
    let mut ci_sl_tj = [0u64; 2]; // S_l + T_j as common interlacer
    let mut ci_total_p = [0u64; 2]; // P = sum_k A_k as common interlacer

    // Track hard cases
    let mut hard_cases: Vec<String> = Vec::new();
    let mut no_interlacer_cases: Vec<String> = Vec::new();

    // Detailed output for first few failures
    let max_detail = 5usize;
    let mut detail_count = 0usize;

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

            let bd = BoardData::compute(board);
            let m = bd.m;

            // Total peak polynomial P = sum_k A_k
            let mut total_peak = vec![0i64];
            for k in 1..=m {
                total_peak = pa(&total_peak, &bd.a[k]);
            }

            for j in 1..=m {
                for l in 1..j {
                    // j > l: this is the gap case

                    let aj_plus = &bd.a_plus[j];
                    let al_plus = &bd.a_plus[l];
                    let wl_plus = &bd.w_plus[l];
                    let sl = &bd.s[l];
                    let sj = &bd.s[j];
                    let tl = &bd.t[l];
                    let tj = &bd.t[j];

                    // Skip trivial cases where key polynomials are zero
                    if pz(aj_plus) || pz(wl_plus) {
                        continue;
                    }

                    total_pairs += 1;

                    let deg_aj_plus = pdeg(aj_plus);
                    let deg_al_plus = pdeg(al_plus);
                    let deg_wl_plus = pdeg(wl_plus);
                    let deg_sl = pdeg(sl);
                    let deg_tl = pdeg(tl);
                    let deg_sj = pdeg(sj);
                    let deg_tj = pdeg(tj);

                    // ── Section 1: Degree analysis ──
                    let sl_lt_tl = match (deg_sl, deg_tl) {
                        (Some(ds), Some(dt)) => ds < dt,
                        (None, Some(_)) => true,
                        _ => false,
                    };
                    if sl_lt_tl {
                        deg_sl_lt_tl += 1;
                    } else {
                        deg_sl_ge_tl += 1;
                    }

                    // ── Section 2: Transitivity tests ──

                    // Chain 1: A_j^+ << A_l^+ << W_l^+
                    let aa_ok = if !pz(al_plus) {
                        chain1_aa[0] += 1;
                        let r = interlaces(aj_plus, al_plus);
                        if !r {
                            chain1_aa[1] += 1;
                        }
                        r
                    } else {
                        false
                    };

                    let aw_diag = if !pz(al_plus) {
                        chain1_aw[0] += 1;
                        let r = interlaces(al_plus, wl_plus);
                        if !r {
                            chain1_aw[1] += 1;
                        }
                        r
                    } else {
                        false
                    };

                    // Direct check A_j^+ << W_l^+
                    chain1_direct[0] += 1;
                    let direct_ok = interlaces(aj_plus, wl_plus);
                    if !direct_ok {
                        chain1_direct[1] += 1;
                    }

                    // Transitivity analysis for chain 1
                    if aa_ok && aw_diag {
                        match (deg_al_plus, deg_wl_plus) {
                            (Some(da), Some(dw)) if da == dw => {
                                chain1_trans_eq[0] += 1;
                                if !direct_ok {
                                    chain1_trans_eq[1] += 1;
                                }
                            }
                            (Some(da), Some(dw)) if da + 1 == dw => {
                                chain1_trans_p1[0] += 1;
                                if !direct_ok {
                                    chain1_trans_p1[1] += 1;
                                }
                            }
                            _ => {}
                        }
                    }

                    // Chain 2: A_j^+ << T_l << W_l^+
                    let at_ok = if !pz(tl) {
                        chain2_at[0] += 1;
                        let r = interlaces(aj_plus, tl);
                        if !r {
                            chain2_at[1] += 1;
                        }
                        r
                    } else {
                        false
                    };

                    let tw_ok = if !pz(tl) {
                        chain2_tw[0] += 1;
                        let r = interlaces(tl, wl_plus);
                        if !r {
                            chain2_tw[1] += 1;
                        }
                        r
                    } else {
                        false
                    };

                    if at_ok && tw_ok {
                        match (deg_tl, deg_wl_plus) {
                            (Some(dt), Some(dw)) if dt == dw => {
                                chain2_trans_eq[0] += 1;
                                if !direct_ok {
                                    chain2_trans_eq[1] += 1;
                                }
                            }
                            (Some(dt), Some(dw)) if dt + 1 == dw => {
                                chain2_trans_p1[0] += 1;
                                if !direct_ok {
                                    chain2_trans_p1[1] += 1;
                                }
                            }
                            _ => {}
                        }
                    }

                    // ── Section 3: Cone decompositions ──
                    if !pz(sj) {
                        sj_wl[0] += 1;
                        if !interlaces(sj, wl_plus) {
                            sj_wl[1] += 1;
                        }
                    }
                    if !pz(tj) {
                        tj_wl[0] += 1;
                        if !interlaces(tj, wl_plus) {
                            tj_wl[1] += 1;
                        }
                    }
                    // S_j << tS_l means S_j << t*S_l. By Wagner reversal, this
                    // is equivalent to S_l << S_j (if S_l has smaller degree).
                    let tsl = pmt(sl);
                    if !pz(sj) && !pz(&tsl) {
                        sj_tsl[0] += 1;
                        if !interlaces(sj, &tsl) {
                            sj_tsl[1] += 1;
                        }
                    }
                    if !pz(sj) && !pz(tl) {
                        sj_tl[0] += 1;
                        if !interlaces(sj, tl) {
                            sj_tl[1] += 1;
                        }
                    }
                    if !pz(tj) && !pz(&tsl) {
                        tj_tsl[0] += 1;
                        if !interlaces(tj, &tsl) {
                            tj_tsl[1] += 1;
                        }
                    }
                    if !pz(tj) && !pz(tl) {
                        tj_tl[0] += 1;
                        if !interlaces(tj, tl) {
                            tj_tl[1] += 1;
                        }
                    }

                    // ── Section 4: (t-1) structure ──
                    // W_l^+ = A_j^+ + (t-1)*h where h = sum_{i<j} D_i + sum_{i<l} U_i
                    // Actually: W_l^+ - A_j^+ = tS_l + T_l - S_j - T_j
                    //   = tS_l - S_j + T_l - T_j
                    //   = tS_l - S_j + (sum_{l<=i<j} W_i)  [since T_l - T_j = sum_{l<=i<j} W_i for l < j]
                    //   Wait, let me recompute.
                    // Actually: W_l^+ - A_j^+ = (tS_l + T_l) - (S_j + T_j)
                    // T_l - T_j = sum_{l <= i < j} W_i  (for l < j, since T_k = sum_{i>=k} W_i)
                    // S_j - S_l = sum_{l <= i < j} A_i  (partial sum difference)
                    // tS_l - S_l = (t-1)*S_l
                    // So: W_l^+ - A_j^+ = (t-1)*S_l + T_l - T_j - (S_j - S_l)
                    //    Hmm, let me just do it: tS_l + T_l - S_j - T_j
                    //    = tS_l - S_l + S_l + T_l - S_j - T_j
                    //    = (t-1)S_l + (S_l + T_l) - (S_j + T_j)
                    //    = (t-1)S_l + A_l^+ - A_j^+
                    //
                    // Alternative direct: W_l^+ - A_j^+ computed directly
                    let diff = ps(wl_plus, aj_plus); // W_l^+ - A_j^+
                                                     // Check: is diff = (t-1)*h for some h?
                                                     // (t-1)*h means coefficients: -h[0], h[0]-h[1], h[1]-h[2], ...
                                                     // Factor out (t-1): diff / (t-1)
                                                     // If diff = c_0 + c_1*t + c_2*t^2 + ...
                                                     // then h = -c_0 + (-c_0-c_1)*t + ... (cumulative sum with negation... or)
                                                     // Actually (t-1) * h(t) = t*h(t) - h(t)
                                                     // If h = h_0 + h_1*t + ..., then (t-1)*h = -h_0 + (h_0-h_1)*t + (h_1-h_2)*t^2 + ...
                                                     // So h_0 = -diff[0], h_1 = h_0 - diff[1], h_2 = h_1 - diff[2], ...
                                                     // h_k = -sum_{i=0}^{k} diff[i]
                    let diff_v = pt(&diff);
                    let mut h = vec![0i64; diff_v.len()];
                    let mut cumsum = 0i64;
                    for i in 0..diff_v.len() {
                        cumsum += diff_v[i];
                        h[i] = -cumsum;
                    }
                    // Verify: the last entry of h should give 0 when we compute (t-1)*h
                    // Actually: (t-1)*h has degree = deg(h) + 1, but diff might have degree = deg(h).
                    // The cumsum of diff should be 0 for the factorization to work.
                    let h = pt(&h);
                    // Verify factorization
                    let tm1_h = ps(&pmt(&h), &h); // t*h - h = (t-1)*h
                    let tm1_h = pt(&tm1_h);
                    let diff_v = pt(&diff);
                    let factorization_ok = tm1_h == diff_v;

                    if factorization_ok && !pz(&h) {
                        h_rr[0] += 1;
                        if !is_real_rooted(&h) {
                            h_rr[1] += 1;
                        }

                        h_direct[0] += 1;
                        // A_j^+ << A_j^+ + (t-1)*h should be same as A_j^+ << W_l^+
                        let ajp_plus_tm1h = pa(aj_plus, &tm1_h);
                        let r = interlaces(aj_plus, &ajp_plus_tm1h);
                        if !r {
                            h_direct[1] += 1;
                        }
                    }

                    // ── Section 5: Common interlacers ──

                    // Candidate 1: A_l^+
                    if !pz(al_plus) {
                        let leg1 = interlaces(aj_plus, al_plus);
                        let leg2 = interlaces(al_plus, wl_plus);
                        ci_al_plus[0] += 1;
                        if !(leg1 && leg2) {
                            ci_al_plus[1] += 1;
                        }
                    }

                    // Candidate 2: T_l
                    if !pz(tl) {
                        let leg1 = interlaces(aj_plus, tl);
                        let leg2 = interlaces(tl, wl_plus);
                        ci_tl[0] += 1;
                        if !(leg1 && leg2) {
                            ci_tl[1] += 1;
                        }
                    }

                    // Candidate 3: S_l + T_j
                    let sl_tj = pa(sl, tj);
                    if !pz(&sl_tj) {
                        let leg1 = interlaces(aj_plus, &sl_tj);
                        let leg2 = interlaces(&sl_tj, wl_plus);
                        ci_sl_tj[0] += 1;
                        if !(leg1 && leg2) {
                            ci_sl_tj[1] += 1;
                        }
                    }

                    // Candidate 4: P = total peak polynomial
                    if !pz(&total_peak) {
                        let leg1 = interlaces(aj_plus, &total_peak);
                        let leg2 = interlaces(&total_peak, wl_plus);
                        ci_total_p[0] += 1;
                        if !(leg1 && leg2) {
                            ci_total_p[1] += 1;
                        }
                    }

                    // Track hard cases and report
                    if !sl_lt_tl {
                        let case_desc = format!(
                            "n={} board={:?} j={} l={} deg(S_l)={:?} deg(T_l)={:?}",
                            n, board, j, l, deg_sl, deg_tl
                        );

                        // Check if any common interlacer works
                        let al_works = !pz(al_plus)
                            && interlaces(aj_plus, al_plus)
                            && interlaces(al_plus, wl_plus);
                        let tl_works =
                            !pz(tl) && interlaces(aj_plus, tl) && interlaces(tl, wl_plus);
                        let sl_tj_works = !pz(&sl_tj)
                            && interlaces(aj_plus, &sl_tj)
                            && interlaces(&sl_tj, wl_plus);
                        let p_works = !pz(&total_peak)
                            && interlaces(aj_plus, &total_peak)
                            && interlaces(&total_peak, wl_plus);

                        let any_works = al_works || tl_works || sl_tj_works || p_works || direct_ok;

                        if !any_works {
                            no_interlacer_cases.push(case_desc.clone());
                        }

                        if hard_cases.len() < 50 {
                            hard_cases.push(format!(
                                "{} direct={} A_l+={} T_l={} S_l+T_j={} P={}",
                                case_desc, direct_ok, al_works, tl_works, sl_tj_works, p_works
                            ));
                        }
                    }

                    // Print details for first few direct failures
                    if !direct_ok && detail_count < max_detail {
                        detail_count += 1;
                        println!(
                            "--- DETAIL #{}: direct A_j^+ << W_l^+ FAILS ---",
                            detail_count
                        );
                        println!("  board={:?} n={} m={} j={} l={}", board, n, m, j, l);
                        println!("  A_j^+ = {}", pf(aj_plus));
                        println!("  W_l^+ = {}", pf(wl_plus));
                        println!(
                            "  deg(A_j^+)={:?} deg(W_l^+)={:?}",
                            deg_aj_plus, deg_wl_plus
                        );
                        println!("  A_l^+ = {}", pf(al_plus));
                        println!("  T_l   = {}", pf(tl));
                        println!("  S_l   = {}", pf(sl));
                        println!("  S_j   = {}", pf(sj));
                        println!("  T_j   = {}", pf(tj));
                        println!(
                            "  deg(S_l)={:?} deg(T_l)={:?} sl<tl={}",
                            deg_sl, deg_tl, sl_lt_tl
                        );

                        // Check all candidates
                        let aa = interlaces(aj_plus, al_plus);
                        let aw = interlaces(al_plus, wl_plus);
                        let at = interlaces(aj_plus, tl);
                        let tw = interlaces(tl, wl_plus);
                        println!("  Chain1: A_j^+<<A_l^+={} A_l^+<<W_l^+={}", aa, aw);
                        println!("  Chain2: A_j^+<<T_l={} T_l<<W_l^+={}", at, tw);
                        println!();
                    }
                }
            }
        }
        println!(
            "n={}: {} valid 312-avoiding boards (cumulative: {})",
            n, n_valid, valid_boards
        );
    }

    // ── SUMMARY OUTPUT ──

    println!("\n================================================================");
    println!("  SUMMARY RESULTS");
    println!("================================================================\n");

    let show = |name: &str, c: [u64; 2]| {
        if c[0] == 0 {
            println!("  {}: (no data)", name);
        } else if c[1] == 0 {
            println!("  {}: {}/{} ALL PASS", name, c[0], c[0]);
        } else {
            println!("  {}: {}/{} pass ({} FAIL)", name, c[0] - c[1], c[0], c[1]);
        }
    };

    println!("Total (j,l) pairs with j > l tested: {}\n", total_pairs);

    println!("=== SECTION 1: Degree analysis ===");
    println!(
        "  deg(S_l) < deg(T_l):  {} pairs (transitivity-friendly)",
        deg_sl_lt_tl
    );
    println!(
        "  deg(S_l) >= deg(T_l): {} pairs (hard cases)",
        deg_sl_ge_tl
    );
    println!();

    println!("=== SECTION 2: Transitivity tests ===");
    println!("  Chain 1: A_j^+ << A_l^+ << W_l^+ (reversed AA + diagonal AW)");
    show("    A_j^+ << A_l^+  (reversed AA)", chain1_aa);
    show("    A_l^+ << W_l^+  (diagonal AW)", chain1_aw);
    show("    A_j^+ << W_l^+  (direct check)", chain1_direct);
    show("    Transitivity deg(A_l^+)=deg(W_l^+)", chain1_trans_eq);
    show("    Transitivity deg(A_l^+)+1=deg(W_l^+)", chain1_trans_p1);
    println!();
    println!("  Chain 2: A_j^+ << T_l << W_l^+ (from the paper)");
    show("    A_j^+ << T_l", chain2_at);
    show("    T_l << W_l^+", chain2_tw);
    show("    Transitivity deg(T_l)=deg(W_l^+)", chain2_trans_eq);
    show("    Transitivity deg(T_l)+1=deg(W_l^+)", chain2_trans_p1);
    println!();

    println!("=== SECTION 3: Alternative cone decompositions ===");
    show("  S_j << W_l^+", sj_wl);
    show("  T_j << W_l^+", tj_wl);
    show("  S_j << tS_l", sj_tsl);
    show("  S_j << T_l", sj_tl);
    show("  T_j << tS_l (S_l << T_j)", tj_tsl);
    show("  T_j << T_l (reverse of T_l << T_j)", tj_tl);
    println!();

    println!("=== SECTION 4: (t-1) structure ===");
    println!("  W_l^+ - A_j^+ = (t-1)*h analysis:");
    show("  h real-rooted", h_rr);
    show("  A_j^+ << A_j^+ + (t-1)*h", h_direct);
    println!();

    println!("=== SECTION 5: Alternative common interlacers ===");
    println!("  For each (j,l) pair, checking candidate C: A_j^+ << C << W_l^+");
    show("  C = A_l^+", ci_al_plus);
    show("  C = T_l", ci_tl);
    show("  C = S_l + T_j", ci_sl_tj);
    show("  C = P (total peak poly)", ci_total_p);
    println!();

    println!("=== HARD CASES: deg(S_l) >= deg(T_l) ===");
    if hard_cases.is_empty() {
        println!("  None!");
    } else {
        for (i, c) in hard_cases.iter().enumerate() {
            println!("  [{}] {}", i + 1, c);
        }
    }
    println!();

    println!("=== CRITICAL: pairs where NO tested common interlacer works ===");
    if no_interlacer_cases.is_empty() {
        println!("  NONE - every hard case has at least one working strategy!");
    } else {
        println!("  {} cases found:", no_interlacer_cases.len());
        for c in &no_interlacer_cases {
            println!("    {}", c);
        }
    }
    println!();

    println!("================================================================");
    println!("  FINAL TALLY");
    println!("================================================================");
    println!("  Total boards: {}", valid_boards);
    println!("  Total (j>l) pairs: {}", total_pairs);
    println!(
        "  Direct A_j^+ << W_l^+ pass: {}/{}",
        chain1_direct[0] - chain1_direct[1],
        chain1_direct[0]
    );
    println!("  deg(S_l) < deg(T_l) (easy): {}", deg_sl_lt_tl);
    println!("  deg(S_l) >= deg(T_l) (hard): {}", deg_sl_ge_tl);
    println!(
        "  Hard cases with no interlacer: {}",
        no_interlacer_cases.len()
    );
    println!("================================================================");
}
