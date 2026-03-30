//! Exploration of partial-sum UU condition, weighted partial sums, recursion
//! structure, and Brändén matrix for 312-avoiding Ferrers boards.
//!
//! Part 1: Partial sum UU -- Σ_{i<j} U_{i,l} ≪ Σ_{i≥j} U_{i,l}
//! Part 2: Weighted partial sums -- Σ_{i<j} U_i ≪ Σ_{i≥j} (U_i + α·D_i)
//! Part 3: Connection to recursion structure (weighted sums of W^μ)
//! Part 4: Brändén matrix structure
//!
//! Usage: cargo run --release --bin peak_partial_uu -- 7

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
/// Scalar multiply
fn ps(c: i64, p: &[i64]) -> Vec<i64> {
    pt(&p.iter().map(|&v| c * v).collect::<Vec<_>>())
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

/// For a board λ = (λ_1, ..., λ_n), extract the "previous board" μ s.t.
/// λ = (m, μ_1+1, ..., μ_{n-1}+1) where m = λ_1.
/// Returns (m, μ) or None if n=1 (base case).
fn extract_prev_board(board: &[u8]) -> Option<(usize, Vec<u8>)> {
    if board.len() <= 1 { return None; }
    let m = board[0] as usize;
    // μ_i = λ_{i+1} - 1 for i = 1, ..., n-1
    let mu: Vec<u8> = board[1..].iter().map(|&v| v - 1).collect();
    // Check μ is a valid board (all entries >= 1, non-decreasing starting from position index+1)
    if mu.iter().any(|&v| v < 1) { return None; }
    Some((m, mu))
}

// ══════════════════════════════════════════════════════════════════════
//  MAIN
// ══════════════════════════════════════════════════════════════════════

fn main() {
    let max_n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(7);
    println!("================================================================");
    println!("  Partial-sum UU exploration");
    println!("  312-avoiding Ferrers boards, n <= {}", max_n);
    println!("================================================================\n");

    // Part 1 counters: partial sum UU
    // For each cutoff j: Σ_{i<j} U_{i,l} ≪ Σ_{i≥j} U_{i,l}
    let mut puu = [0u64; 2];

    // Part 2 counters: weighted partial sums
    // α = 0: Σ_{i<j} U_i ≪ Σ_{i≥j} U_i  (same as Part 1)
    // α = 1: Σ_{i<j} U_i ≪ Σ_{i≥j} (U_i + D_i) = Σ_{i≥j} A_i
    // α = t: Σ_{i<j} U_i ≪ Σ_{i≥j} (U_i + t·D_i) = Σ_{i≥j} W_i
    let mut wu_a0 = [0u64; 2]; // α = 0 (redundant with Part 1, but tracked separately)
    let mut wu_a1 = [0u64; 2]; // α = 1
    let mut wu_at = [0u64; 2]; // α = t

    // Also test: Σ_{i<j} D_i ≪ Σ_{i≥j} U_i (partial sum DU)
    let mut pdu = [0u64; 2];
    // And: Σ_{i<j} A_i ≪ Σ_{i≥j} W_i (partial sum AW)
    let mut paw = [0u64; 2];
    // And: Σ_{i<j} W_i ≪ Σ_{i≥j} W_i (partial sum WW)
    let mut pww = [0u64; 2];

    // Part 3: recursion structure verification
    let mut rec_verified = 0u64;
    let mut rec_checked = 0u64;
    let mut rec_mismatch = 0u64;

    // Part 4: Brändén matrix construction & verification
    let mut mat_verified = 0u64;
    let mut mat_checked = 0u64;
    let mut mat_mismatch = 0u64;

    let mut valid_boards = 0u64;

    // Collect all boards indexed by n for Part 3 lookups
    let mut all_boards_by_n: Vec<Vec<Vec<u8>>> = vec![vec![]; max_n + 1];
    for n in 1..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let perm = board_to_perm(board);
            if is_312_avoiding(&perm) {
                all_boards_by_n[n].push(board.clone());
            }
        }
    }

    // Cache: (board) -> (d_jl, u_jl, w_jl, a_jl) for Part 3 lookups
    // Key = board as Vec<u8>
    use std::collections::HashMap;
    struct BoardData {
        m: usize,
        d_jl: Vec<Vec<Vec<i64>>>,
        u_jl: Vec<Vec<Vec<i64>>>,
        w_jl: Vec<Vec<Vec<i64>>>,
        a_jl: Vec<Vec<Vec<i64>>>,
    }
    let mut board_cache: HashMap<Vec<u8>, BoardData> = HashMap::new();

    for n in 1..=max_n {
        for board in &all_boards_by_n[n] {
            let perm = board_to_perm(board);
            valid_boards += 1;

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

            // Store in cache
            board_cache.insert(board.clone(), BoardData {
                m,
                d_jl: d_jl.clone(),
                u_jl: u_jl.clone(),
                w_jl: w_jl.clone(),
                a_jl: a_jl.clone(),
            });

            // ══════════════════════════════════════════════════════
            //  PART 1 & 2: Partial sum conditions
            // ══════════════════════════════════════════════════════

            for l in 1..=m {
                // Build prefix sums of U, D, A, W for column l
                // prefix_u[j] = Σ_{i=1}^{j-1} U_{i,l}  (i.e., Σ_{i<j} U_i)
                // suffix_u[j] = Σ_{i=j}^{m} U_{i,l}     (i.e., Σ_{i≥j} U_i)
                let mut prefix_u = vec![vec![0i64]; m + 2];
                let mut prefix_d = vec![vec![0i64]; m + 2];
                let mut prefix_a = vec![vec![0i64]; m + 2];
                let mut prefix_w = vec![vec![0i64]; m + 2];

                for i in 1..=m {
                    prefix_u[i + 1] = pa(&prefix_u[i], &u_jl[i][l]);
                    prefix_d[i + 1] = pa(&prefix_d[i], &d_jl[i][l]);
                    prefix_a[i + 1] = pa(&prefix_a[i], &a_jl[i][l]);
                    prefix_w[i + 1] = pa(&prefix_w[i], &w_jl[i][l]);
                }

                let total_u = prefix_u[m + 1].clone();
                let total_d = prefix_d[m + 1].clone();
                let total_a = prefix_a[m + 1].clone();
                let total_w = prefix_w[m + 1].clone();

                // suffix[j] = total - prefix[j]
                let suffix_u = |j: usize| -> Vec<i64> {
                    pa(&total_u, &ps(-1, &prefix_u[j]))
                };
                let suffix_a = |j: usize| -> Vec<i64> {
                    pa(&total_a, &ps(-1, &prefix_a[j]))
                };
                let suffix_w = |j: usize| -> Vec<i64> {
                    pa(&total_w, &ps(-1, &prefix_w[j]))
                };

                for j in 1..=m {
                    let pu = &prefix_u[j]; // Σ_{i<j} U_i
                    let su = suffix_u(j);  // Σ_{i≥j} U_i

                    // Part 1: Σ_{i<j} U_i ≪ Σ_{i≥j} U_i
                    if !pz(pu) && !pz(&su) {
                        puu[0] += 1;
                        if !interlaces(pu, &su) { puu[1] += 1; }
                    }

                    // Part 2 α=0: same test (redundant but tracked)
                    if !pz(pu) && !pz(&su) {
                        wu_a0[0] += 1;
                        if !interlaces(pu, &su) { wu_a0[1] += 1; }
                    }

                    // Part 2 α=1: Σ_{i<j} U_i ≪ Σ_{i≥j} A_i
                    let sa = suffix_a(j);
                    if !pz(pu) && !pz(&sa) {
                        wu_a1[0] += 1;
                        if !interlaces(pu, &sa) { wu_a1[1] += 1; }
                    }

                    // Part 2 α=t: Σ_{i<j} U_i ≪ Σ_{i≥j} W_i
                    let sw = suffix_w(j);
                    if !pz(pu) && !pz(&sw) {
                        wu_at[0] += 1;
                        if !interlaces(pu, &sw) { wu_at[1] += 1; }
                    }

                    // Bonus: Σ_{i<j} D_i ≪ Σ_{i≥j} U_i
                    let pd = &prefix_d[j];
                    if !pz(pd) && !pz(&su) {
                        pdu[0] += 1;
                        if !interlaces(pd, &su) { pdu[1] += 1; }
                    }

                    // Bonus: Σ_{i<j} A_i ≪ Σ_{i≥j} W_i
                    let pa_ = &prefix_a[j];
                    if !pz(pa_) && !pz(&sw) {
                        paw[0] += 1;
                        if !interlaces(pa_, &sw) { paw[1] += 1; }
                    }

                    // Bonus: Σ_{i<j} W_i ≪ Σ_{i≥j} W_i
                    let pw = &prefix_w[j];
                    if !pz(pw) && !pz(&sw) {
                        pww[0] += 1;
                        if !interlaces(pw, &sw) { pww[1] += 1; }
                    }
                }
            }
        }

        println!("n={}: {} boards (cumulative: {})", n,
                 all_boards_by_n[n].len(), valid_boards);
    }

    // ══════════════════════════════════════════════════════
    //  PART 3: Recursion structure (n=3,4,5)
    // ══════════════════════════════════════════════════════

    println!("\n================================================================");
    println!("  PART 3: Recursion structure verification");
    println!("  For boards lam with n=3,4,5, verify U_{{k,l}}^lam = Sum_{{j>=k}} W_{{j,l'}}^mu");
    println!("  and Sum_{{i<cutoff}} U_i^lam = Sum_j min(j,cutoff-1) * W_j^mu");
    println!("================================================================\n");

    let part3_max = max_n.min(5);
    for n in 2..=part3_max {
        for board in &all_boards_by_n[n] {
            let m = board[0] as usize;
            if let Some((m_val, mu)) = extract_prev_board(board) {
                assert_eq!(m_val, m);
                let mu_n = mu.len();
                if mu_n == 0 { continue; }

                // Check if mu is a valid 312-avoiding Ferrers board
                let mu_perm = board_to_perm(&mu);
                if !is_312_avoiding(&mu_perm) { continue; }

                // Get data for mu from cache
                let mu_data = match board_cache.get(&mu) {
                    Some(d) => d,
                    None => continue,
                };

                // Get data for board from cache
                let board_data = match board_cache.get(board) {
                    Some(d) => d,
                    None => continue,
                };

                let m_mu = mu_data.m;

                // For each column l of λ (from 1 to m):
                // When k > l: the "input column" for D^+/U^+ is l' = l
                // When k < l: the "input column" for D^+/U^+ is l' = l - 1
                //
                // D_{k,l}^λ = S_k^{(l')} = Σ_{j<k} A_{j,l'}^μ
                // U_{k,l}^λ = T_k^{(l')} = Σ_{j≥k} W_{j,l'}^μ

                let print_detail = n <= 4; // Print details for small boards

                if print_detail {
                    println!("  Board λ = {:?} (n={}), prev μ = {:?} (n={}), m={}, m_μ={}",
                             board, n, mu, mu_n, m, m_mu);
                }

                for l in 1..=m {
                    for k in 1..=m {
                        if k == l { continue; }
                        let lp = if k > l { l } else { l - 1 };
                        if lp < 1 || lp > m_mu { continue; }

                        // Verify U_{k,l}^λ = T_k^{(lp)} = Σ_{j≥k} W_{j,lp}^μ
                        // But note: indices in μ go from 1 to m_mu
                        let mut tail_w = vec![0i64];
                        for j in k..=m_mu {
                            tail_w = pa(&tail_w, &mu_data.w_jl[j][lp]);
                        }

                        rec_checked += 1;
                        let u_kl = &board_data.u_jl[k][l];
                        if pt(u_kl) != pt(&tail_w) {
                            rec_mismatch += 1;
                            if print_detail {
                                println!("    MISMATCH U_{{{},{}}}^λ: got {} vs expected {} (tail W from μ, lp={})",
                                         k, l, format_poly(u_kl), format_poly(&tail_w), lp);
                            }
                        } else {
                            rec_verified += 1;
                            if print_detail && !pz(u_kl) {
                                println!("    OK U_{{{},{}}}^λ = Σ_{{j≥{}}} W_{{j,{}}}^μ = {}",
                                         k, l, k, lp, format_poly(u_kl));
                            }
                        }

                        // Verify D_{k,l}^λ = S_k^{(lp)} = Σ_{j<k} A_{j,lp}^μ
                        let mut prefix_a_mu = vec![0i64];
                        for j in 1..k {
                            if j <= m_mu {
                                prefix_a_mu = pa(&prefix_a_mu, &mu_data.a_jl[j][lp]);
                            }
                        }

                        rec_checked += 1;
                        let d_kl = &board_data.d_jl[k][l];
                        if pt(d_kl) != pt(&prefix_a_mu) {
                            rec_mismatch += 1;
                            if print_detail {
                                println!("    MISMATCH D_{{{},{}}}^λ: got {} vs expected {} (prefix A from μ, lp={})",
                                         k, l, format_poly(d_kl), format_poly(&prefix_a_mu), lp);
                            }
                        } else {
                            rec_verified += 1;
                            if print_detail && !pz(d_kl) {
                                println!("    OK D_{{{},{}}}^λ = Σ_{{j<{}}} A_{{j,{}}}^μ = {}",
                                         k, l, k, lp, format_poly(d_kl));
                            }
                        }
                    }
                }

                // Now verify: Σ_{i<cutoff} U_{i,l}^λ = Σ_j min(j, cutoff-1) · W_{j,lp}^μ
                // Here we sum over i in the "k > l" group, so lp = l.
                // Actually we need to be careful: U_{i,l} for different i may come from
                // different groups. Let's just verify the identity for a fixed group.
                //
                // For the group where k > l (so lp = l):
                //   Σ_{i=l+1}^{m} U_{i,l}^λ = Σ_{i=l+1}^{m} Σ_{j≥i} W_{j,l}^μ
                //     = Σ_j W_{j,l}^μ · #{i : l+1 ≤ i ≤ j}  (only i from l+1..m, j from l+1..m_mu)
                //     = Σ_j W_{j,l}^μ · (j - l) if j ≤ m, else · (m - l)
                // Wait, more precisely: we need j ≥ i and l+1 ≤ i, so:
                //     = Σ_{j=l+1}^{m_mu} W_{j,l}^μ · min(j - l, m - l)
                //
                // For partial sum Σ_{i < cutoff, i > l} U_{i,l}^λ = Σ_{i=l+1}^{cutoff-1} Σ_{j≥i} W_{j,l}^μ
                //     = Σ_j W_{j,l}^μ · #{i : l+1 ≤ i ≤ min(j, cutoff-1)}
                //     = Σ_j W_{j,l}^μ · max(0, min(j, cutoff-1) - l)

                // Weighted sum check: for each output column l of λ,
                // pick the group k > l with input column lp = l.
                for l_out in 1..m {
                    let lp = l_out;
                    if lp > m_mu { continue; }
                    if !print_detail { continue; }

                    println!("    --- Weighted sum check for group k>l={} (lp={}) ---", l_out, lp);
                    for cutoff in (l_out + 1)..=m {
                        // Direct: Σ_{i=l_out+1}^{cutoff-1} U_{i,l_out}^λ
                        let mut direct_sum = vec![0i64];
                        for i in (l_out + 1)..cutoff {
                            if i <= m {
                                direct_sum = pa(&direct_sum, &board_data.u_jl[i][l_out]);
                            }
                        }

                        // Via weighted W: Σ_j weight_j · W_{j,lp}^μ
                        // weight_j = max(0, min(j, cutoff-1) - l_out)
                        let mut weighted_sum = vec![0i64];
                        for j in 1..=m_mu {
                            let weight = (j.min(cutoff - 1) as i64 - l_out as i64).max(0);
                            if weight > 0 {
                                weighted_sum = pa(&weighted_sum, &ps(weight, &mu_data.w_jl[j][lp]));
                            }
                        }

                        if pt(&direct_sum) != pt(&weighted_sum) {
                            println!("      cutoff={}: MISMATCH direct={} vs weighted={}",
                                     cutoff, format_poly(&direct_sum), format_poly(&weighted_sum));
                        } else if !pz(&direct_sum) {
                            println!("      cutoff={}: OK Sum_{{i<{}}} U_{{i,{}}} = weighted sum = {}",
                                     cutoff, cutoff, l_out, format_poly(&direct_sum));
                        }
                    }
                }

                if print_detail {
                    println!();
                }
            }
        }
    }

    println!("  Recursion structure: checked={}, verified={}, mismatch={}",
             rec_checked, rec_verified, rec_mismatch);

    // ══════════════════════════════════════════════════════
    //  PART 4: Brändén matrix
    // ══════════════════════════════════════════════════════

    println!("\n================================================================");
    println!("  PART 4: Brändén matrix construction");
    println!("  For n=3,4: construct M s.t. M·(D_1,...,D_m,U_1,...,U_m)^T = (D_1^+,...,U_m'^+)^T");
    println!("================================================================\n");

    let part4_max = max_n.min(4);
    for n in 2..=part4_max {
        for board in &all_boards_by_n[n] {
            let m = board[0] as usize;
            if let Some((_m_val, mu)) = extract_prev_board(board) {
                let mu_n = mu.len();
                if mu_n == 0 { continue; }
                let mu_perm = board_to_perm(&mu);
                if !is_312_avoiding(&mu_perm) { continue; }

                let mu_data = match board_cache.get(&mu) {
                    Some(d) => d,
                    None => continue,
                };
                let board_data = match board_cache.get(board) {
                    Some(d) => d,
                    None => continue,
                };

                let m_mu = mu_data.m;

                println!("  Board λ = {:?} (n={}, m={}), prev μ = {:?} (m_μ={})", board, n, m, mu, m_mu);

                // For a fixed input column lp (of μ), the recursion gives:
                // Output D_{k}^+ and U_{k}^+ for k in the group using column lp.
                //
                // If lp is used by k > l (where l is the output last-entry column),
                // then l = lp, and k ranges over {lp+1, ..., m}.
                //
                // D_k^+ = Σ_{j<k} D_j^μ + Σ_{j<k} U_j^μ  = Σ_{j<k} A_j^μ  (prefix sum of A)
                // U_k^+ = t·Σ_{j≥k} D_j^μ + Σ_{j≥k} U_j^μ = Σ_{j≥k} W_j^μ (suffix sum of W)
                //
                // So the matrix M maps the 2m_μ-vector (D_1,...,D_{m_μ}, U_1,...,U_{m_μ}) in column lp
                // to a vector of (D_{k}^+, U_{k}^+) entries.

                for lp in 1..=m_mu {
                    // Determine which output k values use this column
                    // Group k > l=lp: k = lp+1, ..., m  (output column l = lp)
                    // Group k < l=lp+1: k = 1, ..., lp  (output column l = lp+1, and lp = l-1)
                    // Actually, let's just build the matrix for the "k > l" group where l = lp.
                    // Output k values: lp+1, ..., m

                    let out_ks_above: Vec<usize> = ((lp + 1)..=m).collect();
                    // And for the "k < l" group with l = lp+1 (so input col = l-1 = lp):
                    // Only valid if lp+1 <= m (there is a valid output column)
                    let out_ks_below: Vec<usize> = if lp + 1 <= m {
                        (1..=lp).collect()
                    } else {
                        vec![]
                    };
                    // Both groups use the same input column lp.

                    let all_out_ks: Vec<usize> = out_ks_below.iter().chain(out_ks_above.iter()).cloned().collect();

                    if all_out_ks.is_empty() { continue; }

                    println!("    Input column lp={}, m_μ={}", lp, m_mu);
                    println!("    Output k values (k<l={}): {:?}", lp + 1, out_ks_below);
                    println!("    Output k values (k>l={}): {:?}", lp, out_ks_above);

                    // Build the 2*|output| x 2*m_mu matrix M
                    // Row ordering: D_{k1}^+, D_{k2}^+, ..., U_{k1}^+, U_{k2}^+, ...
                    // Column ordering: D_1^μ, D_2^μ, ..., D_{m_μ}^μ, U_1^μ, U_2^μ, ..., U_{m_μ}^μ
                    //
                    // For D_k^+ = Σ_{j<k} A_j = Σ_{j<k} D_j + Σ_{j<k} U_j:
                    //   coeff of D_j = 1 if j < k, else 0
                    //   coeff of U_j = 1 if j < k, else 0
                    //
                    // For U_k^+ = Σ_{j≥k} W_j = t·Σ_{j≥k} D_j + Σ_{j≥k} U_j:
                    //   coeff of D_j = t if j >= k, else 0
                    //   coeff of U_j = 1 if j >= k, else 0

                    let num_out = all_out_ks.len();
                    let num_in = 2 * m_mu;

                    // Matrix entries: "0", "1", or "t"
                    let mut matrix: Vec<Vec<&str>> = vec![vec!["0"; num_in]; 2 * num_out];

                    for (idx, &k) in all_out_ks.iter().enumerate() {
                        // D_k^+ row (index idx)
                        for j in 1..=m_mu {
                            if j < k {
                                matrix[idx][j - 1] = "1";        // coeff of D_j
                                matrix[idx][m_mu + j - 1] = "1"; // coeff of U_j
                            }
                        }
                        // U_k^+ row (index num_out + idx)
                        for j in 1..=m_mu {
                            if j >= k {
                                matrix[num_out + idx][j - 1] = "t";        // coeff of D_j
                                matrix[num_out + idx][m_mu + j - 1] = "1"; // coeff of U_j
                            }
                        }
                    }

                    // Print the matrix
                    println!("    M ({} x {}):", 2 * num_out, num_in);
                    let d_labels: Vec<String> = (1..=m_mu).map(|j| format!("D{}", j)).collect();
                    let u_labels: Vec<String> = (1..=m_mu).map(|j| format!("U{}", j)).collect();
                    let col_labels: Vec<String> = d_labels.iter().chain(u_labels.iter()).cloned().collect();
                    println!("      Cols: {:?}", col_labels);
                    for (idx, &k) in all_out_ks.iter().enumerate() {
                        println!("      D{}^+: {:?}", k, matrix[idx]);
                    }
                    for (idx, &k) in all_out_ks.iter().enumerate() {
                        println!("      U{}^+: {:?}", k, matrix[num_out + idx]);
                    }

                    // Verify M · input = output
                    for (idx, &k) in all_out_ks.iter().enumerate() {
                        // Determine output column l for this k
                        let l = if out_ks_below.contains(&k) { lp + 1 } else { lp };

                        // Compute D_k^+ from matrix
                        let mut d_from_mat = vec![0i64];
                        for j in 1..=m_mu {
                            if matrix[idx][j - 1] == "1" {
                                d_from_mat = pa(&d_from_mat, &mu_data.d_jl[j][lp]);
                            }
                            if matrix[idx][m_mu + j - 1] == "1" {
                                d_from_mat = pa(&d_from_mat, &mu_data.u_jl[j][lp]);
                            }
                        }

                        mat_checked += 1;
                        if pt(&d_from_mat) != pt(&board_data.d_jl[k][l]) {
                            mat_mismatch += 1;
                            println!("      MISMATCH D_{{{}}}^+ (l={}): mat gives {} vs actual {}",
                                     k, l, format_poly(&d_from_mat), format_poly(&board_data.d_jl[k][l]));
                        } else {
                            mat_verified += 1;
                        }

                        // Compute U_k^+ from matrix
                        let mut u_from_mat = vec![0i64];
                        for j in 1..=m_mu {
                            if matrix[num_out + idx][j - 1] == "t" {
                                u_from_mat = pa(&u_from_mat, &pmt(&mu_data.d_jl[j][lp]));
                            }
                            if matrix[num_out + idx][m_mu + j - 1] == "1" {
                                u_from_mat = pa(&u_from_mat, &mu_data.u_jl[j][lp]);
                            }
                        }

                        mat_checked += 1;
                        if pt(&u_from_mat) != pt(&board_data.u_jl[k][l]) {
                            mat_mismatch += 1;
                            println!("      MISMATCH U_{{{}}}^+ (l={}): mat gives {} vs actual {}",
                                     k, l, format_poly(&u_from_mat), format_poly(&board_data.u_jl[k][l]));
                        } else {
                            mat_verified += 1;
                        }
                    }

                    println!();
                }
            }
        }
    }

    println!("  Brändén matrix: checked={}, verified={}, mismatch={}",
             mat_checked, mat_verified, mat_mismatch);

    // ══════════════════════════════════════════════════════
    //  RESULTS SUMMARY
    // ══════════════════════════════════════════════════════

    println!("\n================================================================");
    println!("  PART 1: Partial sum UU");
    println!("  Σ_{{i<j}} U_{{i,l}} ≪ Σ_{{i≥j}} U_{{i,l}}");
    println!("================================================================");
    show("Partial UU: Σ U_{<j} ≪ Σ U_{≥j}", puu);

    println!("\n================================================================");
    println!("  PART 2: Weighted partial sums");
    println!("================================================================");
    show("α=0: Σ U_{<j} ≪ Σ U_{≥j}          ", wu_a0);
    show("α=1: Σ U_{<j} ≪ Σ A_{≥j}          ", wu_a1);
    show("α=t: Σ U_{<j} ≪ Σ W_{≥j}          ", wu_at);
    show("     Σ D_{<j} ≪ Σ U_{≥j}  (partial DU)", pdu);
    show("     Σ A_{<j} ≪ Σ W_{≥j}  (partial AW)", paw);
    show("     Σ W_{<j} ≪ Σ W_{≥j}  (partial WW)", pww);

    println!("\n================================================================");
    println!("  Total boards tested: {}", valid_boards);
    println!("================================================================");
}
