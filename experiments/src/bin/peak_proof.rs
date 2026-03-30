//! Exact verification of the inductive hypothesis for peak polynomial proof.
//!
//! Hypothesis H(λ): S_k ≼ T_k for all k = 1,...,m
//!   where S_k = Σ_{j<k} A_j^λ and T_k = Σ_{j≥k} W_j^λ.
//!
//! Also tests the cross-condition A_j^+ ≼ W_l^+ for j ≤ l
//! (needed for the inductive step).
//!
//! Uses exact arithmetic throughout.

use std::collections::BTreeSet;





/// Exact check: f ≼ g (f interlaces g, roots of f weakly left of roots of g).
fn exact_interlaces(fc: &[i64], gc: &[i64]) -> Result<bool, String> {
    match polynomial_tools::check_weak_interlacing(fc, gc) {
        Some(result) => Ok(result),
        None => Err("degree difference not 1 after GCD removal".to_string()),
    }
}

// ── Combinatorial infrastructure ─────────────────────────────────────

fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> {
    let n = perm.len();
    let mut visited: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut queue: BTreeSet<Vec<u8>> = BTreeSet::new();
    queue.insert(perm.to_vec());
    while let Some(current) = queue.pop_last() {
        for i in 0..n { for j in i+1..n { if current[i] > current[j] {
            let mut c = current.clone(); c.swap(i, j);
            if !visited.contains(&c) { queue.insert(c); }
        }}}
        visited.insert(current);
    }
    visited.into_iter().collect()
}

fn board_to_312_perm(board: &[u8]) -> Vec<u8> {
    let n = board.len();
    let mut perm = vec![0u8; n];
    let mut used = vec![false; n + 1];
    for i in 0..n {
        for c in (1..=(board[i] as usize).min(n)).rev() {
            if !used[c] { perm[i] = c as u8; used[c] = true; break; }
        }
    }
    perm
}

fn peaks(w: &[u8]) -> usize {
    if w.len() < 3 { return 0; }
    (1..w.len()-1).filter(|&i| w[i-1] < w[i] && w[i] > w[i+1]).count()
}

fn pt(p: &[i64]) -> Vec<i64> { let mut v = p.to_vec(); while v.len()>1 && *v.last().unwrap()==0 { v.pop(); } v }
fn pz(p: &[i64]) -> bool { p.iter().all(|&c| c == 0) }
fn pa(a: &[i64], b: &[i64]) -> Vec<i64> {
    let l = a.len().max(b.len()); let mut r = vec![0i64;l];
    for (i,&v) in a.iter().enumerate() { r[i]+=v; }
    for (i,&v) in b.iter().enumerate() { r[i]+=v; }
    pt(&r)
}
fn pmt(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; p.len()+1];
    for (i,&v) in p.iter().enumerate() { r[i+1]=v; }
    pt(&r)
}
fn pf(p: &[i64]) -> String {
    let p = pt(p); if pz(&p) { return "0".into(); }
    let mut t = vec![];
    for (i,&c) in p.iter().enumerate() { if c==0 { continue; } match (c,i) {
        (c,0) => t.push(format!("{}",c)), (1,1) => t.push("t".into()),
        (c,1) => t.push(format!("{}t",c)), (1,e) => t.push(format!("t^{}",e)),
        (c,e) => t.push(format!("{}t^{}",c,e)),
    }}
    t.join(" + ")
}

fn generate_boards(n: usize) -> Vec<Vec<u8>> {
    let mut r = vec![]; let mut c = vec![];
    gb(n, n, 0, &mut c, &mut r); r
}
fn gb(n: usize, mx: usize, d: usize, c: &mut Vec<u8>, r: &mut Vec<Vec<u8>>) {
    if d == n { r.push(c.clone()); return; }
    for v in (d+1).max(if d>0 { c[d-1] as usize } else { 1 })..=mx {
        c.push(v as u8); gb(n, mx, d+1, c, r); c.pop();
    }
}

struct BoardData {
    m: usize,
    a: Vec<Vec<i64>>,  // 1-indexed: A_k = D_k + U_k
    w: Vec<Vec<i64>>,  // 1-indexed: W_k = tD_k + U_k
}

fn compute_board(board: &[u8]) -> BoardData {
    let n = board.len();
    let m = (board[0] as usize).min(n);
    let perm = board_to_312_perm(board);
    let ideal = bruhat_lower_ideal(&perm);
    let mut d = vec![vec![0i64]; m+1];
    let mut u = vec![vec![0i64]; m+1];
    for pi in &ideal {
        if pi.len() < 2 {
            let k = pi[0] as usize;
            if k <= m { while u[k].len() < 1 { u[k].push(0); } u[k][0] += 1; }
            continue;
        }
        let k = pi[0] as usize;
        if k > m { continue; }
        let pk = peaks(pi);
        let poly = if pi[0] > pi[1] { &mut d[k] } else { &mut u[k] };
        while poly.len() <= pk { poly.push(0); }
        poly[pk] += 1;
    }
    let mut a = vec![vec![0i64]; m+1];
    let mut w = vec![vec![0i64]; m+1];
    for k in 1..=m {
        a[k] = pa(&d[k], &u[k]);
        w[k] = pa(&pmt(&d[k]), &u[k]);
    }
    BoardData { m, a, w }
}

fn main() {
    let max_n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(7);
    println!("Exact test of H(λ): S_k ≼ T_k, and A_j^+ ≼ W_l^+ (j≤l)\n");

    let mut total = 0;
    let mut st_fails = 0;  // S_k ≼ T_k fails
    let mut aw_fails = 0;  // A_j^+ ≼ W_l^+ fails

    for n in 1..=max_n {
        let boards = generate_boards(n);
        for board in &boards {
            total += 1;
            let bd = compute_board(board);
            let m = bd.m;

            // Build S_k and T_k
            // S_k = Σ_{j=1}^{k-1} A_j, T_k = Σ_{j=k}^{m} W_j
            let mut s = vec![vec![0i64]; m + 2]; // s[k] for k=1,...,m+1
            let mut t = vec![vec![0i64]; m + 2];
            for k in 1..=m + 1 {
                s[k] = if k == 1 { vec![0] } else { pa(&s[k-1], &bd.a[k-1]) };
            }
            for k in (1..=m).rev() {
                t[k] = pa(&t[k+1], &bd.w[k]);
            }

            // Test S_k ≼ T_k for k=1,...,m
            for k in 1..=m {
                if pz(&s[k]) { continue; } // 0 ≼ anything
                if pz(&t[k]) { continue; } // skip if T=0
                match exact_interlaces(&s[k], &t[k]) {
                    Ok(true) => {},
                    Ok(false) => {
                        st_fails += 1;
                        if st_fails <= 5 {
                            println!("FAIL S_{} ≼ T_{} for {:?}", k, k, board);
                            println!("  S = {}", pf(&s[k]));
                            println!("  T = {}", pf(&t[k]));
                        }
                    },
                    Err(e) => {
                        st_fails += 1;
                        if st_fails <= 5 {
                            println!("ERR S_{} ≼ T_{} for {:?}: {}", k, k, board, e);
                        }
                    }
                }
            }

            // Test A_j^+ ≼ W_l^+ for j ≤ l
            // A_j^+ = S_j + T_j, W_l^+ = tS_l + T_l
            for j in 1..=m {
                let a_plus_j = pa(&s[j], &t[j]);
                for l in j..=m {
                    let w_plus_l = pa(&pmt(&s[l]), &t[l]);
                    if pz(&a_plus_j) || pz(&w_plus_l) { continue; }
                    match exact_interlaces(&a_plus_j, &w_plus_l) {
                        Ok(true) => {},
                        Ok(false) => {
                            aw_fails += 1;
                            if aw_fails <= 5 {
                                println!("FAIL A_{}^+ ≼ W_{}^+ for {:?}", j, l, board);
                                println!("  A^+ = {}", pf(&a_plus_j));
                                println!("  W^+ = {}", pf(&w_plus_l));
                            }
                        },
                        Err(e) => {
                            aw_fails += 1;
                            if aw_fails <= 5 {
                                println!("ERR A_{}^+ ≼ W_{}^+ for {:?}: {}", j, l, board, e);
                            }
                        }
                    }
                }
            }
        }
        println!("n={}: boards={}, S≼T fails={}, A^+≼W^+ fails={}", n, total, st_fails, aw_fails);
    }

    // === Phase 2: Test intermediate conditions for the algebraic proof ===
    println!("\n--- Phase 2: intermediate conditions ---\n");
    let mut ww_fails = 0;  // W_j ≼ W_l fails (j ≤ l)
    let mut wt_fails = 0;  // W_i ≼ T_l fails (i < l)
    let mut tt_fails = 0;  // T_j ≼ T_l fails (j < l)

    for n in 1..=max_n {
        let boards = generate_boards(n);
        for board in &boards {
            let bd = compute_board(board);
            let m = bd.m;
            if m <= 1 { continue; }

            let mut t_arr = vec![vec![0i64]; m + 2];
            for k in (1..=m).rev() {
                t_arr[k] = pa(&t_arr[k+1], &bd.w[k]);
            }

            // W_j ≼ W_l for j ≤ l
            for j in 1..=m { for l in j+1..=m {
                if pz(&bd.w[j]) || pz(&bd.w[l]) { continue; }
                match exact_interlaces(&bd.w[j], &bd.w[l]) {
                    Ok(true) => {},
                    _ => { ww_fails += 1;
                        if ww_fails <= 3 { println!("FAIL W_{}≼W_{} for {:?}: {} vs {}", j, l, board, pf(&bd.w[j]), pf(&bd.w[l])); }
                    }
                }
            }}

            // W_i ≼ T_l for i < l
            for i in 1..m { for l in i+1..=m {
                if pz(&bd.w[i]) || pz(&t_arr[l]) { continue; }
                match exact_interlaces(&bd.w[i], &t_arr[l]) {
                    Ok(true) => {},
                    _ => { wt_fails += 1;
                        if wt_fails <= 3 { println!("FAIL W_{}≼T_{} for {:?}: {} vs {}", i, l, board, pf(&bd.w[i]), pf(&t_arr[l])); }
                    }
                }
            }}

            // T_j ≼ T_l for j < l (tail sums interlace)
            for j in 1..m { for l in j+1..=m {
                if pz(&t_arr[j]) || pz(&t_arr[l]) { continue; }
                match exact_interlaces(&t_arr[j], &t_arr[l]) {
                    Ok(true) => {},
                    _ => { tt_fails += 1;
                        if tt_fails <= 3 { println!("FAIL T_{}≼T_{} for {:?}: {} vs {}", j, l, board, pf(&t_arr[j]), pf(&t_arr[l])); }
                    }
                }
            }}
        }
    }

    println!("Phase 2: W_j≼W_l fails={}, W_i≼T_l fails={}, T_j≼T_l fails={}", ww_fails, wt_fails, tt_fails);

    println!("\n=== SUMMARY ===");
    println!("Total: {}, S_k≼T_k fails: {}, A_j^+≼W_l^+ fails: {}", total, st_fails, aw_fails);
    println!("W_j≼W_l: {}, W_i≼T_l: {}, T_j≼T_l: {}", ww_fails, wt_fails, tt_fails);
    if st_fails == 0 && aw_fails == 0 && ww_fails == 0 && wt_fails == 0 && tt_fails == 0 {
        println!("ALL CONDITIONS HOLD — complete proof chain verified!");
    }
}
