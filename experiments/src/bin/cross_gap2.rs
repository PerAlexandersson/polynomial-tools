//! Test: U_{1,l'} ≪ W_{i,l'-1} (cross-column U1≪W) and
//! D_{k'}^+(col l') ≪ W_{i,l'-1} (cross-column D+≪W).
//! If BOTH hold → cross-group DU by right-cone + transitivity.
use polynomial_tools::real_rootedness::check_weak_interlacing;
use std::collections::BTreeSet;
fn pt(p: &[i64]) -> Vec<i64> { let mut v = p.to_vec(); while v.len() > 1 && *v.last().unwrap() == 0 { v.pop(); } v }
fn pz(p: &[i64]) -> bool { p.iter().all(|&c| c == 0) }
fn pa(a: &[i64], b: &[i64]) -> Vec<i64> { let l = a.len().max(b.len()); let mut r = vec![0i64; l]; for (i, &v) in a.iter().enumerate() { r[i] += v; } for (i, &v) in b.iter().enumerate() { r[i] += v; } pt(&r) }
fn pmt(p: &[i64]) -> Vec<i64> { let mut r = vec![0i64; p.len() + 1]; for (i, &v) in p.iter().enumerate() { r[i + 1] = v; } pt(&r) }
fn pdeg(p: &[i64]) -> Option<usize> { let v = pt(p); if pz(&v) { None } else { Some(v.len() - 1) } }
fn interlaces(f: &[i64], g: &[i64]) -> bool {
    let f = pt(f); let g = pt(g);
    if pz(&f) { return true; }
    if pz(&g) { return false; }
    check_weak_interlacing(&f, &g).unwrap_or(false)
}
fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> { let n = perm.len(); let mut vis: BTreeSet<Vec<u8>> = BTreeSet::new(); let mut q: BTreeSet<Vec<u8>> = BTreeSet::new(); q.insert(perm.to_vec()); while let Some(cur) = q.pop_last() { for i in 0..n { for j in i+1..n { if cur[i] > cur[j] { let mut c = cur.clone(); c.swap(i, j); if !vis.contains(&c) { q.insert(c); } } } } vis.insert(cur); } vis.into_iter().collect() }
fn board_to_perm(b: &[u8]) -> Vec<u8> { let n = b.len(); let mut p = vec![0u8; n]; let mut u = vec![false; n+1]; for i in 0..n { for c in (1..=(b[i] as usize).min(n)).rev() { if !u[c] { p[i] = c as u8; u[c] = true; break; } } } p }
fn is_312_avoiding(perm: &[u8]) -> bool { let n = perm.len(); for i in 0..n { for j in i+1..n { for k in j+1..n { if perm[k] < perm[i] && perm[i] < perm[j] { return false; } } } } true }
fn peaks(w: &[u8]) -> usize { if w.len() < 3 { return 0; } (1..w.len()-1).filter(|&i| w[i-1] < w[i] && w[i] > w[i+1]).count() }
fn gen_boards(n: usize) -> Vec<Vec<u8>> { let mut r = vec![]; let mut c = vec![]; gb(n, n, 0, &mut c, &mut r); r }
fn gb(n: usize, mx: usize, d: usize, c: &mut Vec<u8>, r: &mut Vec<Vec<u8>>) { if d == n { r.push(c.clone()); return; } for v in (d+1).max(if d > 0 { c[d-1] as usize } else { 1 })..=mx { c.push(v as u8); gb(n, mx, d+1, c, r); c.pop(); } }
fn main() {
    let max_n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(7);
    // Strategy: D_k^+(col1) ≪ U_1(col1) ≪ W_i(col2) → D_k^+ ≪ W_i → D_k^+ ≪ U_j^+(col2)
    let mut u1_w_cross = [0u64; 2]; // U_{1,l'} ≪ W_{i,l'-1} for all i
    let mut dk_w_cross = [0u64; 2]; // D_k^+(col1) ≪ W_{i,col2} for all cross-group pairs
    let mut dk_u1_same = [0u64; 2]; // D_k^+(col1) ≪ U_{1,col1} (= A_{1,col1})
    // Also: A_{1,l'} ≪ A_{i,l'-1} for all i (cross-col A1 vs A)
    let mut a1_a_cross = [0u64; 2];
    // A_{1,l'} ≪ W_{i,l'-1} for all i
    let mut a1_w_cross = [0u64; 2];
    // Also test the OTHER cross direction: D_j^+(col2) ≪ U_1(col2) ≪ W_i(col1)?
    // For the other cross-group: k' < l', j' > l'. D uses col l'-1, U uses col l'.
    let mut u1_w_cross2 = [0u64; 2]; // U_{1,l'-1} ≪ W_{i,l'} for all i
    for n in 1..=max_n {
        for board in &gen_boards(n) {
            let perm = board_to_perm(board); if !is_312_avoiding(&perm) { continue; }
            let m = board[0] as usize; if n <= 2 { continue; }
            let ideal = bruhat_lower_ideal(&perm);
            let mut d_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m+1]; m+1];
            let mut u_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m+1]; m+1];
            for pi in &ideal { let j = pi[0] as usize; if j > m { continue; } let l = *pi.last().unwrap() as usize; if l > m { continue; } let pk = peaks(pi);
                let poly = if pi.len() >= 2 && pi[0] > pi[1] { &mut d_jl[j][l] } else { &mut u_jl[j][l] }; while poly.len() <= pk { poly.push(0); } poly[pk] += 1; }
            let mut a_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m+1]; m+1];
            let mut w_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m+1]; m+1];
            for j in 1..=m { for l in 1..=m { a_jl[j][l] = pa(&d_jl[j][l], &u_jl[j][l]); w_jl[j][l] = pa(&pmt(&d_jl[j][l]), &u_jl[j][l]); } }
            // Test cross-column conditions between adjacent columns
            for l in 2..=m {
                let col1 = l;     // "upper" column
                let col2 = l - 1; // "lower" column
                // U_{1,col1} ≪ W_{i,col2} for all i
                for i in 1..=m {
                    if !pz(&u_jl[1][col1]) && !pz(&w_jl[i][col2]) {
                        u1_w_cross[0] += 1;
                        if !interlaces(&u_jl[1][col1], &w_jl[i][col2]) { u1_w_cross[1] += 1; }
                    }
                    // A_{1,col1} ≪ A_{i,col2}
                    if !pz(&a_jl[1][col1]) && !pz(&a_jl[i][col2]) {
                        a1_a_cross[0] += 1;
                        if !interlaces(&a_jl[1][col1], &a_jl[i][col2]) { a1_a_cross[1] += 1; }
                    }
                    // A_{1,col1} ≪ W_{i,col2}
                    if !pz(&a_jl[1][col1]) && !pz(&w_jl[i][col2]) {
                        a1_w_cross[0] += 1;
                        if !interlaces(&a_jl[1][col1], &w_jl[i][col2]) { a1_w_cross[1] += 1; }
                    }
                    // U_{1,col2} ≪ W_{i,col1}
                    if !pz(&u_jl[1][col2]) && !pz(&w_jl[i][col1]) {
                        u1_w_cross2[0] += 1;
                        if !interlaces(&u_jl[1][col2], &w_jl[i][col1]) { u1_w_cross2[1] += 1; }
                    }
                }
                // D_k^+(col1) ≪ W_{i,col2} for cross-group pairs
                let mp = m + 1;
                for kp in (l+1)..=mp {
                    // D_k^+(col1) = Σ_{j<kp} A_{j,col1}
                    let mut dk = vec![0i64];
                    for j in 1..kp.min(m+1) { dk = pa(&dk, &a_jl[j][col1]); }
                    for i in 1..=m {
                        if !pz(&dk) && !pz(&w_jl[i][col2]) {
                            dk_w_cross[0] += 1;
                            if !interlaces(&dk, &w_jl[i][col2]) { dk_w_cross[1] += 1; }
                        }
                    }
                    // D_k^+(col1) ≪ U_{1,col1}
                    if !pz(&dk) && !pz(&u_jl[1][col1]) {
                        dk_u1_same[0] += 1;
                        if !interlaces(&dk, &u_jl[1][col1]) { dk_u1_same[1] += 1; }
                    }
                }
            }
        }
    }
    println!("=== Cross-group strategy: transitivity via U_1 ===");
    let show = |name: &str, c: [u64;2]| { if c[0]==0 { println!("  {}: (no data)", name); } else if c[1]==0 { println!("  {}: {}/{} ALL PASS <<<", name, c[0], c[0]); } else { println!("  {}: {}/{} pass ({} FAIL)", name, c[0]-c[1], c[0], c[1]); } };
    show("D_k^+(col1) ≪ U_{1,col1}", dk_u1_same);
    show("U_{1,l'} ≪ W_{i,l'-1}", u1_w_cross);
    show("A_{1,l'} ≪ W_{i,l'-1}", a1_w_cross);
    show("A_{1,l'} ≪ A_{i,l'-1}", a1_a_cross);
    show("D_k^+(col1) ≪ W_{i,col2}", dk_w_cross);
    show("U_{1,l'-1} ≪ W_{i,l'} (other dir)", u1_w_cross2);
}
