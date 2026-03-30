use polynomial_tools::real_rootedness::check_weak_interlacing;
use std::collections::BTreeSet;
fn pt(p: &[i64]) -> Vec<i64> { let mut v = p.to_vec(); while v.len() > 1 && *v.last().unwrap() == 0 { v.pop(); } v }
fn pz(p: &[i64]) -> bool { p.iter().all(|&c| c == 0) }
fn pa(a: &[i64], b: &[i64]) -> Vec<i64> { let l = a.len().max(b.len()); let mut r = vec![0i64; l]; for (i, &v) in a.iter().enumerate() { r[i] += v; } for (i, &v) in b.iter().enumerate() { r[i] += v; } pt(&r) }
fn pmt(p: &[i64]) -> Vec<i64> { let mut r = vec![0i64; p.len() + 1]; for (i, &v) in p.iter().enumerate() { r[i + 1] = v; } pt(&r) }
fn pdeg(p: &[i64]) -> Option<usize> { let v = pt(p); if pz(&v) { None } else { Some(v.len() - 1) } }
fn interlaces(f: &[i64], g: &[i64]) -> bool {
    let f = pt(f); let g = pt(g); if pz(&f) { return true; } if pz(&g) { return false; }
    match check_weak_interlacing(&f, &g) { Some(true) => true, Some(false) => false,
        None => { match (pdeg(&f), pdeg(&g)) { (Some(df), Some(dg)) if df == dg => { let tf = pmt(&f); check_weak_interlacing(&g, &tf).unwrap_or(false) }, _ => false, } } }
}
fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> { let n = perm.len(); let mut vis: BTreeSet<Vec<u8>> = BTreeSet::new(); let mut q: BTreeSet<Vec<u8>> = BTreeSet::new(); q.insert(perm.to_vec()); while let Some(cur) = q.pop_last() { for i in 0..n { for j in i+1..n { if cur[i] > cur[j] { let mut c = cur.clone(); c.swap(i, j); if !vis.contains(&c) { q.insert(c); } } } } vis.insert(cur); } vis.into_iter().collect() }
fn board_to_perm(b: &[u8]) -> Vec<u8> { let n = b.len(); let mut p = vec![0u8; n]; let mut u = vec![false; n+1]; for i in 0..n { for c in (1..=(b[i] as usize).min(n)).rev() { if !u[c] { p[i] = c as u8; u[c] = true; break; } } } p }
fn is_312_avoiding(perm: &[u8]) -> bool { let n = perm.len(); for i in 0..n { for j in i+1..n { for k in j+1..n { if perm[k] < perm[i] && perm[i] < perm[j] { return false; } } } } true }
fn peaks(w: &[u8]) -> usize { if w.len() < 3 { return 0; } (1..w.len()-1).filter(|&i| w[i-1] < w[i] && w[i] > w[i+1]).count() }
fn gen_boards(n: usize) -> Vec<Vec<u8>> { let mut r = vec![]; let mut c = vec![]; gb(n, n, 0, &mut c, &mut r); r }
fn gb(n: usize, mx: usize, d: usize, c: &mut Vec<u8>, r: &mut Vec<Vec<u8>>) { if d == n { r.push(c.clone()); return; } for v in (d+1).max(if d > 0 { c[d-1] as usize } else { 1 })..=mx { c.push(v as u8); gb(n, mx, d+1, c, r); c.pop(); } }
fn main() {
    let max_n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(7);
    let mut td_u_same = [0u64; 2]; // tD_{j,l} << U_{j,l} (same index)
    let mut td_u_all = [0u64; 2]; // tD_{j,l} << U_{j',l} (all j, j')
    let mut td_w = [0u64; 2]; // tD_{j,l} << W_{j',l} (all j, j')
    // Forward AW via transitivity: A_j << tD_{j'} << W_{j'} for j <= j'
    // Need: A_j << tD_{j'} for j <= j'
    let mut a_td = [0u64; 2]; // A_{j,l} << tD_{j',l} for j <= j'
    for n in 1..=max_n {
        for board in &gen_boards(n) {
            let perm = board_to_perm(board); if !is_312_avoiding(&perm) { continue; }
            let m = board[0] as usize; if n <= 1 { continue; }
            let ideal = bruhat_lower_ideal(&perm);
            let mut d_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m+1]; m+1];
            let mut u_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m+1]; m+1];
            for pi in &ideal { let j = pi[0] as usize; if j > m { continue; } let l = *pi.last().unwrap() as usize; if l > m { continue; } let pk = peaks(pi);
                let poly = if pi.len() >= 2 && pi[0] > pi[1] { &mut d_jl[j][l] } else { &mut u_jl[j][l] }; while poly.len() <= pk { poly.push(0); } poly[pk] += 1; }
            let mut w_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m+1]; m+1];
            for j in 1..=m { for l in 1..=m { w_jl[j][l] = pa(&pmt(&d_jl[j][l]), &u_jl[j][l]); } }
            for l in 1..=m {
                for j in 1..=m {
                    let td = pmt(&d_jl[j][l]);
                    // tD_j << U_j (same index)
                    if !pz(&td) && !pz(&u_jl[j][l]) { td_u_same[0] += 1; if !interlaces(&td, &u_jl[j][l]) { td_u_same[1] += 1; } }
                    for jp in 1..=m {
                        // tD_j << U_{j'} (all)
                        if !pz(&td) && !pz(&u_jl[jp][l]) { td_u_all[0] += 1; if !interlaces(&td, &u_jl[jp][l]) { td_u_all[1] += 1; } }
                        // tD_j << W_{j'} (all)
                        if !pz(&td) && !pz(&w_jl[jp][l]) { td_w[0] += 1; if !interlaces(&td, &w_jl[jp][l]) { td_w[1] += 1; } }
                        // A_j << tD_{j'} for j <= j'
                        let a = pa(&d_jl[j][l], &u_jl[j][l]);
                        let tdp = pmt(&d_jl[jp][l]);
                        if j <= jp && !pz(&a) && !pz(&tdp) { a_td[0] += 1; if !interlaces(&a, &tdp) { a_td[1] += 1; } }
                    }
                }
            }
        }
    }
    println!("tD_j << U_j (same j): {}/{} pass ({} FAIL)", td_u_same[0]-td_u_same[1], td_u_same[0], td_u_same[1]);
    println!("tD_j << U_{{j'}} (all): {}/{} pass ({} FAIL)", td_u_all[0]-td_u_all[1], td_u_all[0], td_u_all[1]);
    println!("tD_j << W_{{j'}} (all): {}/{} pass ({} FAIL)", td_w[0]-td_w[1], td_w[0], td_w[1]);
    println!("A_j << tD_{{j'}} (j<=j'): {}/{} pass ({} FAIL)", a_td[0]-a_td[1], a_td[0], a_td[1]);
}
