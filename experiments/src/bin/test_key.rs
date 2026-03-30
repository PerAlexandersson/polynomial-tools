// The KEY identity: x*D_k^+ + y*U_l^+ = Σ b_j A_j + y(t-1)·Σ_{j≥l} D_j
// where b_j = x*[j<k] + y*[j≥l].
// So DU at λ^+ iff for all k,l: Σ_{j≥l} D_j ≪ Σ b_j A_j and boundary holds.
// Test: does Σ_{j≥l} D_j ≪ Σ b_j A_j hold for ALL k, l, x=1, y=1?
// (This is the shift lemma condition for DU.)
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
    // Key identity: D_k^+ + U_l^+ = Σ b_j A_j + (t-1)·Σ_{j≥l} D_j
    // where b_j = [j<k] + [j≥l].
    // Shift lemma: D_k^+ ≪ D_k^+ + U_l^+ (i.e., DU) iff
    //   (1) Σ_{j≥l} D_j ≪ D_k^+ + U_l^+ and (2) boundary holds.
    // But we can also write: D_k^+ + U_l^+ = P + (t-1)h where
    //   P = D_k^+ + U_l^+ - (t-1)h = Σ b_j A_j and h = Σ_{j≥l} D_j.
    // Actually: for the DU condition D_k^+ ≪ U_l^+, we consider
    //   U_l^+ = D_k^+ + [(t-1) partial D terms + rest].
    // Let me instead directly verify: h = Σ_{j≥l} D_j ≪ Σ b_j A_j 
    // where b_j = [j<k] + [j≥l].
    let mut h_ba = [0u64; 2]; // h ≪ Σ b_j A_j for all k, l
    let mut h_ba_overlap = [0u64; 2]; // same but only k > l (overlapping case)
    // Also test D_j ≪ Σ b_j A_j for each individual D_j (j ≥ l)
    let mut each_d_ba = [0u64; 2];
    for n in 1..=max_n {
        for board in &gen_boards(n) {
            let perm = board_to_perm(board); if !is_312_avoiding(&perm) { continue; }
            let m = board[0] as usize; if n <= 1 { continue; }
            let ideal = bruhat_lower_ideal(&perm);
            let mut d_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m+1]; m+1];
            let mut u_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m+1]; m+1];
            for pi in &ideal { let j = pi[0] as usize; if j > m { continue; } let l = *pi.last().unwrap() as usize; if l > m { continue; } let pk = peaks(pi);
                let poly = if pi.len() >= 2 && pi[0] > pi[1] { &mut d_jl[j][l] } else { &mut u_jl[j][l] }; while poly.len() <= pk { poly.push(0); } poly[pk] += 1; }
            let mut a_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m+1]; m+1];
            for j in 1..=m { for l in 1..=m { a_jl[j][l] = pa(&d_jl[j][l], &u_jl[j][l]); } }
            for col in 1..=m {
                for k in 1..=m+1 { // k ranges 1..=m+1 for the output
                    for l in 1..=m+1 {
                        if k == l { continue; }
                        // b_j = [j<k] + [j≥l] for j in 1..=m
                        let mut ba = vec![0i64]; // Σ b_j A_j
                        for j in 1..=m {
                            let b = (if j < k { 1 } else { 0 }) + (if j >= l { 1 } else { 0 });
                            if b > 0 {
                                let scaled: Vec<i64> = a_jl[j][col].iter().map(|&c| c * b as i64).collect();
                                ba = pa(&ba, &scaled);
                            }
                        }
                        // h = Σ_{j≥l} D_j
                        let mut h = vec![0i64];
                        for j in l.max(1)..=m { h = pa(&h, &d_jl[j][col]); }
                        if pz(&h) || pz(&ba) { continue; }
                        h_ba[0] += 1;
                        if !interlaces(&h, &ba) { h_ba[1] += 1; }
                        if k > l {
                            h_ba_overlap[0] += 1;
                            if !interlaces(&h, &ba) { h_ba_overlap[1] += 1; }
                        }
                        // Each D_j (j >= l) ≪ Σ b_j A_j
                        for j in l.max(1)..=m {
                            if !pz(&d_jl[j][col]) && !pz(&ba) {
                                each_d_ba[0] += 1;
                                if !interlaces(&d_jl[j][col], &ba) { each_d_ba[1] += 1; }
                            }
                        }
                    }
                }
            }
        }
    }
    println!("h = Σ D_{{j≥l}} ≪ Σ b_j A_j (ALL k,l):    {}/{} pass ({} FAIL)", h_ba[0]-h_ba[1], h_ba[0], h_ba[1]);
    println!("h ≪ Σ b_j A_j (overlap k>l only):   {}/{} pass ({} FAIL)", h_ba_overlap[0]-h_ba_overlap[1], h_ba_overlap[0], h_ba_overlap[1]);
    println!("each D_j ≪ Σ b_j A_j:               {}/{} pass ({} FAIL)", each_d_ba[0]-each_d_ba[1], each_d_ba[0], each_d_ba[1]);
}
