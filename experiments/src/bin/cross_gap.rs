//! Test the exact shift-lemma condition for cross-group DU.
//! h = Σ_{i≥j'} D_{i,col2} ≪ f = Σ_{i<k'} A_{i,col1} + Σ_{i≥j'} A_{i,col2}
//! where col1 = l', col2 = l'-1 (adjacent columns).
//! Also test: h ≪ Σ_{i<k'} A_{i,col1} alone (cross-column part).
use polynomial_tools::real_rootedness::check_weak_interlacing;
use std::collections::BTreeSet;
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
        None => match (pdeg(&f), pdeg(&g)) {
            (Some(df), Some(dg)) if df == dg => {
                let tf = pmt(&f);
                check_weak_interlacing(&g, &tf).unwrap_or(false)
            }
            _ => false,
        },
    }
}
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
fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);
    let mut ground = [0u64; 2]; // cross-group DU at lambda+
    let mut h_f = [0u64; 2]; // h ≪ f_s (full shift-lemma condition, s=1)
    let mut h_cross_a = [0u64; 2]; // h ≪ Σ A_{i,col1} (cross-col part alone)
    let mut h_same_a = [0u64; 2]; // h ≪ Σ A_{i,col2} (same-col part alone)
    let mut h_u1_col1 = [0u64; 2]; // h ≪ U_{1,col1} (= A_{1,col1})
    let mut dl_u1 = [0u64; 2]; // D_{j',col2} ≪ U_{1,col1} (single cross-col DU at j=1)
    for n in 1..=max_n {
        for board in &gen_boards(n) {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            let m = board[0] as usize;
            if n <= 2 {
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
                let pk = peaks(pi);
                let poly = if pi.len() >= 2 && pi[0] > pi[1] {
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
            // Precompute prefix/suffix sums for each column
            let mut s_col: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; mp + 2]; m + 1]; // S_k^(col) = Σ_{j<k} A_j
            let mut t_col: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; mp + 2]; m + 1]; // T_k^(col) = Σ_{j≥k} W_j
            for col in 1..=m {
                for k in 1..=m {
                    while s_col[col].len() <= k + 1 {
                        s_col[col].push(vec![0i64]);
                    }
                    s_col[col][k + 1] = pa(&s_col[col][k], &a_jl[k][col]);
                }
                while t_col[col].len() <= mp + 1 {
                    t_col[col].push(vec![0i64]);
                }
                for k in (1..=m).rev() {
                    t_col[col][k] = pa(&t_col[col][k + 1], &w_jl[k][col]);
                }
            }
            // Test cross-group pairs for each l' with 2 ≤ l' ≤ m
            for lp in 2..=m {
                // l' is the last entry at lambda+ level
                let col1 = lp; // upper group column
                let col2 = lp - 1; // lower group column
                if col1 > m || col2 < 1 {
                    continue;
                }
                // Upper group: k' > l', so k' = l'+1, ..., m+1
                // Lower group: j' < l', so j' = 1, ..., l'-1
                for kp in (lp + 1)..=mp {
                    for jp in 1..lp {
                        // D_{k'}^+(col1) = S_{k'}^(col1), U_{j'}^+(col2) = T_{j'}^(col2)
                        let d_plus = if kp < s_col[col1].len() {
                            &s_col[col1][kp]
                        } else {
                            &s_col[col1][m + 1]
                        };
                        let u_plus = if jp < t_col[col2].len() {
                            &t_col[col2][jp]
                        } else {
                            &vec![0i64]
                        };
                        // Ground truth
                        if !pz(d_plus) && !pz(u_plus) {
                            ground[0] += 1;
                            if !interlaces(d_plus, u_plus) {
                                ground[1] += 1;
                            }
                        }
                        // h = Σ_{i≥j'} D_{i,col2}
                        let mut h = vec![0i64];
                        for i in jp..=m {
                            h = pa(&h, &d_jl[i][col2]);
                        }
                        if pz(&h) {
                            continue;
                        }
                        // f = Σ_{i<k'} A_{i,col1} + Σ_{i≥j'} A_{i,col2}
                        let cross_part = if kp < s_col[col1].len() {
                            s_col[col1][kp].clone()
                        } else {
                            s_col[col1][m + 1].clone()
                        };
                        let mut same_part = vec![0i64];
                        for i in jp..=m {
                            same_part = pa(&same_part, &a_jl[i][col2]);
                        }
                        let f = pa(&cross_part, &same_part);
                        // h ≪ f (full condition)
                        if !pz(&f) {
                            h_f[0] += 1;
                            if !interlaces(&h, &f) {
                                h_f[1] += 1;
                            }
                        }
                        // h ≪ cross_part (= Σ A_{i,col1})
                        if !pz(&cross_part) {
                            h_cross_a[0] += 1;
                            if !interlaces(&h, &cross_part) {
                                h_cross_a[1] += 1;
                            }
                        }
                        // h ≪ same_part (= Σ A_{i,col2})
                        if !pz(&same_part) {
                            h_same_a[0] += 1;
                            if !interlaces(&h, &same_part) {
                                h_same_a[1] += 1;
                            }
                        }
                        // h ≪ U_{1,col1} = A_{1,col1}
                        if !pz(&a_jl[1][col1]) {
                            h_u1_col1[0] += 1;
                            if !interlaces(&h, &a_jl[1][col1]) {
                                h_u1_col1[1] += 1;
                            }
                        }
                    }
                }
                // Also test D_{j',col2} ≪ U_{1,col1} for all j'
                for jp in 1..=m {
                    if !pz(&d_jl[jp][col2]) && !pz(&u_jl[1][col1]) {
                        dl_u1[0] += 1;
                        if !interlaces(&d_jl[jp][col2], &u_jl[1][col1]) {
                            dl_u1[1] += 1;
                        }
                    }
                }
            }
        }
    }
    println!("=== Cross-group gap analysis (n <= {}) ===", max_n);
    let show = |name: &str, c: [u64; 2]| {
        if c[0] == 0 {
            println!("  {}: (no data)", name);
        } else if c[1] == 0 {
            println!("  {}: {}/{} ALL PASS <<<", name, c[0], c[0]);
        } else {
            println!("  {}: {}/{} pass ({} FAIL)", name, c[0] - c[1], c[0], c[1]);
        }
    };
    show("Ground truth: D^+(col1) ≪ U^+(col2)", ground);
    show("Shift lemma: h ≪ f_s (full, s=1)", h_f);
    show("h ≪ Σ A_{col1} (cross-col part)", h_cross_a);
    show("h ≪ Σ A_{col2} (same-col part)", h_same_a);
    show("h ≪ U_{1,col1} (= A_{1,col1})", h_u1_col1);
    show("D_{j',col2} ≪ U_{1,col1}", dl_u1);
}
