//! Corrected DU proof gap investigation using float-based root checking.
//!
//! The exact_interlaces function in other binaries has a bug: it uses
//! f0.division(&h) which computes the REMAINDER (= 0 when h | f),
//! so coprime polynomials always "pass" trivially.
//!
//! This version computes actual roots and checks interlacing directly.

use std::collections::BTreeSet;

fn find_real_roots(coeffs: &[i64]) -> Option<Vec<f64>> {
    let roots = polynomial_tools::real_roots(coeffs)?;
    let mut float_roots: Vec<f64> = roots.iter().map(|r| {
        let n = r.numer().to_str_radix(10).parse::<f64>().unwrap();
        let d = r.denom().to_str_radix(10).parse::<f64>().unwrap();
        n / d
    }).collect();
    float_roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Some(float_roots)
}

/// Check f ≪ g (paper convention): roots of f are more negative than roots of g.
/// For same degree d: f[0] ≤ g[0] ≤ f[1] ≤ g[1] ≤ ... ≤ f[d-1] ≤ g[d-1]
/// For deg(f)+1 = deg(g): g[0] ≤ f[0] ≤ g[1] ≤ f[1] ≤ ... ≤ f[d-1] ≤ g[d]
fn interlaces(f: &[i64], g: &[i64]) -> bool {
    let f = poly_trim(f);
    let g = poly_trim(g);
    if pz(&f) { return true; }
    if pz(&g) { return false; }
    let fr = match find_real_roots(&f) { Some(r) => r, None => return false };
    let gr = match find_real_roots(&g) { Some(r) => r, None => return false };
    let df = fr.len();
    let dg = gr.len();
    if dg != df && dg != df + 1 { return false; }
    let eps = 1e-6;
    if df == 0 { return true; }
    if dg == df + 1 {
        // g has one more root: g[0] ≤ f[0] ≤ g[1] ≤ ... ≤ f[df-1] ≤ g[dg-1]
        for i in 0..df {
            if gr[i] > fr[i] + eps || fr[i] > gr[i + 1] + eps { return false; }
        }
        true
    } else {
        // Same degree: f[0] ≤ g[0] ≤ f[1] ≤ g[1] ≤ ... ≤ f[df-1] ≤ g[df-1]
        // (paper convention: f roots more negative)
        for i in 0..df {
            if fr[i] > gr[i] + eps { return false; }
            if i + 1 < df && gr[i] > fr[i + 1] + eps { return false; }
        }
        true
    }
}

fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> {
    let n = perm.len();
    let mut vis: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut q: BTreeSet<Vec<u8>> = BTreeSet::new();
    q.insert(perm.to_vec());
    while let Some(cur) = q.pop_last() {
        for i in 0..n { for j in i+1..n { if cur[i]>cur[j] {
            let mut c=cur.clone(); c.swap(i,j);
            if !vis.contains(&c) { q.insert(c); }
        }}}
        vis.insert(cur);
    }
    vis.into_iter().collect()
}

fn board_to_perm(b: &[u8]) -> Vec<u8> {
    let n=b.len(); let mut p=vec![0u8;n]; let mut u=vec![false;n+1];
    for i in 0..n { for c in (1..=(b[i] as usize).min(n)).rev() {
        if !u[c] { p[i]=c as u8; u[c]=true; break; }
    }} p
}

fn peaks(w: &[u8]) -> usize {
    if w.len()<3 { return 0; }
    (1..w.len()-1).filter(|&i| w[i-1]<w[i] && w[i]>w[i+1]).count()
}

fn poly_trim(p: &[i64]) -> Vec<i64> {
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
    poly_trim(&r)
}
fn pmt(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; p.len() + 1];
    for (i, &v) in p.iter().enumerate() { r[i + 1] = v; }
    poly_trim(&r)
}
fn pf(p: &[i64]) -> String {
    let p = poly_trim(p);
    if pz(&p) { return "0".into(); }
    let mut t = vec![];
    for (i, &c) in p.iter().enumerate() {
        if c == 0 { continue; }
        match (c, i) {
            (c, 0) => t.push(format!("{}", c)),
            (1, 1) => t.push("t".into()),
            (-1, 1) => t.push("-t".into()),
            (c, 1) => t.push(format!("{}t", c)),
            (1, e) => t.push(format!("t^{}", e)),
            (-1, e) => t.push(format!("-t^{}", e)),
            (c, e) => t.push(format!("{}t^{}", c, e)),
        }
    }
    t.join(" + ")
}

fn gen_boards(n: usize) -> Vec<Vec<u8>> {
    let mut r = vec![]; let mut c = vec![];
    gb(n, n, 0, &mut c, &mut r); r
}
fn gb(n: usize, mx: usize, d: usize, c: &mut Vec<u8>, r: &mut Vec<Vec<u8>>) {
    if d == n { r.push(c.clone()); return; }
    for v in (d+1).max(if d > 0 { c[d-1] as usize } else { 1 })..=mx {
        c.push(v as u8); gb(n, mx, d+1, c, r); c.pop();
    }
}

struct BD {
    m: usize,
    d: Vec<Vec<i64>>,
    u: Vec<Vec<i64>>,
    a: Vec<Vec<i64>>,
    w: Vec<Vec<i64>>,
}

fn compute(board: &[u8]) -> BD {
    let n = board.len();
    let m = (board[0] as usize).min(n);
    let perm = board_to_perm(board);
    let ideal = bruhat_lower_ideal(&perm);
    let mut d = vec![vec![0i64]; m+1];
    let mut u = vec![vec![0i64]; m+1];
    for pi in &ideal {
        if pi.len() < 2 {
            let k = pi[0] as usize;
            if k <= m { while u[k].len() < 1 { u[k].push(0); } u[k][0] += 1; }
            continue;
        }
        let k = pi[0] as usize; if k > m { continue; }
        let pk = peaks(pi);
        let poly = if pi[0] > pi[1] { &mut d[k] } else { &mut u[k] };
        while poly.len() <= pk { poly.push(0); } poly[pk] += 1;
    }
    let mut a = vec![vec![0i64]; m+1];
    let mut w = vec![vec![0i64]; m+1];
    for k in 1..=m { a[k] = pa(&d[k], &u[k]); w[k] = pa(&pmt(&d[k]), &u[k]); }
    BD { m, d, u, a, w }
}

fn main() {
    let max_n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(6);
    println!("=== DU proof gap (corrected float-based checker), n ≤ {} ===\n", max_n);

    // First: convention sanity check
    println!("--- Convention check ---");
    let f = vec![1, 1]; // 1+t, root -1
    let g = vec![0, 1]; // t, root 0
    println!("f=1+t (root -1), g=t (root 0)");
    println!("  interlaces(f,g) = {} (paper f≪g: should be TRUE)", interlaces(&f, &g));
    println!("  interlaces(g,f) = {} (paper g≪f: should be FALSE)", interlaces(&g, &f));
    let f2 = vec![3, 4, 1]; // (t+1)(t+3), roots -3,-1
    let g2 = vec![8, 6, 1]; // (t+2)(t+4), roots -4,-2
    println!("f=(t+1)(t+3) roots [-3,-1], g=(t+2)(t+4) roots [-4,-2]");
    println!("  interlaces(f,g) = {} (paper f≪g: should be FALSE)", interlaces(&f2, &g2));
    println!("  interlaces(g,f) = {} (paper g≪f: should be TRUE)", interlaces(&g2, &f2));
    println!();

    // Counters: [tests, fails]
    let mut aw_ord = [0u64; 2];  // A_j ≪ W_l for j ≤ l (paper AW)
    let mut ww_ord = [0u64; 2];  // W_j ≪ W_l for j ≤ l (paper WW)
    let mut du_all = [0u64; 2];  // D_j ≪ U_l for ALL j,l (paper DU)
    let mut uu_ord = [0u64; 2];  // U_j ≪ U_l for j ≤ l (paper UU)
    let mut aw_all = [0u64; 2];  // A_j ≪ W_l for ALL j,l (strengthened AW)
    let mut dw_all = [0u64; 2];  // D_j ≪ W_l for ALL j,l
    let mut dd_fwd = [0u64; 2];  // D_j ≪ D_l for j ≤ l
    let mut dd_rev = [0u64; 2];  // D_l ≪ D_j for j < l
    let mut uu_rev = [0u64; 2];  // U_j ≪ U_l for j > l (reverse UU)
    let mut st_all = [0u64; 2];  // S_j ≪ T_l for ALL j,l (DU at λ+)
    let mut aw_plus = [0u64; 2]; // A_j^+ ≪ W_l^+ for ALL j,l

    for n in 1..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let bd = compute(board);
            let m = bd.m;

            // Build S_k, T_k
            let mut s = vec![vec![0i64]; m+2];
            let mut t = vec![vec![0i64]; m+2];
            for k in 2..=m+1 { s[k] = pa(&s[k-1], &bd.a[k-1]); }
            for k in (1..=m).rev() { t[k] = pa(&t[k+1], &bd.w[k]); }

            // Paper's theorem conditions
            for j in 1..=m { for l in j..=m {
                if !pz(&bd.a[j]) && !pz(&bd.w[l]) {
                    aw_ord[0] += 1;
                    if !interlaces(&bd.a[j], &bd.w[l]) { aw_ord[1] += 1;
                        if aw_ord[1] <= 3 { println!("FAIL AW: A_{}≪W_{} board {:?}: {} vs {}", j, l, board, pf(&bd.a[j]), pf(&bd.w[l])); }
                    }
                }
                if j < l && !pz(&bd.w[j]) && !pz(&bd.w[l]) {
                    ww_ord[0] += 1;
                    if !interlaces(&bd.w[j], &bd.w[l]) { ww_ord[1] += 1;
                        if ww_ord[1] <= 3 { println!("FAIL WW: W_{}≪W_{} board {:?}", j, l, board); }
                    }
                }
                if !pz(&bd.u[j]) && !pz(&bd.u[l]) {
                    uu_ord[0] += 1;
                    if !interlaces(&bd.u[j], &bd.u[l]) { uu_ord[1] += 1;
                        if uu_ord[1] <= 3 { println!("FAIL UU: U_{}≪U_{} board {:?}: {} vs {}", j, l, board, pf(&bd.u[j]), pf(&bd.u[l])); }
                    }
                }
            }}
            // DU for all j,l
            for j in 1..=m { for l in 1..=m {
                if !pz(&bd.d[j]) && !pz(&bd.u[l]) {
                    du_all[0] += 1;
                    if !interlaces(&bd.d[j], &bd.u[l]) { du_all[1] += 1;
                        if du_all[1] <= 3 { println!("FAIL DU: D_{}≪U_{} board {:?}: {} vs {}", j, l, board, pf(&bd.d[j]), pf(&bd.u[l])); }
                    }
                }
            }}

            // Strengthened conditions
            for j in 1..=m { for l in 1..=m {
                // AW all
                if !pz(&bd.a[j]) && !pz(&bd.w[l]) {
                    if j > l { // only count the NEW cases
                        aw_all[0] += 1;
                        if !interlaces(&bd.a[j], &bd.w[l]) { aw_all[1] += 1;
                            if aw_all[1] <= 3 { println!("FAIL AW_all: A_{}≪W_{} board {:?}: {} vs {}", j, l, board, pf(&bd.a[j]), pf(&bd.w[l])); }
                        }
                    }
                }
                // DW all
                if !pz(&bd.d[j]) && !pz(&bd.w[l]) {
                    dw_all[0] += 1;
                    if !interlaces(&bd.d[j], &bd.w[l]) { dw_all[1] += 1;
                        if dw_all[1] <= 3 { println!("FAIL DW: D_{}≪W_{} board {:?}: {} vs {}", j, l, board, pf(&bd.d[j]), pf(&bd.w[l])); }
                    }
                }
            }}
            // DD forward and reverse
            for j in 1..=m { for l in j+1..=m {
                if !pz(&bd.d[j]) && !pz(&bd.d[l]) {
                    dd_fwd[0] += 1;
                    if !interlaces(&bd.d[j], &bd.d[l]) { dd_fwd[1] += 1;
                        if dd_fwd[1] <= 3 { println!("FAIL DD_fwd: D_{}≪D_{} board {:?}: {} vs {}", j, l, board, pf(&bd.d[j]), pf(&bd.d[l])); }
                    }
                    dd_rev[0] += 1;
                    if !interlaces(&bd.d[l], &bd.d[j]) { dd_rev[1] += 1;
                        if dd_rev[1] <= 3 { println!("FAIL DD_rev: D_{}≪D_{} board {:?}: {} vs {}", l, j, board, pf(&bd.d[l]), pf(&bd.d[j])); }
                    }
                }
            }}
            // UU reverse
            for j in 1..=m { for l in 1..j {
                if !pz(&bd.u[j]) && !pz(&bd.u[l]) {
                    uu_rev[0] += 1;
                    if !interlaces(&bd.u[j], &bd.u[l]) { uu_rev[1] += 1;
                        if uu_rev[1] <= 3 { println!("FAIL UU_rev: U_{}≪U_{} board {:?}: {} vs {}", j, l, board, pf(&bd.u[j]), pf(&bd.u[l])); }
                    }
                }
            }}

            // S_j ≪ T_l for ALL j,l
            for j in 1..=m { for l in 1..=m {
                if !pz(&s[j]) && !pz(&t[l]) {
                    st_all[0] += 1;
                    if !interlaces(&s[j], &t[l]) { st_all[1] += 1;
                        if st_all[1] <= 3 { println!("FAIL S_{}≪T_{} board {:?}: {} vs {}", j, l, board, pf(&s[j]), pf(&t[l])); }
                    }
                }
            }}

            // A_j^+ ≪ W_l^+ for ALL j,l
            for mp in 1..=m {
                for j in 1..=mp { for l in 1..=mp {
                    let a_plus = pa(&s[j], &t[j]);
                    let w_plus = pa(&pmt(&s[l]), &t[l]);
                    if !pz(&a_plus) && !pz(&w_plus) {
                        aw_plus[0] += 1;
                        if !interlaces(&a_plus, &w_plus) { aw_plus[1] += 1;
                            if aw_plus[1] <= 3 { println!("FAIL AW+: A_{}^+≪W_{}^+ board {:?} m'={}", j, l, board, mp); }
                        }
                    }
                }}
            }
        }
        println!("n={} done", n);
    }

    println!("\n=== RESULTS (paper convention: f≪g means f roots more negative) ===\n");
    let show = |name: &str, c: [u64; 2]| {
        if c[0] == 0 { println!("  {}: (no tests)", name); }
        else { println!("  {}: {}/{} pass {}", name, c[0]-c[1], c[0],
            if c[1]==0 {"✓"} else {"✗ FAIL"}); }
    };

    println!("Paper's theorem conditions:");
    show("(a) A_j ≪ W_l (j≤l)", aw_ord);
    show("(b) W_j ≪ W_l (j≤l)", ww_ord);
    show("(c) D_j ≪ U_l (all j,l)", du_all);
    show("(d) U_j ≪ U_l (j≤l)", uu_ord);

    println!("\nPotential strengthenings:");
    show("A_j ≪ W_l (j>l only, NEW)", aw_all);
    show("D_j ≪ W_l (all j,l)", dw_all);
    show("D_j ≪ D_l (j≤l, fwd)", dd_fwd);
    show("D_l ≪ D_j (j<l, rev)", dd_rev);
    show("U_j ≪ U_l (j>l, rev UU)", uu_rev);

    println!("\nAt level λ+:");
    show("S_j ≪ T_l (all j,l) = DU at λ+", st_all);
    show("A_j^+ ≪ W_l^+ (all j,l)", aw_plus);
}
