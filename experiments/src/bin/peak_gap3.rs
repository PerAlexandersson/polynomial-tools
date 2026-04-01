//! Corrected DU proof gap investigation.
//! Only tests boards where π_λ is 312-avoiding (as required by the theorem).
//! Uses float-based root checking with correct paper convention.

use std::collections::BTreeSet;

fn find_real_roots(coeffs: &[i64]) -> Option<Vec<f64>> {
    let roots = polynomial_tools::real_roots(coeffs)?;
    let mut float_roots: Vec<f64> = roots
        .iter()
        .map(|r| {
            let n = r.numer().to_str_radix(10).parse::<f64>().unwrap();
            let d = r.denom().to_str_radix(10).parse::<f64>().unwrap();
            n / d
        })
        .collect();
    float_roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Some(float_roots)
}

/// Paper convention: f ≪ g means f roots more negative than g roots.
/// Same degree: f[0] ≤ g[0] ≤ f[1] ≤ g[1] ≤ ...
/// deg(f)+1 = deg(g): g[0] ≤ f[0] ≤ g[1] ≤ ... ≤ f[d-1] ≤ g[d]
fn interlaces(f: &[i64], g: &[i64]) -> bool {
    let f = pt(f);
    let g = pt(g);
    if pz(&f) {
        return true;
    }
    if pz(&g) {
        return false;
    }
    let fr = match find_real_roots(&f) {
        Some(r) => r,
        None => return false,
    };
    let gr = match find_real_roots(&g) {
        Some(r) => r,
        None => return false,
    };
    let df = fr.len();
    let dg = gr.len();
    if dg != df && dg != df + 1 {
        return false;
    }
    let eps = 1e-6;
    if df == 0 {
        return true;
    }
    if dg == df + 1 {
        for i in 0..df {
            if gr[i] > fr[i] + eps || fr[i] > gr[i + 1] + eps {
                return false;
            }
        }
        true
    } else {
        for i in 0..df {
            if fr[i] > gr[i] + eps {
                return false;
            }
            if i + 1 < df && gr[i] > fr[i + 1] + eps {
                return false;
            }
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
fn pf(p: &[i64]) -> String {
    let p = pt(p);
    if pz(&p) {
        return "0".into();
    }
    let mut t = vec![];
    for (i, &c) in p.iter().enumerate() {
        if c == 0 {
            continue;
        }
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
    println!("=== DU gap (312-avoiding boards only), n ≤ {} ===\n", max_n);

    let mut total_boards = 0u64;
    let mut valid_boards = 0u64;

    // [tests, fails]
    let mut aw_ord = [0u64; 2];
    let mut ww_ord = [0u64; 2];
    let mut du_all = [0u64; 2];
    let mut uu_ord = [0u64; 2];
    let mut aw_ext = [0u64; 2]; // A_j ≪ W_l for j > l (the extension)
    let mut dw_all = [0u64; 2];
    let mut dd_fwd = [0u64; 2];
    let mut dd_rev = [0u64; 2];
    let mut uu_rev = [0u64; 2];
    let mut st_all = [0u64; 2];
    let mut st_gap = [0u64; 2]; // S_j ≪ T_l for j > l only (THE GAP)
    let mut aw_plus_all = [0u64; 2];

    for n in 1..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            total_boards += 1;
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            valid_boards += 1;

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

            let mut s = vec![vec![0i64]; m + 2];
            let mut t_arr = vec![vec![0i64]; m + 2];
            for k in 2..=m + 1 {
                s[k] = pa(&s[k - 1], &a[k - 1]);
            }
            for k in (1..=m).rev() {
                t_arr[k] = pa(&t_arr[k + 1], &w[k]);
            }

            // Theorem conditions
            for j in 1..=m {
                for l in j..=m {
                    if !pz(&a[j]) && !pz(&w[l]) {
                        aw_ord[0] += 1;
                        if !interlaces(&a[j], &w[l]) {
                            aw_ord[1] += 1;
                            if aw_ord[1] <= 3 {
                                println!(
                                    "FAIL AW({},{}) {:?}: {} vs {}",
                                    j,
                                    l,
                                    board,
                                    pf(&a[j]),
                                    pf(&w[l])
                                );
                            }
                        }
                    }
                    if j < l && !pz(&w[j]) && !pz(&w[l]) {
                        ww_ord[0] += 1;
                        if !interlaces(&w[j], &w[l]) {
                            ww_ord[1] += 1;
                            if ww_ord[1] <= 3 {
                                println!("FAIL WW({},{}) {:?}", j, l, board);
                            }
                        }
                    }
                    if !pz(&u[j]) && !pz(&u[l]) {
                        uu_ord[0] += 1;
                        if !interlaces(&u[j], &u[l]) {
                            uu_ord[1] += 1;
                            if uu_ord[1] <= 3 {
                                println!(
                                    "FAIL UU({},{}) {:?}: {} vs {}",
                                    j,
                                    l,
                                    board,
                                    pf(&u[j]),
                                    pf(&u[l])
                                );
                            }
                        }
                    }
                }
            }
            for j in 1..=m {
                for l in 1..=m {
                    if !pz(&d[j]) && !pz(&u[l]) {
                        du_all[0] += 1;
                        if !interlaces(&d[j], &u[l]) {
                            du_all[1] += 1;
                            if du_all[1] <= 3 {
                                println!(
                                    "FAIL DU({},{}) {:?}: {} vs {}",
                                    j,
                                    l,
                                    board,
                                    pf(&d[j]),
                                    pf(&u[l])
                                );
                            }
                        }
                    }
                }
            }

            // Potential strengthenings
            for j in 1..=m {
                for l in 1..j {
                    if !pz(&a[j]) && !pz(&w[l]) {
                        aw_ext[0] += 1;
                        if !interlaces(&a[j], &w[l]) {
                            aw_ext[1] += 1;
                            if aw_ext[1] <= 5 {
                                println!(
                                    "FAIL AW_ext({},{}) {:?}: {} vs {}",
                                    j,
                                    l,
                                    board,
                                    pf(&a[j]),
                                    pf(&w[l])
                                );
                            }
                        }
                    }
                }
            }
            for j in 1..=m {
                for l in 1..=m {
                    if !pz(&d[j]) && !pz(&w[l]) {
                        dw_all[0] += 1;
                        if !interlaces(&d[j], &w[l]) {
                            dw_all[1] += 1;
                            if dw_all[1] <= 3 {
                                println!("FAIL DW({},{}) {:?}", j, l, board);
                            }
                        }
                    }
                }
            }
            for j in 1..=m {
                for l in j + 1..=m {
                    if !pz(&d[j]) && !pz(&d[l]) {
                        dd_fwd[0] += 1;
                        if !interlaces(&d[j], &d[l]) {
                            dd_fwd[1] += 1;
                        }
                        dd_rev[0] += 1;
                        if !interlaces(&d[l], &d[j]) {
                            dd_rev[1] += 1;
                        }
                    }
                    if !pz(&u[j]) && !pz(&u[l]) {
                        // we already tested j<l in UU, now test reverse
                    }
                }
            }
            for j in 1..=m {
                for l in 1..j {
                    if !pz(&u[j]) && !pz(&u[l]) {
                        uu_rev[0] += 1;
                        if !interlaces(&u[j], &u[l]) {
                            uu_rev[1] += 1;
                            if uu_rev[1] <= 3 {
                                println!(
                                    "FAIL UU_rev({},{}) {:?}: {} vs {}",
                                    j,
                                    l,
                                    board,
                                    pf(&u[j]),
                                    pf(&u[l])
                                );
                            }
                        }
                    }
                }
            }

            // S_j ≪ T_l for all j,l
            for j in 1..=m {
                for l in 1..=m {
                    if !pz(&s[j]) && !pz(&t_arr[l]) {
                        st_all[0] += 1;
                        if !interlaces(&s[j], &t_arr[l]) {
                            st_all[1] += 1;
                            if st_all[1] <= 3 {
                                println!(
                                    "FAIL S≪T({},{}) {:?}: {} vs {}",
                                    j,
                                    l,
                                    board,
                                    pf(&s[j]),
                                    pf(&t_arr[l])
                                );
                            }
                        }
                        if j > l {
                            st_gap[0] += 1;
                            if !interlaces(&s[j], &t_arr[l]) {
                                st_gap[1] += 1;
                            }
                        }
                    }
                }
            }

            // A_j^+ ≪ W_l^+ for ALL j,l
            for mp in 1..=m {
                for j in 1..=mp {
                    for l in 1..=mp {
                        let ap = pa(&s[j], &t_arr[j]);
                        let wp = pa(&pmt(&s[l]), &t_arr[l]);
                        if !pz(&ap) && !pz(&wp) {
                            aw_plus_all[0] += 1;
                            if !interlaces(&ap, &wp) {
                                aw_plus_all[1] += 1;
                                if aw_plus_all[1] <= 3 {
                                    println!("FAIL AW+({},{}) {:?} m'={}", j, l, board, mp);
                                }
                            }
                        }
                    }
                }
            }
        }
        println!(
            "n={}: {}/{} boards 312-avoiding",
            n, valid_boards, total_boards
        );
    }

    println!("\n=== RESULTS (312-avoiding boards only) ===\n");
    let show = |name: &str, c: [u64; 2]| {
        if c[0] == 0 {
            println!("  {}: (no tests)", name);
        } else {
            println!(
                "  {}: {}/{} pass {}",
                name,
                c[0] - c[1],
                c[0],
                if c[1] == 0 { "✓" } else { "✗ FAIL" }
            );
        }
    };

    println!("Paper's theorem (should all pass):");
    show("(a) A_j ≪ W_l (j≤l)", aw_ord);
    show("(b) W_j ≪ W_l (j≤l)", ww_ord);
    show("(c) D_j ≪ U_l (all j,l)", du_all);
    show("(d) U_j ≪ U_l (j≤l)", uu_ord);

    println!("\nPotential strengthenings:");
    show("A_j ≪ W_l (j>l, extension)", aw_ext);
    show("D_j ≪ W_l (all j,l)", dw_all);
    show("D_j ≪ D_l (j≤l, fwd)", dd_fwd);
    show("D_l ≪ D_j (j<l, rev)", dd_rev);
    show("U_j ≪ U_l (j>l, rev)", uu_rev);

    println!("\nAt level λ+:");
    show("S_j ≪ T_l (all j,l)", st_all);
    show("S_j ≪ T_l (j>l, THE GAP)", st_gap);
    show("A_j^+ ≪ W_l^+ (all j,l)", aw_plus_all);
}
