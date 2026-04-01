//! DU proof gap investigation with CORRECT interlacing check.
//! Handles shared roots via GCD removal (fixed quotient computation).
//! Only tests 312-avoiding boards (as required by the theorem).

use num::{BigInt, BigRational, Zero};
use num_rational::Ratio;
use polynomial_tools::polynomial::{FieldRing, Polynomial};
use std::collections::BTreeSet;

type Poly = Polynomial<Ratio<BigInt>>;
type BR = Ratio<BigInt>;
fn br(n: i64) -> BR {
    BR::from_integer(BigInt::from(n))
}
fn to_poly(c: &[i64]) -> Poly {
    Polynomial::new(
        c.iter()
            .map(|&x| BR::from_integer(BigInt::from(x)))
            .collect(),
    )
}

/// Proper polynomial quotient: returns f / g (assumes g | f).
fn poly_quotient(f: &Poly, g: &Poly) -> Poly {
    f.exact_div(g)
}

fn poly_long_div(f: &Poly, g: &Poly) -> Poly {
    let (q, _) = f.div_rem(g);
    q
}

fn lagrange_coeffs(points: &[BR], values: &[BR]) -> Vec<BR> {
    let p = Polynomial::lagrange_interpolation(points, values);
    let d = p.degree().unwrap_or(0);
    (0..=d).map(|i| p.coeff(i)).collect()
}

fn polynomial_gcd(f: &Poly, g: &Poly) -> Poly {
    f.gcd(g)
}

/// Check f ≪ g (paper convention) using exact arithmetic with GCD removal.
/// Uses Bézout matrix for degree-diff-1, Sturm chains for same-degree.
fn exact_interlaces(fc: &[i64], gc: &[i64]) -> bool {
    // Try weak interlacing (handles shared roots via GCD, degree diff = 1)
    if let Some(result) = polynomial_tools::check_weak_interlacing(fc, gc) {
        return result;
    }
    // Fallback to Sturm-based check (handles same-degree interlacing)
    polynomial_tools::check_interlacing_sturm(fc, gc).unwrap_or(false)
}

// ── Combinatorial infrastructure ───────────────────────────────────

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
            (c, 1) => t.push(format!("{}t", c)),
            (1, e) => t.push(format!("t^{}", e)),
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
    println!(
        "=== DU gap (exact, GCD-fixed, paper convention), n ≤ {} ===\n",
        max_n
    );

    // Sanity check
    println!("--- Sanity ---");
    // f=1+t ≪ g=t: paper TRUE (root -1 ≤ 0)
    println!("(1+t) ≪ t: {}", exact_interlaces(&[1, 1], &[0, 1]));
    // g=t ≪ f=1+t: paper FALSE
    println!("t ≪ (1+t): {}", exact_interlaces(&[0, 1], &[1, 1]));
    // f=5t^2+6t+1 ≪ g=10t^2+2t (share root -1/5): paper TRUE
    println!(
        "(5t^2+6t+1) ≪ (10t^2+2t): {}",
        exact_interlaces(&[1, 6, 5], &[0, 2, 10])
    );
    println!();

    let mut total = 0u64;
    let mut valid = 0u64;
    let mut aw_ord = [0u64; 2];
    let mut ww_ord = [0u64; 2];
    let mut du_all = [0u64; 2];
    let mut uu_ord = [0u64; 2];
    let mut aw_ext = [0u64; 2];
    let mut dd_fwd = [0u64; 2];
    let mut dd_rev = [0u64; 2];
    let mut uu_rev = [0u64; 2];
    let mut st_all = [0u64; 2];
    let mut st_gap = [0u64; 2];

    for n in 1..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            total += 1;
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            valid += 1;
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
            let mut t_a = vec![vec![0i64]; m + 2];
            for k in 2..=m + 1 {
                s[k] = pa(&s[k - 1], &a[k - 1]);
            }
            for k in (1..=m).rev() {
                t_a[k] = pa(&t_a[k + 1], &w[k]);
            }

            for j in 1..=m {
                for l in j..=m {
                    if !pz(&a[j]) && !pz(&w[l]) {
                        aw_ord[0] += 1;
                        if !exact_interlaces(&a[j], &w[l]) {
                            aw_ord[1] += 1;
                            if aw_ord[1] <= 3 {
                                println!("FAIL AW({},{}) {:?}", j, l, board);
                            }
                        }
                    }
                    if j < l && !pz(&w[j]) && !pz(&w[l]) {
                        ww_ord[0] += 1;
                        if !exact_interlaces(&w[j], &w[l]) {
                            ww_ord[1] += 1;
                            if ww_ord[1] <= 3 {
                                println!("FAIL WW({},{}) {:?}", j, l, board);
                            }
                        }
                    }
                    if !pz(&u[j]) && !pz(&u[l]) {
                        uu_ord[0] += 1;
                        if !exact_interlaces(&u[j], &u[l]) {
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
                        if !exact_interlaces(&d[j], &u[l]) {
                            du_all[1] += 1;
                            if du_all[1] <= 3 {
                                println!("FAIL DU({},{}) {:?}", j, l, board);
                            }
                        }
                    }
                }
            }
            for j in 1..=m {
                for l in 1..j {
                    if !pz(&a[j]) && !pz(&w[l]) {
                        aw_ext[0] += 1;
                        if !exact_interlaces(&a[j], &w[l]) {
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
                    if !pz(&u[j]) && !pz(&u[l]) {
                        uu_rev[0] += 1;
                        if !exact_interlaces(&u[j], &u[l]) {
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
            for j in 1..=m {
                for l in j + 1..=m {
                    if !pz(&d[j]) && !pz(&d[l]) {
                        dd_fwd[0] += 1;
                        if !exact_interlaces(&d[j], &d[l]) {
                            dd_fwd[1] += 1;
                        }
                        dd_rev[0] += 1;
                        if !exact_interlaces(&d[l], &d[j]) {
                            dd_rev[1] += 1;
                        }
                    }
                }
            }
            for j in 1..=m {
                for l in 1..=m {
                    if !pz(&s[j]) && !pz(&t_a[l]) {
                        st_all[0] += 1;
                        if !exact_interlaces(&s[j], &t_a[l]) {
                            st_all[1] += 1;
                            if st_all[1] <= 3 {
                                println!("FAIL S≪T({},{}) {:?}", j, l, board);
                            }
                        }
                        if j > l {
                            st_gap[0] += 1;
                            if !exact_interlaces(&s[j], &t_a[l]) {
                                st_gap[1] += 1;
                            }
                        }
                    }
                }
            }
        }
        println!("n={}: {}/{} valid", n, valid, total);
    }

    println!("\n=== RESULTS ===\n");
    let show = |n: &str, c: [u64; 2]| {
        if c[0] == 0 {
            println!("  {}: (none)", n);
        } else {
            println!(
                "  {}: {}/{} {}",
                n,
                c[0] - c[1],
                c[0],
                if c[1] == 0 { "✓" } else { "✗" }
            );
        }
    };
    println!("Theorem:");
    show("(a) A_j≪W_l j≤l", aw_ord);
    show("(b) W_j≪W_l j≤l", ww_ord);
    show("(c) D_j≪U_l all", du_all);
    show("(d) U_j≪U_l j≤l", uu_ord);
    println!("\nStrengthenings:");
    show("A_j≪W_l j>l", aw_ext);
    show("D_j≪D_l fwd", dd_fwd);
    show("D_l≪D_j rev", dd_rev);
    show("U_j≪U_l j>l rev", uu_rev);
    println!("\nλ+ level:");
    show("S_j≪T_l all", st_all);
    show("S_j≪T_l j>l GAP", st_gap);
}
