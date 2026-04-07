//! Continuous deformation approach: does B(α) ≪ A(α) hold for all α ∈ [0,1]?
//!
//! A(α) = G(c+1, s) + α·t·R,  B(α) = G(c, s) + α·t·R
//! where R = R_{mu'_{<s+1>}}.
//!
//! At α=0: G(c+1,s) ≪ G(c,s) (known by IH on s).
//! At α=1: G(c+1,s+1) ≪ G(c,s+1) (what we want to prove).
//!
//! Key insight: B(α) - A(α) = D is CONSTANT (independent of α)!
//! D = G(c,s) - G(c+1,s) = f_c + (t-1)·f_{c+1}, and D(0) = 0.
//!
//! For interlacing to break, a root of A(α) and B(α) must collide at a root of D.
//! But D(0) = 0, and A(α)(0) = B(α)(0) = 1 > 0 for all α.
//! So the collision must happen at a NON-ZERO root of D.
//!
//! Tests:
//! A. Verify B(α) ≪ A(α) for many α ∈ [0,1]
//! B. Compute roots of D and check they don't interfere
//! C. Check: is D/t always real-rooted? Are roots of D/t outside
//!    the range of roots of A(α) and B(α)?

use combpoly::rook_placements::non_nesting_rook_polynomial;
use polynomial_tools::{check_weak_interlacing, format_poly, is_real_rooted};

fn strip_columns(mu: &[usize], c: usize) -> Vec<usize> {
    mu.iter()
        .map(|&x| x.saturating_sub(c))
        .filter(|&x| x > 0)
        .collect()
}

fn partitions(n: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut buf = Vec::new();
    part_rec(n, n, &mut buf, &mut result);
    result
}

fn part_rec(n: usize, max: usize, buf: &mut Vec<usize>, result: &mut Vec<Vec<usize>>) {
    if n == 0 {
        result.push(buf.clone());
        return;
    }
    for k in (1..=n.min(max)).rev() {
        buf.push(k);
        part_rec(n - k, k, buf, result);
        buf.pop();
    }
}

fn poly_degree(p: &[i64]) -> usize {
    p.iter().rposition(|&c| c != 0).unwrap_or(0)
}

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let n = a.len().max(b.len());
    let mut r = vec![0i64; n];
    for (i, &c) in a.iter().enumerate() { r[i] += c; }
    for (i, &c) in b.iter().enumerate() { r[i] += c; }
    r
}

fn poly_sub(a: &[i64], b: &[i64]) -> Vec<i64> {
    let n = a.len().max(b.len());
    let mut r = vec![0i64; n];
    for (i, &c) in a.iter().enumerate() { r[i] += c; }
    for (i, &c) in b.iter().enumerate() { r[i] -= c; }
    r
}

fn poly_mul_t(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; p.len() + 1];
    for (i, &c) in p.iter().enumerate() { r[i + 1] = c; }
    r
}

fn poly_scale(p: &[i64], s: i64) -> Vec<i64> {
    p.iter().map(|&c| c * s).collect()
}

/// Divide by t (remove constant term which must be 0)
fn poly_div_t(p: &[i64]) -> Vec<i64> {
    assert_eq!(p.first().copied().unwrap_or(0), 0, "poly_div_t: constant term must be 0");
    if p.len() <= 1 { return vec![0]; }
    p[1..].to_vec()
}

fn interlaces(f: &[i64], g: &[i64]) -> Option<bool> {
    let df = poly_degree(f);
    let dg = poly_degree(g);
    if df == 0 && dg == 0 { return Some(true); }
    if df > dg + 1 || dg > df + 1 { return None; }
    if df <= dg {
        check_weak_interlacing(f, g)
    } else {
        check_weak_interlacing(g, f)
    }
}

fn main() {
    let max_n = 14;

    // =========================================================================
    // TEST A: Verify B(α) ≪ A(α) along the path α = 0, 0.01, ..., 1
    // =========================================================================
    println!("=== TEST A: Continuous deformation α ∈ [0,1] ===\n");

    // Use scaled version: A(k) = scale * G(c+1,s) + k * t * R, with scale = 100
    // B(k) = scale * G(c,s) + k * t * R
    // Then α = k/scale.
    let scale = 100i64;
    let mut a_pass = 0;
    let mut a_total = 0;
    let mut a_fail: Vec<String> = Vec::new();

    for n in 2..=max_n {
        for mu in partitions(n) {
            if mu.len() < 2 { continue; }
            let ell = mu.len();
            let m = mu[ell - 1];
            let mp = &mu[..ell - 1];

            for c in 0..m {
                for s in (c + 1)..m {
                    // Build G(c+1, s) and G(c, s)
                    let mut g_c1 = non_nesting_rook_polynomial(&strip_columns(mp, c + 1));
                    let mut g_c = non_nesting_rook_polynomial(&strip_columns(mp, c));
                    for j in (c + 2)..=s {
                        let t_fj = poly_mul_t(&non_nesting_rook_polynomial(&strip_columns(mp, j)));
                        g_c1 = poly_add(&g_c1, &t_fj);
                    }
                    for j in (c + 1)..=s {
                        let t_fj = poly_mul_t(&non_nesting_rook_polynomial(&strip_columns(mp, j)));
                        g_c = poly_add(&g_c, &t_fj);
                    }

                    let r = non_nesting_rook_polynomial(&strip_columns(mp, s + 1));
                    let tr = poly_mul_t(&r);

                    // Test at α = k/scale for k = 0, 1, ..., scale
                    for k in 0..=scale {
                        // A(k) = scale * g_c1 + k * tr
                        // B(k) = scale * g_c + k * tr
                        let a_k = poly_add(&poly_scale(&g_c1, scale), &poly_scale(&tr, k));
                        let b_k = poly_add(&poly_scale(&g_c, scale), &poly_scale(&tr, k));

                        a_total += 1;
                        let res = interlaces(&a_k, &b_k);
                        if res == Some(true) {
                            a_pass += 1;
                        } else {
                            if a_fail.len() < 10 {
                                a_fail.push(format!(
                                    "mu={:?} c={} s={} α={}/{}: {:?}",
                                    mu, c, s, k, scale, res,
                                ));
                            }
                        }
                    }
                }
            }
        }
    }

    println!("Deformation check: {a_pass}/{a_total}");
    if !a_fail.is_empty() {
        println!("\nFAILURES:");
        for ex in &a_fail {
            println!("  {ex}");
        }
    } else {
        println!("Interlacing maintained along entire path! ✓");
    }

    // =========================================================================
    // TEST B: Compute D = G(c,s) - G(c+1,s) and check D/t properties
    // =========================================================================
    println!("\n=== TEST B: Properties of D = G(c,s) - G(c+1,s) ===\n");

    let mut d_rr_pass = 0;
    let mut d_rr_total = 0;
    let mut d_fail: Vec<String> = Vec::new();

    for n in 2..=max_n {
        for mu in partitions(n) {
            if mu.len() < 2 { continue; }
            let ell = mu.len();
            let m = mu[ell - 1];
            let mp = &mu[..ell - 1];

            for c in 0..m {
                let f_c = non_nesting_rook_polynomial(&strip_columns(mp, c));
                let f_c1 = non_nesting_rook_polynomial(&strip_columns(mp, c + 1));

                // D = f_c + (t-1)·f_{c+1} = f_c - f_{c+1} + t·f_{c+1}
                let d_poly = poly_add(&poly_sub(&f_c, &f_c1), &poly_mul_t(&f_c1));

                // Verify D(0) = 0
                assert_eq!(d_poly.first().copied().unwrap_or(0), 0,
                    "D(0) should be 0 for mu={:?} c={}", mu, c);

                // D/t
                let d_div_t = poly_div_t(&d_poly);

                d_rr_total += 1;
                if is_real_rooted(&d_div_t) {
                    d_rr_pass += 1;
                } else {
                    if d_fail.len() < 10 {
                        d_fail.push(format!(
                            "mu={:?} c={}: D/t = {} NOT real-rooted",
                            mu, c, format_poly(&d_div_t),
                        ));
                    }
                }
            }
        }
    }

    println!("D/t real-rooted: {d_rr_pass}/{d_rr_total}");
    if !d_fail.is_empty() {
        println!("\nFAILURES:");
        for ex in &d_fail {
            println!("  {ex}");
        }
    } else {
        println!("D/t always real-rooted! ✓");
    }

    // =========================================================================
    // TEST C: D/t nonneg coefficients? (Would mean D/t has all roots ≤ 0)
    // =========================================================================
    println!("\n=== TEST C: D/t nonneg coefficients ===\n");

    let mut c_pass = 0;
    let mut c_total = 0;
    let mut c_fail: Vec<String> = Vec::new();

    for n in 2..=max_n {
        for mu in partitions(n) {
            if mu.len() < 2 { continue; }
            let ell = mu.len();
            let m = mu[ell - 1];
            let mp = &mu[..ell - 1];

            for c in 0..m {
                let f_c = non_nesting_rook_polynomial(&strip_columns(mp, c));
                let f_c1 = non_nesting_rook_polynomial(&strip_columns(mp, c + 1));

                let d_poly = poly_add(&poly_sub(&f_c, &f_c1), &poly_mul_t(&f_c1));
                let d_div_t = poly_div_t(&d_poly);

                c_total += 1;
                if d_div_t.iter().all(|&v| v >= 0) {
                    c_pass += 1;
                } else {
                    if c_fail.len() < 10 {
                        c_fail.push(format!(
                            "mu={:?} c={}: D/t = {}",
                            mu, c, format_poly(&d_div_t),
                        ));
                    }
                }
            }
        }
    }

    println!("D/t nonneg: {c_pass}/{c_total}");
    if !c_fail.is_empty() {
        println!("\nNEGATIVE COEFFICIENTS:");
        for ex in &c_fail {
            println!("  {ex}");
        }
    } else {
        println!("D/t always has nonneg coefficients! ✓");
    }

    // =========================================================================
    // TEST D: Print D/t for small examples to see structure
    // =========================================================================
    println!("\n=== TEST D: D/t for small examples ===\n");

    let test_shapes = vec![
        vec![2, 1], vec![3, 2, 1], vec![3, 2], vec![3, 3],
        vec![4, 3, 2, 1], vec![3, 3, 3], vec![4, 2],
    ];

    for mu in &test_shapes {
        let ell = mu.len();
        if ell < 2 { continue; }
        let m = mu[ell - 1];
        let mp = &mu[..ell - 1];

        println!("mu = {:?}, mu' = {:?}, m = {}:", mu, mp, m);
        for c in 0..m {
            let f_c = non_nesting_rook_polynomial(&strip_columns(mp, c));
            let f_c1 = non_nesting_rook_polynomial(&strip_columns(mp, c + 1));

            let d_poly = poly_add(&poly_sub(&f_c, &f_c1), &poly_mul_t(&f_c1));
            let d_div_t = poly_div_t(&d_poly);

            println!("  c={}: f_c={}, f_{{c+1}}={}, D/t={}",
                c, format_poly(&f_c), format_poly(&f_c1), format_poly(&d_div_t));

            // Also show: D/t = (f_c - f_{c+1})/t + f_{c+1}
            let diff = poly_sub(&f_c, &f_c1);
            let diff_div_t = if diff.first().copied().unwrap_or(0) == 0 && diff.len() > 1 {
                Some(poly_div_t(&diff))
            } else {
                None
            };

            if let Some(ddt) = &diff_div_t {
                println!("       (f_c - f_{{c+1}})/t = {}, so D/t = {} + {}",
                    format_poly(ddt), format_poly(ddt), format_poly(&f_c1));
            }
        }
        println!();
    }

    // =========================================================================
    // TEST E: Does D/t ≪ A(α) and D/t ≪ B(α) for all α?
    // If so, roots of D = t·(D/t) can never collide with roots of A(α) or B(α)
    // because the roots of D/t are "outside" those of A(α) and B(α).
    // =========================================================================
    println!("=== TEST E: Does D/t interlace G(c+1,s) and G(c,s)? ===\n");

    let mut e_pass = 0;
    let mut e_total = 0;
    let mut e_fail: Vec<String> = Vec::new();

    for n in 2..=max_n {
        for mu in partitions(n) {
            if mu.len() < 2 { continue; }
            let ell = mu.len();
            let m = mu[ell - 1];
            let mp = &mu[..ell - 1];

            for c in 0..m {
                let f_c = non_nesting_rook_polynomial(&strip_columns(mp, c));
                let f_c1 = non_nesting_rook_polynomial(&strip_columns(mp, c + 1));
                let d_poly = poly_add(&poly_sub(&f_c, &f_c1), &poly_mul_t(&f_c1));
                let d_div_t = poly_div_t(&d_poly);

                // Check D/t ≪ G(c+1, s) for s = c+1, ..., m
                let mut g_c1 = f_c1.clone();
                for s in (c + 1)..=m {
                    if s > c + 1 {
                        let t_fs = poly_mul_t(&non_nesting_rook_polynomial(&strip_columns(mp, s)));
                        g_c1 = poly_add(&g_c1, &t_fs);
                    }

                    e_total += 1;
                    let res = interlaces(&d_div_t, &g_c1);
                    if res == Some(true) {
                        e_pass += 1;
                    } else {
                        if e_fail.len() < 15 {
                            e_fail.push(format!(
                                "mu={:?} c={} s={}: D/t ≪ G({},{}) = {:?} | {} ≪ {}",
                                mu, c, s, c + 1, s, res,
                                format_poly(&d_div_t), format_poly(&g_c1),
                            ));
                        }
                    }
                }
            }
        }
    }

    println!("D/t ≪ G(c+1,s): {e_pass}/{e_total}");
    if !e_fail.is_empty() {
        println!("\nFAILURES:");
        for ex in &e_fail {
            println!("  {ex}");
        }
    } else {
        println!("D/t always interlaces G(c+1,s)! ✓");
    }

    // =========================================================================
    // TEST F: The KEY claim — D/t ≪ G(c,s) for all s?
    // =========================================================================
    println!("\n=== TEST F: D/t ≪ G(c,s) for all s ===\n");

    let mut f_pass = 0;
    let mut f_total = 0;
    let mut f_fail: Vec<String> = Vec::new();

    for n in 2..=max_n {
        for mu in partitions(n) {
            if mu.len() < 2 { continue; }
            let ell = mu.len();
            let m = mu[ell - 1];
            let mp = &mu[..ell - 1];

            for c in 0..m {
                let f_c = non_nesting_rook_polynomial(&strip_columns(mp, c));
                let f_c1 = non_nesting_rook_polynomial(&strip_columns(mp, c + 1));
                let d_poly = poly_add(&poly_sub(&f_c, &f_c1), &poly_mul_t(&f_c1));
                let d_div_t = poly_div_t(&d_poly);

                let mut g_c = non_nesting_rook_polynomial(&strip_columns(mp, c));
                for s in (c + 1)..=m {
                    let t_fs = poly_mul_t(&non_nesting_rook_polynomial(&strip_columns(mp, s)));
                    g_c = poly_add(&g_c, &t_fs);

                    f_total += 1;
                    let res = interlaces(&d_div_t, &g_c);
                    if res == Some(true) {
                        f_pass += 1;
                    } else {
                        if f_fail.len() < 15 {
                            f_fail.push(format!(
                                "mu={:?} c={} s={}: D/t ≪ G({},{}) = {:?} | {} ≪ {}",
                                mu, c, s, c, s, res,
                                format_poly(&d_div_t), format_poly(&g_c),
                            ));
                        }
                    }
                }
            }
        }
    }

    println!("D/t ≪ G(c,s): {f_pass}/{f_total}");
    if !f_fail.is_empty() {
        println!("\nFAILURES:");
        for ex in &f_fail {
            println!("  {ex}");
        }
    } else {
        println!("D/t always interlaces G(c,s)! ✓");
    }
}
