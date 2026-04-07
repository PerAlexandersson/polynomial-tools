//! Deeper exploration of q-deformation and refinement approaches.
//!
//! Key findings from qdeform:
//! - R_mu(q, t) is real-rooted in t for all q in [0,1]
//! - Column-strip adjacent interlacing holds for all q
//! - Refinement pieces nearly pairwise compatible (245/246)
//!
//! New tests:
//! A. Reverse-order partial sums (adding pieces from high c to low)
//! B. Brändén input pairwise for q=1 (standard rook)
//! C. Does the matrix N step give pairwise interlacing for q=1?
//! D. Alternative q-weight: q^{crossings} where crossing = pair with
//!    (r1 < r2, c1 > c2), so q=0 gives crossing-free = standard rook order,
//!    q=1 gives all placements

use combpoly::rook_placements::{non_nesting_rook_polynomial, rook_polynomial};
use polynomial_tools::{check_weak_interlacing, format_poly, is_real_rooted};

type Placement = Vec<(usize, usize)>;

fn all_placements(mu: &[usize]) -> Vec<Placement> {
    let mut result = Vec::new();
    let mut current = Vec::new();
    let max_col = mu.iter().copied().max().unwrap_or(0);
    let mut used_cols = vec![false; max_col];
    gen_placements(mu, 0, &mut current, &mut used_cols, &mut result);
    result
}

fn gen_placements(
    mu: &[usize],
    row: usize,
    current: &mut Placement,
    used_cols: &mut Vec<bool>,
    result: &mut Vec<Placement>,
) {
    if row == mu.len() {
        result.push(current.clone());
        return;
    }
    gen_placements(mu, row + 1, current, used_cols, result);
    for c in 0..mu[row] {
        if !used_cols[c] {
            used_cols[c] = true;
            current.push((row, c));
            gen_placements(mu, row + 1, current, used_cols, result);
            current.pop();
            used_cols[c] = false;
        }
    }
}

fn count_nestings(placement: &Placement) -> usize {
    let mut count = 0;
    for i in 0..placement.len() {
        for j in i + 1..placement.len() {
            let (r1, c1) = placement[i];
            let (r2, c2) = placement[j];
            if r1 < r2 && c1 < c2 {
                count += 1;
            }
        }
    }
    count
}

fn count_crossings(placement: &Placement) -> usize {
    let mut count = 0;
    for i in 0..placement.len() {
        for j in i + 1..placement.len() {
            let (r1, c1) = placement[i];
            let (r2, c2) = placement[j];
            if r1 < r2 && c1 > c2 {
                count += 1;
            }
        }
    }
    count
}

/// q-rook polynomial with q^{nest} weight
fn q_rook_nest(mu: &[usize]) -> Vec<Vec<i64>> {
    q_rook_generic(mu, count_nestings)
}

/// q-rook polynomial with q^{crossing} weight
fn q_rook_cross(mu: &[usize]) -> Vec<Vec<i64>> {
    q_rook_generic(mu, count_crossings)
}

fn q_rook_generic(mu: &[usize], stat: fn(&Placement) -> usize) -> Vec<Vec<i64>> {
    let placements = all_placements(mu);
    let max_rooks = mu.len().min(mu.iter().copied().max().unwrap_or(0));
    let max_stat = max_rooks * (max_rooks.saturating_sub(1)) / 2;
    let mut coeffs = vec![vec![0i64; max_stat + 1]; max_rooks + 1];
    for p in &placements {
        let k = p.len();
        let s = stat(p);
        if s < coeffs[k].len() {
            coeffs[k][s] += 1;
        } else {
            // Extend
            coeffs[k].resize(s + 1, 0);
            coeffs[k][s] += 1;
        }
    }
    for row in coeffs.iter_mut() {
        while row.len() > 1 && *row.last().unwrap() == 0 {
            row.pop();
        }
    }
    coeffs
}

fn specialize_q_rational(qpoly: &[Vec<i64>], num: i64, den: i64) -> Vec<i64> {
    let max_q_deg = qpoly.iter().map(|r| r.len().saturating_sub(1)).max().unwrap_or(0);
    let mut num_pow = vec![1i128; max_q_deg + 1];
    let mut den_pow = vec![1i128; max_q_deg + 1];
    for i in 1..=max_q_deg {
        num_pow[i] = num_pow[i - 1] * num as i128;
        den_pow[i] = den_pow[i - 1] * den as i128;
    }
    let common = den_pow[max_q_deg];
    let mut result = Vec::new();
    for row in qpoly {
        let mut val = 0i128;
        for (j, &c) in row.iter().enumerate() {
            val += c as i128 * num_pow[j] * den_pow[max_q_deg - j];
        }
        result.push(val);
    }
    let g = result
        .iter()
        .fold(common.unsigned_abs(), |acc, &x| gcd128(acc, x.unsigned_abs()));
    if g > 1 {
        for v in result.iter_mut() {
            *v /= g as i128;
        }
    }
    result.into_iter().map(|x| x as i64).collect()
}

fn gcd128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

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

fn poly_mul_t(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; p.len() + 1];
    for (i, &c) in p.iter().enumerate() { r[i + 1] = c; }
    r
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

fn check_compatible_pair(f: &[i64], g: &[i64]) -> bool {
    let test_alphas: Vec<(i64, i64)> = vec![
        (1, 0), (0, 1), (1, 1), (1, 2), (2, 1), (1, 3), (3, 1),
        (1, 5), (5, 1), (1, 10), (10, 1), (2, 3), (3, 2),
        (7, 11), (11, 7), (1, 100), (100, 1),
    ];
    for &(a, b) in &test_alphas {
        let n = f.len().max(g.len());
        let mut combo = vec![0i64; n];
        for (i, &c) in f.iter().enumerate() { combo[i] += a * c; }
        for (i, &c) in g.iter().enumerate() { combo[i] += b * c; }
        if !is_real_rooted(&combo) {
            return false;
        }
    }
    true
}

fn main() {
    let max_n = 10;

    // =========================================================================
    // TEST A: Reverse-order partial sums (adding pieces from high c to low)
    // S_m = piece[m], S_{m-1} = piece[m] + piece[m-1], ..., S_0 = R_mu
    // =========================================================================
    println!("=== TEST A: Reverse-order partial sums ===\n");

    let mut ra_pass = 0;
    let mut ra_total = 0;
    let mut ra_fail: Vec<String> = Vec::new();

    for n in 2..=max_n {
        for mu in partitions(n) {
            if mu.len() < 2 { continue; }
            let ell = mu.len();
            let last = mu[ell - 1];
            let mp = &mu[..ell - 1];

            // Build pieces: piece[0] = R_{mu'}, piece[c] = t · R_{mu'_{<c>}}
            let mut pieces = vec![non_nesting_rook_polynomial(mp)];
            for c in 1..=last {
                let sub = strip_columns(mp, c);
                pieces.push(poly_mul_t(&non_nesting_rook_polynomial(&sub)));
            }

            // Reverse partial sums: start from piece[last], add down to piece[0]
            let mut sums = Vec::new();
            let mut acc = vec![0i64; 1];
            for c in (0..=last).rev() {
                acc = poly_add(&acc, &pieces[c]);
                sums.push(acc.clone());
            }
            // sums[0] = piece[last]
            // sums[1] = piece[last] + piece[last-1]
            // ...
            // sums[last] = R_mu

            for i in 0..sums.len() - 1 {
                ra_total += 1;
                let di = poly_degree(&sums[i]);
                let di1 = poly_degree(&sums[i + 1]);
                let res = if di <= di1 {
                    check_weak_interlacing(&sums[i], &sums[i + 1])
                } else {
                    check_weak_interlacing(&sums[i + 1], &sums[i])
                };
                if res == Some(true) {
                    ra_pass += 1;
                } else {
                    if ra_fail.len() < 10 {
                        ra_fail.push(format!(
                            "mu={:?}: S[{}] ≪ S[{}] = {:?} | {} ≪ {}",
                            mu, i, i + 1, res,
                            format_poly(&sums[i]), format_poly(&sums[i + 1]),
                        ));
                    }
                }
            }
        }
    }

    println!("Reverse partial sums adjacent interlacing: {ra_pass}/{ra_total}");
    if !ra_fail.is_empty() {
        println!("\nFAILURES:");
        for ex in &ra_fail {
            println!("  {ex}");
        }
    } else {
        println!("All pass! ✓");
    }

    // =========================================================================
    // TEST B: Brändén input pairwise for q=1 (standard rook)
    // =========================================================================
    println!("\n=== TEST B: Brändén input pairwise for q=1 ===\n");

    let mut bp_pass = 0;
    let mut bp_total = 0;
    let mut bp_fail: Vec<String> = Vec::new();

    for n in 2..=max_n {
        for mu in partitions(n) {
            if mu.len() < 2 { continue; }
            let ell = mu.len();
            let m = mu[ell - 1];
            let mp = &mu[..ell - 1];

            let f: Vec<Vec<i64>> = (0..=m)
                .rev()
                .map(|c| rook_polynomial(&strip_columns(mp, c)))
                .collect();

            for i in 0..f.len() {
                for j in i + 1..f.len() {
                    bp_total += 1;
                    if check_compatible_pair(&f[i], &f[j]) {
                        bp_pass += 1;
                    } else {
                        if bp_fail.len() < 10 {
                            bp_fail.push(format!(
                                "mu={:?}: f[{}],f[{}] | {} and {}",
                                mu, i, j,
                                format_poly(&f[i]), format_poly(&f[j]),
                            ));
                        }
                    }
                }
            }
        }
    }

    println!("Brändén input pairwise (q=1): {bp_pass}/{bp_total}");
    if !bp_fail.is_empty() {
        println!("\nINCOMPATIBLE:");
        for ex in &bp_fail {
            println!("  {ex}");
        }
    } else {
        println!("All compatible! ✓");
    }

    // =========================================================================
    // TEST C: Matrix N output pairwise for q=1
    // g_k = f_k + t·(f_1 + ... + f_{k-1}) using standard rook polys
    // =========================================================================
    println!("\n=== TEST C: Matrix N output pairwise for q=1 ===\n");

    let mut mn_pass = 0;
    let mut mn_total = 0;
    let mut mn_fail: Vec<String> = Vec::new();

    for n in 2..=max_n {
        for mu in partitions(n) {
            if mu.len() < 2 { continue; }
            let ell = mu.len();
            let m = mu[ell - 1];
            let mp = &mu[..ell - 1];

            let f: Vec<Vec<i64>> = (0..=m)
                .rev()
                .map(|c| rook_polynomial(&strip_columns(mp, c)))
                .collect();

            let mut g: Vec<Vec<i64>> = Vec::new();
            let mut prefix_sum = vec![0i64; 1];
            for k in 0..f.len() {
                let gk = poly_add(&f[k], &poly_mul_t(&prefix_sum));
                g.push(gk);
                prefix_sum = poly_add(&prefix_sum, &f[k]);
            }

            for i in 0..g.len() {
                for j in i + 1..g.len() {
                    mn_total += 1;
                    if check_compatible_pair(&g[i], &g[j]) {
                        mn_pass += 1;
                    } else {
                        if mn_fail.len() < 10 {
                            mn_fail.push(format!(
                                "mu={:?}: g[{}],g[{}] | {} and {}",
                                mu, i, j,
                                format_poly(&g[i]), format_poly(&g[j]),
                            ));
                        }
                    }
                }
            }
        }
    }

    println!("Matrix N output pairwise (q=1): {mn_pass}/{mn_total}");
    if !mn_fail.is_empty() {
        println!("\nINCOMPATIBLE:");
        for ex in &mn_fail {
            println!("  {ex}");
        }
    } else {
        println!("All compatible! ✓");
    }

    // =========================================================================
    // TEST D: q^{crossing} weight — properties
    // =========================================================================
    println!("\n=== TEST D: q^{{crossing}} weight ===\n");

    println!("q^{{crossing}} polynomials for small shapes:");
    let test_shapes: Vec<Vec<usize>> = vec![
        vec![2, 1], vec![3, 2, 1], vec![2, 2], vec![3, 3],
    ];
    for mu in &test_shapes {
        let qpoly = q_rook_cross(mu);
        println!("  mu={:?}:", mu);
        for (k, row) in qpoly.iter().enumerate() {
            if row.iter().all(|&c| c == 0) { continue; }
            let mut terms = Vec::new();
            for (j, &c) in row.iter().enumerate() {
                if c == 0 { continue; }
                match (c, j) {
                    (c, 0) => terms.push(format!("{c}")),
                    (1, 1) => terms.push("q".to_string()),
                    (c, 1) => terms.push(format!("{c}q")),
                    (1, j) => terms.push(format!("q^{j}")),
                    (c, j) => terms.push(format!("{c}q^{j}")),
                }
            }
            println!("    t^{k}: {}", terms.join(" + "));
        }
    }

    // Check nest + cross = inv for all placements
    println!("\n  Verifying nest + cross = C(k,2) for k-rook placements...");
    let mut inv_check_ok = true;
    for n in 1..=8 {
        for mu in partitions(n) {
            let placements = all_placements(&mu);
            for p in &placements {
                let k = p.len();
                let expected = k * (k.saturating_sub(1)) / 2;
                let actual = count_nestings(p) + count_crossings(p);
                if actual != expected {
                    println!("    FAIL: mu={:?}, placement={:?}, nest+cross={}, expected={}",
                        mu, p, actual, expected);
                    inv_check_ok = false;
                }
            }
        }
    }
    if inv_check_ok {
        println!("  nest + cross = C(k,2) ✓ (so q^cross = q^{{C(k,2)-nest}})\n");
    }

    // Real-rootedness for q^{crossing} weight
    let q_values: Vec<(i64, i64, &str)> = vec![
        (0, 1, "0"), (1, 4, "1/4"), (1, 2, "1/2"), (3, 4, "3/4"), (1, 1, "1"),
    ];

    println!("  Real-rootedness in t for q^{{crossing}} weight:");
    for &(p, d, label) in &q_values {
        let mut rr_pass = 0;
        let mut rr_total = 0;
        for n in 1..=max_n {
            for mu in partitions(n) {
                let qpoly = q_rook_cross(&mu);
                let specialized = specialize_q_rational(&qpoly, p, d);
                rr_total += 1;
                if is_real_rooted(&specialized) {
                    rr_pass += 1;
                }
            }
        }
        println!("    q={label}: {rr_pass}/{rr_total} real-rooted");
    }

    // =========================================================================
    // TEST E: Nesting+crossing bivariate — R(q, p, t) = sum q^nest p^cross t^k
    // Then R(q, 1, t) = our q-deformation.
    // Restrict to q = p for a one-parameter family?
    // Actually, since nest + cross = C(k,2), R(q,p,t) = R(q/p, 1, p·t) up to scaling?
    // No, that doesn't work because nest+cross depends on k.
    //
    // But: R(q, q, t) = sum q^{C(k,2)} * (#k-rook placements) * t^k
    //    = sum C(k,2)_q * c_k t^k  where c_k = total k-rook placements
    //
    // Actually R(q, q, t) = sum_{placements} q^{nest+cross} t^k
    //                     = sum_{k} (#{k-rook pl.}) q^{C(k,2)} t^k
    //
    // So at q = p: sum c_k q^{C(k,2)} t^k. This IS interesting!
    // =========================================================================
    println!("\n=== TEST E: R_mu(q, q, t) = sum c_k q^{{C(k,2)}} t^k ===\n");

    println!("  Verifying this identity...");
    let mut id_ok = true;
    for n in 1..=8 {
        for mu in partitions(n) {
            let placements = all_placements(&mu);
            let max_k = mu.len().min(mu.iter().copied().max().unwrap_or(0));
            let mut c_k = vec![0i64; max_k + 1];
            for p in &placements {
                c_k[p.len()] += 1;
            }

            let _nest_poly = q_rook_nest(&mu);
            let _cross_poly = q_rook_cross(&mu);

            // Compute R(q,q,t) = sum_{k,j1,j2} [q^{j1+j2}] nest[k][j1] * ... no.
            // R(q,q,t) = sum_{placement} q^{nest+cross} t^k = sum_k c_k q^{C(k,2)} t^k
            // Let's verify by brute force
            let mut brute: Vec<Vec<i64>> = vec![vec![]; max_k + 1]; // brute[k][j] = coefficient of t^k q^j in R(q,q,t)
            for p in &placements {
                let k = p.len();
                let total_stat = count_nestings(p) + count_crossings(p);
                if total_stat >= brute[k].len() {
                    brute[k].resize(total_stat + 1, 0);
                }
                brute[k][total_stat] += 1;
            }

            // Should be: brute[k] = [0, ..., 0, c_k] at position C(k,2)
            for k in 0..=max_k {
                let binom = k * k.saturating_sub(1) / 2;
                let expected_len = binom + 1;
                let mut expected_row = vec![0i64; expected_len];
                if c_k[k] > 0 {
                    expected_row[binom] = c_k[k];
                }
                let actual = &brute[k];
                let trimmed: Vec<i64> = if actual.is_empty() {
                    vec![0]
                } else {
                    actual.clone()
                };
                let mut exp_trimmed = expected_row.clone();
                while exp_trimmed.len() > 1 && *exp_trimmed.last().unwrap() == 0 {
                    exp_trimmed.pop();
                }
                let mut act_trimmed = trimmed.clone();
                while act_trimmed.len() > 1 && *act_trimmed.last().unwrap() == 0 {
                    act_trimmed.pop();
                }
                if act_trimmed != exp_trimmed {
                    // This should never happen since nest+cross = C(k,2) is constant
                    println!("    UNEXPECTED: mu={:?}, k={}", mu, k);
                    id_ok = false;
                }
            }
        }
    }
    if id_ok {
        println!("  Identity verified ✓");
    }

    // =========================================================================
    // TEST F: Does the G(c+1,s) ≪ G(c,s) step hold for q=1 (standard rook)?
    // G(c, s) = StdRook(mu'_<c>) + t * sum_{j=c+1}^{s} StdRook(mu'_<j>)
    // =========================================================================
    println!("\n=== TEST F: G(c+1,s) ≪ G(c,s) for standard rook (q=1) ===\n");

    let mut gcs_pass = 0;
    let mut gcs_total = 0;
    let mut gcs_fail: Vec<String> = Vec::new();

    for n in 2..=max_n {
        for mu in partitions(n) {
            if mu.len() < 2 { continue; }
            let ell = mu.len();
            let m = mu[ell - 1];
            let mp = &mu[..ell - 1];

            for c in 0..m {
                let base_c = rook_polynomial(&strip_columns(mp, c));
                let base_c1 = rook_polynomial(&strip_columns(mp, c + 1));

                let term_c1 = poly_mul_t(&rook_polynomial(&strip_columns(mp, c + 1)));
                let mut g_c = poly_add(&base_c, &term_c1);
                let mut g_c1 = base_c1.clone();

                // s = c+1 (before adding more)
                gcs_total += 1;
                let res = interlaces(&g_c1, &g_c);
                if res == Some(true) {
                    gcs_pass += 1;
                } else {
                    if gcs_fail.len() < 10 {
                        gcs_fail.push(format!(
                            "mu={:?}: G({},{}) ≪ G({},{}) = {:?}",
                            mu, c + 1, c + 1, c, c + 1, res,
                        ));
                    }
                }

                for s in (c + 2)..=m {
                    let term = poly_mul_t(&rook_polynomial(&strip_columns(mp, s)));
                    g_c = poly_add(&g_c, &term);
                    g_c1 = poly_add(&g_c1, &term);

                    gcs_total += 1;
                    let res = interlaces(&g_c1, &g_c);
                    if res == Some(true) {
                        gcs_pass += 1;
                    } else {
                        if gcs_fail.len() < 10 {
                            gcs_fail.push(format!(
                                "mu={:?}: G({},{}) ≪ G({},{}) = {:?}",
                                mu, c + 1, s, c, s, res,
                            ));
                        }
                    }
                }
            }
        }
    }

    println!("G(c+1,s) ≪ G(c,s) for std rook: {gcs_pass}/{gcs_total}");
    if !gcs_fail.is_empty() {
        println!("\nFAILURES:");
        for ex in &gcs_fail {
            println!("  {ex}");
        }
    } else {
        println!("All pass! ✓");
    }

    // =========================================================================
    // TEST G: For general q, check G(c+1,s) ≪ G(c,s) using q-rook polys
    // =========================================================================
    println!("\n=== TEST G: G(c+1,s) ≪ G(c,s) for various q (nesting weight) ===\n");

    for &(p, d, label) in &q_values {
        let mut g_pass = 0;
        let mut g_total = 0;

        for n in 2..=max_n {
            for mu in partitions(n) {
                if mu.len() < 2 { continue; }
                let ell = mu.len();
                let m = mu[ell - 1];
                let mp = &mu[..ell - 1];

                for c in 0..m {
                    let base_c = specialize_q_rational(&q_rook_nest(&strip_columns(mp, c)), p, d);
                    let base_c1 = specialize_q_rational(&q_rook_nest(&strip_columns(mp, c + 1)), p, d);

                    let term_c1 = poly_mul_t(&specialize_q_rational(&q_rook_nest(&strip_columns(mp, c + 1)), p, d));
                    let mut g_c = poly_add(&base_c, &term_c1);
                    let mut g_c1 = base_c1.clone();

                    g_total += 1;
                    let res = interlaces(&g_c1, &g_c);
                    if res == Some(true) { g_pass += 1; }

                    for s in (c + 2)..=m {
                        let term = poly_mul_t(&specialize_q_rational(&q_rook_nest(&strip_columns(mp, s)), p, d));
                        g_c = poly_add(&g_c, &term);
                        g_c1 = poly_add(&g_c1, &term);

                        g_total += 1;
                        let res = interlaces(&g_c1, &g_c);
                        if res == Some(true) { g_pass += 1; }
                    }
                }
            }
        }
        println!("  q={label}: {g_pass}/{g_total} G(c+1,s) ≪ G(c,s)");
    }
}
