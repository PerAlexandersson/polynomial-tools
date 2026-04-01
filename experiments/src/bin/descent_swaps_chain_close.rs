/// Test the "chain close" interlacing strategy for consecutive valid positions.
///
/// For each consecutive (p_a, p_b) with 1 not in S:
///
///   Sa(t) = unshifted source for p_a (zones M+R relative to (p_a, p_b))
///         = sum_{q >= p_a-1} N_a(s, q)
///
///   Fb(t) = full source for p_b (all zones, unshifted)
///         = Tb(t) + Rb(t) = sum_q N_b(s, q)
///
///   L^{(p_b)} = t*Tb + Rb  (the output polynomial at position p_b)
///     where Tb = sum_{q <= p_b-2} N_b(s,q)  and  Rb = sum_{q >= p_b-1} N_b(s,q)
///
///   Sb(t) = zones M+R of source B = sum_{q >= p_a-1} N_b(s, q)
///
/// Checks:
///   1. Sa << Fb    (source-to-source interlacing across nearby descent sets)
///   2. Fb << L^{(p_b)}  (by Lemma: L^{(p_b)} = Fb + (t-1)*Tb)
///   3. Sa << Sb    (partial source to partial source)
///   4. Sb << Fb    (partial to full source)
///   5. If Sa << Fb << L^{(p_b)}: Sa << L^{(p_b)} by transitivity
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly};
use std::collections::BTreeMap;

// ============================================================
// Helpers
// ============================================================

fn check_interl(f: &[i64], g: &[i64]) -> bool {
    if f.len() <= 1 || g.len() <= 1 {
        return true;
    }
    let (s, l) = if f.len() <= g.len() { (f, g) } else { (g, f) };
    check_weak_interlacing(s, l).unwrap_or(false)
}

fn build_poly(vals: &[usize]) -> Vec<i64> {
    if vals.is_empty() {
        return vec![0];
    }
    let max_s = *vals.iter().max().unwrap();
    let mut coeffs = vec![0i64; max_s + 1];
    for &s in vals {
        coeffs[s] += 1;
    }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
        coeffs.pop();
    }
    coeffs
}

fn valid_positions(s_mask: u64, n: u8) -> Vec<u8> {
    let mut positions = Vec::new();
    for p in 1..n {
        if (s_mask >> (p - 1)) & 1 == 1 && (p < 2 || (s_mask >> (p - 2)) & 1 == 0) {
            positions.push(p);
        }
    }
    if n >= 2 && (s_mask >> (n - 2)) & 1 == 0 {
        positions.push(n);
    }
    positions
}

fn source_asc(s_mask: u64, p: u8, n: u8) -> u64 {
    if n <= 2 {
        return 0;
    }
    if p == n {
        return s_mask;
    }
    let mut sp = 0u64;
    if p == 1 {
        for j in 2..n {
            if (s_mask >> (j - 1)) & 1 == 1 {
                sp |= 1 << (j - 2);
            }
        }
    } else {
        for pos in 1..=(p.saturating_sub(2)) {
            if (s_mask >> (pos - 1)) & 1 == 1 {
                sp |= 1 << (pos - 1);
            }
        }
        for j in (p + 1)..n {
            if (s_mask >> (j - 1)) & 1 == 1 {
                sp |= 1 << (j - 2);
            }
        }
    }
    sp
}

fn source_desc(s_mask: u64, p: u8, n: u8) -> Option<u64> {
    if p <= 1 || p >= n {
        return None;
    }
    Some(source_asc(s_mask, p, n) | (1 << (p - 2)))
}

fn epsilon1(pi: &[u8], p: u8) -> bool {
    let n = pi.len() as u8 + 1;
    if p <= 1 || p >= n {
        return false;
    }
    pi[(p - 2) as usize] + 1 == pi[(p - 1) as usize]
}

fn eps2(p: u8, q: u8) -> usize {
    if p >= 2 && q <= p.saturating_sub(2) {
        1
    } else {
        0
    }
}

fn descent_set_to_string(mask: u64, n: u8) -> String {
    let mut s = String::from("{");
    let mut first = true;
    for i in 1..n {
        if (mask >> (i - 1)) & 1 == 1 {
            if !first {
                s.push(',');
            }
            s.push_str(&i.to_string());
            first = false;
        }
    }
    s.push('}');
    s
}

fn poly_add(f: &[i64], g: &[i64]) -> Vec<i64> {
    let len = f.len().max(g.len());
    let mut r = vec![0i64; len];
    for (i, &c) in f.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in g.iter().enumerate() {
        r[i] += c;
    }
    while r.len() > 1 && *r.last().unwrap() == 0 {
        r.pop();
    }
    r
}

fn poly_mul_t(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; p.len() + 1];
    for (i, &c) in p.iter().enumerate() {
        r[i + 1] = c;
    }
    r
}

fn is_zero_poly(f: &[i64]) -> bool {
    f.iter().all(|&c| c == 0)
}

/// Build source distribution: map from (modified_swaps, q) -> count
/// for a given target descent set s_mask, insertion position p, and size n.
/// Uses BOTH S'_p (ascent source) and S''_p (descent source), with
/// modified swaps = swaps(pi) + eps1(pi, p) for ascent source,
/// and swaps(pi) for descent source.
fn build_source_dist(
    s_mask: u64,
    p: u8,
    n: u8,
    by_descent_prev: &BTreeMap<u64, Vec<usize>>,
    perm_prev_list: &[&Vec<u8>],
) -> BTreeMap<(usize, u8), usize> {
    let sp_asc = source_asc(s_mask, p, n);
    let sp_desc = source_desc(s_mask, p, n);
    let mut dist: BTreeMap<(usize, u8), usize> = BTreeMap::new();
    if let Some(indices) = by_descent_prev.get(&sp_asc) {
        for &idx in indices {
            let pi = perm_prev_list[idx];
            let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
            let sw = compute(pi, Stat::Swaps) + if epsilon1(pi, p) { 1 } else { 0 };
            *dist.entry((sw, q)).or_default() += 1;
        }
    }
    if let Some(sp) = sp_desc {
        if let Some(indices) = by_descent_prev.get(&sp) {
            for &idx in indices {
                let pi = perm_prev_list[idx];
                let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                let sw = compute(pi, Stat::Swaps);
                *dist.entry((sw, q)).or_default() += 1;
            }
        }
    }
    dist
}

/// Sa(t) = unshifted source for position p, relative to (p_a, p_b) zones.
/// = zones M+R of source for p: sum_{q >= p_a-1} N_p(s, q)
/// When p = p_a, this is {q >= p_a - 1}, which are exactly the q with eps2(p_a, q) = 0.
fn build_sa(dist: &BTreeMap<(usize, u8), usize>, p_a: u8) -> Vec<i64> {
    let mut poly = vec![0i64; 1];
    for (&(s, q), &c) in dist {
        // Zones M+R: q >= p_a - 1
        if q + 1 >= p_a {
            // q >= p_a - 1 means q+1 > p_a - 1, i.e. q >= p_a - 1
            // Since q is 1-indexed position, q >= p_a - 1 means eps2(p_a, q) = 0
            if s >= poly.len() {
                poly.resize(s + 1, 0);
            }
            poly[s] += c as i64;
        }
    }
    while poly.len() > 1 && *poly.last().unwrap() == 0 {
        poly.pop();
    }
    poly
}

/// Sb(t) = zones M+R of source B relative to (p_a, p_b).
/// = sum_{q >= p_a-1} N_b(s, q)
fn build_sb_mr(dist: &BTreeMap<(usize, u8), usize>, p_a: u8) -> Vec<i64> {
    build_sa(dist, p_a)
}

/// Full source polynomial for position p: sum over all q of N_p(s, q).
fn build_full_source(dist: &BTreeMap<(usize, u8), usize>) -> Vec<i64> {
    let mut poly = vec![0i64; 1];
    for (&(s, _q), &c) in dist {
        if s >= poly.len() {
            poly.resize(s + 1, 0);
        }
        poly[s] += c as i64;
    }
    while poly.len() > 1 && *poly.last().unwrap() == 0 {
        poly.pop();
    }
    poly
}

/// Tb(t) = sum_{q <= p_b - 2} N_b(s, q) = "shifted" part of source B
/// (these are the q where eps2(p_b, q) = 1)
fn build_tb(dist: &BTreeMap<(usize, u8), usize>, p_b: u8) -> Vec<i64> {
    let mut poly = vec![0i64; 1];
    for (&(s, q), &c) in dist {
        if eps2(p_b, q) == 1 {
            if s >= poly.len() {
                poly.resize(s + 1, 0);
            }
            poly[s] += c as i64;
        }
    }
    while poly.len() > 1 && *poly.last().unwrap() == 0 {
        poly.pop();
    }
    poly
}

/// Rb(t) = sum_{q >= p_b - 1} N_b(s, q) = "unshifted" part of source B
/// (these are the q where eps2(p_b, q) = 0)
fn build_rb(dist: &BTreeMap<(usize, u8), usize>, p_b: u8) -> Vec<i64> {
    let mut poly = vec![0i64; 1];
    for (&(s, q), &c) in dist {
        if eps2(p_b, q) == 0 {
            if s >= poly.len() {
                poly.resize(s + 1, 0);
            }
            poly[s] += c as i64;
        }
    }
    while poly.len() > 1 && *poly.last().unwrap() == 0 {
        poly.pop();
    }
    poly
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    println!("=== Chain close interlacing experiment ===");
    println!("=== Sa << Fb << L^{{(pb)}} => Sa << L^{{(pb)}} ===\n");

    // Counters
    let mut total = 0u64;
    let mut sa_fb_pass = 0u64;
    let mut sa_fb_fail = 0u64;
    let mut fb_lpb_pass = 0u64;
    let mut fb_lpb_fail = 0u64;
    let mut sa_sb_pass = 0u64;
    let mut sa_sb_fail = 0u64;
    let mut sb_fb_pass = 0u64;
    let mut sb_fb_fail = 0u64;
    let mut chain_complete = 0u64; // Sa << Fb AND Fb << L^{(pb)}
    let mut sa_lpb_direct_pass = 0u64;
    let mut sa_lpb_direct_fail = 0u64;
    let mut tb_fb_pass = 0u64;
    let mut tb_fb_fail = 0u64;

    for n in 5..=max_n {
        println!("--- n = {} ---", n);
        let perms_prev = all_permutations(n - 1);
        let mut by_descent_prev: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
        let mut perm_prev_list: Vec<&Vec<u8>> = Vec::new();
        for pi in &perms_prev {
            let idx = perm_prev_list.len();
            perm_prev_list.push(pi);
            by_descent_prev
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(idx);
        }

        // Also compute output polynomials from actual permutations of size n
        let perms_n = all_permutations(n);
        let mut by_descent_n: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_n {
            by_descent_n
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        for s_mask in 0u64..(1 << (n - 1)) {
            // Filter: 1 not in S
            if s_mask & 1 != 0 {
                continue;
            }

            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            // Pre-compute source distributions
            let mut source_dists: BTreeMap<u8, BTreeMap<(usize, u8), usize>> = BTreeMap::new();
            for &p in &vp {
                source_dists.insert(
                    p,
                    build_source_dist(s_mask, p, n, &by_descent_prev, &perm_prev_list),
                );
            }

            // Pre-compute output polynomials L^{(p)} from actual permutations
            let mut lp_polys: BTreeMap<u8, Vec<i64>> = BTreeMap::new();
            if let Some(class) = by_descent_n.get(&s_mask) {
                for &p in &vp {
                    let vals: Vec<usize> = class
                        .iter()
                        .filter(|pi| pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                        .map(|pi| compute(pi, Stat::Swaps))
                        .collect();
                    if !vals.is_empty() {
                        lp_polys.insert(p, build_poly(&vals));
                    }
                }
            }

            // Process consecutive pairs
            for idx in 0..vp.len() - 1 {
                let p_a = vp[idx];
                let p_b = vp[idx + 1];

                let dist_a = match source_dists.get(&p_a) {
                    Some(d) => d,
                    None => continue,
                };
                let dist_b = match source_dists.get(&p_b) {
                    Some(d) => d,
                    None => continue,
                };

                let lpb = match lp_polys.get(&p_b) {
                    Some(p) => p.clone(),
                    None => continue,
                };

                // Sa = zones M+R of source A (q >= p_a - 1)
                let sa = build_sa(dist_a, p_a);
                // Fb = full source for p_b
                let fb = build_full_source(dist_b);
                // Tb = shifted part of source B (q <= p_b - 2)
                let tb = build_tb(dist_b, p_b);
                // Rb = unshifted part of source B (q >= p_b - 1)
                let rb = build_rb(dist_b, p_b);
                // Sb = zones M+R of source B (q >= p_a - 1)
                let sb = build_sb_mr(dist_b, p_a);

                // Verify: Fb = Tb + Rb
                let tb_plus_rb = poly_add(&tb, &rb);
                assert_eq!(
                    fb,
                    tb_plus_rb,
                    "Fb != Tb + Rb for n={} S={} pa={} pb={}",
                    n,
                    descent_set_to_string(s_mask, n),
                    p_a,
                    p_b
                );

                // Verify: L^{(pb)} = t*Tb + Rb
                let t_tb = poly_mul_t(&tb);
                let lpb_check = poly_add(&t_tb, &rb);
                assert_eq!(
                    lpb,
                    lpb_check,
                    "L^(pb) != t*Tb + Rb for n={} S={} pa={} pb={}",
                    n,
                    descent_set_to_string(s_mask, n),
                    p_a,
                    p_b
                );

                if is_zero_poly(&sa) || is_zero_poly(&fb) {
                    continue;
                }

                total += 1;

                // --- Check 1: Sa << Fb ---
                let c1 = check_interl(&sa, &fb);
                if c1 {
                    sa_fb_pass += 1;
                } else {
                    sa_fb_fail += 1;
                    println!(
                        "  FAIL Sa<<Fb: n={} S={} pa={} pb={}",
                        n,
                        descent_set_to_string(s_mask, n),
                        p_a,
                        p_b
                    );
                    println!("    Sa = {}", format_poly(&sa));
                    println!("    Fb = {}", format_poly(&fb));
                }

                // --- Check 2: Fb << L^{(pb)} ---
                let c2 = check_interl(&fb, &lpb);
                if c2 {
                    fb_lpb_pass += 1;
                } else {
                    fb_lpb_fail += 1;
                    println!(
                        "  FAIL Fb<<Lpb: n={} S={} pa={} pb={}",
                        n,
                        descent_set_to_string(s_mask, n),
                        p_a,
                        p_b
                    );
                    println!("    Fb  = {}", format_poly(&fb));
                    println!("    Lpb = {}", format_poly(&lpb));
                }

                // --- Check 3: Sa << Sb ---
                if !is_zero_poly(&sb) {
                    let c3 = check_interl(&sa, &sb);
                    if c3 {
                        sa_sb_pass += 1;
                    } else {
                        sa_sb_fail += 1;
                        println!(
                            "  FAIL Sa<<Sb: n={} S={} pa={} pb={}",
                            n,
                            descent_set_to_string(s_mask, n),
                            p_a,
                            p_b
                        );
                        println!("    Sa = {}", format_poly(&sa));
                        println!("    Sb = {}", format_poly(&sb));
                    }
                }

                // --- Check 4: Sb << Fb ---
                if !is_zero_poly(&sb) {
                    let c4 = check_interl(&sb, &fb);
                    if c4 {
                        sb_fb_pass += 1;
                    } else {
                        sb_fb_fail += 1;
                        println!(
                            "  FAIL Sb<<Fb: n={} S={} pa={} pb={}",
                            n,
                            descent_set_to_string(s_mask, n),
                            p_a,
                            p_b
                        );
                        println!("    Sb = {}", format_poly(&sb));
                        println!("    Fb = {}", format_poly(&fb));
                    }
                }

                // --- Check: Tb << Fb (source chain: shifted part interlaces full source) ---
                if !is_zero_poly(&tb) {
                    let ct = check_interl(&tb, &fb);
                    if ct {
                        tb_fb_pass += 1;
                    } else {
                        tb_fb_fail += 1;
                        println!(
                            "  FAIL Tb<<Fb: n={} S={} pa={} pb={}",
                            n,
                            descent_set_to_string(s_mask, n),
                            p_a,
                            p_b
                        );
                        println!("    Tb = {}", format_poly(&tb));
                        println!("    Fb = {}", format_poly(&fb));
                    }
                }

                // --- Chain completeness ---
                if c1 && c2 {
                    chain_complete += 1;
                }

                // --- Direct check: Sa << L^{(pb)} ---
                let cd = check_interl(&sa, &lpb);
                if cd {
                    sa_lpb_direct_pass += 1;
                } else {
                    sa_lpb_direct_fail += 1;
                    println!(
                        "  FAIL Sa<<Lpb (direct): n={} S={} pa={} pb={}",
                        n,
                        descent_set_to_string(s_mask, n),
                        p_a,
                        p_b
                    );
                    println!("    Sa  = {}", format_poly(&sa));
                    println!("    Lpb = {}", format_poly(&lpb));
                }
            }
        }
    }

    println!("\n=== SUMMARY (n=5..{}) ===", max_n);
    println!(
        "Total consecutive pairs (1 not in S, nontrivial): {}",
        total
    );
    println!();
    println!(
        "(1) Sa << Fb:        pass={} fail={}",
        sa_fb_pass, sa_fb_fail
    );
    println!(
        "(2) Fb << L^(pb):    pass={} fail={}",
        fb_lpb_pass, fb_lpb_fail
    );
    println!(
        "(3) Sa << Sb:        pass={} fail={}",
        sa_sb_pass, sa_sb_fail
    );
    println!(
        "(4) Sb << Fb:        pass={} fail={}",
        sb_fb_pass, sb_fb_fail
    );
    println!(
        "    Tb << Fb:        pass={} fail={}",
        tb_fb_pass, tb_fb_fail
    );
    println!();
    println!(
        "Chain Sa<<Fb<<Lpb:   complete={} / {}",
        chain_complete, total
    );
    println!(
        "Direct Sa<<Lpb:      pass={} fail={}",
        sa_lpb_direct_pass, sa_lpb_direct_fail
    );
}
