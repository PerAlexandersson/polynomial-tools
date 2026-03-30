/// Detailed investigation of the interlacing relationship Sa(t) << L^{(p_b)}(t).
///
/// Sa(t) = sum_{q >= p_a-1} f~_q^{(p_a)}(t) is the unshifted source polynomial for position p_a.
/// L^{(p)}(t) = sum_q t^{eps2(p,q)} f~_q^{(p)}(t) is the output polynomial at position p.
///
/// Key decomposition:
///   L^{(p_a)} = t * La + Sa  where La = sum_{q < p_a-1} f~_q^{(p_a)}(t) (the "shifted" source)
///
/// Tests:
///   1. Degree relationship: deg(L^{(p_b)}) vs deg(Sa)
///   2. Differences L^{(p_b)} - t*Sa and L^{(p_b)} - Sa: real-rootedness?
///   3. Does L^{(p_b)} = t*Sa + something real-rooted?
///   4. Residuals from zone decomposition
///   5. Full interlacing chain: L^{(p_a)} << L^{(p_b)}
///   6. KEY: Sa << L^{(p_a)}  (since L^{(p_a)} = Sa + t*La)
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{
    check_weak_interlacing, format_poly, is_real_rooted, real_roots,
};
use std::collections::BTreeMap;

// ============================================================
// Helpers
// ============================================================

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

fn coeff(poly: &[i64], k: usize) -> i64 {
    if k < poly.len() {
        poly[k]
    } else {
        0
    }
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

fn poly_sub(f: &[i64], g: &[i64]) -> Vec<i64> {
    let len = f.len().max(g.len());
    let mut r = vec![0i64; len];
    for (i, &c) in f.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in g.iter().enumerate() {
        r[i] -= c;
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

fn poly_deg(f: &[i64]) -> Option<usize> {
    for i in (0..f.len()).rev() {
        if f[i] != 0 {
            return Some(i);
        }
    }
    None
}

fn is_zero_poly(f: &[i64]) -> bool {
    f.iter().all(|&c| c == 0)
}

fn to_f64(p: &[i64]) -> Vec<f64> {
    p.iter().map(|&c| c as f64).collect()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum Zone {
    L,
    M,
    R,
}

fn classify_zone(q: u8, p_a: u8, p_b: u8) -> Zone {
    if p_a >= 2 && q <= p_a - 2 {
        Zone::L
    } else if p_b >= 2 && q <= p_b - 2 {
        Zone::M
    } else {
        Zone::R
    }
}

/// Build the source distribution: map from (modified_swaps, q) -> count
/// for a given target descent set s_mask, insertion position p, and size n.
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

/// Build the raw zone polynomial: sum of N(s,q) for q in the given zone, as polynomial in s.
fn build_zone_raw(
    dist: &BTreeMap<(usize, u8), usize>,
    p_a: u8,
    p_b: u8,
    zone: Zone,
) -> Vec<i64> {
    let mut poly = vec![0i64; 1];
    for (&(s, q), &c) in dist {
        if classify_zone(q, p_a, p_b) == zone {
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

/// Build the f~_q polynomial for a given q: sum of N(s,q) as polynomial in s.
fn build_fq_poly(dist: &BTreeMap<(usize, u8), usize>, q_val: u8) -> Vec<i64> {
    let mut poly = vec![0i64; 1];
    for (&(s, q), &c) in dist {
        if q == q_val {
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

/// Compute L^{(p)} = sum_q t^{eps2(p,q)} f~_q(t) from source distribution
fn build_output_from_source(
    dist: &BTreeMap<(usize, u8), usize>,
    p: u8,
) -> Vec<i64> {
    let mut result = vec![0i64; 1];
    for (&(s, q), &c) in dist {
        let e2 = eps2(p, q);
        let idx = s + e2;
        if idx >= result.len() {
            result.resize(idx + 1, 0);
        }
        result[idx] += c as i64;
    }
    while result.len() > 1 && *result.last().unwrap() == 0 {
        result.pop();
    }
    result
}

/// Compute Sa = sum_{q >= p-1} f~_q^{(p)}(t) (unshifted source)
/// These are all q where eps2(p, q) = 0, i.e., q >= p-1.
fn build_sa(dist: &BTreeMap<(usize, u8), usize>, p: u8) -> Vec<i64> {
    let mut poly = vec![0i64; 1];
    for (&(s, q), &c) in dist {
        if eps2(p, q) == 0 {
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

/// Compute La = sum_{q < p-1} f~_q^{(p)}(t) (shifted source)
/// These are all q where eps2(p, q) = 1, i.e., q <= p-2.
fn build_la(dist: &BTreeMap<(usize, u8), usize>, p: u8) -> Vec<i64> {
    let mut poly = vec![0i64; 1];
    for (&(s, q), &c) in dist {
        if eps2(p, q) == 1 {
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
        .unwrap_or(7);

    println!("=== Sa(t) << L^{{(p_b)}}(t) detailed analysis ===");
    println!("=== KEY TEST: Sa << L^{{(p_a)}} ===\n");

    // Counters for all tests
    let mut total_consec = 0u64;
    let mut deg_lpb_eq_sa_plus1 = 0u64;
    let mut deg_lpb_eq_sa = 0u64;
    let mut deg_other = 0u64;

    // Test 2: differences
    let mut lpb_minus_tsa_rr = 0u64;
    let mut lpb_minus_tsa_not_rr = 0u64;
    let mut lpb_minus_sa_rr = 0u64;
    let mut lpb_minus_sa_not_rr = 0u64;

    // Test 5: L^{(p_a)} << L^{(p_b)}
    let mut lpa_interlace_lpb = 0u64;
    let mut lpa_interlace_lpb_fail = 0u64;

    // Test 6: Sa << L^{(p_a)}
    let mut sa_interlace_lpa = 0u64;
    let mut sa_interlace_lpa_fail = 0u64;
    let mut sa_interlace_lpa_trivial = 0u64; // Sa = 0 or La = 0

    // Additional: La << Sa
    let mut la_interlace_sa = 0u64;
    let mut la_interlace_sa_fail = 0u64;

    // Additional: Sa << L^{(p_b)} (the main target)
    let mut sa_interlace_lpb = 0u64;
    let mut sa_interlace_lpb_fail = 0u64;

    // Additional: t*La << Sa
    let mut tla_interlace_sa = 0u64;
    let mut tla_interlace_sa_fail = 0u64;

    // Test: L^{(p_a)} = Sa + t*La identity verification
    let mut identity_ok = 0u64;
    let mut identity_fail = 0u64;

    // Track failures for detailed reporting
    let mut sa_lpa_failures: Vec<String> = Vec::new();
    let mut sa_lpb_failures: Vec<String> = Vec::new();

    for n in 4..=max_n {
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

        let perms_n = all_permutations(n);
        let mut by_descent_n: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_n {
            by_descent_n
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        for s_mask in 0u64..(1 << (n - 1)) {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            // Pre-compute source distributions for all valid positions
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
                        .filter(|pi| {
                            pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p
                        })
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

                let lpa = match lp_polys.get(&p_a) {
                    Some(p) => p.clone(),
                    None => continue,
                };
                let lpb = match lp_polys.get(&p_b) {
                    Some(p) => p.clone(),
                    None => continue,
                };

                // Compute Sa and La for source A
                let sa = build_sa(dist_a, p_a);
                let la = build_la(dist_a, p_a);

                // Verify identity: L^{(p_a)} = Sa + t*La
                let tla = poly_mul_t(&la);
                let sa_plus_tla = poly_add(&sa, &tla);
                if sa_plus_tla == lpa {
                    identity_ok += 1;
                } else {
                    identity_fail += 1;
                    // Also try building output from source
                    let lpa_from_src = build_output_from_source(dist_a, p_a);
                    if lpa_from_src != lpa {
                        println!(
                            "  IDENTITY MISMATCH: n={} S={} pa={} pb={}",
                            n,
                            descent_set_to_string(s_mask, n),
                            p_a,
                            p_b
                        );
                        println!("    L^(pa) = {}", format_poly(&lpa));
                        println!("    Sa+tLa = {}", format_poly(&sa_plus_tla));
                        println!("    from_src = {}", format_poly(&lpa_from_src));
                    }
                }

                total_consec += 1;

                // === Test 1: Degree relationship ===
                let deg_sa = poly_deg(&sa);
                let deg_lpb = poly_deg(&lpb);
                match (deg_sa, deg_lpb) {
                    (Some(ds), Some(dl)) => {
                        if dl == ds + 1 {
                            deg_lpb_eq_sa_plus1 += 1;
                        } else if dl == ds {
                            deg_lpb_eq_sa += 1;
                        } else {
                            deg_other += 1;
                        }
                    }
                    _ => {
                        deg_other += 1;
                    }
                }

                // === Test 2: Differences ===
                if !is_zero_poly(&sa) {
                    let tsa = poly_mul_t(&sa);
                    let diff1 = poly_sub(&lpb, &tsa);
                    if !is_zero_poly(&diff1) {
                        if diff1.len() <= 2 || is_real_rooted(&diff1) {
                            lpb_minus_tsa_rr += 1;
                        } else {
                            lpb_minus_tsa_not_rr += 1;
                        }
                    } else {
                        lpb_minus_tsa_rr += 1;
                    }

                    let diff2 = poly_sub(&lpb, &sa);
                    if !is_zero_poly(&diff2) {
                        if diff2.len() <= 2 || is_real_rooted(&diff2) {
                            lpb_minus_sa_rr += 1;
                        } else {
                            lpb_minus_sa_not_rr += 1;
                        }
                    } else {
                        lpb_minus_sa_rr += 1;
                    }
                }

                // === Test 5: L^{(p_a)} << L^{(p_b)} ===
                if !is_zero_poly(&lpa) && !is_zero_poly(&lpb) {
                    match check_weak_interlacing(&lpa, &lpb) {
                        Some(true) => {
                            lpa_interlace_lpb += 1;
                        }
                        _ => {
                            lpa_interlace_lpb_fail += 1;
                            println!(
                                "  FAIL L^(pa) << L^(pb): n={} S={} pa={} pb={}",
                                n,
                                descent_set_to_string(s_mask, n),
                                p_a,
                                p_b
                            );
                            println!("    L^(pa) = {}", format_poly(&lpa));
                            println!("    L^(pb) = {}", format_poly(&lpb));
                        }
                    }
                }

                // === Test 6: Sa << L^{(p_a)} === KEY TEST ===
                if is_zero_poly(&sa) || is_zero_poly(&la) {
                    // If La = 0, then L^{(p_a)} = Sa, so Sa << L^{(p_a)} trivially
                    // (same-degree interlacing = self-interlacing, which is true for real-rooted)
                    // If Sa = 0, then 0 << L^{(p_a)} trivially
                    sa_interlace_lpa_trivial += 1;
                    sa_interlace_lpa += 1;
                } else if !is_zero_poly(&lpa) {
                    match check_weak_interlacing(&sa, &lpa) {
                        Some(true) => {
                            sa_interlace_lpa += 1;
                        }
                        _ => {
                            sa_interlace_lpa_fail += 1;
                            let msg = format!(
                                "  FAIL Sa << L^(pa): n={} S={} pa={} pb={}\n    Sa     = {}\n    L^(pa) = {}\n    La     = {}\n    t*La   = {}",
                                n,
                                descent_set_to_string(s_mask, n),
                                p_a,
                                p_b,
                                format_poly(&sa),
                                format_poly(&lpa),
                                format_poly(&la),
                                format_poly(&tla),
                            );
                            println!("{}", msg);
                            sa_lpa_failures.push(msg);
                        }
                    }
                }

                // === Additional: La << Sa ===
                if !is_zero_poly(&la) && !is_zero_poly(&sa) {
                    match check_weak_interlacing(&la, &sa) {
                        Some(true) => {
                            la_interlace_sa += 1;
                        }
                        _ => {
                            la_interlace_sa_fail += 1;
                        }
                    }
                }

                // === Additional: t*La << Sa ===
                if !is_zero_poly(&la) && !is_zero_poly(&sa) {
                    match check_weak_interlacing(&tla, &sa) {
                        Some(true) => {
                            tla_interlace_sa += 1;
                        }
                        _ => {
                            tla_interlace_sa_fail += 1;
                        }
                    }
                }

                // === Sa << L^{(p_b)} === THE MAIN TARGET ===
                if !is_zero_poly(&sa) && !is_zero_poly(&lpb) {
                    match check_weak_interlacing(&sa, &lpb) {
                        Some(true) => {
                            sa_interlace_lpb += 1;
                        }
                        _ => {
                            sa_interlace_lpb_fail += 1;
                            let msg = format!(
                                "  FAIL Sa << L^(pb): n={} S={} pa={} pb={}\n    Sa     = {}\n    L^(pb) = {}",
                                n,
                                descent_set_to_string(s_mask, n),
                                p_a,
                                p_b,
                                format_poly(&sa),
                                format_poly(&lpb),
                            );
                            if sa_lpb_failures.len() < 10 {
                                println!("{}", msg);
                            }
                            sa_lpb_failures.push(msg);
                        }
                    }
                }
            }
        }

        println!("  n={}: {} consecutive pairs processed", n, total_consec);
    }

    println!();
    println!("{}", "=".repeat(80));
    println!("SUMMARY (n=4..{})", max_n);
    println!("{}", "=".repeat(80));

    println!(
        "\nIdentity L^(pa) = Sa + t*La:  OK={}/{}  {}",
        identity_ok,
        identity_ok + identity_fail,
        if identity_fail == 0 {
            "ALWAYS HOLDS"
        } else {
            "HAS FAILURES"
        }
    );

    println!("\n--- Degree relationship: deg(L^(pb)) vs deg(Sa) ---");
    println!(
        "  deg(L^(pb)) = deg(Sa)+1:  {}  ({:.1}%)",
        deg_lpb_eq_sa_plus1,
        100.0 * deg_lpb_eq_sa_plus1 as f64 / total_consec as f64
    );
    println!(
        "  deg(L^(pb)) = deg(Sa):    {}  ({:.1}%)",
        deg_lpb_eq_sa,
        100.0 * deg_lpb_eq_sa as f64 / total_consec as f64
    );
    println!("  Other:                    {}", deg_other);

    println!("\n--- Differences ---");
    println!(
        "  L^(pb) - t*Sa real-rooted: {}/{}",
        lpb_minus_tsa_rr,
        lpb_minus_tsa_rr + lpb_minus_tsa_not_rr
    );
    println!(
        "  L^(pb) - Sa real-rooted:   {}/{}",
        lpb_minus_sa_rr,
        lpb_minus_sa_rr + lpb_minus_sa_not_rr
    );

    println!("\n--- Interlacing tests ---");
    let total5 = lpa_interlace_lpb + lpa_interlace_lpb_fail;
    println!(
        "  L^(pa) << L^(pb):  {}/{}  {}",
        lpa_interlace_lpb,
        total5,
        if lpa_interlace_lpb_fail == 0 {
            "ALWAYS"
        } else {
            "FAILS"
        }
    );

    let total6 = sa_interlace_lpa + sa_interlace_lpa_fail;
    println!(
        "\n  *** Sa << L^(pa):  {}/{}  {} ***  (trivial: {})",
        sa_interlace_lpa,
        total6,
        if sa_interlace_lpa_fail == 0 {
            "ALWAYS"
        } else {
            "FAILS"
        },
        sa_interlace_lpa_trivial,
    );

    let total_la_sa = la_interlace_sa + la_interlace_sa_fail;
    println!(
        "  La << Sa:          {}/{}  {}",
        la_interlace_sa,
        total_la_sa,
        if la_interlace_sa_fail == 0 {
            "ALWAYS"
        } else {
            "FAILS"
        }
    );

    let total_tla_sa = tla_interlace_sa + tla_interlace_sa_fail;
    println!(
        "  t*La << Sa:        {}/{}  {}",
        tla_interlace_sa,
        total_tla_sa,
        if tla_interlace_sa_fail == 0 {
            "ALWAYS"
        } else {
            "FAILS"
        }
    );

    let total_main = sa_interlace_lpb + sa_interlace_lpb_fail;
    println!(
        "\n  *** Sa << L^(pb):  {}/{}  {} ***  (TARGET)",
        sa_interlace_lpb,
        total_main,
        if sa_interlace_lpb_fail == 0 {
            "ALWAYS"
        } else {
            "FAILS"
        }
    );

    // Proof strategy summary
    println!("\n{}", "=".repeat(80));
    println!("PROOF STRATEGY ANALYSIS");
    println!("{}", "=".repeat(80));

    if sa_interlace_lpa_fail == 0 && lpa_interlace_lpb_fail == 0 {
        println!("\nBOTH Sa << L^(pa) and L^(pa) << L^(pb) hold universally!");
        println!("By transitivity of interlacing: Sa << L^(pb) follows!");
        println!();
        println!("Proof path:");
        println!("  1. Sa << L^(pa)  [from: L^(pa) = Sa + t*La with La << Sa]");
        println!("  2. L^(pa) << L^(pb)  [this is the interlacing chain at level n]");
        println!("  3. Sa << L^(pb)  [by transitivity of (1) and (2)]");
        println!();
        println!("For step 1: L^(pa) = Sa + t*La, and if La << Sa (from source chain),");
        println!("  then by Lemma 2.4(iii) or similar, Sa << Sa + t*La = L^(pa).");
    } else {
        println!("\nSome test failed. Detailed failures above.");
    }

    if !sa_lpa_failures.is_empty() {
        println!("\n--- Sa << L^(pa) FAILURES ({}) ---", sa_lpa_failures.len());
        for f in sa_lpa_failures.iter().take(5) {
            println!("{}", f);
        }
    }

    if !sa_lpb_failures.is_empty() {
        println!("\n--- Sa << L^(pb) FAILURES ({}) ---", sa_lpb_failures.len());
        for f in sa_lpb_failures.iter().take(5) {
            println!("{}", f);
        }
    }

    // Extra detail: print examples for small n
    println!("\n{}", "=".repeat(80));
    println!("DETAILED EXAMPLES (n=5)");
    println!("{}", "=".repeat(80));
    {
        let n: u8 = 5;
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

        let perms_n = all_permutations(n);
        let mut by_descent_n: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_n {
            by_descent_n
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        let mut example_count = 0;
        for s_mask in 0u64..(1 << (n - 1)) {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }
            if example_count >= 4 {
                break;
            }

            let mut source_dists: BTreeMap<u8, BTreeMap<(usize, u8), usize>> = BTreeMap::new();
            for &p in &vp {
                source_dists.insert(
                    p,
                    build_source_dist(s_mask, p, n, &by_descent_prev, &perm_prev_list),
                );
            }

            let mut lp_polys: BTreeMap<u8, Vec<i64>> = BTreeMap::new();
            if let Some(class) = by_descent_n.get(&s_mask) {
                for &p in &vp {
                    let vals: Vec<usize> = class
                        .iter()
                        .filter(|pi| {
                            pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p
                        })
                        .map(|pi| compute(pi, Stat::Swaps))
                        .collect();
                    if !vals.is_empty() {
                        lp_polys.insert(p, build_poly(&vals));
                    }
                }
            }

            println!(
                "\nS = {} positions = {:?}",
                descent_set_to_string(s_mask, n),
                vp
            );

            for idx in 0..vp.len() - 1 {
                let p_a = vp[idx];
                let p_b = vp[idx + 1];

                let dist_a = match source_dists.get(&p_a) {
                    Some(d) => d,
                    None => continue,
                };

                let lpa = match lp_polys.get(&p_a) {
                    Some(p) => p.clone(),
                    None => continue,
                };
                let lpb = match lp_polys.get(&p_b) {
                    Some(p) => p.clone(),
                    None => continue,
                };

                let sa = build_sa(dist_a, p_a);
                let la = build_la(dist_a, p_a);
                let tla = poly_mul_t(&la);

                println!("  (pa={}, pb={})", p_a, p_b);
                println!("    La     = {}", format_poly(&la));
                println!("    Sa     = {}", format_poly(&sa));
                println!("    t*La   = {}", format_poly(&tla));
                println!("    L^(pa) = {}", format_poly(&lpa));
                println!("    L^(pb) = {}", format_poly(&lpb));

                let sa_plus_tla = poly_add(&sa, &tla);
                println!(
                    "    Sa+tLa = {} (== L^(pa)? {})",
                    format_poly(&sa_plus_tla),
                    sa_plus_tla == lpa
                );

                if !is_zero_poly(&la) && !is_zero_poly(&sa) {
                    let la_ll_sa = check_weak_interlacing(&la, &sa) == Some(true);
                    let tla_ll_sa = check_weak_interlacing(&tla, &sa) == Some(true);
                    let sa_ll_lpa = check_weak_interlacing(&sa, &lpa) == Some(true);
                    println!(
                        "    La<<Sa: {}  t*La<<Sa: {}  Sa<<L^(pa): {}",
                        la_ll_sa, tla_ll_sa, sa_ll_lpa
                    );
                }

                let sa_ll_lpb = check_weak_interlacing(&sa, &lpb) == Some(true);
                let lpa_ll_lpb = check_weak_interlacing(&lpa, &lpb) == Some(true);
                println!(
                    "    L^(pa)<<L^(pb): {}  Sa<<L^(pb): {}",
                    lpa_ll_lpb, sa_ll_lpb
                );

                // Show roots
                if !is_zero_poly(&sa) {
                    if let Some(roots) = real_roots(&sa) {
                        let root_strs: Vec<String> = roots.iter().map(|r| format!("{}", r)).collect();
                        println!("    roots(Sa) = [{}]", root_strs.join(", "));
                    }
                }
                if let Some(roots) = real_roots(&lpa) {
                    let root_strs: Vec<String> = roots.iter().map(|r| format!("{}", r)).collect();
                    println!("    roots(L^(pa)) = [{}]", root_strs.join(", "));
                }
                if let Some(roots) = real_roots(&lpb) {
                    let root_strs: Vec<String> = roots.iter().map(|r| format!("{}", r)).collect();
                    println!("    roots(L^(pb)) = [{}]", root_strs.join(", "));
                }
            }
            example_count += 1;
        }
    }
}
