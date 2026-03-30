/// Generalized shift lemma: investigate weaker sufficient conditions for
/// L = f + (t-1)h to be real-rooted when h does NOT weakly interlace f.
///
/// For descent sets S with 1 in S (the failure cases), compute f_p, h_p
/// for each valid position p where the shift lemma fails, and look for
/// structural relationships between f and h that could serve as alternative
/// sufficient conditions.
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{
    check_weak_interlacing, format_poly, is_real_rooted, real_roots,
};
use polynomial_tools::{Polynomial, FieldRing};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::ToPrimitive;
use std::collections::BTreeMap;

type Q = Ratio<BigInt>;

fn build_poly(perms: &[&Vec<u8>]) -> Vec<i64> {
    if perms.is_empty() {
        return vec![0];
    }
    let max_s = perms.iter().map(|s| compute(s, Stat::Swaps)).max().unwrap_or(0);
    let mut coeffs = vec![0i64; max_s + 1];
    for s in perms {
        coeffs[compute(s, Stat::Swaps)] += 1;
    }
    trim(&mut coeffs);
    coeffs
}

fn trim(p: &mut Vec<i64>) {
    while p.len() > 1 && *p.last().unwrap() == 0 {
        p.pop();
    }
}

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut c = vec![0i64; len];
    for (i, &v) in a.iter().enumerate() { c[i] += v; }
    for (i, &v) in b.iter().enumerate() { c[i] += v; }
    trim(&mut c);
    c
}

fn poly_mul_t(a: &[i64]) -> Vec<i64> {
    if a.is_empty() || (a.len() == 1 && a[0] == 0) {
        return vec![0];
    }
    let mut c = vec![0i64; a.len() + 1];
    for (i, &v) in a.iter().enumerate() { c[i + 1] = v; }
    c
}

fn poly_mul_t_minus_1(a: &[i64]) -> Vec<i64> {
    let ta = poly_mul_t(a);
    let len = ta.len().max(a.len());
    let mut c = vec![0i64; len];
    for (i, &v) in ta.iter().enumerate() { c[i] += v; }
    for (i, &v) in a.iter().enumerate() { c[i] -= v; }
    trim(&mut c);
    c
}

fn weakly_interlace(f: &[i64], g: &[i64]) -> bool {
    if f.len() <= 1 || g.len() <= 1 {
        return true;
    }
    let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
    check_weak_interlacing(small, large) == Some(true)
}

fn descent_set_str(s: u64, n: u8) -> String {
    let positions: Vec<String> = (0..n - 1)
        .filter(|&i| s & (1 << i) != 0)
        .map(|i| (i + 1).to_string())
        .collect();
    if positions.is_empty() { "{}".to_string() }
    else { format!("{{{}}}", positions.join(",")) }
}

fn epsilon1(pi: &[u8], p: usize, n: u8) -> usize {
    if p <= 1 || p >= n as usize { return 0; }
    if pi[p - 2] + 1 == pi[p - 1] { 1 } else { 0 }
}

fn epsilon2(q: usize, p: usize) -> usize {
    if q + 2 <= p { 1 } else { 0 }
}

fn roots_as_floats(coeffs: &[i64]) -> String {
    if let Some(roots) = real_roots(coeffs) {
        let mut sorted: Vec<f64> = roots.iter().map(|r| {
            let num = r.numer().to_f64().unwrap_or(f64::NAN);
            let den = r.denom().to_f64().unwrap_or(1.0);
            num / den
        }).collect();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        format!("{:?}", sorted.iter().map(|x| format!("{:.4}", x)).collect::<Vec<_>>())
    } else {
        "NOT REAL-ROOTED".to_string()
    }
}

fn poly_degree(p: &[i64]) -> Option<usize> {
    p.iter().rposition(|&c| c != 0)
}

/// Compute GCD of two polynomials over Q using polynomial-tools Polynomial type.
fn poly_gcd_i64(f: &[i64], g: &[i64]) -> Vec<i64> {
    let fq: Polynomial<Q> = Polynomial::new(
        f.iter().map(|&c| Q::from_integer(BigInt::from(c))).collect()
    );
    let gq: Polynomial<Q> = Polynomial::new(
        g.iter().map(|&c| Q::from_integer(BigInt::from(c))).collect()
    );
    let d = fq.gcd(&gq);
    // Convert back to i64 by clearing denominators
    let coeffs = d.coeffs();
    if coeffs.is_empty() {
        return vec![0];
    }
    // Find LCM of denominators
    let mut lcm = BigInt::from(1);
    for c in coeffs {
        let den = c.denom().clone();
        lcm = num::integer::lcm(lcm, den);
    }
    let lcm_q = Q::from_integer(lcm);
    coeffs.iter().map(|c| {
        let scaled = c * &lcm_q;
        scaled.to_integer().to_i64().unwrap_or(0)
    }).collect()
}

/// Try to divide f by g exactly. Returns Some(quotient) if exact, None otherwise.
fn try_exact_div(f: &[i64], g: &[i64]) -> Option<Vec<i64>> {
    let fq: Polynomial<Q> = Polynomial::new(
        f.iter().map(|&c| Q::from_integer(BigInt::from(c))).collect()
    );
    let gq: Polynomial<Q> = Polynomial::new(
        g.iter().map(|&c| Q::from_integer(BigInt::from(c))).collect()
    );
    if gq.is_zero() {
        return None;
    }
    let (quot, rem) = fq.div_rem(&gq);
    if rem.is_zero() {
        let coeffs = quot.coeffs();
        // Check all rational coefficients are actually integers
        let mut result = Vec::new();
        for c in coeffs {
            if *c.denom() != BigInt::from(1) {
                // Not integer coefficients; clear denominators
                let mut lcm = BigInt::from(1);
                for cc in coeffs {
                    lcm = num::integer::lcm(lcm, cc.denom().clone());
                }
                let lcm_q = Q::from_integer(lcm);
                return Some(coeffs.iter().map(|cc| {
                    let scaled = cc * &lcm_q;
                    scaled.to_integer().to_i64().unwrap_or(0)
                }).collect());
            }
            result.push(c.to_integer().to_i64().unwrap_or(0));
        }
        Some(result)
    } else {
        None
    }
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    println!("================================================================");
    println!("  Generalized shift lemma: analyzing failure cases (1 in S)");
    println!("================================================================\n");

    let mut total_failures = 0u64;
    let mut l_always_rr = 0u64;
    let mut f_divides_h = 0u64;
    let mut h_divides_f = 0u64;
    let mut nontrivial_gcd = 0u64;
    let mut ratio_rr = 0u64;
    let mut deg_h_le_half_deg_f = 0u64;

    for n in 4..=max_n {
        let all_nm1 = all_permutations(n - 1);
        let all_n = all_permutations(n);

        // Group target by descent set
        let mut target_by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for s in &all_n {
            target_by_des.entry(descent_set_bitmask(s)).or_default().push(s);
        }

        let mut n_failures = 0u32;
        let mut n_l_rr = 0u32;

        for (&target_ds, _target_perms) in &target_by_des {
            // Only examine descent sets with 1 in S
            if target_ds & 1 == 0 {
                continue;
            }

            let valid_pos: Vec<usize> = (1..=n as usize)
                .filter(|&p| {
                    let ok_descent = p >= n as usize || (target_ds & (1 << (p - 1))) != 0;
                    let ok_ascent = p <= 1 || (target_ds & (1 << (p - 2))) == 0;
                    ok_descent && ok_ascent
                })
                .collect();

            if valid_pos.len() < 2 {
                continue;
            }

            for &p in &valid_pos {
                let mut f_p = vec![0i64; 1];
                let mut h_p = vec![0i64; 1];

                for pi in &all_nm1 {
                    let mut sigma = Vec::with_capacity(n as usize);
                    for i in 0..p - 1 { sigma.push(pi[i]); }
                    sigma.push(n);
                    for i in p - 1..pi.len() { sigma.push(pi[i]); }

                    if descent_set_bitmask(&sigma) != target_ds {
                        continue;
                    }

                    let sw = compute(pi, Stat::Swaps);
                    let q = pi.iter().position(|&v| v == n - 1).unwrap() + 1;
                    let e1 = epsilon1(pi, p, n);
                    let e2 = epsilon2(q, p);

                    let f_deg = e2 + sw;
                    if f_deg >= f_p.len() { f_p.resize(f_deg + 1, 0); }
                    f_p[f_deg] += 1;

                    if e1 == 1 {
                        let h_deg = e2 + sw;
                        if h_deg >= h_p.len() { h_p.resize(h_deg + 1, 0); }
                        h_p[h_deg] += 1;
                    }
                }

                trim(&mut f_p);
                trim(&mut h_p);

                let h_int_f = weakly_interlace(&h_p, &f_p);
                if h_int_f {
                    continue; // shift lemma works here, skip
                }

                // This is a FAILURE case for the shift lemma
                total_failures += 1;
                n_failures += 1;

                let correction = poly_mul_t_minus_1(&h_p);
                let l_p = poly_add(&f_p, &correction);
                let l_rr = l_p.len() <= 2 || is_real_rooted(&l_p);
                if l_rr { l_always_rr += 1; n_l_rr += 1; }

                let deg_f = poly_degree(&f_p);
                let deg_h = poly_degree(&h_p);

                // Check: does h divide f?
                let h_div_f = if deg_h.is_some() && deg_h.unwrap() > 0 {
                    try_exact_div(&f_p, &h_p).is_some()
                } else {
                    false
                };
                if h_div_f { h_divides_f += 1; }

                // Check: does f divide h?
                let f_div_h = if deg_f.is_some() && deg_f.unwrap() > 0 {
                    try_exact_div(&h_p, &f_p).is_some()
                } else {
                    false
                };
                if f_div_h { f_divides_h += 1; }

                // Compute gcd(f, h)
                let gcd = poly_gcd_i64(&f_p, &h_p);
                let gcd_deg = poly_degree(&gcd).unwrap_or(0);
                let has_nontrivial_gcd = gcd_deg > 0;
                if has_nontrivial_gcd { nontrivial_gcd += 1; }

                // Check if f/h is a polynomial with all real roots
                // (only if h divides f)
                let ratio_is_rr = if h_div_f {
                    let q = try_exact_div(&f_p, &h_p).unwrap();
                    q.len() <= 2 || is_real_rooted(&q)
                } else {
                    false
                };
                if ratio_is_rr { ratio_rr += 1; }

                // Check degree relationship
                let deg_h_val = deg_h.unwrap_or(0);
                let deg_f_val = deg_f.unwrap_or(0);
                let deg_check = deg_f_val > 0 && deg_h_val <= deg_f_val / 2;
                if deg_check { deg_h_le_half_deg_f += 1; }

                // Print details for first few failures per n
                if n_failures <= 8 {
                    println!("--- n={} Des={} p={} ---", n, descent_set_str(target_ds, n), p);
                    println!("  f  = {}  (deg {})", format_poly(&f_p), deg_f.unwrap_or(0));
                    println!("  h  = {}  (deg {})", format_poly(&h_p), deg_h.unwrap_or(0));
                    println!("  L  = {}  real-rooted: {}", format_poly(&l_p), l_rr);
                    println!("  h interlaces f: {}", h_int_f);
                    println!("  h divides f: {}", h_div_f);
                    if h_div_f {
                        let q = try_exact_div(&f_p, &h_p).unwrap();
                        println!("  f/h = {}  real-rooted: {}", format_poly(&q),
                            q.len() <= 2 || is_real_rooted(&q));
                    }
                    println!("  f divides h: {}", f_div_h);
                    println!("  gcd(f,h) = {}  (deg {})", format_poly(&gcd), gcd_deg);
                    if has_nontrivial_gcd && !h_div_f {
                        let f_red = try_exact_div(&f_p, &gcd);
                        let h_red = try_exact_div(&h_p, &gcd);
                        if let (Some(fr), Some(hr)) = (f_red, h_red) {
                            println!("  f/gcd = {}", format_poly(&fr));
                            println!("  h/gcd = {}", format_poly(&hr));
                            let red_interlace = weakly_interlace(&hr, &fr);
                            println!("  (h/gcd) interlaces (f/gcd): {}", red_interlace);
                        }
                    }
                    println!("  deg(h) <= deg(f)/2: {}", deg_check);
                    println!("  roots of f: {}", roots_as_floats(&f_p));
                    println!("  roots of h: {}", roots_as_floats(&h_p));
                    println!("  roots of L: {}", roots_as_floats(&l_p));

                    // Also check: after removing common factor of t^k, do they interlace?
                    let f_min_pow = f_p.iter().position(|&c| c != 0).unwrap_or(0);
                    let h_min_pow = h_p.iter().position(|&c| c != 0).unwrap_or(0);
                    let common_t_pow = f_min_pow.min(h_min_pow);
                    if common_t_pow > 0 {
                        let f_shifted: Vec<i64> = f_p[common_t_pow..].to_vec();
                        let h_shifted: Vec<i64> = h_p[common_t_pow..].to_vec();
                        let int_after_shift = weakly_interlace(&h_shifted, &f_shifted);
                        println!("  After removing t^{}: h' interlaces f': {}", common_t_pow, int_after_shift);
                    }
                    println!();
                }
            }
        }

        if n_failures > 0 {
            println!("n={}: {} shift lemma failures (1 in S), {} of {} have L real-rooted",
                n, n_failures, n_l_rr, n_failures);
        } else {
            println!("n={}: no shift lemma failures with 1 in S", n);
        }
    }

    println!("\n================================================================");
    println!("  Summary across all n");
    println!("================================================================");
    println!("Total failure cases (h does NOT interlace f): {}", total_failures);
    println!("L = f + (t-1)h is always real-rooted: {} / {}", l_always_rr, total_failures);
    println!();
    println!("Structural checks:");
    println!("  h divides f:                {} / {}", h_divides_f, total_failures);
    println!("  f divides h:                {} / {}", f_divides_h, total_failures);
    println!("  nontrivial gcd(f,h):        {} / {}", nontrivial_gcd, total_failures);
    println!("  f/h is real-rooted (when h|f): {} / {}", ratio_rr, total_failures);
    println!("  deg(h) <= deg(f)/2:         {} / {}", deg_h_le_half_deg_f, total_failures);
}
