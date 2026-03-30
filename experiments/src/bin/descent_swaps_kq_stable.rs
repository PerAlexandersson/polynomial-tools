/// Stability analysis of the G_p(t) structure of Phi_S(x,t).
///
/// Phi_{n,S}(x,t) = sum_{p in P(S)} x^p * G_p(t)
///
/// Each G_p(t) decomposes by the eps2 correction:
///   G_p(t) = B_p(t) + t * A_p(t)
/// where:
///   B_p(t) = sum_{pi: source for p, q >= p-1} t^{s_tilde(pi)}
///   A_p(t) = sum_{pi: source for p, q <= p-2} t^{s_tilde(pi)}
/// (here q = pi^{-1}(n-1), s_tilde = swaps(pi) + eps1(pi,p))
///
/// QUESTIONS:
///   1. Are A_p and B_p individually real-rooted?
///   2. Does B_p interlace A_p (or vice versa)? This would make G_p = B_p + t*A_p stable in t.
///   3. Since Phi = sum_p x^p (B_p + t*A_p), can we write Phi = B(x,t) + t*A(x,t)?
///      Then: Phi stable iff Im(B(x,t)/A(x,t)) >= 0 in the upper half-plane...
///      But Phi is degree > 1 in t, so this isn't directly linear in t.
///   4. What does the FULL bivariate source decomposition look like?
///      f_{p,q}(t) = source poly for insertion pos p, pos(n-1)=q.
///      Phi = sum_{p,q} x^p * t^{eps2(p,q)} * f_{p,q}(t)
///   5. MATRIX STRUCTURE: Think of M(t) as a matrix with rows indexed by q
///      and columns by p, with entry f_{p,q}(t) * t^{eps2(p,q)}.
///      Then Phi = sum_q [sum_p x^p * M_{q,p}(t)] = x^T M(t) 1?
///      No, more like: Phi = sum_p x^p sum_q M_{q,p}(t) = sum_p x^p G_p(t).
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly, is_real_rooted};
use std::collections::BTreeMap;

// ── Complex arithmetic ───────────────────────────────────────────────
fn cmul(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}
fn cadd(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 + b.0, a.1 + b.1)
}
fn cpow(base: (f64, f64), exp: u32) -> (f64, f64) {
    let mut result = (1.0, 0.0);
    for _ in 0..exp {
        result = cmul(result, base);
    }
    result
}
fn cabs(a: (f64, f64)) -> f64 {
    (a.0 * a.0 + a.1 * a.1).sqrt()
}
fn cdiv(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    let denom = b.0 * b.0 + b.1 * b.1;
    (
        (a.0 * b.0 + a.1 * b.1) / denom,
        (a.1 * b.0 - a.0 * b.1) / denom,
    )
}

// ── Simple PRNG ──────────────────────────────────────────────────────
struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed.wrapping_add(1))
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn next_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

// ── Helpers ──────────────────────────────────────────────────────────

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

fn epsilon2(p: u8, q: u8) -> usize {
    if q + 2 <= p {
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

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        r[i] += c;
    }
    while r.len() > 1 && *r.last().unwrap() == 0 {
        r.pop();
    }
    r
}

/// Multiply polynomial by t (shift coefficients by 1)
fn poly_mul_t(a: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; a.len() + 1];
    for (i, &c) in a.iter().enumerate() {
        r[i + 1] = c;
    }
    r
}

fn eval_poly_complex(coeffs: &[i64], z: (f64, f64)) -> (f64, f64) {
    let mut result = (0.0, 0.0);
    for (k, &c) in coeffs.iter().enumerate() {
        if c != 0 {
            result = cadd(result, cmul((c as f64, 0.0), cpow(z, k as u32)));
        }
    }
    result
}

fn eval_phi_direct(class: &[&Vec<u8>], n: u8, x: (f64, f64), t: (f64, f64)) -> (f64, f64) {
    let mut result = (0.0, 0.0);
    for sigma in class {
        let pos = sigma.iter().position(|&v| v == n).unwrap() as u8 + 1;
        let swaps = compute(*sigma, Stat::Swaps);
        result = cadd(result, cmul(cpow(x, pos as u32), cpow(t, swaps as u32)));
    }
    result
}

fn check_stability<F>(eval: F, num_tests: usize, seed: u64) -> (bool, f64)
where
    F: Fn((f64, f64), (f64, f64)) -> (f64, f64),
{
    let mut rng = Rng::new(seed);
    let mut min_abs = f64::MAX;
    for _ in 0..num_tests {
        let a = rng.next_range(-3.0, 3.0);
        let b = rng.next_range(0.01, 3.0);
        let c = rng.next_range(-3.0, 3.0);
        let d = rng.next_range(0.01, 3.0);
        let val = eval((a, b), (c, d));
        let abs_val = cabs(val);
        if abs_val < min_abs {
            min_abs = abs_val;
        }
        if abs_val < 1e-8 {
            return (false, min_abs);
        }
    }
    for _ in 0..num_tests {
        let a = rng.next_range(-5.0, 5.0);
        let b = rng.next_range(0.001, 0.1);
        let c = rng.next_range(-5.0, 5.0);
        let d = rng.next_range(0.001, 0.1);
        let val = eval((a, b), (c, d));
        let abs_val = cabs(val);
        if abs_val < min_abs {
            min_abs = abs_val;
        }
        if abs_val < 1e-8 {
            return (false, min_abs);
        }
    }
    (true, min_abs)
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    let num_tests = 200;

    println!("=== Phi_S(x,t) = sum_p x^p G_p(t) stability decomposition ===");
    println!("=== G_p(t) = B_p(t) + t*A_p(t) where A_p from eps2=1 sources, B_p from eps2=0 ===");
    println!("n = 5..{}\n", max_n);

    // Global counters
    let mut total_cases = 0u32;
    let mut total_phi_stable = 0u32;
    let mut total_ap_rr = 0u32;
    let mut total_bp_rr = 0u32;
    let mut total_ap_bp_tested = 0u32;
    let mut total_bp_int_ap = 0u32; // B_p interlaces A_p
    let mut total_ap_int_bp = 0u32; // A_p interlaces B_p
    let mut total_gp_rr = 0u32;

    // Track f_{p,q} individual stability
    let mut total_fpq_rr = 0u32;
    let mut total_fpq_tested = 0u32;
    let mut total_compat = 0u32;
    let mut total_chain = 0u32;

    for n in 5..=max_n {
        let perms = all_permutations(n);
        let perms_prev = all_permutations(n - 1);

        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        let mut prev_by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        println!("{}", "=".repeat(70));
        println!("  n = {}", n);
        println!("{}", "=".repeat(70));

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            total_cases += 1;
            let s_str = descent_set_to_string(mask, n);

            // Build the full (p, q, s_tilde) data
            // f_{p,q}(t) = generating polynomial of modified swaps for source at (p, q)
            let mut fpq: BTreeMap<(u8, u8), Vec<usize>> = BTreeMap::new();

            for &p in &vp {
                let sp_a = source_asc(mask, p, n);
                let sp_d = source_desc(mask, p, n);

                let source_masks: Vec<u64> = {
                    let mut v = vec![sp_a];
                    if let Some(sd) = sp_d {
                        if sd != sp_a {
                            v.push(sd);
                        }
                    }
                    v
                };

                for &src_mask in &source_masks {
                    if let Some(cls) = prev_by_des.get(&src_mask) {
                        for pi in cls {
                            let e1 = if epsilon1(pi, p) { 1 } else { 0 };
                            let ms = compute(pi, Stat::Swaps) + e1;
                            let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            fpq.entry((p, q)).or_default().push(ms);
                        }
                    }
                }
            }

            // Build A_p(t) and B_p(t) for each p
            let mut ap_polys: BTreeMap<u8, Vec<i64>> = BTreeMap::new();
            let mut bp_polys: BTreeMap<u8, Vec<i64>> = BTreeMap::new();
            let mut gp_polys: BTreeMap<u8, Vec<i64>> = BTreeMap::new();

            for &p in &vp {
                let mut a_vals: Vec<usize> = Vec::new();
                let mut b_vals: Vec<usize> = Vec::new();
                for (&(pp, q), vals) in &fpq {
                    if pp != p {
                        continue;
                    }
                    if epsilon2(p, q) == 1 {
                        a_vals.extend(vals);
                    } else {
                        b_vals.extend(vals);
                    }
                }
                let ap = build_poly(&a_vals);
                let bp = build_poly(&b_vals);
                let gp = poly_add(&bp, &poly_mul_t(&ap));
                ap_polys.insert(p, ap);
                bp_polys.insert(p, bp);
                gp_polys.insert(p, gp);
            }

            // Verify G_p matches direct computation
            let mut gp_direct: BTreeMap<u8, Vec<i64>> = BTreeMap::new();
            for &p in &vp {
                let vals: Vec<usize> = class
                    .iter()
                    .filter(|sigma| sigma.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                    .map(|sigma| compute(*sigma, Stat::Swaps))
                    .collect();
                if !vals.is_empty() {
                    gp_direct.insert(p, build_poly(&vals));
                }
            }

            let gp_match = vp.iter().all(|&p| {
                gp_polys.get(&p).map(|v| v.as_slice()) == gp_direct.get(&p).map(|v| v.as_slice())
            });

            // ── Print ──
            println!(
                "\n  S={}, P(S)={:?}, G_p match: {}",
                s_str, vp, gp_match
            );

            if !gp_match {
                println!("    *** G_p MISMATCH - bug in source computation! ***");
                for &p in &vp {
                    println!("      G_{}(computed) = {}", p, format_poly(gp_polys.get(&p).unwrap_or(&vec![0])));
                    println!("      G_{}(direct)   = {}", p, format_poly(gp_direct.get(&p).unwrap_or(&vec![0])));
                }
                continue;
            }

            for &p in &vp {
                let ap = ap_polys.get(&p).unwrap();
                let bp = bp_polys.get(&p).unwrap();
                let gp = gp_polys.get(&p).unwrap();
                let ap_rr = if ap.len() <= 2 { true } else { is_real_rooted(ap) };
                let bp_rr = if bp.len() <= 2 { true } else { is_real_rooted(bp) };
                let gp_rr = if gp.len() <= 2 { true } else { is_real_rooted(gp) };

                if ap_rr {
                    total_ap_rr += 1;
                }
                if bp_rr {
                    total_bp_rr += 1;
                }
                if gp_rr {
                    total_gp_rr += 1;
                }
                total_ap_bp_tested += 1;

                // Check interlacing: does B_p interlace A_p?
                // B_p interlaces A_p means the roots of B_p separate the roots of A_p.
                // For G_p = B_p + t*A_p to be stable in (variable t with parameter in coeffs),
                // we need B_p to interlace A_p (Hermite-Biehler).
                let bp_int_ap = if bp_rr && ap_rr && bp.len() >= 2 && ap.len() >= 2 {
                    // B_p should interlace A_p if deg(B_p) = deg(A_p) or deg(B_p) = deg(A_p) + 1
                    // Use check_weak_interlacing
                    let (smaller, larger) = if bp.len() <= ap.len() {
                        (bp.as_slice(), ap.as_slice())
                    } else {
                        (ap.as_slice(), bp.as_slice())
                    };
                    match check_weak_interlacing(smaller, larger) {
                        Some(true) => true,
                        _ => false,
                    }
                } else if ap == &vec![0i64] || ap.len() <= 1 {
                    // A_p is constant or zero; G_p = B_p is trivially "stable in t" direction
                    true
                } else {
                    false
                };

                if bp_int_ap {
                    total_bp_int_ap += 1;
                }

                // Also check A_p interlaces B_p
                let ap_int_bp = if bp_rr && ap_rr && bp.len() >= 2 && ap.len() >= 2 {
                    let (smaller, larger) = if ap.len() <= bp.len() {
                        (ap.as_slice(), bp.as_slice())
                    } else {
                        (bp.as_slice(), ap.as_slice())
                    };
                    match check_weak_interlacing(smaller, larger) {
                        Some(true) => true,
                        _ => false,
                    }
                } else {
                    false
                };
                if ap_int_bp {
                    total_ap_int_bp += 1;
                }

                println!(
                    "    p={}: A_p={} [rr:{}]  B_p={} [rr:{}]  G_p={} [rr:{}]  B<=A:{} A<=B:{}",
                    p,
                    format_poly(ap),
                    ap_rr,
                    format_poly(bp),
                    bp_rr,
                    format_poly(gp),
                    gp_rr,
                    bp_int_ap,
                    ap_int_bp
                );
            }

            // Check f_{p,q} individual real-rootedness
            let mut all_fpq_rr = true;
            let mut fpq_details: Vec<String> = Vec::new();
            for (&(p, q), vals) in &fpq {
                total_fpq_tested += 1;
                let poly = build_poly(vals);
                let rr = if poly.len() <= 2 {
                    true
                } else {
                    is_real_rooted(&poly)
                };
                if rr {
                    total_fpq_rr += 1;
                } else {
                    all_fpq_rr = false;
                    fpq_details.push(format!("f_({},{})={} NOT rr", p, q, format_poly(&poly)));
                }
            }
            if !all_fpq_rr {
                for d in &fpq_details {
                    println!("    {}", d);
                }
            } else {
                println!("    All f_{{p,q}} real-rooted: YES");
            }

            // Check COMPATIBLE FAMILY: are G_{p_1}, G_{p_2} compatible for all pairs?
            // Compatible means alpha*G_{p_1} + beta*G_{p_2} is real-rooted for all alpha, beta > 0.
            // Equivalently, G_{p_1} and G_{p_2} have interlacing roots.
            let gp_list: Vec<(u8, Vec<i64>)> = vp
                .iter()
                .filter_map(|&p| gp_polys.get(&p).map(|g| (p, g.clone())))
                .collect();
            let mut all_compat = true;
            if gp_list.len() >= 2 {
                for i in 0..gp_list.len() {
                    for j in (i + 1)..gp_list.len() {
                        let (pi, ref gi) = gp_list[i];
                        let (pj, ref gj) = gp_list[j];
                        if gi.len() <= 1 || gj.len() <= 1 {
                            continue;
                        }
                        // Check interlacing (either direction)
                        let (smaller, larger) = if gi.len() <= gj.len() {
                            (gi.as_slice(), gj.as_slice())
                        } else {
                            (gj.as_slice(), gi.as_slice())
                        };
                        let interlace = match check_weak_interlacing(smaller, larger) {
                            Some(true) => true,
                            _ => false,
                        };
                        if !interlace {
                            all_compat = false;
                            // Also test with random positive linear combinations
                            let mut combo_rr = true;
                            let weights = [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (1, 5), (5, 1)];
                            for &(a, b) in &weights {
                                let maxlen = gi.len().max(gj.len());
                                let mut c = vec![0i64; maxlen];
                                for (k, &v) in gi.iter().enumerate() {
                                    c[k] += a * v;
                                }
                                for (k, &v) in gj.iter().enumerate() {
                                    c[k] += b * v;
                                }
                                while c.len() > 1 && *c.last().unwrap() == 0 {
                                    c.pop();
                                }
                                if c.len() > 2 && !is_real_rooted(&c) {
                                    combo_rr = false;
                                    break;
                                }
                            }
                            println!(
                                "    COMPAT FAIL: G_{} vs G_{} (interlace: false, combos rr: {})",
                                pi, pj, combo_rr
                            );
                        }
                    }
                }
                if all_compat {
                    total_compat += 1;
                    println!("    G_p compatible family: YES");
                }
            }

            // Check interlacing CHAIN: G_{p_1} <= G_{p_2} <= ... <= G_{p_m}
            let mut chain_ok = true;
            if gp_list.len() >= 2 {
                for i in 0..gp_list.len() - 1 {
                    let (_, ref gi) = gp_list[i];
                    let (_, ref gj) = gp_list[i + 1];
                    if gi.len() <= 1 || gj.len() <= 1 {
                        continue;
                    }
                    let (smaller, larger) = if gi.len() <= gj.len() {
                        (gi.as_slice(), gj.as_slice())
                    } else {
                        (gj.as_slice(), gi.as_slice())
                    };
                    match check_weak_interlacing(smaller, larger) {
                        Some(true) => {}
                        _ => {
                            chain_ok = false;
                        }
                    }
                }
                if chain_ok {
                    total_chain += 1;
                }
                println!(
                    "    G_p interlacing chain: {}",
                    if chain_ok { "YES" } else { "NO" }
                );
            }

            // Check: for fixed p, do the f_{p,q} form an interlacing chain in q?
            // Also: for fixed q, do the f_{p,q} form an interlacing chain in p?
            {
                let mut chain_q_ok = true;
                let mut chain_p_ok = true;

                // For each p: collect f_{p,q} sorted by q
                for &p in &vp {
                    let mut q_polys: Vec<(u8, Vec<i64>)> = fpq
                        .iter()
                        .filter(|&(&(pp, _), _)| pp == p)
                        .map(|(&(_, q), vals)| (q, build_poly(vals)))
                        .filter(|(_, poly)| poly.len() > 1)
                        .collect();
                    q_polys.sort_by_key(|&(q, _)| q);
                    for i in 0..q_polys.len().saturating_sub(1) {
                        let (_, ref a) = q_polys[i];
                        let (_, ref b) = q_polys[i + 1];
                        if a.len() > 1 && b.len() > 1 {
                            let (s, l) = if a.len() <= b.len() { (a.as_slice(), b.as_slice()) } else { (b.as_slice(), a.as_slice()) };
                            match check_weak_interlacing(s, l) {
                                Some(true) => {}
                                _ => { chain_q_ok = false; }
                            }
                        }
                    }
                }

                // For each q: collect f_{p,q} sorted by p
                let mut all_qs: Vec<u8> = fpq.keys().map(|&(_, q)| q).collect();
                all_qs.sort();
                all_qs.dedup();
                for &q in &all_qs {
                    let mut p_polys: Vec<(u8, Vec<i64>)> = fpq
                        .iter()
                        .filter(|&(&(_, qq), _)| qq == q)
                        .map(|(&(p, _), vals)| (p, build_poly(vals)))
                        .filter(|(_, poly)| poly.len() > 1)
                        .collect();
                    p_polys.sort_by_key(|&(p, _)| p);
                    for i in 0..p_polys.len().saturating_sub(1) {
                        let (_, ref a) = p_polys[i];
                        let (_, ref b) = p_polys[i + 1];
                        if a.len() > 1 && b.len() > 1 {
                            let (s, l) = if a.len() <= b.len() { (a.as_slice(), b.as_slice()) } else { (b.as_slice(), a.as_slice()) };
                            match check_weak_interlacing(s, l) {
                                Some(true) => {}
                                _ => { chain_p_ok = false; }
                            }
                        }
                    }
                }

                println!(
                    "    f_{{p,q}} chain in q (fixed p): {}  chain in p (fixed q): {}",
                    if chain_q_ok { "YES" } else { "NO" },
                    if chain_p_ok { "YES" } else { "NO" }
                );
            }

            // Full Phi stability
            let (phi_stable, phi_min) = check_stability(
                |x, t| eval_phi_direct(&class, n, x, t),
                num_tests,
                mask.wrapping_mul(111) + n as u64,
            );
            if phi_stable {
                total_phi_stable += 1;
            }
            println!(
                "    Phi: {} (min={:.2e})",
                if phi_stable { "STABLE" } else { "UNSTABLE" },
                phi_min
            );

            // For the full matrix f_{p,q}(t): display its structure for small cases
            if vp.len() <= 5 {
                // Get all q values that appear
                let mut all_qs: Vec<u8> = fpq.keys().map(|&(_, q)| q).collect();
                all_qs.sort();
                all_qs.dedup();

                println!("    f_{{p,q}} matrix (rows=q, cols=p):");
                for &q in &all_qs {
                    let row: Vec<String> = vp
                        .iter()
                        .map(|&p| {
                            fpq.get(&(p, q))
                                .map(|vals| format_poly(&build_poly(vals)))
                                .unwrap_or_else(|| "0".to_string())
                        })
                        .collect();
                    let e2_row: Vec<String> = vp
                        .iter()
                        .map(|&p| format!("e2={}", epsilon2(p, q)))
                        .collect();
                    println!("      q={}: [{}]  eps2=[{}]", q, row.join(" | "), e2_row.join(","));
                }
            }
        }
        println!();
    }

    println!("{}", "=".repeat(70));
    println!("SUMMARY (n=5..{})", max_n);
    println!("{}", "=".repeat(70));
    println!("Total descent sets tested:     {}", total_cases);
    println!(
        "A_p real-rooted:               {}/{}",
        total_ap_rr, total_ap_bp_tested
    );
    println!(
        "B_p real-rooted:               {}/{}",
        total_bp_rr, total_ap_bp_tested
    );
    println!(
        "G_p real-rooted:               {}/{}",
        total_gp_rr, total_ap_bp_tested
    );
    println!(
        "B_p interlaces A_p:            {}/{}",
        total_bp_int_ap, total_ap_bp_tested
    );
    println!(
        "A_p interlaces B_p:            {}/{}",
        total_ap_int_bp, total_ap_bp_tested
    );
    println!(
        "f_{{p,q}} real-rooted:            {}/{}",
        total_fpq_rr, total_fpq_tested
    );
    println!(
        "G_p compatible family:         {}/{}",
        total_compat, total_cases
    );
    println!(
        "G_p interlacing chain:         {}/{}",
        total_chain, total_cases
    );
    println!(
        "Full Phi stable:               {}/{}",
        total_phi_stable, total_cases
    );
    println!();
    if total_ap_rr == total_ap_bp_tested && total_bp_rr == total_ap_bp_tested {
        println!("*** Both A_p and B_p are ALWAYS real-rooted. ***");
    }
    if total_fpq_rr == total_fpq_tested {
        println!("*** All source polys f_{{p,q}}(t) are real-rooted. ***");
    }
    if total_compat == total_cases {
        println!("*** G_p form a COMPATIBLE FAMILY for all cases! ***");
        println!("    This directly implies Phi(x,t) is stable by the Chudnovsky-Seymour theorem.");
    }
    if total_chain == total_cases {
        println!("*** G_p form an interlacing CHAIN for all cases! ***");
    }
    if total_phi_stable == total_cases {
        println!("*** Full Phi is STABLE for all tested cases. ***");
    }
}
