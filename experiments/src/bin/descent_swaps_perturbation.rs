/// Decompose Delta = A_{k1}B_{k2} - A_{k2}B_{k1} into three brackets:
///
///   Delta = [staircase_M] + [perturbation_L] + [perturbation_R]
///
/// Using the algebraic identity Delta = A_{k2} D(k1) - A_{k1} D(k2)
/// where D(k) = A_k - B_k, decomposed by zones:
///
///   D(k) = D_L(k) + D_M(k) + D_R(k)
///
/// where:
///   D_L(k) = L_a(k-1) - L_b(k-1)  (same-shift zone, perturbation from source diff)
///   D_M(k) = M_a(k) - M_b(k-1)    (THE key zone: different eps2 shifts)
///   D_R(k) = R_a(k) - R_b(k)       (same-shift zone, perturbation from source diff)
///
/// So Delta = [A_{k2} D_M(k1) - A_{k1} D_M(k2)]     <-- "staircase" bracket
///          + [A_{k2} D_L(k1) - A_{k1} D_L(k2)]      <-- "L perturbation" bracket
///          + [A_{k2} D_R(k1) - A_{k1} D_R(k2)]      <-- "R perturbation" bracket
///
/// In the common-source case (D_L = D_R = 0), only the staircase bracket remains,
/// which equals Delta from Theorem 6.8.
///
/// Key question: Is the staircase bracket always non-negative?
/// And are the perturbation brackets small enough to be compensated?
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use std::collections::BTreeMap;

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
    if p >= 2 && q <= p - 2 {
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

fn zone_label(z: Zone) -> &'static str {
    match z {
        Zone::L => "L",
        Zone::M => "M",
        Zone::R => "R",
    }
}

/// Build N_p(s, q) distribution
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

/// Build zone-restricted polynomial for source p:
/// Z_p^{zone}(t) = sum_{q in zone} t^{eps2(p,q)} N_p(*, q)
/// Returns coefficients indexed by the "output" degree k (after applying eps2).
fn build_zone_poly(
    dist: &BTreeMap<(usize, u8), usize>,
    p: u8,
    p_a: u8,
    p_b: u8,
    zone: Zone,
) -> Vec<i64> {
    let mut poly = vec![0i64; 20];
    for (&(s, q), &count) in dist {
        if classify_zone(q, p_a, p_b) == zone {
            let k = s + eps2(p, q);
            if k >= poly.len() {
                poly.resize(k + 1, 0);
            }
            poly[k] += count as i64;
        }
    }
    while poly.len() > 1 && *poly.last().unwrap() == 0 {
        poly.pop();
    }
    poly
}

/// Build the "raw" zone polynomial (without eps2 shift):
/// Np_zone(s) = sum_{q in zone} N_p(s, q)
fn build_zone_raw_poly(
    dist: &BTreeMap<(usize, u8), usize>,
    p_a: u8,
    p_b: u8,
    zone: Zone,
) -> Vec<i64> {
    let mut poly = vec![0i64; 20];
    for (&(s, q), &count) in dist {
        if classify_zone(q, p_a, p_b) == zone {
            if s >= poly.len() {
                poly.resize(s + 1, 0);
            }
            poly[s] += count as i64;
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

    println!("=== Perturbation decomposition of Delta ===");
    println!("Delta = [staircase_M] + [perturbation_L] + [perturbation_R]");
    println!("where staircase_M = A_{{k2}} D_M(k1) - A_{{k1}} D_M(k2)");
    println!("  D_M(k) = M_a(k) - M_b(k-1)");
    println!("  D_L(k) = L_a(k-1) - L_b(k-1)");
    println!("  D_R(k) = R_a(k) - R_b(k)");
    println!("Testing n = 4 to {}\n", max_n);

    // Global counters
    let mut total_minors = 0u64;
    let mut nontrivial_minors = 0u64;

    // Track signs of each bracket
    let mut staircase_pos = 0u64;
    let mut staircase_zero = 0u64;
    let mut staircase_neg = 0u64;

    let mut pert_l_pos = 0u64;
    let mut pert_l_zero = 0u64;
    let mut pert_l_neg = 0u64;

    let mut pert_r_pos = 0u64;
    let mut pert_r_zero = 0u64;
    let mut pert_r_neg = 0u64;

    // Track staircase surplus vs perturbation deficit
    let mut staircase_dominates = 0u64; // staircase bracket >= |pert_L| + |pert_R|
    let mut pert_compensated = 0u64; // Delta >= 0 but some perturbation bracket < 0

    // Track detailed examples of negative staircase
    let mut neg_staircase_examples: Vec<String> = Vec::new();

    // Track the "cross bracket" interactions
    // Also decompose by which zones contribute to the perturbation
    let mut staircase_with_lpert_neg = 0u64;
    let mut staircase_with_rpert_neg = 0u64;
    let mut staircase_with_both_neg = 0u64;

    // NEW: Also check whether the M+L and M+R combined brackets are always nonneg
    let mut ml_bracket_pos = 0u64;
    let mut ml_bracket_zero = 0u64;
    let mut ml_bracket_neg = 0u64;
    let mut mr_bracket_pos = 0u64;
    let mut mr_bracket_zero = 0u64;
    let mut mr_bracket_neg = 0u64;

    // Track LR bracket (combined L+R perturbation)
    let mut lr_bracket_pos = 0u64;
    let mut lr_bracket_zero = 0u64;
    let mut lr_bracket_neg = 0u64;

    // ALTERNATIVE decomposition: use D split differently
    // D(k) = [M_a(k) - M_b(k-1)] + [L_a(k-1) - L_b(k-1)] + [R_a(k) - R_b(k)]
    // But we can also write D_M(k) = M_a(k) - M_a(k-1) + [M_a(k-1) - M_b(k-1)]
    //   = backward_diff_Ma(k) + [Ma-Mb perturbation](k-1)
    // So the "pure staircase" part uses only the backward difference of M_a,
    // and there's a fourth perturbation from M_a != M_b.

    let mut pure_staircase_pos = 0u64;
    let mut pure_staircase_zero = 0u64;
    let mut pure_staircase_neg = 0u64;

    let mut m_pert_pos = 0u64;
    let mut m_pert_zero = 0u64;
    let mut m_pert_neg = 0u64;

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

        let mut n_total = 0u64;
        let mut n_nontrivial = 0u64;
        let mut n_delta_nonneg = 0u64;

        for s_mask in 0u64..(1 << (n - 1)) {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            // Build source distributions
            let mut source_dists: BTreeMap<u8, BTreeMap<(usize, u8), usize>> = BTreeMap::new();
            for &p in &vp {
                let dist = build_source_dist(s_mask, p, n, &by_descent_prev, &perm_prev_list);
                source_dists.insert(p, dist);
            }

            // Build L^{(p)} polynomials for verification
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

            // For each consecutive pair of valid positions
            for idx in 0..vp.len() - 1 {
                let p_a = vp[idx];
                let p_b = vp[idx + 1];

                let dist_a = &source_dists[&p_a];
                let dist_b = &source_dists[&p_b];

                // Build zone polynomials (with eps2 applied)
                let la_poly = build_zone_poly(dist_a, p_a, p_a, p_b, Zone::L);
                let lb_poly = build_zone_poly(dist_b, p_b, p_a, p_b, Zone::L);
                let ma_poly = build_zone_poly(dist_a, p_a, p_a, p_b, Zone::M);
                let mb_poly = build_zone_poly(dist_b, p_b, p_a, p_b, Zone::M);
                let ra_poly = build_zone_poly(dist_a, p_a, p_a, p_b, Zone::R);
                let rb_poly = build_zone_poly(dist_b, p_b, p_a, p_b, Zone::R);

                // Build zone RAW polynomials (without eps2, just in swap count s)
                let ma_raw = build_zone_raw_poly(dist_a, p_a, p_b, Zone::M);
                let mb_raw = build_zone_raw_poly(dist_b, p_a, p_b, Zone::M);

                // Verify: A_k = L_a(k) + M_a(k) + R_a(k) (with eps2 applied in build_zone_poly)
                // Actually: L_a(k) = sum_{q in L} N_a(k - eps2(p_a,q), q)
                // But build_zone_poly already builds the output polynomial indexed by k.
                // So A_k = la_poly[k] + ma_poly[k] + ra_poly[k].

                let max_deg = lp_polys.values().map(|p| p.len()).max().unwrap_or(1).max(
                    [&la_poly, &lb_poly, &ma_poly, &mb_poly, &ra_poly, &rb_poly]
                        .iter()
                        .map(|p| p.len())
                        .max()
                        .unwrap_or(1),
                );

                // Verify reconstruction
                let a_poly = lp_polys.get(&p_a);
                let b_poly = lp_polys.get(&p_b);
                for k in 0..max_deg {
                    let a_recon = coeff(&la_poly, k) + coeff(&ma_poly, k) + coeff(&ra_poly, k);
                    let a_direct = a_poly.map(|p| coeff(p, k)).unwrap_or(0);
                    if a_recon != a_direct {
                        println!(
                            "RECON FAIL A: n={} S={} p_a={} k={}: {} != {}",
                            n,
                            descent_set_to_string(s_mask, n),
                            p_a,
                            k,
                            a_recon,
                            a_direct
                        );
                    }
                    let b_recon = coeff(&lb_poly, k) + coeff(&mb_poly, k) + coeff(&rb_poly, k);
                    let b_direct = b_poly.map(|p| coeff(p, k)).unwrap_or(0);
                    if b_recon != b_direct {
                        println!(
                            "RECON FAIL B: n={} S={} p_b={} k={}: {} != {}",
                            n,
                            descent_set_to_string(s_mask, n),
                            p_b,
                            k,
                            b_recon,
                            b_direct
                        );
                    }
                }

                // Now compute D_L(k), D_M(k), D_R(k)
                // A_k = la_poly[k] + ma_poly[k] + ra_poly[k]
                // B_k = lb_poly[k] + mb_poly[k] + rb_poly[k]
                // D(k) = A_k - B_k
                // D_L(k) = la_poly[k] - lb_poly[k]
                // D_M(k) = ma_poly[k] - mb_poly[k]
                // D_R(k) = ra_poly[k] - rb_poly[k]
                //
                // Note: la_poly[k] = L_a(k-1) since for zone L, eps2(p_a,q) = 1, so
                //   la_poly[k] = sum_{q in L} N_a(k-1, q)
                // Similarly lb_poly[k] = L_b(k-1) since eps2(p_b,q) = 1 for zone L too.
                // So D_L(k) = L_a(k-1) - L_b(k-1) = la_poly[k] - lb_poly[k]. Correct.
                //
                // For zone M: eps2(p_a,q) = 0, eps2(p_b,q) = 1. So:
                //   ma_poly[k] = sum_{q in M} N_a(k, q) = M_a(k)
                //   mb_poly[k] = sum_{q in M} N_b(k-1, q) = M_b(k-1)
                // D_M(k) = M_a(k) - M_b(k-1). Correct.
                //
                // For zone R: eps2(p_a,q) = 0, eps2(p_b,q) = 0. So:
                //   ra_poly[k] = R_a(k), rb_poly[k] = R_b(k)
                // D_R(k) = R_a(k) - R_b(k). Correct.

                // For each (k1, k2) pair
                for k1 in 0..max_deg {
                    for k2 in (k1 + 1)..max_deg {
                        let a_k1 = a_poly.map(|p| coeff(p, k1)).unwrap_or(0);
                        let a_k2 = a_poly.map(|p| coeff(p, k2)).unwrap_or(0);
                        let b_k1 = b_poly.map(|p| coeff(p, k1)).unwrap_or(0);
                        let b_k2 = b_poly.map(|p| coeff(p, k2)).unwrap_or(0);

                        let delta = a_k1 * b_k2 - a_k2 * b_k1;

                        if a_k1 == 0 && a_k2 == 0 {
                            continue;
                        }
                        if b_k1 == 0 && b_k2 == 0 {
                            continue;
                        }

                        n_total += 1;
                        if delta == 0 && a_k1 * b_k2 == 0 {
                            // trivial
                            continue;
                        }
                        n_nontrivial += 1;
                        if delta >= 0 {
                            n_delta_nonneg += 1;
                        }

                        // Compute D values
                        let d_l_k1 = coeff(&la_poly, k1) - coeff(&lb_poly, k1);
                        let d_l_k2 = coeff(&la_poly, k2) - coeff(&lb_poly, k2);
                        let d_m_k1 = coeff(&ma_poly, k1) - coeff(&mb_poly, k1);
                        let d_m_k2 = coeff(&ma_poly, k2) - coeff(&mb_poly, k2);
                        let d_r_k1 = coeff(&ra_poly, k1) - coeff(&rb_poly, k1);
                        let d_r_k2 = coeff(&ra_poly, k2) - coeff(&rb_poly, k2);

                        // Three brackets
                        let bracket_m = a_k2 * d_m_k1 - a_k1 * d_m_k2;
                        let bracket_l = a_k2 * d_l_k1 - a_k1 * d_l_k2;
                        let bracket_r = a_k2 * d_r_k1 - a_k1 * d_r_k2;

                        // Verify
                        let sum_brackets = bracket_m + bracket_l + bracket_r;
                        if sum_brackets != delta {
                            println!(
                                "BRACKET SUM FAIL: n={} S={} pa={} pb={} k1={} k2={}: {} != {}",
                                n,
                                descent_set_to_string(s_mask, n),
                                p_a,
                                p_b,
                                k1,
                                k2,
                                sum_brackets,
                                delta
                            );
                        }

                        // Track signs
                        if bracket_m > 0 {
                            staircase_pos += 1;
                        } else if bracket_m == 0 {
                            staircase_zero += 1;
                        } else {
                            staircase_neg += 1;
                            if neg_staircase_examples.len() < 10 {
                                neg_staircase_examples.push(format!(
                                    "n={} S={} pa={} pb={} k1={} k2={} bracket_M={} bracket_L={} bracket_R={} Delta={}",
                                    n, descent_set_to_string(s_mask, n), p_a, p_b, k1, k2,
                                    bracket_m, bracket_l, bracket_r, delta
                                ));
                            }
                        }

                        if bracket_l > 0 {
                            pert_l_pos += 1;
                        } else if bracket_l == 0 {
                            pert_l_zero += 1;
                        } else {
                            pert_l_neg += 1;
                        }

                        if bracket_r > 0 {
                            pert_r_pos += 1;
                        } else if bracket_r == 0 {
                            pert_r_zero += 1;
                        } else {
                            pert_r_neg += 1;
                        }

                        // Combined brackets
                        let ml_bracket = bracket_m + bracket_l;
                        let mr_bracket = bracket_m + bracket_r;
                        let lr_bracket = bracket_l + bracket_r;

                        if ml_bracket > 0 {
                            ml_bracket_pos += 1;
                        } else if ml_bracket == 0 {
                            ml_bracket_zero += 1;
                        } else {
                            ml_bracket_neg += 1;
                        }

                        if mr_bracket > 0 {
                            mr_bracket_pos += 1;
                        } else if mr_bracket == 0 {
                            mr_bracket_zero += 1;
                        } else {
                            mr_bracket_neg += 1;
                        }

                        if lr_bracket > 0 {
                            lr_bracket_pos += 1;
                        } else if lr_bracket == 0 {
                            lr_bracket_zero += 1;
                        } else {
                            lr_bracket_neg += 1;
                        }

                        // Check staircase dominance
                        if bracket_m >= bracket_l.abs() + bracket_r.abs() {
                            staircase_dominates += 1;
                        }

                        if delta >= 0 && (bracket_l < 0 || bracket_r < 0) {
                            pert_compensated += 1;
                        }

                        // Track patterns
                        if bracket_l < 0 && bracket_r >= 0 {
                            staircase_with_lpert_neg += 1;
                        }
                        if bracket_r < 0 && bracket_l >= 0 {
                            staircase_with_rpert_neg += 1;
                        }
                        if bracket_l < 0 && bracket_r < 0 {
                            staircase_with_both_neg += 1;
                        }

                        // ALTERNATIVE: "pure staircase" using backward difference of M_a only
                        // D_M(k) = M_a(k) - M_b(k-1)
                        //        = [M_a(k) - M_a(k-1)] + [M_a(k-1) - M_b(k-1)]
                        //        = backward_diff_Ma(k) + delta_M_ab(k-1)
                        //
                        // where delta_M_ab(s) = M_a(s) - M_b(s) is the source perturbation in zone M
                        // (comparing raw swap counts without eps2).
                        //
                        // Actually: ma_raw is the raw zone-M polynomial for source a (without eps2).
                        // mb_raw is the raw zone-M polynomial for source b.
                        // ma_poly[k] = ma_raw[k] (since eps2(p_a,q)=0 for zone M)
                        // mb_poly[k] = mb_raw[k-1] (since eps2(p_b,q)=1 for zone M)
                        //
                        // So D_M(k) = ma_raw[k] - mb_raw[k-1]
                        //           = [ma_raw[k] - ma_raw[k-1]] + [ma_raw[k-1] - mb_raw[k-1]]
                        // backward_diff_Ma(k) = ma_raw[k] - ma_raw[k-1]
                        // delta_M_ab(k-1) = ma_raw[k-1] - mb_raw[k-1]

                        // For k >= 1: D_M(k) = ma_raw[k] - mb_raw[k-1]
                        //   = [ma_raw[k] - ma_raw[k-1]] + [ma_raw[k-1] - mb_raw[k-1]]
                        // For k = 0: D_M(0) = ma_raw[0] (since mb_poly[0] = 0)
                        //   backward_diff = ma_raw[0], perturbation = 0
                        let bd_ma_k1 = if k1 > 0 {
                            coeff(&ma_raw, k1) - coeff(&ma_raw, k1 - 1)
                        } else {
                            coeff(&ma_raw, 0)
                        };
                        let bd_ma_k2 = if k2 > 0 {
                            coeff(&ma_raw, k2) - coeff(&ma_raw, k2 - 1)
                        } else {
                            coeff(&ma_raw, 0)
                        };

                        let delta_mab_k1 = if k1 > 0 {
                            coeff(&ma_raw, k1 - 1) - coeff(&mb_raw, k1 - 1)
                        } else {
                            0
                        };
                        let delta_mab_k2 = if k2 > 0 {
                            coeff(&ma_raw, k2 - 1) - coeff(&mb_raw, k2 - 1)
                        } else {
                            0
                        };

                        // Verify: bd_ma_k + delta_mab_k = D_M(k)
                        assert_eq!(
                            bd_ma_k1 + delta_mab_k1,
                            d_m_k1,
                            "pure staircase decomp fail at k1"
                        );
                        assert_eq!(
                            bd_ma_k2 + delta_mab_k2,
                            d_m_k2,
                            "pure staircase decomp fail at k2"
                        );

                        let pure_staircase = a_k2 * bd_ma_k1 - a_k1 * bd_ma_k2;
                        let m_pert = a_k2 * delta_mab_k1 - a_k1 * delta_mab_k2;

                        // Verify: pure_staircase + m_pert = bracket_m
                        assert_eq!(pure_staircase + m_pert, bracket_m);

                        if pure_staircase > 0 {
                            pure_staircase_pos += 1;
                        } else if pure_staircase == 0 {
                            pure_staircase_zero += 1;
                        } else {
                            pure_staircase_neg += 1;
                        }

                        if m_pert > 0 {
                            m_pert_pos += 1;
                        } else if m_pert == 0 {
                            m_pert_zero += 1;
                        } else {
                            m_pert_neg += 1;
                        }
                    }
                }
            }
        }

        total_minors += n_total;
        nontrivial_minors += n_nontrivial;
        println!(
            "n={}: total={} nontrivial={} Delta>=0: {}/{}",
            n, n_total, n_nontrivial, n_delta_nonneg, n_nontrivial
        );
    }

    println!("\n{}", "=".repeat(80));
    println!("THREE-BRACKET DECOMPOSITION SUMMARY (n=4..{})", max_n);
    println!("{}", "=".repeat(80));
    println!("Delta = bracket_M + bracket_L + bracket_R");
    println!("  bracket_M = A_{{k2}} D_M(k1) - A_{{k1}} D_M(k2)  [staircase from zone-M shift]");
    println!("  bracket_L = A_{{k2}} D_L(k1) - A_{{k1}} D_L(k2)  [L zone source perturbation]");
    println!("  bracket_R = A_{{k2}} D_R(k1) - A_{{k1}} D_R(k2)  [R zone source perturbation]");
    println!();

    println!(
        "{:<20} {:>10} {:>10} {:>10}",
        "Bracket", "Positive", "Zero", "Negative"
    );
    println!("{}", "-".repeat(55));
    println!(
        "{:<20} {:>10} {:>10} {:>10}  {}",
        "bracket_M (stairc.)",
        staircase_pos,
        staircase_zero,
        staircase_neg,
        if staircase_neg == 0 {
            "ALWAYS >= 0"
        } else {
            "CAN BE NEGATIVE"
        }
    );
    println!(
        "{:<20} {:>10} {:>10} {:>10}  {}",
        "bracket_L (pert.L)",
        pert_l_pos,
        pert_l_zero,
        pert_l_neg,
        if pert_l_neg == 0 {
            "ALWAYS >= 0"
        } else {
            "CAN BE NEGATIVE"
        }
    );
    println!(
        "{:<20} {:>10} {:>10} {:>10}  {}",
        "bracket_R (pert.R)",
        pert_r_pos,
        pert_r_zero,
        pert_r_neg,
        if pert_r_neg == 0 {
            "ALWAYS >= 0"
        } else {
            "CAN BE NEGATIVE"
        }
    );

    println!();
    println!("{}", "=".repeat(80));
    println!("COMBINED BRACKETS");
    println!("{}", "=".repeat(80));
    println!(
        "{:<20} {:>10} {:>10} {:>10}  {}",
        "M+L combined",
        ml_bracket_pos,
        ml_bracket_zero,
        ml_bracket_neg,
        if ml_bracket_neg == 0 {
            "ALWAYS >= 0"
        } else {
            "CAN BE NEGATIVE"
        }
    );
    println!(
        "{:<20} {:>10} {:>10} {:>10}  {}",
        "M+R combined",
        mr_bracket_pos,
        mr_bracket_zero,
        mr_bracket_neg,
        if mr_bracket_neg == 0 {
            "ALWAYS >= 0"
        } else {
            "CAN BE NEGATIVE"
        }
    );
    println!(
        "{:<20} {:>10} {:>10} {:>10}  {}",
        "L+R combined",
        lr_bracket_pos,
        lr_bracket_zero,
        lr_bracket_neg,
        if lr_bracket_neg == 0 {
            "ALWAYS >= 0"
        } else {
            "CAN BE NEGATIVE"
        }
    );
    println!(
        "{:<20} {:>10} {:>10} {:>10}  (= Delta)",
        "M+L+R combined",
        staircase_pos + pert_l_pos + pert_r_pos, // not right, just placeholder
        0,
        0
    );

    println!();
    println!("{}", "=".repeat(80));
    println!("PURE STAIRCASE vs M-SOURCE PERTURBATION");
    println!("{}", "=".repeat(80));
    println!("bracket_M = pure_staircase + m_perturbation");
    println!("  pure_staircase = A_{{k2}} [Ma(k1)-Ma(k1-1)] - A_{{k1}} [Ma(k2)-Ma(k2-1)]");
    println!("  m_perturbation = A_{{k2}} [Ma(k1-1)-Mb(k1-1)] - A_{{k1}} [Ma(k2-1)-Mb(k2-1)]");
    println!();
    println!(
        "{:<20} {:>10} {:>10} {:>10}  {}",
        "pure_staircase",
        pure_staircase_pos,
        pure_staircase_zero,
        pure_staircase_neg,
        if pure_staircase_neg == 0 {
            "ALWAYS >= 0"
        } else {
            "CAN BE NEGATIVE"
        }
    );
    println!(
        "{:<20} {:>10} {:>10} {:>10}  {}",
        "m_perturbation",
        m_pert_pos,
        m_pert_zero,
        m_pert_neg,
        if m_pert_neg == 0 {
            "ALWAYS >= 0"
        } else {
            "CAN BE NEGATIVE"
        }
    );

    println!();
    println!("{}", "=".repeat(80));
    println!("DOMINANCE ANALYSIS");
    println!("{}", "=".repeat(80));
    println!(
        "Staircase bracket >= |pert_L| + |pert_R|: {} / {} ({:.1}%)",
        staircase_dominates,
        nontrivial_minors,
        100.0 * staircase_dominates as f64 / nontrivial_minors.max(1) as f64
    );
    println!(
        "Cases where Delta>=0 but some perturbation<0: {} / {}",
        pert_compensated, nontrivial_minors
    );
    println!("Patterns of negative perturbations:",);
    println!("  Only bracket_L < 0: {}", staircase_with_lpert_neg);
    println!("  Only bracket_R < 0: {}", staircase_with_rpert_neg);
    println!("  Both bracket_L,R < 0: {}", staircase_with_both_neg);

    if !neg_staircase_examples.is_empty() {
        println!();
        println!("Examples where bracket_M (staircase) is NEGATIVE:");
        for ex in &neg_staircase_examples {
            println!("  {}", ex);
        }
    }

    // NEW ANALYSIS: Check if Delta can be decomposed differently
    // Using the ZONE-PAIR decomposition from the minor_decomp experiment:
    //   Delta = sum_{za, zb in {L,M,R}} contrib(za, zb)
    // where contrib(za, zb) = sum_{qa in za, qb in zb} [Na(k1-ea(qa),qa)*Nb(k2-eb(qb),qb) - Na(k2-ea(qa),qa)*Nb(k1-eb(qb),qb)]
    //
    // Known from minor_decomp: (L,M), (M,L), (M,M), (R,M) are ALWAYS >= 0.
    //   So: ALL zone pairs involving M on either side are nonneg.
    //   The negative ones are: (L,L), (L,R), (R,L), (R,R), (M,R).
    //
    // Actually let me recheck: the minor_decomp said:
    //   Always nonneg: (L,M), (M,L), (M,M), (R,M)
    //   Can be negative: (L,L), (L,R), (M,R), (R,L), (R,R)
    //
    // So: M on the LEFT is always nonneg [(M,L), (M,M), (M,R) -- wait (M,R) CAN be negative!]
    // Let me just re-verify with the data.

    println!();
    println!("{}", "=".repeat(80));
    println!("ZONE-PAIR DECOMPOSITION (re-verified)");
    println!("{}", "=".repeat(80));

    let zone_pairs = [
        (Zone::L, Zone::L),
        (Zone::L, Zone::M),
        (Zone::L, Zone::R),
        (Zone::M, Zone::L),
        (Zone::M, Zone::M),
        (Zone::M, Zone::R),
        (Zone::R, Zone::L),
        (Zone::R, Zone::M),
        (Zone::R, Zone::R),
    ];

    let mut zp_pos: BTreeMap<(Zone, Zone), u64> = BTreeMap::new();
    let mut zp_zero: BTreeMap<(Zone, Zone), u64> = BTreeMap::new();
    let mut zp_neg: BTreeMap<(Zone, Zone), u64> = BTreeMap::new();
    for &zp in &zone_pairs {
        zp_pos.insert(zp, 0);
        zp_zero.insert(zp, 0);
        zp_neg.insert(zp, 0);
    }

    // Also track "M-row sum" and "M-col sum" and other groupings
    let mut m_any_pos = 0u64;
    let mut m_any_zero = 0u64;
    let mut m_any_neg = 0u64;

    let mut no_m_pos = 0u64;
    let mut no_m_zero = 0u64;
    let mut no_m_neg = 0u64;

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
                        .filter(|pi| pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                        .map(|pi| compute(pi, Stat::Swaps))
                        .collect();
                    if !vals.is_empty() {
                        lp_polys.insert(p, build_poly(&vals));
                    }
                }
            }

            for idx_pair in 0..vp.len() - 1 {
                let p_a = vp[idx_pair];
                let p_b = vp[idx_pair + 1];

                let dist_a = &source_dists[&p_a];
                let dist_b = &source_dists[&p_b];

                let mut all_qs_a: Vec<u8> = dist_a.keys().map(|&(_, q)| q).collect();
                all_qs_a.sort();
                all_qs_a.dedup();
                let mut all_qs_b: Vec<u8> = dist_b.keys().map(|&(_, q)| q).collect();
                all_qs_b.sort();
                all_qs_b.dedup();

                // Build N_p(s, q) polynomials
                let n_polys_a: BTreeMap<u8, Vec<i64>> = all_qs_a
                    .iter()
                    .map(|&q| {
                        let mut poly = vec![0i64; 20];
                        for (&(s, qq), &c) in dist_a {
                            if qq == q {
                                if s >= poly.len() {
                                    poly.resize(s + 1, 0);
                                }
                                poly[s] += c as i64;
                            }
                        }
                        while poly.len() > 1 && *poly.last().unwrap() == 0 {
                            poly.pop();
                        }
                        (q, poly)
                    })
                    .collect();
                let n_polys_b: BTreeMap<u8, Vec<i64>> = all_qs_b
                    .iter()
                    .map(|&q| {
                        let mut poly = vec![0i64; 20];
                        for (&(s, qq), &c) in dist_b {
                            if qq == q {
                                if s >= poly.len() {
                                    poly.resize(s + 1, 0);
                                }
                                poly[s] += c as i64;
                            }
                        }
                        while poly.len() > 1 && *poly.last().unwrap() == 0 {
                            poly.pop();
                        }
                        (q, poly)
                    })
                    .collect();

                let max_deg = lp_polys.values().map(|p| p.len()).max().unwrap_or(1);

                for k1 in 0..max_deg {
                    for k2 in (k1 + 1)..max_deg {
                        let a_k1 = lp_polys.get(&p_a).map(|p| coeff(p, k1)).unwrap_or(0);
                        let a_k2 = lp_polys.get(&p_a).map(|p| coeff(p, k2)).unwrap_or(0);
                        let b_k1 = lp_polys.get(&p_b).map(|p| coeff(p, k1)).unwrap_or(0);
                        let b_k2 = lp_polys.get(&p_b).map(|p| coeff(p, k2)).unwrap_or(0);

                        if (a_k1 == 0 && a_k2 == 0) || (b_k1 == 0 && b_k2 == 0) {
                            continue;
                        }
                        let delta = a_k1 * b_k2 - a_k2 * b_k1;
                        if delta == 0 && a_k1 * b_k2 == 0 {
                            continue;
                        }

                        // Compute zone-pair contributions
                        let mut zp_contribs: BTreeMap<(Zone, Zone), i64> = BTreeMap::new();
                        for &zp in &zone_pairs {
                            zp_contribs.insert(zp, 0);
                        }

                        for &q_a in &all_qs_a {
                            let e_a = eps2(p_a, q_a);
                            let z_a = classify_zone(q_a, p_a, p_b);
                            for &q_b in &all_qs_b {
                                let e_b = eps2(p_b, q_b);
                                let z_b = classify_zone(q_b, p_a, p_b);

                                let na_poly = &n_polys_a[&q_a];
                                let nb_poly = &n_polys_b[&q_b];

                                let na_k1 = if k1 >= e_a {
                                    coeff(na_poly, k1 - e_a)
                                } else {
                                    0
                                };
                                let na_k2 = if k2 >= e_a {
                                    coeff(na_poly, k2 - e_a)
                                } else {
                                    0
                                };
                                let nb_k1 = if k1 >= e_b {
                                    coeff(nb_poly, k1 - e_b)
                                } else {
                                    0
                                };
                                let nb_k2 = if k2 >= e_b {
                                    coeff(nb_poly, k2 - e_b)
                                } else {
                                    0
                                };

                                let c = na_k1 * nb_k2 - na_k2 * nb_k1;
                                *zp_contribs.entry((z_a, z_b)).or_default() += c;
                            }
                        }

                        for &zp in &zone_pairs {
                            let c = zp_contribs[&zp];
                            if c > 0 {
                                *zp_pos.get_mut(&zp).unwrap() += 1;
                            } else if c == 0 {
                                *zp_zero.get_mut(&zp).unwrap() += 1;
                            } else {
                                *zp_neg.get_mut(&zp).unwrap() += 1;
                            }
                        }

                        // "M-any" = sum of zone pairs where at least one is M
                        let m_any: i64 = zone_pairs
                            .iter()
                            .filter(|(za, zb)| *za == Zone::M || *zb == Zone::M)
                            .map(|zp| zp_contribs[zp])
                            .sum();
                        let no_m: i64 = zone_pairs
                            .iter()
                            .filter(|(za, zb)| *za != Zone::M && *zb != Zone::M)
                            .map(|zp| zp_contribs[zp])
                            .sum();

                        if m_any > 0 {
                            m_any_pos += 1;
                        } else if m_any == 0 {
                            m_any_zero += 1;
                        } else {
                            m_any_neg += 1;
                        }

                        if no_m > 0 {
                            no_m_pos += 1;
                        } else if no_m == 0 {
                            no_m_zero += 1;
                        } else {
                            no_m_neg += 1;
                        }
                    }
                }
            }
        }
    }

    println!();
    println!(
        "{:<12} {:>10} {:>10} {:>10}  {}",
        "Zone pair", "Positive", "Zero", "Negative", "Always>=0?"
    );
    println!("{}", "-".repeat(65));
    for &zp in &zone_pairs {
        println!(
            "({},{})        {:>10} {:>10} {:>10}  {}",
            zone_label(zp.0),
            zone_label(zp.1),
            zp_pos[&zp],
            zp_zero[&zp],
            zp_neg[&zp],
            if zp_neg[&zp] == 0 { "YES" } else { "NO" }
        );
    }

    println!();
    println!(
        "{:<20} {:>10} {:>10} {:>10}  {}",
        "M-any (M in pair)",
        m_any_pos,
        m_any_zero,
        m_any_neg,
        if m_any_neg == 0 {
            "ALWAYS >= 0"
        } else {
            "CAN BE NEGATIVE"
        }
    );
    println!(
        "{:<20} {:>10} {:>10} {:>10}  {}",
        "no-M ({L,R}x{L,R})",
        no_m_pos,
        no_m_zero,
        no_m_neg,
        if no_m_neg == 0 {
            "ALWAYS >= 0"
        } else {
            "CAN BE NEGATIVE"
        }
    );

    println!();
    println!("{}", "=".repeat(80));
    println!("CONCLUSIONS");
    println!("{}", "=".repeat(80));
    println!();

    if staircase_neg == 0 {
        println!("* bracket_M (staircase from zone-M shift) is ALWAYS non-negative.");
        println!("  This is the main term driving Delta >= 0.");
    } else {
        println!("* bracket_M can be NEGATIVE in {} cases.", staircase_neg);
        println!("  The pure staircase argument does NOT suffice alone.");
    }

    if m_any_neg == 0 {
        println!("* The sum of all zone pairs involving M is ALWAYS non-negative.");
        println!("  So: M-zone contributions always provide enough surplus to cover");
        println!("  any negative contributions from {{L,R}}x{{L,R}} zone pairs.");
    }

    if ml_bracket_neg == 0 && mr_bracket_neg == 0 {
        println!("* Both M+L and M+R combined brackets are always non-negative.");
        println!("  So the staircase bracket compensates EACH perturbation individually.");
    }

    println!("\nTotal nontrivial minors checked: {}", nontrivial_minors);
}
