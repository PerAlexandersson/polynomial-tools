/// Deep investigation of the two always-nonneg terms:
///
///   (1) staircase(A, t*Tb) >= 0: A_{k1} * Tb(k2-1) >= A_{k2} * Tb(k1-1)
///       where Tb(s) = sum_{q in L+M} N_b(s, q) = zone (L+M) of source B.
///       Equivalently: A >=_{staircase LR} t*Tb.
///
///   (2) standard(Sa, B) >= 0: Sa(k1) * B_{k2} >= Sa(k2) * B_{k1}
///       where Sa(s) = sum_{q in M+R} N_a(s, q) = zone (M+R) of source A.
///       Equivalently: Sa >=_{LR} B.
///
/// Investigation:
///   - Print Sa and B for small cases to see the structure.
///   - Check if Sa is a "left shift" of B or a "left part" of some larger polynomial.
///   - Check if the relationship follows from some structural property.
///
/// For (2), note:
///   B = L^{(pb)} = sum_q t^{eps2(pb,q)} N_b(*, q)
///     = t * [sum_{q in L+M} N_b(*, q)] + [sum_{q in R} N_b(*, q)]
///     = t * Tb + Rb   (as polynomials)
///
///   Sa = Ma + Ra = sum_{q in M+R} N_a(*, q)
///
///   So B_k = Tb(k-1) + Rb(k) and Sa(k) = Ma(k) + Ra(k).
///
///   standard(Sa, B) >= 0 means:
///     Sa(k1) * [Tb(k2-1) + Rb(k2)] >= Sa(k2) * [Tb(k1-1) + Rb(k1)]
///
///   = Sa(k1)*Tb(k2-1) - Sa(k2)*Tb(k1-1) + Sa(k1)*Rb(k2) - Sa(k2)*Rb(k1)
///   = staircase(Sa, t*Tb) + standard(Sa, Rb)
///
///   So standard(Sa, B) = staircase(Sa, t*Tb) + standard(Sa, Rb).
///
///   The staircase(Sa, t*Tb) is related to (but not identical to) the zone-pair
///   staircase (0,1) group. And standard(Sa, Rb) is related to the (M+R, R) zone pairs.
///
///   Can we show Sa >=_{LR} B directly? Sa is the "unshifted part" of A
///   (the part that doesn't get the t*... shift). B is the full output for pb.
///   If Sa has roots <= those of B (in a suitable sense), this would follow.
///
///   Actually, let me think about what Sa and B represent combinatorially:
///   Sa(s) = #{source pi_a with swaps~ = s and pos(n-1) >= pa-1}
///   B_k   = #{sigma_b in Des(S) with pos(n) = pb and swaps(sigma_b) = k}
///
///   The LR ordering Sa >=_{LR} B means: Sa is "more concentrated to the left"
///   than B. In other words, the source permutations for pa with pos(n-1) >= pa-1
///   have fewer swaps (on average) than the output permutations for pb.
///
///   This makes intuitive sense because:
///   - Sa counts SOURCE permutations (n-1 elements), which have fewer swaps.
///   - B counts OUTPUT permutations (n elements), which have swaps >= source swaps.
///   - The eps2 shift in B pushes mass to the right.
///
///   More precisely: B_k = sum_q N_b(k - eps2(pb, q), q)
///   For q in L+M (zone L+M for pb): eps2(pb, q) = 1, so B_k gets N_b(k-1, q).
///   For q in R: eps2 = 0, so B_k gets N_b(k, q).
///
///   Total B: B_k = Tb(k-1) + Rb(k).
///   Total source b: Tb(s) + Rb(s) = total_b(s).
///   So B_k = total_b(k) - Rb(k) + Rb(k) + Tb(k-1) - Tb(k) + Tb(k)... hmm.
///   Simply: B_k = Tb(k-1) + Rb(k).
///
///   Sa(k) = Ma(k) + Ra(k).
///
///   Sa >=_{LR} B:
///   Sa(k1)/Sa(k2) >= B(k1)/B(k2) for k1 < k2.
///   [Ma(k1)+Ra(k1)] / [Ma(k2)+Ra(k2)] >= [Tb(k1-1)+Rb(k1)] / [Tb(k2-1)+Rb(k2)].
///
///   Since Sa is a "raw" polynomial (no shift) and B has a shift component (Tb at k-1),
///   the shift in B pushes B rightward, making B/Sa increase for large k.
///   This means Sa >=_{LR} B.
///
///   But this heuristic needs to be made rigorous!
///
///   KEY QUESTION: Does Sa >=_{LR} B follow from the inductive hypothesis
///   (that the source coefficient matrices are TN_2)?
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use std::collections::BTreeMap;

fn build_poly(vals: &[usize]) -> Vec<i64> {
    if vals.is_empty() { return vec![0]; }
    let max_s = *vals.iter().max().unwrap();
    let mut coeffs = vec![0i64; max_s + 1];
    for &s in vals { coeffs[s] += 1; }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 { coeffs.pop(); }
    coeffs
}

fn valid_positions(s_mask: u64, n: u8) -> Vec<u8> {
    let mut positions = Vec::new();
    for p in 1..n {
        if (s_mask >> (p - 1)) & 1 == 1 && (p < 2 || (s_mask >> (p - 2)) & 1 == 0) {
            positions.push(p);
        }
    }
    if n >= 2 && (s_mask >> (n - 2)) & 1 == 0 { positions.push(n); }
    positions
}

fn source_asc(s_mask: u64, p: u8, n: u8) -> u64 {
    if n <= 2 { return 0; }
    if p == n { return s_mask; }
    let mut sp = 0u64;
    if p == 1 {
        for j in 2..n { if (s_mask >> (j - 1)) & 1 == 1 { sp |= 1 << (j - 2); } }
    } else {
        for pos in 1..=(p.saturating_sub(2)) { if (s_mask >> (pos - 1)) & 1 == 1 { sp |= 1 << (pos - 1); } }
        for j in (p + 1)..n { if (s_mask >> (j - 1)) & 1 == 1 { sp |= 1 << (j - 2); } }
    }
    sp
}

fn source_desc(s_mask: u64, p: u8, n: u8) -> Option<u64> {
    if p <= 1 || p >= n { return None; }
    Some(source_asc(s_mask, p, n) | (1 << (p - 2)))
}

fn epsilon1(pi: &[u8], p: u8) -> bool {
    let n = pi.len() as u8 + 1;
    if p <= 1 || p >= n { return false; }
    pi[(p - 2) as usize] + 1 == pi[(p - 1) as usize]
}

fn eps2(p: u8, q: u8) -> usize {
    if p >= 2 && q <= p - 2 { 1 } else { 0 }
}

fn descent_set_to_string(mask: u64, n: u8) -> String {
    let mut s = String::from("{");
    let mut first = true;
    for i in 1..n {
        if (mask >> (i - 1)) & 1 == 1 {
            if !first { s.push(','); }
            s.push_str(&i.to_string());
            first = false;
        }
    }
    s.push('}');
    s
}

fn coeff(poly: &[i64], k: usize) -> i64 {
    if k < poly.len() { poly[k] } else { 0 }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum Zone { L, M, R }

fn classify_zone(q: u8, p_a: u8, p_b: u8) -> Zone {
    if p_a >= 2 && q <= p_a - 2 { Zone::L }
    else if p_b >= 2 && q <= p_b - 2 { Zone::M }
    else { Zone::R }
}

fn build_source_dist(
    s_mask: u64, p: u8, n: u8,
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

fn build_zone_raw(dist: &BTreeMap<(usize, u8), usize>, p_a: u8, p_b: u8, zone: Zone) -> Vec<i64> {
    let mut poly = vec![0i64; 20];
    for (&(s, q), &c) in dist {
        if classify_zone(q, p_a, p_b) == zone {
            if s >= poly.len() { poly.resize(s + 1, 0); }
            poly[s] += c as i64;
        }
    }
    while poly.len() > 1 && *poly.last().unwrap() == 0 { poly.pop(); }
    poly
}

fn poly_sum2(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &c) in a.iter().enumerate() { r[i] += c; }
    for (i, &c) in b.iter().enumerate() { r[i] += c; }
    while r.len() > 1 && *r.last().unwrap() == 0 { r.pop(); }
    r
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    println!("=== Two-term decomposition: detailed analysis ===\n");

    // Print Sa and B for small cases
    for n in 4u8..=max_n.min(7) {
        let perms_prev = all_permutations(n - 1);
        let mut by_descent_prev: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
        let mut perm_prev_list: Vec<&Vec<u8>> = Vec::new();
        for pi in &perms_prev {
            let idx = perm_prev_list.len();
            perm_prev_list.push(pi);
            by_descent_prev.entry(descent_set_bitmask(pi)).or_default().push(idx);
        }

        let perms_n = all_permutations(n);
        let mut by_descent_n: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_n {
            by_descent_n.entry(descent_set_bitmask(pi)).or_default().push(pi);
        }

        for s_mask in 0u64..(1 << (n - 1)) {
            if s_mask & 1 != 0 { continue; }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 { continue; }

            let mut source_dists: BTreeMap<u8, BTreeMap<(usize, u8), usize>> = BTreeMap::new();
            for &p in &vp {
                source_dists.insert(p,
                    build_source_dist(s_mask, p, n, &by_descent_prev, &perm_prev_list));
            }

            let mut lp_polys: BTreeMap<u8, Vec<i64>> = BTreeMap::new();
            if let Some(class) = by_descent_n.get(&s_mask) {
                for &p in &vp {
                    let vals: Vec<usize> = class.iter()
                        .filter(|pi| pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                        .map(|pi| compute(pi, Stat::Swaps))
                        .collect();
                    if !vals.is_empty() { lp_polys.insert(p, build_poly(&vals)); }
                }
            }

            for idx_pair in 0..vp.len() - 1 {
                let p_a = vp[idx_pair];
                let p_b = vp[idx_pair + 1];

                let dist_a = &source_dists[&p_a];
                let dist_b = &source_dists[&p_b];

                let la_raw = build_zone_raw(dist_a, p_a, p_b, Zone::L);
                let ma_raw = build_zone_raw(dist_a, p_a, p_b, Zone::M);
                let ra_raw = build_zone_raw(dist_a, p_a, p_b, Zone::R);
                let _lb_raw = build_zone_raw(dist_b, p_a, p_b, Zone::L);
                let _mb_raw = build_zone_raw(dist_b, p_a, p_b, Zone::M);
                let _rb_raw = build_zone_raw(dist_b, p_a, p_b, Zone::R);

                let a_poly_ref = lp_polys.get(&p_a);
                let b_poly_ref = lp_polys.get(&p_b);
                if a_poly_ref.is_none() || b_poly_ref.is_none() { continue; }
                let a_poly = a_poly_ref.unwrap();
                let b_poly = b_poly_ref.unwrap();

                let sa_raw = poly_sum2(&ma_raw, &ra_raw);

                // Print comparison
                println!("n={} S={} pa={} pb={}", n, descent_set_to_string(s_mask, n), p_a, p_b);
                println!("  A  = L^(pa) = {:?}", a_poly);
                println!("  B  = L^(pb) = {:?}", b_poly);
                println!("  Sa = Ma+Ra  = {:?}  (zone M+R raw for source A)", sa_raw);
                println!("  La = {:?}  (zone L raw for source A)", la_raw);
                println!("  Note: A_k = La(k-1) + Sa(k)");

                // Check: A_k = La(k-1) + Sa(k)
                let max_deg = a_poly.len().max(b_poly.len());
                for k in 0..max_deg {
                    let la_km1 = if k > 0 { coeff(&la_raw, k-1) } else { 0 };
                    let a_check = la_km1 + coeff(&sa_raw, k);
                    assert_eq!(a_check, coeff(a_poly, k), "A_k check failed");
                }

                // Show the LR ratios
                println!("  LR check Sa >=_LR B:");
                let mut ratios_ok = true;
                for k in 0..max_deg {
                    let sa_k = coeff(&sa_raw, k);
                    let b_k = coeff(b_poly, k);
                    if sa_k > 0 || b_k > 0 {
                        print!("    k={}: Sa={}, B={}", k, sa_k, b_k);
                        if sa_k > 0 && b_k > 0 {
                            print!("  ratio Sa/B = {:.4}", sa_k as f64 / b_k as f64);
                        }
                        println!();
                    }
                }

                // Check all minors
                for k1 in 0..max_deg {
                    for k2 in (k1+1)..max_deg {
                        let sa_k1 = coeff(&sa_raw, k1);
                        let sa_k2 = coeff(&sa_raw, k2);
                        let b_k1 = coeff(b_poly, k1);
                        let b_k2 = coeff(b_poly, k2);
                        let minor = sa_k1 * b_k2 - sa_k2 * b_k1;
                        if minor < 0 {
                            println!("  *** NEGATIVE: Sa({})B({}) - Sa({})B({}) = {} < 0!",
                                k1, k2, k2, k1, minor);
                            ratios_ok = false;
                        }
                    }
                }
                if ratios_ok {
                    println!("  Sa >=_LR B: VERIFIED");
                }
                println!();
            }
        }
    }

    // Deeper analysis: WHY does Sa >=_LR B hold?
    //
    // Sa(s) = sum_{q >= pa-1} N_a(s, q)  (zone M+R for source a, at raw degree s)
    // B_k = sum_q N_b(k - eps2(pb, q), q)
    //     = sum_{q <= pb-2} N_b(k-1, q) + sum_{q >= pb-1} N_b(k, q)
    //     = Tb(k-1) + Rb(k)
    //
    // The key observation:
    //   - Sa is a "partial row sum" of the source matrix for p_a, summing over q >= pa-1.
    //   - B is the FULL output polynomial for p_b, which includes a rightward shift.
    //   - Sa involves source permutations (n-1 elements) at raw swap count s.
    //   - B involves output permutations (n elements) at swap count k = s + eps.
    //   - The rightward shift in B (from eps2 applied to L+M zones) makes B
    //     stochastically larger (rightward) compared to raw source polynomials.
    //   - Sa is a raw source polynomial (no shift), hence lighter left tail.
    //   - Therefore Sa >=_LR B should hold because Sa is a raw source (no shift)
    //     and B is an output that has been shifted right.
    //
    // To make this precise:
    //   Total source b = Tb(s) + Rb(s)  (sum over all q of N_b(s,q))
    //   B_k = Tb(k-1) + Rb(k)
    //
    //   If Sa were equal to total_b, then:
    //     Sa(k1) * B(k2) - Sa(k2) * B(k1)
    //     = [Tb(k1)+Rb(k1)] * [Tb(k2-1)+Rb(k2)] - [Tb(k2)+Rb(k2)] * [Tb(k1-1)+Rb(k1)]
    //     = Tb(k1)*Tb(k2-1) - Tb(k2)*Tb(k1-1)  +  Tb(k1)*Rb(k2) - Tb(k2)*Rb(k1)
    //       + Rb(k1)*Tb(k2-1) - Rb(k2)*Tb(k1-1) + Rb(k1)*Rb(k2) - Rb(k2)*Rb(k1)
    //     = staircase(Tb, Tb) + standard(Tb, Rb) - standard(Rb, Tb) + staircase(Rb, Tb) + 0
    //     ... this gets messy.
    //
    // Better approach: think about it as a single polynomial comparison.
    //   Sa >=_LR B means Sa / B is non-increasing as a ratio.
    //   B has an extra +1 to certain terms (the t-shift), making B larger at higher k.
    //   So B_k / Sa(k) should be non-decreasing, which is what >=_LR means.
    //
    // This is really just saying: the output B has MORE mass at high k than the
    // source Sa, because B includes the eps2-shift bonus.
    //
    // Can we prove this formally?
    //
    // Lemma: If f(s) is a polynomial with nonneg coefficients and
    //   g(k) = sum_{i=0}^{m} c_i * f(k - d_i) where c_i > 0, d_i >= 0,
    //   and not all d_i = 0, then f >=_LR g (f has heavier left tail than g).
    //
    // Proof sketch: g is a "shifted version" of f (shifted right by positive amounts),
    //   so g has heavier right tail, hence f >=_LR g.
    //
    // For our setting:
    //   B_k = Tb(k-1) + Rb(k) = [lb + mb](k-1) + rb(k)
    //   Sa = ma + ra
    //
    // These are NOT directly "f" and "shifted f" because Sa and total_b are different
    // polynomials (different sources). So the lemma doesn't directly apply.
    //
    // However, if we could show:
    //   (i)  Sa >=_LR total_b  (Sa dominates total source b in LR)
    //   (ii) total_b >=_LR B   (source b dominates its own output in LR)
    // then by transitivity, Sa >=_LR B.
    //
    // (ii) is automatic from the eps2-shift: B = t*Tb + Rb, and total_b = Tb + Rb.
    //   So B - total_b = (t-1)*Tb, which shifts mass to the right. Hence total_b >=_LR B.
    //   (This requires Tb to have all nonneg coefficients, which it does.)
    //
    // Actually, (ii) needs more care. B_k = Tb(k-1) + Rb(k) and total_b(k) = Tb(k) + Rb(k).
    //   total_b(k) - B_k = Tb(k) - Tb(k-1) = backward difference of Tb.
    //   This is positive for small k (increasing part of Tb) and negative for large k.
    //   This means B has MORE mass at the tails (both left and right)? No, just at the right.
    //   Actually, sum_k B_k = sum_k [Tb(k-1) + Rb(k)] = sum_k Tb(k) + sum_k Rb(k) = sum_k total_b(k).
    //   Same total. And the shift of Tb to the right makes B stochastically larger.
    //   So total_b >=_LR B should hold.
    //
    // For (i), Sa >=_LR total_b means:
    //   Sa(k1)/Sa(k2) >= total_b(k1)/total_b(k2) for k1 < k2.
    //
    // Sa = Ma + Ra = sum_{q >= pa-1} N_a(s, q)
    // total_b = sum_q N_b(s, q)
    //
    // These involve different sources (a and b) summed over different q ranges.
    // This is the hard part.

    println!("{}", "=".repeat(80));
    println!("CHECKING INTERMEDIATE LR ORDERINGS");
    println!("{}", "=".repeat(80));
    println!();

    // Check: total_b >=_LR B (source dominates output in LR)
    let mut tb_lr_b_pass = 0u64;
    let mut tb_lr_b_fail = 0u64;

    // Check: Sa >=_LR total_b
    let mut sa_lr_totalb_pass = 0u64;
    let mut sa_lr_totalb_fail = 0u64;

    // Check: total_a >=_LR total_b
    let mut totala_lr_totalb_pass = 0u64;
    let mut totala_lr_totalb_fail = 0u64;

    // Check: Sa >=_LR Sb (zone M+R of a vs zone M+R of b)
    let mut sa_lr_sb_pass = 0u64;
    let mut sa_lr_sb_fail = 0u64;

    // Check: La >=_LR Lb (zone L of a vs zone L of b)
    let mut la_lr_lb_pass = 0u64;
    let mut la_lr_lb_fail = 0u64;

    for n in 4..=max_n {
        let perms_prev = all_permutations(n - 1);
        let mut by_descent_prev: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
        let mut perm_prev_list: Vec<&Vec<u8>> = Vec::new();
        for pi in &perms_prev {
            let idx = perm_prev_list.len();
            perm_prev_list.push(pi);
            by_descent_prev.entry(descent_set_bitmask(pi)).or_default().push(idx);
        }

        let perms_n = all_permutations(n);
        let mut by_descent_n: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_n {
            by_descent_n.entry(descent_set_bitmask(pi)).or_default().push(pi);
        }

        for s_mask in 0u64..(1 << (n - 1)) {
            if s_mask & 1 != 0 { continue; }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 { continue; }

            let mut source_dists: BTreeMap<u8, BTreeMap<(usize, u8), usize>> = BTreeMap::new();
            for &p in &vp {
                source_dists.insert(p,
                    build_source_dist(s_mask, p, n, &by_descent_prev, &perm_prev_list));
            }

            let mut lp_polys: BTreeMap<u8, Vec<i64>> = BTreeMap::new();
            if let Some(class) = by_descent_n.get(&s_mask) {
                for &p in &vp {
                    let vals: Vec<usize> = class.iter()
                        .filter(|pi| pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                        .map(|pi| compute(pi, Stat::Swaps))
                        .collect();
                    if !vals.is_empty() { lp_polys.insert(p, build_poly(&vals)); }
                }
            }

            for idx_pair in 0..vp.len() - 1 {
                let p_a = vp[idx_pair];
                let p_b = vp[idx_pair + 1];

                let dist_a = &source_dists[&p_a];
                let dist_b = &source_dists[&p_b];

                let la_raw = build_zone_raw(dist_a, p_a, p_b, Zone::L);
                let ma_raw = build_zone_raw(dist_a, p_a, p_b, Zone::M);
                let ra_raw = build_zone_raw(dist_a, p_a, p_b, Zone::R);
                let lb_raw = build_zone_raw(dist_b, p_a, p_b, Zone::L);
                let mb_raw = build_zone_raw(dist_b, p_a, p_b, Zone::M);
                let rb_raw = build_zone_raw(dist_b, p_a, p_b, Zone::R);

                let b_poly_ref = lp_polys.get(&p_b);
                if b_poly_ref.is_none() { continue; }
                let b_poly = b_poly_ref.unwrap();

                let sa_raw = poly_sum2(&ma_raw, &ra_raw);
                let sb_raw = poly_sum2(&mb_raw, &rb_raw);
                let total_a = poly_sum2(&la_raw, &sa_raw);
                let total_b = poly_sum2(&lb_raw, &sb_raw);

                let max_s = [&sa_raw, &total_b, b_poly, &total_a, &sb_raw, &la_raw, &lb_raw]
                    .iter().map(|p| p.len()).max().unwrap_or(1);

                // Check total_b >=_LR B
                let mut ok = true;
                for k1 in 0..max_s {
                    for k2 in (k1+1)..max_s {
                        let tb_k1 = coeff(&total_b, k1);
                        let tb_k2 = coeff(&total_b, k2);
                        let b_k1 = coeff(b_poly, k1);
                        let b_k2 = coeff(b_poly, k2);
                        if tb_k1 * b_k2 < tb_k2 * b_k1 {
                            ok = false;
                        }
                    }
                }
                if ok { tb_lr_b_pass += 1; } else { tb_lr_b_fail += 1; }

                // Check Sa >=_LR total_b
                ok = true;
                for k1 in 0..max_s {
                    for k2 in (k1+1)..max_s {
                        let sa_k1 = coeff(&sa_raw, k1);
                        let sa_k2 = coeff(&sa_raw, k2);
                        let tb_k1 = coeff(&total_b, k1);
                        let tb_k2 = coeff(&total_b, k2);
                        if sa_k1 * tb_k2 < sa_k2 * tb_k1 {
                            ok = false;
                        }
                    }
                }
                if ok { sa_lr_totalb_pass += 1; } else { sa_lr_totalb_fail += 1; }

                // Check total_a >=_LR total_b
                ok = true;
                for k1 in 0..max_s {
                    for k2 in (k1+1)..max_s {
                        if coeff(&total_a, k1) * coeff(&total_b, k2) < coeff(&total_a, k2) * coeff(&total_b, k1) {
                            ok = false;
                        }
                    }
                }
                if ok { totala_lr_totalb_pass += 1; } else { totala_lr_totalb_fail += 1; }

                // Check Sa >=_LR Sb
                ok = true;
                for k1 in 0..max_s {
                    for k2 in (k1+1)..max_s {
                        if coeff(&sa_raw, k1) * coeff(&sb_raw, k2) < coeff(&sa_raw, k2) * coeff(&sb_raw, k1) {
                            ok = false;
                        }
                    }
                }
                if ok { sa_lr_sb_pass += 1; } else { sa_lr_sb_fail += 1; }

                // Check La >=_LR Lb
                ok = true;
                for k1 in 0..max_s {
                    for k2 in (k1+1)..max_s {
                        if coeff(&la_raw, k1) * coeff(&lb_raw, k2) < coeff(&la_raw, k2) * coeff(&lb_raw, k1) {
                            ok = false;
                        }
                    }
                }
                if ok { la_lr_lb_pass += 1; } else { la_lr_lb_fail += 1; }
            }
        }
    }

    let total_pairs = tb_lr_b_pass + tb_lr_b_fail;
    println!("{:<30} {:>6}/{:<6} {}", "total_b >=_LR B", tb_lr_b_pass, total_pairs,
        if tb_lr_b_fail == 0 { "ALWAYS" } else { "FAILS" });
    println!("{:<30} {:>6}/{:<6} {}", "Sa >=_LR total_b", sa_lr_totalb_pass, total_pairs,
        if sa_lr_totalb_fail == 0 { "ALWAYS" } else { "FAILS" });
    println!("{:<30} {:>6}/{:<6} {}", "total_a >=_LR total_b", totala_lr_totalb_pass, total_pairs,
        if totala_lr_totalb_fail == 0 { "ALWAYS" } else { "FAILS" });
    println!("{:<30} {:>6}/{:<6} {}", "Sa >=_LR Sb", sa_lr_sb_pass, total_pairs,
        if sa_lr_sb_fail == 0 { "ALWAYS" } else { "FAILS" });
    println!("{:<30} {:>6}/{:<6} {}", "La >=_LR Lb", la_lr_lb_pass, total_pairs,
        if la_lr_lb_fail == 0 { "ALWAYS" } else { "FAILS" });

    println!();
    if tb_lr_b_fail == 0 {
        println!("total_b >=_LR B: PROVED by the eps2-shift structure.");
        println!("  B = t*Tb + Rb is a right-shifted version of total_b = Tb + Rb.");
    }
    if sa_lr_totalb_fail == 0 {
        println!("Sa >=_LR total_b: HOLDS! This is the key non-trivial claim.");
        println!("  Combined with total_b >=_LR B, this gives Sa >=_LR B by transitivity.");
    } else {
        println!("Sa >=_LR total_b: FAILS in {} cases. Cannot use transitivity.", sa_lr_totalb_fail);
        if totala_lr_totalb_fail == 0 {
            println!("But total_a >=_LR total_b holds, which is a weaker statement.");
        }
    }
}
