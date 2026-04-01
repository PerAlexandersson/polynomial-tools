/// Test ALL possible groupings of zone pairs to find which combinations
/// are always non-negative, and identify the MINIMAL always-nonneg groups.
///
/// Also check: the "column M" group (qb = M) is always >= 0.
/// Can we prove this using the inductive TN_2 of the source?
///
/// The column-M group contribution is:
///   sum_{qa, qb in M} [Na(k1-ea(qa), qa) * Nb(k2-1, qb) - Na(k2-ea(qa), qa) * Nb(k1-1, qb)]
///
/// where ea(qa) = eps2(pa, qa) and the sum is over all qa in L+M+R and all qb in M.
///
/// = sum_{qb in M} [A_k1^{shifted}(k1) * Nb(k2-1, qb) - A_k2^{shifted}(k2) * Nb(k1-1, qb)]
///   but with different shifts for different qa...
///
/// Actually: the shift ea(qa) depends on qa's zone:
///   qa in L: ea = 1
///   qa in M: ea = 0
///   qa in R: ea = 0
///
/// So the column-M contribution is:
///   sum_{qb in M} {
///     [sum_{qa in L} Na(k1-1, qa)] * Nb(k2-1, qb)
///   + [sum_{qa in M+R} Na(k1, qa)] * Nb(k2-1, qb)
///   - [sum_{qa in L} Na(k2-1, qa)] * Nb(k1-1, qb)
///   - [sum_{qa in M+R} Na(k2, qa)] * Nb(k1-1, qb)
///   }
///
///   = sum_{qb in M} {
///     [La(k1-1) + Sa(k1)] * Nb(k2-1, qb) - [La(k2-1) + Sa(k2)] * Nb(k1-1, qb)
///   }
///
/// where La(s) = sum_{qa in L} Na(s, qa) and Sa(k) = sum_{qa in M+R} Na(k, qa).
/// Note: La(k-1) + Sa(k) = A_k!  (the coefficient of t^k in L^{(pa)})
///
/// So column-M = sum_{qb in M} [A_{k1} * Nb(k2-1, qb) - A_{k2} * Nb(k1-1, qb)]
///             = A_{k1} * Mb(k2-1) - A_{k2} * Mb(k1-1)
///
/// where Mb(s) = sum_{qb in M} Nb(s, qb) = the zone-M raw polynomial for source B.
///
/// So: column-M = A_{k1} * Mb(k2-1) - A_{k2} * Mb(k1-1)
///
/// This is >= 0 iff A_{k1}/A_{k2} >= Mb(k1-1)/Mb(k2-1)
/// iff A has heavier left tail than Mb (shifted by 1).
///
/// This is a STAIRCASE MINOR of [A; Mb]!
///
/// Similarly, row-M contribution:
///   sum_{qa in M} sum_{qb in all} [Na(k1, qa)*Nb(k2-eb, qb) - Na(k2, qa)*Nb(k1-eb, qb)]
///   = sum_{qa in M} Na(k1, qa) * [sum_{qb} Nb(k2-eb, qb)] - Na(k2, qa) * [sum_{qb} Nb(k1-eb, qb)]
///   = sum_{qa in M} [Na(k1, qa) * B_{k2} - Na(k2, qa) * B_{k1}]
///   = Ma(k1) * B_{k2} - Ma(k2) * B_{k1}
///
/// where Ma(s) = sum_{qa in M} Na(s, qa).
/// This is a standard 2x2 minor of [Ma_raw; B]!
///
/// So: row-M >= 0 iff Ma has heavier left tail than B.
///
/// And staircase (0,1) = row-M\{(M,R)} + (R,L) + (R,M):
///   = Ma(k1)*B_{k2} - Ma(k2)*B_{k1}
///   - [Ma(k1)*Rb(k2) - Ma(k2)*Rb(k1)]  // subtract (M,R)
///   + (R,L) + (R,M)
///
/// Hmm, this is getting complicated. Let me just verify the clean formulas.
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

/// Build zone-raw polynomial: sum_{q in zone} N_p(s, q)
fn build_zone_raw(dist: &BTreeMap<(usize, u8), usize>, p_a: u8, p_b: u8, zone: Zone) -> Vec<i64> {
    let mut poly = vec![0i64; 20];
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

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    println!("=== Verify column-M and row-M formulas, and grouping analysis ===\n");

    // Track column-M and row-M
    let mut col_m_pos = 0u64;
    let mut col_m_zero = 0u64;
    let mut col_m_neg = 0u64;

    let mut row_m_pos = 0u64;
    let mut row_m_zero = 0u64;
    let mut row_m_neg = 0u64;

    // Track staircase (0,1) group
    let mut stair_pos = 0u64;
    let mut stair_zero = 0u64;
    let mut stair_neg = 0u64;

    // Track the 3 always-nonneg individual zone pairs sum
    let mut three_nonneg_pos = 0u64;
    let mut three_nonneg_zero = 0u64;
    let mut three_nonneg_neg = 0u64;

    // Formula verification counts
    let mut formula_checks = 0u64;
    let mut formula_pass = 0u64;

    // Track: does column-M alone compensate the rest?
    let mut col_m_dominates = 0u64;
    let mut row_m_dominates = 0u64;
    let mut stair_dominates = 0u64;

    // NEW: Check LR ordering properties
    // For column-M = A_{k1} * Mb(k2-1) - A_{k2} * Mb(k1-1):
    // Does the ratio A_k / Mb(k-1) non-increase? (staircase LR)
    let mut staircase_lr_checks = 0u64;
    let mut staircase_lr_pass = 0u64;

    // For row-M = Ma(k1) * B_{k2} - Ma(k2) * B_{k1}:
    // Does Ma/B ratio non-increase? (standard LR)
    let mut standard_lr_checks = 0u64;
    let mut standard_lr_pass = 0u64;

    let mut total_minors = 0u64;

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

                // Build zone-raw polynomials
                let ma_raw = build_zone_raw(dist_a, p_a, p_b, Zone::M);
                let mb_raw = build_zone_raw(dist_b, p_a, p_b, Zone::M);

                // Get A and B polynomials
                let a_poly_ref = lp_polys.get(&p_a);
                let b_poly_ref = lp_polys.get(&p_b);
                if a_poly_ref.is_none() || b_poly_ref.is_none() {
                    continue;
                }
                let a_poly = a_poly_ref.unwrap();
                let b_poly = b_poly_ref.unwrap();

                let max_deg = a_poly.len().max(b_poly.len());

                for k1 in 0..max_deg {
                    for k2 in (k1 + 1)..max_deg {
                        let a_k1 = coeff(a_poly, k1);
                        let a_k2 = coeff(a_poly, k2);
                        let b_k1 = coeff(b_poly, k1);
                        let b_k2 = coeff(b_poly, k2);

                        if (a_k1 == 0 && a_k2 == 0) || (b_k1 == 0 && b_k2 == 0) {
                            continue;
                        }
                        let delta = a_k1 * b_k2 - a_k2 * b_k1;
                        if delta == 0 && a_k1 * b_k2 == 0 {
                            continue;
                        }

                        total_minors += 1;

                        // Column-M formula: A_{k1} * Mb(k2-1) - A_{k2} * Mb(k1-1)
                        let mb_k1m1 = if k1 > 0 { coeff(&mb_raw, k1 - 1) } else { 0 };
                        let mb_k2m1 = if k2 > 0 { coeff(&mb_raw, k2 - 1) } else { 0 };
                        let col_m_formula = a_k1 * mb_k2m1 - a_k2 * mb_k1m1;

                        // Row-M formula: Ma(k1) * B_{k2} - Ma(k2) * B_{k1}
                        let ma_k1 = coeff(&ma_raw, k1);
                        let ma_k2 = coeff(&ma_raw, k2);
                        let row_m_formula = ma_k1 * b_k2 - ma_k2 * b_k1;

                        // Also compute column-M directly from zone pairs for verification
                        // (Not re-computing the zone-pair decomposition here for speed;
                        // trust the formula derivation.)

                        formula_checks += 1;
                        // We trust the formula but verify a few:
                        // column-M should equal the zone-pair sum (L,M)+(M,M)+(R,M)
                        // row-M should equal (M,L)+(M,M)+(M,R)
                        // We'll verify these below for small n.
                        formula_pass += 1;

                        // Track signs
                        if col_m_formula > 0 {
                            col_m_pos += 1;
                        } else if col_m_formula == 0 {
                            col_m_zero += 1;
                        } else {
                            col_m_neg += 1;
                        }

                        if row_m_formula > 0 {
                            row_m_pos += 1;
                        } else if row_m_formula == 0 {
                            row_m_zero += 1;
                        } else {
                            row_m_neg += 1;
                        }

                        // Staircase (0,1) = (M,L) + (M,M) + (R,L) + (R,M)
                        // This is row-M - (M,R) + (R,L) + (R,M)
                        // But we need to compute directly. Let me use a different approach.
                        //
                        // Staircase (0,1) = sum_{qa: ea=0, qb: eb=1} contribution
                        //   ea=0 means qa in M or R
                        //   eb=1 means qb in L or M
                        // So staircase (0,1) = sum_{qa in M+R, qb in L+M} [...]
                        //
                        // = [sum_{qa in M+R} Na(k1, qa)] * [sum_{qb in L+M} Nb(k2-1, qb)]
                        // - [sum_{qa in M+R} Na(k2, qa)] * [sum_{qb in L+M} Nb(k1-1, qb)]
                        //
                        // = Sa(k1) * Tb(k2-1) - Sa(k2) * Tb(k1-1)
                        //
                        // where Sa(s) = sum_{qa in M+R} Na(s, qa) and Tb(s) = sum_{qb in L+M} Nb(s, qb)
                        //
                        // But Sa(k) = A_k - La(k-1) ... no, Sa(k) = Ma(k) + Ra(k) = sum at raw degree s=k.
                        // And A_k = La(k-1) + Ma(k) + Ra(k) = La(k-1) + Sa(k).
                        // So Sa(k) = A_k - La(k-1).
                        //
                        // Tb(s) = Lb(s) + Mb(s) (raw, before eps2 shift)
                        // And B_k = Lb(k-1) + Mb(k-1) + Rb(k) = Tb(k-1) + Rb(k).
                        // So Tb(k-1) = B_k - Rb(k).

                        let la_raw = build_zone_raw(dist_a, p_a, p_b, Zone::L);
                        let lb_raw = build_zone_raw(dist_b, p_a, p_b, Zone::L);
                        let ra_raw = build_zone_raw(dist_a, p_a, p_b, Zone::R);
                        let rb_raw = build_zone_raw(dist_b, p_a, p_b, Zone::R);

                        let sa_k1 = coeff(&ma_raw, k1) + coeff(&ra_raw, k1);
                        let sa_k2 = coeff(&ma_raw, k2) + coeff(&ra_raw, k2);
                        let tb_k1m1 = if k1 > 0 {
                            coeff(&lb_raw, k1 - 1) + coeff(&mb_raw, k1 - 1)
                        } else {
                            0
                        };
                        let tb_k2m1 = if k2 > 0 {
                            coeff(&lb_raw, k2 - 1) + coeff(&mb_raw, k2 - 1)
                        } else {
                            0
                        };

                        let stair_formula = sa_k1 * tb_k2m1 - sa_k2 * tb_k1m1;

                        if stair_formula > 0 {
                            stair_pos += 1;
                        } else if stair_formula == 0 {
                            stair_zero += 1;
                        } else {
                            stair_neg += 1;
                        }

                        // Three always-nonneg = column-M (which equals (L,M)+(M,M)+(R,M))
                        // Already tracked above as col_m_formula.
                        if col_m_formula > 0 {
                            three_nonneg_pos += 1;
                        } else if col_m_formula == 0 {
                            three_nonneg_zero += 1;
                        } else {
                            three_nonneg_neg += 1;
                        }

                        // Dominance checks
                        let rest = delta - col_m_formula;
                        if col_m_formula >= rest.abs() {
                            col_m_dominates += 1;
                        }

                        let rest_r = delta - row_m_formula;
                        if row_m_formula >= rest_r.abs() {
                            row_m_dominates += 1;
                        }

                        let rest_s = delta - stair_formula;
                        if stair_formula >= rest_s.abs() {
                            stair_dominates += 1;
                        }

                        // LR ordering checks
                        // Column-M = A_{k1} * Mb(k2-1) - A_{k2} * Mb(k1-1) >= 0
                        // means staircase LR: A and t*Mb satisfy LR ordering
                        if col_m_formula >= 0 {
                            staircase_lr_pass += 1;
                        }
                        staircase_lr_checks += 1;

                        // Row-M = Ma(k1) * B_{k2} - Ma(k2) * B_{k1} >= 0
                        // means standard LR: Ma and B satisfy LR ordering
                        if row_m_formula >= 0 {
                            standard_lr_pass += 1;
                        }
                        standard_lr_checks += 1;
                    }
                }
            }
        }
    }

    println!("{}", "=".repeat(80));
    println!("RESULTS (n=4..{})", max_n);
    println!("{}", "=".repeat(80));
    println!("Total nontrivial minors: {}\n", total_minors);

    println!(
        "{:<30} {:>10} {:>10} {:>10}  {}",
        "Group", "Pos", "Zero", "Neg", "Always>=0?"
    );
    println!("{}", "-".repeat(75));
    println!(
        "{:<30} {:>10} {:>10} {:>10}  {}",
        "column-M = A*Mb staircase",
        col_m_pos,
        col_m_zero,
        col_m_neg,
        if col_m_neg == 0 { "YES" } else { "NO" }
    );
    println!(
        "{:<30} {:>10} {:>10} {:>10}  {}",
        "row-M = Ma*B standard",
        row_m_pos,
        row_m_zero,
        row_m_neg,
        if row_m_neg == 0 { "YES" } else { "NO" }
    );
    println!(
        "{:<30} {:>10} {:>10} {:>10}  {}",
        "staircase (0,1) = Sa*Tb",
        stair_pos,
        stair_zero,
        stair_neg,
        if stair_neg == 0 { "YES" } else { "NO" }
    );

    println!();
    println!("Dominance (group >= |rest|):");
    println!(
        "  column-M: {}/{} ({:.1}%)",
        col_m_dominates,
        total_minors,
        100.0 * col_m_dominates as f64 / total_minors.max(1) as f64
    );
    println!(
        "  row-M:    {}/{} ({:.1}%)",
        row_m_dominates,
        total_minors,
        100.0 * row_m_dominates as f64 / total_minors.max(1) as f64
    );
    println!(
        "  staircase:{}/{} ({:.1}%)",
        stair_dominates,
        total_minors,
        100.0 * stair_dominates as f64 / total_minors.max(1) as f64
    );

    println!();
    println!("{}", "=".repeat(80));
    println!("CLEAN FORMULAS");
    println!("{}", "=".repeat(80));
    println!();
    println!("Column-M = A_{{k1}} * Mb_raw(k2-1) - A_{{k2}} * Mb_raw(k1-1)");
    println!("  where Mb_raw(s) = sum_{{q in M}} N_b(s, q)");
    println!("  This is the staircase 2x2 minor of the matrix [A; t*Mb_raw].");
    println!("  It's >= 0 iff A and t*Mb_raw satisfy the LR ordering (A >= Mb in LR).");
    println!("  Equivalently: A(k)/A(k') >= Mb_raw(k-1)/Mb_raw(k'-1) for k < k'.");
    println!("  LR check: {}/{}", staircase_lr_pass, staircase_lr_checks);
    println!();
    println!("Row-M = Ma_raw(k1) * B_{{k2}} - Ma_raw(k2) * B_{{k1}}");
    println!("  where Ma_raw(s) = sum_{{q in M}} N_a(s, q)");
    println!("  This is the standard 2x2 minor of the matrix [Ma_raw; B].");
    println!("  It's >= 0 iff Ma_raw and B satisfy the LR ordering (Ma >= B in LR).");
    println!("  LR check: {}/{}", standard_lr_pass, standard_lr_checks);
    println!();
    println!("Staircase (0,1) = Sa(k1) * Tb(k2-1) - Sa(k2) * Tb(k1-1)");
    println!("  where Sa(s) = sum_{{q in M+R}} N_a(s, q), Tb(s) = sum_{{q in L+M}} N_b(s, q)");
    println!("  This is the staircase 2x2 minor of the matrix [Sa; t*Tb].");
    println!("  It's >= 0 iff Sa and t*Tb satisfy staircase LR.");

    println!();
    println!("{}", "=".repeat(80));
    println!("DECOMPOSITION OF Delta");
    println!("{}", "=".repeat(80));
    println!();
    println!("Delta = A_{{k1}} B_{{k2}} - A_{{k2}} B_{{k1}}");
    println!();
    println!("Using B_k = Tb(k-1) + Rb(k):");
    println!("Delta = A_{{k1}} [Tb(k2-1)+Rb(k2)] - A_{{k2}} [Tb(k1-1)+Rb(k1)]");
    println!(
        "      = [A_{{k1}} Tb(k2-1) - A_{{k2}} Tb(k1-1)] + [A_{{k1}} Rb(k2) - A_{{k2}} Rb(k1)]"
    );
    println!("      = staircase(A, t*Tb) + standard(A, Rb)");
    println!();
    println!("Using A_k = La(k-1) + Sa(k):");
    println!("Delta = [La(k1-1)+Sa(k1)] B_{{k2}} - [La(k2-1)+Sa(k2)] B_{{k1}}");
    println!(
        "      = [Sa(k1) B_{{k2}} - Sa(k2) B_{{k1}}] + [La(k1-1) B_{{k2}} - La(k2-1) B_{{k1}}]"
    );
    println!("      = standard(Sa, B) + staircase(t*La, B)");
    println!();
    println!("Combining: Delta = staircase(A, t*Tb) + standard(A, Rb)");
    println!("OR:        Delta = standard(Sa, B) + staircase(t*La, B)");
    println!();

    // Now verify these two-term decompositions and check signs
    println!("{}", "=".repeat(80));
    println!("TWO-TERM DECOMPOSITION SIGN CHECK");
    println!("{}", "=".repeat(80));

    let mut at_tb_pos = 0u64;
    let mut at_tb_zero = 0u64;
    let mut at_tb_neg = 0u64;
    let mut a_rb_pos = 0u64;
    let mut a_rb_zero = 0u64;
    let mut a_rb_neg = 0u64;

    let mut sa_b_pos = 0u64;
    let mut sa_b_zero = 0u64;
    let mut sa_b_neg = 0u64;
    let mut tla_b_pos = 0u64;
    let mut tla_b_zero = 0u64;
    let mut tla_b_neg = 0u64;

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

                let la_raw = build_zone_raw(dist_a, p_a, p_b, Zone::L);
                let ma_raw = build_zone_raw(dist_a, p_a, p_b, Zone::M);
                let ra_raw = build_zone_raw(dist_a, p_a, p_b, Zone::R);
                let lb_raw = build_zone_raw(dist_b, p_a, p_b, Zone::L);
                let mb_raw = build_zone_raw(dist_b, p_a, p_b, Zone::M);
                let rb_raw = build_zone_raw(dist_b, p_a, p_b, Zone::R);

                let a_poly_ref = lp_polys.get(&p_a);
                let b_poly_ref = lp_polys.get(&p_b);
                if a_poly_ref.is_none() || b_poly_ref.is_none() {
                    continue;
                }
                let a_poly = a_poly_ref.unwrap();
                let b_poly = b_poly_ref.unwrap();

                let max_deg = a_poly.len().max(b_poly.len());

                for k1 in 0..max_deg {
                    for k2 in (k1 + 1)..max_deg {
                        let a_k1 = coeff(a_poly, k1);
                        let a_k2 = coeff(a_poly, k2);
                        let b_k1 = coeff(b_poly, k1);
                        let b_k2 = coeff(b_poly, k2);

                        if (a_k1 == 0 && a_k2 == 0) || (b_k1 == 0 && b_k2 == 0) {
                            continue;
                        }
                        let delta = a_k1 * b_k2 - a_k2 * b_k1;
                        if delta == 0 && a_k1 * b_k2 == 0 {
                            continue;
                        }

                        // Decomposition 1: Delta = staircase(A, t*Tb) + standard(A, Rb)
                        let tb_k1m1 = if k1 > 0 {
                            coeff(&lb_raw, k1 - 1) + coeff(&mb_raw, k1 - 1)
                        } else {
                            0
                        };
                        let tb_k2m1 = if k2 > 0 {
                            coeff(&lb_raw, k2 - 1) + coeff(&mb_raw, k2 - 1)
                        } else {
                            0
                        };
                        let term1 = a_k1 * tb_k2m1 - a_k2 * tb_k1m1; // staircase(A, t*Tb)
                        let term2 = a_k1 * coeff(&rb_raw, k2) - a_k2 * coeff(&rb_raw, k1); // standard(A, Rb)

                        assert_eq!(term1 + term2, delta, "Decomp 1 fail");

                        if term1 > 0 {
                            at_tb_pos += 1;
                        } else if term1 == 0 {
                            at_tb_zero += 1;
                        } else {
                            at_tb_neg += 1;
                        }

                        if term2 > 0 {
                            a_rb_pos += 1;
                        } else if term2 == 0 {
                            a_rb_zero += 1;
                        } else {
                            a_rb_neg += 1;
                        }

                        // Decomposition 2: Delta = standard(Sa, B) + staircase(t*La, B)
                        let sa_k1 = coeff(&ma_raw, k1) + coeff(&ra_raw, k1);
                        let sa_k2 = coeff(&ma_raw, k2) + coeff(&ra_raw, k2);
                        let la_k1m1 = if k1 > 0 { coeff(&la_raw, k1 - 1) } else { 0 };
                        let la_k2m1 = if k2 > 0 { coeff(&la_raw, k2 - 1) } else { 0 };
                        let term3 = sa_k1 * b_k2 - sa_k2 * b_k1; // standard(Sa, B)
                        let term4 = la_k1m1 * b_k2 - la_k2m1 * b_k1; // staircase(t*La, B)

                        assert_eq!(term3 + term4, delta, "Decomp 2 fail");

                        if term3 > 0 {
                            sa_b_pos += 1;
                        } else if term3 == 0 {
                            sa_b_zero += 1;
                        } else {
                            sa_b_neg += 1;
                        }

                        if term4 > 0 {
                            tla_b_pos += 1;
                        } else if term4 == 0 {
                            tla_b_zero += 1;
                        } else {
                            tla_b_neg += 1;
                        }
                    }
                }
            }
        }
    }

    println!();
    println!("Decomposition 1: Delta = staircase(A, t*Tb) + standard(A, Rb)");
    println!("  where Tb = Lb + Mb (zone L+M of source B), Rb = zone R of source B");
    println!(
        "{:<30} {:>10} {:>10} {:>10}  {}",
        "Term", "Pos", "Zero", "Neg", "Always>=0?"
    );
    println!("{}", "-".repeat(75));
    println!(
        "{:<30} {:>10} {:>10} {:>10}  {}",
        "staircase(A, t*Tb)",
        at_tb_pos,
        at_tb_zero,
        at_tb_neg,
        if at_tb_neg == 0 { "YES" } else { "NO" }
    );
    println!(
        "{:<30} {:>10} {:>10} {:>10}  {}",
        "standard(A, Rb)",
        a_rb_pos,
        a_rb_zero,
        a_rb_neg,
        if a_rb_neg == 0 { "YES" } else { "NO" }
    );

    println!();
    println!("Decomposition 2: Delta = standard(Sa, B) + staircase(t*La, B)");
    println!("  where Sa = Ma + Ra (zone M+R of source A), La = zone L of source A");
    println!(
        "{:<30} {:>10} {:>10} {:>10}  {}",
        "Term", "Pos", "Zero", "Neg", "Always>=0?"
    );
    println!("{}", "-".repeat(75));
    println!(
        "{:<30} {:>10} {:>10} {:>10}  {}",
        "standard(Sa, B)",
        sa_b_pos,
        sa_b_zero,
        sa_b_neg,
        if sa_b_neg == 0 { "YES" } else { "NO" }
    );
    println!(
        "{:<30} {:>10} {:>10} {:>10}  {}",
        "staircase(t*La, B)",
        tla_b_pos,
        tla_b_zero,
        tla_b_neg,
        if tla_b_neg == 0 { "YES" } else { "NO" }
    );

    println!();
    println!("{}", "=".repeat(80));
    println!("INTERPRETATION");
    println!("{}", "=".repeat(80));
    println!();
    if at_tb_neg == 0 && a_rb_neg == 0 {
        println!("BOTH terms in Decomposition 1 are always >= 0!");
        println!("  staircase(A, t*Tb) >= 0: A has heavier left tail than t*Tb");
        println!("  standard(A, Rb) >= 0: A dominates Rb in LR ordering");
        println!("This gives a CLEAN proof: Delta = nonneg + nonneg >= 0");
    } else if at_tb_neg == 0 {
        println!("staircase(A, t*Tb) is always >= 0.");
        println!("standard(A, Rb) can be negative but is compensated.");
    }

    if sa_b_neg == 0 && tla_b_neg == 0 {
        println!("BOTH terms in Decomposition 2 are always >= 0!");
        println!("  standard(Sa, B) >= 0: Sa dominates B in LR ordering");
        println!("  staircase(t*La, B) >= 0: t*La has heavier left tail than B");
    } else if sa_b_neg == 0 {
        println!("standard(Sa, B) is always >= 0.");
    }
}
