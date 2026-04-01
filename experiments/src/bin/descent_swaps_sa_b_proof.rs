/// Investigate WHY Sa >=_LR B always holds.
///
/// Sa(s) = sum_{q >= pa-1} N_a(s, q)
/// B_k = L^{(pb)}_k = sum_q N_b(k - eps2(pb, q), q)
///
/// Expand both using the insertion structure:
///   N_p(s, q) counts source permutations pi in source(p) with:
///     - modified swaps = swaps(pi) + eps1(pi, p) = s
///     - pos(n-1, pi) = q
///
/// For consecutive valid positions pa < pb:
///   source(pa) comes from descent sets S'_pa (ascent) and S''_pa (descent)
///   source(pb) comes from descent sets S'_pb (ascent) and S''_pb (descent)
///
/// The sources overlap but differ at position pa-1 (for S'_pa/S'_pb difference)
/// and possibly at position pb-1.
///
/// KEY STRUCTURAL OBSERVATION:
/// Sa sums over q >= pa-1. This includes ALL q values except q <= pa-2.
/// In the source for pa: q = pos(n-1) in pi. The eps2 threshold for pa is pa-2.
/// So Sa sums over exactly the q values where eps2(pa, q) = 0.
/// These are the "unshifted" contributions to A.
/// In other words: A_k = La(k-1) + Sa(k), where La(k-1) is the "shifted" part.
///
/// Sa(k) = A_k - La(k-1).
///
/// So Sa >=_LR B means: A - t*La >=_LR B (as formal power series / polynomials).
///
/// Alternatively: Sa(k1)*B(k2) >= Sa(k2)*B(k1) for k1 < k2.
///   [A(k1) - La(k1-1)] * B(k2) >= [A(k2) - La(k2-1)] * B(k1)
///   A(k1)*B(k2) - A(k2)*B(k1) >= La(k1-1)*B(k2) - La(k2-1)*B(k1)
///   Delta >= staircase(t*La, B)
///
/// So Sa >=_LR B is EQUIVALENT to Delta >= staircase(t*La, B).
/// This is Decomposition 2: Delta = standard(Sa, B) + staircase(t*La, B).
///
/// And we know staircase(t*La, B) can be negative. So the claim is that
/// Delta (which is always >= 0) is at least as large as staircase(t*La, B).
///
/// This is a STRONGER version of Delta >= 0!
///
/// Similarly, staircase(A, t*Tb) >= 0 means:
///   A(k1)*Tb(k2-1) >= A(k2)*Tb(k1-1)
///   Using B_k = Tb(k-1) + Rb(k):
///   A(k1)*[B(k2) - Rb(k2)] >= A(k2)*[B(k1) - Rb(k1)]
///   A(k1)*B(k2) - A(k2)*B(k1) >= A(k1)*Rb(k2) - A(k2)*Rb(k1)
///   Delta >= standard(A, Rb)
///
/// So staircase(A, t*Tb) >= 0 is EQUIVALENT to Delta >= standard(A, Rb).
///
/// Both say: Delta is large enough to compensate a certain negative term.
///
/// APPROACH TO PROOF:
/// Perhaps we should prove BOTH Sa >=_LR B AND staircase(A, t*Tb) >= 0
/// simultaneously, as they are both statements about Delta being "large enough."
///
/// Or better: since Sa >=_LR B is equivalent to Delta >= staircase(t*La, B),
/// and staircase(t*La, B) involves La (zone L raw for source A) and B,
/// perhaps we can show |staircase(t*La, B)| is bounded by Delta using
/// properties of La.
///
/// When pa = 1 (leftmost valid position): La = 0 (no q <= pa-2 = -1).
/// Then Sa = total source for A, and Sa >=_LR B is the full claim.
///
/// When pa > 1: La counts source permutations with q <= pa-2.
/// The size of La relative to Sa determines how "tight" the bound is.
///
/// ALTERNATIVE APPROACH (bypassing the two-term decomposition):
/// Go back to the zone-pair decomposition and notice that:
///   The "column M" group (L,M)+(M,M)+(R,M) = A(k1)*Mb(k2-1) - A(k2)*Mb(k1-1)
/// is always >= 0. This is a staircase minor of [A; t*Mb].
///
/// Can we prove A >=_{staircase LR} t*Mb?
///   A(k1)/A(k2) >= Mb(k1-1)/Mb(k2-1) for k1 < k2.
///
/// Mb(s) is the zone-M raw polynomial for source B:
///   Mb(s) = sum_{pa-1 <= q <= pb-2} N_b(s, q)
///
/// This is a PART of the source for B, summing over the "M zone" positions of n-1.
///
/// A is the full output for pa. A includes both "shifted" (La) and "unshifted" (Sa) parts.
/// Mb is a raw source polynomial (no shift). So Mb is "lighter" than A (fewer swaps on average).
///
/// The staircase LR A >=_staircase Mb means the "leftward pull" of A
/// (from the raw source) dominates the "rightward shift" applied to Mb.
///
/// SIMPLEST POSSIBLE PROOF ROUTE:
/// Show that A(k)/A(k+1) >= Mb(k-1)/Mb(k) for all k.
/// This is the "adjacent minors" version of staircase LR.
///
/// If A and Mb are both log-concave (real-rooted), then A(k)/A(k+1) and Mb(k)/Mb(k+1)
/// are both non-increasing. But we need a CROSS comparison.
///
/// Actually, let me check: are A and t*Mb interlacing? Or compatible?
/// If they interlace or are compatible, the staircase LR follows.
///
use combpoly::fixed_descent::{
    permutations_grouped_by_descent_set_bitmask,
    q_refined_modified_source_distribution_from_grouped_source_permutations,
    valid_insertion_positions_for_target_descent_set as valid_positions,
};
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::is_real_rooted;
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

fn poly_sum2(a: &[i64], b: &[i64]) -> Vec<i64> {
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

/// Multiply polynomial by t (shift coefficients right by 1)
fn poly_mul_t(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; p.len() + 1];
    for (i, &c) in p.iter().enumerate() {
        r[i + 1] = c;
    }
    r
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    println!("=== Investigating Sa >=_LR B and compatibility/interlacing ===\n");

    // Check: are A and t*Mb compatible? Interlacing?
    // If f and g are compatible (every alpha*f + beta*g is real-rooted),
    // does that imply staircase LR?
    //
    // Compatibility: f and g have interlacing roots (Hermite-Biehler).
    // Staircase LR for (f, t*g): f(k)/f(k+1) >= g(k-1)/g(k).
    //
    // If f and g are compatible, then f and t*g... well, t*g has an extra root at 0.
    // f and t*g may or may not be compatible.

    let mut a_tmb_compat = 0u64;
    let mut a_tmb_not_compat = 0u64;
    let mut sa_b_compat = 0u64;
    let mut sa_b_not_compat = 0u64;

    // Also check: does Sa interlace B? Or vice versa?
    // Interlacing means one polynomial's roots separate the other's.

    // And check: is A + alpha * t*Mb real-rooted for all alpha >= 0?
    // This is one-sided compatibility.
    let mut a_tmb_one_sided = 0u64;
    let mut a_tmb_one_sided_fail = 0u64;

    // Check: Sa >=_LR B for ALL (pa, pb) pairs, not just consecutive
    let mut sa_lr_b_all_pass = 0u64;
    let mut sa_lr_b_all_fail = 0u64;

    // Additional: check column-M formula for ALL position pairs
    let mut col_m_all_pass = 0u64;
    let mut col_m_all_fail = 0u64;

    for n in 4..=max_n {
        let grouped_source_permutations = permutations_grouped_by_descent_set_bitmask(n - 1);

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
                    q_refined_modified_source_distribution_from_grouped_source_permutations(
                        s_mask,
                        p,
                        n,
                        &grouped_source_permutations,
                    )
                    .counts_by_modified_source_swaps_and_previous_maximum_position,
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

            // Check ALL pairs (pa, pb) with pa < pb (not just consecutive)
            for i in 0..vp.len() {
                for j in (i + 1)..vp.len() {
                    let p_a = vp[i];
                    let p_b = vp[j];

                    let dist_a = &source_dists[&p_a];
                    let dist_b = &source_dists[&p_b];

                    let ma_raw = build_zone_raw(dist_a, p_a, p_b, Zone::M);
                    let ra_raw = build_zone_raw(dist_a, p_a, p_b, Zone::R);
                    let mb_raw = build_zone_raw(dist_b, p_a, p_b, Zone::M);

                    let a_poly_ref = lp_polys.get(&p_a);
                    let b_poly_ref = lp_polys.get(&p_b);
                    if a_poly_ref.is_none() || b_poly_ref.is_none() {
                        continue;
                    }
                    let a_poly = a_poly_ref.unwrap();
                    let b_poly = b_poly_ref.unwrap();

                    let sa_raw = poly_sum2(&ma_raw, &ra_raw);
                    let t_mb = poly_mul_t(&mb_raw);

                    // Check Sa >=_LR B
                    let max_deg = a_poly.len().max(b_poly.len()).max(sa_raw.len());
                    let mut sa_lr_b = true;
                    for k1 in 0..max_deg {
                        for k2 in (k1 + 1)..max_deg {
                            let sa_k1 = coeff(&sa_raw, k1);
                            let sa_k2 = coeff(&sa_raw, k2);
                            let b_k1 = coeff(b_poly, k1);
                            let b_k2 = coeff(b_poly, k2);
                            if sa_k1 * b_k2 < sa_k2 * b_k1 {
                                sa_lr_b = false;
                            }
                        }
                    }
                    if sa_lr_b {
                        sa_lr_b_all_pass += 1;
                    } else {
                        sa_lr_b_all_fail += 1;
                    }

                    // Check column-M: A(k1)*Mb(k2-1) >= A(k2)*Mb(k1-1) for all k1 < k2
                    let mut col_m_ok = true;
                    for k1 in 0..max_deg {
                        for k2 in (k1 + 1)..max_deg {
                            let a_k1 = coeff(a_poly, k1);
                            let a_k2 = coeff(a_poly, k2);
                            let mb_k1m1 = if k1 > 0 { coeff(&mb_raw, k1 - 1) } else { 0 };
                            let mb_k2m1 = if k2 > 0 { coeff(&mb_raw, k2 - 1) } else { 0 };
                            if a_k1 * mb_k2m1 < a_k2 * mb_k1m1 {
                                col_m_ok = false;
                            }
                        }
                    }
                    if col_m_ok {
                        col_m_all_pass += 1;
                    } else {
                        col_m_all_fail += 1;
                    }

                    // Only for consecutive pairs: check compatibility
                    if j == i + 1 {
                        // Check interlacing of A and t*Mb
                        if a_poly.len() >= 2 && t_mb.len() >= 2 {
                            // Compatible = every convex combination is real-rooted
                            // Check by testing A + alpha * t*Mb for several alpha values
                            let mut compat = true;
                            for &alpha in &[1i64, 2, 5, 10, 100] {
                                let max_len = a_poly.len().max(t_mb.len());
                                let mut sum = vec![0i64; max_len];
                                for (i, &c) in a_poly.iter().enumerate() {
                                    sum[i] += c;
                                }
                                for (i, &c) in t_mb.iter().enumerate() {
                                    sum[i] += alpha * c;
                                }
                                while sum.len() > 1 && *sum.last().unwrap() == 0 {
                                    sum.pop();
                                }
                                if sum.len() >= 2 && !is_real_rooted(&sum) {
                                    compat = false;
                                    break;
                                }
                            }
                            if compat {
                                a_tmb_compat += 1;
                            } else {
                                a_tmb_not_compat += 1;
                            }
                        }

                        if sa_raw.len() >= 2 && b_poly.len() >= 2 {
                            let mut compat = true;
                            for &alpha in &[1i64, 2, 5, 10, 100] {
                                let max_len = sa_raw.len().max(b_poly.len());
                                let mut sum = vec![0i64; max_len];
                                for (i, &c) in sa_raw.iter().enumerate() {
                                    sum[i] += c;
                                }
                                for (i, &c) in b_poly.iter().enumerate() {
                                    sum[i] += alpha * c;
                                }
                                while sum.len() > 1 && *sum.last().unwrap() == 0 {
                                    sum.pop();
                                }
                                if sum.len() >= 2 && !is_real_rooted(&sum) {
                                    compat = false;
                                    break;
                                }
                            }
                            if compat {
                                sa_b_compat += 1;
                            } else {
                                sa_b_not_compat += 1;
                            }
                        }

                        // Check one-sided: A + alpha*t*Mb real-rooted for several alpha
                        let mut one_sided_ok = true;
                        for &alpha in &[1i64, 2, 5, 10, 100] {
                            let max_len = a_poly.len().max(t_mb.len());
                            let mut sum = vec![0i64; max_len];
                            for (i, &c) in a_poly.iter().enumerate() {
                                sum[i] += c;
                            }
                            for (i, &c) in t_mb.iter().enumerate() {
                                sum[i] += alpha * c;
                            }
                            while sum.len() > 1 && *sum.last().unwrap() == 0 {
                                sum.pop();
                            }
                            if sum.len() >= 2 && !is_real_rooted(&sum) {
                                one_sided_ok = false;
                            }
                        }
                        if one_sided_ok {
                            a_tmb_one_sided += 1;
                        } else {
                            a_tmb_one_sided_fail += 1;
                        }
                    }
                }
            }
        }
    }

    println!();
    println!("{}", "=".repeat(80));
    println!("RESULTS (n=4..{})", max_n);
    println!("{}", "=".repeat(80));

    let all_pairs = sa_lr_b_all_pass + sa_lr_b_all_fail;
    println!("\nSa >=_LR B (ALL position pairs, not just consecutive):");
    println!(
        "  Pass: {}/{} {}",
        sa_lr_b_all_pass,
        all_pairs,
        if sa_lr_b_all_fail == 0 {
            "ALWAYS"
        } else {
            "FAILS"
        }
    );

    println!("\nColumn-M = A*Mb staircase (ALL position pairs):");
    println!(
        "  Pass: {}/{} {}",
        col_m_all_pass,
        all_pairs,
        if col_m_all_fail == 0 {
            "ALWAYS"
        } else {
            "FAILS"
        }
    );

    let consec = a_tmb_compat + a_tmb_not_compat;
    println!("\nCompatibility checks (consecutive pairs only):");
    println!("  A and t*Mb compatible: {}/{}", a_tmb_compat, consec);
    println!(
        "  Sa and B compatible:   {}/{}",
        sa_b_compat,
        sa_b_compat + sa_b_not_compat
    );

    println!("\nOne-sided compatibility A + alpha*t*Mb:");
    println!(
        "  Pass: {}/{}",
        a_tmb_one_sided,
        a_tmb_one_sided + a_tmb_one_sided_fail
    );

    if sa_lr_b_all_fail == 0 {
        println!("\n*** Sa >=_LR B holds for ALL position pairs, not just consecutive! ***");
        println!("This is a STRONGER statement that might be easier to prove by induction.");
    }

    if col_m_all_fail == 0 {
        println!("\n*** Column-M (staircase A, t*Mb) >= 0 holds for ALL position pairs! ***");
    }
}
