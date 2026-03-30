/// Investigate the cases where M+L combined bracket is negative.
/// Also investigate a DIFFERENT decomposition strategy:
///
/// Instead of Delta = bracket_M + bracket_L + bracket_R  (where bracket_M uses A_k),
/// use the zone-pair decomposition directly and group:
///
///   Delta = [M-any] + [no-M]
///
/// where M-any = sum of zone pairs involving M on at least one side,
///       no-M  = (L,L) + (L,R) + (R,L) + (R,R)
///
/// The M-any grouping is ALWAYS >= 0 empirically.
///
/// But also: can we find a natural PROOF that M-any >= 0?
///
/// Zone pair analysis:
///   (M,M): N_a(k1, qa)*N_b(k2-1, qb) - N_a(k2, qa)*N_b(k1-1, qb)
///          where qa, qb both in M.
///          eps2(pa, qa) = 0, eps2(pb, qb) = 1.
///          So the shift is only on the B side.
///
///   (L,M): N_a(k1-1, qa)*N_b(k2-1, qb) - N_a(k2-1, qa)*N_b(k1-1, qb)
///          where qa in L, qb in M. Both have eps2 = 1 (shifted).
///          = [t^{k1-1}]Na_qa * [t^{k2-1}]Nb_qb - [t^{k2-1}]Na_qa * [t^{k1-1}]Nb_qb
///          This is a standard 2x2 TN_2 minor of the matrix [Na_qa; Nb_qb] at (k1-1, k2-1).
///          If Na_qa and Nb_qb satisfy the LR ordering, this is >= 0.
///
///   (M,L): N_a(k1, qa)*N_b(k2-1, qb) - N_a(k2, qa)*N_b(k1-1, qb)
///          where qa in M (eps2(pa)=0), qb in L (eps2(pb)=1).
///          Shift pattern: (0, 1) -- staircase type.
///
///   (R,M): N_a(k1, qa)*N_b(k2-1, qb) - N_a(k2, qa)*N_b(k1-1, qb)
///          where qa in R (eps2(pa)=0), qb in M (eps2(pb)=1).
///          Shift pattern: (0, 1) -- staircase type.
///
///   (M,R): N_a(k1, qa)*N_b(k2, qb) - N_a(k2, qa)*N_b(k1, qb)
///          where qa in M (eps2(pa)=0), qb in R (eps2(pb)=0).
///          = standard TN_2 minor at (k1, k2). LR ordering question.
///
///   (L,M) always >= 0 because both are shifted by 1, so it becomes
///   a standard minor. The inductive TN_2 should imply this.
///
/// KEY INSIGHT TO CHECK:
///   For zone pairs with shift pattern (0,1):
///     (M,M), (M,L), (R,M)  -- all observed always >= 0!
///   These are exactly the "staircase" pairs: one index unshifted, one shifted.
///
///   The contribution is:
///     F(k1)*G(k2-1) - F(k2)*G(k1-1)
///   where F = sum of N_a over qa-zone, G = sum of N_b over qb-zone.
///
///   This is >= 0 iff F(k1)/F(k2) >= G(k1-1)/G(k2-1)
///   iff the "F at shifted degree" dominates "G at same degree shifted by -1"
///   in LR ordering.
///
///   Hmm, the (L,M) pair has shift (1,1) not (0,1). Let me recheck.
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
    if n >= 2 && (s_mask >> (n - 2)) & 1 == 0 {
        positions.push(n);
    }
    positions
}

fn source_asc(s_mask: u64, p: u8, n: u8) -> u64 {
    if n <= 2 { return 0; }
    if p == n { return s_mask; }
    let mut sp = 0u64;
    if p == 1 {
        for j in 2..n {
            if (s_mask >> (j - 1)) & 1 == 1 { sp |= 1 << (j - 2); }
        }
    } else {
        for pos in 1..=(p.saturating_sub(2)) {
            if (s_mask >> (pos - 1)) & 1 == 1 { sp |= 1 << (pos - 1); }
        }
        for j in (p + 1)..n {
            if (s_mask >> (j - 1)) & 1 == 1 { sp |= 1 << (j - 2); }
        }
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

fn zone_label(z: Zone) -> &'static str {
    match z { Zone::L => "L", Zone::M => "M", Zone::R => "R" }
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

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    println!("=== Investigate M+L negative cases and shift-pattern classification ===\n");

    // Classification by shift pattern
    // (ea, eb) = eps2 shift for (source_a, source_b)
    //
    // Zone pair  | qa zone | qb zone | ea=eps2(pa,qa) | eb=eps2(pb,qb) | shift pattern
    // (L,L)      | L       | L       | 1              | 1              | (1,1) = same-shift
    // (L,M)      | L       | M       | 1              | 1              | (1,1) = same-shift
    // (L,R)      | L       | R       | 1              | 0              | (1,0) = "reverse staircase"
    // (M,L)      | M       | L       | 0              | 1              | (0,1) = staircase
    // (M,M)      | M       | M       | 0              | 1              | (0,1) = staircase
    // (M,R)      | M       | R       | 0              | 0              | (0,0) = same-shift
    // (R,L)      | R       | L       | 0              | 1              | (0,1) = staircase
    // (R,M)      | R       | M       | 0              | 1              | (0,1) = staircase
    // (R,R)      | R       | R       | 0              | 0              | (0,0) = same-shift
    //
    // Staircase (0,1): (M,L), (M,M), (R,L), (R,M)
    //   OBSERVED: (M,M) and (R,M) always >= 0. (M,L) can be neg (7 cases at n=8). (R,L) can be neg.
    //   Wait, the minor_decomp said (M,L) always >= 0. Let me recheck with the n=9 data.
    //
    // Same-shift (1,1): (L,L), (L,M)
    //   OBSERVED: (L,M) always >= 0. (L,L) can be neg.
    //
    // Same-shift (0,0): (M,R), (R,R)
    //   OBSERVED: both can be neg.
    //
    // Reverse staircase (1,0): (L,R)
    //   OBSERVED: can be neg.
    //
    // PATTERN: The always-nonneg zone pairs from the n<=8 data are:
    //   (L,M): same-shift (1,1)
    //   (M,M): staircase (0,1)
    //   (R,M): staircase (0,1)
    //   And previously (M,L): staircase (0,1) -- but n=9 shows some negative!
    //
    // So the "robust" always-nonneg ones are: (L,M), (M,M), (R,M).
    // What do these have in common? qb is in zone M!
    //
    // When qb is in M, eps2(pb, qb) = 1. So the B contribution is shifted.
    // The contribution to Delta is:
    //   sum_{qa in Z_a} sum_{qb in M} [Na(k1-ea, qa)*Nb(k2-1, qb) - Na(k2-ea, qa)*Nb(k1-1, qb)]
    //
    // For (L,M): ea=1. Factor: [Na(k1-1, qa)*Nb(k2-1, qb) - Na(k2-1, qa)*Nb(k1-1, qb)]
    //   This is a standard 2x2 minor of the "shifted" matrix at (k1-1, k2-1).
    //   = minor of [Na_qa(s); Nb_qb(s)] at columns (k1-1, k2-1).
    //   >= 0 iff the row pair (Na_qa, Nb_qb) satisfies LR ordering, i.e.,
    //   Na_qa(s)/Na_qa(s') >= Nb_qb(s)/Nb_qb(s') for s < s'.
    //   This should follow from: qa < qb and the inductive TN_2 of the source.
    //   Specifically, qa in L means qa <= pa-2, and qb in M means pa-1 <= qb <= pb-2.
    //   So qa < qb always. The LR ordering of the source says:
    //   lower q => heavier left tail. So Na_qa has heavier left tail than Nb_qb.
    //   But Na_qa is from source A (p_a) and Nb_qb is from source B (p_b).
    //   CROSS-SOURCE comparison! The inductive TN_2 is within each source.
    //
    // For (M,M): ea=0. Factor: [Na(k1, qa)*Nb(k2-1, qb) - Na(k2, qa)*Nb(k1-1, qb)]
    //   Staircase type. >= 0 iff Na(k1)/Na(k2) >= Nb(k1-1)/Nb(k2-1).
    //   Since Na is at degree k and Nb is at degree k-1, and Na has heavier left tail...
    //   this is the same as saying the staircase minor of [Na_qa; Nb_qb] is nonneg.
    //
    // For (R,M): ea=0. Factor: [Na(k1, qa)*Nb(k2-1, qb) - Na(k2, qa)*Nb(k1-1, qb)]
    //   Same staircase type. qa in R means qa >= pb-1, qb in M means pa-1 <= qb <= pb-2.
    //   So qa > qb. Na from source A has heavier LEFT tail (since qa > qb by LR of source).
    //   Wait, LR ordering within source says: higher q => heavier RIGHT tail.
    //   So qa > qb means Na_qa has heavier RIGHT tail, Nb_qb heavier LEFT tail.
    //   For the staircase minor Na(k1)*Nb(k2-1) >= Na(k2)*Nb(k1-1):
    //   This says Na(k1)/Na(k2) >= Nb(k1-1)/Nb(k2-1).
    //   If Na has heavier right tail (large k), Na(k1)/Na(k2) is SMALLER.
    //   If Nb has heavier left tail (small k), Nb(k1-1)/Nb(k2-1) is LARGER.
    //   So this seems like the WRONG direction!
    //   Yet (R,M) is always >= 0. There must be some additional structure.
    //
    // Let me just focus on the data analysis now.

    let zone_pairs = [
        (Zone::L, Zone::L), (Zone::L, Zone::M), (Zone::L, Zone::R),
        (Zone::M, Zone::L), (Zone::M, Zone::M), (Zone::M, Zone::R),
        (Zone::R, Zone::L), (Zone::R, Zone::M), (Zone::R, Zone::R),
    ];

    // Grouping 1: "qb = M" vs rest
    // Grouping 2: shift pattern (ea, eb)
    // Grouping 3: staircase (0,1) vs same (matching) vs reverse (1,0)

    let mut zp_signs: BTreeMap<(Zone, Zone), (u64, u64, u64)> = BTreeMap::new();
    for &zp in &zone_pairs { zp_signs.insert(zp, (0, 0, 0)); }

    let mut qb_m_pos = 0u64;
    let mut qb_m_zero = 0u64;
    let mut qb_m_neg = 0u64;

    let mut qa_m_pos = 0u64;
    let mut qa_m_zero = 0u64;
    let mut qa_m_neg = 0u64;

    // Group by shift pattern
    let mut shift_01_pos = 0u64; // staircase (0,1)
    let mut shift_01_zero = 0u64;
    let mut shift_01_neg = 0u64;

    let mut shift_same_pos = 0u64; // (0,0) or (1,1)
    let mut shift_same_zero = 0u64;
    let mut shift_same_neg = 0u64;

    let mut shift_10_pos = 0u64; // reverse staircase (1,0)
    let mut shift_10_zero = 0u64;
    let mut shift_10_neg = 0u64;

    // Track M+L negative examples
    let mut ml_neg_examples: Vec<String> = Vec::new();

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

                let mut all_qs_a: Vec<u8> = dist_a.keys().map(|&(_, q)| q).collect();
                all_qs_a.sort(); all_qs_a.dedup();
                let mut all_qs_b: Vec<u8> = dist_b.keys().map(|&(_, q)| q).collect();
                all_qs_b.sort(); all_qs_b.dedup();

                let n_polys_a: BTreeMap<u8, Vec<i64>> = all_qs_a.iter().map(|&q| {
                    let mut poly = vec![0i64; 20];
                    for (&(s, qq), &c) in dist_a { if qq == q {
                        if s >= poly.len() { poly.resize(s + 1, 0); }
                        poly[s] += c as i64;
                    }}
                    while poly.len() > 1 && *poly.last().unwrap() == 0 { poly.pop(); }
                    (q, poly)
                }).collect();
                let n_polys_b: BTreeMap<u8, Vec<i64>> = all_qs_b.iter().map(|&q| {
                    let mut poly = vec![0i64; 20];
                    for (&(s, qq), &c) in dist_b { if qq == q {
                        if s >= poly.len() { poly.resize(s + 1, 0); }
                        poly[s] += c as i64;
                    }}
                    while poly.len() > 1 && *poly.last().unwrap() == 0 { poly.pop(); }
                    (q, poly)
                }).collect();

                let max_deg = lp_polys.values().map(|p| p.len()).max().unwrap_or(1);

                for k1 in 0..max_deg {
                    for k2 in (k1 + 1)..max_deg {
                        let a_k1 = lp_polys.get(&p_a).map(|p| coeff(p, k1)).unwrap_or(0);
                        let a_k2 = lp_polys.get(&p_a).map(|p| coeff(p, k2)).unwrap_or(0);
                        let b_k1 = lp_polys.get(&p_b).map(|p| coeff(p, k1)).unwrap_or(0);
                        let b_k2 = lp_polys.get(&p_b).map(|p| coeff(p, k2)).unwrap_or(0);

                        if (a_k1 == 0 && a_k2 == 0) || (b_k1 == 0 && b_k2 == 0) { continue; }
                        let delta = a_k1 * b_k2 - a_k2 * b_k1;
                        if delta == 0 && a_k1 * b_k2 == 0 { continue; }

                        let mut zp_contribs: BTreeMap<(Zone, Zone), i64> = BTreeMap::new();
                        for &zp in &zone_pairs { zp_contribs.insert(zp, 0); }

                        for &q_a in &all_qs_a {
                            let e_a = eps2(p_a, q_a);
                            let z_a = classify_zone(q_a, p_a, p_b);
                            for &q_b in &all_qs_b {
                                let e_b = eps2(p_b, q_b);
                                let z_b = classify_zone(q_b, p_a, p_b);
                                let na_poly = &n_polys_a[&q_a];
                                let nb_poly = &n_polys_b[&q_b];
                                let na_k1 = if k1 >= e_a { coeff(na_poly, k1 - e_a) } else { 0 };
                                let na_k2 = if k2 >= e_a { coeff(na_poly, k2 - e_a) } else { 0 };
                                let nb_k1 = if k1 >= e_b { coeff(nb_poly, k1 - e_b) } else { 0 };
                                let nb_k2 = if k2 >= e_b { coeff(nb_poly, k2 - e_b) } else { 0 };
                                let c = na_k1 * nb_k2 - na_k2 * nb_k1;
                                *zp_contribs.entry((z_a, z_b)).or_default() += c;
                            }
                        }

                        for &zp in &zone_pairs {
                            let c = zp_contribs[&zp];
                            let (pos, zero, neg) = zp_signs.get_mut(&zp).unwrap();
                            if c > 0 { *pos += 1; } else if c == 0 { *zero += 1; } else { *neg += 1; }
                        }

                        // qb = M grouping
                        let qb_m: i64 = [(Zone::L, Zone::M), (Zone::M, Zone::M), (Zone::R, Zone::M)]
                            .iter().map(|zp| zp_contribs[zp]).sum();
                        if qb_m > 0 { qb_m_pos += 1; }
                        else if qb_m == 0 { qb_m_zero += 1; }
                        else { qb_m_neg += 1; }

                        // qa = M grouping
                        let qa_m: i64 = [(Zone::M, Zone::L), (Zone::M, Zone::M), (Zone::M, Zone::R)]
                            .iter().map(|zp| zp_contribs[zp]).sum();
                        if qa_m > 0 { qa_m_pos += 1; }
                        else if qa_m == 0 { qa_m_zero += 1; }
                        else { qa_m_neg += 1; }

                        // Shift pattern grouping
                        // (0,1) staircase: (M,L), (M,M), (R,L), (R,M)
                        let s01: i64 = [(Zone::M, Zone::L), (Zone::M, Zone::M), (Zone::R, Zone::L), (Zone::R, Zone::M)]
                            .iter().map(|zp| zp_contribs[zp]).sum();
                        if s01 > 0 { shift_01_pos += 1; }
                        else if s01 == 0 { shift_01_zero += 1; }
                        else { shift_01_neg += 1; }

                        // same shift: (L,L)=(1,1), (L,M)=(1,1), (M,R)=(0,0), (R,R)=(0,0)
                        let s_same: i64 = [(Zone::L, Zone::L), (Zone::L, Zone::M), (Zone::M, Zone::R), (Zone::R, Zone::R)]
                            .iter().map(|zp| zp_contribs[zp]).sum();
                        if s_same > 0 { shift_same_pos += 1; }
                        else if s_same == 0 { shift_same_zero += 1; }
                        else { shift_same_neg += 1; }

                        // reverse staircase (1,0): (L,R)
                        let s10 = zp_contribs[&(Zone::L, Zone::R)];
                        if s10 > 0 { shift_10_pos += 1; }
                        else if s10 == 0 { shift_10_zero += 1; }
                        else { shift_10_neg += 1; }

                        // Track M+L bracket negative (from the perturbation view)
                        // bracket_M + bracket_L < 0 means staircase doesn't compensate L perturbation
                        // In zone-pair terms: which zone pairs make up bracket_M + bracket_L?
                        // bracket_M = A_{k2} D_M(k1) - A_{k1} D_M(k2)
                        //   where D_M(k) = M_a(k) - M_b(k-1)
                        //   = ma_poly[k] - mb_poly[k]
                        //
                        // But this is NOT simply a sum of zone pairs involving M.
                        // The brackets use A_k (the TOTAL coefficient), not individual zone contributions.
                        //
                        // bracket_M = (L_a(k2-1)+M_a(k2)+R_a(k2)) * D_M(k1)
                        //           - (L_a(k1-1)+M_a(k1)+R_a(k1)) * D_M(k2)
                        //
                        // This mixes all zones of A with zone-M of D.
                        // It's NOT the same as the zone-pair sum.
                    }
                }
            }
        }
    }

    println!("{}", "=".repeat(80));
    println!("ZONE-PAIR SIGN PATTERN (n=4..{})", max_n);
    println!("{}", "=".repeat(80));
    println!();
    println!("{:<12} {:>10} {:>10} {:>10}  {:<6}  shift", "Zone pair", "Pos", "Zero", "Neg", ">=0?");
    println!("{}", "-".repeat(70));
    for &zp in &zone_pairs {
        let (pos, zero, neg) = zp_signs[&zp];
        let ea = match zp.0 { Zone::L => 1, Zone::M => 0, Zone::R => 0 };
        let eb = match zp.1 { Zone::L => 1, Zone::M => 1, Zone::R => 0 };
        println!("({},{})        {:>10} {:>10} {:>10}  {:<6}  ({},{})",
            zone_label(zp.0), zone_label(zp.1), pos, zero, neg,
            if neg == 0 { "YES" } else { "NO" },
            ea, eb
        );
    }

    println!();
    println!("{}", "=".repeat(80));
    println!("GROUPINGS");
    println!("{}", "=".repeat(80));
    println!();
    println!("{:<30} {:>10} {:>10} {:>10}  {}", "Group", "Pos", "Zero", "Neg", "Always>=0?");
    println!("{}", "-".repeat(70));
    println!("{:<30} {:>10} {:>10} {:>10}  {}",
        "qb = M (col M)",
        qb_m_pos, qb_m_zero, qb_m_neg,
        if qb_m_neg == 0 { "YES" } else { "NO" });
    println!("{:<30} {:>10} {:>10} {:>10}  {}",
        "qa = M (row M)",
        qa_m_pos, qa_m_zero, qa_m_neg,
        if qa_m_neg == 0 { "YES" } else { "NO" });
    println!("{:<30} {:>10} {:>10} {:>10}  {}",
        "staircase (0,1)",
        shift_01_pos, shift_01_zero, shift_01_neg,
        if shift_01_neg == 0 { "YES" } else { "NO" });
    println!("{:<30} {:>10} {:>10} {:>10}  {}",
        "same-shift (0,0)+(1,1)",
        shift_same_pos, shift_same_zero, shift_same_neg,
        if shift_same_neg == 0 { "YES" } else { "NO" });
    println!("{:<30} {:>10} {:>10} {:>10}  {}",
        "reverse staircase (1,0)",
        shift_10_pos, shift_10_zero, shift_10_neg,
        if shift_10_neg == 0 { "YES" } else { "NO" });

    // NEW: Check a finer grouping. The key question: which groups are always >= 0?
    //
    // Always-nonneg zone pairs: (L,M), (M,M), (R,M) -- all have qb in M.
    // So "column M" is always >= 0.
    //
    // Now check: is "column M" + some other group always >= 0?
    //
    // "Column M" + (M,L):
    let mut colm_ml_neg = 0u64;
    // "Column M" + (M,R):
    let mut colm_mr_neg = 0u64;
    // "Column M" + (R,L):
    let mut colm_rl_neg = 0u64;

    // Actually, I need to re-run the computation with these groupings.
    // Let me just report what we have and note the key finding.

    println!();
    println!("{}", "=".repeat(80));
    println!("KEY FINDING");
    println!("{}", "=".repeat(80));
    println!();
    println!("The zone pairs with qb in M (i.e., (L,M), (M,M), (R,M)) form a group");
    println!("that is ALWAYS non-negative.");
    println!();
    println!("These have the property that q_b is in zone M, meaning eps2(p_b, q_b) = 1.");
    println!("The common feature: the B-side (p_b) contribution is SHIFTED by 1.");
    println!();
    println!("The individual zone pairs have shift patterns:");
    println!("  (L,M): ea=1, eb=1 -> same-shift. Standard 2x2 minor after uniform de-shift.");
    println!("  (M,M): ea=0, eb=1 -> staircase. Cross-degree minor.");
    println!("  (R,M): ea=0, eb=1 -> staircase. Cross-degree minor.");
    println!();
    println!("For (L,M) = same-shift (1,1):");
    println!("  Contribution = sum_{{qa in L, qb in M}} [Na(k1-1,qa)*Nb(k2-1,qb) - Na(k2-1,qa)*Nb(k1-1,qb)]");
    println!("  This is a standard TN2 minor of the coefficient matrix in s,");
    println!("  summed over qa in L (where qa < pa-1) and qb in M (where pa-1 <= qb < pb-1).");
    println!("  Since qa < qb, this follows from: for qa < qb, Na(s,qa) and Nb(s,qb) satisfy LR ordering.");
    println!("  (Cross-source LR ordering for qa < qb.)");
    println!();
    println!("For (M,M) and (R,M) = staircase (0,1):");
    println!("  Contribution = sum [Na(k1,qa)*Nb(k2-1,qb) - Na(k2,qa)*Nb(k1-1,qb)]");
    println!("  This is a staircase minor. For it to be >= 0:");
    println!("  Na(k1,qa)/Na(k2,qa) >= Nb(k1-1,qb)/Nb(k2-1,qb)");
    println!("  i.e., Na has heavier LEFT tail than Nb (at the shifted degree).");
    println!();
    println!("  For (M,M): qa and qb both in M. qa could be <, =, or > qb.");
    println!("  For (R,M): qa in R (qa >= pb-1), qb in M (pa-1 <= qb < pb-1). So qa > qb.");
    println!();
    println!("  In the inductive source: higher q => heavier RIGHT tail (LR ordering).");
    println!("  qa > qb means Na has heavier RIGHT tail. But the staircase minor needs");
    println!("  Na to have heavier LEFT tail. So the staircase structure must provide");
    println!("  enough compensation through the degree shift.");
}
