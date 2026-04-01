/// Investigate the LR (likelihood ratio) ordering for TN_2 proof of the
/// position-refined swaps coefficient matrix.
///
/// For S subset [n-1] with 1 not in S, order valid positions p_1 < p_2 < ... < p_m.
/// The coefficient matrix M_{i,k} = [t^k] L_{n,S}^{(p_i)}(t) is conjectured TN_2.
///
/// Equivalent: for p_a < p_b and all k:
///   A_k * B_{k+1} >= A_{k+1} * B_k
///
/// Via insertion: sigma with sigma^{-1}(n) = p arises from inserting n at position p
/// in some pi in S_{n-1}. swaps(sigma) = swaps(pi) + eps1(pi,p) + eps2(p,q)
/// where q = pi^{-1}(n-1).
///
/// For consecutive valid positions (p_a, p_b), we analyze:
/// 1. Source descent set comparison (S'_{p_a} vs S'_{p_b}, S''_{p_a} vs S''_{p_b})
/// 2. Shared vs exclusive source permutations
/// 3. How eps2 thresholds differ
/// 4. Whether a natural injection exists for the LR inequality
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use std::collections::{BTreeMap, BTreeSet};

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

fn mask_to_string(mask: u64, n: u8) -> String {
    descent_set_to_string(mask, n)
}

fn coeff(poly: &[i64], k: usize) -> i64 {
    if k < poly.len() {
        poly[k]
    } else {
        0
    }
}

/// Check TN_2 (all 2x2 minors >= 0) for a matrix given as rows of polynomials.
fn check_tn2(rows: &[Vec<i64>]) -> bool {
    let max_deg = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    for i in 0..rows.len() {
        for j in (i + 1)..rows.len() {
            for k1 in 0..max_deg {
                for k2 in (k1 + 1)..max_deg {
                    let minor = coeff(&rows[i], k1) * coeff(&rows[j], k2)
                        - coeff(&rows[i], k2) * coeff(&rows[j], k1);
                    if minor < 0 {
                        return false;
                    }
                }
            }
        }
    }
    true
}

/// Represent a source permutation with all relevant data.
#[derive(Clone, Debug)]
struct SourceInfo {
    /// The (n-1)-permutation
    perm: Vec<u8>,
    /// swaps(pi) before modification
    base_swaps: usize,
    /// eps1 contribution at the given insertion position p
    eps1: bool,
    /// Position of n-1 in pi (1-indexed)
    q: u8,
    /// Which source descent set this came from ('a' for ascent, 'd' for descent)
    source_type: char,
    /// The source descent set mask
    source_mask: u64,
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    println!("=== LR Ordering Analysis for TN_2 Proof ===");
    println!("Checking n = 4 to {}\n", max_n);

    let mut total_pairs = 0u64;
    let mut same_source_sets = 0u64;
    let mut lr_holds = 0u64;

    for n in 4..=max_n {
        let perms_prev = all_permutations(n - 1);
        let perms = all_permutations(n);

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

        let mut n_pairs = 0u32;
        let mut n_same_src = 0u32;
        let mut n_lr_holds = 0u32;

        for (&mask, _class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            // Gather full source info for each valid position
            let mut sources_by_pos: BTreeMap<u8, Vec<SourceInfo>> = BTreeMap::new();

            for &p in &vp {
                let sp_a = source_asc(mask, p, n);
                let sp_d = source_desc(mask, p, n);

                let mut sources = Vec::new();

                // Ascent sources
                if let Some(cls) = prev_by_des.get(&sp_a) {
                    for pi in cls {
                        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                        let e1 = epsilon1(pi, p);
                        sources.push(SourceInfo {
                            perm: (*pi).clone(),
                            base_swaps: compute(pi, Stat::Swaps),
                            eps1: e1,
                            q,
                            source_type: 'a',
                            source_mask: sp_a,
                        });
                    }
                }

                // Descent sources (only for 1 < p < n)
                if let Some(sp_d_val) = sp_d {
                    if sp_d_val != sp_a {
                        if let Some(cls) = prev_by_des.get(&sp_d_val) {
                            for pi in cls {
                                let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                                let e1 = epsilon1(pi, p);
                                sources.push(SourceInfo {
                                    perm: (*pi).clone(),
                                    base_swaps: compute(pi, Stat::Swaps),
                                    eps1: e1,
                                    q,
                                    source_type: 'd',
                                    source_mask: sp_d_val,
                                });
                            }
                        }
                    }
                }

                sources_by_pos.insert(p, sources);
            }

            // Analyze consecutive pairs
            for w in vp.windows(2) {
                let p_a = w[0];
                let p_b = w[1];
                n_pairs += 1;

                let src_a = sources_by_pos.get(&p_a).unwrap();
                let src_b = sources_by_pos.get(&p_b).unwrap();

                // --- Analysis 1: Compare source descent sets ---
                let sp_a_asc = source_asc(mask, p_a, n);
                let sp_a_desc = source_desc(mask, p_a, n);
                let sp_b_asc = source_asc(mask, p_b, n);
                let sp_b_desc = source_desc(mask, p_b, n);

                let same_asc = sp_a_asc == sp_b_asc;
                let same_desc = sp_a_desc == sp_b_desc;
                let same_both = same_asc && same_desc;

                if same_both {
                    n_same_src += 1;
                }

                // --- Analysis 2: Shared source permutations ---
                // A permutation pi is a source for both p_a and p_b if it appears
                // in the source set for both.
                let perm_set_a: BTreeSet<Vec<u8>> = src_a.iter().map(|s| s.perm.clone()).collect();
                let perm_set_b: BTreeSet<Vec<u8>> = src_b.iter().map(|s| s.perm.clone()).collect();

                let shared: BTreeSet<Vec<u8>> =
                    perm_set_a.intersection(&perm_set_b).cloned().collect();
                let excl_a: BTreeSet<Vec<u8>> =
                    perm_set_a.difference(&perm_set_b).cloned().collect();
                let excl_b: BTreeSet<Vec<u8>> =
                    perm_set_b.difference(&perm_set_a).cloned().collect();

                // --- Analysis 3: eps2 threshold difference ---
                // eps2(p,q) = 1[q <= p-2]
                // For p_a: eps2=1 when q <= p_a-2
                // For p_b: eps2=1 when q <= p_b-2
                // Transition zone: q in [p_a-1, p_b-2] gets eps2=0 for p_a, eps2=1 for p_b

                // --- Analysis 4: LR inequality check ---
                // For the SHARED permutations, check whether the eps2 shift
                // alone implies LR ordering.
                // For shared pi, swaps_a(pi) = swaps(pi) + eps1(pi,p_a) + eps2(p_a,q)
                //                swaps_b(pi) = swaps(pi) + eps1(pi,p_b) + eps2(p_b,q)
                // Note: eps1 can differ between p_a and p_b for the same pi!

                // Compute A_k and B_k from full sigma-level counts
                let mut a_vals: Vec<usize> = Vec::new();
                for si in src_a {
                    let eps2 = if si.q <= p_a.saturating_sub(2) { 1 } else { 0 };
                    let total = si.base_swaps + if si.eps1 { 1 } else { 0 } + eps2;
                    a_vals.push(total);
                }
                let mut b_vals: Vec<usize> = Vec::new();
                for si in src_b {
                    let eps2 = if si.q <= p_b.saturating_sub(2) { 1 } else { 0 };
                    let total = si.base_swaps + if si.eps1 { 1 } else { 0 } + eps2;
                    b_vals.push(total);
                }

                let a_poly = build_poly(&a_vals);
                let b_poly = build_poly(&b_vals);

                // Check LR: A_k * B_{k+1} >= A_{k+1} * B_k for all k
                let max_k = a_poly.len().max(b_poly.len());
                let mut lr_ok = true;
                for k in 0..max_k {
                    let ak = coeff(&a_poly, k);
                    let bk1 = coeff(&b_poly, k + 1);
                    let ak1 = coeff(&a_poly, k + 1);
                    let bk = coeff(&b_poly, k);
                    if ak * bk1 < ak1 * bk {
                        lr_ok = false;
                        break;
                    }
                }
                if lr_ok {
                    n_lr_holds += 1;
                }

                // Detailed output for small n
                if n <= 6 {
                    let s_str = descent_set_to_string(mask, n);
                    println!(
                        "n={} S={} pair ({},{}): same_src_sets={} shared={} excl_a={} excl_b={} LR={}",
                        n, s_str, p_a, p_b, same_both, shared.len(), excl_a.len(), excl_b.len(),
                        if lr_ok { "OK" } else { "FAIL" }
                    );

                    // Show source descent sets
                    println!(
                        "  S'_asc({})={} S'_asc({})={}  same={}",
                        p_a,
                        mask_to_string(sp_a_asc, n - 1),
                        p_b,
                        mask_to_string(sp_b_asc, n - 1),
                        same_asc
                    );
                    if let (Some(da), Some(db)) = (sp_a_desc, sp_b_desc) {
                        println!(
                            "  S'_desc({})={} S'_desc({})={}  same={}",
                            p_a,
                            mask_to_string(da, n - 1),
                            p_b,
                            mask_to_string(db, n - 1),
                            same_desc
                        );
                    } else {
                        println!(
                            "  S'_desc({})={:?} S'_desc({})={:?}",
                            p_a, sp_a_desc, p_b, sp_b_desc
                        );
                    }

                    // eps2 transition zone
                    let tz_lo = p_a.saturating_sub(1);
                    let tz_hi = p_b.saturating_sub(2);
                    println!(
                        "  eps2 transition zone: q in [{}..{}] (size {})",
                        tz_lo,
                        tz_hi,
                        if tz_hi >= tz_lo { tz_hi - tz_lo + 1 } else { 0 }
                    );

                    // For shared permutations, show swaps difference
                    if !shared.is_empty() && shared.len() <= 20 {
                        println!("  Shared permutations (swaps_a, swaps_b, diff):");
                        for pi in &shared {
                            let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            let base_sw = compute(pi, Stat::Swaps);

                            let e1_a = epsilon1(pi, p_a);
                            let eps2_a = if q <= p_a.saturating_sub(2) { 1 } else { 0 };
                            let sw_a = base_sw + if e1_a { 1 } else { 0 } + eps2_a;

                            let e1_b = epsilon1(pi, p_b);
                            let eps2_b = if q <= p_b.saturating_sub(2) { 1 } else { 0 };
                            let sw_b = base_sw + if e1_b { 1 } else { 0 } + eps2_b;

                            let in_tz = q >= tz_lo && q <= tz_hi;
                            println!(
                                "    pi={:?} q={} base_sw={} e1_a={} e1_b={} eps2_a={} eps2_b={} => sw_a={} sw_b={} diff={} in_tz={}",
                                pi, q, base_sw, e1_a as u8, e1_b as u8, eps2_a, eps2_b,
                                sw_a, sw_b, sw_b as i32 - sw_a as i32, in_tz
                            );
                        }
                    }

                    // Show the A and B polynomials
                    println!("  A(t) = {:?}  B(t) = {:?}", a_poly, b_poly);

                    // Decompose: shared vs exclusive contributions
                    let mut shared_a_vals = Vec::new();
                    let mut shared_b_vals = Vec::new();
                    let mut excl_a_vals = Vec::new();
                    let mut excl_b_vals = Vec::new();

                    for si in src_a {
                        let eps2 = if si.q <= p_a.saturating_sub(2) { 1 } else { 0 };
                        let total = si.base_swaps + if si.eps1 { 1 } else { 0 } + eps2;
                        if shared.contains(&si.perm) {
                            shared_a_vals.push(total);
                        } else {
                            excl_a_vals.push(total);
                        }
                    }
                    for si in src_b {
                        let eps2 = if si.q <= p_b.saturating_sub(2) { 1 } else { 0 };
                        let total = si.base_swaps + if si.eps1 { 1 } else { 0 } + eps2;
                        if shared.contains(&si.perm) {
                            shared_b_vals.push(total);
                        } else {
                            excl_b_vals.push(total);
                        }
                    }

                    if !shared_a_vals.is_empty() {
                        let sa_poly = build_poly(&shared_a_vals);
                        let sb_poly = build_poly(&shared_b_vals);
                        println!("  Shared: A_sh(t) = {:?}  B_sh(t) = {:?}", sa_poly, sb_poly);

                        // Check LR just on shared part
                        let max_k2 = sa_poly.len().max(sb_poly.len());
                        let mut shared_lr_ok = true;
                        for k in 0..max_k2 {
                            let ak = coeff(&sa_poly, k);
                            let bk1 = coeff(&sb_poly, k + 1);
                            let ak1 = coeff(&sa_poly, k + 1);
                            let bk = coeff(&sb_poly, k);
                            if ak * bk1 < ak1 * bk {
                                shared_lr_ok = false;
                                break;
                            }
                        }
                        println!("  Shared LR: {}", if shared_lr_ok { "OK" } else { "FAIL" });
                    }
                    if !excl_a_vals.is_empty() || !excl_b_vals.is_empty() {
                        let ea_poly = if excl_a_vals.is_empty() {
                            vec![0]
                        } else {
                            build_poly(&excl_a_vals)
                        };
                        let eb_poly = if excl_b_vals.is_empty() {
                            vec![0]
                        } else {
                            build_poly(&excl_b_vals)
                        };
                        println!("  Excl:   A_ex(t) = {:?}  B_ex(t) = {:?}", ea_poly, eb_poly);
                    }

                    println!();
                }
            }
        }

        total_pairs += n_pairs as u64;
        same_source_sets += n_same_src as u64;
        lr_holds += n_lr_holds as u64;

        println!(
            "--- n={}: {} pairs, {} same source sets, {} LR holds ---\n",
            n, n_pairs, n_same_src, n_lr_holds
        );
    }

    println!("========================================");
    println!("TOTALS (n=4..{}):", max_n);
    println!("  Total consecutive pairs:   {}", total_pairs);
    println!("  Same source descent sets:  {}", same_source_sets);
    println!("  LR ordering holds:         {}", lr_holds);
    println!();

    // =====================================================================
    // PART 2: Deeper analysis of the eps2 shift mechanism
    // =====================================================================
    println!("=== PART 2: Decomposition by q-position ===\n");

    for n in 4..=max_n {
        let perms_prev = all_permutations(n - 1);

        let mut prev_by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        let perms = all_permutations(n);
        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        for (&mask, _class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            for w in vp.windows(2) {
                let p_a = w[0];
                let p_b = w[1];

                // Get source info for both positions
                let mut gather_sources = |p: u8| -> Vec<SourceInfo> {
                    let sp_a_val = source_asc(mask, p, n);
                    let sp_d_val = source_desc(mask, p, n);
                    let mut sources = Vec::new();

                    if let Some(cls) = prev_by_des.get(&sp_a_val) {
                        for pi in cls {
                            let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            sources.push(SourceInfo {
                                perm: (*pi).clone(),
                                base_swaps: compute(pi, Stat::Swaps),
                                eps1: epsilon1(pi, p),
                                q,
                                source_type: 'a',
                                source_mask: sp_a_val,
                            });
                        }
                    }
                    if let Some(sp_d) = sp_d_val {
                        if sp_d != sp_a_val {
                            if let Some(cls) = prev_by_des.get(&sp_d) {
                                for pi in cls {
                                    let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                                    sources.push(SourceInfo {
                                        perm: (*pi).clone(),
                                        base_swaps: compute(pi, Stat::Swaps),
                                        eps1: epsilon1(pi, p),
                                        q,
                                        source_type: 'd',
                                        source_mask: sp_d,
                                    });
                                }
                            }
                        }
                    }
                    sources
                };

                let src_a = gather_sources(p_a);
                let src_b = gather_sources(p_b);

                // Decompose by q-value
                // For p_a: swaps = base + eps1(pi,p_a) + eps2(p_a,q) = base + eps1 + [q<=p_a-2]
                // For p_b: swaps = base + eps1(pi,p_b) + eps2(p_b,q) = base + eps1 + [q<=p_b-2]
                //
                // Three q-ranges:
                //   Zone L: q <= p_a - 2  => eps2=1 for both
                //   Zone M: p_a - 1 <= q <= p_b - 2  => eps2=0 for p_a, eps2=1 for p_b
                //   Zone R: q >= p_b - 1  => eps2=0 for both

                let mut a_by_zone: [Vec<usize>; 3] = [vec![], vec![], vec![]];
                let mut b_by_zone: [Vec<usize>; 3] = [vec![], vec![], vec![]];

                let zone_of = |q: u8, p: u8| -> usize {
                    if q <= p_a.saturating_sub(2) {
                        0 // Zone L
                    } else if q >= p_b.saturating_sub(1) {
                        2 // Zone R
                    } else {
                        1 // Zone M
                    }
                };

                for si in &src_a {
                    let eps2 = if si.q <= p_a.saturating_sub(2) { 1 } else { 0 };
                    let total = si.base_swaps + if si.eps1 { 1 } else { 0 } + eps2;
                    let z = zone_of(si.q, p_a);
                    a_by_zone[z].push(total);
                }
                for si in &src_b {
                    let eps2 = if si.q <= p_b.saturating_sub(2) { 1 } else { 0 };
                    let total = si.base_swaps + if si.eps1 { 1 } else { 0 } + eps2;
                    let z = zone_of(si.q, p_b);
                    b_by_zone[z].push(total);
                }

                let s_str = descent_set_to_string(mask, n);
                println!("n={} S={} pair ({},{}): zones L/M/R", n, s_str, p_a, p_b);
                let zone_names = ["L (q<=pa-2)", "M (pa-1<=q<=pb-2)", "R (q>=pb-1)"];
                for z in 0..3 {
                    let ap = if a_by_zone[z].is_empty() {
                        vec![0]
                    } else {
                        build_poly(&a_by_zone[z])
                    };
                    let bp = if b_by_zone[z].is_empty() {
                        vec![0]
                    } else {
                        build_poly(&b_by_zone[z])
                    };
                    println!(
                        "  {}: #a={} #b={} A={:?} B={:?}",
                        zone_names[z],
                        a_by_zone[z].len(),
                        b_by_zone[z].len(),
                        ap,
                        bp,
                    );
                }

                // KEY OBSERVATION: In zone M, p_b gets an extra +1 from eps2.
                // So if pi is shared and in zone M:
                //   For p_a: total = base + eps1(pi,p_a) + 0
                //   For p_b: total = base + eps1(pi,p_b) + 1
                // The "+1" shift in zone M is exactly what creates the LR ordering!

                // Check: for shared permutations in zone M,
                // does eps1 change between p_a and p_b?
                let perm_set_a: BTreeSet<Vec<u8>> = src_a.iter().map(|s| s.perm.clone()).collect();
                let perm_set_b: BTreeSet<Vec<u8>> = src_b.iter().map(|s| s.perm.clone()).collect();
                let shared: BTreeSet<Vec<u8>> =
                    perm_set_a.intersection(&perm_set_b).cloned().collect();

                let mut shared_in_m = 0u32;
                let mut shared_in_m_eps1_same = 0u32;
                let mut shared_in_m_eps1_diff = 0u32;

                for pi in &shared {
                    let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                    let z = zone_of(q, p_a); // zone is same for both
                    if z == 1 {
                        shared_in_m += 1;
                        let e1a = epsilon1(pi, p_a);
                        let e1b = epsilon1(pi, p_b);
                        if e1a == e1b {
                            shared_in_m_eps1_same += 1;
                        } else {
                            shared_in_m_eps1_diff += 1;
                        }
                    }
                }

                if shared_in_m > 0 {
                    println!(
                        "  Zone M shared: {} total, {} eps1 same, {} eps1 diff",
                        shared_in_m, shared_in_m_eps1_same, shared_in_m_eps1_diff
                    );
                }

                println!();
            }
        }
    }

    // =====================================================================
    // PART 3: Check if injection exists
    // =====================================================================
    println!("=== PART 3: Injection check for LR inequality ===\n");
    println!("For each (n,S,p_a,p_b,k), need A_k*B_{{k+1}} >= A_{{k+1}}*B_k.");
    println!("Equivalently: inject from 'bad' pairs to 'good' pairs.\n");

    for n in 4..=max_n.min(7) {
        let perms = all_permutations(n);

        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            for w in vp.windows(2) {
                let p_a = w[0];
                let p_b = w[1];

                // Collect (sigma, swaps) for each position
                let mut sigma_a: Vec<(&Vec<u8>, usize)> = Vec::new();
                let mut sigma_b: Vec<(&Vec<u8>, usize)> = Vec::new();
                for sig in class {
                    let pos = sig.iter().position(|&v| v == n).unwrap() as u8 + 1;
                    let sw = compute(sig, Stat::Swaps);
                    if pos == p_a {
                        sigma_a.push((sig, sw));
                    } else if pos == p_b {
                        sigma_b.push((sig, sw));
                    }
                }

                // Compute A and B polynomials
                let a_vals: Vec<usize> = sigma_a.iter().map(|&(_, s)| s).collect();
                let b_vals: Vec<usize> = sigma_b.iter().map(|&(_, s)| s).collect();

                if a_vals.is_empty() || b_vals.is_empty() {
                    continue;
                }

                let a_poly = build_poly(&a_vals);
                let b_poly = build_poly(&b_vals);

                let s_str = descent_set_to_string(mask, n);
                let max_k = a_poly.len().max(b_poly.len());

                // Check each k for LR and report margins
                let mut all_ok = true;
                let mut details = Vec::new();
                for k in 0..max_k {
                    let ak = coeff(&a_poly, k);
                    let bk1 = coeff(&b_poly, k + 1);
                    let ak1 = coeff(&a_poly, k + 1);
                    let bk = coeff(&b_poly, k);
                    let lhs = ak * bk1;
                    let rhs = ak1 * bk;
                    let margin = lhs - rhs;
                    if margin < 0 {
                        all_ok = false;
                    }
                    if lhs != 0 || rhs != 0 {
                        details.push((k, ak, ak1, bk, bk1, margin));
                    }
                }

                if !all_ok || n <= 5 {
                    println!(
                        "n={} S={} ({},{}): A={:?} B={:?}",
                        n, s_str, p_a, p_b, a_poly, b_poly
                    );
                    for &(k, ak, ak1, bk, bk1, margin) in &details {
                        println!(
                            "  k={}: A_k={} B_k+1={} - A_k+1={} B_k={} = {} {}",
                            k,
                            ak,
                            bk1,
                            ak1,
                            bk,
                            margin,
                            if margin < 0 { " <<< FAIL" } else { "" }
                        );
                    }
                    println!();
                }
            }
        }
    }

    // =====================================================================
    // PART 4: Summary statistics on source set relationships
    // =====================================================================
    println!("=== PART 4: Source set relationship summary ===\n");

    for n in 4..=max_n {
        let perms_prev = all_permutations(n - 1);
        let perms = all_permutations(n);

        let mut prev_by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        let mut n_total = 0u32;
        let mut n_same_asc = 0u32;
        let mut n_same_desc = 0u32;
        let mut n_same_both = 0u32;
        let mut n_a_subset_b = 0u32; // sources(p_a) subset of sources(p_b)
        let mut n_b_subset_a = 0u32;
        let mut n_disjoint = 0u32;

        for (&mask, _class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            for w in vp.windows(2) {
                let p_a = w[0];
                let p_b = w[1];
                n_total += 1;

                let sp_a_asc = source_asc(mask, p_a, n);
                let sp_a_desc = source_desc(mask, p_a, n);
                let sp_b_asc = source_asc(mask, p_b, n);
                let sp_b_desc = source_desc(mask, p_b, n);

                if sp_a_asc == sp_b_asc {
                    n_same_asc += 1;
                }
                if sp_a_desc == sp_b_desc {
                    n_same_desc += 1;
                }
                if sp_a_asc == sp_b_asc && sp_a_desc == sp_b_desc {
                    n_same_both += 1;
                }

                // Check permutation-level containment
                let mut gather = |p: u8| -> BTreeSet<Vec<u8>> {
                    let sp_a_val = source_asc(mask, p, n);
                    let sp_d_val = source_desc(mask, p, n);
                    let mut set = BTreeSet::new();
                    if let Some(cls) = prev_by_des.get(&sp_a_val) {
                        for pi in cls {
                            set.insert((*pi).clone());
                        }
                    }
                    if let Some(sp_d) = sp_d_val {
                        if sp_d != sp_a_val {
                            if let Some(cls) = prev_by_des.get(&sp_d) {
                                for pi in cls {
                                    set.insert((*pi).clone());
                                }
                            }
                        }
                    }
                    set
                };

                let set_a = gather(p_a);
                let set_b = gather(p_b);

                if set_a.is_subset(&set_b) {
                    n_a_subset_b += 1;
                }
                if set_b.is_subset(&set_a) {
                    n_b_subset_a += 1;
                }
                if set_a.is_disjoint(&set_b) {
                    n_disjoint += 1;
                }
            }
        }

        println!(
            "n={}: {} pairs | same_asc={} same_desc={} same_both={} | A<=B={} B<=A={} disjoint={}",
            n,
            n_total,
            n_same_asc,
            n_same_desc,
            n_same_both,
            n_a_subset_b,
            n_b_subset_a,
            n_disjoint
        );
    }

    println!("\nDone.");
}
