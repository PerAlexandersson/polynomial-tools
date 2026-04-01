/// Decompose 2x2 minors of the L^{(p)} coefficient matrix by (q_a, q_b) zone pairs.
///
/// For consecutive valid positions p_a < p_b and k_1 < k_2:
///   Delta = A_{k1} B_{k2} - A_{k2} B_{k1} >= 0  (to prove)
///
/// where A_k = |{sigma: Des=S, pos(n)=p_a, swaps=k}|
///       B_k = |{sigma: Des=S, pos(n)=p_b, swaps=k}|
///
/// Each sigma comes from inserting n at position p into a source permutation pi.
/// The key decomposition uses N_p(s, q) = |{pi in source(p): swaps_tilde=s, pos(n-1)=q}|
/// and eps2(p, q) = 1 if q <= p-2, else 0.
///
/// Then A_k = sum_q N_a(k - e_a(q), q) and B_k = sum_q N_b(k - e_b(q), q),
/// and the minor decomposes by (q_a, q_b) pairs into zone contributions:
///
///   Zone L: q <= p_a - 2       => e_a = 1, e_b = 1
///   Zone M: p_a - 1 <= q <= p_b - 2  => e_a = 0, e_b = 1
///   Zone R: q >= p_b - 1       => e_a = 0, e_b = 0
///
use combpoly::fixed_descent::{
    permutations_grouped_by_descent_set_bitmask,
    q_refined_modified_source_distribution_from_grouped_source_permutations,
    valid_insertion_positions_for_target_descent_set as valid_positions,
};
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

/// eps2(p, q) = 1 if q <= p-2, else 0
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

/// Classify q into zone relative to (p_a, p_b):
///   Zone L: q <= p_a - 2  (both eps2 = 1)
///   Zone M: p_a - 1 <= q <= p_b - 2  (eps2(p_a,q)=0, eps2(p_b,q)=1)
///   Zone R: q >= p_b - 1  (both eps2 = 0)
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

/// Build the polynomial N_p(*, q) for a fixed q: coefficients in swaps_tilde.
fn build_n_poly(dist: &BTreeMap<(usize, u8), usize>, q: u8) -> Vec<i64> {
    let relevant: Vec<usize> = dist
        .iter()
        .filter(|&(&(_, qq), _)| qq == q)
        .map(|(&(s, _), &c)| {
            // Expand: c copies of value s
            std::iter::repeat(s).take(c).collect::<Vec<_>>()
        })
        .flatten()
        .collect();

    build_poly(&relevant)
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    println!("=== Minor decomposition by (q_a, q_b) zone pairs ===");
    println!("Testing n = 5 to {}\n", max_n);

    // Zone pair labels
    let zone_pairs: [(Zone, Zone); 9] = [
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

    // Global sign tracking per zone pair
    // For each zone pair: count of (positive, zero, negative) contributions
    let mut global_signs: BTreeMap<(Zone, Zone), (u64, u64, u64)> = BTreeMap::new();
    for &zp in &zone_pairs {
        global_signs.insert(zp, (0, 0, 0));
    }

    // Track whether full minor is always nonneg
    let mut total_minors = 0u64;
    let mut nonneg_minors = 0u64;

    // Track the worst negative zone-pair contributions
    let mut worst_negative: BTreeMap<(Zone, Zone), i64> = BTreeMap::new();
    for &zp in &zone_pairs {
        worst_negative.insert(zp, 0);
    }

    // Track cases where (R,M) or (M,R) is negative
    let mut rm_negative_examples: Vec<String> = Vec::new();
    let mut mr_negative_examples: Vec<String> = Vec::new();

    for n in 5..=max_n {
        let grouped_source_permutations = permutations_grouped_by_descent_set_bitmask(n - 1);

        // Also build full permutations of n for verification
        let perms_n = all_permutations(n);
        let mut by_descent_n: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_n {
            by_descent_n
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        let mut n_minors = 0u64;
        let mut n_nonneg = 0u64;

        for s_mask in 0u64..(1 << (n - 1)) {
            if s_mask & 1 != 0 {
                continue; // 1 not in S
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            // Build source distributions for each valid position
            let mut source_dists: BTreeMap<u8, BTreeMap<(usize, u8), usize>> = BTreeMap::new();
            for &p in &vp {
                let dist = q_refined_modified_source_distribution_from_grouped_source_permutations(
                    s_mask,
                    p,
                    n,
                    &grouped_source_permutations,
                )
                .counts_by_modified_source_swaps_and_previous_maximum_position;
                source_dists.insert(p, dist);
            }

            // Build L^{(p)} polynomials directly for verification
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

            // Also verify L^{(p)} from source distributions
            for &p in &vp {
                let dist = &source_dists[&p];
                let mut reconstructed = vec![0i64; 20];
                for (&(s, q), &count) in dist {
                    let k = s + eps2(p, q);
                    if k >= reconstructed.len() {
                        reconstructed.resize(k + 1, 0);
                    }
                    reconstructed[k] += count as i64;
                }
                while reconstructed.len() > 1 && *reconstructed.last().unwrap() == 0 {
                    reconstructed.pop();
                }

                if let Some(direct) = lp_polys.get(&p) {
                    if reconstructed != *direct {
                        println!(
                            "VERIFICATION FAIL: n={} S={} p={}: reconstructed {:?} != direct {:?}",
                            n,
                            descent_set_to_string(s_mask, n),
                            p,
                            reconstructed,
                            direct,
                        );
                    }
                }
            }

            // For each consecutive pair of valid positions
            for idx in 0..vp.len() - 1 {
                let p_a = vp[idx];
                let p_b = vp[idx + 1];

                let dist_a = &source_dists[&p_a];
                let dist_b = &source_dists[&p_b];

                // Find all q values that appear in either distribution
                let mut all_qs_a: Vec<u8> = dist_a.keys().map(|&(_, q)| q).collect();
                all_qs_a.sort();
                all_qs_a.dedup();
                let mut all_qs_b: Vec<u8> = dist_b.keys().map(|&(_, q)| q).collect();
                all_qs_b.sort();
                all_qs_b.dedup();

                // Build N_p(s, q) polynomials
                let n_polys_a: BTreeMap<u8, Vec<i64>> = all_qs_a
                    .iter()
                    .map(|&q| (q, build_n_poly(dist_a, q)))
                    .collect();
                let n_polys_b: BTreeMap<u8, Vec<i64>> = all_qs_b
                    .iter()
                    .map(|&q| (q, build_n_poly(dist_b, q)))
                    .collect();

                // Determine max degree
                let max_deg = lp_polys.values().map(|p| p.len()).max().unwrap_or(1);

                // For each (k_1 < k_2) pair, decompose the minor
                for k1 in 0..max_deg {
                    for k2 in (k1 + 1)..max_deg {
                        let a_k1 = lp_polys.get(&p_a).map(|p| coeff(p, k1)).unwrap_or(0);
                        let a_k2 = lp_polys.get(&p_a).map(|p| coeff(p, k2)).unwrap_or(0);
                        let b_k1 = lp_polys.get(&p_b).map(|p| coeff(p, k1)).unwrap_or(0);
                        let b_k2 = lp_polys.get(&p_b).map(|p| coeff(p, k2)).unwrap_or(0);

                        let minor = a_k1 * b_k2 - a_k2 * b_k1;

                        // Skip trivial minors
                        if a_k1 == 0 && a_k2 == 0 {
                            continue;
                        }
                        if b_k1 == 0 && b_k2 == 0 {
                            continue;
                        }

                        n_minors += 1;
                        if minor >= 0 {
                            n_nonneg += 1;
                        }

                        // Now decompose into zone-pair contributions
                        // Delta = sum_{q_a, q_b} [N_a(k1 - e_a(q_a), q_a) * N_b(k2 - e_b(q_b), q_b)
                        //                       - N_a(k2 - e_a(q_a), q_a) * N_b(k1 - e_b(q_b), q_b)]
                        let mut zone_contribs: BTreeMap<(Zone, Zone), i64> = BTreeMap::new();
                        for &zp in &zone_pairs {
                            zone_contribs.insert(zp, 0);
                        }

                        for &q_a in &all_qs_a {
                            let e_a = eps2(p_a, q_a);
                            let z_a = classify_zone(q_a, p_a, p_b);

                            for &q_b in &all_qs_b {
                                let e_b = eps2(p_b, q_b);
                                let z_b = classify_zone(q_b, p_a, p_b);

                                let zp = (z_a, z_b);

                                // N_a(k1 - e_a, q_a) * N_b(k2 - e_b, q_b)
                                // - N_a(k2 - e_a, q_a) * N_b(k1 - e_b, q_b)
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

                                let contrib = na_k1 * nb_k2 - na_k2 * nb_k1;
                                *zone_contribs.entry(zp).or_default() += contrib;
                            }
                        }

                        // Verify: sum of zone contributions = minor
                        let sum_contribs: i64 = zone_contribs.values().sum();
                        if sum_contribs != minor {
                            println!(
                                "DECOMPOSITION ERROR: n={} S={} p_a={} p_b={} k1={} k2={}: sum={} minor={}",
                                n,
                                descent_set_to_string(s_mask, n),
                                p_a, p_b, k1, k2, sum_contribs, minor,
                            );
                        }

                        // Record signs for each zone pair
                        for &zp in &zone_pairs {
                            let c = zone_contribs[&zp];
                            let (pos, zero, neg) = global_signs.get_mut(&zp).unwrap();
                            if c > 0 {
                                *pos += 1;
                            } else if c == 0 {
                                *zero += 1;
                            } else {
                                *neg += 1;
                            }

                            let w = worst_negative.get_mut(&zp).unwrap();
                            if c < *w {
                                *w = c;
                            }

                            // Track (R,M) and (M,R) negative examples
                            if c < 0 {
                                let example = format!(
                                    "n={} S={} p_a={} p_b={} k1={} k2={} contrib={}",
                                    n,
                                    descent_set_to_string(s_mask, n),
                                    p_a,
                                    p_b,
                                    k1,
                                    k2,
                                    c,
                                );
                                if zp == (Zone::R, Zone::M) && rm_negative_examples.len() < 5 {
                                    rm_negative_examples.push(example.clone());
                                }
                                if zp == (Zone::M, Zone::R) && mr_negative_examples.len() < 5 {
                                    mr_negative_examples.push(example);
                                }
                            }
                        }
                    }
                }
            }
        }

        total_minors += n_minors;
        nonneg_minors += n_nonneg;

        println!(
            "n={}: minors={}, nonneg={}/{}",
            n, n_minors, n_nonneg, n_minors,
        );
    }

    // Report zone-pair sign patterns
    println!("\n{}", "=".repeat(80));
    println!("ZONE-PAIR SIGN PATTERN SUMMARY (n=5..{})", max_n);
    println!("{}", "=".repeat(80));
    println!();
    println!("Zone definitions relative to (p_a, p_b) with p_a < p_b:");
    println!("  L: q <= p_a - 2       [both eps2 = 1]");
    println!("  M: p_a-1 <= q <= p_b-2  [eps2(p_a)=0, eps2(p_b)=1]");
    println!("  R: q >= p_b - 1       [both eps2 = 0]");
    println!();

    println!(
        "{:<10} {:>10} {:>10} {:>10}  {:>12}  {}",
        "Zone pair", "Positive", "Zero", "Negative", "Worst neg.", "Always >= 0?"
    );
    println!("{}", "-".repeat(75));

    for &zp in &zone_pairs {
        let (pos, zero, neg) = global_signs[&zp];
        let w = worst_negative[&zp];
        let always_nonneg = neg == 0;
        println!(
            "({},{})      {:>10} {:>10} {:>10}  {:>12}  {}",
            zone_label(zp.0),
            zone_label(zp.1),
            pos,
            zero,
            neg,
            w,
            if always_nonneg { "YES" } else { "NO" },
        );
    }

    println!();
    println!(
        "Total minors: {}, non-negative: {}/{}",
        total_minors, nonneg_minors, total_minors
    );

    // Report cross-zone negative examples
    if !rm_negative_examples.is_empty() {
        println!("\n--- (R,M) negative examples (first 5) ---");
        for ex in &rm_negative_examples {
            println!("  {}", ex);
        }
    }
    if !mr_negative_examples.is_empty() {
        println!("\n--- (M,R) negative examples (first 5) ---");
        for ex in &mr_negative_examples {
            println!("  {}", ex);
        }
    }

    // Additional analysis: within each zone pair, check if the contribution
    // can be expressed as a TN_2 minor
    println!("\n{}", "=".repeat(80));
    println!("DETAILED ZONE-PAIR ANALYSIS");
    println!("{}", "=".repeat(80));

    // For same-zone pairs (L,L), (M,M), (R,R), the eps2 shifts are the same.
    // So the contribution is sum_{q_a in Z_a, q_b in Z_b} [N_a(k1-e, q_a)*N_b(k2-e, q_b) - N_a(k2-e, q_a)*N_b(k1-e, q_b)]
    // = (sum_{q_a} N_a(k1-e, q_a)) * (sum_{q_b} N_b(k2-e, q_b)) - (sum N_a(k2-e, q_a)) * (sum N_b(k1-e, q_b))
    // This is a single 2x2 minor of the zone-restricted polynomials!
    println!();
    println!("For same-eps zone pairs (L,L), (R,R): both have same shift.");
    println!("The zone contribution is a SINGLE 2x2 minor of zone-restricted L^(p) polynomials.");
    println!("For (M,M): e_a=0, e_b=1, so the shift DIFFERS -- NOT a single minor.");
    println!();

    // Analyze whether same-shift zone pairs always give nonneg contributions
    // by computing zone-restricted L^(p) polynomials and checking their TN_2
    println!("Checking TN_2 of zone-restricted polynomials...\n");

    for n in 5..=max_n {
        let grouped_source_permutations = permutations_grouped_by_descent_set_bitmask(n - 1);

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
                let dist = q_refined_modified_source_distribution_from_grouped_source_permutations(
                    s_mask,
                    p,
                    n,
                    &grouped_source_permutations,
                )
                .counts_by_modified_source_swaps_and_previous_maximum_position;
                source_dists.insert(p, dist);
            }

            for idx in 0..vp.len() - 1 {
                let p_a = vp[idx];
                let p_b = vp[idx + 1];

                // Build zone-restricted L^(p) polys
                // For each zone Z, define L^(p)_Z(t) = sum_{q in Z} t^{eps2(p,q)} N_p(*, q)
                for zone in &[Zone::L, Zone::M, Zone::R] {
                    let mut la_zone = vec![0i64; 20];
                    let mut lb_zone = vec![0i64; 20];

                    let dist_a = &source_dists[&p_a];
                    let dist_b = &source_dists[&p_b];

                    for (&(s, q), &count) in dist_a {
                        if classify_zone(q, p_a, p_b) == *zone {
                            let k = s + eps2(p_a, q);
                            if k >= la_zone.len() {
                                la_zone.resize(k + 1, 0);
                            }
                            la_zone[k] += count as i64;
                        }
                    }
                    for (&(s, q), &count) in dist_b {
                        if classify_zone(q, p_a, p_b) == *zone {
                            let k = s + eps2(p_b, q);
                            if k >= lb_zone.len() {
                                lb_zone.resize(k + 1, 0);
                            }
                            lb_zone[k] += count as i64;
                        }
                    }

                    // Trim
                    while la_zone.len() > 1 && *la_zone.last().unwrap() == 0 {
                        la_zone.pop();
                    }
                    while lb_zone.len() > 1 && *lb_zone.last().unwrap() == 0 {
                        lb_zone.pop();
                    }

                    // Check TN_2 of the 2-row matrix [la_zone; lb_zone]
                    let max_deg = la_zone.len().max(lb_zone.len());
                    let mut zone_tn2 = true;
                    for k1 in 0..max_deg {
                        for k2 in (k1 + 1)..max_deg {
                            let minor = coeff(&la_zone, k1) * coeff(&lb_zone, k2)
                                - coeff(&la_zone, k2) * coeff(&lb_zone, k1);
                            if minor < 0 {
                                zone_tn2 = false;
                            }
                        }
                    }

                    if !zone_tn2 && n <= 7 {
                        println!(
                            "  Zone {} restricted polys NOT TN_2: n={} S={} p_a={} p_b={}",
                            zone_label(*zone),
                            n,
                            descent_set_to_string(s_mask, n),
                            p_a,
                            p_b,
                        );
                        println!("    L^(p_a)_{} = {:?}", zone_label(*zone), la_zone);
                        println!("    L^(p_b)_{} = {:?}", zone_label(*zone), lb_zone);
                    }
                }
            }
        }
    }

    // Final analysis: check individual (q_a, q_b) pair contributions
    // to understand finer structure
    println!("\n{}", "=".repeat(80));
    println!("FINE-GRAINED: Individual (q_a, q_b) contributions that are negative");
    println!("{}", "=".repeat(80));
    println!();

    let mut neg_individual_count = 0u64;
    let mut neg_individual_examples: Vec<String> = Vec::new();

    for n in 5..=max_n.min(6) {
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

            for idx in 0..vp.len() - 1 {
                let p_a = vp[idx];
                let p_b = vp[idx + 1];

                let dist_a = &source_dists[&p_a];
                let dist_b = &source_dists[&p_b];

                let mut all_qs_a: Vec<u8> = dist_a.keys().map(|&(_, q)| q).collect();
                all_qs_a.sort();
                all_qs_a.dedup();
                let mut all_qs_b: Vec<u8> = dist_b.keys().map(|&(_, q)| q).collect();
                all_qs_b.sort();
                all_qs_b.dedup();

                let n_polys_a: BTreeMap<u8, Vec<i64>> = all_qs_a
                    .iter()
                    .map(|&q| (q, build_n_poly(dist_a, q)))
                    .collect();
                let n_polys_b: BTreeMap<u8, Vec<i64>> = all_qs_b
                    .iter()
                    .map(|&q| (q, build_n_poly(dist_b, q)))
                    .collect();

                let max_deg = lp_polys.values().map(|p| p.len()).max().unwrap_or(1);

                for k1 in 0..max_deg {
                    for k2 in (k1 + 1)..max_deg {
                        // Look at individual (q_a, q_b) contributions
                        for &q_a in &all_qs_a {
                            let e_a = eps2(p_a, q_a);
                            for &q_b in &all_qs_b {
                                let e_b = eps2(p_b, q_b);

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

                                let contrib = na_k1 * nb_k2 - na_k2 * nb_k1;

                                if contrib < 0 {
                                    neg_individual_count += 1;
                                    let z_a = classify_zone(q_a, p_a, p_b);
                                    let z_b = classify_zone(q_b, p_a, p_b);
                                    if neg_individual_examples.len() < 20 {
                                        neg_individual_examples.push(format!(
                                            "n={} S={} pa={} pb={} k1={} k2={} qa={} qb={} zone=({},{}) e_a={} e_b={} c={}  Na={:?} Nb={:?}",
                                            n,
                                            descent_set_to_string(s_mask, n),
                                            p_a, p_b, k1, k2, q_a, q_b,
                                            zone_label(z_a), zone_label(z_b),
                                            e_a, e_b, contrib,
                                            na_poly, nb_poly,
                                        ));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    println!(
        "Total negative individual (q_a, q_b) contributions: {}",
        neg_individual_count
    );
    if !neg_individual_examples.is_empty() {
        println!("Examples:");
        for ex in &neg_individual_examples {
            println!("  {}", ex);
        }
    }

    // Summary of what we learned
    println!("\n{}", "=".repeat(80));
    println!("CONCLUSIONS");
    println!("{}", "=".repeat(80));

    let always_nonneg_zones: Vec<String> = zone_pairs
        .iter()
        .filter(|zp| global_signs[zp].2 == 0)
        .map(|zp| format!("({},{})", zone_label(zp.0), zone_label(zp.1)))
        .collect();
    let sometimes_neg_zones: Vec<String> = zone_pairs
        .iter()
        .filter(|zp| global_signs[zp].2 > 0)
        .map(|zp| format!("({},{})", zone_label(zp.0), zone_label(zp.1)))
        .collect();

    println!(
        "Zone pairs that are ALWAYS non-negative: {:?}",
        always_nonneg_zones
    );
    println!(
        "Zone pairs that can be NEGATIVE:         {:?}",
        sometimes_neg_zones
    );

    if nonneg_minors == total_minors {
        println!("\nAll full minors are non-negative (TN_2 verified).");
        if !sometimes_neg_zones.is_empty() {
            println!("BUT some zone-pair contributions are negative -- they cancel!");
            println!("This means the proof CANNOT proceed by showing each zone pair is nonneg.");
            println!("Instead, compensation between zone pairs is needed.");
        }
    }
}
