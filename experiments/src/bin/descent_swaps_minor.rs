/// Decompose 2x2 minors of the coefficient matrix of {L_{n,S}^{(p)}} into
/// contributions from individual permutations.
///
/// For each S with 1 not in S (n = 5..max_n):
///   - Build L^{(p)}(t) = sum_{sigma: Des(sigma)=S, pos(n,sigma)=p} t^{swaps(sigma)}
///   - Form the coefficient matrix M with M_{p,k} = [t^k] L^{(p)}
///   - For each 2x2 minor (rows p_a < p_b, cols k_1 < k_2):
///       minor = A_{k1}*B_{k2} - A_{k2}*B_{k1}
///     where A_k = #{sigma: Des=S, pos(n)=p_a, swaps=k}
///           B_k = #{sigma: Des=S, pos(n)=p_b, swaps=k}
///
/// Tests:
/// 1. Verify all minors are >= 0 (TN_2)
/// 2. Check ratio monotonicity: A_k/A_{k+1} >= B_k/B_{k+1} for p_a < p_b
///    (equivalent to log-concavity shifting / stochastic dominance)
/// 3. Look at the structure of the minor decomposition to find injection patterns
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::format_poly;
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

/// Compute the position of value `val` (1-indexed) in permutation pi.
fn pos_of(pi: &[u8], val: u8) -> u8 {
    pi.iter().position(|&v| v == val).unwrap() as u8 + 1
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    println!("=== 2x2 Minor Decomposition for L^{{(p)}} Coefficient Matrix ===");
    println!("Checking n = 5 to {}\n", max_n);

    let mut total_minors = 0u64;
    let mut nonneg_minors = 0u64;
    let mut total_families = 0u64;
    let mut tn2_families = 0u64;

    // Track ratio monotonicity
    let mut ratio_mono_pairs_total = 0u64;
    let mut ratio_mono_pairs_ok = 0u64;

    // Track minor value distribution
    let mut zero_minors = 0u64;
    let mut positive_minors = 0u64;

    // Track injection patterns
    let mut injection_analysis_count = 0u64;

    for n in 5..=max_n {
        let perms = all_permutations(n);

        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        let mut n_total_minors = 0u64;
        let mut n_nonneg = 0u64;
        let mut n_families = 0u64;
        let mut n_tn2 = 0u64;
        let mut n_ratio_total = 0u64;
        let mut n_ratio_ok = 0u64;
        let mut n_zero = 0u64;
        let mut n_positive = 0u64;

        for (&mask, class) in &by_descent {
            // Only consider S with 1 not in S
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            // Build L^{(p)} polynomials and store per-permutation data
            // For each valid position p, collect swaps values of all permutations with pos(n)=p
            let mut lp_polys: Vec<(u8, Vec<i64>)> = Vec::new();
            let mut lp_perms: BTreeMap<u8, Vec<(&Vec<u8>, usize)>> = BTreeMap::new(); // p -> [(perm, swaps)]

            for &p in &vp {
                let mut vals: Vec<usize> = Vec::new();
                let mut perm_data: Vec<(&Vec<u8>, usize)> = Vec::new();
                for pi in class.iter() {
                    if pos_of(pi, n) == p {
                        let sw = compute(pi, Stat::Swaps);
                        vals.push(sw);
                        perm_data.push((pi, sw));
                    }
                }
                if !vals.is_empty() {
                    lp_polys.push((p, build_poly(&vals)));
                    lp_perms.insert(p, perm_data);
                }
            }

            if lp_polys.len() < 2 {
                continue;
            }

            n_families += 1;
            let s_str = descent_set_to_string(mask, n);

            // Check TP_2 and decompose each minor
            let max_deg = lp_polys
                .iter()
                .map(|(_, poly)| poly.len())
                .max()
                .unwrap_or(0);
            let mut family_tn2 = true;
            let mut family_minors: Vec<(u8, u8, usize, usize, i64)> = Vec::new();

            for i in 0..lp_polys.len() {
                for j in (i + 1)..lp_polys.len() {
                    let (p_a, ref poly_a) = lp_polys[i];
                    let (p_b, ref poly_b) = lp_polys[j];

                    for k1 in 0..max_deg {
                        for k2 in (k1 + 1)..max_deg {
                            let a_k1 = coeff(poly_a, k1);
                            let a_k2 = coeff(poly_a, k2);
                            let b_k1 = coeff(poly_b, k1);
                            let b_k2 = coeff(poly_b, k2);

                            let minor = a_k1 * b_k2 - a_k2 * b_k1;
                            n_total_minors += 1;

                            if minor >= 0 {
                                n_nonneg += 1;
                                if minor == 0 {
                                    n_zero += 1;
                                } else {
                                    n_positive += 1;
                                }
                            } else {
                                family_tn2 = false;
                            }

                            family_minors.push((p_a, p_b, k1, k2, minor));
                        }
                    }
                }
            }

            if family_tn2 {
                n_tn2 += 1;
            } else {
                println!("  TN_2 FAIL: n={} S={} positions={:?}", n, s_str, vp);
                for (p, poly) in &lp_polys {
                    println!("    L^({}) = {}", p, format_poly(poly));
                }
                for &(p_a, p_b, k1, k2, minor) in &family_minors {
                    if minor < 0 {
                        println!(
                            "    Negative minor: rows ({},{}) cols ({},{}) = {}",
                            p_a, p_b, k1, k2, minor
                        );
                    }
                }
            }

            // ===== RATIO MONOTONICITY CHECK =====
            // For each pair p_a < p_b, check: A_k/A_{k+1} >= B_k/B_{k+1}
            // Equivalently: A_k * B_{k+1} >= A_{k+1} * B_k
            // This is exactly the consecutive-column minor >= 0 condition.
            // But let's check the full ratio sequence.
            for i in 0..lp_polys.len() {
                for j in (i + 1)..lp_polys.len() {
                    let (p_a, ref poly_a) = lp_polys[i];
                    let (p_b, ref poly_b) = lp_polys[j];

                    // For ratio analysis, look at consecutive degrees where both are nonzero
                    let max_k = poly_a.len().max(poly_b.len());
                    let mut ratios_a: Vec<(usize, f64)> = Vec::new(); // (k, A_k/A_{k+1})
                    let mut ratios_b: Vec<(usize, f64)> = Vec::new();

                    for k in 0..max_k.saturating_sub(1) {
                        let a_k = coeff(poly_a, k);
                        let a_k1 = coeff(poly_a, k + 1);
                        let b_k = coeff(poly_b, k);
                        let b_k1 = coeff(poly_b, k + 1);

                        if a_k1 > 0 {
                            ratios_a.push((k, a_k as f64 / a_k1 as f64));
                        }
                        if b_k1 > 0 {
                            ratios_b.push((k, b_k as f64 / b_k1 as f64));
                        }

                        // Check the consecutive-column minor (this IS ratio monotonicity)
                        if a_k1 > 0 && b_k1 > 0 {
                            n_ratio_total += 1;
                            if a_k * b_k1 >= a_k1 * b_k {
                                n_ratio_ok += 1;
                            }
                        }
                    }
                }
            }

            // ===== INJECTION / STRUCTURAL ANALYSIS =====
            // For small cases, print detailed minor decomposition
            if n <= 6 && lp_polys.len() >= 2 {
                injection_analysis_count += 1;
                println!("\n  --- Detailed analysis: n={} S={} ---", n, s_str);
                for (p, poly) in &lp_polys {
                    println!(
                        "    L^({}) = {} (total={})",
                        p,
                        format_poly(poly),
                        poly.iter().sum::<i64>()
                    );
                }

                // For each pair of positions, show the ratio sequence
                for i in 0..lp_polys.len() {
                    for j in (i + 1)..lp_polys.len() {
                        let (p_a, ref poly_a) = lp_polys[i];
                        let (p_b, ref poly_b) = lp_polys[j];

                        println!("    Pair (p_a={}, p_b={}):", p_a, p_b);

                        // Print coefficients side by side
                        let max_k = poly_a.len().max(poly_b.len());
                        print!("      degree:  ");
                        for k in 0..max_k {
                            print!("{:>6}", k);
                        }
                        println!();
                        print!("      A (p={}): ", p_a);
                        for k in 0..max_k {
                            print!("{:>6}", coeff(poly_a, k));
                        }
                        println!();
                        print!("      B (p={}): ", p_b);
                        for k in 0..max_k {
                            print!("{:>6}", coeff(poly_b, k));
                        }
                        println!();

                        // Print ratios A_k/A_{k+1} vs B_k/B_{k+1}
                        print!("      A_k/A_{{k+1}}: ");
                        for k in 0..max_k.saturating_sub(1) {
                            let a_k = coeff(poly_a, k);
                            let a_k1 = coeff(poly_a, k + 1);
                            if a_k1 > 0 {
                                print!("{:>8.3}", a_k as f64 / a_k1 as f64);
                            } else if a_k > 0 {
                                print!("{:>8}", "inf");
                            } else {
                                print!("{:>8}", "-");
                            }
                        }
                        println!();
                        print!("      B_k/B_{{k+1}}: ");
                        for k in 0..max_k.saturating_sub(1) {
                            let b_k = coeff(poly_b, k);
                            let b_k1 = coeff(poly_b, k + 1);
                            if b_k1 > 0 {
                                print!("{:>8.3}", b_k as f64 / b_k1 as f64);
                            } else if b_k > 0 {
                                print!("{:>8}", "inf");
                            } else {
                                print!("{:>8}", "-");
                            }
                        }
                        println!();

                        // Show all 2x2 minors for this pair
                        print!("      minors: ");
                        for k1 in 0..max_k {
                            for k2 in (k1 + 1)..max_k {
                                let minor = coeff(poly_a, k1) * coeff(poly_b, k2)
                                    - coeff(poly_a, k2) * coeff(poly_b, k1);
                                if minor != 0 {
                                    print!("({},{})={} ", k1, k2, minor);
                                }
                            }
                        }
                        println!();

                        // ===== INJECTION ANALYSIS =====
                        // Check: can we find a structural explanation?
                        // The minor A_{k1}*B_{k2} - A_{k2}*B_{k1} >= 0
                        // means: among perms with pos(n)=p_a having k1 swaps
                        //        paired with perms with pos(n)=p_b having k2 swaps,
                        //        there are at least as many pairs as the "swapped" pairing.
                        //
                        // Look at the distribution of swaps conditioned on pos(n):
                        // Is it true that the distribution for smaller p stochastically
                        // dominates (in reversed likelihood ratio order) the one for larger p?

                        // Check LIKELIHOOD RATIO ORDER:
                        // p(k|p_a)/p(k|p_b) should be non-increasing in k
                        // (where p(k|p) = A_k/sum(A) or B_k/sum(B))
                        let sum_a: f64 = (0..max_k).map(|k| coeff(poly_a, k) as f64).sum();
                        let sum_b: f64 = (0..max_k).map(|k| coeff(poly_b, k) as f64).sum();

                        if sum_a > 0.0 && sum_b > 0.0 {
                            print!("      LR ratio A/B: ");
                            let mut lr_ratios: Vec<(usize, f64)> = Vec::new();
                            for k in 0..max_k {
                                let a_norm = coeff(poly_a, k) as f64 / sum_a;
                                let b_norm = coeff(poly_b, k) as f64 / sum_b;
                                if b_norm > 0.0 {
                                    let lr = a_norm / b_norm;
                                    print!("{:>8.3}", lr);
                                    lr_ratios.push((k, lr));
                                } else if a_norm > 0.0 {
                                    print!("{:>8}", "inf");
                                } else {
                                    print!("{:>8}", "-");
                                }
                            }

                            // Check if LR ratios are non-increasing
                            let lr_monotone =
                                lr_ratios.windows(2).all(|w| w[0].1 >= w[1].1 - 1e-12);
                            print!("  [{}]", if lr_monotone { "MONOTONE" } else { "NOT mono" });
                            println!();
                        }

                        // ===== SOURCE DECOMPOSITION =====
                        // For each permutation in A_{k} (pos(n)=p_a, swaps=k),
                        // look at the source permutation (remove n) and its properties.
                        // Group by: (source descent set, pos(n-1) in source, eps1, eps2)
                        if let (Some(perms_a), Some(perms_b)) =
                            (lp_perms.get(&p_a), lp_perms.get(&p_b))
                        {
                            // Track (swaps, pos_of_n-1, pos_of_1) for injection search
                            let mut a_data: Vec<(usize, u8, u8, Vec<u8>)> = Vec::new(); // (swaps, pos(n-1), pos(1), reduced_perm)
                            for &(pi, sw) in perms_a.iter() {
                                let pos_nm1 = if n >= 2 { pos_of(pi, n - 1) } else { 0 };
                                let pos_1 = pos_of(pi, 1);
                                // Reduced permutation: remove n, standardize
                                let reduced: Vec<u8> =
                                    pi.iter().filter(|&&v| v != n).copied().collect();
                                a_data.push((sw, pos_nm1, pos_1, reduced));
                            }

                            let mut b_data: Vec<(usize, u8, u8, Vec<u8>)> = Vec::new();
                            for &(pi, sw) in perms_b.iter() {
                                let pos_nm1 = if n >= 2 { pos_of(pi, n - 1) } else { 0 };
                                let pos_1 = pos_of(pi, 1);
                                let reduced: Vec<u8> =
                                    pi.iter().filter(|&&v| v != n).copied().collect();
                                b_data.push((sw, pos_nm1, pos_1, reduced));
                            }

                            // Show distribution of pos(n-1) conditioned on swaps count
                            if n <= 6 {
                                println!("      pos(n-1) distribution by swaps:");
                                let max_sw = a_data
                                    .iter()
                                    .chain(b_data.iter())
                                    .map(|d| d.0)
                                    .max()
                                    .unwrap_or(0);
                                for sw in 0..=max_sw {
                                    let a_posns: Vec<u8> =
                                        a_data.iter().filter(|d| d.0 == sw).map(|d| d.1).collect();
                                    let b_posns: Vec<u8> =
                                        b_data.iter().filter(|d| d.0 == sw).map(|d| d.1).collect();
                                    if !a_posns.is_empty() || !b_posns.is_empty() {
                                        // Count by pos(n-1)
                                        let mut a_counts: BTreeMap<u8, usize> = BTreeMap::new();
                                        let mut b_counts: BTreeMap<u8, usize> = BTreeMap::new();
                                        for &q in &a_posns {
                                            *a_counts.entry(q).or_default() += 1;
                                        }
                                        for &q in &b_posns {
                                            *b_counts.entry(q).or_default() += 1;
                                        }
                                        println!(
                                            "        sw={}: A(p={}) {:?}  B(p={}) {:?}",
                                            sw, p_a, a_counts, p_b, b_counts
                                        );
                                    }
                                }

                                // Check: do shared source permutations contribute?
                                // Two perms sigma (pos(n)=p_a) and tau (pos(n)=p_b)
                                // share the same reduced perm iff removing n gives the same thing.
                                let mut shared_sources = 0usize;
                                for ad in &a_data {
                                    for bd in &b_data {
                                        if ad.3 == bd.3 {
                                            shared_sources += 1;
                                        }
                                    }
                                }
                                println!(
                                    "      Shared source pairs (same reduced perm): {}",
                                    shared_sources
                                );
                            }
                        }
                    }
                }
            }
        }

        total_minors += n_total_minors;
        nonneg_minors += n_nonneg;
        total_families += n_families;
        tn2_families += n_tn2;
        ratio_mono_pairs_total += n_ratio_total;
        ratio_mono_pairs_ok += n_ratio_ok;
        zero_minors += n_zero;
        positive_minors += n_positive;

        println!(
            "\nn={}: families={}, TN_2={}/{}, minors={} (pos={}, zero={}, neg={}), ratio_mono={}/{}",
            n,
            n_families,
            n_tn2,
            n_families,
            n_total_minors,
            n_positive,
            n_zero,
            n_total_minors - n_nonneg,
            n_ratio_ok,
            n_ratio_total,
        );
    }

    println!("\n=== SUMMARY (n=5..{}) ===", max_n);
    println!("Total families:           {}", total_families);
    println!(
        "TN_2 families:            {}/{}",
        tn2_families, total_families
    );
    println!("Total 2x2 minors:         {}", total_minors);
    println!(
        "Non-negative minors:      {}/{}",
        nonneg_minors, total_minors
    );
    println!("  Positive:               {}", positive_minors);
    println!("  Zero:                   {}", zero_minors);
    println!("  Negative:               {}", total_minors - nonneg_minors);
    println!(
        "Ratio monotonicity:       {}/{}",
        ratio_mono_pairs_ok, ratio_mono_pairs_total
    );

    // ===== STOCHASTIC DOMINANCE CHECK =====
    // Check: does p_a < p_b imply L^{(p_a)} stochastically dominates L^{(p_b)}
    // in the likelihood ratio order?
    //
    // LR order: P(X=k|p_a)/P(X=k|p_b) is non-increasing in k.
    // Equivalently: A_k * B_{k'} >= A_{k'} * B_k for all k < k'.
    // This is exactly TN_2!
    //
    // So TN_2 <=> LR order between the swap distributions conditioned on pos(n).

    println!("\n=== KEY INSIGHT ===");
    println!("TN_2 of the coefficient matrix is EQUIVALENT to:");
    println!("  The swap distribution conditioned on pos(n)=p_a");
    println!("  dominates the distribution conditioned on pos(n)=p_b");
    println!("  in the LIKELIHOOD RATIO ORDER, for all p_a < p_b.");
    println!();
    println!("Ratio monotonicity (A_k/A_{{k+1}} >= B_k/B_{{k+1}})");
    println!("is equivalent to consecutive-column TN_2 minors,");
    println!("which implies all TN_2 minors (since TN_2 = all 2x2 minors >= 0).");
    println!();
    if ratio_mono_pairs_ok == ratio_mono_pairs_total {
        println!("RESULT: Consecutive ratio monotonicity HOLDS for all tested cases.");
        println!("This means: smaller pos(n) => stochastically fewer swaps (in LR order).");
    } else {
        println!(
            "RESULT: Consecutive ratio monotonicity FAILS for {}/{} pairs.",
            ratio_mono_pairs_total - ratio_mono_pairs_ok,
            ratio_mono_pairs_total
        );
    }

    if tn2_families == total_families {
        println!(
            "\nTN_2 VERIFIED for all {} families across n=5..{}.",
            total_families, max_n
        );
    } else {
        println!(
            "\nTN_2 FAILS for {} out of {} families.",
            total_families - tn2_families,
            total_families
        );
    }

    // ===== EXTRA: Check if the normalized coefficient vectors
    //       show a clean pattern (e.g., componentwise dominance) =====
    println!("\n=== NORMALIZED COEFFICIENT ANALYSIS ===");
    println!(
        "Checking whether L^{{(p_a)}}/|L^{{(p_a)}}| >= L^{{(p_b)}}/|L^{{(p_b)}}| componentwise"
    );
    println!("for the lower-degree coefficients (stochastic dominance in CDF sense).\n");

    let mut cdf_dom_total = 0u64;
    let mut cdf_dom_ok = 0u64;

    for n in 5..=max_n {
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

            let mut lp_polys: Vec<(u8, Vec<i64>)> = Vec::new();
            for &p in &vp {
                let vals: Vec<usize> = class
                    .iter()
                    .filter(|pi| pos_of(pi, n) == p)
                    .map(|pi| compute(pi, Stat::Swaps))
                    .collect();
                if !vals.is_empty() {
                    lp_polys.push((p, build_poly(&vals)));
                }
            }

            if lp_polys.len() < 2 {
                continue;
            }

            // Check CDF dominance: For p_a < p_b,
            // sum_{k=0}^{K} A_k / sum(A) >= sum_{k=0}^{K} B_k / sum(B)  for all K
            // (i.e., smaller p has stochastically smaller swaps in the usual sense)
            for i in 0..lp_polys.len() {
                for j in (i + 1)..lp_polys.len() {
                    let (p_a, ref poly_a) = lp_polys[i];
                    let (p_b, ref poly_b) = lp_polys[j];

                    let sum_a: i64 = poly_a.iter().sum();
                    let sum_b: i64 = poly_b.iter().sum();
                    if sum_a == 0 || sum_b == 0 {
                        continue;
                    }

                    let max_k = poly_a.len().max(poly_b.len());
                    let mut cum_a = 0i64;
                    let mut cum_b = 0i64;
                    let mut cdf_ok = true;
                    for k in 0..max_k {
                        cum_a += coeff(poly_a, k);
                        cum_b += coeff(poly_b, k);
                        // Check: cum_a / sum_a >= cum_b / sum_b
                        // <=> cum_a * sum_b >= cum_b * sum_a
                        if cum_a * sum_b < cum_b * sum_a {
                            cdf_ok = false;
                            if n <= 7 {
                                let s_str = descent_set_to_string(mask, n);
                                println!(
                                    "  CDF dominance FAIL: n={} S={} p_a={} p_b={} at k={}",
                                    n, s_str, p_a, p_b, k
                                );
                                println!(
                                    "    cum_a/sum_a = {}/{} = {:.4}, cum_b/sum_b = {}/{} = {:.4}",
                                    cum_a,
                                    sum_a,
                                    cum_a as f64 / sum_a as f64,
                                    cum_b,
                                    sum_b,
                                    cum_b as f64 / sum_b as f64
                                );
                            }
                            break;
                        }
                    }
                    cdf_dom_total += 1;
                    if cdf_ok {
                        cdf_dom_ok += 1;
                    }
                }
            }
        }
    }

    println!(
        "First-order stochastic dominance: {}/{}",
        cdf_dom_ok, cdf_dom_total
    );
    if cdf_dom_ok == cdf_dom_total {
        println!("=> Smaller pos(n) gives stochastically fewer swaps (in CDF sense too).");
        println!("   Combined with LR order, this is a strong structural result.");
    }
}
