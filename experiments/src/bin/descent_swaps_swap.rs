/// Investigate the "source swap" operation for TN_2 proof of the
/// position-refined swaps coefficient matrix.
///
/// For consecutive valid positions p_a < p_b in P(S):
/// - sigma has n at position p_a, tau has n at position p_b
/// - Remove n from sigma to get pi_a, remove n from tau to get pi_b
/// - "Source swap": sigma' = insert(n @ p_a in pi_b), tau' = insert(n @ p_b in pi_a)
///
/// Questions:
/// 1. What are the source descent sets S'_{p_a}, S''_{p_a}, S'_{p_b}, S''_{p_b}?
/// 2. Where exactly do they differ?
/// 3. For each source pi_a (for p_a), how does Des(pi_a) differ from the valid sources for p_b?
/// 4. If we insert n at p_b into pi_a (ignoring descent constraints), what Des does the result have?
/// 5. Does the source swap preserve total swaps: swaps(sigma) + swaps(tau) = swaps(sigma') + swaps(tau')?
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

/// Insert value n at 1-indexed position p into pi in S_{n-1}.
/// Returns sigma in S_n.
fn insert(pi: &[u8], p: usize) -> Vec<u8> {
    let n = pi.len() as u8 + 1;
    let mut sigma = Vec::with_capacity(n as usize);
    for i in 0..p - 1 {
        sigma.push(pi[i]);
    }
    sigma.push(n);
    for i in p - 1..pi.len() {
        sigma.push(pi[i]);
    }
    sigma
}

/// Remove the value n from sigma in S_n, returning pi in S_{n-1}
/// and the 1-indexed position where n was.
fn remove_max(sigma: &[u8]) -> (Vec<u8>, usize) {
    let n = sigma.len() as u8;
    let pos = sigma.iter().position(|&v| v == n).unwrap();
    let mut pi = Vec::with_capacity(sigma.len() - 1);
    for (i, &v) in sigma.iter().enumerate() {
        if i != pos {
            pi.push(v);
        }
    }
    (pi, pos + 1) // 1-indexed position
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

fn mask_bits(mask: u64, n: u8) -> Vec<u8> {
    (1..n).filter(|&i| (mask >> (i - 1)) & 1 == 1).collect()
}

fn mask_diff_positions(a: u64, b: u64, n: u8) -> Vec<(u8, bool, bool)> {
    // Returns positions where a and b differ, with (pos, in_a, in_b)
    let mut diffs = Vec::new();
    for i in 1..n {
        let in_a = (a >> (i - 1)) & 1 == 1;
        let in_b = (b >> (i - 1)) & 1 == 1;
        if in_a != in_b {
            diffs.push((i, in_a, in_b));
        }
    }
    diffs
}

fn coeff(poly: &[i64], k: usize) -> i64 {
    if k < poly.len() { poly[k] } else { 0 }
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    // =========================================================================
    // PART 1: Source descent set comparison for consecutive valid positions
    // =========================================================================
    println!("=== PART 1: Source Descent Sets for Consecutive Valid Positions ===\n");

    for n in 5..=max_n.min(6) {
        let perms = all_permutations(n);
        let perms_prev = all_permutations(n - 1);

        let mut by_descent: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        let mut prev_by_des: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        println!("--- n = {} ---\n", n);

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            let s_str = descent_set_to_string(mask, n);
            println!("S = {} (n={}), valid positions P(S) = {:?}", s_str, n, vp);

            for w in vp.windows(2) {
                let pa = w[0];
                let pb = w[1];

                let spa_asc = source_asc(mask, pa, n);
                let spa_desc = source_desc(mask, pa, n);
                let spb_asc = source_asc(mask, pb, n);
                let spb_desc = source_desc(mask, pb, n);

                println!("\n  Consecutive pair: p_a={}, p_b={}", pa, pb);
                println!(
                    "    S'_asc(p_a={}) = {}",
                    pa,
                    descent_set_to_string(spa_asc, n - 1)
                );
                if let Some(sd) = spa_desc {
                    println!(
                        "    S''_desc(p_a={}) = {}",
                        pa,
                        descent_set_to_string(sd, n - 1)
                    );
                } else {
                    println!("    S''_desc(p_a={}) = None", pa);
                }
                println!(
                    "    S'_asc(p_b={}) = {}",
                    pb,
                    descent_set_to_string(spb_asc, n - 1)
                );
                if let Some(sd) = spb_desc {
                    println!(
                        "    S''_desc(p_b={}) = {}",
                        pb,
                        descent_set_to_string(sd, n - 1)
                    );
                } else {
                    println!("    S''_desc(p_b={}) = None", pb);
                }

                // Differences in source sets
                let diff_asc = mask_diff_positions(spa_asc, spb_asc, n - 1);
                println!("    Diff in S'_asc: {:?}", diff_asc);

                if let (Some(da), Some(db)) = (spa_desc, spb_desc) {
                    let diff_desc = mask_diff_positions(da, db, n - 1);
                    println!("    Diff in S''_desc: {:?}", diff_desc);
                }

                // Source permutation sets
                let mut src_perms_a: BTreeSet<Vec<u8>> = BTreeSet::new();
                if let Some(cls) = prev_by_des.get(&spa_asc) {
                    for pi in cls {
                        src_perms_a.insert(pi.clone());
                    }
                }
                if let Some(sd) = spa_desc {
                    if sd != spa_asc {
                        if let Some(cls) = prev_by_des.get(&sd) {
                            for pi in cls {
                                src_perms_a.insert(pi.clone());
                            }
                        }
                    }
                }

                let mut src_perms_b: BTreeSet<Vec<u8>> = BTreeSet::new();
                if let Some(cls) = prev_by_des.get(&spb_asc) {
                    for pi in cls {
                        src_perms_b.insert(pi.clone());
                    }
                }
                if let Some(sd) = spb_desc {
                    if sd != spb_asc {
                        if let Some(cls) = prev_by_des.get(&sd) {
                            for pi in cls {
                                src_perms_b.insert(pi.clone());
                            }
                        }
                    }
                }

                let shared: BTreeSet<Vec<u8>> =
                    src_perms_a.intersection(&src_perms_b).cloned().collect();
                let excl_a: BTreeSet<Vec<u8>> =
                    src_perms_a.difference(&src_perms_b).cloned().collect();
                let excl_b: BTreeSet<Vec<u8>> =
                    src_perms_b.difference(&src_perms_a).cloned().collect();

                println!(
                    "    Source perms: |A|={}, |B|={}, shared={}, excl_a={}, excl_b={}",
                    src_perms_a.len(),
                    src_perms_b.len(),
                    shared.len(),
                    excl_a.len(),
                    excl_b.len()
                );
            }
            println!();
        }
    }

    // =========================================================================
    // PART 2: Source swap analysis
    // For each (sigma, tau) in the "bad" set, compute the source swap and check:
    //   (a) What descent set does sigma' = insert(n @ p_a in pi_b) have?
    //   (b) What descent set does tau' = insert(n @ p_b in pi_a) have?
    //   (c) Is swaps(sigma) + swaps(tau) = swaps(sigma') + swaps(tau')?
    // =========================================================================
    println!("\n=== PART 2: Source Swap Operation Analysis ===\n");

    for n in 5..=max_n.min(7) {
        let perms = all_permutations(n);

        let mut by_descent: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        println!("--- n = {} ---\n", n);

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            let s_str = descent_set_to_string(mask, n);

            for w in vp.windows(2) {
                let pa = w[0];
                let pb = w[1];

                // Collect permutations at each position
                let sigmas_a: Vec<&Vec<u8>> = class
                    .iter()
                    .filter(|sig| sig.iter().position(|&v| v == n).unwrap() + 1 == pa as usize)
                    .collect();
                let sigmas_b: Vec<&Vec<u8>> = class
                    .iter()
                    .filter(|sig| sig.iter().position(|&v| v == n).unwrap() + 1 == pb as usize)
                    .collect();

                if sigmas_a.is_empty() || sigmas_b.is_empty() {
                    continue;
                }

                // For the "bad" set at each k: (sigma, tau) with swaps(sigma)=k+1, swaps(tau)=k
                // For the "good" set: (sigma', tau') with swaps(sigma')=k, swaps(tau')=k+1

                // First, let's check the source swap for ALL pairs, not just bad/good
                // For each (sigma, tau) with sigma^{-1}(n)=pa, tau^{-1}(n)=pb:
                //   Remove n: pi_a = sigma \ n,  pi_b = tau \ n
                //   sigma' = insert(n @ pa, pi_b), tau' = insert(n @ pb, pi_a)
                //   Check: Des(sigma') =? S, Des(tau') =? S
                //   Check: swaps(sigma) + swaps(tau) =? swaps(sigma') + swaps(tau')

                let mut total_pairs = 0u64;
                let mut des_sigma_prime_ok = 0u64;
                let mut des_tau_prime_ok = 0u64;
                let mut both_des_ok = 0u64;
                let mut total_swaps_preserved = 0u64;
                let mut swaps_diff_counts: BTreeMap<i32, u64> = BTreeMap::new();

                // Also track descent set differences
                let mut des_sigma_prime_diffs: BTreeMap<u64, u64> = BTreeMap::new();
                let mut des_tau_prime_diffs: BTreeMap<u64, u64> = BTreeMap::new();

                // Limit detailed output
                let show_detail = n <= 5;
                let mut detail_count = 0u64;

                for sigma in &sigmas_a {
                    for tau in &sigmas_b {
                        total_pairs += 1;
                        let (pi_a, _pos_a) = remove_max(sigma);
                        let (pi_b, _pos_b) = remove_max(tau);

                        // Source swap
                        let sigma_prime = insert(&pi_b, pa as usize);
                        let tau_prime = insert(&pi_a, pb as usize);

                        let des_sp = descent_set_bitmask(&sigma_prime);
                        let des_tp = descent_set_bitmask(&tau_prime);

                        if des_sp == mask {
                            des_sigma_prime_ok += 1;
                        }
                        if des_tp == mask {
                            des_tau_prime_ok += 1;
                        }
                        if des_sp == mask && des_tp == mask {
                            both_des_ok += 1;
                        }

                        *des_sigma_prime_diffs.entry(des_sp).or_insert(0) += 1;
                        *des_tau_prime_diffs.entry(des_tp).or_insert(0) += 1;

                        let sw_orig =
                            compute(sigma, Stat::Swaps) + compute(tau, Stat::Swaps);
                        let sw_new =
                            compute(&sigma_prime, Stat::Swaps) + compute(&tau_prime, Stat::Swaps);
                        let diff = sw_new as i32 - sw_orig as i32;

                        if diff == 0 {
                            total_swaps_preserved += 1;
                        }
                        *swaps_diff_counts.entry(diff).or_insert(0) += 1;

                        if show_detail && detail_count < 8 {
                            println!(
                                "  S={} pa={} pb={}: sigma={:?} tau={:?}",
                                s_str, pa, pb, sigma, tau
                            );
                            println!(
                                "    pi_a={:?} (Des={}) pi_b={:?} (Des={})",
                                pi_a,
                                descent_set_to_string(descent_set_bitmask(&pi_a), n - 1),
                                pi_b,
                                descent_set_to_string(descent_set_bitmask(&pi_b), n - 1),
                            );
                            println!(
                                "    sigma'={:?} (Des={}) tau'={:?} (Des={})",
                                sigma_prime,
                                descent_set_to_string(des_sp, n),
                                tau_prime,
                                descent_set_to_string(des_tp, n),
                            );
                            println!(
                                "    swaps: sigma={} tau={} sum={} | sigma'={} tau'={} sum={} | diff={}",
                                compute(sigma, Stat::Swaps),
                                compute(tau, Stat::Swaps),
                                sw_orig,
                                compute(&sigma_prime, Stat::Swaps),
                                compute(&tau_prime, Stat::Swaps),
                                sw_new,
                                diff,
                            );
                            println!(
                                "    Des match: sigma'={} tau'={}",
                                if des_sp == mask { "YES" } else { "NO" },
                                if des_tp == mask { "YES" } else { "NO" },
                            );
                            detail_count += 1;
                        }
                    }
                }

                println!(
                    "  S={} (pa={}, pb={}): {} total pairs",
                    s_str, pa, pb, total_pairs
                );
                println!(
                    "    Des(sigma')=S: {}/{} ({:.1}%)",
                    des_sigma_prime_ok,
                    total_pairs,
                    100.0 * des_sigma_prime_ok as f64 / total_pairs as f64
                );
                println!(
                    "    Des(tau')=S:   {}/{} ({:.1}%)",
                    des_tau_prime_ok,
                    total_pairs,
                    100.0 * des_tau_prime_ok as f64 / total_pairs as f64
                );
                println!(
                    "    Both Des=S:    {}/{} ({:.1}%)",
                    both_des_ok,
                    total_pairs,
                    100.0 * both_des_ok as f64 / total_pairs as f64
                );
                println!(
                    "    Total swaps preserved (diff=0): {}/{} ({:.1}%)",
                    total_swaps_preserved,
                    total_pairs,
                    100.0 * total_swaps_preserved as f64 / total_pairs as f64
                );
                println!("    Swaps diff distribution: {:?}", swaps_diff_counts);

                // Show distinct descent sets that appear
                if des_sigma_prime_diffs.len() <= 5 {
                    println!("    Des(sigma') values:");
                    for (&d, &cnt) in &des_sigma_prime_diffs {
                        println!(
                            "      {} : {} times",
                            descent_set_to_string(d, n),
                            cnt
                        );
                    }
                }
                if des_tau_prime_diffs.len() <= 5 {
                    println!("    Des(tau') values:");
                    for (&d, &cnt) in &des_tau_prime_diffs {
                        println!(
                            "      {} : {} times",
                            descent_set_to_string(d, n),
                            cnt
                        );
                    }
                }
                println!();
            }
        }
    }

    // =========================================================================
    // PART 3: Focus on pairs where both descent sets ARE preserved.
    // For those pairs, is the total swaps always preserved?
    // And does the swap restrict to a bijection on the "bad" set to "good" set?
    // =========================================================================
    println!("\n=== PART 3: Source Swap Restricted to Valid Pairs ===\n");
    println!("When both Des(sigma')=S AND Des(tau')=S, check swaps preservation.\n");

    for n in 5..=max_n {
        let perms = all_permutations(n);

        let mut by_descent: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        println!("--- n = {} ---", n);

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            let s_str = descent_set_to_string(mask, n);

            for w in vp.windows(2) {
                let pa = w[0];
                let pb = w[1];

                let sigmas_a: Vec<&Vec<u8>> = class
                    .iter()
                    .filter(|sig| sig.iter().position(|&v| v == n).unwrap() + 1 == pa as usize)
                    .collect();
                let sigmas_b: Vec<&Vec<u8>> = class
                    .iter()
                    .filter(|sig| sig.iter().position(|&v| v == n).unwrap() + 1 == pb as usize)
                    .collect();

                if sigmas_a.is_empty() || sigmas_b.is_empty() {
                    continue;
                }

                let mut valid_pairs = 0u64;
                let mut swaps_preserved_in_valid = 0u64;
                let mut valid_bad_to_good = 0u64; // bad pair maps to good pair
                let mut valid_good_to_bad = 0u64;
                let mut valid_same_k = 0u64;

                for sigma in &sigmas_a {
                    for tau in &sigmas_b {
                        let (pi_a, _) = remove_max(sigma);
                        let (pi_b, _) = remove_max(tau);

                        let sigma_prime = insert(&pi_b, pa as usize);
                        let tau_prime = insert(&pi_a, pb as usize);

                        let des_sp = descent_set_bitmask(&sigma_prime);
                        let des_tp = descent_set_bitmask(&tau_prime);

                        if des_sp != mask || des_tp != mask {
                            continue;
                        }

                        valid_pairs += 1;

                        let sw_s = compute(sigma, Stat::Swaps);
                        let sw_t = compute(tau, Stat::Swaps);
                        let sw_sp = compute(&sigma_prime, Stat::Swaps);
                        let sw_tp = compute(&tau_prime, Stat::Swaps);

                        if sw_s + sw_t == sw_sp + sw_tp {
                            swaps_preserved_in_valid += 1;
                        }

                        // Check if "bad" pair maps to "good" pair
                        // Bad: swaps(sigma) = k+1, swaps(tau) = k  (for some k)
                        // Good: swaps(sigma') = k, swaps(tau') = k+1
                        if sw_s > sw_t && sw_sp < sw_tp {
                            // sigma had more swaps, but sigma' has fewer
                            if sw_s == sw_t + 1 && sw_sp + 1 == sw_tp {
                                valid_bad_to_good += 1;
                            }
                        }
                        if sw_s < sw_t && sw_sp > sw_tp {
                            if sw_s + 1 == sw_t && sw_sp == sw_tp + 1 {
                                valid_good_to_bad += 1;
                            }
                        }
                        if sw_s == sw_t {
                            valid_same_k += 1;
                        }
                    }
                }

                if valid_pairs > 0 {
                    println!(
                        "  S={} (pa={}, pb={}): valid_pairs={} swaps_preserved={} bad->good={} good->bad={} same_k={}",
                        s_str, pa, pb, valid_pairs, swaps_preserved_in_valid,
                        valid_bad_to_good, valid_good_to_bad, valid_same_k
                    );
                }
            }
        }
        println!();
    }

    // =========================================================================
    // PART 4: Direct test - does the source swap preserve swaps ignoring descent?
    // For ALL pairs (sigma, tau), compute swaps(sigma') + swaps(tau') vs
    // swaps(sigma) + swaps(tau), even when Des is wrong.
    // =========================================================================
    println!("\n=== PART 4: Total Swaps Preservation (Ignoring Descent Constraints) ===\n");

    for n in 5..=max_n {
        let perms = all_permutations(n);

        let mut by_descent: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        let mut total_all = 0u64;
        let mut total_preserved = 0u64;

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            let s_str = descent_set_to_string(mask, n);

            for w in vp.windows(2) {
                let pa = w[0];
                let pb = w[1];

                let sigmas_a: Vec<&Vec<u8>> = class
                    .iter()
                    .filter(|sig| sig.iter().position(|&v| v == n).unwrap() + 1 == pa as usize)
                    .collect();
                let sigmas_b: Vec<&Vec<u8>> = class
                    .iter()
                    .filter(|sig| sig.iter().position(|&v| v == n).unwrap() + 1 == pb as usize)
                    .collect();

                if sigmas_a.is_empty() || sigmas_b.is_empty() {
                    continue;
                }

                let mut preserved = 0u64;
                let mut total = 0u64;
                let mut diff_dist: BTreeMap<i32, u64> = BTreeMap::new();

                for sigma in &sigmas_a {
                    for tau in &sigmas_b {
                        total += 1;
                        let (pi_a, _) = remove_max(sigma);
                        let (pi_b, _) = remove_max(tau);

                        let sigma_prime = insert(&pi_b, pa as usize);
                        let tau_prime = insert(&pi_a, pb as usize);

                        let sw_orig =
                            compute(sigma, Stat::Swaps) + compute(tau, Stat::Swaps);
                        let sw_new =
                            compute(&sigma_prime, Stat::Swaps) + compute(&tau_prime, Stat::Swaps);
                        let diff = sw_new as i32 - sw_orig as i32;

                        if diff == 0 {
                            preserved += 1;
                        }
                        *diff_dist.entry(diff).or_insert(0) += 1;
                    }
                }

                total_all += total;
                total_preserved += preserved;

                if n <= 6 || preserved != total {
                    println!(
                        "  n={} S={} (pa={}, pb={}): preserved {}/{} diff_dist={:?}",
                        n, s_str, pa, pb, preserved, total, diff_dist
                    );
                }
            }
        }

        println!(
            "  n={} TOTAL: preserved {}/{} ({:.1}%)\n",
            n,
            total_preserved,
            total_all,
            100.0 * total_preserved as f64 / total_all.max(1) as f64
        );
    }

    // =========================================================================
    // PART 5: Detailed analysis of what insert(n @ p_b, pi_a) produces.
    // For a source pi_a valid for p_a, inserting n at p_b gives tau'.
    // What is Des(tau')? How does it differ from S?
    // Key: the positions in [p_a, p_b-1] that are in S create forced descent
    // positions. When we insert at p_b instead of p_a, the descent at position
    // p_a is no longer forced, and the descent at p_b is now forced.
    // =========================================================================
    println!("\n=== PART 5: Detailed Descent Analysis of Cross-Insertions ===\n");

    for n in 5..=max_n.min(6) {
        let perms = all_permutations(n);
        let perms_prev = all_permutations(n - 1);

        let mut by_descent: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        let mut prev_by_des: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        println!("--- n = {} ---\n", n);

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            let s_str = descent_set_to_string(mask, n);

            for w in vp.windows(2) {
                let pa = w[0];
                let pb = w[1];

                // Get source permutations for p_a
                let spa_asc = source_asc(mask, pa, n);
                let spa_desc = source_desc(mask, pa, n);

                let mut sources_a: Vec<(Vec<u8>, char)> = Vec::new(); // (perm, type)
                if let Some(cls) = prev_by_des.get(&spa_asc) {
                    for pi in cls {
                        sources_a.push((pi.clone(), 'a'));
                    }
                }
                if let Some(sd) = spa_desc {
                    if sd != spa_asc {
                        if let Some(cls) = prev_by_des.get(&sd) {
                            for pi in cls {
                                sources_a.push((pi.clone(), 'd'));
                            }
                        }
                    }
                }

                if sources_a.is_empty() {
                    continue;
                }

                println!(
                    "S={} pa={} pb={}: inserting n={} at p_b={} into sources of p_a={}",
                    s_str, pa, pb, n, pb, pa
                );

                // For positions in [pa..pb-1], which are in S?
                let descent_window: Vec<u8> = (pa..pb)
                    .filter(|&i| (mask >> (i - 1)) & 1 == 1)
                    .map(|i| i as u8)
                    .collect();
                println!(
                    "  Descent positions in [{}, {}]: {:?}",
                    pa,
                    pb - 1,
                    descent_window
                );

                let mut show_count = 0;
                for (pi_a, src_type) in &sources_a {
                    if show_count >= 6 {
                        println!("  ... (truncated)");
                        break;
                    }

                    // Insert n at position p_b into pi_a
                    let tau_prime = insert(pi_a, pb as usize);
                    let des_tp = descent_set_bitmask(&tau_prime);

                    // How does des_tp differ from mask (= S)?
                    let diff_positions = mask_diff_positions(mask, des_tp, n);

                    // Also insert at p_a for comparison (this should give Des = S)
                    let sigma = insert(pi_a, pa as usize);
                    let des_s = descent_set_bitmask(&sigma);

                    println!(
                        "  pi_a={:?} (type={}, Des={})",
                        pi_a,
                        src_type,
                        descent_set_to_string(descent_set_bitmask(pi_a), n - 1)
                    );
                    println!(
                        "    insert@p_a={}: sigma={:?} Des={} swaps={}  [Des=S? {}]",
                        pa,
                        sigma,
                        descent_set_to_string(des_s, n),
                        compute(&sigma, Stat::Swaps),
                        if des_s == mask { "YES" } else { "NO" }
                    );
                    println!(
                        "    insert@p_b={}: tau'={:?} Des={} swaps={}  [Des=S? {}]",
                        pb,
                        tau_prime,
                        descent_set_to_string(des_tp, n),
                        compute(&tau_prime, Stat::Swaps),
                        if des_tp == mask { "YES" } else { "NO" }
                    );
                    if !diff_positions.is_empty() {
                        println!(
                            "    Des diff from S: {:?} (pos, in_S, in_Des(tau'))",
                            diff_positions
                        );
                    }

                    show_count += 1;
                }

                // Summary: what fraction of cross-insertions preserve the descent set?
                let total = sources_a.len();
                let preserved = sources_a
                    .iter()
                    .filter(|(pi_a, _)| {
                        let tau_prime = insert(pi_a, pb as usize);
                        descent_set_bitmask(&tau_prime) == mask
                    })
                    .count();
                println!(
                    "  Summary: {}/{} sources for p_a give Des=S when inserted at p_b\n",
                    preserved, total
                );
            }
        }
    }

    // =========================================================================
    // PART 6: Check the LR inequality directly and compute margins.
    // For each consecutive (pa, pb), for each k, compute A_k*B_{k+1} - A_{k+1}*B_k
    // where A is the polynomial for pa and B for pb.
    // =========================================================================
    println!("\n=== PART 6: LR Inequality Margins ===\n");

    for n in 5..=max_n {
        let perms = all_permutations(n);

        let mut by_descent: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        let mut all_ok = true;

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            let s_str = descent_set_to_string(mask, n);

            // Build polynomials for each position
            let mut polys: BTreeMap<u8, Vec<i64>> = BTreeMap::new();
            for &p in &vp {
                let vals: Vec<usize> = class
                    .iter()
                    .filter(|sig| sig.iter().position(|&v| v == n).unwrap() + 1 == p as usize)
                    .map(|sig| compute(sig, Stat::Swaps))
                    .collect();
                if !vals.is_empty() {
                    polys.insert(p, build_poly(&vals));
                }
            }

            // Check ALL pairs (not just consecutive) for TN_2
            let positions: Vec<u8> = polys.keys().copied().collect();
            for i in 0..positions.len() {
                for j in (i + 1)..positions.len() {
                    let pa = positions[i];
                    let pb = positions[j];
                    let a_poly = polys.get(&pa).unwrap();
                    let b_poly = polys.get(&pb).unwrap();

                    let max_k = a_poly.len().max(b_poly.len());
                    for k in 0..max_k {
                        let ak = coeff(a_poly, k);
                        let bk1 = coeff(b_poly, k + 1);
                        let ak1 = coeff(a_poly, k + 1);
                        let bk = coeff(b_poly, k);
                        let margin = ak * bk1 - ak1 * bk;
                        if margin < 0 {
                            all_ok = false;
                            println!(
                                "  FAIL n={} S={} pa={} pb={} k={}: {}*{} - {}*{} = {}",
                                n, s_str, pa, pb, k, ak, bk1, ak1, bk, margin
                            );
                        }
                    }
                }
            }
        }

        if all_ok {
            println!("  n={}: ALL TN_2 inequalities hold.", n);
        }
    }

    println!("\nDone.");
}
