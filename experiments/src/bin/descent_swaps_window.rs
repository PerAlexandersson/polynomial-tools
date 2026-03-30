/// Investigate the WINDOW structure for the varying-source MTP_2 proof.
///
/// For consecutive valid positions p_a < p_b in P(S):
///   Source(p_a) = permutations with Des in {S'_{p_a}, S''_{p_a}}
///   Source(p_b) = permutations with Des in {S'_{p_b}, S''_{p_b}}
///
/// These descent sets agree outside a "window" W = [p_a-1, p_b-1] (1-indexed in [n-2]).
///
/// Key questions:
///   1. What is |W|? (Always 1 or 2? Or sometimes larger?)
///   2. Fix values outside W. How do inside-W arrangements relate between sources?
///   3. For each outside-W config, compute Count_a(s, q) and Count_b(s, q)
///      and look for a simple relationship.
///
/// If W has size 1, then the source sets are related by "toggling a descent at
/// a single position" -- which means there's a near-bijection between them.
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use std::collections::{BTreeMap, BTreeSet};

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

fn mask_to_set(mask: u64, n: u8) -> BTreeSet<u8> {
    let mut s = BTreeSet::new();
    for i in 1..n {
        if (mask >> (i - 1)) & 1 == 1 {
            s.insert(i);
        }
    }
    s
}

fn set_to_string(s: &BTreeSet<u8>) -> String {
    if s.is_empty() {
        return "{}".to_string();
    }
    let parts: Vec<String> = s.iter().map(|x| x.to_string()).collect();
    format!("{{{}}}", parts.join(","))
}

fn mask_to_string(mask: u64, n: u8) -> String {
    set_to_string(&mask_to_set(mask, n))
}

fn sym_diff(a: u64, b: u64) -> u64 {
    a ^ b
}

/// Compute the window W where source descent sets can differ.
/// Returns 1-indexed positions in [1..n-2] where S'_{p_a} and S'_{p_b} differ.
fn compute_window(s_mask: u64, p_a: u8, p_b: u8, n: u8) -> BTreeSet<u8> {
    let sp_a = source_asc(s_mask, p_a, n);
    let sp_b = source_asc(s_mask, p_b, n);
    let sd = sym_diff(sp_a, sp_b);
    let mut w = BTreeSet::new();
    for i in 1..n - 1 {
        if (sd >> (i - 1)) & 1 == 1 {
            w.insert(i);
        }
    }
    w
}

/// Compute FULL window: union of positions where S'_pa differs from S'_pb
/// AND where S''_pa differs from S''_pb.
fn compute_full_window(s_mask: u64, p_a: u8, p_b: u8, n: u8) -> BTreeSet<u8> {
    let sp_a = source_asc(s_mask, p_a, n);
    let sp_b = source_asc(s_mask, p_b, n);
    let sd1 = sym_diff(sp_a, sp_b);

    let spp_a = source_desc(s_mask, p_a, n);
    let spp_b = source_desc(s_mask, p_b, n);
    let sd2 = match (spp_a, spp_b) {
        (Some(a), Some(b)) => sym_diff(a, b),
        (Some(a), None) => a, // not quite right but won't matter; p=1 or p=n edge cases
        (None, Some(b)) => b,
        (None, None) => 0,
    };

    let sd = sd1 | sd2;
    let mut w = BTreeSet::new();
    for i in 1..n - 1 {
        if (sd >> (i - 1)) & 1 == 1 {
            w.insert(i);
        }
    }
    w
}

/// Modified swaps: swaps(pi) + eps1(pi, p) + eps2(p, q).
/// q = 1-indexed position of (n-1) in pi.
fn modified_swaps(pi: &[u8], p: u8) -> usize {
    let n = pi.len() as u8 + 1;
    let base = compute(pi, Stat::Swaps);
    let e1 = if epsilon1(pi, p) { 1 } else { 0 };
    let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
    let e2 = if p > 1 && q < p - 1 { 1 } else { 0 };
    base + e1 + e2
}

/// Position of value (n-1) in pi (1-indexed).
fn pos_of_nm1(pi: &[u8]) -> u8 {
    let n = pi.len() as u8 + 1;
    pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1
}

/// Extract the values at positions OUTSIDE the window (0-indexed positions).
/// Returns (outside_config, inside_values, inside_positions_0indexed)
fn split_by_window(
    pi: &[u8],
    window_positions_1indexed: &BTreeSet<u8>,
) -> (Vec<(u8, u8)>, Vec<(u8, u8)>) {
    // window_positions are 1-indexed descent-check positions.
    // Position j in [1..n-2] for a permutation of [n-1] corresponds to
    // the descent between pi[j-1] and pi[j] (0-indexed).
    // The "window positions" in the permutation are indices j-1 and j for each j in W.
    // But to keep it simple: we say position j (1-indexed, in the descent sense)
    // involves the permutation indices j-1 and j.

    // A cleaner approach: the window W = {j1, j2, ...} means the descent pattern
    // at those positions can differ. The permutation values involved are at
    // 0-indexed positions {j-1 : j in W} union {j : j in W} (capped at len-1).
    let len = pi.len();
    let mut inside_indices: BTreeSet<usize> = BTreeSet::new();
    for &j in window_positions_1indexed {
        if (j as usize) >= 1 {
            inside_indices.insert((j as usize) - 1);
        }
        if (j as usize) < len {
            inside_indices.insert(j as usize);
        }
    }

    let mut outside = Vec::new();
    let mut inside = Vec::new();
    for (i, &v) in pi.iter().enumerate() {
        if inside_indices.contains(&i) {
            inside.push((i as u8, v));
        } else {
            outside.push((i as u8, v));
        }
    }
    (outside, inside)
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

fn format_poly(p: &[i64]) -> String {
    if p.is_empty() || (p.len() == 1 && p[0] == 0) {
        return "0".to_string();
    }
    let parts: Vec<String> = p
        .iter()
        .enumerate()
        .filter(|(_, &c)| c != 0)
        .map(|(k, &c)| {
            if k == 0 {
                format!("{}", c)
            } else if k == 1 {
                if c == 1 {
                    "t".to_string()
                } else {
                    format!("{}t", c)
                }
            } else {
                if c == 1 {
                    format!("t^{}", k)
                } else {
                    format!("{}t^{}", c, k)
                }
            }
        })
        .collect();
    parts.join(" + ")
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    println!("=== Window analysis for varying-source MTP_2 proof ===\n");

    // ========================================================================
    // Part 1: Window size statistics + structural classification
    // ========================================================================
    println!("--- Part 1: Window sizes and structural classification ---\n");

    let mut window_size_counts: BTreeMap<usize, usize> = BTreeMap::new();
    let mut full_window_size_counts: BTreeMap<usize, usize> = BTreeMap::new();

    // Classification of consecutive (p_a, p_b) into cases:
    //   Case A: p_a >= 2 and p_b < n  (both S' and S'' exist for both)
    //   Case B: p_a = 1  (S''_pa doesn't exist, only S'_pa)
    //   Case C: p_b = n  (S''_pb doesn't exist, only S'_pb)
    //   Case D: p_a = 1 and p_b = n (only S' for both, no S'')
    let mut case_counts: BTreeMap<String, usize> = BTreeMap::new();

    for n in 4..=max_n {
        println!("  n = {}", n);
        let max_mask = 1u64 << (n - 1);
        for s_mask in 0..max_mask {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }
            for w in vp.windows(2) {
                let p_a = w[0];
                let p_b = w[1];
                let window = compute_window(s_mask, p_a, p_b, n);
                let full_window = compute_full_window(s_mask, p_a, p_b, n);
                *window_size_counts.entry(window.len()).or_insert(0) += 1;
                *full_window_size_counts.entry(full_window.len()).or_insert(0) += 1;

                let case = if p_a == 1 && p_b == n {
                    "D(pa=1,pb=n)"
                } else if p_a == 1 {
                    "B(pa=1)"
                } else if p_b == n {
                    "C(pb=n)"
                } else {
                    "A(interior)"
                };
                *case_counts.entry(case.to_string()).or_insert(0) += 1;

                let sp_a = source_asc(s_mask, p_a, n);
                let sp_b = source_asc(s_mask, p_b, n);
                let spp_a = source_desc(s_mask, p_a, n);
                let spp_b = source_desc(s_mask, p_b, n);

                if n <= 7 {
                    println!(
                        "    S={} p_a={} p_b={} [{}]: W(S')={} (size {})  W(full)={} (size {})",
                        mask_to_string(s_mask, n),
                        p_a, p_b, case,
                        set_to_string(&window), window.len(),
                        set_to_string(&full_window), full_window.len(),
                    );
                    println!(
                        "      S'_pa={} S'_pb={}",
                        mask_to_string(sp_a, n - 1),
                        mask_to_string(sp_b, n - 1),
                    );
                    if let (Some(a), Some(b)) = (spp_a, spp_b) {
                        println!(
                            "      S''_pa={} S''_pb={}",
                            mask_to_string(a, n - 1),
                            mask_to_string(b, n - 1),
                        );
                    }

                    // New: analyze exactly which 4 source descent sets are involved
                    // and their pairwise relationships
                    let mut all_source_sets: Vec<(String, u64)> = Vec::new();
                    all_source_sets.push(("S'_pa".to_string(), sp_a));
                    if let Some(da) = spp_a {
                        all_source_sets.push(("S''_pa".to_string(), da));
                    }
                    all_source_sets.push(("S'_pb".to_string(), sp_b));
                    if let Some(db) = spp_b {
                        all_source_sets.push(("S''_pb".to_string(), db));
                    }

                    // Check which pairs are equal
                    for i in 0..all_source_sets.len() {
                        for j in (i+1)..all_source_sets.len() {
                            if all_source_sets[i].1 == all_source_sets[j].1 {
                                println!(
                                    "      {} = {} = {}",
                                    all_source_sets[i].0,
                                    all_source_sets[j].0,
                                    mask_to_string(all_source_sets[i].1, n-1),
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    println!("\n  Window size distribution (S' only): {:?}", window_size_counts);
    println!("  Full window size distribution: {:?}", full_window_size_counts);
    println!("  Case distribution: {:?}", case_counts);

    // ========================================================================
    // Part 2: Outside-W configuration analysis
    // ========================================================================
    println!("\n--- Part 2: Outside-W configuration matching ---\n");

    for n in 5..=max_n {
        println!("  ========== n = {} ==========", n);
        let perms_prev = all_permutations(n - 1);
        let mut prev_by_des: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        let max_mask = 1u64 << (n - 1);
        for s_mask in 0..max_mask {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            for w_idx in 0..vp.len() - 1 {
                let p_a = vp[w_idx];
                let p_b = vp[w_idx + 1];

                // Gather sources for p_a and p_b
                let mut src_a: Vec<Vec<u8>> = Vec::new();
                let mut src_b: Vec<Vec<u8>> = Vec::new();

                let sp_a = source_asc(s_mask, p_a, n);
                let spp_a = source_desc(s_mask, p_a, n);
                let sp_b = source_asc(s_mask, p_b, n);
                let spp_b = source_desc(s_mask, p_b, n);

                if let Some(cls) = prev_by_des.get(&sp_a) {
                    src_a.extend(cls.iter().cloned());
                }
                if let Some(sd) = spp_a {
                    if sd != sp_a {
                        if let Some(cls) = prev_by_des.get(&sd) {
                            src_a.extend(cls.iter().cloned());
                        }
                    }
                }

                if let Some(cls) = prev_by_des.get(&sp_b) {
                    src_b.extend(cls.iter().cloned());
                }
                if let Some(sd) = spp_b {
                    if sd != sp_b {
                        if let Some(cls) = prev_by_des.get(&sd) {
                            src_b.extend(cls.iter().cloned());
                        }
                    }
                }

                // Compute the full window
                let full_window = compute_full_window(s_mask, p_a, p_b, n);
                let window_size = full_window.len();

                // Determine inside indices (0-indexed positions in pi involved in window)
                let len = (n - 1) as usize;
                let mut inside_indices: BTreeSet<usize> = BTreeSet::new();
                for &j in &full_window {
                    if (j as usize) >= 1 {
                        inside_indices.insert((j as usize) - 1);
                    }
                    if (j as usize) < len {
                        inside_indices.insert(j as usize);
                    }
                }

                // Outside indices
                let outside_indices: Vec<usize> = (0..len)
                    .filter(|i| !inside_indices.contains(i))
                    .collect();
                let inside_idx_vec: Vec<usize> =
                    inside_indices.iter().copied().collect();

                // Group source permutations by outside configuration
                // outside_config = values at outside positions (in order)
                let extract_outside = |pi: &[u8]| -> Vec<u8> {
                    outside_indices.iter().map(|&i| pi[i]).collect()
                };
                let extract_inside = |pi: &[u8]| -> Vec<u8> {
                    inside_idx_vec.iter().map(|&i| pi[i]).collect()
                };

                // Group A: config -> [(inside_vals, swaps, q)]
                let mut groups_a: BTreeMap<Vec<u8>, Vec<(Vec<u8>, usize, u8)>> =
                    BTreeMap::new();
                for pi in &src_a {
                    let cfg = extract_outside(pi);
                    let s = modified_swaps(pi, p_a);
                    let q = pos_of_nm1(pi);
                    let inside = extract_inside(pi);
                    groups_a
                        .entry(cfg)
                        .or_default()
                        .push((inside, s, q));
                }

                let mut groups_b: BTreeMap<Vec<u8>, Vec<(Vec<u8>, usize, u8)>> =
                    BTreeMap::new();
                for pi in &src_b {
                    let cfg = extract_outside(pi);
                    let s = modified_swaps(pi, p_b);
                    let q = pos_of_nm1(pi);
                    let inside = extract_inside(pi);
                    groups_b
                        .entry(cfg)
                        .or_default()
                        .push((inside, s, q));
                }

                // Analyze relationship
                let all_configs: BTreeSet<Vec<u8>> = groups_a
                    .keys()
                    .chain(groups_b.keys())
                    .cloned()
                    .collect();

                let only_in_a = all_configs
                    .iter()
                    .filter(|c| groups_a.contains_key(*c) && !groups_b.contains_key(*c))
                    .count();
                let only_in_b = all_configs
                    .iter()
                    .filter(|c| !groups_a.contains_key(*c) && groups_b.contains_key(*c))
                    .count();
                let in_both = all_configs
                    .iter()
                    .filter(|c| groups_a.contains_key(*c) && groups_b.contains_key(*c))
                    .count();

                println!(
                    "  S={} p_a={} p_b={} | W={} (size {}, inside_indices={:?})",
                    mask_to_string(s_mask, n),
                    p_a,
                    p_b,
                    set_to_string(&full_window),
                    window_size,
                    inside_idx_vec,
                );
                println!(
                    "    |src_a|={} |src_b|={} | configs: both={} only_a={} only_b={}",
                    src_a.len(),
                    src_b.len(),
                    in_both,
                    only_in_a,
                    only_in_b,
                );

                // For shared configs, show the inside arrangements and (s,q) distributions
                if n <= 6 || window_size <= 2 {
                    let mut mismatch_count = 0;
                    for cfg in all_configs.iter().filter(|c| {
                        groups_a.contains_key(*c) && groups_b.contains_key(*c)
                    }) {
                        let entries_a = &groups_a[cfg];
                        let entries_b = &groups_b[cfg];

                        // Build distributions (s, q) -> count
                        let mut dist_a: BTreeMap<(usize, u8), usize> = BTreeMap::new();
                        let mut dist_b: BTreeMap<(usize, u8), usize> = BTreeMap::new();
                        for (_, s, q) in entries_a {
                            *dist_a.entry((*s, *q)).or_insert(0) += 1;
                        }
                        for (_, s, q) in entries_b {
                            *dist_b.entry((*s, *q)).or_insert(0) += 1;
                        }

                        // Compare
                        let all_keys: BTreeSet<(usize, u8)> = dist_a
                            .keys()
                            .chain(dist_b.keys())
                            .copied()
                            .collect();

                        let same = dist_a == dist_b;

                        if !same || (n <= 6 && in_both <= 10) {
                            println!(
                                "    config={:?} | inside_a entries={} inside_b entries={}",
                                cfg,
                                entries_a.len(),
                                entries_b.len(),
                            );

                            // Show inside arrangements
                            if entries_a.len() <= 8 {
                                for (inside, s, q) in entries_a {
                                    println!(
                                        "      src_a: inside={:?} swaps={} q={}",
                                        inside, s, q
                                    );
                                }
                            }
                            if entries_b.len() <= 8 {
                                for (inside, s, q) in entries_b {
                                    println!(
                                        "      src_b: inside={:?} swaps={} q={}",
                                        inside, s, q
                                    );
                                }
                            }

                            // Show (s,q) distribution comparison
                            if !same {
                                mismatch_count += 1;
                                println!("      ** (s,q) distributions DIFFER:");
                                for key in &all_keys {
                                    let ca = dist_a.get(key).copied().unwrap_or(0);
                                    let cb = dist_b.get(key).copied().unwrap_or(0);
                                    if ca != cb {
                                        println!(
                                            "        (s={}, q={}): a={} b={}",
                                            key.0, key.1, ca, cb
                                        );
                                    }
                                }
                            } else {
                                println!("      (s,q) distributions MATCH");
                            }
                        }
                    }
                    if mismatch_count > 0 {
                        println!(
                            "    ==> {} configs with different (s,q) distributions",
                            mismatch_count
                        );
                    }
                }

                // ============================================================
                // Part 3: For window size 1, check the "toggle descent" bijection idea
                // ============================================================
                if window_size == 1 {
                    let j_star = *full_window.iter().next().unwrap(); // the single window position
                    println!("    >>> Window size 1: j* = {}", j_star);

                    // The inside indices are {j*-1, j*} (0-indexed).
                    // For source(p_a): at position j* we have either a specific descent/ascent constraint.
                    // For source(p_b): the constraint at j* flips.
                    //
                    // Check: can we biject src_a <-> src_b by swapping values at positions j*-1 and j*?
                    let pos0 = (j_star - 1) as usize; // 0-indexed
                    let pos1 = j_star as usize;

                    let mut swap_maps_to_b = 0;
                    let mut swap_maps_to_a = 0;
                    let mut swap_is_identity = 0;
                    let mut swap_total = 0;
                    let src_a_set: BTreeSet<Vec<u8>> = src_a.iter().cloned().collect();
                    let src_b_set: BTreeSet<Vec<u8>> = src_b.iter().cloned().collect();

                    // Also track effect on (swaps, q)
                    let mut swap_delta_s: BTreeMap<i32, usize> = BTreeMap::new();
                    let mut swap_delta_q: BTreeMap<i32, usize> = BTreeMap::new();

                    for pi_a in &src_a {
                        swap_total += 1;
                        let mut pi_swapped = pi_a.clone();
                        if pos1 < pi_swapped.len() {
                            pi_swapped.swap(pos0, pos1);
                        }

                        if src_b_set.contains(&pi_swapped) {
                            swap_maps_to_b += 1;
                            let s_a = modified_swaps(pi_a, p_a) as i32;
                            let s_b = modified_swaps(&pi_swapped, p_b) as i32;
                            let q_a = pos_of_nm1(pi_a) as i32;
                            let q_b = pos_of_nm1(&pi_swapped) as i32;
                            *swap_delta_s.entry(s_b - s_a).or_insert(0) += 1;
                            *swap_delta_q.entry(q_b - q_a).or_insert(0) += 1;
                        } else if src_a_set.contains(&pi_swapped) {
                            swap_is_identity += 1;
                        }
                    }
                    for pi_b in &src_b {
                        let mut pi_swapped = pi_b.clone();
                        if pos1 < pi_swapped.len() {
                            pi_swapped.swap(pos0, pos1);
                        }
                        if src_a_set.contains(&pi_swapped) {
                            swap_maps_to_a += 1;
                        }
                    }

                    println!(
                        "    swap at ({},{}) : {}/{} src_a -> src_b, {}/{} src_b -> src_a, {} self-maps",
                        pos0,
                        pos1,
                        swap_maps_to_b,
                        src_a.len(),
                        swap_maps_to_a,
                        src_b.len(),
                        swap_is_identity,
                    );
                    if !swap_delta_s.is_empty() {
                        println!(
                            "    delta(modified_swaps): {:?}",
                            swap_delta_s
                        );
                        println!(
                            "    delta(q=pos(n-1)): {:?}",
                            swap_delta_q
                        );
                    }
                }

                // ============================================================
                // Part 4: For window size 2, enumerate inside arrangements
                // ============================================================
                if window_size == 2 {
                    let w_vec: Vec<u8> = full_window.iter().copied().collect();
                    println!(
                        "    >>> Window size 2: positions {:?}",
                        w_vec
                    );

                    // The inside indices might be 2 or 3 positions.
                    println!("    inside_indices = {:?} ({} positions)", inside_idx_vec, inside_idx_vec.len());

                    // Enumerate all distinct inside arrangements seen
                    let mut inside_arrangements_a: BTreeSet<Vec<u8>> = BTreeSet::new();
                    let mut inside_arrangements_b: BTreeSet<Vec<u8>> = BTreeSet::new();
                    for pi in &src_a {
                        inside_arrangements_a.insert(extract_inside(pi));
                    }
                    for pi in &src_b {
                        inside_arrangements_b.insert(extract_inside(pi));
                    }

                    // What we care about is: for a fixed set of values at inside positions,
                    // what are the valid arrangements?
                    // Group by the SET of inside values (sorted)
                    let mut inside_by_values_a: BTreeMap<Vec<u8>, BTreeSet<Vec<u8>>> =
                        BTreeMap::new();
                    let mut inside_by_values_b: BTreeMap<Vec<u8>, BTreeSet<Vec<u8>>> =
                        BTreeMap::new();

                    for pi in &src_a {
                        let inside = extract_inside(pi);
                        let mut vals = inside.clone();
                        vals.sort();
                        inside_by_values_a
                            .entry(vals)
                            .or_default()
                            .insert(inside);
                    }
                    for pi in &src_b {
                        let inside = extract_inside(pi);
                        let mut vals = inside.clone();
                        vals.sort();
                        inside_by_values_b
                            .entry(vals)
                            .or_default()
                            .insert(inside);
                    }

                    let all_val_sets: BTreeSet<Vec<u8>> = inside_by_values_a
                        .keys()
                        .chain(inside_by_values_b.keys())
                        .cloned()
                        .collect();

                    if all_val_sets.len() <= 20 {
                        for vals in &all_val_sets {
                            let arrs_a = inside_by_values_a
                                .get(vals)
                                .map(|s| s.iter().cloned().collect::<Vec<_>>())
                                .unwrap_or_default();
                            let arrs_b = inside_by_values_b
                                .get(vals)
                                .map(|s| s.iter().cloned().collect::<Vec<_>>())
                                .unwrap_or_default();
                            println!(
                                "    vals={:?}: arrangements_a={:?} arrangements_b={:?}",
                                vals, arrs_a, arrs_b
                            );
                        }
                    }
                }
            }
        }
    }

    // ========================================================================
    // Part 5: Summary statistics
    // ========================================================================
    println!("\n\n=== SUMMARY ===");
    println!("Window size (S' sym diff) distribution: {:?}", window_size_counts);
    println!("Full window size distribution: {:?}", full_window_size_counts);

    // ========================================================================
    // Part 6: For window size 1 cases, deeper analysis of the bijection
    // ========================================================================
    println!("\n\n--- Part 6: Deeper window-1 bijection analysis ---\n");
    println!("For each window-1 case, classify which descent constraint flips.");
    println!("Check if the swap at j* gives a PERFECT bijection.\n");

    let mut perfect_bij_count = 0u64;
    let mut imperfect_bij_count = 0u64;
    let mut total_w1_cases = 0u64;

    for n in 5..=max_n {
        let perms_prev = all_permutations(n - 1);
        let mut prev_by_des: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        let max_mask = 1u64 << (n - 1);
        for s_mask in 0..max_mask {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            for w_idx in 0..vp.len() - 1 {
                let p_a = vp[w_idx];
                let p_b = vp[w_idx + 1];
                let full_window = compute_full_window(s_mask, p_a, p_b, n);
                if full_window.len() != 1 {
                    continue;
                }
                total_w1_cases += 1;

                let j_star = *full_window.iter().next().unwrap();

                // Gather sources
                let mut src_a: Vec<Vec<u8>> = Vec::new();
                let mut src_b: Vec<Vec<u8>> = Vec::new();
                let sp_a = source_asc(s_mask, p_a, n);
                let spp_a = source_desc(s_mask, p_a, n);
                let sp_b = source_asc(s_mask, p_b, n);
                let spp_b = source_desc(s_mask, p_b, n);

                if let Some(cls) = prev_by_des.get(&sp_a) {
                    src_a.extend(cls.iter().cloned());
                }
                if let Some(sd) = spp_a {
                    if sd != sp_a {
                        if let Some(cls) = prev_by_des.get(&sd) {
                            src_a.extend(cls.iter().cloned());
                        }
                    }
                }
                if let Some(cls) = prev_by_des.get(&sp_b) {
                    src_b.extend(cls.iter().cloned());
                }
                if let Some(sd) = spp_b {
                    if sd != sp_b {
                        if let Some(cls) = prev_by_des.get(&sd) {
                            src_b.extend(cls.iter().cloned());
                        }
                    }
                }

                // Check swap-at-j* bijection
                let pos0 = (j_star - 1) as usize;
                let pos1 = j_star as usize;
                let src_b_set: BTreeSet<Vec<u8>> = src_b.iter().cloned().collect();

                let mut hit_count = 0;
                for pi_a in &src_a {
                    let mut pi_swapped = pi_a.clone();
                    if pos1 < pi_swapped.len() {
                        pi_swapped.swap(pos0, pos1);
                    }
                    if src_b_set.contains(&pi_swapped) {
                        hit_count += 1;
                    }
                }

                let is_perfect = hit_count == src_a.len() && src_a.len() == src_b.len();
                if is_perfect {
                    perfect_bij_count += 1;
                } else {
                    imperfect_bij_count += 1;
                    // Which positions in the descent set differ?
                    let sp_a_set = mask_to_set(sp_a, n - 1);
                    let sp_b_set = mask_to_set(sp_b, n - 1);
                    let in_a_not_b: BTreeSet<u8> =
                        sp_a_set.difference(&sp_b_set).copied().collect();
                    let in_b_not_a: BTreeSet<u8> =
                        sp_b_set.difference(&sp_a_set).copied().collect();

                    println!(
                        "  IMPERFECT: n={} S={} pa={} pb={} j*={} | hit={}/{} |src_b|={}",
                        n,
                        mask_to_string(s_mask, n),
                        p_a,
                        p_b,
                        j_star,
                        hit_count,
                        src_a.len(),
                        src_b.len(),
                    );
                    println!(
                        "    S'_pa \\ S'_pb = {}  S'_pb \\ S'_pa = {}",
                        set_to_string(&in_a_not_b),
                        set_to_string(&in_b_not_a),
                    );
                }
            }
        }
    }

    println!(
        "\nWindow-1 bijection (swap at j*): perfect={} imperfect={} total={}",
        perfect_bij_count, imperfect_bij_count, total_w1_cases,
    );

    // ========================================================================
    // Part 7: For ALL consecutive pairs, what is the precise effect on
    //         swaps and q when we "toggle" inside the window?
    // ========================================================================
    println!("\n\n--- Part 7: Effect of window on (swaps, q) ---\n");
    println!("For each (n, S, p_a, p_b), compute marginal swaps polynomials for source(p_a) and source(p_b),");
    println!("grouped by outside-config, and report whether the marginals are LR-ordered.\n");

    let mut lr_via_config_ok = 0u64;
    let mut lr_via_config_fail = 0u64;

    for n in 5..=max_n {
        let perms_prev = all_permutations(n - 1);
        let mut prev_by_des: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        let max_mask = 1u64 << (n - 1);
        for s_mask in 0..max_mask {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            for w_idx in 0..vp.len() - 1 {
                let p_a = vp[w_idx];
                let p_b = vp[w_idx + 1];

                let mut src_a: Vec<Vec<u8>> = Vec::new();
                let mut src_b: Vec<Vec<u8>> = Vec::new();
                let sp_a = source_asc(s_mask, p_a, n);
                let spp_a = source_desc(s_mask, p_a, n);
                let sp_b = source_asc(s_mask, p_b, n);
                let spp_b = source_desc(s_mask, p_b, n);

                if let Some(cls) = prev_by_des.get(&sp_a) {
                    src_a.extend(cls.iter().cloned());
                }
                if let Some(sd) = spp_a {
                    if sd != sp_a {
                        if let Some(cls) = prev_by_des.get(&sd) {
                            src_a.extend(cls.iter().cloned());
                        }
                    }
                }
                if let Some(cls) = prev_by_des.get(&sp_b) {
                    src_b.extend(cls.iter().cloned());
                }
                if let Some(sd) = spp_b {
                    if sd != sp_b {
                        if let Some(cls) = prev_by_des.get(&sd) {
                            src_b.extend(cls.iter().cloned());
                        }
                    }
                }

                // Compute overall L^{pa} and L^{pb} polynomials
                let vals_a: Vec<usize> = src_a.iter().map(|pi| modified_swaps(pi, p_a)).collect();
                let vals_b: Vec<usize> = src_b.iter().map(|pi| modified_swaps(pi, p_b)).collect();
                let poly_a = build_poly(&vals_a);
                let poly_b = build_poly(&vals_b);

                // Also compute the joint (q, s) matrix for each:
                // M_a[q][s] = #{pi in src_a : pos(n-1)=q, modified_swaps=s}
                let mut joint_a: BTreeMap<u8, BTreeMap<usize, usize>> = BTreeMap::new();
                let mut joint_b: BTreeMap<u8, BTreeMap<usize, usize>> = BTreeMap::new();
                for pi in &src_a {
                    let s = modified_swaps(pi, p_a);
                    let q = pos_of_nm1(pi);
                    *joint_a.entry(q).or_default().entry(s).or_insert(0) += 1;
                }
                for pi in &src_b {
                    let s = modified_swaps(pi, p_b);
                    let q = pos_of_nm1(pi);
                    *joint_b.entry(q).or_default().entry(s).or_insert(0) += 1;
                }

                // Check: for FIXED q, is the swaps polynomial for source(p_a) LR-ordered
                // against the swaps polynomial for source(p_b)?
                let all_qs: BTreeSet<u8> = joint_a
                    .keys()
                    .chain(joint_b.keys())
                    .copied()
                    .collect();

                let mut all_lr_ok = true;
                for &q in &all_qs {
                    let dist_a = joint_a.get(&q);
                    let dist_b = joint_b.get(&q);

                    let poly_qa = match dist_a {
                        Some(d) => {
                            let max_s = *d.keys().max().unwrap_or(&0);
                            let mut p = vec![0i64; max_s + 1];
                            for (&s, &c) in d {
                                p[s] = c as i64;
                            }
                            p
                        }
                        None => vec![0],
                    };
                    let poly_qb = match dist_b {
                        Some(d) => {
                            let max_s = *d.keys().max().unwrap_or(&0);
                            let mut p = vec![0i64; max_s + 1];
                            for (&s, &c) in d {
                                p[s] = c as i64;
                            }
                            p
                        }
                        None => vec![0],
                    };

                    // LR check: for all k, poly_qa[k]*poly_qb[k+1] >= poly_qa[k+1]*poly_qb[k]
                    let max_deg = poly_qa.len().max(poly_qb.len());
                    for k in 0..max_deg {
                        let a_k = if k < poly_qa.len() { poly_qa[k] } else { 0 };
                        let a_k1 = if k + 1 < poly_qa.len() {
                            poly_qa[k + 1]
                        } else {
                            0
                        };
                        let b_k = if k < poly_qb.len() { poly_qb[k] } else { 0 };
                        let b_k1 = if k + 1 < poly_qb.len() {
                            poly_qb[k + 1]
                        } else {
                            0
                        };
                        if a_k * b_k1 < a_k1 * b_k {
                            all_lr_ok = false;
                        }
                    }
                }

                if all_lr_ok {
                    lr_via_config_ok += 1;
                } else {
                    lr_via_config_fail += 1;
                    let full_window = compute_full_window(s_mask, p_a, p_b, n);
                    println!(
                        "  LR-by-q FAIL: n={} S={} pa={} pb={} W={}",
                        n,
                        mask_to_string(s_mask, n),
                        p_a,
                        p_b,
                        set_to_string(&full_window),
                    );
                    // Show the q-refined polynomials
                    for &q in &all_qs {
                        let dist_a = joint_a.get(&q);
                        let dist_b = joint_b.get(&q);
                        let make_poly = |d: Option<&BTreeMap<usize, usize>>| -> Vec<i64> {
                            match d {
                                Some(d) => {
                                    let max_s = *d.keys().max().unwrap_or(&0);
                                    let mut p = vec![0i64; max_s + 1];
                                    for (&s, &c) in d {
                                        p[s] = c as i64;
                                    }
                                    p
                                }
                                None => vec![0],
                            }
                        };
                        let pa_q = make_poly(dist_a);
                        let pb_q = make_poly(dist_b);
                        println!(
                            "    q={}: L^(pa)={} L^(pb)={}",
                            q,
                            format_poly(&pa_q),
                            format_poly(&pb_q),
                        );
                    }
                }
            }
        }
    }

    println!(
        "\nLR-by-q: pass={} fail={} total={}",
        lr_via_config_ok,
        lr_via_config_fail,
        lr_via_config_ok + lr_via_config_fail,
    );

    // ========================================================================
    // Part 8: Four-source decomposition and the S'-only window
    // ========================================================================
    println!("\n\n--- Part 8: Four-source decomposition (S'_pa, S''_pa, S'_pb, S''_pb) ---\n");
    println!("The S' sym diff has size 1 or 2. The S'' = S' union {{p-1}}.");
    println!("For each pair, decompose sources into 4 components and examine relationships.\n");

    // For each consecutive (pa, pb), decompose:
    //   A1 = D(n-1, S'_pa)   with eps1=false (the "non-breaking" sources for pa)
    //   A2 = D(n-1, S''_pa)  with eps1=true  (the "breaking" sources for pa, contribute eps1=1)
    //   B1 = D(n-1, S'_pb)   with eps1=false
    //   B2 = D(n-1, S''_pb)  with eps1=true
    //
    // Question: What is the relationship between A1 and B1?
    // Since S'_pa and S'_pb have sym diff of size 1 or 2, the permutations
    // differ at the descent constraint at those positions.
    //
    // ALSO: is S'_pa = S''_pb or S''_pa = S'_pb ever?
    // If so, a component of source(pa) IS a component of source(pb), giving
    // "shared" permutations (which would help the coupling).

    let mut shared_sp_a_eq_spp_b = 0u64;
    let mut shared_spp_a_eq_sp_b = 0u64;
    let mut total_pairs_p8 = 0u64;

    for n in 5..=max_n {
        let perms_prev = all_permutations(n - 1);
        let mut prev_by_des: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        let max_mask = 1u64 << (n - 1);
        for s_mask in 0..max_mask {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            for w_idx in 0..vp.len() - 1 {
                let p_a = vp[w_idx];
                let p_b = vp[w_idx + 1];
                total_pairs_p8 += 1;

                let sp_a = source_asc(s_mask, p_a, n);
                let spp_a = source_desc(s_mask, p_a, n);
                let sp_b = source_asc(s_mask, p_b, n);
                let spp_b = source_desc(s_mask, p_b, n);

                // Check cross-equalities
                if let Some(db) = spp_b {
                    if sp_a == db {
                        shared_sp_a_eq_spp_b += 1;
                        if n <= 7 {
                            println!(
                                "  S'_pa = S''_pb! n={} S={} pa={} pb={} set={}",
                                n, mask_to_string(s_mask, n), p_a, p_b,
                                mask_to_string(sp_a, n - 1),
                            );
                        }
                    }
                }
                if let (Some(da), _) = (spp_a, sp_b) {
                    if da == sp_b {
                        shared_spp_a_eq_sp_b += 1;
                        if n <= 7 {
                            println!(
                                "  S''_pa = S'_pb! n={} S={} pa={} pb={} set={}",
                                n, mask_to_string(s_mask, n), p_a, p_b,
                                mask_to_string(sp_b, n - 1),
                            );
                        }
                    }
                }

                // Detailed analysis: for S' sym diff of size 1
                let window = compute_window(s_mask, p_a, p_b, n);
                if window.len() == 1 && n <= 7 {
                    let j_star = *window.iter().next().unwrap();

                    // S'_pa and S'_pb differ only at position j_star.
                    // S''_pa = S'_pa ∪ {p_a - 1}
                    // S''_pb = S'_pb ∪ {p_b - 1}
                    //
                    // So S''_pa and S''_pb differ at:
                    //   {j_star} from the S' diff
                    //   PLUS potentially {p_a - 1} and {p_b - 1}
                    //
                    // When p_b = n, S''_pb doesn't exist.
                    // When p_a = 1, S''_pa doesn't exist.

                    // Check: is j_star related to p_a, p_b?
                    println!(
                        "  W(S')={{{}}} n={} S={} pa={} pb={}: j*={} pa-1={} pb-1={} [j*==pa-1? {}] [j*==pb-1? {}]",
                        j_star, n,
                        mask_to_string(s_mask, n),
                        p_a, p_b,
                        j_star, p_a.saturating_sub(1), p_b.saturating_sub(1),
                        j_star == p_a.saturating_sub(1),
                        j_star == p_b.saturating_sub(1),
                    );

                    // If j* = p_b - 1, then:
                    //   S'_pa agrees with S'_pb except at j* = p_b-1
                    //   S'_pa doesn't have p_b-1
                    //   S'_pb has p_b-1
                    //   S''_pa = S'_pa ∪ {p_a - 1}
                    //   S''_pb = S'_pb ∪ {p_b - 1} = S'_pb (since p_b-1 already in S'_pb)
                    //   So S''_pb = S'_pb! This means source for p_b has S'_pb = S''_pb.

                    if j_star == p_b - 1 && spp_b.is_some() {
                        let is_same = sp_b == spp_b.unwrap();
                        println!("    => j*=pb-1: S'_pb == S''_pb? {}", is_same);
                    }

                    // Compute the four source swap polynomials
                    let get_swaps_poly = |des_mask: u64, p: u8| -> Vec<i64> {
                        let mut vals = Vec::new();
                        if let Some(perms) = prev_by_des.get(&des_mask) {
                            for pi in perms {
                                vals.push(modified_swaps(pi, p));
                            }
                        }
                        build_poly(&vals)
                    };

                    let poly_sp_a = get_swaps_poly(sp_a, p_a);
                    let poly_sp_b = get_swaps_poly(sp_b, p_b);
                    println!(
                        "    poly(S'_pa, pa)={} poly(S'_pb, pb)={}",
                        format_poly(&poly_sp_a), format_poly(&poly_sp_b),
                    );
                    if let Some(da) = spp_a {
                        let poly_spp_a = get_swaps_poly(da, p_a);
                        println!("    poly(S''_pa, pa)={}", format_poly(&poly_spp_a));
                    }
                    if let Some(db) = spp_b {
                        let poly_spp_b = get_swaps_poly(db, p_b);
                        println!("    poly(S''_pb, pb)={}", format_poly(&poly_spp_b));
                    }
                }
            }
        }
    }

    println!(
        "\nCross-equalities: S'_pa=S''_pb: {} times, S''_pa=S'_pb: {} times (out of {})",
        shared_sp_a_eq_spp_b, shared_spp_a_eq_sp_b, total_pairs_p8,
    );

    // ========================================================================
    // Part 9: When S' window = {j*} with j* = p_b - 1:
    //         Analyze the "toggle descent" bijection on the S'-only components
    // ========================================================================
    println!("\n\n--- Part 9: Toggle-descent bijection for S'-window of size 1 ---\n");
    println!("When S'_pa and S'_pb differ at exactly one position j*,");
    println!("can we biject D(n-1, S'_pa) <-> D(n-1, S'_pb) by swapping values at j*-1, j* ?");
    println!("And what is the effect on modified_swaps(pi, p_a) vs modified_swaps(swapped, p_b)?\n");

    let mut bij_results: BTreeMap<(String, String), usize> = BTreeMap::new();

    for n in 5..=max_n {
        let perms_prev = all_permutations(n - 1);
        let mut prev_by_des: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        let max_mask = 1u64 << (n - 1);
        for s_mask in 0..max_mask {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            for w_idx in 0..vp.len() - 1 {
                let p_a = vp[w_idx];
                let p_b = vp[w_idx + 1];
                let window = compute_window(s_mask, p_a, p_b, n);
                if window.len() != 1 {
                    continue;
                }

                let j_star = *window.iter().next().unwrap();
                let sp_a = source_asc(s_mask, p_a, n);
                let sp_b = source_asc(s_mask, p_b, n);

                // Check: at j*, which side has the descent?
                let a_has_j = (sp_a >> (j_star - 1)) & 1 == 1;
                let b_has_j = (sp_b >> (j_star - 1)) & 1 == 1;

                let src_a: Vec<Vec<u8>> = prev_by_des.get(&sp_a)
                    .map(|v| v.clone())
                    .unwrap_or_default();
                let src_b: Vec<Vec<u8>> = prev_by_des.get(&sp_b)
                    .map(|v| v.clone())
                    .unwrap_or_default();

                let src_b_set: BTreeSet<Vec<u8>> = src_b.iter().cloned().collect();

                // Try swapping at (j*-1, j*) in 0-indexed
                let pos0 = (j_star - 1) as usize;
                let pos1 = j_star as usize;

                let mut swap_success = 0;
                let mut delta_swaps_counts: BTreeMap<i32, usize> = BTreeMap::new();
                let mut delta_q_counts: BTreeMap<i32, usize> = BTreeMap::new();

                for pi_a in &src_a {
                    let mut pi_swapped = pi_a.clone();
                    if pos1 < pi_swapped.len() {
                        pi_swapped.swap(pos0, pos1);
                    }
                    if src_b_set.contains(&pi_swapped) {
                        swap_success += 1;
                        let s_before = modified_swaps(pi_a, p_a) as i32;
                        let s_after = modified_swaps(&pi_swapped, p_b) as i32;
                        let q_before = pos_of_nm1(pi_a) as i32;
                        let q_after = pos_of_nm1(&pi_swapped) as i32;
                        *delta_swaps_counts.entry(s_after - s_before).or_insert(0) += 1;
                        *delta_q_counts.entry(q_after - q_before).or_insert(0) += 1;
                    }
                }

                let bij_type = if swap_success == src_a.len() && src_a.len() == src_b.len() {
                    "PERFECT"
                } else {
                    "PARTIAL"
                };

                *bij_results
                    .entry((
                        bij_type.to_string(),
                        format!("a_has_j*={}, b_has_j*={}", a_has_j, b_has_j),
                    ))
                    .or_insert(0) += 1;

                println!(
                    "  n={} S={} pa={} pb={} j*={}: {}, hit {}/{}, |src_b|={}, a_has_j*={} b_has_j*={}",
                    n, mask_to_string(s_mask, n), p_a, p_b, j_star,
                    bij_type, swap_success, src_a.len(), src_b.len(),
                    a_has_j, b_has_j,
                );
                if !delta_swaps_counts.is_empty() {
                    println!(
                        "    delta(swaps): {:?}  delta(q): {:?}",
                        delta_swaps_counts, delta_q_counts,
                    );
                }
            }
        }
    }

    println!("\nBijection results summary: {:?}", bij_results);

    // ========================================================================
    // Part 10: S' window of size 2 -- what are the two differing positions?
    // ========================================================================
    println!("\n\n--- Part 10: S' window of size 2 classification ---\n");

    let mut w2_pattern_counts: BTreeMap<String, usize> = BTreeMap::new();

    for n in 5..=max_n {
        let max_mask = 1u64 << (n - 1);
        for s_mask in 0..max_mask {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }
            for w_idx in 0..vp.len() - 1 {
                let p_a = vp[w_idx];
                let p_b = vp[w_idx + 1];
                let window = compute_window(s_mask, p_a, p_b, n);
                if window.len() != 2 {
                    continue;
                }

                let w_vec: Vec<u8> = window.iter().copied().collect();
                let j1 = w_vec[0];
                let j2 = w_vec[1];

                let sp_a = source_asc(s_mask, p_a, n);
                let sp_b = source_asc(s_mask, p_b, n);

                // At j1: who has descent?
                let a_j1 = (sp_a >> (j1 - 1)) & 1 == 1;
                let b_j1 = (sp_b >> (j1 - 1)) & 1 == 1;
                let a_j2 = (sp_a >> (j2 - 1)) & 1 == 1;
                let b_j2 = (sp_b >> (j2 - 1)) & 1 == 1;

                // Classify the pattern
                // The sym diff at j1: a_j1 vs b_j1 must differ, same for j2
                let pattern = format!(
                    "j1: a={} b={}; j2: a={} b={}; gap=j2-j1={}; j1-pa+1={}; j2-pb+1={}",
                    if a_j1 { "D" } else { "A" },
                    if b_j1 { "D" } else { "A" },
                    if a_j2 { "D" } else { "A" },
                    if b_j2 { "D" } else { "A" },
                    j2 - j1,
                    j1 as i32 - p_a as i32 + 1,
                    j2 as i32 - p_b as i32 + 1,
                );

                *w2_pattern_counts.entry(pattern.clone()).or_insert(0) += 1;

                if n <= 7 {
                    println!(
                        "  n={} S={} pa={} pb={}: W(S')={{{},{}}} {}",
                        n, mask_to_string(s_mask, n), p_a, p_b, j1, j2, pattern,
                    );
                }
            }
        }
    }

    println!("\nSize-2 window pattern counts: {:?}", w2_pattern_counts);

    // ========================================================================
    // Part 11: Size-2 window bijection attempts
    // ========================================================================
    println!("\n\n--- Part 11: Size-2 S' window -- bijection via descent exchange ---\n");
    println!("S'_pa has A at j1, D at j2. S'_pb has D at j1, A at j2.");
    println!("Try: swap values at j1 and j2 (0-indexed: j1-1 and j2-1).");
    println!("Also try: swap at j1-1,j1 and then j2-1,j2 (adjacent transpositions).\n");

    let mut w2_bij_results: BTreeMap<String, usize> = BTreeMap::new();

    for n in 5..=max_n {
        let perms_prev = all_permutations(n - 1);
        let mut prev_by_des: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        let max_mask = 1u64 << (n - 1);
        for s_mask in 0..max_mask {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            for w_idx in 0..vp.len() - 1 {
                let p_a = vp[w_idx];
                let p_b = vp[w_idx + 1];
                let window = compute_window(s_mask, p_a, p_b, n);
                if window.len() != 2 {
                    continue;
                }

                let w_vec: Vec<u8> = window.iter().copied().collect();
                let j1 = w_vec[0] as usize; // 1-indexed position
                let j2 = w_vec[1] as usize;

                let sp_a = source_asc(s_mask, p_a, n);
                let sp_b = source_asc(s_mask, p_b, n);
                let spp_a = source_desc(s_mask, p_a, n);
                let spp_b = source_desc(s_mask, p_b, n);

                // Work with S'-only sources first
                let src_a_prime: Vec<Vec<u8>> = prev_by_des.get(&sp_a)
                    .map(|v| v.clone())
                    .unwrap_or_default();
                let src_b_prime: Vec<Vec<u8>> = prev_by_des.get(&sp_b)
                    .map(|v| v.clone())
                    .unwrap_or_default();

                let src_b_prime_set: BTreeSet<Vec<u8>> = src_b_prime.iter().cloned().collect();
                let src_a_prime_set: BTreeSet<Vec<u8>> = src_a_prime.iter().cloned().collect();

                // Strategy 1: swap values at 0-indexed positions j1-1 and j2-1
                let mut s1_hits = 0;
                let mut s1_delta_swaps: BTreeMap<i32, usize> = BTreeMap::new();
                let mut s1_delta_q: BTreeMap<i32, usize> = BTreeMap::new();
                for pi in &src_a_prime {
                    let mut swapped = pi.clone();
                    swapped.swap(j1 - 1, j2 - 1);
                    if src_b_prime_set.contains(&swapped) {
                        s1_hits += 1;
                        let ds = modified_swaps(&swapped, p_b) as i32 - modified_swaps(pi, p_a) as i32;
                        let dq = pos_of_nm1(&swapped) as i32 - pos_of_nm1(pi) as i32;
                        *s1_delta_swaps.entry(ds).or_insert(0) += 1;
                        *s1_delta_q.entry(dq).or_insert(0) += 1;
                    }
                }

                // Strategy 2: swap adjacent at j2-1,j2 (to fix the descent at j2)
                let mut s2_hits = 0;
                for pi in &src_a_prime {
                    let mut swapped = pi.clone();
                    if j2 < swapped.len() {
                        swapped.swap(j2 - 1, j2);
                    }
                    if src_b_prime_set.contains(&swapped) {
                        s2_hits += 1;
                    }
                }

                // Strategy 3: swap adjacent at j1-1,j1 (to fix the descent at j1)
                let mut s3_hits = 0;
                for pi in &src_a_prime {
                    let mut swapped = pi.clone();
                    swapped.swap(j1 - 1, j1);
                    if src_b_prime_set.contains(&swapped) {
                        s3_hits += 1;
                    }
                }

                // Strategy 4: both adjacent swaps: j1-1<->j1 and j2-1<->j2
                let mut s4_hits = 0;
                for pi in &src_a_prime {
                    let mut swapped = pi.clone();
                    swapped.swap(j1 - 1, j1);
                    if j2 < swapped.len() {
                        swapped.swap(j2 - 1, j2);
                    }
                    if src_b_prime_set.contains(&swapped) {
                        s4_hits += 1;
                    }
                }

                // Strategy 5: swap values at 0-indexed j1 and j2-1 (the inner pair)
                let mut s5_hits = 0;
                let mut s5_delta_swaps: BTreeMap<i32, usize> = BTreeMap::new();
                let mut s5_delta_q: BTreeMap<i32, usize> = BTreeMap::new();
                if j1 < j2 - 1 {
                    for pi in &src_a_prime {
                        let mut swapped = pi.clone();
                        swapped.swap(j1, j2 - 1);
                        if src_b_prime_set.contains(&swapped) {
                            s5_hits += 1;
                            let ds = modified_swaps(&swapped, p_b) as i32 - modified_swaps(pi, p_a) as i32;
                            let dq = pos_of_nm1(&swapped) as i32 - pos_of_nm1(pi) as i32;
                            *s5_delta_swaps.entry(ds).or_insert(0) += 1;
                            *s5_delta_q.entry(dq).or_insert(0) += 1;
                        }
                    }
                }

                let is_s1_perfect = s1_hits == src_a_prime.len() && src_a_prime.len() == src_b_prime.len();
                let best = *[s1_hits, s2_hits, s3_hits, s4_hits, s5_hits].iter().max().unwrap();
                let best_name = if best == s1_hits { "S1(j1-1<->j2-1)" }
                    else if best == s2_hits { "S2(adj-j2)" }
                    else if best == s3_hits { "S3(adj-j1)" }
                    else if best == s4_hits { "S4(both-adj)" }
                    else { "S5(j1<->j2-1)" };

                *w2_bij_results
                    .entry(if is_s1_perfect { "S1_PERFECT" } else { "S1_partial" }.to_string())
                    .or_insert(0) += 1;

                println!(
                    "  n={} S={} pa={} pb={} j1={} j2={}: |S'_a|={} |S'_b|={} S1={}/{} S2={} S3={} S4={} S5={} best={}",
                    n, mask_to_string(s_mask, n), p_a, p_b, j1, j2,
                    src_a_prime.len(), src_b_prime.len(),
                    s1_hits, src_a_prime.len(),
                    s2_hits, s3_hits, s4_hits, s5_hits,
                    best_name,
                );
                if s1_hits > 0 {
                    println!(
                        "    S1 delta(swaps): {:?}  delta(q): {:?}",
                        s1_delta_swaps, s1_delta_q,
                    );
                }
                if s5_hits > 0 {
                    println!(
                        "    S5 delta(swaps): {:?}  delta(q): {:?}",
                        s5_delta_swaps, s5_delta_q,
                    );
                }

                // Also check: for the COMBINED source (S' union S''), what are the sizes?
                let mut full_src_a: Vec<Vec<u8>> = src_a_prime.clone();
                let mut full_src_b: Vec<Vec<u8>> = src_b_prime.clone();
                if let Some(da) = spp_a {
                    if da != sp_a {
                        if let Some(cls) = prev_by_des.get(&da) {
                            full_src_a.extend(cls.iter().cloned());
                        }
                    }
                }
                if let Some(db) = spp_b {
                    if db != sp_b {
                        if let Some(cls) = prev_by_des.get(&db) {
                            full_src_b.extend(cls.iter().cloned());
                        }
                    }
                }
                println!(
                    "    full: |src_a|={} |src_b|={} (S''_a adds {}, S''_b adds {})",
                    full_src_a.len(), full_src_b.len(),
                    full_src_a.len() - src_a_prime.len(),
                    full_src_b.len() - src_b_prime.len(),
                );
            }
        }
    }

    println!("\nSize-2 bijection results: {:?}", w2_bij_results);

    // ========================================================================
    // Part 12: Check the PRECISE relationship between j1, j2, p_a, p_b
    // ========================================================================
    println!("\n\n--- Part 12: Precise relationship between j1, j2, p_a, p_b ---\n");
    println!("In size-2 S' window: j2 = p_b - 1 always.");
    println!("What about j1? Is j1 = p_a - 1?  Or j1 = p_a?\n");

    for n in 5..=max_n {
        let max_mask = 1u64 << (n - 1);
        for s_mask in 0..max_mask {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }
            for w_idx in 0..vp.len() - 1 {
                let p_a = vp[w_idx];
                let p_b = vp[w_idx + 1];
                let window = compute_window(s_mask, p_a, p_b, n);
                if window.len() != 2 {
                    continue;
                }

                let w_vec: Vec<u8> = window.iter().copied().collect();
                let j1 = w_vec[0];
                let j2 = w_vec[1];

                // Compute p_a - 1 in source coords
                // After removing position p_a from [n-1], positions p_a+1,...,n-1 shift down
                // S'_pa doesn't include p_a-1 as a descent, it shifts things above p_a down
                // j2 should be p_b - 1 (the shifted version of p_b in the (n-1) world)
                // j1 should be p_a - 1 (the position just before p_a)

                // But wait -- for p_a >= 2, position p_a-1 in the source is
                // the same as position p_a-1 in the original (no shift needed, as p is removed)
                // And for p_b > p_a, position p_b - 1 in the source is the shifted version.

                let expected_j1 = p_a - 1; // but only if p_a >= 2... if p_a = 1, unclear
                let expected_j2 = p_b - 1;

                if j1 != expected_j1 || j2 != expected_j2 {
                    println!(
                        "  UNEXPECTED: n={} S={} pa={} pb={}: j1={} (expected {}) j2={} (expected {})",
                        n, mask_to_string(s_mask, n), p_a, p_b, j1, expected_j1, j2, expected_j2,
                    );
                }
            }
        }
    }

    println!("  (If no UNEXPECTED lines above, then j1 = p_a - 1 and j2 = p_b - 1 always.)");

    // ========================================================================
    // Part 13: For size-1 S' window (pb=n): deeper analysis
    // ========================================================================
    println!("\n\n--- Part 13: Size-1 S' window (pb=n) ---\n");
    println!("j* is the single differing position. What is j* in terms of p_a?\n");

    for n in 5..=max_n {
        let max_mask = 1u64 << (n - 1);
        for s_mask in 0..max_mask {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }
            for w_idx in 0..vp.len() - 1 {
                let p_a = vp[w_idx];
                let p_b = vp[w_idx + 1];
                if p_b != n {
                    continue;
                }
                let window = compute_window(s_mask, p_a, p_b, n);
                assert!(window.len() == 1, "expected size 1 for pb=n");

                let j_star = *window.iter().next().unwrap();

                // When p_b = n, S'_pb = S (the original target descent set)
                // S'_pa = S ∩ [p_a-2] ∪ {j-1 : j ∈ S, j > p_a}
                // The difference is at positions where the "shift" matters.
                // Since p_b = n, the last element of S might be n-1 or less.
                // j* should be related to the highest element of S or something near p_a.

                // After the shift, S'_pa maps S elements above p_a to {j-1 : j > p_a, j ∈ S}.
                // S'_pb = S'_n = S (no shift).
                // The difference is at positions where either:
                //   (1) some j ∈ S with j > p_a maps to j-1 in S'_pa but j in S'_pb, or
                //   (2) p_a itself is in S'_pb but not in S'_pa (since S'_pa skips position p_a-1).

                // Actually S'_pa for p_a < n:
                //   positions 1..p_a-2 of S'_pa = same as S
                //   positions p_a-1..n-2 of S'_pa = positions p_a..n-1 of S, shifted down by 1
                // S'_n = S: positions 1..n-1

                // At position j < p_a - 1: S'_pa[j] = S[j], S'_n[j] = S[j] → same
                // At position j = p_a - 1: S'_pa[p_a-1] = S[p_a] (shifted), S'_n[p_a-1] = S[p_a-1]
                //   These differ iff S[p_a] != S[p_a-1]
                // At position j >= p_a: S'_pa[j] = S[j+1] (shifted), S'_n[j] = S[j]
                //   These differ iff S[j+1] != S[j]

                // So the sym diff is at positions where the shift changes the descent pattern.
                // For the sym diff to have size 1, there must be exactly one position where
                // S[j+1] != S[j] (for j >= p_a-1 in the source indexing).

                // If S has no elements in [p_a-1, n-1], then S'_pa = S'_n (trivially), so size 0.
                // Actually, S always has an element (at least p_a which starts a run).

                // In the valid position structure: p_a starts a run (p_a ∈ S and p_a-1 ∉ S),
                // and the next valid position is n (meaning S has elements p_a, p_a+1, ..., p_a+k-1
                // forming a consecutive block, with p_a+k-1 being the last element of S).

                // The shift only affects position p_a-1:
                //   S'_pa[p_a-1] = S[p_a] = 1 (since p_a ∈ S)
                //   S'_n[p_a-1] = S[p_a-1] = 0 (since p_a-1 ∉ S, as p_a starts a run)
                //   These DIFFER.
                // For j >= p_a: S'_pa[j] = S[j+1], S'_n[j] = S[j]
                //   Within the run: S[j+1] = 1 = S[j] if both in the run → same
                //   At the end of the run: S[last] = 1 but S[last+1] = 0, while S[last-1] = 1.
                //     S'_pa[last-1] = S[last] = 1, S'_n[last-1] = S[last-1] = 1 → same
                //     S'_pa[last] = S[last+1] = 0, S'_n[last] = S[last] = 1 → DIFFER!
                //   But wait, last = p_a + k - 1, and the next position after the run...
                //   Actually if the run goes p_a, p_a+1, ..., m (last element of S above p_a),
                //   then in the source:
                //     position m-1 (shifted from m): S'_pa[m-1] = S[m] = 1, S'_n[m-1] = S[m-1] = 1 → same
                //     position m (shifted from m+1): S'_pa[m] = S[m+1] = 0, S'_n[m] = S[m] = 1 → DIFFER!

                // So we get TWO differing positions: p_a-1 and m (last S-element above p_a).
                // But the experiment says size 1! Let me re-examine...

                // Hmm wait. If S = {p_a, p_a+1, ..., n-1} (the run goes all the way to n-1),
                // then m = n-1, and position m = n-1 is BEYOND the index range [1..n-2],
                // so only position p_a-1 differs. That gives size 1!
                // And if pb = n, the run from p_a must go all the way to n-1.
                // Because: the valid positions are starts of maximal consecutive blocks in S.
                // If pb = n (the rightmost valid position), it means n-1 ∉ S.
                // But wait, pb=n means there IS no S element at n-1...
                // Actually pb=n is always a valid position when (n-1) ∉ S (the last position is always valid if there's no descent at n-1).

                // Let me just record j* vs p_a.
                println!(
                    "  n={} S={} pa={} pb=n: j*={} (j*==pa-1? {}, j*==pa? {})",
                    n, mask_to_string(s_mask, n), p_a,
                    j_star, j_star == p_a - 1, j_star == p_a as u8,
                );
            }
        }
    }
}
