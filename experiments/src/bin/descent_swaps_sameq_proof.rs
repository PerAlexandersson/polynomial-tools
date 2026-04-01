/// Same-q cross-chain compatibility proof exploration.
///
/// For consecutive valid positions p_a < p_b and each q, examine the modified source
/// polynomials F~_q^{(p_a)} and F~_q^{(p_b)} in detail.
///
/// For n=5 alternating S={2,4}: p_a=2, p_b=4:
///   Source for p_a=2: D(4, S'_2) ∪ D(4, S''_2) = D(4, {3}) ∪ D(4, {1,3})
///   Source for p_b=4: D(4, S'_4) ∪ D(4, S''_4) = D(4, {2}) ∪ D(4, {2,3})
///
/// Lists ALL permutations, shows window analysis, and tests for bijection structure.
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{format_poly, is_real_rooted};
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

fn mask_to_string(mask: u64, n: u8) -> String {
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

fn perm_to_string(pi: &[u8]) -> String {
    let parts: Vec<String> = pi.iter().map(|x| x.to_string()).collect();
    format!("[{}]", parts.join(","))
}

/// Check if f and g are pairwise compatible (all nonnegative combinations are real-rooted).
fn check_pairwise_compatible(f: &[i64], g: &[i64]) -> bool {
    if f.iter().all(|&c| c == 0) || g.iter().all(|&c| c == 0) {
        return true;
    }
    if f.len() <= 1 && g.len() <= 1 {
        return true;
    }
    if !is_real_rooted(f) || !is_real_rooted(g) {
        return false;
    }
    let weights: Vec<(i64, i64)> = vec![
        (1, 1),
        (1, 2),
        (2, 1),
        (1, 3),
        (3, 1),
        (1, 5),
        (5, 1),
        (2, 3),
        (3, 2),
        (1, 7),
        (7, 1),
        (3, 5),
        (5, 3),
    ];
    let maxlen = f.len().max(g.len());
    for (a, b) in &weights {
        let mut combo = vec![0i64; maxlen];
        for (i, &c) in f.iter().enumerate() {
            combo[i] += a * c;
        }
        for (i, &c) in g.iter().enumerate() {
            combo[i] += b * c;
        }
        while combo.len() > 1 && *combo.last().unwrap() == 0 {
            combo.pop();
        }
        if combo.len() > 1 && !is_real_rooted(&combo) {
            return false;
        }
    }
    true
}

/// Information about a single permutation's contribution.
#[derive(Clone, Debug)]
struct PermInfo {
    pi: Vec<u8>,
    descent_mask: u64,
    source_type: &'static str, // "asc" or "desc"
    q: u8,                     // position of value n-1
    swaps: usize,
    eps1: bool,
    modified_swaps: usize,
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    println!("=== Same-q cross-chain compatibility: detailed proof exploration ===\n");

    for n in 5..=max_n {
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

        // Only look at alternating S for detailed analysis
        // S = {2,4,...} (alternating descent set)
        let mut s_masks: Vec<u64> = Vec::new();
        // Alternating S: descents at even positions
        let alt_mask = {
            let mut m = 0u64;
            for i in (2..n).step_by(2) {
                m |= 1 << (i - 1);
            }
            m
        };
        s_masks.push(alt_mask);

        // For larger n, also include non-alternating sets for completeness testing
        if n <= 6 {
            for s_mask in 0u64..(1 << (n - 1)) {
                if s_mask & 1 != 0 {
                    continue;
                }
                if s_mask != alt_mask {
                    s_masks.push(s_mask);
                }
            }
        }

        for &s_mask in &s_masks {
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            let is_alt = s_mask == alt_mask;
            if !is_alt && n > 5 {
                continue; // Skip non-alternating for n > 5 except in summary
            }

            println!("{}", "=".repeat(80));
            println!(
                "n={}, S={}, P(S)={:?}{}",
                n,
                mask_to_string(s_mask, n),
                vp,
                if is_alt { " [ALTERNATING]" } else { "" }
            );
            println!("{}", "=".repeat(80));

            // For each consecutive pair of valid positions
            for idx in 0..vp.len() - 1 {
                let p_a = vp[idx];
                let p_b = vp[idx + 1];

                let sp_a_asc = source_asc(s_mask, p_a, n);
                let sp_a_desc = source_desc(s_mask, p_a, n);
                let sp_b_asc = source_asc(s_mask, p_b, n);
                let sp_b_desc = source_desc(s_mask, p_b, n);

                println!("\n--- p_a={}, p_b={} ---", p_a, p_b);
                println!(
                    "  Source for p_a={}: S'_{{p_a}} = {}  S''_{{p_a}} = {}",
                    p_a,
                    mask_to_string(sp_a_asc, n - 1),
                    sp_a_desc
                        .map(|m| mask_to_string(m, n - 1))
                        .unwrap_or_else(|| "N/A".to_string())
                );
                println!(
                    "  Source for p_b={}: S'_{{p_b}} = {}  S''_{{p_b}} = {}",
                    p_b,
                    mask_to_string(sp_b_asc, n - 1),
                    sp_b_desc
                        .map(|m| mask_to_string(m, n - 1))
                        .unwrap_or_else(|| "N/A".to_string())
                );
                println!(
                    "  Window: [{}, {}]",
                    p_a.saturating_sub(1),
                    p_b.saturating_sub(1)
                );

                // Collect all permutations for p_a
                let mut perms_a: Vec<PermInfo> = Vec::new();
                if let Some(indices) = by_descent_prev.get(&sp_a_asc) {
                    for &i in indices {
                        let pi = perm_prev_list[i];
                        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                        let sw = compute(pi, Stat::Swaps);
                        let e1 = epsilon1(pi, p_a);
                        perms_a.push(PermInfo {
                            pi: pi.clone(),
                            descent_mask: sp_a_asc,
                            source_type: "asc",
                            q,
                            swaps: sw,
                            eps1: e1,
                            modified_swaps: sw + if e1 { 1 } else { 0 },
                        });
                    }
                }
                if let Some(sp) = sp_a_desc {
                    if let Some(indices) = by_descent_prev.get(&sp) {
                        for &i in indices {
                            let pi = perm_prev_list[i];
                            let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            let sw = compute(pi, Stat::Swaps);
                            // eps1 is 0 for descent source (always)
                            perms_a.push(PermInfo {
                                pi: pi.clone(),
                                descent_mask: sp,
                                source_type: "desc",
                                q,
                                swaps: sw,
                                eps1: false,
                                modified_swaps: sw,
                            });
                        }
                    }
                }

                // Collect all permutations for p_b
                let mut perms_b: Vec<PermInfo> = Vec::new();
                if let Some(indices) = by_descent_prev.get(&sp_b_asc) {
                    for &i in indices {
                        let pi = perm_prev_list[i];
                        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                        let sw = compute(pi, Stat::Swaps);
                        let e1 = epsilon1(pi, p_b);
                        perms_b.push(PermInfo {
                            pi: pi.clone(),
                            descent_mask: sp_b_asc,
                            source_type: "asc",
                            q,
                            swaps: sw,
                            eps1: e1,
                            modified_swaps: sw + if e1 { 1 } else { 0 },
                        });
                    }
                }
                if let Some(sp) = sp_b_desc {
                    if let Some(indices) = by_descent_prev.get(&sp) {
                        for &i in indices {
                            let pi = perm_prev_list[i];
                            let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            let sw = compute(pi, Stat::Swaps);
                            perms_b.push(PermInfo {
                                pi: pi.clone(),
                                descent_mask: sp,
                                source_type: "desc",
                                q,
                                swaps: sw,
                                eps1: false,
                                modified_swaps: sw,
                            });
                        }
                    }
                }

                // Group by q
                let mut by_q_a: BTreeMap<u8, Vec<&PermInfo>> = BTreeMap::new();
                let mut by_q_b: BTreeMap<u8, Vec<&PermInfo>> = BTreeMap::new();
                for p in &perms_a {
                    by_q_a.entry(p.q).or_default().push(p);
                }
                for p in &perms_b {
                    by_q_b.entry(p.q).or_default().push(p);
                }

                // All q values present in either
                let mut all_qs: Vec<u8> = by_q_a.keys().chain(by_q_b.keys()).copied().collect();
                all_qs.sort();
                all_qs.dedup();

                // Only show full permutations for small n
                let show_perms = n <= 6;

                for &q in &all_qs {
                    let group_a = by_q_a.get(&q).cloned().unwrap_or_default();
                    let group_b = by_q_b.get(&q).cloned().unwrap_or_default();

                    if group_a.is_empty() && group_b.is_empty() {
                        continue;
                    }

                    // Build polynomials
                    let vals_a: Vec<usize> = group_a.iter().map(|p| p.modified_swaps).collect();
                    let vals_b: Vec<usize> = group_b.iter().map(|p| p.modified_swaps).collect();
                    let poly_a = build_poly(&vals_a);
                    let poly_b = build_poly(&vals_b);

                    let compatible = check_pairwise_compatible(&poly_a, &poly_b);

                    println!(
                        "\n  q={}: F~_q^{{(p_a={})}} = {}  |  F~_q^{{(p_b={})}} = {}  | compatible: {}",
                        q,
                        p_a,
                        format_poly(&poly_a),
                        p_b,
                        format_poly(&poly_b),
                        if compatible { "YES" } else { "**NO**" }
                    );

                    if show_perms {
                        println!("    Permutations for p_a={}:", p_a);
                        for info in &group_a {
                            println!(
                                "      pi={}  Des={}  src={}  swaps={}  eps1={}  mod_swaps={}",
                                perm_to_string(&info.pi),
                                mask_to_string(info.descent_mask, n - 1),
                                info.source_type,
                                info.swaps,
                                info.eps1 as u8,
                                info.modified_swaps,
                            );
                        }
                        println!("    Permutations for p_b={}:", p_b);
                        for info in &group_b {
                            println!(
                                "      pi={}  Des={}  src={}  swaps={}  eps1={}  mod_swaps={}",
                                perm_to_string(&info.pi),
                                mask_to_string(info.descent_mask, n - 1),
                                info.source_type,
                                info.swaps,
                                info.eps1 as u8,
                                info.modified_swaps,
                            );
                        }

                        // Window analysis: compare permutations outside the window
                        let window_lo = p_a.saturating_sub(1) as usize; // 0-indexed position p_a-1
                        let window_hi = (p_b.saturating_sub(1)) as usize; // 0-indexed position p_b-1
                                                                          // Positions OUTSIDE the window (0-indexed)
                        println!(
                            "    Window analysis (0-indexed positions [{}, {}]):",
                            window_lo, window_hi
                        );

                        // For each pair (a,b), check if they agree outside the window
                        let mut bijection_candidates: Vec<(usize, usize, i64)> = Vec::new();
                        for (ia, info_a) in group_a.iter().enumerate() {
                            for (ib, info_b) in group_b.iter().enumerate() {
                                let mut agree_outside = true;
                                for pos in 0..info_a.pi.len() {
                                    if pos >= window_lo && pos <= window_hi {
                                        continue;
                                    }
                                    if info_a.pi[pos] != info_b.pi[pos] {
                                        agree_outside = false;
                                        break;
                                    }
                                }
                                let delta =
                                    info_b.modified_swaps as i64 - info_a.modified_swaps as i64;
                                if agree_outside {
                                    bijection_candidates.push((ia, ib, delta));
                                    println!(
                                        "      MATCH: {} <-> {}  delta_mod_swaps={}  window_a={:?} window_b={:?}",
                                        perm_to_string(&info_a.pi),
                                        perm_to_string(&info_b.pi),
                                        delta,
                                        &info_a.pi[window_lo..=window_hi.min(info_a.pi.len()-1)],
                                        &info_b.pi[window_lo..=window_hi.min(info_b.pi.len()-1)],
                                    );
                                }
                            }
                        }

                        if bijection_candidates.is_empty()
                            && !group_a.is_empty()
                            && !group_b.is_empty()
                        {
                            println!("      No permutations agree outside the window!");
                            // Try a relaxed version: agree outside a LARGER window
                            // Maybe some values get rearranged further?
                            println!("      Trying value-based comparison...");
                            for (ia, info_a) in group_a.iter().enumerate() {
                                for (ib, info_b) in group_b.iter().enumerate() {
                                    // Which positions differ?
                                    let diffs: Vec<usize> = (0..info_a.pi.len())
                                        .filter(|&pos| info_a.pi[pos] != info_b.pi[pos])
                                        .collect();
                                    let delta =
                                        info_b.modified_swaps as i64 - info_a.modified_swaps as i64;
                                    if diffs.len() <= 4 {
                                        println!(
                                            "        {} <-> {}  diff_positions={:?}  delta_mod_swaps={}",
                                            perm_to_string(&info_a.pi),
                                            perm_to_string(&info_b.pi),
                                            diffs,
                                            delta,
                                        );
                                    }
                                }
                            }
                        }

                        // Test 4: For perms that agree outside window, is |delta| <= 1?
                        if !bijection_candidates.is_empty() {
                            let max_delta = bijection_candidates
                                .iter()
                                .map(|(_, _, d)| d.abs())
                                .max()
                                .unwrap();
                            println!(
                                "      Max |delta_mod_swaps| for window-matching pairs: {}",
                                max_delta
                            );
                        }
                    }
                }

                // Test 5: Look for a bijection with controlled delta
                println!(
                    "\n  --- Bijection analysis for (p_a={}, p_b={}) ---",
                    p_a, p_b
                );
                for &q in &all_qs {
                    let group_a = by_q_a.get(&q).cloned().unwrap_or_default();
                    let group_b = by_q_b.get(&q).cloned().unwrap_or_default();

                    if group_a.is_empty() || group_b.is_empty() {
                        if !group_a.is_empty() || !group_b.is_empty() {
                            println!(
                                "  q={}: |group_a|={}, |group_b|={} -- sizes differ, one empty",
                                q,
                                group_a.len(),
                                group_b.len()
                            );
                        }
                        continue;
                    }

                    println!(
                        "  q={}: |group_a|={}, |group_b|={}",
                        q,
                        group_a.len(),
                        group_b.len()
                    );

                    // Sorted mod_swaps values
                    let mut ms_a: Vec<usize> = group_a.iter().map(|p| p.modified_swaps).collect();
                    let mut ms_b: Vec<usize> = group_b.iter().map(|p| p.modified_swaps).collect();
                    ms_a.sort();
                    ms_b.sort();
                    println!("    mod_swaps_a (sorted): {:?}", ms_a);
                    println!("    mod_swaps_b (sorted): {:?}", ms_b);

                    // If same size, try sorted bijection (match sorted lists)
                    if ms_a.len() == ms_b.len() {
                        let deltas: Vec<i64> = ms_a
                            .iter()
                            .zip(ms_b.iter())
                            .map(|(&a, &b)| b as i64 - a as i64)
                            .collect();
                        let max_d = deltas.iter().map(|d| d.abs()).max().unwrap();
                        println!(
                            "    Sorted bijection deltas: {:?}  max_|delta|={}",
                            deltas, max_d
                        );
                        if max_d <= 1 {
                            println!("    ==> SORTED BIJECTION with |delta| <= 1 EXISTS!");
                        }
                    }

                    // Try greedy matching: for each a, find closest unmatched b
                    if ms_a.len() == ms_b.len() {
                        let mut used = vec![false; ms_b.len()];
                        let mut greedy_deltas = Vec::new();
                        let mut greedy_ok = true;
                        for &va in &ms_a {
                            let mut best_j = None;
                            let mut best_d = i64::MAX;
                            for (j, &vb) in ms_b.iter().enumerate() {
                                if used[j] {
                                    continue;
                                }
                                let d = (vb as i64 - va as i64).abs();
                                if d < best_d {
                                    best_d = d;
                                    best_j = Some(j);
                                }
                            }
                            if let Some(j) = best_j {
                                used[j] = true;
                                greedy_deltas.push(ms_b[j] as i64 - va as i64);
                                if best_d > 1 {
                                    greedy_ok = false;
                                }
                            }
                        }
                        if greedy_ok {
                            println!(
                                "    Greedy bijection deltas: {:?}  ALL |delta| <= 1",
                                greedy_deltas
                            );
                        }
                    }
                }

                // Overall compatibility check
                println!("\n  --- Compatibility summary ---");
                let mut all_compatible = true;
                for &q in &all_qs {
                    let group_a = by_q_a.get(&q).cloned().unwrap_or_default();
                    let group_b = by_q_b.get(&q).cloned().unwrap_or_default();
                    let vals_a: Vec<usize> = group_a.iter().map(|p| p.modified_swaps).collect();
                    let vals_b: Vec<usize> = group_b.iter().map(|p| p.modified_swaps).collect();
                    let poly_a = build_poly(&vals_a);
                    let poly_b = build_poly(&vals_b);
                    let compat = check_pairwise_compatible(&poly_a, &poly_b);
                    if !compat {
                        println!("  q={}: INCOMPATIBLE!", q);
                        all_compatible = false;
                    }
                }
                if all_compatible {
                    println!("  All q-values: COMPATIBLE for p_a={}, p_b={}", p_a, p_b);
                }
            }
        }

        // Summary: test ALL descent sets for this n
        println!("\n{}", "#".repeat(80));
        println!("# n={} FULL SUMMARY (all valid S with 1 not in S)", n);
        println!("{}", "#".repeat(80));

        let mut total_tested = 0u64;
        let mut total_compatible = 0u64;
        let mut total_same_size = 0u64;
        let mut total_sorted_bij_ok = 0u64;

        for s_mask in 0u64..(1 << (n - 1)) {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            for idx in 0..vp.len() - 1 {
                let p_a = vp[idx];
                let p_b = vp[idx + 1];

                let sp_a_asc = source_asc(s_mask, p_a, n);
                let sp_a_desc = source_desc(s_mask, p_a, n);
                let sp_b_asc = source_asc(s_mask, p_b, n);
                let sp_b_desc = source_desc(s_mask, p_b, n);

                // Collect permutations
                let mut by_q_a: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
                let mut by_q_b: BTreeMap<u8, Vec<usize>> = BTreeMap::new();

                // p_a sources
                if let Some(indices) = by_descent_prev.get(&sp_a_asc) {
                    for &i in indices {
                        let pi = perm_prev_list[i];
                        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                        let sw = compute(pi, Stat::Swaps);
                        let e1 = epsilon1(pi, p_a);
                        let ms = sw + if e1 { 1 } else { 0 };
                        by_q_a.entry(q).or_default().push(ms);
                    }
                }
                if let Some(sp) = sp_a_desc {
                    if let Some(indices) = by_descent_prev.get(&sp) {
                        for &i in indices {
                            let pi = perm_prev_list[i];
                            let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            let sw = compute(pi, Stat::Swaps);
                            by_q_a.entry(q).or_default().push(sw);
                        }
                    }
                }

                // p_b sources
                if let Some(indices) = by_descent_prev.get(&sp_b_asc) {
                    for &i in indices {
                        let pi = perm_prev_list[i];
                        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                        let sw = compute(pi, Stat::Swaps);
                        let e1 = epsilon1(pi, p_b);
                        let ms = sw + if e1 { 1 } else { 0 };
                        by_q_b.entry(q).or_default().push(ms);
                    }
                }
                if let Some(sp) = sp_b_desc {
                    if let Some(indices) = by_descent_prev.get(&sp) {
                        for &i in indices {
                            let pi = perm_prev_list[i];
                            let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            let sw = compute(pi, Stat::Swaps);
                            by_q_b.entry(q).or_default().push(sw);
                        }
                    }
                }

                let mut all_qs: Vec<u8> = by_q_a.keys().chain(by_q_b.keys()).copied().collect();
                all_qs.sort();
                all_qs.dedup();

                for &q in &all_qs {
                    let vals_a = by_q_a.get(&q).cloned().unwrap_or_default();
                    let vals_b = by_q_b.get(&q).cloned().unwrap_or_default();

                    if vals_a.is_empty() && vals_b.is_empty() {
                        continue;
                    }

                    let poly_a = build_poly(&vals_a);
                    let poly_b = build_poly(&vals_b);

                    total_tested += 1;
                    if check_pairwise_compatible(&poly_a, &poly_b) {
                        total_compatible += 1;
                    } else {
                        println!(
                            "  FAIL: S={} p_a={} p_b={} q={}  F_a={} F_b={}",
                            mask_to_string(s_mask, n),
                            p_a,
                            p_b,
                            q,
                            format_poly(&poly_a),
                            format_poly(&poly_b),
                        );
                    }

                    // Check sizes and sorted bijection
                    if vals_a.len() == vals_b.len() {
                        total_same_size += 1;
                        let mut sa = vals_a.clone();
                        let mut sb = vals_b.clone();
                        sa.sort();
                        sb.sort();
                        let max_d = sa
                            .iter()
                            .zip(sb.iter())
                            .map(|(&a, &b)| (b as i64 - a as i64).abs())
                            .max()
                            .unwrap_or(0);
                        if max_d <= 1 {
                            total_sorted_bij_ok += 1;
                        }
                    }
                }
            }
        }

        println!(
            "  Compatibility: {}/{} same-q pairs compatible",
            total_compatible, total_tested
        );
        println!(
            "  Same-size groups: {}/{} pairs have |group_a|=|group_b|",
            total_same_size, total_tested
        );
        println!(
            "  Sorted bijection |delta|<=1: {}/{}",
            total_sorted_bij_ok, total_same_size
        );
        println!();
    }
}
