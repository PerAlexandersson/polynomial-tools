/// Coupling / bijection approach for proving LR ordering of position-refined
/// swaps polynomials.
///
/// For S ⊆ [n-1] with 1 ∉ S and consecutive valid positions p_a < p_b:
///   [t^k] L^{(p_a)} · [t^{k+1}] L^{(p_b)} >= [t^{k+1}] L^{(p_a)} · [t^k] L^{(p_b)}
///
/// Strategy: try to find an explicit injection from sources for p_a to sources
/// for p_b that increases the swaps statistic. The source descent sets S'_{p_a}
/// and S'_{p_b} differ at specific positions. We try local operations (value swaps
/// at those positions) to convert a source for p_a into a valid source for p_b,
/// and track how the total swaps statistic changes.
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use std::collections::{BTreeMap, BTreeSet, HashMap};

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

/// Source descent set S'_p (ascent at p-1): the descent set of pi in S_{n-1}
/// such that inserting n at position p gives a permutation in D(n, S) with
/// an ascent at position p-1 (i.e., pi(p-1) < n).
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

/// Source descent set S''_p = S'_p ∪ {p-1} (descent at p-1).
fn source_desc(s_mask: u64, p: u8, n: u8) -> Option<u64> {
    if p <= 1 || p >= n {
        return None;
    }
    Some(source_asc(s_mask, p, n) | (1 << (p - 2)))
}

/// epsilon_1: does inserting n between pi[p-2] and pi[p-1] break the
/// consecutive pair (pi[p-2], pi[p-2]+1)?
fn epsilon1(pi: &[u8], p: u8) -> bool {
    let n = pi.len() as u8 + 1;
    if p <= 1 || p >= n {
        return false;
    }
    pi[(p - 2) as usize] + 1 == pi[(p - 1) as usize]
}

/// epsilon_2: does q < p-1 (position of n-1 is to the left of insertion point)?
/// When true, inserting n at position p creates a new long swap for (n-1, n).
fn epsilon2(p: u8, q: u8) -> bool {
    // q is 1-indexed position of n-1 in pi (of length n-1)
    // We need q < p-1 (with appropriate adjustment)
    q < p - 1
}

/// Insert value n at 1-indexed position p into pi in S_{n-1}.
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

/// Check LR ordering: for all k, a[k]*b[k+1] >= a[k+1]*b[k]
fn check_lr(a: &[i64], b: &[i64]) -> bool {
    let max_deg = a.len().max(b.len());
    for k in 0..max_deg {
        if coeff(a, k) * coeff(b, k + 1) < coeff(a, k + 1) * coeff(b, k) {
            return false;
        }
    }
    true
}

/// Find positions where two bitmasks differ (1-indexed positions in [1..n-1]).
fn diff_positions(a: u64, b: u64, n: u8) -> Vec<(u8, bool, bool)> {
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

/// Compute the "modified swaps" for source pi at insertion position p:
///   swaps(insert(n, p, pi)) = swaps(pi) + eps1(pi, p) + eps2(p, pos_{n-1}(pi))
fn modified_swaps(pi: &[u8], p: u8) -> usize {
    let n = pi.len() as u8 + 1;
    let base = compute(pi, Stat::Swaps);
    let e1 = if epsilon1(pi, p) { 1 } else { 0 };
    // Position of n-1 in pi (1-indexed)
    let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
    let e2 = if p > 1 && q < p - 1 { 1 } else { 0 };
    base + e1 + e2
}

/// Try swapping values at 0-indexed positions i and j in pi.
fn swap_values_at(pi: &[u8], i: usize, j: usize) -> Vec<u8> {
    let mut result = pi.to_vec();
    result.swap(i, j);
    result
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    println!("=== Coupling experiment for LR ordering ===\n");

    for n in 5..=max_n {
        let perms_prev = all_permutations(n - 1);

        let mut prev_by_des: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        // Enumerate all S with 1 not in S
        let max_mask = 1u64 << (n - 1);
        let mut total_pairs = 0u64;
        let mut lr_ok = 0u64;

        let mut total_coupling_attempts = 0u64;
        let mut naive_swap_works = 0u64;
        let mut naive_swap_valid = 0u64;
        let mut naive_swap_monotone = 0u64;
        let mut best_local_works = 0u64;

        for s_mask in 0..max_mask {
            // 1 not in S
            if s_mask & 1 != 0 {
                continue;
            }

            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            // Collect sources for each valid position
            let mut sources: BTreeMap<u8, Vec<Vec<u8>>> = BTreeMap::new();
            for &p in &vp {
                let spa = source_asc(s_mask, p, n);
                let spd = source_desc(s_mask, p, n);
                let mut src = Vec::new();
                if let Some(cls) = prev_by_des.get(&spa) {
                    src.extend(cls.iter().cloned());
                }
                if let Some(sd) = spd {
                    if sd != spa {
                        if let Some(cls) = prev_by_des.get(&sd) {
                            src.extend(cls.iter().cloned());
                        }
                    }
                }
                sources.insert(p, src);
            }

            // Build L^{(p)} polynomials
            let mut polys: BTreeMap<u8, Vec<i64>> = BTreeMap::new();
            for &p in &vp {
                let vals: Vec<usize> = sources[&p].iter().map(|pi| modified_swaps(pi, p)).collect();
                polys.insert(p, build_poly(&vals));
            }

            // Check LR for consecutive pairs
            for w in vp.windows(2) {
                let pa = w[0];
                let pb = w[1];
                total_pairs += 1;

                let poly_a = &polys[&pa];
                let poly_b = &polys[&pb];

                if check_lr(poly_a, poly_b) {
                    lr_ok += 1;
                } else {
                    println!(
                        "  LR FAIL: n={} S={} pa={} pb={}",
                        n,
                        descent_set_to_string(s_mask, n),
                        pa,
                        pb
                    );
                }

                // ============================================================
                // Now try the coupling: map sources for p_a to sources for p_b
                // ============================================================

                // Compute source descent sets and their differences
                let spa_asc = source_asc(s_mask, pa, n);
                let spa_desc = source_desc(s_mask, pa, n);
                let spb_asc = source_asc(s_mask, pb, n);
                let spb_desc = source_desc(s_mask, pb, n);

                // Source sets for pa come from D(n-1, S'_pa) ∪ D(n-1, S''_pa)
                // Source sets for pb come from D(n-1, S'_pb) ∪ D(n-1, S''_pb)
                // Find where S'_pa and S'_pb differ:
                let diff_asc = diff_positions(spa_asc, spb_asc, n - 1);

                let diff_desc = match (spa_desc, spb_desc) {
                    (Some(a), Some(b)) => diff_positions(a, b, n - 1),
                    _ => vec![],
                };

                let src_a = &sources[&pa];
                let src_b_set: BTreeSet<Vec<u8>> = sources[&pb].iter().cloned().collect();

                // For each source pi_a for p_a, try to map it to a source pi_b for p_b
                let mut pair_attempts = 0u64;
                let mut pair_naive_works = 0u64;
                let mut pair_naive_valid = 0u64;
                let mut pair_naive_monotone = 0u64;
                let mut pair_best_local = 0u64;

                // Track the image multiset for injectivity checking
                let mut naive_images: Vec<Option<Vec<u8>>> = Vec::new();
                let mut best_images: Vec<Option<(Vec<u8>, String)>> = Vec::new();

                for pi_a in src_a {
                    pair_attempts += 1;
                    let swaps_a = modified_swaps(pi_a, pa);

                    // --- Strategy 1: Naive swap at each differing position ---
                    // The source descent sets differ at specific positions.
                    // At each differing position j* (0-indexed j*-1):
                    //   if j* is in S'_pb but not in S'_pa, pi_a has ASCENT at j*,
                    //   we need DESCENT at j* => swap values at positions j*-1 and j*.

                    let des_a = descent_set_bitmask(pi_a);

                    // Try: apply all diffs from S'_asc simultaneously
                    let mut candidate = pi_a.clone();
                    let mut all_diffs_applicable = true;

                    // Collect the diff positions between all source sets
                    // Combined: diffs between (S'_pa ∪ S''_pa) and (S'_pb ∪ S''_pb)
                    // But actually, a single pi_a belongs to exactly one of S'_pa or S''_pa.
                    let pi_a_des = descent_set_bitmask(pi_a);
                    let pi_a_source_type = if pi_a_des == spa_asc { "asc" } else { "desc" };

                    // Determine which source descent set pi_a belongs to
                    let pi_a_source_mask = pi_a_des;

                    // Target: which source descent sets for pb could pi_b belong to?
                    // It could be S'_pb (asc) or S''_pb (desc).

                    // Strategy: find diff positions between pi_a's descent set and each target,
                    // and try the simplest operation at each.

                    // First, try: swap values at positions (j*-1, j*) for each diff position j*
                    // (0-indexed) between pi_a's descent set and spb_asc.
                    let target_masks: Vec<(u64, &str)> = {
                        let mut v = vec![(spb_asc, "asc")];
                        if let Some(sd) = spb_desc {
                            v.push((sd, "desc"));
                        }
                        v
                    };

                    let mut found_naive = false;
                    let mut found_any = false;
                    let mut best_op_desc = String::new();
                    let mut best_candidate: Option<Vec<u8>> = None;
                    let mut best_swaps_b = 0usize;

                    for &(target_mask, _ttype) in &target_masks {
                        // Compute diff between pi_a's actual descent set and target
                        let diffs = diff_positions(pi_a_des, target_mask, n - 1);
                        if diffs.is_empty() {
                            // pi_a already has the target descent set!
                            if src_b_set.contains(pi_a) {
                                let swaps_b = modified_swaps(pi_a, pb);
                                if !found_any || swaps_b > best_swaps_b {
                                    best_candidate = Some(pi_a.clone());
                                    best_swaps_b = swaps_b;
                                    best_op_desc = "identity".to_string();
                                    found_any = true;
                                }
                                found_naive = true;
                            }
                            continue;
                        }

                        // Try naive: swap at each diff position
                        let mut c = pi_a.clone();
                        let mut naive_ok = true;
                        for &(pos, in_a_des, in_target) in &diffs {
                            // pos is 1-indexed position in [1..n-2]
                            // descent at position pos means c[pos-1] > c[pos] (0-indexed)
                            let idx = (pos - 1) as usize;
                            if idx + 1 >= c.len() {
                                naive_ok = false;
                                break;
                            }
                            if in_a_des && !in_target {
                                // Need to remove descent at pos: make c[idx] < c[idx+1]
                                // Swap if c[idx] > c[idx+1]
                                if c[idx] > c[idx + 1] {
                                    c.swap(idx, idx + 1);
                                } else {
                                    naive_ok = false;
                                    break;
                                }
                            } else if !in_a_des && in_target {
                                // Need to create descent at pos: make c[idx] > c[idx+1]
                                if c[idx] < c[idx + 1] {
                                    c.swap(idx, idx + 1);
                                } else {
                                    naive_ok = false;
                                    break;
                                }
                            }
                        }

                        if naive_ok {
                            let c_des = descent_set_bitmask(&c);
                            if c_des == target_mask && src_b_set.contains(&c) {
                                let swaps_b = modified_swaps(&c, pb);
                                found_naive = true;
                                if !found_any || swaps_b > best_swaps_b {
                                    best_candidate = Some(c.clone());
                                    best_swaps_b = swaps_b;
                                    best_op_desc = format!("naive_swap(diffs={:?})", diffs);
                                    found_any = true;
                                }
                            }
                        }
                    }

                    if found_naive {
                        pair_naive_valid += 1;
                        let swaps_b = best_swaps_b;
                        if swaps_b >= swaps_a {
                            pair_naive_monotone += 1;
                        }
                        pair_naive_works += 1;
                        naive_images.push(best_candidate.clone());
                    } else {
                        naive_images.push(None);
                    }

                    // --- Strategy 2: Try ALL local operations ---
                    // For each single swap at positions (i, i+1) for i in 0..n-3,
                    // check if it gives a valid source for pb with swaps >= swaps_a.
                    if !found_any {
                        let nn = pi_a.len();
                        // Try single adjacent swaps
                        for i in 0..nn.saturating_sub(1) {
                            let c = swap_values_at(pi_a, i, i + 1);
                            let c_des = descent_set_bitmask(&c);
                            if src_b_set.contains(&c) {
                                let swaps_b = modified_swaps(&c, pb);
                                if !found_any || swaps_b > best_swaps_b {
                                    best_candidate = Some(c.clone());
                                    best_swaps_b = swaps_b;
                                    best_op_desc = format!("adj_swap({})", i);
                                    found_any = true;
                                }
                            }
                        }
                        // Try single non-adjacent swaps (i, j)
                        if !found_any {
                            for i in 0..nn {
                                for j in i + 2..nn {
                                    let c = swap_values_at(pi_a, i, j);
                                    let c_des = descent_set_bitmask(&c);
                                    if src_b_set.contains(&c) {
                                        let swaps_b = modified_swaps(&c, pb);
                                        if !found_any || swaps_b > best_swaps_b {
                                            best_candidate = Some(c.clone());
                                            best_swaps_b = swaps_b;
                                            best_op_desc = format!("nonadj_swap({},{})", i, j);
                                            found_any = true;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if found_any {
                        pair_best_local += 1;
                        best_images.push(Some((best_candidate.unwrap(), best_op_desc.clone())));
                    } else {
                        best_images.push(None);
                    }
                }

                total_coupling_attempts += pair_attempts;
                naive_swap_works += pair_naive_works;
                naive_swap_valid += pair_naive_valid;
                naive_swap_monotone += pair_naive_monotone;
                best_local_works += pair_best_local;

                // Check injectivity of the naive map
                let naive_injective = {
                    let images: Vec<&Vec<u8>> =
                        naive_images.iter().filter_map(|x| x.as_ref()).collect();
                    let unique: BTreeSet<&Vec<u8>> = images.iter().copied().collect();
                    images.len() == unique.len()
                };

                // Check injectivity of best local map
                let best_injective = {
                    let images: Vec<&Vec<u8>> = best_images
                        .iter()
                        .filter_map(|x| x.as_ref().map(|(v, _)| v))
                        .collect();
                    let unique: BTreeSet<&Vec<u8>> = images.iter().copied().collect();
                    images.len() == unique.len()
                };

                // Check monotonicity: does every mapped source have swaps_b >= swaps_a?
                let mut all_monotone = true;
                let mut monotone_count = 0u64;
                let mut strict_increase = 0u64;
                let mut delta_histogram: BTreeMap<i64, u64> = BTreeMap::new();

                for (idx, pi_a) in src_a.iter().enumerate() {
                    let swaps_a = modified_swaps(pi_a, pa);
                    if let Some((ref pi_b, _)) = best_images[idx] {
                        let swaps_b = modified_swaps(pi_b, pb);
                        let delta = swaps_b as i64 - swaps_a as i64;
                        *delta_histogram.entry(delta).or_insert(0) += 1;
                        if swaps_b >= swaps_a {
                            monotone_count += 1;
                        } else {
                            all_monotone = false;
                        }
                        if swaps_b > swaps_a {
                            strict_increase += 1;
                        }
                    }
                }

                // Report for this pair (only if interesting or n small)
                if n <= 6 || !all_monotone || pair_best_local < pair_attempts {
                    let s_str = descent_set_to_string(s_mask, n);
                    println!(
                        "n={} S={} (pa={}, pb={}): sources_a={}, sources_b={}",
                        n,
                        s_str,
                        pa,
                        pb,
                        src_a.len(),
                        sources[&pb].len()
                    );

                    // Show source set diffs
                    if !diff_asc.is_empty() {
                        println!("  Diff S'_asc: {:?}", diff_asc);
                    }
                    if !diff_desc.is_empty() {
                        println!("  Diff S''_desc: {:?}", diff_desc);
                    }

                    println!(
                        "  Naive swap: {}/{} valid, {}/{} monotone(swaps_b>=swaps_a), injective={}",
                        pair_naive_valid,
                        pair_attempts,
                        pair_naive_monotone,
                        pair_naive_valid,
                        naive_injective
                    );
                    println!(
                        "  Best local op: {}/{} found, injective={}, all_monotone={}",
                        pair_best_local, pair_attempts, best_injective, all_monotone
                    );
                    if !delta_histogram.is_empty() {
                        println!(
                            "  Delta histogram (swaps_b - swaps_a): {:?}",
                            delta_histogram
                        );
                    }

                    // Show unmapped sources
                    let unmapped: Vec<usize> = (0..src_a.len())
                        .filter(|&i| best_images[i].is_none())
                        .collect();
                    if !unmapped.is_empty() && unmapped.len() <= 5 {
                        println!("  Unmapped sources:");
                        for &idx in &unmapped {
                            let pi = &src_a[idx];
                            let sw = modified_swaps(pi, pa);
                            println!(
                                "    pi={:?} swaps_a={} des={}",
                                pi,
                                sw,
                                descent_set_to_string(descent_set_bitmask(pi), n - 1)
                            );
                        }
                    } else if !unmapped.is_empty() {
                        println!("  Unmapped: {} sources", unmapped.len());
                    }

                    // Show collisions if not injective
                    if !best_injective {
                        let mut image_count: BTreeMap<&Vec<u8>, Vec<usize>> = BTreeMap::new();
                        for (idx, img) in best_images.iter().enumerate() {
                            if let Some((ref v, _)) = img {
                                image_count.entry(v).or_default().push(idx);
                            }
                        }
                        let collisions: Vec<_> = image_count
                            .iter()
                            .filter(|(_, idxs)| idxs.len() > 1)
                            .collect();
                        if collisions.len() <= 3 {
                            for (img, idxs) in &collisions {
                                println!(
                                    "  Collision: pi_b={:?} from sources {:?}",
                                    img,
                                    idxs.iter()
                                        .map(|&i| format!("{:?}", src_a[i]))
                                        .collect::<Vec<_>>()
                                );
                            }
                        } else {
                            println!("  {} collisions", collisions.len());
                        }
                    }

                    println!();
                }
            }
        }

        println!(
            "--- n={}: LR ok {}/{}, naive_swap {}/{}, best_local {}/{} ---\n",
            n,
            lr_ok,
            total_pairs,
            naive_swap_works,
            total_coupling_attempts,
            best_local_works,
            total_coupling_attempts
        );
    }

    // ========================================================================
    // Part 2: More detailed analysis for small cases
    // Show the actual permutation-level mapping for n=5
    // ========================================================================
    println!("\n=== Part 2: Detailed mapping for n=5 ===\n");

    let n: u8 = 5;
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

        let mut sources: BTreeMap<u8, Vec<Vec<u8>>> = BTreeMap::new();
        for &p in &vp {
            let spa = source_asc(s_mask, p, n);
            let spd = source_desc(s_mask, p, n);
            let mut src = Vec::new();
            if let Some(cls) = prev_by_des.get(&spa) {
                src.extend(cls.iter().cloned());
            }
            if let Some(sd) = spd {
                if sd != spa {
                    if let Some(cls) = prev_by_des.get(&sd) {
                        src.extend(cls.iter().cloned());
                    }
                }
            }
            sources.insert(p, src);
        }

        let s_str = descent_set_to_string(s_mask, n);

        for w in vp.windows(2) {
            let pa = w[0];
            let pb = w[1];

            let src_a = &sources[&pa];
            let src_b = &sources[&pb];
            let src_b_set: BTreeSet<Vec<u8>> = src_b.iter().cloned().collect();

            println!(
                "S={} pa={} pb={}: |src_a|={} |src_b|={}",
                s_str,
                pa,
                pb,
                src_a.len(),
                src_b.len()
            );

            // Show the full insertion picture
            for pi_a in src_a {
                let swaps_a = modified_swaps(pi_a, pa);
                let sigma_a = insert(pi_a, pa as usize);
                let sigma_a_swaps = compute(&sigma_a, Stat::Swaps);

                // What happens if we insert n at pb into pi_a?
                let sigma_cross = insert(pi_a, pb as usize);
                let sigma_cross_des = descent_set_bitmask(&sigma_cross);
                let sigma_cross_swaps = compute(&sigma_cross, Stat::Swaps);
                let cross_valid = sigma_cross_des == s_mask;

                // Try: find pi_b in src_b such that swaps(insert(n, pb, pi_b)) >= swaps(sigma_a)
                // using the nearest-swaps heuristic
                let mut best_match: Option<(&Vec<u8>, usize, i64)> = None;
                for pi_b in src_b.iter() {
                    let swaps_b = modified_swaps(pi_b, pb);
                    let delta = swaps_b as i64 - swaps_a as i64;
                    if best_match.is_none()
                        || delta.abs() < best_match.unwrap().2.abs()
                        || (delta.abs() == best_match.unwrap().2.abs() && delta >= 0)
                    {
                        best_match = Some((pi_b, swaps_b, delta));
                    }
                }

                if let Some((pi_b, swaps_b, delta)) = best_match {
                    println!(
                        "  pi_a={:?} sw_a={} -> cross_valid={} | best pi_b={:?} sw_b={} delta={}",
                        pi_a, swaps_a, cross_valid, pi_b, swaps_b, delta
                    );
                }
            }
            println!();
        }
    }

    // ========================================================================
    // Part 3: Analysis of the "cross-insertion" approach
    // For each sigma in D(n,S) with pos(n)=pa, move n to position pb.
    // Does the result land in D(n,S)? If not, what goes wrong?
    // ========================================================================
    println!("\n=== Part 3: Cross-insertion analysis ===\n");

    for n in 5..=max_n.min(7) {
        let perms = all_permutations(n);
        let perms_prev = all_permutations(n - 1);

        let mut prev_by_des: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        let max_mask = 1u64 << (n - 1);

        let mut cross_total = 0u64;
        let mut cross_valid = 0u64;
        let mut cross_monotone = 0u64;

        for s_mask in 0..max_mask {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            let mut sources: BTreeMap<u8, Vec<Vec<u8>>> = BTreeMap::new();
            for &p in &vp {
                let spa = source_asc(s_mask, p, n);
                let spd = source_desc(s_mask, p, n);
                let mut src = Vec::new();
                if let Some(cls) = prev_by_des.get(&spa) {
                    src.extend(cls.iter().cloned());
                }
                if let Some(sd) = spd {
                    if sd != spa {
                        if let Some(cls) = prev_by_des.get(&sd) {
                            src.extend(cls.iter().cloned());
                        }
                    }
                }
                sources.insert(p, src);
            }

            for w in vp.windows(2) {
                let pa = w[0];
                let pb = w[1];

                for pi_a in &sources[&pa] {
                    cross_total += 1;
                    let swaps_a = modified_swaps(pi_a, pa);

                    // Direct cross-insertion: insert n at pb into pi_a
                    let sigma_cross = insert(pi_a, pb as usize);
                    let sigma_cross_des = descent_set_bitmask(&sigma_cross);

                    if sigma_cross_des == s_mask {
                        cross_valid += 1;
                        let swaps_cross = compute(&sigma_cross, Stat::Swaps);
                        if swaps_cross >= swaps_a {
                            cross_monotone += 1;
                        }
                    }
                }
            }
        }

        println!(
            "n={}: cross-insertion valid {}/{} ({:.1}%), monotone {}/{}",
            n,
            cross_valid,
            cross_total,
            100.0 * cross_valid as f64 / cross_total.max(1) as f64,
            cross_monotone,
            cross_valid
        );
    }

    // ========================================================================
    // Part 4: Analysis of how source sets differ and what operations fix them
    // For each (n, S, pa, pb), catalog the exact diff positions and try
    // systematic operations.
    // ========================================================================
    println!("\n=== Part 4: Source set diff structure ===\n");

    for n in 5..=max_n.min(7) {
        let perms_prev = all_permutations(n - 1);
        let mut prev_by_des: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        let max_mask = 1u64 << (n - 1);
        let mut diff_pattern_counts: BTreeMap<String, u64> = BTreeMap::new();

        for s_mask in 0..max_mask {
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            for w in vp.windows(2) {
                let pa = w[0];
                let pb = w[1];

                let spa_asc = source_asc(s_mask, pa, n);
                let spb_asc = source_asc(s_mask, pb, n);
                let spa_desc = source_desc(s_mask, pa, n);
                let spb_desc = source_desc(s_mask, pb, n);

                let diff_a = diff_positions(spa_asc, spb_asc, n - 1);

                let diff_d = match (spa_desc, spb_desc) {
                    (Some(a), Some(b)) => diff_positions(a, b, n - 1),
                    _ => vec![],
                };

                // Characterize the diff pattern
                let pattern = format!(
                    "asc_diffs={}, desc_diffs={}, pa={}, pb={}, gap={}",
                    diff_a.len(),
                    diff_d.len(),
                    pa,
                    pb,
                    pb - pa
                );
                *diff_pattern_counts.entry(pattern).or_insert(0) += 1;

                // Print detailed diff for small n
                if n <= 6 {
                    let s_str = descent_set_to_string(s_mask, n);
                    println!("n={} S={} pa={} pb={}:", n, s_str, pa, pb);
                    println!(
                        "  S'_asc(pa)={} S'_asc(pb)={}",
                        descent_set_to_string(spa_asc, n - 1),
                        descent_set_to_string(spb_asc, n - 1)
                    );
                    if let (Some(a), Some(b)) = (spa_desc, spb_desc) {
                        println!(
                            "  S''_desc(pa)={} S''_desc(pb)={}",
                            descent_set_to_string(a, n - 1),
                            descent_set_to_string(b, n - 1)
                        );
                    }
                    for (pos, in_a, in_b) in &diff_a {
                        let dir = if *in_a { "asc_only_a" } else { "asc_only_b" };
                        println!("  Diff at pos {}: {}", pos, dir);
                    }
                    println!();
                }
            }
        }

        println!("n={} diff pattern summary:", n);
        for (pat, cnt) in &diff_pattern_counts {
            println!("  {} : {} cases", pat, cnt);
        }
        println!();
    }
}
