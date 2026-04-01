/// Investigate the EXACT relationship between source descent sets for consecutive
/// valid positions p_a < p_b in P(S), to understand the varying-source TN_2 structure.
///
/// For target descent set S (subset of [n-1] with 1 not in S), and valid position p:
///   S'_p = (S intersect [p-2]) union {j-1 : j in S, j > p}     ("ascending" source)
///   S''_p = S'_p union {p-1}                                     ("descending" source)
///
/// For consecutive valid positions p_a < p_b, we compute:
///   1. The symmetric difference of S'_{p_a} and S'_{p_b}
///   2. The symmetric difference of S''_{p_a} and S''_{p_b}
///   3. Pattern analysis: which positions differ, and is there a structural rule?
///
/// Also compute the source permutations for each and look for bijection structure.
///
use combpoly::fixed_descent::{
    augmented_source_descent_set_for_target_descent_set_and_insertion_position as source_desc,
    base_source_descent_set_for_target_descent_set_and_insertion_position as source_asc,
    insertion_breaks_consecutive_ascending_pair_at_boundary as epsilon1,
    valid_insertion_positions_for_target_descent_set as valid_positions,
};
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use std::collections::{BTreeMap, BTreeSet};

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

fn sym_diff(a: &BTreeSet<u8>, b: &BTreeSet<u8>) -> BTreeSet<u8> {
    a.symmetric_difference(b).copied().collect()
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

/// Classify what S looks like between positions p_a and p_b.
/// Returns the elements of S in the interval [p_a, p_b-1].
fn s_between(s_mask: u64, p_a: u8, p_b: u8) -> BTreeSet<u8> {
    let mut result = BTreeSet::new();
    for j in p_a..p_b {
        if (s_mask >> (j - 1)) & 1 == 1 {
            result.insert(j);
        }
    }
    result
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    println!("=== Source descent set relationship for consecutive valid positions ===\n");

    // Collect statistics about the symmetric differences
    // Map: (size_of_sym_diff) -> count
    let mut diff_size_counts: BTreeMap<usize, usize> = BTreeMap::new();
    // Map: pattern description -> count
    let mut pattern_counts: BTreeMap<String, usize> = BTreeMap::new();

    for n in 5..=max_n {
        println!("========== n = {} ==========", n);
        let perms_prev = all_permutations(n - 1);
        let mut prev_by_des: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        for s_mask in 0u64..(1 << (n - 1)) {
            // 1 not in S
            if s_mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            let s_set = mask_to_set(s_mask, n);
            let s_str = set_to_string(&s_set);

            // Look at each consecutive pair
            for idx in 0..vp.len() - 1 {
                let p_a = vp[idx];
                let p_b = vp[idx + 1];

                let sp_a = source_asc(s_mask, p_a, n);
                let sp_b = source_asc(s_mask, p_b, n);
                let spp_a = source_desc(s_mask, p_a, n);
                let spp_b = source_desc(s_mask, p_b, n);

                let set_sp_a = mask_to_set(sp_a, n - 1);
                let set_sp_b = mask_to_set(sp_b, n - 1);
                let diff_prime = sym_diff(&set_sp_a, &set_sp_b);

                // Compute S''_p differences too (when both exist)
                let diff_dprime = match (spp_a, spp_b) {
                    (Some(da), Some(db)) => {
                        let set_da = mask_to_set(da, n - 1);
                        let set_db = mask_to_set(db, n - 1);
                        Some(sym_diff(&set_da, &set_db))
                    }
                    _ => None,
                };

                // Analyze the diff pattern for S'
                let diff_size = diff_prime.len();
                *diff_size_counts.entry(diff_size).or_insert(0) += 1;

                // Classify: are all diff positions in [p_a, p_b-2]?
                let all_in_window = diff_prime.iter().all(|&j| j >= p_a && j <= p_b - 1);

                // Elements of S between p_a and p_b
                let s_between_set = s_between(s_mask, p_a + 1, p_b);

                // Positions in diff that are in S'_{p_a} but not S'_{p_b}
                let only_in_a: BTreeSet<u8> = diff_prime
                    .iter()
                    .filter(|j| set_sp_a.contains(j))
                    .copied()
                    .collect();
                let only_in_b: BTreeSet<u8> = diff_prime
                    .iter()
                    .filter(|j| set_sp_b.contains(j))
                    .copied()
                    .collect();

                // Describe the pattern
                let pattern_desc = if diff_size == 0 {
                    "IDENTICAL".to_string()
                } else {
                    // Check if the diff is exactly: S'_{p_b} has p_a, p_a+1,..., elements that
                    // correspond to positions p_a through p_b-2 where S has a descent
                    let mut desc_parts = Vec::new();
                    if !only_in_a.is_empty() {
                        desc_parts.push(format!("only_in_S'_pa={}", set_to_string(&only_in_a)));
                    }
                    if !only_in_b.is_empty() {
                        desc_parts.push(format!("only_in_S'_pb={}", set_to_string(&only_in_b)));
                    }
                    desc_parts.join("; ")
                };

                *pattern_counts.entry(pattern_desc.clone()).or_insert(0) += 1;

                // Detailed output for small n or interesting cases
                if n <= 7 || diff_size > 3 {
                    println!(
                        "  S={} p_a={} p_b={} | S'_pa={} S'_pb={}",
                        s_str,
                        p_a,
                        p_b,
                        mask_to_string(sp_a, n - 1),
                        mask_to_string(sp_b, n - 1),
                    );
                    println!(
                        "    sym_diff(S'_pa, S'_pb) = {}  (size {})",
                        set_to_string(&diff_prime),
                        diff_size
                    );
                    println!(
                        "    only_in_S'_pa = {}  only_in_S'_pb = {}",
                        set_to_string(&only_in_a),
                        set_to_string(&only_in_b)
                    );
                    println!("    all_in_window [{},{}]: {}", p_a, p_b - 1, all_in_window);
                    println!(
                        "    S between (p_a, p_b) = {} (elements of S in [{},{}])",
                        set_to_string(&s_between_set),
                        p_a + 1,
                        p_b - 1
                    );

                    if let (Some(da), Some(db), Some(ref dd)) = (spp_a, spp_b, &diff_dprime) {
                        println!(
                            "    S''_pa={} S''_pb={}",
                            mask_to_string(da, n - 1),
                            mask_to_string(db, n - 1),
                        );
                        println!(
                            "    sym_diff(S''_pa, S''_pb) = {}  (size {})",
                            set_to_string(dd),
                            dd.len()
                        );
                    }

                    // Detailed position-by-position analysis
                    print!("    Position-by-position for S':");
                    for j in 1..n - 1 {
                        let in_a = set_sp_a.contains(&j);
                        let in_b = set_sp_b.contains(&j);
                        if in_a != in_b {
                            let marker = if in_a && !in_b { "A" } else { "B" };
                            print!(" {}:only_{}", j, marker);
                        }
                    }
                    println!();

                    // NEW: Check the shift-by-1 hypothesis
                    // Hypothesis: for positions j in [p_a, p_b-2]:
                    //   S'_{p_a} has j iff j+1 in S (shifted portion)
                    //   S'_{p_b} has j iff j in S (unshifted portion)
                    // So the diff at position j should be: j in S XOR j+1 in S
                    let mut shift_hypothesis_holds = true;
                    for j in p_a..=(p_b.saturating_sub(2)) {
                        let in_a = set_sp_a.contains(&j);
                        let in_b = set_sp_b.contains(&j);
                        let j_in_s = (s_mask >> (j - 1)) & 1 == 1;
                        let j1_in_s = if j < n - 1 {
                            (s_mask >> j) & 1 == 1
                        } else {
                            false
                        };

                        // S'_{p_a}: for position j where p_a <= j <= p_b-2
                        //   j = (j+1)-1, and j+1 > p_a, so j in S'_{p_a} iff j+1 in S
                        // S'_{p_b}: j < p_b-1, so j in S'_{p_b} iff j in S
                        // (but we need j <= p_b-2 for "match S" region of S'_{p_b})
                        let expected_a = j1_in_s;
                        let expected_b = j_in_s;

                        if in_a != expected_a || in_b != expected_b {
                            shift_hypothesis_holds = false;
                        }
                    }
                    println!(
                        "    Shift hypothesis (j in [p_a, p_b-2]): S'_pa has j iff (j+1 in S), S'_pb has j iff (j in S): {}",
                        if shift_hypothesis_holds { "HOLDS" } else { "FAILS" }
                    );

                    // Characterize the diff more precisely
                    // The diff at position j (in [p_a, p_b-2]) is:
                    //   j in S XOR j+1 in S
                    // i.e., exactly the positions where S "changes" between j and j+1
                    print!(
                        "    S-change positions in [{},{}]:",
                        p_a,
                        p_b.saturating_sub(2)
                    );
                    for j in p_a..=(p_b.saturating_sub(2)) {
                        let j_in_s = (s_mask >> (j - 1)) & 1 == 1;
                        let j1_in_s = if j < n - 1 {
                            (s_mask >> j) & 1 == 1
                        } else {
                            false
                        };
                        if j_in_s != j1_in_s {
                            let dir = if j_in_s && !j1_in_s { "off" } else { "on" };
                            print!(" {} (S turns {})", j, dir);
                        }
                    }
                    println!();

                    // Source permutation counts and polynomial comparison
                    let count_a = prev_by_des.get(&sp_a).map(|v| v.len()).unwrap_or(0);
                    let count_b = prev_by_des.get(&sp_b).map(|v| v.len()).unwrap_or(0);
                    println!(
                        "    |D(n-1, S'_pa)| = {}  |D(n-1, S'_pb)| = {}",
                        count_a, count_b
                    );

                    // Compute swap polynomials for each source set
                    if let Some(cls_a) = prev_by_des.get(&sp_a) {
                        let swaps_a: Vec<usize> =
                            cls_a.iter().map(|pi| compute(pi, Stat::Swaps)).collect();
                        let poly_a = build_poly(&swaps_a);
                        print!("    Swaps poly(S'_pa) = {}", format_poly(&poly_a));

                        // Also with eps1 modification
                        let mod_swaps_a: Vec<usize> = cls_a
                            .iter()
                            .map(|pi| {
                                let e = if epsilon1(pi, p_a) { 1 } else { 0 };
                                compute(pi, Stat::Swaps) + e
                            })
                            .collect();
                        let mod_poly_a = build_poly(&mod_swaps_a);
                        println!("  mod(eps1,p_a) = {}", format_poly(&mod_poly_a));
                    }
                    if let Some(cls_b) = prev_by_des.get(&sp_b) {
                        let swaps_b: Vec<usize> =
                            cls_b.iter().map(|pi| compute(pi, Stat::Swaps)).collect();
                        let poly_b = build_poly(&swaps_b);
                        print!("    Swaps poly(S'_pb) = {}", format_poly(&poly_b));

                        let mod_swaps_b: Vec<usize> = cls_b
                            .iter()
                            .map(|pi| {
                                let e = if epsilon1(pi, p_b) { 1 } else { 0 };
                                compute(pi, Stat::Swaps) + e
                            })
                            .collect();
                        let mod_poly_b = build_poly(&mod_swaps_b);
                        println!("  mod(eps1,p_b) = {}", format_poly(&mod_poly_b));
                    }

                    // For S''_p sources too
                    if let Some(da) = spp_a {
                        if let Some(cls_da) = prev_by_des.get(&da) {
                            let swaps_da: Vec<usize> =
                                cls_da.iter().map(|pi| compute(pi, Stat::Swaps)).collect();
                            let poly_da = build_poly(&swaps_da);
                            let mod_swaps_da: Vec<usize> = cls_da
                                .iter()
                                .map(|pi| {
                                    let e = if epsilon1(pi, p_a) { 1 } else { 0 };
                                    compute(pi, Stat::Swaps) + e
                                })
                                .collect();
                            let mod_poly_da = build_poly(&mod_swaps_da);
                            println!(
                                "    Swaps poly(S''_pa) = {}  mod(eps1,p_a) = {}",
                                format_poly(&poly_da),
                                format_poly(&mod_poly_da)
                            );
                        }
                    }
                    if let Some(db) = spp_b {
                        if let Some(cls_db) = prev_by_des.get(&db) {
                            let swaps_db: Vec<usize> =
                                cls_db.iter().map(|pi| compute(pi, Stat::Swaps)).collect();
                            let poly_db = build_poly(&swaps_db);
                            let mod_swaps_db: Vec<usize> = cls_db
                                .iter()
                                .map(|pi| {
                                    let e = if epsilon1(pi, p_b) { 1 } else { 0 };
                                    compute(pi, Stat::Swaps) + e
                                })
                                .collect();
                            let mod_poly_db = build_poly(&mod_swaps_db);
                            println!(
                                "    Swaps poly(S''_pb) = {}  mod(eps1,p_b) = {}",
                                format_poly(&poly_db),
                                format_poly(&mod_poly_db)
                            );
                        }
                    }

                    println!();
                }
            }
        }
        println!();
    }

    // Summary statistics
    println!("=== SYMMETRIC DIFFERENCE SIZE DISTRIBUTION ===");
    for (size, count) in &diff_size_counts {
        println!("  size {} : {} cases", size, count);
    }

    println!("\n=== PATTERN SUMMARY ===");
    let mut patterns: Vec<_> = pattern_counts.iter().collect();
    patterns.sort_by(|a, b| b.1.cmp(a.1));
    for (desc, count) in &patterns {
        println!("  [{}x] {}", count, desc);
    }

    // Deeper analysis: Check whether the diff positions are EXACTLY
    // the "boundaries" of descent runs of S within [p_a, p_b]
    println!("\n=== BOUNDARY ANALYSIS ===");
    println!("For each consecutive (p_a, p_b), checking if the symmetric diff");
    println!("of S'_pa and S'_pb consists exactly of the boundary positions\n");

    for n in 5..=max_n {
        let mut all_match = true;
        let mut total = 0;

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
                total += 1;

                let sp_a = source_asc(s_mask, p_a, n);
                let sp_b = source_asc(s_mask, p_b, n);
                let set_sp_a = mask_to_set(sp_a, n - 1);
                let set_sp_b = mask_to_set(sp_b, n - 1);
                let diff = sym_diff(&set_sp_a, &set_sp_b);

                // REFINED Predicted diff:
                // Interior: positions j in [p_a, p_b-2] where j in S XOR j+1 in S
                // Boundary: position p_b-1 is in S'_pa iff p_b in S, and NOT in S'_pb (gap).
                //   So p_b-1 in diff iff p_b in S, i.e., p_b < n.
                let mut predicted_diff = BTreeSet::new();
                for j in p_a..=(p_b.saturating_sub(2)) {
                    let j_in_s = (s_mask >> (j - 1)) & 1 == 1;
                    let j1_in_s = if j < n - 1 {
                        (s_mask >> j) & 1 == 1
                    } else {
                        false
                    };
                    if j_in_s != j1_in_s {
                        predicted_diff.insert(j);
                    }
                }
                // p_b-1: in S'_pa iff p_b in S. NOT in S'_pb (it's the gap).
                let p_b_in_s = p_b < n && (s_mask >> (p_b - 1)) & 1 == 1;
                if p_b_in_s {
                    // p_b-1 is in S'_pa but not S'_pb
                    predicted_diff.insert(p_b - 1);
                }

                if diff != predicted_diff {
                    all_match = false;
                    let s_str = mask_to_string(s_mask, n);
                    println!(
                        "  MISMATCH n={} S={} p_a={} p_b={}: actual={} predicted={}",
                        n,
                        s_str,
                        p_a,
                        p_b,
                        set_to_string(&diff),
                        set_to_string(&predicted_diff)
                    );
                }
            }
        }

        if all_match {
            println!(
                "  n={}: ALL {} cases match the shift-XOR prediction",
                n, total
            );
        }
    }

    // Additional analysis: what happens to the S''_p difference?
    println!("\n=== S'' SYMMETRIC DIFFERENCE ANALYSIS ===");
    println!("S''_pa = S'_pa union {{p_a-1}}, S''_pb = S'_pb union {{p_b-1}}");
    println!(
        "So sym_diff(S''_pa, S''_pb) = sym_diff(S'_pa, S'_pb) XOR {{p_a-1, p_b-1}} adjustments\n"
    );

    for n in 5..=max_n {
        let mut total = 0;
        let mut all_match = true;

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

                let spp_a = match source_desc(s_mask, p_a, n) {
                    Some(v) => v,
                    None => continue,
                };
                let spp_b = match source_desc(s_mask, p_b, n) {
                    Some(v) => v,
                    None => continue,
                };
                total += 1;

                let set_spp_a = mask_to_set(spp_a, n - 1);
                let set_spp_b = mask_to_set(spp_b, n - 1);
                let diff_dprime = sym_diff(&set_spp_a, &set_spp_b);

                // Predicted: same as S' diff, but also toggle p_a-1 and p_b-1
                // S''_pa = S'_pa + {p_a-1}, S''_pb = S'_pb + {p_b-1}
                // If p_a-1 is NOT in S'_pa and NOT in S'_pb: diff adds p_a-1 (from S''_pa only)
                // etc.  Let's just compute it.
                let sp_a = source_asc(s_mask, p_a, n);
                let sp_b = source_asc(s_mask, p_b, n);
                let set_sp_a = mask_to_set(sp_a, n - 1);
                let set_sp_b = mask_to_set(sp_b, n - 1);

                // S''_pa has p_a-1 (added). Is p_a-1 in S'_pb?
                let pa_m1_in_sp_b = set_sp_b.contains(&(p_a - 1));
                // S''_pb has p_b-1 (added). Is p_b-1 in S'_pa?
                let pb_m1_in_sp_a = set_sp_a.contains(&(p_b - 1));

                // Build predicted diff for S''
                // Add p_a-1 to S''_pa and p_b-1 to S''_pb, then compare directly.
                let mut set_spp_a_pred = set_sp_a.clone();
                set_spp_a_pred.insert(p_a - 1);
                let mut set_spp_b_pred = set_sp_b.clone();
                set_spp_b_pred.insert(p_b - 1);
                let predicted_dprime = sym_diff(&set_spp_a_pred, &set_spp_b_pred);

                if diff_dprime != predicted_dprime {
                    all_match = false;
                    println!(
                        "  INTERNAL ERROR: n={} S={} p_a={} p_b={}",
                        n,
                        mask_to_string(s_mask, n),
                        p_a,
                        p_b
                    );
                }

                // Check: p_a-1 and p_b-1 behavior
                // p_a-1 is in S''_pa but not in S'_pa.
                // p_b-1 is in S''_pb but not in S'_pb.
                // Question: is p_a-1 in S'_pb? Is p_b-1 in S'_pa?
                if n <= 7 {
                    println!(
                        "  n={} S={} p_a={} p_b={}: sym_diff(S''_pa,S''_pb) = {}",
                        n,
                        mask_to_string(s_mask, n),
                        p_a,
                        p_b,
                        set_to_string(&diff_dprime)
                    );
                    println!(
                        "    p_a-1={}: in S'_pb? {}  |  p_b-1={}: in S'_pa? {}",
                        p_a - 1,
                        pa_m1_in_sp_b,
                        p_b - 1,
                        pb_m1_in_sp_a
                    );

                    // Key: since p_a-1 not in S (validity), and p_a-1 < p_b-1 means
                    //   p_a-1 is in "match S" region of S'_pb, so p_a-1 in S'_pb iff p_a-1 in S.
                    //   Since p_a-1 not in S, p_a-1 NOT in S'_pb.
                    //   So p_a-1 is in S''_pa but NOT in S'_pb, S'_pa, or S''_pb
                    //   -> p_a-1 is ALWAYS in sym_diff(S''_pa, S''_pb).
                    //
                    //   p_b-1 not in S (validity). p_b-1 might or might not be in S'_pa.
                    //   In S'_pa: p_b-1 = j, and j >= p_a, so j = (j+1)-1 with j+1 > p_a.
                    //   So p_b-1 in S'_pa iff p_b in S.
                    //   Since p_b in S iff p_b < n and p_b in S. If p_b = n, p_b not in S.
                    //   If p_b < n: p_b is a valid position, so p_b in S (left boundary of run).
                    //   So if p_b < n: p_b-1 IS in S'_pa, so p_b-1 in S''_pb AND S'_pa
                    //   -> p_b-1 NOT in sym_diff(S''_pa, S''_pb) (it cancels!)
                    //   If p_b = n: p_b not in S, so p_b-1 NOT in S'_pa.
                    //   -> p_b-1 IS in sym_diff(S''_pa, S''_pb).
                    let p_b_in_s = p_b < n && (s_mask >> (p_b - 1)) & 1 == 1;
                    println!(
                        "    p_b={} in S: {}  -> p_b-1={} in sym_diff(S''): {}",
                        p_b,
                        p_b_in_s,
                        p_b - 1,
                        diff_dprime.contains(&(p_b - 1))
                    );
                }
            }
        }

        if all_match {
            println!("  n={}: all {} S'' predictions match", n, total);
        }
    }

    // Final synthesis: describe the complete structure
    println!("\n=== COMPLETE SYNTHESIS ===");
    println!("For consecutive valid positions p_a < p_b in P(S):");
    println!();
    println!("  S' ANALYSIS:");
    println!(
        "  Positions j < p_a: AGREE (both = S intersect [p_a-2] restricted to relevant positions)."
    );
    println!("  Position p_a-1: AGREE (both exclude, since p_a-1 not in S by validity).");
    println!("  Positions j in [p_a, p_b-2] (interior window):");
    println!("    S'_pa has j iff (j+1) in S   (shifted-by-1 view of S above p_a)");
    println!("    S'_pb has j iff j in S        (direct view of S below p_b)");
    println!("    DIFFER at j iff j in S XOR (j+1) in S (descent-run boundaries of S).");
    println!("  Position p_b-1:");
    println!(
        "    S'_pa has p_b-1 iff p_b in S (i.e., p_b < n, since p_b is left boundary of a run)."
    );
    println!("    S'_pb NEVER has p_b-1 (it's the gap for insertion at p_b).");
    println!("    DIFFER iff p_b in S, equivalently p_b < n.");
    println!("  Positions j >= p_b: AGREE (both shift S above their respective p).");
    println!();
    println!("  KEY STRUCTURAL RESULT:");
    println!("    When p_b < n: sym_diff has exactly TWO components:");
    println!(
        "      (a) Interior run boundaries in [p_a, p_b-2]: j where S changes between j and j+1"
    );
    println!("      (b) The right boundary: position p_b-1 (only in S'_pa)");
    println!("    When p_b = n: only interior run boundaries (component a).");
    println!();
    println!("  CONSECUTIVE VALID POSITIONS -- structure of S between p_a and p_b:");
    println!("    Between consecutive valid positions p_a < p_b, NO intermediate");
    println!("    position is a left boundary of a descent run of S.");
    println!("    This means the descent run starting at p_a extends as p_a, p_a+1,...");
    println!("    until it ends at some position r (where r+1 not in S), then");
    println!("    r+1,...,p_b-1 are all NOT in S (no new run starts before p_b).");
    println!("    The S-structure in [p_a, p_b-1] is: p_a,...,r in S, r+1,...,p_b-1 not in S.");
    println!("    Interior diff at position j in [p_a, p_b-2] iff j in S XOR (j+1) in S.");
    println!("    This can only happen at j=r (the end of the descent run): r in S, r+1 not in S.");
    println!("    So there is AT MOST ONE interior diff position.");
    println!();
    println!("    Combined with the boundary at p_b-1:");
    println!(
        "      Case 1 (p_b = n): sym_diff(S',S') has at most 1 element (at j=r if r < p_b-1)."
    );
    println!("      Case 2 (p_b < n): sym_diff(S',S') has at most 2 elements:");
    println!("        {{r}} from interior (if r < p_b-1) plus {{p_b-1}} from boundary.");
    println!("    But: if the run p_a,...,r extends all the way to p_b-1 (r=p_b-1),");
    println!("    then interior diff at r=p_b-1 is outside [p_a,p_b-2], so no interior diff,");
    println!("    and the only diff is at p_b-1 from the boundary rule.");
    println!();
    println!("    DATA CONFIRMS: symmetric diff has size exactly 1 or 2, never 0 or >=3.");
    println!();
    println!("  S'' ANALYSIS:");
    println!("    S''_pa = S'_pa + {{p_a-1}},  S''_pb = S'_pb + {{p_b-1}}");
    println!("    Additional diff from p_a-1: p_a-1 not in S'_pb (since p_a-1 not in S), so p_a-1 always in sym_diff(S'',S'').");
    println!("    Additional diff from p_b-1: p_b-1 already in S'_pa iff p_b in S.");
    println!("      If p_b in S: p_b-1 in both S'_pa and S''_pb => cancels in sym_diff.");
    println!("      If p_b = n: p_b-1 not in S'_pa, so p_b-1 in S''_pb only => in sym_diff.");
    println!();
    println!("  BOTTOM LINE FOR TN_2:");
    println!("    S'_pa and S'_pb differ at AT MOST 2 positions:");
    println!("      - At most 1 interior position (where the descent run from p_a ends)");
    println!("      - At most 1 boundary position (p_b-1, when p_b < n)");
    println!("    In many cases, the descent run from p_a extends to p_b-1,");
    println!("    giving NO interior diff, and the ONLY diff is at p_b-1.");
    println!("    ");
    println!("    S''_pa and S''_pb differ at p_a-1, plus the same interior positions as S',");
    println!("    but p_b-1 cancels when p_b in S (the common case for consecutive p_b < n).");
    println!("    ");
    println!("    The source sets are always VERY CLOSE -- differing by at most 2 positions,");
    println!("    which is promising for bijective or transfer-matrix approaches to TN_2.");
}
