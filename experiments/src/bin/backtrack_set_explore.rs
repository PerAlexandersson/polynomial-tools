/// Explore set-valued statistics under the backtrack map.
///
/// For each (pattern, set-statistic) pair, checks whether
///   SetStat(backtrack_π(Av_pat(n))) = SetStat(π^{-1})   or   SetStat(π)
/// holds for all π ∈ S_n, and reports:
///   - which set-stat identities hold (and up to which n),
///   - equidistribution: whether the multiset of set-values over S_n matches,
///   - the number of distinct set-values attained by the backtrack image,
///   - cross-pattern coincidences (same multiset of sets for different patterns).
///
/// Also investigates finer structure:
///   - descent-bottom and descent-top preservation,
///   - joint preservation of pairs of set-statistics,
///   - the fiber structure: for each descent set D, how many π map to D?
use combpoly::permutation::{
    all_permutations, avoiding_permutations, backtrack_image, inverse, parse_sequence,
};
use combpoly::statistics::{compute_set, SetStat};
use std::collections::{BTreeMap, BTreeSet};
use std::env;

type Set = BTreeSet<usize>;

/// Complement: S -> [n] \ S  (positions not in S)
fn set_complement(s: &Set, n: usize) -> Set {
    (1..=n).filter(|x| !s.contains(x)).collect()
}

/// Flip: S -> {n+1-x : x in S}
fn set_flip(s: &Set, n: usize) -> Set {
    s.iter().map(|&x| n + 1 - x).collect()
}

/// Both: complement then flip
fn set_comp_flip(s: &Set, n: usize) -> Set {
    set_flip(&set_complement(s, n), n)
}

fn main() {
    let max_n: u8 = env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    let patterns = ["123", "132", "213", "231", "312", "321"];
    let set_stats: &[(&str, SetStat)] = &[
        ("DesSet", SetStat::DesSet),
        ("AscSet", SetStat::AscSet),
        ("DesBottomSet", SetStat::DesBottomSet),
        ("DesTopSet", SetStat::DesTopSet),
        ("ExcSet", SetStat::ExcSet),
        ("FixSet", SetStat::FixSet),
        ("PeakSet", SetStat::PeakSet),
        ("ValleySet", SetStat::ValleySet),
        ("LrminSet", SetStat::LrminSet),
        ("LrmaxSet", SetStat::LrmaxSet),
        ("RlminSet", SetStat::RlminSet),
        ("RlmaxSet", SetStat::RlmaxSet),
    ];

    println!("========================================");
    println!("  Set-valued backtrack exploration");
    println!("  max_n = {}", max_n);
    println!("========================================\n");

    // Part 1: Check set-stat identities
    println!("=== PART 1: Identity checks ===");
    println!("Testing: SetStat(backtrack_π(Av_pat)) = SetStat(π^{{-1}}) or SetStat(π)\n");

    for pat_str in &patterns {
        let pat = parse_sequence(pat_str);
        println!("--- Av_{} ---", pat_str);

        for &(stat_name, stat) in set_stats {
            let mut matches_inverse = true;
            let mut matches_pi = true;
            let mut fail_inv_n = 0u8;
            let mut fail_pi_n = 0u8;

            for n in 1..=max_n {
                let all_pi = all_permutations(n);
                let s_set = avoiding_permutations(n, &[pat.clone()]);

                for pi in &all_pi {
                    let sigma = backtrack_image(pi, &s_set);
                    let bt_set = compute_set(&sigma, stat);
                    if matches_inverse {
                        let inv = inverse(pi);
                        let inv_set = compute_set(&inv, stat);
                        if bt_set != inv_set {
                            matches_inverse = false;
                            fail_inv_n = n;
                        }
                    }
                    if matches_pi {
                        let pi_set = compute_set(pi, stat);
                        if bt_set != pi_set {
                            matches_pi = false;
                            fail_pi_n = n;
                        }
                    }
                }
                if !matches_inverse && !matches_pi {
                    break;
                }
            }

            let mut notes = Vec::new();
            if matches_inverse {
                notes.push(format!("= inv(π) ✓ (n≤{})", max_n));
            } else if fail_inv_n > 2 {
                notes.push(format!("= inv(π) fails at n={}", fail_inv_n));
            }
            if matches_pi {
                notes.push(format!("= π ✓ (n≤{})", max_n));
            } else if fail_pi_n > 2 {
                notes.push(format!("= π fails at n={}", fail_pi_n));
            }
            if !notes.is_empty() {
                println!("  {:14} {}", stat_name, notes.join(", "));
            }
        }
        println!();
    }

    // Part 2: Equidistribution of set-values
    println!("\n=== PART 2: Equidistribution ===");
    println!("Comparing multisets of set-values: backtrack vs π vs π^{{-1}}\n");

    for n in [4, 5, 6, 7].iter().copied().filter(|&n| n <= max_n) {
        println!("--- n = {} ---", n);
        let all_pi = all_permutations(n);

        for pat_str in &patterns {
            let pat = parse_sequence(pat_str);
            let s_set = avoiding_permutations(n, &[pat.clone()]);

            for &(stat_name, stat) in set_stats {
                // Collect multisets: map from set-value to count
                let mut bt_dist: BTreeMap<Set, usize> = BTreeMap::new();
                let mut pi_dist: BTreeMap<Set, usize> = BTreeMap::new();
                let mut inv_dist: BTreeMap<Set, usize> = BTreeMap::new();

                for pi in &all_pi {
                    let sigma = backtrack_image(pi, &s_set);
                    *bt_dist.entry(compute_set(&sigma, stat)).or_default() += 1;
                    *pi_dist.entry(compute_set(pi, stat)).or_default() += 1;
                    let inv = inverse(pi);
                    *inv_dist.entry(compute_set(&inv, stat)).or_default() += 1;
                }

                let eq_pi = bt_dist == pi_dist;
                let eq_inv = bt_dist == inv_dist;
                if eq_pi || eq_inv {
                    let tag = match (eq_pi, eq_inv) {
                        (true, true) => "≡ π ≡ inv(π)",
                        (true, false) => "≡ π",
                        (false, true) => "≡ inv(π)",
                        _ => unreachable!(),
                    };
                    println!(
                        "  Av_{} + {:14} {} ({} distinct sets)",
                        pat_str,
                        stat_name,
                        tag,
                        bt_dist.len()
                    );
                }
            }
        }
        println!();
    }

    // Part 3: Cross-pattern coincidences
    println!("\n=== PART 3: Cross-pattern coincidences ===");
    println!("Pairs of patterns giving identical set-value multisets\n");

    let n_cross: u8 = max_n.min(6);
    for &(stat_name, stat) in set_stats {
        let all_pi = all_permutations(n_cross);
        let mut dists: Vec<(&str, BTreeMap<Set, usize>)> = Vec::new();

        for pat_str in &patterns {
            let pat = parse_sequence(pat_str);
            let s_set = avoiding_permutations(n_cross, &[pat.clone()]);

            let mut dist: BTreeMap<Set, usize> = BTreeMap::new();
            for pi in &all_pi {
                let sigma = backtrack_image(pi, &s_set);
                *dist.entry(compute_set(&sigma, stat)).or_default() += 1;
            }
            dists.push((pat_str, dist));
        }

        let mut coincidences: Vec<String> = Vec::new();
        for i in 0..dists.len() {
            for j in i + 1..dists.len() {
                if dists[i].1 == dists[j].1 {
                    coincidences.push(format!("Av_{} = Av_{}", dists[i].0, dists[j].0));
                }
            }
        }
        if !coincidences.is_empty() {
            println!(
                "  {:14} (n={}): {}",
                stat_name,
                n_cross,
                coincidences.join(", ")
            );
        }
    }

    // Part 4: Joint preservation of pairs
    println!("\n\n=== PART 4: Joint preservation of set-stat pairs ===");
    println!(
        "Testing: (A(bt), B(bt)) = (A(inv), B(inv)) for all π, n≤{}\n",
        max_n
    );

    let key_stats: &[(&str, SetStat)] = &[
        ("DesSet", SetStat::DesSet),
        ("DesBottomSet", SetStat::DesBottomSet),
        ("DesTopSet", SetStat::DesTopSet),
        ("PeakSet", SetStat::PeakSet),
        ("ExcSet", SetStat::ExcSet),
        ("ValleySet", SetStat::ValleySet),
    ];

    for pat_str in &["213", "312"] {
        let pat = parse_sequence(pat_str);
        println!("--- Av_{} ---", pat_str);

        for i in 0..key_stats.len() {
            for j in i + 1..key_stats.len() {
                let (name_a, stat_a) = key_stats[i];
                let (name_b, stat_b) = key_stats[j];
                let mut holds = true;
                let mut fail_n = 0u8;

                'outer: for n in 1..=max_n {
                    let all_pi = all_permutations(n);
                    let s_set = avoiding_permutations(n, &[pat.clone()]);

                    for pi in &all_pi {
                        let sigma = backtrack_image(pi, &s_set);
                        let inv = inverse(pi);
                        let bt_a = compute_set(&sigma, stat_a);
                        let bt_b = compute_set(&sigma, stat_b);
                        let inv_a = compute_set(&inv, stat_a);
                        let inv_b = compute_set(&inv, stat_b);
                        if bt_a != inv_a || bt_b != inv_b {
                            holds = false;
                            fail_n = n;
                            break 'outer;
                        }
                    }
                }

                if holds {
                    println!(
                        "  ({}, {}) = on inv(π)  ✓ (n≤{})",
                        name_a, name_b, max_n
                    );
                } else if fail_n > 3 {
                    println!(
                        "  ({}, {}) = on inv(π)  fails at n={}",
                        name_a, name_b, fail_n
                    );
                }
            }
        }
        println!();
    }

    // Part 5: Fiber structure — for Av_213 + DesSet
    println!("\n=== PART 5: Fiber structure for Av_213 + DesSet ===");
    println!("For each descent set D, how many π give backtrack with Des = D?\n");

    for n in (3..=max_n.min(6)).rev().take(1) {
        let pat = parse_sequence("213");
        let all_pi = all_permutations(n);
        let s_set = avoiding_permutations(n, &[pat]);

        let mut fibers: BTreeMap<Set, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &all_pi {
            let sigma = backtrack_image(pi, &s_set);
            let d = compute_set(&sigma, SetStat::DesSet);
            fibers.entry(d).or_default().push(pi.clone());
        }

        println!("n = {} ({} distinct descent sets):", n, fibers.len());
        for (d, pis) in &fibers {
            let d_vec: Vec<usize> = d.iter().copied().collect();
            // Also show the descent set of inv(π) for the first few fibers
            let sample_inv_des: Vec<Set> = pis
                .iter()
                .take(3)
                .map(|pi| compute_set(&inverse(pi), SetStat::DesSet))
                .collect();
            println!(
                "  Des = {:?}  |fiber| = {}  (sample inv Des: {:?})",
                d_vec,
                pis.len(),
                sample_inv_des
            );
        }
    }

    // Part 6: Descent-bottom exploration across all patterns
    println!("\n\n=== PART 6: Descent-bottom/top preservation ===");
    println!("For each pattern: does DesBottomSet(bt) determine or equal DesBottomSet(inv(π))?\n");

    for pat_str in &patterns {
        let pat = parse_sequence(pat_str);
        let mut bottom_holds = true;
        let mut top_holds = true;
        let mut bottom_fail = 0u8;
        let mut top_fail = 0u8;

        for n in 1..=max_n {
            let all_pi = all_permutations(n);
            let s_set = avoiding_permutations(n, &[pat.clone()]);

            for pi in &all_pi {
                let sigma = backtrack_image(pi, &s_set);
                let inv = inverse(pi);
                if bottom_holds
                    && compute_set(&sigma, SetStat::DesBottomSet)
                        != compute_set(&inv, SetStat::DesBottomSet)
                {
                    bottom_holds = false;
                    bottom_fail = n;
                }
                if top_holds
                    && compute_set(&sigma, SetStat::DesTopSet)
                        != compute_set(&inv, SetStat::DesTopSet)
                {
                    top_holds = false;
                    top_fail = n;
                }
            }
            if !bottom_holds && !top_holds {
                break;
            }
        }

        let mut notes = Vec::new();
        if bottom_holds {
            notes.push(format!("DesBottom = inv(π) ✓ (n≤{})", max_n));
        } else {
            notes.push(format!("DesBottom fails n={}", bottom_fail));
        }
        if top_holds {
            notes.push(format!("DesTop = inv(π) ✓ (n≤{})", max_n));
        } else {
            notes.push(format!("DesTop fails n={}", top_fail));
        }
        println!("  Av_{}: {}", pat_str, notes.join(", "));
    }

    // Part 7: What about SetStat(backtrack) = SetStat(foata(π)) or other maps?
    println!("\n\n=== PART 7: Alternative target maps ===");
    println!("Testing: SetStat(bt) = SetStat(foata(π)) or SetStat(reverse(π))\n");

    for pat_str in &patterns {
        let pat = parse_sequence(pat_str);
        let mut notes: Vec<String> = Vec::new();

        for &(stat_name, stat) in &[
            ("DesSet", SetStat::DesSet),
            ("PeakSet", SetStat::PeakSet),
            ("DesBottomSet", SetStat::DesBottomSet),
        ] {
            let mut foata_ok = true;
            let mut rev_ok = true;
            let mut comp_ok = true;

            for n in 1..=max_n {
                let all_pi = all_permutations(n);
                let s_set = avoiding_permutations(n, &[pat.clone()]);

                for pi in &all_pi {
                    let sigma = backtrack_image(pi, &s_set);
                    let bt_set = compute_set(&sigma, stat);

                    if foata_ok {
                        let foata = combpoly::permutation::foata_map(pi);
                        if compute_set(&foata, stat) != bt_set {
                            foata_ok = false;
                        }
                    }
                    if rev_ok {
                        let rev: Vec<u8> = pi.iter().rev().copied().collect();
                        if compute_set(&rev, stat) != bt_set {
                            rev_ok = false;
                        }
                    }
                    if comp_ok {
                        let n_val = n as u8;
                        let comp: Vec<u8> = pi.iter().map(|&x| n_val + 1 - x).collect();
                        if compute_set(&comp, stat) != bt_set {
                            comp_ok = false;
                        }
                    }
                }
                if !foata_ok && !rev_ok && !comp_ok {
                    break;
                }
            }

            if foata_ok {
                notes.push(format!("{}(bt) = {}(foata(π)) ✓", stat_name, stat_name));
            }
            if rev_ok {
                notes.push(format!("{}(bt) = {}(rev(π)) ✓", stat_name, stat_name));
            }
            if comp_ok {
                notes.push(format!("{}(bt) = {}(comp(π)) ✓", stat_name, stat_name));
            }
        }

        if !notes.is_empty() {
            println!("  Av_{}:", pat_str);
            for note in &notes {
                println!("    {}", note);
            }
        }
    }

    // Part 8: Twisted identities using set operations
    println!("\n\n=== PART 8: Twisted identities via set operations ===");
    println!("Testing: SetStat(bt) = op(SetStat(target))  where op ∈ {{complement, flip, comp∘flip}}");
    println!("         target ∈ {{π, inv(π)}}\n");

    let transforms: &[(&str, fn(&Set, usize) -> Set)] = &[
        ("complement", set_complement),
        ("flip", set_flip),
        ("comp_flip", set_comp_flip),
    ];
    let targets: &[&str] = &["pi", "inv"];

    for pat_str in &patterns {
        let pat = parse_sequence(pat_str);
        let mut found_any = false;

        for &(stat_name, stat) in set_stats {
            for &(op_name, op_fn) in transforms {
                for &target_name in targets {
                    let mut holds = true;
                    let mut fail_n = 0u8;

                    'check: for n in 1..=max_n {
                        let all_pi = all_permutations(n);
                        let s_set = avoiding_permutations(n, &[pat.clone()]);
                        let nn = n as usize;

                        for pi in &all_pi {
                            let sigma = backtrack_image(pi, &s_set);
                            let bt_set = compute_set(&sigma, stat);

                            let target_perm = match target_name {
                                "pi" => pi.clone(),
                                "inv" => inverse(pi),
                                _ => unreachable!(),
                            };
                            let target_set = compute_set(&target_perm, stat);
                            let transformed = op_fn(&target_set, nn);

                            if bt_set != transformed {
                                holds = false;
                                fail_n = n;
                                break 'check;
                            }
                        }
                    }

                    if holds {
                        if !found_any {
                            println!("  Av_{}:", pat_str);
                            found_any = true;
                        }
                        println!(
                            "    {}(bt) = {}({}({}))  ✓ (n≤{})",
                            stat_name, op_name, stat_name, target_name, max_n
                        );
                    }
                }
            }
        }
        if !found_any {
            println!("  Av_{}: (no twisted identities found)", pat_str);
        }
        println!();
    }

    // Part 9: Equidistribution with set transforms
    println!("\n=== PART 9: Twisted equidistribution ===");
    println!("Cases where multiset of SetStat(bt) ≡ multiset of op(SetStat(target))\n");

    let n_eq: u8 = max_n.min(6);
    let all_pi = all_permutations(n_eq);

    for pat_str in &patterns {
        let pat = parse_sequence(pat_str);
        let s_set = avoiding_permutations(n_eq, &[pat.clone()]);
        let nn = n_eq as usize;
        let mut found_any = false;

        for &(stat_name, stat) in set_stats {
            // Compute backtrack distribution
            let mut bt_dist: BTreeMap<Set, usize> = BTreeMap::new();
            for pi in &all_pi {
                let sigma = backtrack_image(pi, &s_set);
                *bt_dist.entry(compute_set(&sigma, stat)).or_default() += 1;
            }

            // Already equidist with plain inv? Skip (reported in Part 2)
            let mut inv_dist: BTreeMap<Set, usize> = BTreeMap::new();
            for pi in &all_pi {
                let inv = inverse(pi);
                *inv_dist.entry(compute_set(&inv, stat)).or_default() += 1;
            }
            if bt_dist == inv_dist {
                continue;
            }

            // Already equidist with plain π? Skip
            let mut pi_dist: BTreeMap<Set, usize> = BTreeMap::new();
            for pi in &all_pi {
                *pi_dist.entry(compute_set(pi, stat)).or_default() += 1;
            }
            if bt_dist == pi_dist {
                continue;
            }

            // Try transforms
            for &(op_name, op_fn) in transforms {
                for &target_name in targets {
                    let mut transformed_dist: BTreeMap<Set, usize> = BTreeMap::new();
                    for pi in &all_pi {
                        let target_perm = match target_name {
                            "pi" => pi.clone(),
                            "inv" => inverse(pi),
                            _ => unreachable!(),
                        };
                        let target_set = compute_set(&target_perm, stat);
                        let transformed = op_fn(&target_set, nn);
                        *transformed_dist.entry(transformed).or_default() += 1;
                    }

                    if bt_dist == transformed_dist {
                        if !found_any {
                            println!("  Av_{} (n={}):", pat_str, n_eq);
                            found_any = true;
                        }
                        println!(
                            "    {}(bt) ≡ {}({}({}))  ({} distinct)",
                            stat_name,
                            op_name,
                            stat_name,
                            target_name,
                            bt_dist.len()
                        );
                    }
                }
            }
        }
        if found_any {
            println!();
        }
    }

    // Part 10: Cross-statistic identities with transforms
    println!("\n=== PART 10: Cross-statistic identities ===");
    println!("Testing: StatA(bt) = op(StatB(target))  for A ≠ B, target ∈ {{π, inv(π)}}\n");

    let cross_stats: &[(&str, SetStat)] = &[
        ("DesSet", SetStat::DesSet),
        ("AscSet", SetStat::AscSet),
        ("DesBottomSet", SetStat::DesBottomSet),
        ("DesTopSet", SetStat::DesTopSet),
        ("ExcSet", SetStat::ExcSet),
        ("FixSet", SetStat::FixSet),
        ("PeakSet", SetStat::PeakSet),
        ("ValleySet", SetStat::ValleySet),
        ("LrminSet", SetStat::LrminSet),
        ("LrmaxSet", SetStat::LrmaxSet),
        ("RlminSet", SetStat::RlminSet),
        ("RlmaxSet", SetStat::RlmaxSet),
    ];

    let all_transforms: &[(&str, Option<fn(&Set, usize) -> Set>)] = &[
        ("id", None),
        ("complement", Some(set_complement as fn(&Set, usize) -> Set)),
        ("flip", Some(set_flip as fn(&Set, usize) -> Set)),
        ("comp_flip", Some(set_comp_flip as fn(&Set, usize) -> Set)),
    ];

    for pat_str in &patterns {
        let pat = parse_sequence(pat_str);
        let mut found = Vec::new();

        for &(name_a, stat_a) in cross_stats {
            for &(name_b, stat_b) in cross_stats {
                if name_a == name_b {
                    continue;
                }
                for &(op_name, op_fn) in all_transforms {
                    for &target_name in targets {
                        let mut holds = true;

                        'outer2: for n in 1..=max_n {
                            let all_pi_n = all_permutations(n);
                            let s_set = avoiding_permutations(n, &[pat.clone()]);
                            let nn = n as usize;

                            for pi in &all_pi_n {
                                let sigma = backtrack_image(pi, &s_set);
                                let bt_a = compute_set(&sigma, stat_a);
                                let target_perm = match target_name {
                                    "pi" => pi.clone(),
                                    "inv" => inverse(pi),
                                    _ => unreachable!(),
                                };
                                let target_set = compute_set(&target_perm, stat_b);
                                let transformed = match op_fn {
                                    Some(f) => f(&target_set, nn),
                                    None => target_set,
                                };
                                if bt_a != transformed {
                                    holds = false;
                                    break 'outer2;
                                }
                            }
                        }

                        if holds {
                            found.push(format!(
                                "{}(bt) = {}({}({}))  ✓",
                                name_a, op_name, name_b, target_name
                            ));
                        }
                    }
                }
            }
        }

        if !found.is_empty() {
            println!("  Av_{}:", pat_str);
            for f in &found {
                println!("    {}", f);
            }
            println!();
        }
    }

    println!("\n========================================");
    println!("  Done.");
    println!("========================================");
}
