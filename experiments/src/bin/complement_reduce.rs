/// Complement reduction: map descent sets S with 1 in S to complement descent
/// sets S_bar = {n-i : i in S} and check if S_bar has 1 not in S_bar.
///
/// The complement map on S_n sends sigma to sigma_bar where
/// sigma_bar(i) = n+1 - sigma(n+1-i). This sends Des(sigma) = S to
/// Des(sigma_bar) = {n-i : i in S}, and swaps is invariant under complement.
///
/// So L_{n,S}(t) = L_{n,S_bar}(t).
///
/// If S_bar has 1 not in S_bar, then the shift lemma / compatible family approach
/// works for S_bar, giving us real-rootedness of L_{n,S_bar} = L_{n,S} for free.
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{
    check_weak_interlacing, format_poly, is_real_rooted,
};
use std::collections::BTreeMap;

fn build_poly(perms: &[&Vec<u8>]) -> Vec<i64> {
    if perms.is_empty() {
        return vec![0];
    }
    let max_s = perms.iter().map(|s| compute(s, Stat::Swaps)).max().unwrap_or(0);
    let mut coeffs = vec![0i64; max_s + 1];
    for s in perms {
        coeffs[compute(s, Stat::Swaps)] += 1;
    }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
        coeffs.pop();
    }
    coeffs
}

fn weakly_interlace(f: &[i64], g: &[i64]) -> bool {
    if f.len() <= 1 || g.len() <= 1 {
        return true;
    }
    let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
    check_weak_interlacing(small, large) == Some(true)
}

fn descent_set_str(s: u64, n: u8) -> String {
    let positions: Vec<String> = (0..n - 1)
        .filter(|&i| s & (1 << i) != 0)
        .map(|i| (i + 1).to_string())
        .collect();
    if positions.is_empty() { "{}".to_string() }
    else { format!("{{{}}}", positions.join(",")) }
}

/// Compute complement of descent set: S_bar = {n-i : i in S} where S is a subset of [n-1].
/// Input and output are bitmasks where bit j (0-indexed) corresponds to position j+1.
fn complement_descent_set(ds: u64, n: u8) -> u64 {
    let mut result = 0u64;
    for i in 0..(n - 1) {
        if ds & (1 << i) != 0 {
            // Position i+1 is in S, so n-(i+1) = n-i-1 is in S_bar
            // That is position (n-i-1) which is bit (n-i-2) in the bitmask
            let new_pos = (n - 1) - (i + 1); // 0-indexed bit for position (n-i-1)
            result |= 1 << new_pos;
        }
    }
    result
}

/// Apply complement map to a permutation: sigma_bar(i) = n+1 - sigma(n+1-i).
/// Input is 1-indexed (values 1..n), output same.
fn complement_perm(sigma: &[u8]) -> Vec<u8> {
    let n = sigma.len() as u8;
    let mut result = vec![0u8; n as usize];
    for i in 0..n as usize {
        // sigma_bar(i+1) = n+1 - sigma(n+1-(i+1)) = n+1 - sigma(n-i)
        // 0-indexed: result[i] = n+1 - sigma[n-1-i]
        result[i] = n + 1 - sigma[n as usize - 1 - i];
    }
    result
}

/// Check if position-refined pieces form a compatible (pairwise interlacing) family.
fn check_pos_compatible(perms: &[&Vec<u8>], n: u8) -> bool {
    let mut by_pos: BTreeMap<usize, Vec<&Vec<u8>>> = BTreeMap::new();
    for s in perms {
        let pos = s.iter().position(|&v| v == n).unwrap() + 1;
        by_pos.entry(pos).or_default().push(s);
    }

    let pos_polys: Vec<Vec<i64>> = by_pos
        .values()
        .map(|grp| build_poly(grp))
        .collect();

    for i in 0..pos_polys.len() {
        for j in i + 1..pos_polys.len() {
            if !weakly_interlace(&pos_polys[i], &pos_polys[j]) {
                return false;
            }
        }
    }
    true
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    println!("================================================================");
    println!("  Complement reduction: mapping '1 in S' cases to '1 not in S'");
    println!("================================================================\n");

    let mut grand_total_bad = 0u64;
    let mut grand_reduced = 0u64;
    let mut grand_not_reduced = 0u64;

    for n in 3..=max_n {
        let all = all_permutations(n);

        // Group by descent set
        let mut by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for s in &all {
            by_des.entry(descent_set_bitmask(s)).or_default().push(s);
        }

        // First, identify which descent sets are "bad" (1 in S AND pos-of-max incompatible)
        let mut bad_sets: Vec<u64> = Vec::new();
        let mut good_sets: Vec<u64> = Vec::new();

        for (&ds, perms) in &by_des {
            let compatible = check_pos_compatible(perms, n);
            if compatible {
                good_sets.push(ds);
            } else {
                bad_sets.push(ds);
            }
        }

        let total_1_in_s: Vec<u64> = by_des.keys()
            .filter(|&&ds| ds & 1 != 0)
            .copied()
            .collect();

        let bad_with_1: Vec<u64> = bad_sets.iter()
            .filter(|&&ds| ds & 1 != 0)
            .copied()
            .collect();

        println!("--- n = {} ---", n);
        println!("  Total descent sets: {}", by_des.len());
        println!("  Descent sets with 1 in S: {}", total_1_in_s.len());
        println!("  Pos-compatible (good): {}", good_sets.len());
        println!("  Pos-incompatible (bad): {}", bad_sets.len());
        println!("  Bad with 1 in S: {}", bad_with_1.len());

        // For each bad descent set, compute complement and check
        let mut n_reduced = 0u32;
        let mut n_not_reduced = 0u32;
        let mut not_reduced_details: Vec<String> = Vec::new();

        for &ds in &bad_sets {
            let ds_bar = complement_descent_set(ds, n);
            let bar_has_1 = ds_bar & 1 != 0;
            let bar_is_good = good_sets.contains(&ds_bar);

            // Verify: L_{n,S} = L_{n,S_bar}
            let poly_s = build_poly(&by_des[&ds]);
            let poly_s_bar = if let Some(perms_bar) = by_des.get(&ds_bar) {
                build_poly(perms_bar)
            } else {
                vec![0]
            };

            let polys_match = poly_s == poly_s_bar;

            if !bar_has_1 && bar_is_good {
                n_reduced += 1;
            } else if !bar_has_1 && polys_match {
                // S_bar has 1 not in it but is also bad for pos-of-max
                // Still, maybe it's handled by shift lemma on S_bar's side?
                n_reduced += 1; // 1 not in S_bar means shift lemma works
            } else if !bar_has_1 {
                n_reduced += 1; // 1 not in S_bar
            } else {
                n_not_reduced += 1;
                if not_reduced_details.len() < 5 {
                    not_reduced_details.push(format!(
                        "  S={} -> S_bar={} (1_in_S_bar:{}, bar_good:{}, match:{})",
                        descent_set_str(ds, n), descent_set_str(ds_bar, n),
                        bar_has_1, bar_is_good, polys_match,
                    ));
                }
            }
        }

        // Also verify the complement identity computationally for ALL descent sets
        let mut complement_verified = 0u32;
        let mut complement_failed = 0u32;
        for (&ds, perms) in &by_des {
            let ds_bar = complement_descent_set(ds, n);
            let poly_s = build_poly(perms);

            // Alternative: compute L_{n,S_bar} by complementing permutations
            let compl_perms: Vec<Vec<u8>> = perms.iter()
                .map(|s| complement_perm(s))
                .collect();

            // Verify each complemented perm has descent set ds_bar
            let mut all_correct_des = true;
            for cp in &compl_perms {
                if descent_set_bitmask(cp) != ds_bar {
                    all_correct_des = false;
                    break;
                }
            }

            // Verify swaps is invariant
            let mut swaps_match = true;
            for (orig, compl) in perms.iter().zip(compl_perms.iter()) {
                if compute(orig, Stat::Swaps) != compute(compl, Stat::Swaps) {
                    swaps_match = false;
                    break;
                }
            }

            if all_correct_des && swaps_match {
                // The complement map gives a bijection D(n,S) -> D(n,S_bar)
                // preserving swaps, so L_{n,S} = L_{n,S_bar}
                let poly_bar = if let Some(perms_bar) = by_des.get(&ds_bar) {
                    build_poly(perms_bar)
                } else {
                    vec![0]
                };
                if poly_s == poly_bar {
                    complement_verified += 1;
                } else {
                    complement_failed += 1;
                    println!("  MISMATCH: S={} poly={}, S_bar={} poly={}",
                        descent_set_str(ds, n), format_poly(&poly_s),
                        descent_set_str(ds_bar, n), format_poly(&poly_bar));
                }
            } else {
                complement_failed += 1;
                if !all_correct_des {
                    println!("  Des set mismatch for complement of S={}!", descent_set_str(ds, n));
                }
                if !swaps_match {
                    println!("  Swaps NOT invariant for complement of S={}!", descent_set_str(ds, n));
                }
            }
        }

        println!("  Complement identity verified: {}/{}", complement_verified, by_des.len());
        if complement_failed > 0 {
            println!("  WARNING: complement identity failed for {} sets!", complement_failed);
        }

        println!("  Bad sets reduced via complement (1 not in S_bar): {}/{}", n_reduced, bad_sets.len());
        if n_not_reduced > 0 {
            println!("  NOT reduced (both S and S_bar have 1): {}", n_not_reduced);
            for detail in &not_reduced_details {
                println!("{}", detail);
            }
        }
        println!();

        grand_total_bad += bad_sets.len() as u64;
        grand_reduced += n_reduced as u64;
        grand_not_reduced += n_not_reduced as u64;
    }

    println!("================================================================");
    println!("  Grand summary");
    println!("================================================================");
    println!("Total bad (pos-incompatible) descent sets: {}", grand_total_bad);
    println!("Reduced via complement (1 not in S_bar): {}", grand_reduced);
    println!("Not reducible (both S and S_bar have 1): {}", grand_not_reduced);
    if grand_not_reduced > 0 {
        println!("\n--- Characterizing irreducible cases ---");
        println!("These S satisfy: 1 in S AND n-1 in S (so both S and S_bar have 1).\n");

        // Re-run to characterize the remaining cases
        for n in 3..=max_n {
            let all = all_permutations(n);
            let mut by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
            for s in &all {
                by_des.entry(descent_set_bitmask(s)).or_default().push(s);
            }

            for (&ds, perms) in &by_des {
                if ds & 1 == 0 { continue; } // skip if 1 not in S
                let ds_bar = complement_descent_set(ds, n);
                if ds_bar & 1 == 0 { continue; } // reducible

                // This set is irreducible: 1 in S and 1 in S_bar
                // Equivalently: 1 in S and n-1 in S
                let compatible = check_pos_compatible(perms, n);
                if compatible { continue; } // actually good after all

                let poly = build_poly(perms);
                let rr = poly.len() <= 2 || is_real_rooted(&poly);

                // Characterize S
                let elements: Vec<usize> = (0..n-1)
                    .filter(|&i| ds & (1u64 << i) != 0)
                    .map(|i| (i + 1) as usize)
                    .collect();

                println!("  n={} S={} |S|={} poly={} rr={} self-complement={}",
                    n, descent_set_str(ds, n), elements.len(), format_poly(&poly), rr,
                    ds == ds_bar);
            }
        }
    }

    // Part 2: Detailed analysis of the complement map on interlacing
    println!("\n================================================================");
    println!("  Approach detail: Does the complement map resolve ALL bad sets?");
    println!("================================================================\n");

    for n in 3..=max_n {
        let all = all_permutations(n);
        let mut by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for s in &all {
            by_des.entry(descent_set_bitmask(s)).or_default().push(s);
        }

        let mut bad = 0u32;
        let mut resolved_by_complement = 0u32;
        let mut resolved_by_complement_being_good = 0u32;
        let mut unresolved = 0u32;

        for (&ds, perms) in &by_des {
            let compatible = check_pos_compatible(perms, n);
            if compatible { continue; }

            bad += 1;
            let ds_bar = complement_descent_set(ds, n);

            // Check if S_bar is good (pos-compatible)
            if let Some(perms_bar) = by_des.get(&ds_bar) {
                let bar_compat = check_pos_compatible(perms_bar, n);
                if bar_compat {
                    resolved_by_complement_being_good += 1;
                    resolved_by_complement += 1;
                    continue;
                }
            }

            // Even if S_bar is not pos-compatible, if 1 not in S_bar then the
            // shift lemma works for ALL positions of S_bar
            if ds_bar & 1 == 0 {
                resolved_by_complement += 1;
                continue;
            }

            unresolved += 1;
        }

        println!("n={}: {} bad, {} resolved by complement, {} unresolved",
            n, bad, resolved_by_complement, unresolved);
    }
}
