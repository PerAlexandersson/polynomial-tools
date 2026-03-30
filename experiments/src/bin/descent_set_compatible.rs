/// For each descent set S ⊆ [n-1], refine L_{n,S}(t) by position of max value n.
/// Check if the position-refined pieces form a compatible family
/// (pairwise weakly interlacing), which would imply real-rootedness of L_{n,S}.
///
/// Also test the insertion-based decomposition: for each valid insertion
/// position p, the contribution to L_{n,S} comes from L_{n-1,S'_p} with
/// known ε₂ and potentially nonzero ε₁. Check which ε₁ pattern emerges.
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

fn descent_set_str(s: u64, n: u8) -> String {
    let positions: Vec<String> = (0..n - 1)
        .filter(|&i| s & (1 << i) != 0)
        .map(|i| (i + 1).to_string())
        .collect();
    if positions.is_empty() {
        "∅".to_string()
    } else {
        format!("{{{}}}", positions.join(","))
    }
}

/// Check if two polynomials are weakly interlacing (in either direction).
fn weakly_interlace(f: &[i64], g: &[i64]) -> bool {
    if f.len() <= 1 || g.len() <= 1 {
        return true;
    }
    let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
    check_weak_interlacing(small, large) == Some(true)
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    println!("═══ Compatible family test: position-of-max refinement ═══\n");

    let mut total_sets = 0u64;
    let mut compatible_count = 0u64;
    let mut incompatible_count = 0u64;
    let mut rr_count = 0u64;
    let mut not_rr_count = 0u64;

    for n in 3..=max_n {
        let all = all_permutations(n);

        // Group by descent set
        let mut by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for s in &all {
            by_des.entry(descent_set_bitmask(s)).or_default().push(s);
        }

        let mut n_total = 0u32;
        let mut n_compatible = 0u32;
        let mut n_rr = 0u32;
        let mut n_fail_examples: Vec<String> = Vec::new();

        for (&ds, perms) in &by_des {
            let poly = build_poly(&perms);
            let rr = poly.len() <= 2 || is_real_rooted(&poly);
            if rr {
                n_rr += 1;
                rr_count += 1;
            } else {
                not_rr_count += 1;
            }

            // Refine by position of max value n
            let mut by_pos: BTreeMap<usize, Vec<&Vec<u8>>> = BTreeMap::new();
            for s in perms {
                let pos = s.iter().position(|&v| v == n).unwrap() + 1; // 1-indexed
                by_pos.entry(pos).or_default().push(s);
            }

            let pos_polys: Vec<(usize, Vec<i64>)> = by_pos
                .iter()
                .map(|(&p, grp)| (p, build_poly(grp)))
                .collect();

            // Check all pieces real-rooted
            let all_pieces_rr = pos_polys
                .iter()
                .all(|(_, p)| p.len() <= 2 || is_real_rooted(p));

            // Check pairwise weak interlacing
            let mut compatible = true;
            let mut fail_pair = None;
            'outer: for i in 0..pos_polys.len() {
                for j in i + 1..pos_polys.len() {
                    if !weakly_interlace(&pos_polys[i].1, &pos_polys[j].1) {
                        compatible = false;
                        fail_pair = Some((pos_polys[i].0, pos_polys[j].0));
                        break 'outer;
                    }
                }
            }

            n_total += 1;
            total_sets += 1;
            if compatible {
                n_compatible += 1;
                compatible_count += 1;
            } else {
                incompatible_count += 1;
                if n_fail_examples.len() < 3 {
                    let (p1, p2) = fail_pair.unwrap();
                    n_fail_examples.push(format!(
                        "Des={} (pos {},{} fail) pieces_rr:{} L={}",
                        descent_set_str(ds, n),
                        p1,
                        p2,
                        all_pieces_rr,
                        format_poly(&poly),
                    ));
                }
            }
        }

        println!(
            "n={}: {} descent sets, {} rr, {}/{} compatible (pos-of-max)",
            n, n_total, n_rr, n_compatible, n_total,
        );
        for ex in &n_fail_examples {
            println!("  FAIL: {}", ex);
        }
    }

    println!("\n═══ Summary ═══");
    println!(
        "Total: {} descent sets, {} real-rooted, {} not",
        total_sets, rr_count, not_rr_count
    );
    println!(
        "Compatible family (pos-of-max): {} yes, {} no",
        compatible_count, incompatible_count
    );

    if incompatible_count > 0 {
        // Try alternative refinement: position of SECOND largest value
        println!("\n═══ Alternative: refine by (pos_max, pos_second_max) ═══\n");

        for n in 3..=max_n.min(8) {
            let all = all_permutations(n);
            let mut by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
            for s in &all {
                by_des.entry(descent_set_bitmask(s)).or_default().push(s);
            }

            let mut n_total = 0u32;
            let mut n_compatible2 = 0u32;

            for (&ds, perms) in &by_des {
                // Refine by (pos_of_n, pos_of_(n-1))
                let mut by_pos2: BTreeMap<(usize, usize), Vec<&Vec<u8>>> = BTreeMap::new();
                for s in perms {
                    let pos_n = s.iter().position(|&v| v == n).unwrap() + 1;
                    let pos_nm1 = s.iter().position(|&v| v == n - 1).unwrap() + 1;
                    by_pos2.entry((pos_n, pos_nm1)).or_default().push(s);
                }

                let pos2_polys: Vec<((usize, usize), Vec<i64>)> = by_pos2
                    .iter()
                    .map(|(&key, grp)| (key, build_poly(grp)))
                    .collect();

                let mut compatible2 = true;
                for i in 0..pos2_polys.len() {
                    for j in i + 1..pos2_polys.len() {
                        if !weakly_interlace(&pos2_polys[i].1, &pos2_polys[j].1) {
                            compatible2 = false;
                            break;
                        }
                    }
                    if !compatible2 {
                        break;
                    }
                }

                n_total += 1;
                if compatible2 {
                    n_compatible2 += 1;
                }
            }

            println!(
                "n={}: {}/{} compatible (pos_max, pos_second_max)",
                n, n_compatible2, n_total,
            );
        }
    }

    // Part 3: Check insertion-based structure
    println!("\n═══ Insertion structure: ε₁ frequency by descent set ═══\n");
    println!("For each (n, S), count how often ε₁=1 among all insertions.\n");

    for n in 4..=max_n.min(8) {
        let prev_perms = all_permutations(n - 1);
        let all_perms = all_permutations(n);

        // Group target perms by descent set
        let mut by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for s in &all_perms {
            by_des.entry(descent_set_bitmask(s)).or_default().push(s);
        }

        // For each target descent set, analyze insertions
        let mut eps1_always_zero = 0u32;
        let mut eps1_sometimes_one = 0u32;

        for (&ds, _) in &by_des {
            // Determine valid insertion positions for this descent set
            // p is valid if p ∈ S (or p=n) and p-1 ∉ S
            let mut e1_total = 0u64;
            let mut e1_one = 0u64;

            for pi in &prev_perms {
                for p in 1..=n as usize {
                    // Check if inserting at p gives descent set ds
                    let mut sigma = Vec::with_capacity(n as usize);
                    for i in 0..p - 1 {
                        sigma.push(pi[i]);
                    }
                    sigma.push(n);
                    for i in p - 1..pi.len() {
                        sigma.push(pi[i]);
                    }

                    if descent_set_bitmask(&sigma) != ds {
                        continue;
                    }

                    e1_total += 1;

                    // Compute ε₁
                    if p > 1 && p < n as usize {
                        if pi[p - 2] + 1 == pi[p - 1] {
                            e1_one += 1;
                        }
                    }
                }
            }

            if e1_one == 0 {
                eps1_always_zero += 1;
            } else {
                eps1_sometimes_one += 1;
            }
        }

        println!(
            "n={}: ε₁≡0 for {}/{} descent sets, ε₁ sometimes 1 for {}",
            n,
            eps1_always_zero,
            by_des.len(),
            eps1_sometimes_one,
        );
    }
}
