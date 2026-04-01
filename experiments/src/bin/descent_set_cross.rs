/// Cross-compatibility of L_{n,S}(t) across different descent sets S.
/// Also test alternative refinements for the failing descent sets.
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly, is_real_rooted};
use std::collections::BTreeMap;

fn build_poly(perms: &[&Vec<u8>]) -> Vec<i64> {
    if perms.is_empty() {
        return vec![0];
    }
    let max_s = perms
        .iter()
        .map(|s| compute(s, Stat::Swaps))
        .max()
        .unwrap_or(0);
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
    if positions.is_empty() {
        "∅".to_string()
    } else {
        format!("{{{}}}", positions.join(","))
    }
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    // Part 1: Cross-compatibility of L_{n,S} across different S
    println!("═══ Part 1: Are all L_{{n,S}} mutually compatible? ═══\n");

    for n in 3..=max_n {
        let all = all_permutations(n);
        let mut by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for s in &all {
            by_des.entry(descent_set_bitmask(s)).or_default().push(s);
        }

        let polys: Vec<(u64, Vec<i64>)> = by_des
            .iter()
            .map(|(&ds, perms)| (ds, build_poly(perms)))
            .filter(|(_, p)| p.len() > 1) // skip constant polys
            .collect();

        let mut all_compatible = true;
        let mut fail_count = 0u32;
        for i in 0..polys.len() {
            for j in i + 1..polys.len() {
                if !weakly_interlace(&polys[i].1, &polys[j].1) {
                    all_compatible = false;
                    fail_count += 1;
                    if fail_count <= 3 {
                        println!(
                            "  n={}: L_{{{}}} ~ L_{{{}}}: ✗",
                            n,
                            descent_set_str(polys[i].0, n),
                            descent_set_str(polys[j].0, n),
                        );
                    }
                }
            }
        }
        let total_pairs = polys.len() * (polys.len() - 1) / 2;
        if all_compatible {
            println!("n={}: all {} pairs compatible ✓", n, total_pairs,);
        } else {
            println!("n={}: {} / {} pairs fail", n, fail_count, total_pairs,);
        }
    }

    // Part 2: For failing descent sets, try refinement by LAST value σ(n)
    println!("\n═══ Part 2: Refine by last value σ(n) ═══\n");

    for n in 6..=max_n {
        let all = all_permutations(n);
        let mut by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for s in &all {
            by_des.entry(descent_set_bitmask(s)).or_default().push(s);
        }

        let mut n_compat = 0u32;
        let mut n_total = 0u32;

        for (&ds, perms) in &by_des {
            // Only test descent sets where pos-of-max failed (1 ∈ S)
            if ds & 1 == 0 {
                continue;
            }
            n_total += 1;

            let mut by_last: BTreeMap<u8, Vec<&Vec<u8>>> = BTreeMap::new();
            for s in perms {
                by_last.entry(s[n as usize - 1]).or_default().push(s);
            }

            let last_polys: Vec<(u8, Vec<i64>)> = by_last
                .iter()
                .map(|(&v, grp)| (v, build_poly(grp)))
                .collect();

            let mut compatible = true;
            for i in 0..last_polys.len() {
                for j in i + 1..last_polys.len() {
                    if !weakly_interlace(&last_polys[i].1, &last_polys[j].1) {
                        compatible = false;
                        break;
                    }
                }
                if !compatible {
                    break;
                }
            }
            if compatible {
                n_compat += 1;
            }
        }

        println!(
            "n={}: {}/{} descent sets with 1∈S compatible (last-value refine)",
            n, n_compat, n_total,
        );
    }

    // Part 3: Refine by FIRST value σ(1)
    println!("\n═══ Part 3: Refine by first value σ(1) ═══\n");

    for n in 6..=max_n {
        let all = all_permutations(n);
        let mut by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for s in &all {
            by_des.entry(descent_set_bitmask(s)).or_default().push(s);
        }

        let mut n_compat = 0u32;
        let mut n_total = 0u32;

        for (&ds, perms) in &by_des {
            if ds & 1 == 0 {
                continue;
            }
            n_total += 1;

            let mut by_first: BTreeMap<u8, Vec<&Vec<u8>>> = BTreeMap::new();
            for s in perms {
                by_first.entry(s[0]).or_default().push(s);
            }

            let first_polys: Vec<(u8, Vec<i64>)> = by_first
                .iter()
                .map(|(&v, grp)| (v, build_poly(grp)))
                .collect();

            let mut compatible = true;
            for i in 0..first_polys.len() {
                for j in i + 1..first_polys.len() {
                    if !weakly_interlace(&first_polys[i].1, &first_polys[j].1) {
                        compatible = false;
                        break;
                    }
                }
                if !compatible {
                    break;
                }
            }
            if compatible {
                n_compat += 1;
            }
        }

        println!(
            "n={}: {}/{} descent sets with 1∈S compatible (first-value refine)",
            n, n_compat, n_total,
        );
    }

    // Part 4: Refine by number of adjacencies
    println!("\n═══ Part 4: Refine by number of adjacencies ═══\n");

    for n in 3..=max_n {
        let all = all_permutations(n);
        let mut by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for s in &all {
            by_des.entry(descent_set_bitmask(s)).or_default().push(s);
        }

        let mut n_compat = 0u32;
        let mut n_total = 0u32;

        for (_, perms) in &by_des {
            n_total += 1;

            // Count ascending adjacencies: positions i where σ(i)+1 = σ(i+1)
            let adjacencies = |sigma: &[u8]| -> usize {
                (0..sigma.len() - 1)
                    .filter(|&i| sigma[i] + 1 == sigma[i + 1])
                    .count()
            };

            let mut by_adj: BTreeMap<usize, Vec<&Vec<u8>>> = BTreeMap::new();
            for s in perms {
                by_adj.entry(adjacencies(s)).or_default().push(s);
            }

            let adj_polys: Vec<(usize, Vec<i64>)> = by_adj
                .iter()
                .map(|(&a, grp)| (a, build_poly(grp)))
                .collect();

            let mut compatible = true;
            for i in 0..adj_polys.len() {
                for j in i + 1..adj_polys.len() {
                    if !weakly_interlace(&adj_polys[i].1, &adj_polys[j].1) {
                        compatible = false;
                        break;
                    }
                }
                if !compatible {
                    break;
                }
            }
            if compatible {
                n_compat += 1;
            }
        }

        println!(
            "n={}: {}/{} compatible (adjacency-count refine)",
            n, n_compat, n_total,
        );
    }

    // Part 5: Hybrid: refine by (pos_max, number_of_adjacencies)
    println!("\n═══ Part 5: Refine by (pos_max, adj_count) ═══\n");

    for n in 3..=max_n {
        let all = all_permutations(n);
        let mut by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for s in &all {
            by_des.entry(descent_set_bitmask(s)).or_default().push(s);
        }

        let adjacencies = |sigma: &[u8]| -> usize {
            (0..sigma.len() - 1)
                .filter(|&i| sigma[i] + 1 == sigma[i + 1])
                .count()
        };

        let mut n_compat = 0u32;
        let mut n_total = 0u32;

        for (_, perms) in &by_des {
            n_total += 1;

            let mut by_hybrid: BTreeMap<(usize, usize), Vec<&Vec<u8>>> = BTreeMap::new();
            for s in perms {
                let pos = s.iter().position(|&v| v == n).unwrap() + 1;
                let adj = adjacencies(s);
                by_hybrid.entry((pos, adj)).or_default().push(s);
            }

            let hybrid_polys: Vec<((usize, usize), Vec<i64>)> = by_hybrid
                .iter()
                .map(|(&key, grp)| (key, build_poly(grp)))
                .collect();

            let mut compatible = true;
            for i in 0..hybrid_polys.len() {
                for j in i + 1..hybrid_polys.len() {
                    if !weakly_interlace(&hybrid_polys[i].1, &hybrid_polys[j].1) {
                        compatible = false;
                        break;
                    }
                }
                if !compatible {
                    break;
                }
            }
            if compatible {
                n_compat += 1;
            }
        }

        println!(
            "n={}: {}/{} compatible (pos_max, adj_count)",
            n, n_compat, n_total,
        );
    }
}
