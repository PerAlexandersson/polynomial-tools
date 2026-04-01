/// Check whether the "total modified source" T_p = sum_q F~_q^{(p)} is the same
/// polynomial for all valid p.
///
/// For each target descent set S (with 1 not in S, 2 in S) and each valid p in P(S):
///   - The modified source for position p consists of source permutations (from
///     both S'_p and S''_p) with statistic swaps + eps1.
///   - F~_q^{(p)} = sum_{pi in source(p), pi^{-1}(n-1)=q} t^{swaps(pi)+eps1(pi,p)}
///   - T_p = sum_q F~_q^{(p)} is the "total modified source" (unweighted sum).
///
/// Question: is T_p the same for all p? If so, the abstract staircase argument
/// would close the gap for p=2 vs p>=3 compatibility.
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

fn poly_equal(a: &[i64], b: &[i64]) -> bool {
    let la = a.len();
    let lb = b.len();
    let maxl = la.max(lb);
    for i in 0..maxl {
        let va = if i < la { a[i] } else { 0 };
        let vb = if i < lb { b[i] } else { 0 };
        if va != vb {
            return false;
        }
    }
    true
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    println!("=== Check: is T_p = sum_q F~_q^(p) the same for all valid p? ===\n");
    println!("T_p is the total (unweighted) modified source polynomial for position p.");
    println!("If T_p is independent of p, the abstract staircase argument closes the gap.\n");

    let mut total_same = 0u64;
    let mut total_diff = 0u64;

    for n in 5..=max_n {
        let perms_prev = all_permutations(n - 1);

        let mut prev_by_des: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
        let mut perm_prev_list: Vec<&Vec<u8>> = Vec::new();
        for pi in &perms_prev {
            let idx = perm_prev_list.len();
            perm_prev_list.push(pi);
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(idx);
        }

        // We need to iterate over valid descent sets of S_n
        // Collect all descent sets that appear
        let perms = all_permutations(n);
        let mut all_masks: std::collections::BTreeSet<u64> = std::collections::BTreeSet::new();
        for pi in &perms {
            all_masks.insert(descent_set_bitmask(pi));
        }

        let mut n_same = 0u32;
        let mut n_diff = 0u32;

        for &mask in &all_masks {
            // Filter: 1 not in S, 2 in S
            if mask & 1 != 0 {
                continue;
            }
            if mask & 2 == 0 {
                continue;
            }

            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            // For each valid p, compute T_p = sum_q F~_q^{(p)}
            let mut t_polys: Vec<(u8, Vec<i64>)> = Vec::new();

            for &p in &vp {
                let sp_a = source_asc(mask, p, n);
                let sp_d = source_desc(mask, p, n);

                let mut all_vals: Vec<usize> = Vec::new();

                // Ascent source
                if let Some(indices) = prev_by_des.get(&sp_a) {
                    for &idx in indices {
                        let pi = perm_prev_list[idx];
                        let e1 = if epsilon1(pi, p) { 1 } else { 0 };
                        all_vals.push(compute(pi, Stat::Swaps) + e1);
                    }
                }

                // Descent source
                if let Some(sp) = sp_d {
                    if let Some(indices) = prev_by_des.get(&sp) {
                        for &idx in indices {
                            let pi = perm_prev_list[idx];
                            // eps1 is always 0 for descent source
                            all_vals.push(compute(pi, Stat::Swaps));
                        }
                    }
                }

                let t_p = build_poly(&all_vals);
                t_polys.push((p, t_p));
            }

            // Check if all T_p are the same
            if t_polys.len() < 2 {
                continue;
            }

            let ref_poly = &t_polys[0].1;
            let all_same = t_polys[1..]
                .iter()
                .all(|(_, poly)| poly_equal(ref_poly, poly));

            let s_str = descent_set_to_string(mask, n);

            if all_same {
                n_same += 1;
                println!(
                    "  n={} S={} P(S)={:?}: ALL SAME, T = {}",
                    n,
                    s_str,
                    vp,
                    format_poly(ref_poly)
                );
            } else {
                n_diff += 1;
                println!("  n={} S={} P(S)={:?}: DIFFERENT!", n, s_str, vp);
                for (p, poly) in &t_polys {
                    println!("    T_{} = {}", p, format_poly(poly));
                }
            }
        }

        total_same += n_same as u64;
        total_diff += n_diff as u64;

        println!(
            "n={}: same={}, different={} (out of {} descent sets with 1 not in S, 2 in S)",
            n,
            n_same,
            n_diff,
            n_same + n_diff
        );
        println!();
    }

    println!("=== SUMMARY ===");
    println!("Total same:      {}", total_same);
    println!("Total different: {}", total_diff);
    println!();
    if total_diff == 0 {
        println!("RESULT: T_p is always the same for all p. The abstract argument works!");
        println!("Since all T_p = F, and F = g_0 in the staircase family {{g_0,...,g_m}},");
        println!("F is compatible with each g_k by interlacing. QED.");
    } else {
        println!("RESULT: T_p differs across p values in some cases.");
        println!("The abstract staircase argument does NOT directly apply.");
        println!("Need a different approach to close the gap.");
    }
}
