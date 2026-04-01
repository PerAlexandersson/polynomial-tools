/// Approach 8 combined with Approach 6:
/// Can we prove the TP_2 property inductively?
///
/// The recurrence gives L_{n,S}^{(p)} in terms of L_{n-1,T}^{(q)} for various T.
/// Question: if the coefficient matrices at level n-1 are TP_2 for all T,
/// does the staircase + correction transfer preserve TP_2 at level n?
///
/// This experiment:
/// 1. For each (n, S), compute the coefficient matrix M of L^{(p)}.
/// 2. Check TP_2 of M.
/// 3. Check if M can be factored as a product of "simple" TP_2 matrices.
/// 4. Check if the individual 2x2 minors of M have a COMBINATORIAL formula
///    in terms of the source data.
///
/// KEY INSIGHT: If {L^{(p)}} is an interlacing sequence (not just compatible family),
/// then Conjecture 6.7 (compatible family) is automatically satisfied!
/// And TP_2 of the coefficient matrix PROVES interlacing.
///
/// So the new conjecture to prove is:
///   The coefficient matrix of (L^{(p)})_{p in P(S)} is TP_2.
///
/// For the INDUCTIVE STEP, we need:
///   TP_2 at level n-1 (for all source S'_p, S''_p) => TP_2 at level n.
///
/// The staircase matrix maps interlacing sequences to interlacing sequences,
/// but we need to verify that the COMBINED effect (staircase + correction + mixing)
/// preserves TP_2.
///
/// Test: Compute the 2x2 minors explicitly and check if they factor nicely.
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

fn coeff(poly: &[i64], k: usize) -> i64 {
    if k < poly.len() {
        poly[k]
    } else {
        0
    }
}

/// Compute the 2x2 minor M_{i,k1} * M_{j,k2} - M_{i,k2} * M_{j,k1}
fn minor_2x2(rows: &[Vec<i64>], i: usize, j: usize, k1: usize, k2: usize) -> i64 {
    coeff(&rows[i], k1) * coeff(&rows[j], k2) - coeff(&rows[i], k2) * coeff(&rows[j], k1)
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    println!("=== TP_2 Inductive Analysis ===\n");

    // For each (n, S), look at the 2x2 minors and try to express them
    // in terms of source quantities.

    for n in 5..=max_n {
        let perms = all_permutations(n);
        let perms_prev = all_permutations(n - 1);

        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        let mut prev_by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        println!("--- n = {} ---", n);

        // Focus on examples with exactly 2 valid positions including p=2
        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }
            if !vp.contains(&2) {
                continue;
            }

            // Only print detailed analysis for small families
            if vp.len() > 4 {
                continue;
            }

            let s_str = descent_set_to_string(mask, n);

            // Build L^{(p)} for each p
            let mut lp_polys: Vec<(u8, Vec<i64>)> = Vec::new();
            for &p in &vp {
                let vals: Vec<usize> = class
                    .iter()
                    .filter(|pi| pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                    .map(|pi| compute(pi, Stat::Swaps))
                    .collect();
                if !vals.is_empty() {
                    lp_polys.push((p, build_poly(&vals)));
                }
            }

            if lp_polys.len() < 2 {
                continue;
            }

            println!("\nS = {} (positions {:?})", s_str, vp);
            for (p, poly) in &lp_polys {
                println!("  L^({}) = {}", p, format_poly(poly));
            }

            // Print 2x2 minors between p=2 and each p>=3
            let rows: Vec<Vec<i64>> = lp_polys.iter().map(|(_, poly)| poly.clone()).collect();
            let max_deg = rows.iter().map(|r| r.len()).max().unwrap_or(0);

            // Find the index of p=2
            let p2_idx = lp_polys.iter().position(|(p, _)| *p == 2);
            if let Some(p2_idx) = p2_idx {
                for (j_idx, (p_j, _)) in lp_polys.iter().enumerate() {
                    if *p_j <= 2 {
                        continue;
                    }
                    print!("  Minors L^(2) vs L^({}): ", p_j);
                    for k1 in 0..max_deg {
                        for k2 in (k1 + 1)..max_deg {
                            let m = minor_2x2(&rows, p2_idx, j_idx, k1, k2);
                            if m != 0 {
                                print!("M[{},{}]={} ", k1, k2, m);
                            }
                        }
                    }
                    println!();
                }
            }

            // Now decompose: what are the SOURCE polynomials for each p?
            for &p in &vp {
                if p <= 1 {
                    continue;
                }
                let sp_a = source_asc(mask, p, n);
                let sp_d = source_desc(mask, p, n);

                let sp_a_str = descent_set_to_string(sp_a, n - 1);
                let sp_d_str = sp_d
                    .map(|d| descent_set_to_string(d, n - 1))
                    .unwrap_or("N/A".to_string());

                println!("  Source for p={}: S'={}, S''={}", p, sp_a_str, sp_d_str);

                // Check if S'_p is the same for all p
            }

            // Check: do all valid p share the same source descent set?
            let sources: Vec<(u64, Option<u64>)> = vp
                .iter()
                .map(|&p| (source_asc(mask, p, n), source_desc(mask, p, n)))
                .collect();
            let same_asc = sources.iter().all(|(a, _)| *a == sources[0].0);
            println!("  Same ascent source for all p: {}", same_asc);
        }
    }

    // Summary: check TP_2 of coefficient matrix for CUMULATIVE sums
    println!("\n\n=== TP_2 of cumulative sums ===");
    println!("  C_k = L^(p_1) + ... + L^(p_k) for p_1 < ... < p_m\n");

    for n in 4..=max_n {
        let perms = all_permutations(n);
        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        let mut cumul_tp2_ok = 0u32;
        let mut cumul_tp2_total = 0u32;

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            let mut lp_polys: Vec<Vec<i64>> = Vec::new();
            for &p in &vp {
                let vals: Vec<usize> = class
                    .iter()
                    .filter(|pi| pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                    .map(|pi| compute(pi, Stat::Swaps))
                    .collect();
                if !vals.is_empty() {
                    lp_polys.push(build_poly(&vals));
                }
            }

            if lp_polys.len() < 2 {
                continue;
            }

            // Build cumulative sums
            let max_deg = lp_polys.iter().map(|r| r.len()).max().unwrap_or(0);
            let mut cumul_rows: Vec<Vec<i64>> = Vec::new();
            let mut running = vec![0i64; max_deg];
            for poly in &lp_polys {
                for (i, &c) in poly.iter().enumerate() {
                    running[i] += c;
                }
                cumul_rows.push(running.clone());
            }

            cumul_tp2_total += 1;
            // Check TP_2 of cumulative sum matrix
            let mut is_tp2 = true;
            for i in 0..cumul_rows.len() {
                for j in (i + 1)..cumul_rows.len() {
                    for k1 in 0..max_deg {
                        for k2 in (k1 + 1)..max_deg {
                            let minor = coeff(&cumul_rows[i], k1) * coeff(&cumul_rows[j], k2)
                                - coeff(&cumul_rows[i], k2) * coeff(&cumul_rows[j], k1);
                            if minor < 0 {
                                is_tp2 = false;
                            }
                        }
                    }
                }
            }
            if is_tp2 {
                cumul_tp2_ok += 1;
            }
        }

        println!(
            "n={}: cumulative TP_2 {}/{}",
            n, cumul_tp2_ok, cumul_tp2_total
        );
    }
}
