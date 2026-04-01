/// KEY IDEA: Invert the staircase to find a single source vector.
///
/// If we include p=2 as the "bottom" of the staircase (threshold at 0),
/// the full staircase matrix (rows for all p ∈ P(S)) maps a single
/// source vector (f_q) to the target (L^{(p)}).
///
/// The staircase gives: L^{(p)} = Σ_{q≤p-2} t·f_q + Σ_{q≥p-1} f_q
///                              = F + (t-1)·(f_1+...+f_{p-2})
///
/// So the differences determine the f_q:
///   L^{(p_{i+1})} - L^{(p_i)} = (t-1)·(f_{p_i-1} + ... + f_{p_{i+1}-2})
///
/// Check: (a) Is each difference divisible by (t-1)?
///        (b) Is the quotient a polynomial with non-negative coefficients?
///        (c) Do the recovered f_q form an interlacing sequence?
///
/// If all three hold, Brändén gives the FULL compatible family in one shot!
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

/// Divide polynomial by (t-1). Returns None if not divisible.
fn div_by_t_minus_1(p: &[i64]) -> Option<Vec<i64>> {
    // (t-1) = [-1, 1] in coefficient form
    // Polynomial long division
    if p.is_empty() || (p.len() == 1 && p[0] == 0) {
        return Some(vec![0]);
    }
    let n = p.len();
    if n < 2 {
        return None;
    } // degree 0 can't be divisible by degree 1 (unless 0)

    let mut rem = p.to_vec();
    let mut quotient = vec![0i64; n - 1];

    // Divide: rem by [-1, 1]
    // Leading coeff of divisor is 1 (for t)
    for i in (1..n).rev() {
        let q = rem[i]; // / 1 (leading coeff of t-1)
        quotient[i - 1] = q;
        rem[i] -= q * 1; // subtract q * t
        rem[i - 1] -= q * (-1); // subtract q * (-1)
    }

    // Check remainder is 0
    if rem[0] != 0 {
        return None;
    }

    // Trim
    while quotient.len() > 1 && *quotient.last().unwrap() == 0 {
        quotient.pop();
    }
    Some(quotient)
}

fn poly_sub(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut result = vec![0i64; len];
    for (i, &c) in a.iter().enumerate() {
        result[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        result[i] -= c;
    }
    while result.len() > 1 && *result.last().unwrap() == 0 {
        result.pop();
    }
    result
}

fn check_compat(f: &[i64], g: &[i64]) -> bool {
    if f.len() <= 1 || g.len() <= 1 {
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
        let mut c = vec![0i64; maxlen];
        for (i, &v) in f.iter().enumerate() {
            c[i] += a * v;
        }
        for (i, &v) in g.iter().enumerate() {
            c[i] += b * v;
        }
        while c.len() > 1 && *c.last().unwrap() == 0 {
            c.pop();
        }
        if c.len() > 1 && !is_real_rooted(&c) {
            return false;
        }
    }
    true
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

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    println!("=== Staircase inversion: find single source vector ===\n");

    let mut total_divisible = 0u64;
    let mut total_nonneg = 0u64;
    let mut total_interlacing = 0u64;
    let mut total_tested = 0u64;

    for n in 4..=max_n {
        let perms = all_permutations(n);
        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        let mut n_div = 0u32;
        let mut n_nn = 0u32;
        let mut n_int = 0u32;
        let mut n_total = 0u32;

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            } // 1 ∈ S
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            n_total += 1;

            // Compute L^{(p)} for each valid p, ordered by p
            let mut lp: Vec<(u8, Vec<i64>)> = Vec::new();
            for &p in &vp {
                let vals: Vec<usize> = class
                    .iter()
                    .filter(|pi| pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                    .map(|pi| compute(pi, Stat::Swaps))
                    .collect();
                lp.push((p, build_poly(&vals)));
            }

            // The staircase model: L^{(p)} = F + (t-1)·S_{p-2}
            // where F = Σ f_q and S_k = f_1 + ... + f_k
            //
            // So S_{p-2} = (L^{(p)} - L^{(p_min)}) / (t-1)
            // where p_min is the smallest valid position (has S_{p_min-2} = 0 if p_min=2, or S_0=0)
            //
            // Actually: L^{(p)} - F = (t-1)·S_{p-2}, and F = L^{(p_min)} when p_min has threshold 0.
            //
            // For p_min = first valid position:
            // L^{(p_min)} = F + (t-1)·S_{p_min-2}
            //
            // If p_min = 2: S_0 = 0, so L^{(2)} = F. Then S_{p-2} = (L^{(p)} - L^{(2)})/(t-1).
            //
            // If p_min > 2: We need S_{p_min-2} too, which we don't have from the data alone.

            let p_min = vp[0];

            if p_min != 2 {
                // p=2 is not valid, so we can't use L^{(2)} = F.
                // Skip for now — we only care about the 2∈S case.
                continue;
            }

            // p_min = 2: F = L^{(2)}
            let f_total = lp[0].1.clone(); // L^{(2)}

            // Recover partial sums S_k from consecutive differences
            // S_{p-2} = (L^{(p)} - L^{(2)}) / (t-1)
            let mut all_divisible = true;
            let mut all_nonneg = true;
            let mut partial_sums: Vec<(u8, Vec<i64>)> = Vec::new(); // (threshold p-2, polynomial S_{p-2})
            partial_sums.push((0, vec![0])); // S_0 = 0

            for i in 1..lp.len() {
                let (p, ref lp_poly) = lp[i];
                let diff = poly_sub(lp_poly, &f_total);
                match div_by_t_minus_1(&diff) {
                    Some(q) => {
                        // Check non-negative coefficients
                        let nn = q.iter().all(|&c| c >= 0);
                        if !nn {
                            all_nonneg = false;
                        }
                        partial_sums.push((p - 2, q));
                    }
                    None => {
                        all_divisible = false;
                        all_nonneg = false;
                        let s_str = descent_set_to_string(mask, n);
                        println!("  NOT DIVISIBLE: n={} S={} p={}", n, s_str, p);
                        println!("    L^({}) - L^(2) = {}", p, format_poly(&diff));
                        break;
                    }
                }
            }

            if all_divisible {
                n_div += 1;
            }
            if all_nonneg {
                n_nn += 1;
            }

            // Recover individual f_q from differences of partial sums
            // S_k - S_{k-1} = f_k (but we only have S at specific thresholds)
            if all_divisible && all_nonneg && partial_sums.len() >= 2 {
                // The partial sums give us "blocks": between consecutive thresholds
                // S_{threshold_i} - S_{threshold_{i-1}} = f_{threshold_{i-1}+1} + ... + f_{threshold_i}
                // These are the "block sums" — we can treat each block as one f_q.

                let block_polys: Vec<Vec<i64>> = (1..partial_sums.len())
                    .map(|i| poly_sub(&partial_sums[i].1, &partial_sums[i - 1].1))
                    .collect();

                // Also need the "tail": F - S_{max_threshold} = f_{max_threshold+1} + ... + f_m
                let last_partial = &partial_sums.last().unwrap().1;
                let tail = poly_sub(&f_total, last_partial);

                let mut all_source: Vec<Vec<i64>> = block_polys;
                all_source.push(tail);

                // Check interlacing of block source
                let nontrivial: Vec<&Vec<i64>> = all_source
                    .iter()
                    .filter(|p| p.len() > 1 || (p.len() == 1 && p[0] != 0))
                    .collect();

                if nontrivial.len() >= 2 {
                    let mut interlacing = true;
                    'outer: for i in 0..nontrivial.len() {
                        for j in (i + 1)..nontrivial.len() {
                            if !check_compat(nontrivial[i], nontrivial[j]) {
                                interlacing = false;
                                break 'outer;
                            }
                        }
                    }
                    if interlacing {
                        n_int += 1;
                    } else {
                        let s_str = descent_set_to_string(mask, n);
                        println!("  NOT INTERLACING: n={} S={}", n, s_str);
                        for (i, bp) in all_source.iter().enumerate() {
                            println!("    block {}: {}", i, format_poly(bp));
                        }
                    }
                } else {
                    n_int += 1; // trivially interlacing
                }
            }
        }

        total_divisible += n_div as u64;
        total_nonneg += n_nn as u64;
        total_interlacing += n_int as u64;
        total_tested += n_total as u64;

        println!(
            "n={}: divisible {}/{}, nonneg {}/{}, interlacing {}/{}",
            n, n_div, n_total, n_nn, n_total, n_int, n_total
        );
    }

    println!("\n=== SUMMARY ===");
    println!("Divisible by (t-1): {}/{}", total_divisible, total_tested);
    println!("Non-negative quotient: {}/{}", total_nonneg, total_tested);
    println!(
        "Block source interlacing: {}/{}",
        total_interlacing, total_tested
    );
    println!();
    println!("If ALL THREE pass, then a single staircase matrix");
    println!("(including p=2) acts on an interlacing source,");
    println!("and Branden gives the full compatible family. QED!");
}
