/// Compute and display the transfer matrix for the insertion recurrence.
///
/// For a target descent set S in S_n, the transfer matrix maps
/// {L_{n-1,S'_p}^{(q)}} to {L_{n,S}^{(p)}} via:
///
///   L_{n,S}^{(p)} = sum_q M_{p,q} * L_{n-1,S'_p}^{(q)}
///
/// where M_{p,q} encodes the staircase (t vs 1) and correction terms.
///
/// We display:
/// 1. The "pure staircase" matrix (t^{eps2} entries)
/// 2. The actual transfer including corrections
/// 3. Properties of the matrix (total nonnegativity, etc.)
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{format_poly, is_real_rooted};
use std::collections::{BTreeMap, HashMap};

fn build_poly_from_vals(vals: &[usize]) -> Vec<i64> {
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
        let p_in_s = (s_mask >> (p - 1)) & 1 == 1;
        let pm1_in_s = if p >= 2 {
            (s_mask >> (p - 2)) & 1 == 1
        } else {
            false
        };
        if p_in_s && !pm1_in_s {
            positions.push(p);
        }
    }
    if n >= 2 && (s_mask >> (n - 2)) & 1 == 0 {
        positions.push(n);
    }
    positions
}

fn source_descent_set(s_mask: u64, p: u8, n: u8) -> u64 {
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

fn epsilon1(pi: &[u8], p: u8) -> bool {
    let n = pi.len() as u8 + 1;
    if p <= 1 || p >= n {
        return false;
    }
    let val_left = pi[(p - 2) as usize];
    let val_right = pi[(p - 1) as usize];
    val_left + 1 == val_right
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

/// Polynomial addition
fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut result = vec![0i64; len];
    for (i, &c) in a.iter().enumerate() {
        result[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        result[i] += c;
    }
    while result.len() > 1 && *result.last().unwrap() == 0 {
        result.pop();
    }
    result
}

/// Multiply polynomial by t (shift coefficients)
fn poly_mul_t(a: &[i64]) -> Vec<i64> {
    if a.is_empty() || (a.len() == 1 && a[0] == 0) {
        return vec![0];
    }
    let mut result = vec![0i64; a.len() + 1];
    for (i, &c) in a.iter().enumerate() {
        result[i + 1] = c;
    }
    result
}

/// Multiply polynomial by (t-1)
fn poly_mul_tm1(a: &[i64]) -> Vec<i64> {
    if a.is_empty() || (a.len() == 1 && a[0] == 0) {
        return vec![0];
    }
    // (t-1)*a = t*a - a
    let ta = poly_mul_t(a);
    let mut result = vec![0i64; ta.len()];
    for (i, &c) in ta.iter().enumerate() {
        result[i] += c;
    }
    for (i, &c) in a.iter().enumerate() {
        result[i] -= c;
    }
    while result.len() > 1 && *result.last().unwrap() == 0 {
        result.pop();
    }
    result
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    for n in 4..=max_n {
        let perms = all_permutations(n);
        let perms_prev = all_permutations(n - 1);

        // Group current perms by descent set
        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        // Group previous perms by descent set
        let mut by_descent_prev: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            by_descent_prev
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        // Pick interesting descent sets to display
        // Focus on small examples and alternating
        let interesting: Vec<u64> = by_descent
            .keys()
            .filter(|&&mask| {
                // 1 not in S
                (mask & 1) == 0
            })
            .copied()
            .collect();

        for &mask in &interesting {
            let s_str = descent_set_to_string(mask, n);
            let vp = valid_positions(mask, n);

            if vp.len() < 2 {
                continue;
            }

            // For each valid p, compute source descent set and check
            // if all source sets are the same
            let source_sets: Vec<(u8, u64)> = vp
                .iter()
                .map(|&p| (p, source_descent_set(mask, p, n)))
                .collect();

            let all_same_source = source_sets.iter().all(|(_, s)| *s == source_sets[0].1);

            // Compute position-refined polynomials for source
            // We need: for each (p, q), the contribution to L_{n,S}^{(p)}
            // from source permutations with max at q

            println!(
                "=== n={}, S={} (|P(S)|={}{}) ===",
                n,
                s_str,
                vp.len(),
                if all_same_source {
                    ", same source"
                } else {
                    ", VARYING sources"
                }
            );

            for &(p, sp_mask) in &source_sets {
                let sp_str = descent_set_to_string(sp_mask, n - 1);
                let vq = valid_positions(sp_mask, n - 1);

                let source_class = match by_descent_prev.get(&sp_mask) {
                    Some(c) => c,
                    None => {
                        println!("  p={}: empty source S'={}", p, sp_str);
                        continue;
                    }
                };

                println!("  p={}: S'_p={}, valid q={:?}", p, sp_str, vq);

                // For each q, compute:
                // - L_{n-1,S'}^{(q)}: source poly at position q
                // - C_{p,q}: correction (epsilon1=1 subset)
                // - The actual contribution to L_{n,S}^{(p)}
                //   = t^{eps2} * [L^{(q)} + (t-1)*C_{p,q}]

                let mut staircase_row = String::new();
                let mut actual_row = String::new();

                for &q in &vq {
                    let source_at_q: Vec<&Vec<u8>> = source_class
                        .iter()
                        .filter(|pi| {
                            let pos = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            pos == q
                        })
                        .copied()
                        .collect();

                    if source_at_q.is_empty() {
                        continue;
                    }

                    let l_vals: Vec<usize> = source_at_q
                        .iter()
                        .map(|pi| compute(pi, Stat::Swaps))
                        .collect();
                    let l_poly = build_poly_from_vals(&l_vals);

                    let c_vals: Vec<usize> = source_at_q
                        .iter()
                        .filter(|pi| epsilon1(pi, p))
                        .map(|pi| compute(pi, Stat::Swaps))
                        .collect();
                    let c_poly = build_poly_from_vals(&c_vals);

                    let eps2 = if q <= p.saturating_sub(2) { 1 } else { 0 };

                    // Staircase entry (no correction): t^{eps2} * L^{(q)}
                    let staircase_entry = if eps2 == 1 {
                        poly_mul_t(&l_poly)
                    } else {
                        l_poly.clone()
                    };

                    // Actual entry: t^{eps2} * [L^{(q)} + (t-1)*C_{p,q}]
                    let corrected = poly_add(&l_poly, &poly_mul_tm1(&c_poly));
                    let actual_entry = if eps2 == 1 {
                        poly_mul_t(&corrected)
                    } else {
                        corrected
                    };

                    let has_correction = c_vals.iter().any(|_| true);

                    println!(
                        "    q={}: eps2={}, L^(q)={}, C={}{}, contribution={}",
                        q,
                        eps2,
                        format_poly(&l_poly),
                        format_poly(&c_poly),
                        if has_correction { "" } else { " (zero)" },
                        format_poly(&actual_entry),
                    );
                }

                // Also show the total L_{n,S}^{(p)}
                let target_vals: Vec<usize> = by_descent[&mask]
                    .iter()
                    .filter(|pi| {
                        let pos = pi.iter().position(|&v| v == n).unwrap() as u8 + 1;
                        pos == p
                    })
                    .map(|pi| compute(pi, Stat::Swaps))
                    .collect();
                let target_poly = build_poly_from_vals(&target_vals);
                println!("    => L_{{n,S}}^({}) = {}", p, format_poly(&target_poly));
            }

            // Show the "matrix" in compact form
            // For each p, the contribution is: sum_q M_{p,q}(t) * f_q
            // where M_{p,q}(t) is a polynomial multiplier
            println!("  Transfer matrix M_{{p,q}}(t):");
            for &(p, sp_mask) in &source_sets {
                let vq = valid_positions(sp_mask, n - 1);
                let source_class = match by_descent_prev.get(&sp_mask) {
                    Some(c) => c,
                    None => continue,
                };

                let mut row_entries: Vec<String> = Vec::new();
                for &q in &vq {
                    let source_at_q: Vec<&Vec<u8>> = source_class
                        .iter()
                        .filter(|pi| {
                            let pos = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            pos == q
                        })
                        .copied()
                        .collect();
                    if source_at_q.is_empty() {
                        continue;
                    }

                    let eps2 = if q <= p.saturating_sub(2) { 1 } else { 0 };
                    let c_count: usize = source_at_q.iter().filter(|pi| epsilon1(pi, p)).count();
                    let total = source_at_q.len();

                    if c_count == 0 {
                        // Pure staircase
                        row_entries.push(if eps2 == 1 {
                            format!("q={}: t", q)
                        } else {
                            format!("q={}: 1", q)
                        });
                    } else {
                        // Has correction
                        // The multiplier is t^{eps2} * [1 + (t-1)*c/l]
                        // where c/l is the fraction with epsilon1=1
                        // But this isn't a scalar - it depends on the poly
                        row_entries.push(format!(
                            "q={}: {}*(1+(t-1)*[{}/{}])",
                            q,
                            if eps2 == 1 { "t" } else { "1" },
                            c_count,
                            total,
                        ));
                    }
                }
                println!("    p={}: {}", p, row_entries.join(", "));
            }
            println!();
        }
    }
}
