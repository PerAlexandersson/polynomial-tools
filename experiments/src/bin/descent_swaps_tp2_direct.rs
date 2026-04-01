/// Direct TP_2 induction test.
///
/// The strongest version of Approach 8:
/// Prove by induction on n that for all S with 1 not in S:
///   (A) L_{n,S}(t) is real-rooted
///   (B) The coefficient matrix of {L_{n,S}^{(p)}}_{p in P(S)} is TP_2
///
/// At level n, L^{(p)} is computed from source polynomials at level n-1.
/// Question: does (B) at level n-1 (for all source descent sets) imply (B) at level n?
///
/// The recurrence gives:
///   L_{n,S}^{(p)} = sum_q t^{eps_2(p,q)} F~_q^{(p)}
///
/// where F~_q^{(p)} involves permutations from S'_p and S''_p.
///
/// Key insight: the OUTPUT TP_2 property is a statement about the L^{(p)}
/// polynomials themselves, NOT about the source vectors.
///
/// Test: Can we express the 2x2 minors of the output coefficient matrix
/// directly in terms of quantities at level n-1?
///
/// Specifically, for rows p_1 < p_2 and columns k_1 < k_2:
///   M(p_1,p_2; k_1,k_2) = [t^k_1]L^{(p_1)} * [t^k_2]L^{(p_2)} - [t^k_2]L^{(p_1)} * [t^k_1]L^{(p_2)}
///
/// If we can show this is >= 0 using the inductive hypothesis, we're done.
///
/// Special focus: minors involving p=2.
/// Since L^{(2)} = sum_q F~_q (no staircase weighting), and
/// L^{(p)} = sum_q t^{eps_2(p,q)} F~_q for p >= 3,
/// the minor decomposes into a sum over pairs (q1, q2) of source positions.
///
/// Actually: L^{(p)} = F_p + (t-1) * F_{p, <= p-2}
/// where F_p = sum_q F~_q and F_{p, <= k} = sum_{q<=k} F~_q.
///
/// So: L^{(2)} = F_2 (since eps_2(2,q) = 0 for all q)
///     L^{(p)} = F_p + (t-1) * F_{p, <= p-2} for p >= 3
///
/// Now F_2 and F_p are sums of source polynomials from DIFFERENT descent sets.
/// But maybe we can relate them through the inductive hypothesis.
///
/// TEST: At what level do the different source descent sets S'_p start to differ?
/// For the simplest case (alternating permutations), check how many distinct
/// source descent sets arise.
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::format_poly;
use std::collections::{BTreeMap, BTreeSet};

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

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        r[i] += c;
    }
    while r.len() > 1 && *r.last().unwrap() == 0 {
        r.pop();
    }
    r
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    println!("=== Direct TP_2 Induction Analysis ===\n");

    // Focus: for each (n, S), decompose L^{(p)} in terms of the modified source vectors
    // F_q^{(p)} and check how the TP_2 minors relate.

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

        // Look at alternating descent set
        let alt_mask: u64 = (1..n)
            .filter(|&i| i % 2 == 0)
            .map(|i| 1u64 << (i - 1))
            .fold(0, |a, b| a | b);
        // Actually just check {2,4,...} intersected with [n-1]

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if !vp.contains(&2) || vp.len() < 2 {
                continue;
            }
            // Only show cases with exactly 2 or 3 valid positions for clarity
            if vp.len() > 3 {
                continue;
            }

            let s_str = descent_set_to_string(mask, n);

            // Collect all source descent sets
            let mut all_source_sets: BTreeSet<u64> = BTreeSet::new();
            for &p in &vp {
                all_source_sets.insert(source_asc(mask, p, n));
                if let Some(sd) = source_desc(mask, p, n) {
                    all_source_sets.insert(sd);
                }
            }

            // For each p, compute the total modified source F_p = sum_q F~_q^{(p)}
            // and the partial sum F_{p, <= p-2}
            let mut f_total: BTreeMap<u8, Vec<i64>> = BTreeMap::new();
            let mut f_partial: BTreeMap<u8, Vec<i64>> = BTreeMap::new();

            for &p in &vp {
                if p <= 1 {
                    continue;
                }

                let sp_a = source_asc(mask, p, n);
                let sp_d = source_desc(mask, p, n);

                let mut by_q: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
                for &sp in &[Some(sp_a), sp_d]
                    .iter()
                    .filter_map(|x| *x)
                    .collect::<Vec<_>>()
                {
                    if let Some(cls) = prev_by_des.get(&sp) {
                        for pi in cls {
                            let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            let e1 = if epsilon1(pi, p) { 1 } else { 0 };
                            by_q.entry(q)
                                .or_default()
                                .push(compute(pi, Stat::Swaps) + e1);
                        }
                    }
                }

                let mut total = vec![0i64];
                let mut partial = vec![0i64]; // sum for q <= p-2

                for (&q, vals) in &by_q {
                    let poly = build_poly(vals);
                    total = poly_add(&total, &poly);
                    if q <= p.saturating_sub(2) {
                        partial = poly_add(&partial, &poly);
                    }
                }

                f_total.insert(p, total);
                f_partial.insert(p, partial);
            }

            // Now L^{(p)} = F_p + (t-1) * F_{p, <= p-2}
            // = F_p + t*F_{p, <= p-2} - F_{p, <= p-2}
            // For p=2: L^{(2)} = F_2 (since F_{2, <= 0} = 0)

            // Check: do all p share the same total F_p?
            let all_same_total = {
                let first = f_total.values().next().unwrap();
                f_total.values().all(|v| v == first)
            };

            if all_same_total {
                println!(
                    "  S={}: ALL F_p EQUAL! F = {}",
                    s_str,
                    format_poly(f_total.values().next().unwrap())
                );
            }

            // Print the decomposition for small cases
            if vp.len() <= 3 && n <= 7 {
                println!("  S={} (|P(S)|={})", s_str, vp.len());
                for &p in &vp {
                    if p <= 1 {
                        continue;
                    }
                    let ft = f_total.get(&p).unwrap();
                    let fp = f_partial.get(&p).unwrap();
                    println!(
                        "    p={}: F_total={}, F_partial(<= p-2)={}",
                        p,
                        format_poly(ft),
                        format_poly(fp)
                    );

                    // Verify: L^{(p)} should equal F_total + (t-1)*F_partial
                    // = F_total - F_partial + t*F_partial
                    let mut check = vec![0i64; ft.len().max(fp.len() + 1)];
                    for (i, &c) in ft.iter().enumerate() {
                        check[i] += c;
                    }
                    for (i, &c) in fp.iter().enumerate() {
                        check[i] -= c;
                    } // subtract F_partial
                    for (i, &c) in fp.iter().enumerate() {
                        // add t*F_partial
                        if i + 1 >= check.len() {
                            check.resize(i + 2, 0);
                        }
                        check[i + 1] += c;
                    }
                    while check.len() > 1 && *check.last().unwrap() == 0 {
                        check.pop();
                    }

                    // Compare with actual L^{(p)}
                    let lp_vals: Vec<usize> = class
                        .iter()
                        .filter(|pi| pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                        .map(|pi| compute(pi, Stat::Swaps))
                        .collect();
                    let lp = build_poly(&lp_vals);
                    if check == lp {
                        println!("      [verified: F + (t-1)*F_partial = L^({})]", p);
                    } else {
                        println!(
                            "      [MISMATCH: computed={}, actual={}]",
                            format_poly(&check),
                            format_poly(&lp)
                        );
                    }
                }
            }
        }
    }

    // KEY TEST: For the case where there are exactly 2 valid positions {2, p_max}:
    // L^{(2)} = F_2
    // L^{(p_max)} = F_{p_max} + (t-1)*F_{p_max, <= p_max-2}
    //
    // The TP_2 condition between these two requires:
    // For all k1 < k2:
    //   [t^k1]L^(2) * [t^k2]L^(p_max) >= [t^k2]L^(2) * [t^k1]L^(p_max)
    //
    // Since L^(p_max) = F_{p_max} + (t-1)*F_{p_max,<=p_max-2},
    // this involves F_2 and F_{p_max} from DIFFERENT source descent sets.
    //
    // Question: is there a relationship between F_2 and F_{p_max}?

    println!("\n\n=== Relationship between F_2 and F_p for two-position cases ===\n");

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

        for (&mask, _) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() != 2 || !vp.contains(&2) {
                continue;
            }

            let s_str = descent_set_to_string(mask, n);
            let p2 = 2u8;
            let p_max = *vp.iter().filter(|&&p| p > 2).next().unwrap();

            // Compute F_2 and F_{p_max}
            let mut f2_vals: Vec<usize> = Vec::new();
            let mut fp_vals: Vec<usize> = Vec::new();

            for &sp in &[Some(source_asc(mask, p2, n)), source_desc(mask, p2, n)]
                .iter()
                .filter_map(|x| *x)
                .collect::<Vec<_>>()
            {
                if let Some(cls) = prev_by_des.get(&sp) {
                    for pi in cls {
                        let e1 = if epsilon1(pi, p2) { 1 } else { 0 };
                        f2_vals.push(compute(pi, Stat::Swaps) + e1);
                    }
                }
            }

            for &sp in &[
                Some(source_asc(mask, p_max, n)),
                source_desc(mask, p_max, n),
            ]
            .iter()
            .filter_map(|x| *x)
            .collect::<Vec<_>>()
            {
                if let Some(cls) = prev_by_des.get(&sp) {
                    for pi in cls {
                        let e1 = if epsilon1(pi, p_max) { 1 } else { 0 };
                        fp_vals.push(compute(pi, Stat::Swaps) + e1);
                    }
                }
            }

            let f2 = build_poly(&f2_vals);
            let fp = build_poly(&fp_vals);

            // Check: is F_2 = F_p?
            if f2 == fp {
                println!(
                    "n={} S={}: F_2 = F_p (p={})  F={}",
                    n,
                    s_str,
                    p_max,
                    format_poly(&f2)
                );
            } else {
                // Check ratio
                let f2_sum: i64 = f2.iter().sum();
                let fp_sum: i64 = fp.iter().sum();
                println!(
                    "n={} S={}: F_2 != F_p (p={})  |F_2|={} |F_p|={}",
                    n, s_str, p_max, f2_sum, fp_sum
                );
                println!("    F_2 = {}", format_poly(&f2));
                println!("    F_p = {}", format_poly(&fp));
            }
        }
    }
}
