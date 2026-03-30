/// Verify the algebraic identity for the "staircase" TN_2 minor decomposition.
///
/// Setup: Given polynomials f_1, ..., f_m with coefficient matrix A (A_{j,k} = [t^k]f_j),
/// define the "staircase" polynomials:
///   g_i(t) = t * (f_1 + ... + f_{i-1}) + (f_i + ... + f_m)
///
/// This models: L^{(p)}(t) = sum_{q <= p-2} t * F~_q^{(p)} + sum_{q >= p-1} F~_q^{(p)}
/// when the source polynomials F~_q^{(p)} are the same for all p (indexed by q).
///
/// The coefficient matrix B has B_{i,k} = [t^k]g_i = L_i(k-1) + R_i(k)
/// where L_i(k) = sum_{j<i} A_{j,k} and R_i(k) = sum_{j>=i} A_{j,k}.
///
/// Key algebraic identity (to verify):
///   Delta = B_{i1,k2} * D(k1) - B_{i1,k1} * D(k2)
///
/// where Delta is the 2x2 minor B_{i1,k1}*B_{i2,k2} - B_{i1,k2}*B_{i2,k1}
/// and D(k) = sum_{j=i1}^{i2-1} [A_{j,k} - A_{j,k-1}].
///
/// Tests:
/// 1. For n=5,6,7 and various S, compute A from source polynomials (for fixed p = max P(S))
/// 2. Apply the staircase transform to get B
/// 3. Verify the identity Delta = B_{i1,k2}*D(k1) - B_{i1,k1}*D(k2)
/// 4. Analyze the sign pattern of D(k) and whether Delta >= 0 follows
///
/// Also test: does the staircase preserve TN_2 when A is TN_2?
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

fn is_tn2(rows: &[Vec<i64>]) -> bool {
    let nrows = rows.len();
    if nrows < 2 {
        return true;
    }
    let max_deg = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    for i1 in 0..nrows {
        for i2 in (i1 + 1)..nrows {
            for k1 in 0..max_deg {
                for k2 in (k1 + 1)..max_deg {
                    let minor = coeff(&rows[i1], k1) * coeff(&rows[i2], k2)
                        - coeff(&rows[i1], k2) * coeff(&rows[i2], k1);
                    if minor < 0 {
                        return false;
                    }
                }
            }
        }
    }
    true
}

fn first_neg_minor(rows: &[Vec<i64>]) -> Option<(usize, usize, usize, usize, i64)> {
    let nrows = rows.len();
    if nrows < 2 {
        return None;
    }
    let max_deg = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    for i1 in 0..nrows {
        for i2 in (i1 + 1)..nrows {
            for k1 in 0..max_deg {
                for k2 in (k1 + 1)..max_deg {
                    let minor = coeff(&rows[i1], k1) * coeff(&rows[i2], k2)
                        - coeff(&rows[i1], k2) * coeff(&rows[i2], k1);
                    if minor < 0 {
                        return Some((i1, i2, k1, k2, minor));
                    }
                }
            }
        }
    }
    None
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

/// Apply the staircase transform to a matrix A (given as rows).
/// g_i(t) = t * (f_1 + ... + f_{i-1}) + (f_i + ... + f_m)
/// B_{i,k} = sum_{j<i} A_{j,k-1} + sum_{j>=i} A_{j,k}
fn staircase_transform(a_rows: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let m = a_rows.len();
    let max_deg = a_rows.iter().map(|r| r.len()).max().unwrap_or(0);
    // Need one extra degree since t * f shifts degree up by 1
    let out_deg = max_deg + 1;

    let mut b_rows = Vec::new();
    for i in 0..m {
        let mut row = vec![0i64; out_deg];
        for k in 0..out_deg {
            // sum_{j<i} A_{j,k-1}
            if k >= 1 {
                for j in 0..i {
                    row[k] += coeff(&a_rows[j], k - 1);
                }
            }
            // sum_{j>=i} A_{j,k}
            for j in i..m {
                row[k] += coeff(&a_rows[j], k);
            }
        }
        // Trim trailing zeros
        while row.len() > 1 && *row.last().unwrap() == 0 {
            row.pop();
        }
        b_rows.push(row);
    }
    b_rows
}

/// Compute D(k) = sum_{j=i1}^{i2-1} [A_{j,k} - A_{j,k-1}]
/// for the algebraic identity check.
fn compute_d(a_rows: &[Vec<i64>], i1: usize, i2: usize, k: usize) -> i64 {
    let mut d = 0i64;
    for j in i1..i2 {
        let a_jk = coeff(&a_rows[j], k);
        let a_jk_minus1 = if k >= 1 { coeff(&a_rows[j], k - 1) } else { 0 };
        d += a_jk - a_jk_minus1;
    }
    d
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    println!("=== Staircase TN_2 Identity Verification ===");
    println!("Delta = B_{{i1,k2}} * D(k1) - B_{{i1,k1}} * D(k2)");
    println!("Checking n = 4 to {}\n", max_n);

    let mut identity_checks = 0u64;
    let mut identity_pass = 0u64;
    let mut staircase_tn2_total = 0u64;
    let mut staircase_tn2_pass = 0u64;
    let mut source_tn2_total = 0u64;
    let mut source_tn2_pass = 0u64;
    let mut actual_tn2_total = 0u64;
    let mut actual_tn2_pass = 0u64;

    // Track D(k) sign patterns
    let mut d_monotone_cases = 0u64;
    let mut d_total_cases = 0u64;

    for n in 4..=max_n {
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

        let mut n_id_check = 0u64;
        let mut n_id_pass = 0u64;
        let mut n_staircase_total = 0u32;
        let mut n_staircase_pass = 0u32;
        let mut n_source_total = 0u32;
        let mut n_source_pass = 0u32;
        let mut n_actual_total = 0u32;
        let mut n_actual_pass = 0u32;

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            // ============================================================
            // Part 1: Compute source polynomials F~_q for each valid p
            // The staircase model assumes F~_q is the SAME for all p.
            // In reality, sources may differ for different p.
            // We check both the idealized staircase model and the actual L^{(p)}.
            // ============================================================

            // For each p, compute modified source polys refined by q = pos(n-1)
            let mut source_by_p: BTreeMap<u8, BTreeMap<u8, Vec<i64>>> = BTreeMap::new();

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
                            let q =
                                pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            let e1 = if epsilon1(pi, p) { 1 } else { 0 };
                            by_q.entry(q)
                                .or_default()
                                .push(compute(pi, Stat::Swaps) + e1);
                        }
                    }
                }

                let mut qmap: BTreeMap<u8, Vec<i64>> = BTreeMap::new();
                for (&q, vals) in &by_q {
                    qmap.insert(q, build_poly(vals));
                }
                source_by_p.insert(p, qmap);
            }

            // ============================================================
            // Part 2: Use p = max(P(S)) as the "canonical" source
            // (or try p = vp[last])
            // ============================================================

            let p_max = *vp.last().unwrap();
            if !source_by_p.contains_key(&p_max) || source_by_p[&p_max].is_empty() {
                continue;
            }

            // Get all q values that appear for any p
            let mut all_qs: Vec<u8> = Vec::new();
            for qmap in source_by_p.values() {
                for &q in qmap.keys() {
                    if !all_qs.contains(&q) {
                        all_qs.push(q);
                    }
                }
            }
            all_qs.sort();

            if all_qs.len() < 2 {
                continue;
            }

            // Build source coefficient matrix A from p_max
            let qmap = &source_by_p[&p_max];
            let a_qs: Vec<u8> = all_qs.clone(); // rows indexed by q
            let a_rows: Vec<Vec<i64>> = a_qs
                .iter()
                .map(|q| {
                    qmap.get(q).cloned().unwrap_or_else(|| vec![0])
                })
                .collect();

            // Check TN_2 of source matrix A
            let a_is_tn2 = is_tn2(&a_rows);
            n_source_total += 1;
            if a_is_tn2 {
                n_source_pass += 1;
            }

            // Apply the staircase transform
            let b_rows = staircase_transform(&a_rows);
            let max_b_deg = b_rows.iter().map(|r| r.len()).max().unwrap_or(0);

            // Check TN_2 of staircase output B
            let b_is_tn2 = is_tn2(&b_rows);
            n_staircase_total += 1;
            if b_is_tn2 {
                n_staircase_pass += 1;
            }

            // ============================================================
            // Part 3: Verify the algebraic identity
            // Delta = B_{i1,k2} * D(k1) - B_{i1,k1} * D(k2)
            // for every 2x2 minor (i1 < i2, k1 < k2)
            // ============================================================

            let m = b_rows.len();
            let mut local_id_checks = 0u64;
            let mut local_id_pass = 0u64;

            for i1 in 0..m {
                for i2 in (i1 + 1)..m {
                    for k1 in 0..max_b_deg {
                        for k2 in (k1 + 1)..max_b_deg {
                            // Direct minor computation
                            let delta_direct = coeff(&b_rows[i1], k1)
                                * coeff(&b_rows[i2], k2)
                                - coeff(&b_rows[i1], k2)
                                    * coeff(&b_rows[i2], k1);

                            // Via identity
                            let d_k1 = compute_d(&a_rows, i1, i2, k1);
                            let d_k2 = compute_d(&a_rows, i1, i2, k2);
                            let delta_identity =
                                coeff(&b_rows[i1], k2) * d_k1
                                    - coeff(&b_rows[i1], k1) * d_k2;

                            local_id_checks += 1;
                            if delta_direct == delta_identity {
                                local_id_pass += 1;
                            } else {
                                let s_str = descent_set_to_string(mask, n);
                                println!(
                                    "IDENTITY FAIL: n={} S={} i1={} i2={} k1={} k2={} direct={} identity={}",
                                    n, s_str, i1, i2, k1, k2, delta_direct, delta_identity
                                );
                            }
                        }
                    }
                }
            }

            n_id_check += local_id_checks;
            n_id_pass += local_id_pass;

            // ============================================================
            // Part 4: Check D(k) sign pattern for each (i1, i2)
            // ============================================================

            for i1 in 0..m {
                for i2 in (i1 + 1)..m {
                    let mut d_values: Vec<(usize, i64)> = Vec::new();
                    for k in 0..max_b_deg {
                        let d_k = compute_d(&a_rows, i1, i2, k);
                        if d_k != 0 {
                            d_values.push((k, d_k));
                        }
                    }

                    d_total_cases += 1;

                    // Check if D changes sign
                    let pos_count = d_values.iter().filter(|(_, v)| *v > 0).count();
                    let neg_count = d_values.iter().filter(|(_, v)| *v < 0).count();
                    let is_monotone_sign =
                        pos_count == 0 || neg_count == 0; // single sign

                    if is_monotone_sign {
                        d_monotone_cases += 1;
                    }
                }
            }

            // ============================================================
            // Part 5: Compute the ACTUAL L^{(p)} and check TN_2
            // ============================================================

            let mut actual_rows: Vec<Vec<i64>> = Vec::new();
            let mut actual_labels: Vec<u8> = Vec::new();
            for &p in &vp {
                let vals: Vec<usize> = class
                    .iter()
                    .filter(|pi| {
                        pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p
                    })
                    .map(|pi| compute(pi, Stat::Swaps))
                    .collect();
                if !vals.is_empty() {
                    actual_rows.push(build_poly(&vals));
                    actual_labels.push(p);
                }
            }

            if actual_rows.len() >= 2 {
                n_actual_total += 1;
                if is_tn2(&actual_rows) {
                    n_actual_pass += 1;
                } else {
                    let s_str = descent_set_to_string(mask, n);
                    println!("ACTUAL TN_2 FAIL: n={} S={}", n, s_str);
                }
            }

            // ============================================================
            // Part 6: Detailed output for small cases
            // ============================================================

            if n <= 6 && vp.len() <= 4 {
                let s_str = descent_set_to_string(mask, n);
                println!(
                    "\nn={} S={} P(S)={:?}  source_TN2={} staircase_TN2={}",
                    n, s_str, vp, a_is_tn2, b_is_tn2
                );

                println!("  Source matrix A (rows indexed by q={:?}):", a_qs);
                for (idx, q) in a_qs.iter().enumerate() {
                    println!("    q={}: {}", q, format_poly(&a_rows[idx]));
                }

                println!("  Staircase matrix B (rows indexed by q={:?}):", a_qs);
                for (idx, q) in a_qs.iter().enumerate() {
                    println!("    i={} (q={}): {}", idx, q, format_poly(&b_rows[idx]));
                }

                // Show D(k) for each (i1, i2)
                for i1 in 0..m {
                    for i2 in (i1 + 1)..m {
                        let d_vals: Vec<i64> =
                            (0..max_b_deg).map(|k| compute_d(&a_rows, i1, i2, k)).collect();
                        println!(
                            "  D(k) for i1={},i2={}: {:?}",
                            i1, i2, d_vals
                        );

                        // Show 2x2 minors and identity check
                        for k1 in 0..max_b_deg {
                            for k2 in (k1 + 1)..max_b_deg {
                                let delta = coeff(&b_rows[i1], k1)
                                    * coeff(&b_rows[i2], k2)
                                    - coeff(&b_rows[i1], k2)
                                        * coeff(&b_rows[i2], k1);
                                let d_k1 = compute_d(&a_rows, i1, i2, k1);
                                let d_k2 = compute_d(&a_rows, i1, i2, k2);
                                let b_i1_k1 = coeff(&b_rows[i1], k1);
                                let b_i1_k2 = coeff(&b_rows[i1], k2);

                                if delta != 0 || d_k1 != 0 || d_k2 != 0 {
                                    println!(
                                        "    minor({},{};{},{})={:>4}  B_{{i1,k2}}={:>3} D(k1)={:>3} B_{{i1,k1}}={:>3} D(k2)={:>3}  identity={}",
                                        i1, i2, k1, k2, delta,
                                        b_i1_k2, d_k1, b_i1_k1, d_k2,
                                        b_i1_k2 * d_k1 - b_i1_k1 * d_k2
                                    );
                                }
                            }
                        }
                    }
                }

                // Show actual L^{(p)} for comparison
                if actual_rows.len() >= 2 {
                    println!("  Actual L^{{(p)}} (for comparison):");
                    for (idx, p) in actual_labels.iter().enumerate() {
                        println!("    p={}: {}", p, format_poly(&actual_rows[idx]));
                    }
                }
            }
        }

        identity_checks += n_id_check;
        identity_pass += n_id_pass;
        staircase_tn2_total += n_staircase_total as u64;
        staircase_tn2_pass += n_staircase_pass as u64;
        source_tn2_total += n_source_total as u64;
        source_tn2_pass += n_source_pass as u64;
        actual_tn2_total += n_actual_total as u64;
        actual_tn2_pass += n_actual_pass as u64;

        println!(
            "\nn={}: identity {}/{} | source_TN2 {}/{} | staircase_TN2 {}/{} | actual_TN2 {}/{}",
            n, n_id_pass, n_id_check, n_source_pass, n_source_total,
            n_staircase_pass, n_staircase_total, n_actual_pass, n_actual_total
        );
    }

    println!("\n=== OVERALL SUMMARY (n=4..{}) ===", max_n);
    println!("Identity checks:        {}/{}", identity_pass, identity_checks);
    println!("Source matrix A TN_2:   {}/{}", source_tn2_pass, source_tn2_total);
    println!("Staircase matrix B TN_2:{}/{}", staircase_tn2_pass, staircase_tn2_total);
    println!("Actual L^(p) TN_2:      {}/{}", actual_tn2_pass, actual_tn2_total);
    println!(
        "D(k) single-sign:       {}/{} ({:.1}%)",
        d_monotone_cases,
        d_total_cases,
        if d_total_cases > 0 {
            100.0 * d_monotone_cases as f64 / d_total_cases as f64
        } else {
            0.0
        }
    );

    // ============================================================
    // Part 7: Abstract staircase test with constructed TN_2 matrices
    // ============================================================
    println!("\n=== Abstract Staircase Tests ===");
    println!("Does staircase preserve TN_2 for arbitrary TN_2 input matrices?\n");

    // Test 1: Simple 2x2 TN_2 matrix [a b; c d] with ad >= bc
    // Staircase: g_1 = f_1 + f_2, g_2 = t*f_1 + f_2
    // B = [[a+c, b+d], [c, a+d]]  (wait, let's compute properly)
    // Actually g_1 = (f_1 + f_2) (since i=1, no j<1 terms), g_2 = t*f_1 + f_2
    // B_{1,k} = sum_{j>=1} A_{j,k} = A_{1,k} + A_{2,k}
    // B_{2,k} = A_{1,k-1} + A_{2,k}

    // Test with known TN_2 matrices
    let test_matrices: Vec<(&str, Vec<Vec<i64>>)> = vec![
        (
            "Geometric [1,2,4; 1,3,9]",
            vec![vec![1, 2, 4], vec![1, 3, 9]],
        ),
        (
            "Binomial row [1,3,3,1; 1,4,6,4]",
            vec![vec![1, 3, 3, 1], vec![1, 4, 6, 4]],
        ),
        (
            "3 rows [1,2,1; 1,3,3; 1,4,6]",
            vec![vec![1, 2, 1], vec![1, 3, 3], vec![1, 4, 6]],
        ),
        (
            "4 rows [1,1,0,0; 0,1,1,0; 0,0,1,1; 0,0,0,1]",
            vec![
                vec![1, 1, 0, 0],
                vec![0, 1, 1, 0],
                vec![0, 0, 1, 1],
                vec![0, 0, 0, 1],
            ],
        ),
        (
            "3 rows [1,0; 1,1; 0,1]",
            vec![vec![1, 0], vec![1, 1], vec![0, 1]],
        ),
        (
            "4 rows [3,2,1; 2,2,1; 1,2,2; 1,1,2]",
            vec![vec![3, 2, 1], vec![2, 2, 1], vec![1, 2, 2], vec![1, 1, 2]],
        ),
    ];

    for (name, a_rows) in &test_matrices {
        let a_tn2 = is_tn2(a_rows);
        let b_rows = staircase_transform(a_rows);
        let b_tn2 = is_tn2(&b_rows);

        println!("  {}: A TN_2={} -> B TN_2={}", name, a_tn2, b_tn2);

        if a_tn2 && !b_tn2 {
            println!("    *** COUNTEREXAMPLE: staircase does NOT preserve TN_2! ***");
            println!("    A:");
            for row in a_rows {
                println!("      {:?}", row);
            }
            println!("    B:");
            for row in &b_rows {
                println!("      {:?}", row);
            }
            if let Some((i1, i2, k1, k2, val)) = first_neg_minor(&b_rows) {
                println!(
                    "    First negative minor: rows({},{}) cols({},{}) = {}",
                    i1, i2, k1, k2, val
                );
            }
        }

        // Verify identity for this abstract test
        let m = b_rows.len();
        let max_deg = b_rows.iter().map(|r| r.len()).max().unwrap_or(0);
        let mut id_ok = true;
        for i1 in 0..m {
            for i2 in (i1 + 1)..m {
                for k1 in 0..max_deg {
                    for k2 in (k1 + 1)..max_deg {
                        let delta = coeff(&b_rows[i1], k1) * coeff(&b_rows[i2], k2)
                            - coeff(&b_rows[i1], k2) * coeff(&b_rows[i2], k1);
                        let d_k1 = compute_d(a_rows, i1, i2, k1);
                        let d_k2 = compute_d(a_rows, i1, i2, k2);
                        let identity = coeff(&b_rows[i1], k2) * d_k1
                            - coeff(&b_rows[i1], k1) * d_k2;
                        if delta != identity {
                            id_ok = false;
                            println!(
                                "    ID FAIL: i1={} i2={} k1={} k2={} delta={} identity={}",
                                i1, i2, k1, k2, delta, identity
                            );
                        }
                    }
                }
            }
        }
        if id_ok {
            println!("    Identity verified for all minors.");
        }

        // Show D(k) sign pattern
        for i1 in 0..m {
            for i2 in (i1 + 1)..m {
                let d_vals: Vec<i64> =
                    (0..max_deg).map(|k| compute_d(a_rows, i1, i2, k)).collect();
                let has_pos = d_vals.iter().any(|&v| v > 0);
                let has_neg = d_vals.iter().any(|&v| v < 0);
                let sign_str = if has_pos && has_neg {
                    "MIXED"
                } else if has_pos {
                    "ALL >= 0"
                } else if has_neg {
                    "ALL <= 0"
                } else {
                    "ALL ZERO"
                };
                println!(
                    "    D(k) for i1={},i2={}: {:?}  [{}]",
                    i1, i2, d_vals, sign_str
                );
            }
        }
    }

    // ============================================================
    // Part 8: Deeper analysis - when the identity implies Delta >= 0
    // ============================================================
    println!("\n=== Analysis: When Does Delta >= 0 Follow? ===");
    println!("Delta = B_{{i1,k2}} * D(k1) - B_{{i1,k1}} * D(k2)");
    println!("For Delta >= 0 we need: D(k1)/D(k2) <= B_{{i1,k2}}/B_{{i1,k1}} (when D(k2) > 0)");
    println!("Or equivalently: the ratio D(k)/B_{{i1,k}} is non-increasing in k.\n");

    // Check this ratio condition on all the actual data
    let mut ratio_total = 0u64;
    let mut ratio_pass = 0u64;

    for n in 4..=max_n {
        let perms_prev = all_permutations(n - 1);

        let perms = all_permutations(n);
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

        for (&mask, _class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            let p_max = *vp.last().unwrap();
            if p_max <= 1 {
                continue;
            }

            let sp_a = source_asc(mask, p_max, n);
            let sp_d = source_desc(mask, p_max, n);

            let mut by_q: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
            for &sp in &[Some(sp_a), sp_d]
                .iter()
                .filter_map(|x| *x)
                .collect::<Vec<_>>()
            {
                if let Some(cls) = prev_by_des.get(&sp) {
                    for pi in cls {
                        let q =
                            pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                        let e1 = if epsilon1(pi, p_max) { 1 } else { 0 };
                        by_q.entry(q)
                            .or_default()
                            .push(compute(pi, Stat::Swaps) + e1);
                    }
                }
            }

            let a_qs: Vec<u8> = by_q.keys().copied().collect();
            if a_qs.len() < 2 {
                continue;
            }
            let a_rows: Vec<Vec<i64>> =
                a_qs.iter().map(|q| build_poly(&by_q[q])).collect();

            let b_rows = staircase_transform(&a_rows);
            let m = b_rows.len();
            let max_deg = b_rows.iter().map(|r| r.len()).max().unwrap_or(0);

            for i1 in 0..m {
                for i2 in (i1 + 1)..m {
                    // Check ratio D(k)/B_{i1,k} is non-increasing
                    // Use cross-multiplication: D(k1)*B_{i1,k2} <= D(k2)*B_{i1,k1} (when >=0)
                    // Actually Delta >= 0 means B_{i1,k2}*D(k1) >= B_{i1,k1}*D(k2)
                    // i.e., D(k1)/B_{i1,k1} >= D(k2)/B_{i1,k2} when B > 0
                    // i.e., D(k)/B_{i1,k} is NON-INCREASING.

                    let ratios: Vec<(usize, i64, i64)> = (0..max_deg)
                        .map(|k| (k, compute_d(&a_rows, i1, i2, k), coeff(&b_rows[i1], k)))
                        .filter(|&(_, _, b)| b > 0)
                        .collect();

                    // Check non-increasing: D(k1)*B_{i1,k2} >= D(k2)*B_{i1,k1}
                    let mut ok = true;
                    for w in ratios.windows(2) {
                        let (_, d1, b1) = w[0];
                        let (_, d2, b2) = w[1];
                        // D(k1)/B(k1) >= D(k2)/B(k2) <=> D(k1)*B(k2) >= D(k2)*B(k1)
                        if d1 * b2 < d2 * b1 {
                            ok = false;
                            break;
                        }
                    }
                    ratio_total += 1;
                    if ok {
                        ratio_pass += 1;
                    }
                }
            }
        }
    }

    println!(
        "Ratio D(k)/B_{{i1,k}} non-increasing: {}/{} ({:.1}%)",
        ratio_pass,
        ratio_total,
        if ratio_total > 0 {
            100.0 * ratio_pass as f64 / ratio_total as f64
        } else {
            0.0
        }
    );

    println!("\n=== CONCLUSION ===");
    if identity_pass == identity_checks {
        println!("The algebraic identity Delta = B_{{i1,k2}}*D(k1) - B_{{i1,k1}}*D(k2) is VERIFIED.");
    } else {
        println!(
            "The algebraic identity FAILED in {} cases!",
            identity_checks - identity_pass
        );
    }
    if staircase_tn2_pass == staircase_tn2_total {
        println!("The staircase transform preserves TN_2 in ALL tested cases.");
    } else {
        println!(
            "The staircase transform FAILS to preserve TN_2 in {} cases.",
            staircase_tn2_total - staircase_tn2_pass
        );
    }
    if ratio_pass == ratio_total {
        println!("The ratio D(k)/B_{{i1,k}} is non-increasing in ALL cases => proof strategy viable!");
    } else {
        println!(
            "The ratio condition fails in {} cases => need additional analysis.",
            ratio_total - ratio_pass
        );
    }
}
