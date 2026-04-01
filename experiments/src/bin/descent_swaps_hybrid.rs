/// Hybrid approach: split by eps1 for p=2, modified source for p>=3.
/// Check if the SPLIT source at p=2 forms an interlacing sequence.
/// Check if ALL source columns (from both approaches) jointly interlace.
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

fn check_pairwise_compatible(f: &[i64], g: &[i64]) -> bool {
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
        let mut combo = vec![0i64; maxlen];
        for (i, &c) in f.iter().enumerate() {
            combo[i] += a * c;
        }
        for (i, &c) in g.iter().enumerate() {
            combo[i] += b * c;
        }
        while combo.len() > 1 && *combo.last().unwrap() == 0 {
            combo.pop();
        }
        if combo.len() > 1 && !is_real_rooted(&combo) {
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

    println!("=== Hybrid approach: split at p=2, modified for p>=3 ===\n");

    let mut total_split_ok = 0u64;
    let mut total_split_tested = 0u64;
    let mut total_joint_ok = 0u64;
    let mut total_joint_tested = 0u64;

    for n in 4..=max_n {
        let perms_prev = all_permutations(n - 1);
        let perms = all_permutations(n);

        let mut prev_by_des: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
        let mut prev_list: Vec<&Vec<u8>> = Vec::new();
        for pi in &perms_prev {
            let idx = prev_list.len();
            prev_list.push(pi);
            prev_by_des
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(idx);
        }

        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        let mut n_split_ok = 0u32;
        let mut n_split_tested = 0u32;
        let mut n_joint_ok = 0u32;
        let mut n_joint_tested = 0u32;

        for (&mask, _) in &by_descent {
            let one_not_in_s = (mask & 1) == 0;
            if !one_not_in_s {
                continue;
            } // Only 1 not in S
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            // For p=2: build SPLIT source (A_q, C_q, D_q) separately
            // For p>=3: build MODIFIED source (F_q = A_q + t*C_q + D_q)
            // Collect ALL source polynomials, check joint interlacing

            let mut all_source_polys: Vec<(String, Vec<i64>)> = Vec::new(); // label, poly

            let has_p2 = vp.contains(&2);

            for &p in &vp {
                let sp_a = source_asc(mask, p, n);
                let sp_d = source_desc(mask, p, n);

                // Gather source perms grouped by q
                let mut asc_by_q: BTreeMap<u8, Vec<(usize, bool)>> = BTreeMap::new(); // (swaps, eps1)
                if let Some(indices) = prev_by_des.get(&sp_a) {
                    for &idx in indices {
                        let pi = prev_list[idx];
                        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                        let s = compute(pi, Stat::Swaps);
                        let e1 = epsilon1(pi, p);
                        asc_by_q.entry(q).or_default().push((s, e1));
                    }
                }

                let mut desc_by_q: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
                if let Some(sp) = sp_d {
                    if let Some(indices) = prev_by_des.get(&sp) {
                        for &idx in indices {
                            let pi = prev_list[idx];
                            let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            desc_by_q
                                .entry(q)
                                .or_default()
                                .push(compute(pi, Stat::Swaps));
                        }
                    }
                }

                let all_qs: Vec<u8> = {
                    let mut qs: Vec<u8> =
                        asc_by_q.keys().chain(desc_by_q.keys()).copied().collect();
                    qs.sort();
                    qs.dedup();
                    qs
                };

                if p == 2 {
                    // SPLIT approach: separate A_q, C_q, D_q
                    for &q in &all_qs {
                        // A_q: ascent source, eps1=0
                        let a_vals: Vec<usize> = asc_by_q
                            .get(&q)
                            .map(|v| v.iter().filter(|(_, e)| !e).map(|(s, _)| *s).collect())
                            .unwrap_or_default();
                        if !a_vals.is_empty() {
                            all_source_polys.push((format!("p2_A_q{}", q), build_poly(&a_vals)));
                        }

                        // C_q: ascent source, eps1=1
                        let c_vals: Vec<usize> = asc_by_q
                            .get(&q)
                            .map(|v| v.iter().filter(|(_, e)| *e).map(|(s, _)| *s).collect())
                            .unwrap_or_default();
                        if !c_vals.is_empty() {
                            all_source_polys.push((format!("p2_C_q{}", q), build_poly(&c_vals)));
                        }

                        // D_q: descent source
                        let d_vals: Vec<usize> =
                            desc_by_q.get(&q).map(|v| v.clone()).unwrap_or_default();
                        if !d_vals.is_empty() {
                            all_source_polys.push((format!("p2_D_q{}", q), build_poly(&d_vals)));
                        }
                    }
                } else {
                    // MODIFIED approach: F_q = A_q + t*C_q + D_q (using swaps+eps1)
                    for &q in &all_qs {
                        let mut vals: Vec<usize> = Vec::new();
                        if let Some(v) = asc_by_q.get(&q) {
                            for &(s, e1) in v {
                                vals.push(s + if e1 { 1 } else { 0 });
                            }
                        }
                        if let Some(v) = desc_by_q.get(&q) {
                            for &s in v {
                                vals.push(s);
                            }
                        }
                        if !vals.is_empty() {
                            all_source_polys.push((format!("p{}_F_q{}", p, q), build_poly(&vals)));
                        }
                    }
                }
            }

            // Check: do ALL source polynomials jointly form an interlacing sequence?
            n_joint_tested += 1;
            let polys: Vec<&Vec<i64>> = all_source_polys.iter().map(|(_, p)| p).collect();
            let mut joint_ok = true;
            'outer: for i in 0..polys.len() {
                for j in (i + 1)..polys.len() {
                    if !check_pairwise_compatible(polys[i], polys[j]) {
                        joint_ok = false;
                        break 'outer;
                    }
                }
            }
            if joint_ok {
                n_joint_ok += 1;
            } else {
                let s_str = descent_set_to_string(mask, n);
                println!("  JOINT FAIL: n={} S={} ({} polys)", n, s_str, polys.len());
                // Show which pairs fail
                let mut fail_count = 0;
                for i in 0..polys.len() {
                    for j in (i + 1)..polys.len() {
                        if !check_pairwise_compatible(polys[i], polys[j]) {
                            if fail_count < 3 {
                                println!(
                                    "    {} vs {}: {} vs {}",
                                    all_source_polys[i].0,
                                    all_source_polys[j].0,
                                    format_poly(polys[i]),
                                    format_poly(polys[j])
                                );
                            }
                            fail_count += 1;
                        }
                    }
                }
                if fail_count > 3 {
                    println!("    ... and {} more failures", fail_count - 3);
                }
            }

            // Also check split source at p=2 alone
            if has_p2 {
                let p2_polys: Vec<&Vec<i64>> = all_source_polys
                    .iter()
                    .filter(|(label, _)| label.starts_with("p2_"))
                    .map(|(_, p)| p)
                    .collect();
                if p2_polys.len() >= 2 {
                    n_split_tested += 1;
                    let mut split_ok = true;
                    'outer2: for i in 0..p2_polys.len() {
                        for j in (i + 1)..p2_polys.len() {
                            if !check_pairwise_compatible(p2_polys[i], p2_polys[j]) {
                                split_ok = false;
                                break 'outer2;
                            }
                        }
                    }
                    if split_ok {
                        n_split_ok += 1;
                    } else {
                        let s_str = descent_set_to_string(mask, n);
                        println!("  SPLIT@p2 FAIL: n={} S={}", n, s_str);
                    }
                }
            }
        }

        total_split_ok += n_split_ok as u64;
        total_split_tested += n_split_tested as u64;
        total_joint_ok += n_joint_ok as u64;
        total_joint_tested += n_joint_tested as u64;
        println!(
            "n={}: split@p2 {}/{}, joint {}/{}",
            n, n_split_ok, n_split_tested, n_joint_ok, n_joint_tested
        );
    }

    println!("\n=== SUMMARY ===");
    println!(
        "Split source at p=2: {}/{}",
        total_split_ok, total_split_tested
    );
    println!(
        "Joint interlacing (all sources): {}/{}",
        total_joint_ok, total_joint_tested
    );
}
