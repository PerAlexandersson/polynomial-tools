/// Explore root structure and interlacing patterns for L_{n,S}(t).
///
/// Key experiments:
/// 1. Does L_{n-1,S'}(t) interlace L_{n,S}^{(p)}(t) for each valid p?
///    (Source-to-target interlacing)
/// 2. Do position-refined polynomials have a "proper position" property?
/// 3. Explore alternative refinements (position of min, position of value 1, etc.)
/// 4. Check L_{n,S} interlacing across n for fixed descent pattern
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly, is_real_rooted};
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

fn check_pairwise_compatible(f: &[i64], g: &[i64]) -> bool {
    if f.len() <= 1 || g.len() <= 1 {
        return true;
    }
    if !is_real_rooted(f) || !is_real_rooted(g) {
        return false;
    }
    let test_weights: Vec<(i64, i64)> = vec![
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
    for (a, b) in &test_weights {
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

/// Precompute all data for a given n
struct DescentData {
    /// For each descent set mask, the permutations grouped by swaps value
    by_descent: BTreeMap<u64, Vec<usize>>, // mask -> list of swaps values
    /// Position-of-max refined: (mask, position) -> swaps values
    by_descent_pos: HashMap<(u64, u8), Vec<usize>>,
    /// Position-of-value-1 refined: (mask, position) -> swaps values
    by_descent_pos1: HashMap<(u64, u8), Vec<usize>>,
    /// Position-of-min refined: (mask, position) -> swaps values
    by_descent_posmin: HashMap<(u64, u8), Vec<usize>>,
}

fn compute_descent_data(n: u8) -> DescentData {
    let perms = all_permutations(n);
    let mut by_descent: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
    let mut by_descent_pos: HashMap<(u64, u8), Vec<usize>> = HashMap::new();
    let mut by_descent_pos1: HashMap<(u64, u8), Vec<usize>> = HashMap::new();
    let mut by_descent_posmin: HashMap<(u64, u8), Vec<usize>> = HashMap::new();

    for pi in &perms {
        let mask = descent_set_bitmask(pi);
        let swaps = compute(pi, Stat::Swaps);
        by_descent.entry(mask).or_default().push(swaps);

        // Position of max (value n)
        let pos_max = pi.iter().position(|&v| v == n).unwrap() as u8 + 1;
        by_descent_pos
            .entry((mask, pos_max))
            .or_default()
            .push(swaps);

        // Position of value 1
        let pos_1 = pi.iter().position(|&v| v == 1).unwrap() as u8 + 1;
        by_descent_pos1
            .entry((mask, pos_1))
            .or_default()
            .push(swaps);

        // Position of min (value 1) — same as above, just clearer naming
        let pos_min = pos_1;
        by_descent_posmin
            .entry((mask, pos_min))
            .or_default()
            .push(swaps);
    }

    DescentData {
        by_descent,
        by_descent_pos,
        by_descent_pos1,
        by_descent_posmin,
    }
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    // Precompute data for all n
    let mut data: Vec<Option<DescentData>> = (0..=(max_n as usize)).map(|_| None).collect();
    for n in 2..=max_n {
        data[n as usize] = Some(compute_descent_data(n));
    }

    println!("=== EXPERIMENT 1: Source-to-target interlacing ===");
    println!("Does L_{{n-1,S'}}(t) interlace L_{{n,S}}^{{(p)}}(t)?\n");

    let mut src_interlace_ok = 0u64;
    let mut src_interlace_total = 0u64;

    for n in 4..=max_n {
        let d = data[n as usize].as_ref().unwrap();
        let d_prev = data[(n - 1) as usize].as_ref().unwrap();

        let mut n_ok = 0;
        let mut n_total = 0;

        for (&mask, _) in &d.by_descent {
            let vp = valid_positions(mask, n);
            for &p in &vp {
                let sp_mask = source_descent_set(mask, p, n);
                let target_vals = match d.by_descent_pos.get(&(mask, p)) {
                    Some(v) => v,
                    None => continue,
                };
                let source_vals = match d_prev.by_descent.get(&sp_mask) {
                    Some(v) => v,
                    None => continue,
                };

                let target_poly = build_poly_from_vals(target_vals);
                let source_poly = build_poly_from_vals(source_vals);

                if target_poly.len() <= 1 || source_poly.len() <= 1 {
                    n_ok += 1;
                    n_total += 1;
                    continue;
                }

                n_total += 1;
                let (small, large) = if source_poly.len() <= target_poly.len() {
                    (&source_poly, &target_poly)
                } else {
                    (&target_poly, &source_poly)
                };
                match check_weak_interlacing(small, large) {
                    Some(true) => {
                        n_ok += 1;
                    }
                    _ => {
                        // Also try: does source interlace target?
                        // Or: are they compatible?
                        if check_pairwise_compatible(&source_poly, &target_poly) {
                            n_ok += 1;
                        } else {
                            let s_str = descent_set_to_string(mask, n);
                            println!(
                                "  FAIL n={} S={} p={}: src={} tgt={}",
                                n,
                                s_str,
                                p,
                                format_poly(&source_poly),
                                format_poly(&target_poly)
                            );
                        }
                    }
                }
            }
        }

        src_interlace_ok += n_ok;
        src_interlace_total += n_total;
        println!("n={}: {}/{} source-target compatible", n, n_ok, n_total);
    }
    println!("Total: {}/{}\n", src_interlace_ok, src_interlace_total);

    println!("=== EXPERIMENT 2: Alternative refinement by position of value 1 ===");
    println!("For S with 1 not in S, check if {{L_{{n,S}}^{{pos1=q}}}} is compatible.\n");

    let mut alt_compat_ok = 0u64;
    let mut alt_compat_total = 0u64;

    for n in 4..=max_n {
        let d = data[n as usize].as_ref().unwrap();
        let mut n_ok = 0;
        let mut n_total = 0;

        for (&mask, _) in &d.by_descent {
            // Collect position-of-1 refined polynomials
            let mut pos1_polys: Vec<(u8, Vec<i64>)> = Vec::new();
            for q in 1..=n {
                if let Some(vals) = d.by_descent_pos1.get(&(mask, q)) {
                    if !vals.is_empty() {
                        pos1_polys.push((q, build_poly_from_vals(vals)));
                    }
                }
            }

            if pos1_polys.len() < 2 {
                continue;
            }
            n_total += 1;

            // Check pairwise compatibility
            let mut all_compat = true;
            for i in 0..pos1_polys.len() {
                for j in (i + 1)..pos1_polys.len() {
                    if !check_pairwise_compatible(&pos1_polys[i].1, &pos1_polys[j].1) {
                        all_compat = false;
                    }
                }
            }

            if all_compat {
                n_ok += 1;
            } else {
                let one_not_in_s = (mask & 1) == 0;
                let s_str = descent_set_to_string(mask, n);
                println!("  FAIL n={} S={} (1 not in S: {})", n, s_str, one_not_in_s);
            }
        }

        alt_compat_ok += n_ok;
        alt_compat_total += n_total;
        println!("n={}: {}/{} pos-of-1 compatible", n, n_ok, n_total);
    }
    println!("Total: {}/{}\n", alt_compat_ok, alt_compat_total);

    println!("=== EXPERIMENT 3: Cross-n interlacing for fixed descent patterns ===");
    println!("For k-alternating: does H_n^{{(k)}} interlace H_{{n+1}}^{{(k)}}?\n");

    for k in 2..=5u8 {
        println!("k={}", k);
        let mut prev_poly: Option<Vec<i64>> = None;
        for n in 2..=max_n {
            // k-alternating descent set: descents at k, 2k, 3k, ...
            let mut mask = 0u64;
            let mut j = k;
            while j < n {
                mask |= 1 << (j - 1);
                j += k;
            }

            let d = data[n as usize].as_ref().unwrap();
            let vals = match d.by_descent.get(&mask) {
                Some(v) => v,
                None => {
                    prev_poly = None;
                    continue;
                }
            };
            let poly = build_poly_from_vals(vals);

            if let Some(ref pp) = prev_poly {
                if pp.len() > 1 && poly.len() > 1 {
                    let (small, large) = if pp.len() <= poly.len() {
                        (pp, &poly)
                    } else {
                        (&poly, pp)
                    };
                    match check_weak_interlacing(small, large) {
                        Some(true) => print!("  n={}: interlace OK  ", n),
                        Some(false) => print!("  n={}: FAIL  ", n),
                        None => print!("  n={}: ?  ", n),
                    }
                }
            }
            prev_poly = Some(poly);
        }
        println!();
    }

    println!();
    println!("=== EXPERIMENT 4: What if we refine by (pos_max, pos_min) jointly? ===");
    println!("Check compatible family for joint refinement.\n");

    for n in 4..=max_n.min(8) {
        let d = data[n as usize].as_ref().unwrap();
        let mut n_ok = 0;
        let mut n_total = 0;

        for (&mask, _) in &d.by_descent {
            let one_not_in_s = (mask & 1) == 0;
            if !one_not_in_s {
                continue;
            } // Only test 1 not in S

            // Joint refinement: group by (pos_max, pos_min)
            let perms = all_permutations(n);
            let mut joint: HashMap<(u8, u8), Vec<usize>> = HashMap::new();
            for pi in &perms {
                if descent_set_bitmask(pi) != mask {
                    continue;
                }
                let pos_max = pi.iter().position(|&v| v == n).unwrap() as u8 + 1;
                let pos_min = pi.iter().position(|&v| v == 1).unwrap() as u8 + 1;
                joint
                    .entry((pos_max, pos_min))
                    .or_default()
                    .push(compute(pi, Stat::Swaps));
            }

            let polys: Vec<Vec<i64>> = joint
                .values()
                .filter(|v| !v.is_empty())
                .map(|v| build_poly_from_vals(v))
                .collect();

            if polys.len() < 2 {
                continue;
            }
            n_total += 1;

            let mut all_compat = true;
            for i in 0..polys.len() {
                if !all_compat {
                    break;
                }
                for j in (i + 1)..polys.len() {
                    if !check_pairwise_compatible(&polys[i], &polys[j]) {
                        all_compat = false;
                        break;
                    }
                }
            }

            if all_compat {
                n_ok += 1;
            } else {
                let s_str = descent_set_to_string(mask, n);
                println!("  FAIL n={} S={}", n, s_str);
            }
        }
        println!(
            "n={}: {}/{} joint (pos_max, pos_min) compatible",
            n, n_ok, n_total
        );
    }

    println!();
    println!("=== EXPERIMENT 5: Refine by position of second-largest (n-1) ===");
    println!("Check compatible family for pos(n-1) refinement, 1 not in S.\n");

    for n in 4..=max_n {
        let d = data[n as usize].as_ref().unwrap();
        let perms = all_permutations(n);
        let mut n_ok = 0;
        let mut n_total = 0;

        for (&mask, _) in &d.by_descent {
            let one_not_in_s = (mask & 1) == 0;
            if !one_not_in_s {
                continue;
            }

            // Refine by position of value n-1
            let mut by_pos_nm1: HashMap<u8, Vec<usize>> = HashMap::new();
            for pi in &perms {
                if descent_set_bitmask(pi) != mask {
                    continue;
                }
                let pos = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                by_pos_nm1
                    .entry(pos)
                    .or_default()
                    .push(compute(pi, Stat::Swaps));
            }

            let polys: Vec<Vec<i64>> = by_pos_nm1
                .values()
                .filter(|v| !v.is_empty())
                .map(|v| build_poly_from_vals(v))
                .collect();

            if polys.len() < 2 {
                continue;
            }
            n_total += 1;

            let mut all_compat = true;
            for i in 0..polys.len() {
                if !all_compat {
                    break;
                }
                for j in (i + 1)..polys.len() {
                    if !check_pairwise_compatible(&polys[i], &polys[j]) {
                        all_compat = false;
                        break;
                    }
                }
            }

            if all_compat {
                n_ok += 1;
            } else {
                let s_str = descent_set_to_string(mask, n);
                println!("  FAIL n={} S={}", n, s_str);
            }
        }
        println!(
            "n={}: {}/{} pos(n-1) compatible (1 not in S)",
            n, n_ok, n_total
        );
    }
}
