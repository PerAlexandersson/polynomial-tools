/// Deep dive into the p=2 case.
///
/// For p=2, eps2=0 always, so the transfer is flat (all entries 1).
/// The target is L_{n,S}^{(2)} = sum_q F~_q where F~_q = A_q + t*C_q + D_q.
///
/// Questions:
/// 1. What does the p=2 contribution actually look like?
///    Is it always "small" relative to the other L_{n,S}^{(p)}?
/// 2. Can we find a DIFFERENT refinement at p=2 that gives a non-flat matrix?
///    E.g., refine by position of (n-2) instead of (n-1)?
/// 3. Can we refine the p=2 source by the VALUE at position 1 (= pi(1))?
///    Since p-1=1, the eps1 condition is pi(1)+1=pi(2).
/// 4. Does refining by the value pi(1) give interlacing?
/// 5. What if we refine by position of the MINIMUM (value 1)?
/// 6. Can we decompose L_{n,S}^{(2)} = (1+t)*R for some real-rooted R?
///    Or L_{n,S}^{(2)} = R + t*Q with R,Q compatible?
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly, is_real_rooted};
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

fn check_family_compatible(polys: &[Vec<i64>]) -> bool {
    for i in 0..polys.len() {
        for j in (i + 1)..polys.len() {
            if !check_pairwise_compatible(&polys[i], &polys[j]) {
                return false;
            }
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

    println!("=== Deep dive into p=2 case ===\n");

    // For each S with 2 ∈ S and 1 ∉ S, p=2 is valid.
    // At p=2: source perms have n at position 2 in sigma.
    // Source descent sets: S'_2 and S''_2 = S'_2 ∪ {1}.
    // eps2 = 0 always. eps1 = [pi(1)+1 = pi(2)].

    println!("──── Experiment 1: Refine p=2 by value at position 1 ────");
    println!("For pi in the source, let v = pi(1). Refine by v.\n");

    let mut exp1_ok = 0u32;
    let mut exp1_total = 0u32;

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

        let mut n_ok = 0u32;
        let mut n_total = 0u32;

        for (&mask, _) in &by_descent {
            if mask & 1 != 0 {
                continue;
            } // 1 ∈ S
            if mask & 2 == 0 {
                continue;
            } // 2 ∉ S, so p=2 not valid
            let vp = valid_positions(mask, n);
            if !vp.contains(&2) {
                continue;
            }

            // Compute source descent sets for p=2
            // S'_2: positions ≥ 2 of S shifted down by 1
            let sp_a = {
                let mut sp = 0u64;
                for j in 3..n {
                    if (mask >> (j - 1)) & 1 == 1 {
                        sp |= 1 << (j - 2);
                    }
                }
                sp
            };
            let sp_d = sp_a | 1; // S''_2 = S'_2 ∪ {1}

            // Collect all source perms for p=2, refine by v = pi(1)
            let mut by_v: BTreeMap<u8, Vec<usize>> = BTreeMap::new();

            for &sp in &[sp_a, sp_d] {
                if let Some(class) = prev_by_des.get(&sp) {
                    for pi in class {
                        let v = pi[0]; // pi(1) = first element, 1-indexed value
                        let e1 = if pi[0] + 1 == pi[1] { 1 } else { 0 };
                        by_v.entry(v)
                            .or_default()
                            .push(compute(pi, Stat::Swaps) + e1);
                    }
                }
            }

            let mut sorted_vs: Vec<u8> = by_v.keys().copied().collect();
            sorted_vs.sort();

            let v_polys: Vec<Vec<i64>> = sorted_vs.iter().map(|v| build_poly(&by_v[v])).collect();

            if v_polys.len() >= 2 {
                n_total += 1;
                if check_family_compatible(&v_polys) {
                    n_ok += 1;
                } else {
                    let s_str = descent_set_to_string(mask, n);
                    println!("  FAIL n={} S={}", n, s_str);
                }
            }
        }
        exp1_ok += n_ok;
        exp1_total += n_total;
        println!("n={}: refine-by-v {}/{}", n, n_ok, n_total);
    }
    println!("Total refine-by-v: {}/{}\n", exp1_ok, exp1_total);

    println!("──── Experiment 2: Refine p=2 by value at position 2 ────");
    println!("For pi in the source, let w = pi(2). Refine by w.\n");

    let mut exp2_ok = 0u32;
    let mut exp2_total = 0u32;

    for n in 5..=max_n {
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

        let mut n_ok = 0u32;
        let mut n_total = 0u32;

        for (&mask, _) in &by_descent {
            if mask & 1 != 0 || mask & 2 == 0 {
                continue;
            }

            let sp_a = {
                let mut sp = 0u64;
                for j in 3..n {
                    if (mask >> (j - 1)) & 1 == 1 {
                        sp |= 1 << (j - 2);
                    }
                }
                sp
            };
            let sp_d = sp_a | 1;

            let mut by_w: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
            for &sp in &[sp_a, sp_d] {
                if let Some(class) = prev_by_des.get(&sp) {
                    for pi in class {
                        let w = pi[1]; // pi(2)
                        let e1 = if pi[0] + 1 == pi[1] { 1 } else { 0 };
                        by_w.entry(w)
                            .or_default()
                            .push(compute(pi, Stat::Swaps) + e1);
                    }
                }
            }

            let mut sorted_ws: Vec<u8> = by_w.keys().copied().collect();
            sorted_ws.sort();
            let w_polys: Vec<Vec<i64>> = sorted_ws.iter().map(|w| build_poly(&by_w[w])).collect();

            if w_polys.len() >= 2 {
                n_total += 1;
                if check_family_compatible(&w_polys) {
                    n_ok += 1;
                } else {
                    let s_str = descent_set_to_string(mask, n);
                    println!("  FAIL n={} S={}", n, s_str);
                }
            }
        }
        exp2_ok += n_ok;
        exp2_total += n_total;
        println!("n={}: refine-by-w {}/{}", n, n_ok, n_total);
    }
    println!("Total refine-by-w: {}/{}\n", exp2_ok, exp2_total);

    println!("──── Experiment 3: Refine p=2 by (v, eps1) ────");
    println!("Split by value v=pi(1) AND whether eps1=0 or 1.\n");

    let mut exp3_ok = 0u32;
    let mut exp3_total = 0u32;

    for n in 5..=max_n {
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

        let mut n_ok = 0u32;
        let mut n_total = 0u32;

        for (&mask, _) in &by_descent {
            if mask & 1 != 0 || mask & 2 == 0 {
                continue;
            }

            let sp_a = {
                let mut sp = 0u64;
                for j in 3..n {
                    if (mask >> (j - 1)) & 1 == 1 {
                        sp |= 1 << (j - 2);
                    }
                }
                sp
            };
            let sp_d = sp_a | 1;

            // key = (v, eps1_flag)
            let mut by_key: BTreeMap<(u8, bool), Vec<usize>> = BTreeMap::new();
            for &sp in &[sp_a, sp_d] {
                if let Some(class) = prev_by_des.get(&sp) {
                    for pi in class {
                        let v = pi[0];
                        let e1 = pi[0] + 1 == pi[1];
                        // Use swaps (not modified) since we split by eps1
                        by_key
                            .entry((v, e1))
                            .or_default()
                            .push(compute(pi, Stat::Swaps));
                    }
                }
            }

            let mut sorted_keys: Vec<(u8, bool)> = by_key.keys().copied().collect();
            sorted_keys.sort();
            let key_polys: Vec<Vec<i64>> =
                sorted_keys.iter().map(|k| build_poly(&by_key[k])).collect();

            if key_polys.len() >= 2 {
                n_total += 1;
                if check_family_compatible(&key_polys) {
                    n_ok += 1;
                } else {
                    let s_str = descent_set_to_string(mask, n);
                    println!("  FAIL n={} S={}", n, s_str);
                }
            }
        }
        exp3_ok += n_ok;
        exp3_total += n_total;
        println!("n={}: refine-by-(v,eps1) {}/{}", n, n_ok, n_total);
    }
    println!("Total refine-by-(v,eps1): {}/{}\n", exp3_ok, exp3_total);

    println!("──── Experiment 4: Refine p=2 by position of value 1 ────");
    println!("pos(1) = position of the minimum value in pi.\n");

    let mut exp4_ok = 0u32;
    let mut exp4_total = 0u32;

    for n in 5..=max_n {
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

        let mut n_ok = 0u32;
        let mut n_total = 0u32;

        for (&mask, _) in &by_descent {
            if mask & 1 != 0 || mask & 2 == 0 {
                continue;
            }

            let sp_a = {
                let mut sp = 0u64;
                for j in 3..n {
                    if (mask >> (j - 1)) & 1 == 1 {
                        sp |= 1 << (j - 2);
                    }
                }
                sp
            };
            let sp_d = sp_a | 1;

            let mut by_pos1: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
            for &sp in &[sp_a, sp_d] {
                if let Some(class) = prev_by_des.get(&sp) {
                    for pi in class {
                        let pos1 = pi.iter().position(|&v| v == 1).unwrap() as u8 + 1;
                        let e1 = if pi[0] + 1 == pi[1] { 1 } else { 0 };
                        by_pos1
                            .entry(pos1)
                            .or_default()
                            .push(compute(pi, Stat::Swaps) + e1);
                    }
                }
            }

            let mut sorted_ps: Vec<u8> = by_pos1.keys().copied().collect();
            sorted_ps.sort();
            let pos1_polys: Vec<Vec<i64>> =
                sorted_ps.iter().map(|p| build_poly(&by_pos1[p])).collect();

            if pos1_polys.len() >= 2 {
                n_total += 1;
                if check_family_compatible(&pos1_polys) {
                    n_ok += 1;
                } else {
                    let s_str = descent_set_to_string(mask, n);
                    println!("  FAIL n={} S={}", n, s_str);
                }
            }
        }
        exp4_ok += n_ok;
        exp4_total += n_total;
        println!("n={}: refine-by-pos(1) {}/{}", n, n_ok, n_total);
    }
    println!("Total refine-by-pos(1): {}/{}\n", exp4_ok, exp4_total);

    println!("──── Experiment 5: Refine p=2 by position of n-2 ────");
    println!("Instead of pos(n-1), try pos(n-2) as the refinement.\n");

    let mut exp5_ok = 0u32;
    let mut exp5_total = 0u32;

    for n in 5..=max_n {
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

        let mut n_ok = 0u32;
        let mut n_total = 0u32;

        for (&mask, _) in &by_descent {
            if mask & 1 != 0 || mask & 2 == 0 {
                continue;
            }
            if n < 4 {
                continue;
            }

            let sp_a = {
                let mut sp = 0u64;
                for j in 3..n {
                    if (mask >> (j - 1)) & 1 == 1 {
                        sp |= 1 << (j - 2);
                    }
                }
                sp
            };
            let sp_d = sp_a | 1;

            let mut by_pos_nm2: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
            for &sp in &[sp_a, sp_d] {
                if let Some(class) = prev_by_des.get(&sp) {
                    for pi in class {
                        let pos = pi.iter().position(|&v| v == n - 2).unwrap() as u8 + 1;
                        let e1 = if pi[0] + 1 == pi[1] { 1 } else { 0 };
                        by_pos_nm2
                            .entry(pos)
                            .or_default()
                            .push(compute(pi, Stat::Swaps) + e1);
                    }
                }
            }

            let mut sorted_ps: Vec<u8> = by_pos_nm2.keys().copied().collect();
            sorted_ps.sort();
            let polys: Vec<Vec<i64>> = sorted_ps
                .iter()
                .map(|p| build_poly(&by_pos_nm2[p]))
                .collect();

            if polys.len() >= 2 {
                n_total += 1;
                if check_family_compatible(&polys) {
                    n_ok += 1;
                } else {
                    let s_str = descent_set_to_string(mask, n);
                    println!("  FAIL n={} S={}", n, s_str);
                }
            }
        }
        exp5_ok += n_ok;
        exp5_total += n_total;
        println!("n={}: refine-by-pos(n-2) {}/{}", n, n_ok, n_total);
    }
    println!("Total refine-by-pos(n-2): {}/{}\n", exp5_ok, exp5_total);

    println!("──── Experiment 6: Check if L^(2) is compatible with all L^(p>=3) ────");
    println!(
        "Even if the p=2 source isn't interlacing, maybe L^(2) is still\n\
              compatible with the other L^(p) pieces.\n"
    );

    let mut exp6_ok = 0u32;
    let mut exp6_total = 0u32;

    for n in 5..=max_n {
        let perms = all_permutations(n);
        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        let mut n_ok = 0u32;
        let mut n_total = 0u32;

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 || mask & 2 == 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            // Build L^(p) for each valid p
            let mut pos_polys: Vec<(u8, Vec<i64>)> = Vec::new();
            for &p in &vp {
                let vals: Vec<usize> = class
                    .iter()
                    .filter(|pi| pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                    .map(|pi| compute(pi, Stat::Swaps))
                    .collect();
                if !vals.is_empty() {
                    pos_polys.push((p, build_poly(&vals)));
                }
            }

            // Check: is L^(2) compatible with each L^(p) for p >= 3?
            let l2 = pos_polys.iter().find(|(p, _)| *p == 2);
            if l2.is_none() {
                continue;
            }
            let l2_poly = &l2.unwrap().1;

            n_total += 1;
            let mut ok = true;
            for (p, poly) in &pos_polys {
                if *p == 2 {
                    continue;
                }
                if !check_pairwise_compatible(l2_poly, poly) {
                    ok = false;
                    let s_str = descent_set_to_string(mask, n);
                    println!(
                        "  FAIL n={} S={}: L^(2)={} not compat with L^({})={}",
                        n,
                        s_str,
                        format_poly(l2_poly),
                        p,
                        format_poly(poly)
                    );
                    break;
                }
            }
            if ok {
                n_ok += 1;
            }
        }
        exp6_ok += n_ok;
        exp6_total += n_total;
        println!("n={}: L^(2) compat with L^(p>=3): {}/{}", n, n_ok, n_total);
    }
    println!("Total L^(2) compat: {}/{}\n", exp6_ok, exp6_total);

    println!("──── Experiment 7: Matrix at p=2 with entries {{0,1,t}} ────");
    println!("Since eps2=0 at p=2, the eps1 split gives entries 1 (eps1=0) and t (eps1=1).");
    println!("The matrix is [[1,t]] (single row). This trivially satisfies Brändén.");
    println!("Question: does the split source (A_q, C_q separately) form interlacing?\n");

    let mut exp7_ok = 0u32;
    let mut exp7_total = 0u32;

    for n in 5..=max_n {
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

        let mut n_ok = 0u32;
        let mut n_total = 0u32;

        for (&mask, _) in &by_descent {
            if mask & 1 != 0 || mask & 2 == 0 {
                continue;
            }

            let sp_a = {
                let mut sp = 0u64;
                for j in 3..n {
                    if (mask >> (j - 1)) & 1 == 1 {
                        sp |= 1 << (j - 2);
                    }
                }
                sp
            };
            let sp_d = sp_a | 1;

            // Split ALL source perms by (q, eps1_flag, source_type)
            // For Brändén at p=2: A pieces get multiplier 1, C pieces get multiplier t, D pieces get 1
            // All entries in {0, 1, t} ✓
            // Need: the source vector (all A_q, C_q, D_q) forms interlacing
            let mut pieces: Vec<(String, Vec<i64>)> = Vec::new();

            for &(sp, label) in &[(sp_a, "asc"), (sp_d, "desc")] {
                if let Some(class) = prev_by_des.get(&sp) {
                    // Group by q = pos(n-1)
                    let mut by_q: BTreeMap<u8, (Vec<usize>, Vec<usize>)> = BTreeMap::new();
                    for pi in class {
                        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                        let s = compute(pi, Stat::Swaps);
                        let e1 = pi[0] + 1 == pi[1];
                        let entry = by_q.entry(q).or_insert((vec![], vec![]));
                        if label == "desc" || !e1 {
                            entry.0.push(s); // eps1=0 part
                        } else {
                            entry.1.push(s); // eps1=1 part
                        }
                    }

                    for (q, (a_vals, c_vals)) in &by_q {
                        if !a_vals.is_empty() {
                            pieces.push((format!("{}_A_q{}", label, q), build_poly(a_vals)));
                        }
                        if !c_vals.is_empty() {
                            pieces.push((format!("{}_C_q{}", label, q), build_poly(c_vals)));
                        }
                    }
                }
            }

            let polys: Vec<Vec<i64>> = pieces.iter().map(|(_, p)| p.clone()).collect();
            if polys.len() >= 2 {
                n_total += 1;
                if check_family_compatible(&polys) {
                    n_ok += 1;
                } else {
                    let s_str = descent_set_to_string(mask, n);
                    println!("  FAIL n={} S={} ({} pieces)", n, s_str, polys.len());
                }
            }
        }
        exp7_ok += n_ok;
        exp7_total += n_total;
        println!("n={}: split (A,C,D) interlacing {}/{}", n, n_ok, n_total);
    }
    println!("Total split interlacing: {}/{}\n", exp7_ok, exp7_total);

    println!("=== SUMMARY ===");
    println!("Exp 1 (refine by v=pi(1)):     {}/{}", exp1_ok, exp1_total);
    println!("Exp 2 (refine by w=pi(2)):     {}/{}", exp2_ok, exp2_total);
    println!("Exp 3 (refine by (v, eps1)):   {}/{}", exp3_ok, exp3_total);
    println!("Exp 4 (refine by pos(1)):      {}/{}", exp4_ok, exp4_total);
    println!("Exp 5 (refine by pos(n-2)):    {}/{}", exp5_ok, exp5_total);
    println!("Exp 6 (L^(2) compat w/ rest): {}/{}", exp6_ok, exp6_total);
    println!("Exp 7 (split A,C,D interlac): {}/{}", exp7_ok, exp7_total);
}
