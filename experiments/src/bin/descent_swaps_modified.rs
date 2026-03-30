/// Check if the "modified source" approach gives a pure staircase transfer.
///
/// For each target S and valid position p, define the modified source:
///   F_q(t) = sum_{pi: source for p, pi^{-1}(n-1) = q} t^{swaps(pi) + eps1(pi, p)}
///
/// Then L_{n,S}^{(p)} = sum_q t^{eps2(p,q)} F_q
/// which is a pure staircase (entries 1 and t) on the vector (F_q).
///
/// We check: does (F_q)_q form an interlacing sequence for each p?
/// If yes, by Brändén's Thm 7.8.5, the staircase preserves it.
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{format_poly, is_real_rooted};
use std::collections::BTreeMap;

fn build_poly(vals: &[usize]) -> Vec<i64> {
    if vals.is_empty() { return vec![0]; }
    let max_s = *vals.iter().max().unwrap();
    let mut coeffs = vec![0i64; max_s + 1];
    for &s in vals { coeffs[s] += 1; }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 { coeffs.pop(); }
    coeffs
}

fn valid_positions(s_mask: u64, n: u8) -> Vec<u8> {
    let mut positions = Vec::new();
    for p in 1..n {
        if (s_mask >> (p - 1)) & 1 == 1 && (p < 2 || (s_mask >> (p - 2)) & 1 == 0) {
            positions.push(p);
        }
    }
    if n >= 2 && (s_mask >> (n - 2)) & 1 == 0 { positions.push(n); }
    positions
}

fn source_asc(s_mask: u64, p: u8, n: u8) -> u64 {
    if n <= 2 { return 0; }
    if p == n { return s_mask; }
    let mut sp = 0u64;
    if p == 1 {
        for j in 2..n { if (s_mask >> (j - 1)) & 1 == 1 { sp |= 1 << (j - 2); } }
    } else {
        for pos in 1..=(p.saturating_sub(2)) {
            if (s_mask >> (pos - 1)) & 1 == 1 { sp |= 1 << (pos - 1); }
        }
        for j in (p + 1)..n { if (s_mask >> (j - 1)) & 1 == 1 { sp |= 1 << (j - 2); } }
    }
    sp
}

fn source_desc(s_mask: u64, p: u8, n: u8) -> Option<u64> {
    if p <= 1 || p >= n { return None; }
    Some(source_asc(s_mask, p, n) | (1 << (p - 2)))
}

fn epsilon1(pi: &[u8], p: u8) -> bool {
    let n = pi.len() as u8 + 1;
    if p <= 1 || p >= n { return false; }
    pi[(p - 2) as usize] + 1 == pi[(p - 1) as usize]
}

fn check_pairwise_compatible(f: &[i64], g: &[i64]) -> bool {
    if f.len() <= 1 || g.len() <= 1 { return true; }
    if !is_real_rooted(f) || !is_real_rooted(g) { return false; }
    let weights: Vec<(i64, i64)> = vec![
        (1,1),(1,2),(2,1),(1,3),(3,1),(1,5),(5,1),(2,3),(3,2),(1,7),(7,1),(3,5),(5,3),
    ];
    let maxlen = f.len().max(g.len());
    for (a, b) in &weights {
        let mut combo = vec![0i64; maxlen];
        for (i, &c) in f.iter().enumerate() { combo[i] += a * c; }
        for (i, &c) in g.iter().enumerate() { combo[i] += b * c; }
        while combo.len() > 1 && *combo.last().unwrap() == 0 { combo.pop(); }
        if combo.len() > 1 && !is_real_rooted(&combo) { return false; }
    }
    true
}

fn descent_set_to_string(mask: u64, n: u8) -> String {
    let mut s = String::from("{");
    let mut first = true;
    for i in 1..n { if (mask >> (i - 1)) & 1 == 1 { if !first { s.push(','); } s.push_str(&i.to_string()); first = false; } }
    s.push('}'); s
}

fn main() {
    let max_n: u8 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(9);

    println!("=== Modified source approach: absorb eps1 into source statistic ===\n");
    println!("For each target (S, p), compute F_q = sum t^{{swaps+eps1}} over source perms at q.");
    println!("Check if (F_q)_q is interlacing.\n");

    let mut total_ok = 0u64;
    let mut total_tested = 0u64;
    let mut total_target_ok = 0u64;
    let mut total_target_tested = 0u64;

    for n in 4..=max_n {
        let perms_prev = all_permutations(n - 1);
        let perms = all_permutations(n);

        let mut by_descent_prev: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
        let mut perm_prev_list: Vec<&Vec<u8>> = Vec::new();
        for pi in &perms_prev {
            let idx = perm_prev_list.len();
            perm_prev_list.push(pi);
            by_descent_prev.entry(descent_set_bitmask(pi)).or_default().push(idx);
        }

        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms { by_descent.entry(descent_set_bitmask(pi)).or_default().push(pi); }

        let mut n_src_ok = 0u32;
        let mut n_src_tested = 0u32;
        let mut n_tgt_ok = 0u32;
        let mut n_tgt_tested = 0u32;

        for (&mask, _) in &by_descent {
            let one_not_in_s = (mask & 1) == 0;
            let vp = valid_positions(mask, n);
            if vp.len() < 2 { continue; }

            // For each p, compute modified source vector and check interlacing
            let mut target_polys: Vec<(u8, Vec<i64>)> = Vec::new();

            for &p in &vp {
                let sp_a = source_asc(mask, p, n);
                let sp_d = source_desc(mask, p, n);

                // Collect all source perms for this p, grouped by q = pos(n-1)
                let mut by_q: BTreeMap<u8, Vec<usize>> = BTreeMap::new();

                // Ascent source
                if let Some(indices) = by_descent_prev.get(&sp_a) {
                    for &idx in indices {
                        let pi = perm_prev_list[idx];
                        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                        let s = compute(pi, Stat::Swaps) + if epsilon1(pi, p) { 1 } else { 0 };
                        by_q.entry(q).or_default().push(s);
                    }
                }

                // Descent source
                if let Some(sp) = sp_d {
                    if let Some(indices) = by_descent_prev.get(&sp) {
                        for &idx in indices {
                            let pi = perm_prev_list[idx];
                            let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            let s = compute(pi, Stat::Swaps); // eps1=0 always
                            by_q.entry(q).or_default().push(s);
                        }
                    }
                }

                // Build modified source polys, ordered by q
                let mut sorted_qs: Vec<u8> = by_q.keys().copied().collect();
                sorted_qs.sort();

                let f_polys: Vec<Vec<i64>> = sorted_qs.iter()
                    .map(|q| build_poly(&by_q[q]))
                    .collect();

                // Check: is (F_q)_q interlacing?
                if f_polys.len() >= 2 {
                    n_src_tested += 1;
                    let mut all_compat = true;
                    'outer: for i in 0..f_polys.len() {
                        for j in (i+1)..f_polys.len() {
                            if !check_pairwise_compatible(&f_polys[i], &f_polys[j]) {
                                all_compat = false;
                                break 'outer;
                            }
                        }
                    }
                    if all_compat {
                        n_src_ok += 1;
                    } else {
                        if one_not_in_s {
                            let s_str = descent_set_to_string(mask, n);
                            println!("  SRC FAIL (1 not in S): n={} S={} p={}", n, s_str, p);
                            for (qi, q) in sorted_qs.iter().enumerate() {
                                println!("    q={}: {}", q, format_poly(&f_polys[qi]));
                            }
                        }
                    }
                }

                // Also verify the staircase reconstruction
                let mut reconstructed = vec![0i64; 1];
                for (qi, &q) in sorted_qs.iter().enumerate() {
                    let eps2 = if q <= p.saturating_sub(2) { 1usize } else { 0 };
                    for (k, &c) in f_polys[qi].iter().enumerate() {
                        let idx = k + eps2;
                        if idx >= reconstructed.len() { reconstructed.resize(idx + 1, 0); }
                        reconstructed[idx] += c;
                    }
                }
                while reconstructed.len() > 1 && *reconstructed.last().unwrap() == 0 { reconstructed.pop(); }
                target_polys.push((p, reconstructed));
            }

            // Check target interlacing
            if target_polys.len() >= 2 {
                n_tgt_tested += 1;
                let mut all_compat = true;
                for i in 0..target_polys.len() {
                    for j in (i+1)..target_polys.len() {
                        if !check_pairwise_compatible(&target_polys[i].1, &target_polys[j].1) {
                            all_compat = false;
                        }
                    }
                }
                if all_compat { n_tgt_ok += 1; }
                else if one_not_in_s {
                    let s_str = descent_set_to_string(mask, n);
                    println!("  TGT FAIL (1 not in S): n={} S={}", n, s_str);
                }
            }
        }

        total_ok += n_src_ok as u64;
        total_tested += n_src_tested as u64;
        total_target_ok += n_tgt_ok as u64;
        total_target_tested += n_tgt_tested as u64;

        println!("n={}: modified source interlacing {}/{}, target interlacing {}/{}",
            n, n_src_ok, n_src_tested, n_tgt_ok, n_tgt_tested);
    }

    println!("\n=== SUMMARY ===");
    println!("Modified source interlacing: {}/{}", total_ok, total_tested);
    println!("Target interlacing: {}/{}", total_target_ok, total_target_tested);
    println!("\nIf modified source is always interlacing (for 1 not in S),");
    println!("then Branden's Thm 7.8.5 on the staircase matrix proves the conjecture.");
}
