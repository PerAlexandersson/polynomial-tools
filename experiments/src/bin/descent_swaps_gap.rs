/// Attack the gap: is each G_r (pos(1)-piece at p=2)
/// compatible with each L^{(p)} for p >= 3?
///
/// If yes, then {G_1,...,G_k, L^{(p)}} is a compatible family
/// (since {G_r} is already compatible), so L^{(2)} = sum G_r
/// is compatible with L^{(p)}.
///
/// Even stronger: is each G_r compatible with each F~_q^{(p)}
/// (the pos(n-1)-piece at p >= 3)?
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{is_real_rooted, format_poly};
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

fn check_compat(f: &[i64], g: &[i64]) -> bool {
    if f.len() <= 1 || g.len() <= 1 { return true; }
    if !is_real_rooted(f) || !is_real_rooted(g) { return false; }
    let weights: Vec<(i64, i64)> = vec![
        (1,1),(1,2),(2,1),(1,3),(3,1),(1,5),(5,1),(2,3),(3,2),(1,7),(7,1),(3,5),(5,3),
    ];
    let maxlen = f.len().max(g.len());
    for (a, b) in &weights {
        let mut c = vec![0i64; maxlen];
        for (i, &v) in f.iter().enumerate() { c[i] += a * v; }
        for (i, &v) in g.iter().enumerate() { c[i] += b * v; }
        while c.len() > 1 && *c.last().unwrap() == 0 { c.pop(); }
        if c.len() > 1 && !is_real_rooted(&c) { return false; }
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

    println!("=== Attacking the gap: G_r vs L^(p) and G_r vs F~_q^(p) ===\n");

    let mut gr_vs_lp_ok = 0u64;
    let mut gr_vs_lp_total = 0u64;
    let mut gr_vs_fq_ok = 0u64;
    let mut gr_vs_fq_total = 0u64;
    let mut all_gr_fq_family_ok = 0u64;
    let mut all_gr_fq_family_total = 0u64;

    for n in 5..=max_n {
        let perms = all_permutations(n);
        let perms_prev = all_permutations(n - 1);

        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms { by_descent.entry(descent_set_bitmask(pi)).or_default().push(pi); }

        let mut prev_by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev { prev_by_des.entry(descent_set_bitmask(pi)).or_default().push(pi); }

        let mut n_gr_lp_ok = 0u32;
        let mut n_gr_lp_total = 0u32;
        let mut n_gr_fq_ok = 0u32;
        let mut n_gr_fq_total = 0u32;
        let mut n_family_ok = 0u32;
        let mut n_family_total = 0u32;

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 || mask & 2 == 0 { continue; }
            let vp = valid_positions(mask, n);
            if !vp.contains(&2) || vp.len() < 2 { continue; }

            // Compute G_r (pos(1)-refined, p=2)
            let sp_a2 = source_asc(mask, 2, n);
            let sp_d2 = source_desc(mask, 2, n);

            let mut g_by_r: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
            for &sp in &[Some(sp_a2), sp_d2].iter().filter_map(|x| *x).collect::<Vec<_>>() {
                if let Some(cls) = prev_by_des.get(&sp) {
                    for pi in cls {
                        let r = pi.iter().position(|&v| v == 1).unwrap() as u8 + 1;
                        let e1 = if epsilon1(pi, 2) { 1 } else { 0 };
                        g_by_r.entry(r).or_default().push(compute(pi, Stat::Swaps) + e1);
                    }
                }
            }
            let g_polys: Vec<(u8, Vec<i64>)> = g_by_r.iter()
                .map(|(&r, vals)| (r, build_poly(vals)))
                .collect();

            // Compute L^{(p)} and F~_q^{(p)} for p >= 3
            let mut lp_polys: Vec<(u8, Vec<i64>)> = Vec::new();
            let mut fq_polys: Vec<(u8, u8, Vec<i64>)> = Vec::new(); // (p, q, poly)

            for &p in &vp {
                if p <= 2 { continue; }

                // L^{(p)}
                let lp_vals: Vec<usize> = class.iter()
                    .filter(|pi| pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                    .map(|pi| compute(pi, Stat::Swaps))
                    .collect();
                if !lp_vals.is_empty() {
                    lp_polys.push((p, build_poly(&lp_vals)));
                }

                // F~_q^{(p)} (modified source, refined by pos(n-1))
                let sp_a = source_asc(mask, p, n);
                let sp_d = source_desc(mask, p, n);

                let mut fq_by_q: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
                for &sp in &[Some(sp_a), sp_d].iter().filter_map(|x| *x).collect::<Vec<_>>() {
                    if let Some(cls) = prev_by_des.get(&sp) {
                        for pi in cls {
                            let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            let e1 = if epsilon1(pi, p) { 1 } else { 0 };
                            fq_by_q.entry(q).or_default().push(compute(pi, Stat::Swaps) + e1);
                        }
                    }
                }
                for (&q, vals) in &fq_by_q {
                    fq_polys.push((p, q, build_poly(vals)));
                }
            }

            // Test 1: G_r vs L^{(p)}
            for (_, g_poly) in &g_polys {
                for (_, lp_poly) in &lp_polys {
                    n_gr_lp_total += 1;
                    if check_compat(g_poly, lp_poly) {
                        n_gr_lp_ok += 1;
                    } else {
                        let s_str = descent_set_to_string(mask, n);
                        println!("  G_r vs L^(p) FAIL: n={} S={}", n, s_str);
                        println!("    G={} L={}", format_poly(g_poly), format_poly(lp_poly));
                    }
                }
            }

            // Test 2: G_r vs F~_q^{(p)}
            for (_, g_poly) in &g_polys {
                for (_, _, fq_poly) in &fq_polys {
                    n_gr_fq_total += 1;
                    if check_compat(g_poly, fq_poly) {
                        n_gr_fq_ok += 1;
                    }
                }
            }

            // Test 3: Full family {G_r, F~_q^{(p)} for all p>=3} compatible?
            n_family_total += 1;
            let mut all_polys: Vec<&Vec<i64>> = Vec::new();
            for (_, poly) in &g_polys { all_polys.push(poly); }
            for (_, _, poly) in &fq_polys { all_polys.push(poly); }

            let mut family_ok = true;
            'outer: for i in 0..all_polys.len() {
                for j in (i+1)..all_polys.len() {
                    if !check_compat(all_polys[i], all_polys[j]) {
                        family_ok = false;
                        break 'outer;
                    }
                }
            }
            if family_ok { n_family_ok += 1; }
        }

        gr_vs_lp_ok += n_gr_lp_ok as u64;
        gr_vs_lp_total += n_gr_lp_total as u64;
        gr_vs_fq_ok += n_gr_fq_ok as u64;
        gr_vs_fq_total += n_gr_fq_total as u64;
        all_gr_fq_family_ok += n_family_ok as u64;
        all_gr_fq_family_total += n_family_total as u64;

        println!("n={}: G_r vs L^(p) {}/{}, G_r vs F~_q {}/{}, full family {}/{}",
            n, n_gr_lp_ok, n_gr_lp_total,
            n_gr_fq_ok, n_gr_fq_total,
            n_family_ok, n_family_total);
    }

    println!("\n=== SUMMARY ===");
    println!("G_r vs L^(p):          {}/{}", gr_vs_lp_ok, gr_vs_lp_total);
    println!("G_r vs F~_q^(p):       {}/{}", gr_vs_fq_ok, gr_vs_fq_total);
    println!("Full {{G_r, F~_q}} family: {}/{}", all_gr_fq_family_ok, all_gr_fq_family_total);
    println!();
    println!("If G_r vs L^(p) = ALL PASS, then:");
    println!("  {{G_r}} compatible (Prop 6.4) + G_r compat w/ L^(p)");
    println!("  => {{G_r, L^(p)}} compatible family");
    println!("  => L^(2) = sum G_r compatible w/ L^(p) (common interlacing thm)");
    println!("  => full {{L^(p)}} compatible family. QED.");
}
