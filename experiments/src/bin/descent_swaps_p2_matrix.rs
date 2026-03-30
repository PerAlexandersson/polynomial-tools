/// Investigate the matrix structure of the p=2 case when refined by pos(1).
///
/// At p=2 with eps2=0 always:
/// - eps1 depends on pi(1)+1 = pi(2)
/// - swaps~(pi) = swaps(pi) + eps1
///
/// If we refine by r = pos(value 1 in pi), what does the "matrix" look like?
/// L_{n,S}^{(2)} = sum_r M_{2,r} * G_r
/// where G_r = sum_{pi: source, pos(1)=r} t^{swaps~(pi)}
///
/// What is M_{2,r}? Since eps2=0, the entry is just 1 for all r.
/// So L^{(2)} = sum_r G_r. This is flat.
///
/// But maybe there's a STAIRCASE if we look at it differently?
/// What if the modified statistic swaps~ introduces a natural ordering
/// where the refinement by pos(1) has a staircase-like structure?
///
/// Let's also check: for p>=3, is the pos(1)-refined source also interlacing?
/// If so, maybe pos(1) works as a UNIVERSAL refinement.
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

fn check_family(polys: &[Vec<i64>]) -> bool {
    for i in 0..polys.len() {
        for j in (i+1)..polys.len() {
            if !check_pairwise_compatible(&polys[i], &polys[j]) { return false; }
        }
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

    println!("=== p=2 matrix structure with pos(1) refinement ===\n");

    println!("──── Check: joint refinement by (pos(n-1) for p>=3, pos(1) for p=2) ────");
    println!("For each S with 2 in S and 1 not in S:");
    println!("  For p>=3: modified source refined by q=pos(n-1) -> interlacing");
    println!("  For p=2: modified source refined by r=pos(1) -> interlacing");
    println!("  ALL these polynomials together: compatible family?\n");

    let mut total_ok = 0u32;
    let mut total_tested = 0u32;

    for n in 5..=max_n {
        let perms = all_permutations(n);
        let perms_prev = all_permutations(n - 1);

        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms { by_descent.entry(descent_set_bitmask(pi)).or_default().push(pi); }

        let mut prev_by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev { prev_by_des.entry(descent_set_bitmask(pi)).or_default().push(pi); }

        let mut n_ok = 0u32;
        let mut n_total = 0u32;

        for (&mask, _) in &by_descent {
            if mask & 1 != 0 { continue; } // 1 ∈ S
            let vp = valid_positions(mask, n);
            if vp.len() < 2 { continue; }
            let has_p2 = vp.contains(&2);

            let mut all_polys: Vec<(String, Vec<i64>)> = Vec::new();

            for &p in &vp {
                let sp_a = source_asc(mask, p, n);
                let sp_d = source_desc(mask, p, n);

                // Collect source perms
                let mut source_perms: Vec<(&Vec<u8>, bool)> = Vec::new(); // (pi, is_desc_source)
                if let Some(class) = prev_by_des.get(&sp_a) {
                    for pi in class { source_perms.push((pi, false)); }
                }
                if let Some(sp) = sp_d {
                    if let Some(class) = prev_by_des.get(&sp) {
                        for pi in class { source_perms.push((pi, true)); }
                    }
                }

                if p == 2 && has_p2 {
                    // Refine by r = pos(value 1)
                    let mut by_r: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
                    for &(pi, _is_desc) in &source_perms {
                        let r = pi.iter().position(|&v| v == 1).unwrap() as u8 + 1;
                        let e1 = if epsilon1(pi, p) { 1 } else { 0 };
                        by_r.entry(r).or_default().push(compute(pi, Stat::Swaps) + e1);
                    }
                    for (r, vals) in &by_r {
                        all_polys.push((format!("p2_r{}", r), build_poly(vals)));
                    }
                } else {
                    // Refine by q = pos(n-1)
                    let mut by_q: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
                    for &(pi, _) in &source_perms {
                        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                        let e1 = if epsilon1(pi, p) { 1 } else { 0 };
                        by_q.entry(q).or_default().push(compute(pi, Stat::Swaps) + e1);
                    }
                    for (q, vals) in &by_q {
                        all_polys.push((format!("p{}_q{}", p, q), build_poly(vals)));
                    }
                }
            }

            let polys: Vec<Vec<i64>> = all_polys.iter().map(|(_, p)| p.clone()).collect();
            if polys.len() >= 2 {
                n_total += 1;
                if check_family(&polys) {
                    n_ok += 1;
                } else {
                    let s_str = descent_set_to_string(mask, n);
                    println!("  FAIL n={} S={} ({} polys)", n, s_str, polys.len());
                    // Find which pairs fail
                    let mut fail_count = 0;
                    for i in 0..polys.len() {
                        for j in (i+1)..polys.len() {
                            if !check_pairwise_compatible(&polys[i], &polys[j]) {
                                if fail_count < 2 {
                                    println!("    {} vs {}", all_polys[i].0, all_polys[j].0);
                                    println!("    {} vs {}", format_poly(&polys[i]), format_poly(&polys[j]));
                                }
                                fail_count += 1;
                            }
                        }
                    }
                    if fail_count > 2 { println!("    ... {} more", fail_count - 2); }
                }
            }
        }
        total_ok += n_ok;
        total_tested += n_total;
        println!("n={}: joint compatible {}/{}", n, n_ok, n_total);
    }
    println!("\nTotal joint: {}/{}", total_ok, total_tested);
}
