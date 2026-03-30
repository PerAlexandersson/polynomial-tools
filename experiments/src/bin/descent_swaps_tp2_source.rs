/// Check TP_2 of the MODIFIED SOURCE coefficient matrix.
///
/// For each target (n, S) and insertion position p:
///   The modified source F~_q^{(p)} = sum_{pi in source(p), pi^{-1}(n-1)=q} t^{swaps(pi)+eps1(pi,p)}
///   is indexed by q (position of n-1 in pi).
///
/// Form the coefficient matrix with rows indexed by q, columns by degree.
/// Check TP_2.
///
/// If TP_2 holds for the source AND the staircase preserves TP_2,
/// then TP_2 of the target follows.
///
/// Also check: does TP_2 hold for the COMBINED source across all p?
/// I.e., stack the source matrices for different p values.
///
/// KEY NEW TEST: Check TP_2 of the "universal staircase output" matrix.
/// The staircase maps: source row q -> output row p via t^{eps_2(p,q)}.
/// If we compute the coefficient matrix of L^{(p)} as M_{p,k} = [t^k] L^{(p)},
/// and the source matrix as N_{q,k} = [t^k] F~_q,
/// then M = Staircase * N (in polynomial multiplication sense).
/// The staircase itself is TP (it has 0/1 entries in staircase pattern).
/// Claim: if N is TP_2, then M is TP_2 (since TP * TP = TP for real matrices,
/// but here the "multiplication" is polynomial, not real).
///
/// Actually, the staircase acts as: L^{(p)} = sum_q t^{eps_2(p,q)} F~_q
/// In coefficient terms: M_{p,k} = sum_q [t^k] (t^{eps_2(p,q)} F~_q)
///                                = sum_q [t^{k-eps_2(p,q)}] F~_q
///                                = sum_q N_{q, k-eps_2(p,q)}
/// This is a "shifted" multiplication, not a standard matrix product.
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::format_poly;
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

fn descent_set_to_string(mask: u64, n: u8) -> String {
    let mut s = String::from("{");
    let mut first = true;
    for i in 1..n { if (mask >> (i - 1)) & 1 == 1 { if !first { s.push(','); } s.push_str(&i.to_string()); first = false; } }
    s.push('}'); s
}

fn coeff(poly: &[i64], k: usize) -> i64 {
    if k < poly.len() { poly[k] } else { 0 }
}

fn check_tp2(rows: &[Vec<i64>]) -> bool {
    let max_deg = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    for i in 0..rows.len() {
        for j in (i+1)..rows.len() {
            for k1 in 0..max_deg {
                for k2 in (k1+1)..max_deg {
                    let minor = coeff(&rows[i], k1) * coeff(&rows[j], k2)
                              - coeff(&rows[i], k2) * coeff(&rows[j], k1);
                    if minor < 0 { return false; }
                }
            }
        }
    }
    true
}

fn main() {
    let max_n: u8 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(9);

    println!("=== TP_2 of source coefficient matrices ===\n");

    let mut src_tp2_ok = 0u64;
    let mut src_tp2_total = 0u64;
    let mut p2_src_tp2_ok = 0u64;
    let mut p2_src_tp2_total = 0u64;

    // Also: check TP_2 of the pos(1)-refined source at p=2
    let mut p2_pos1_tp2_ok = 0u64;
    let mut p2_pos1_tp2_total = 0u64;

    for n in 5..=max_n {
        let perms_prev = all_permutations(n - 1);
        let perms = all_permutations(n);

        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms { by_descent.entry(descent_set_bitmask(pi)).or_default().push(pi); }

        let mut prev_by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev { prev_by_des.entry(descent_set_bitmask(pi)).or_default().push(pi); }

        let mut n_src_ok = 0u32;
        let mut n_src_total = 0u32;
        let mut n_p2_src_ok = 0u32;
        let mut n_p2_src_total = 0u32;
        let mut n_p2_pos1_ok = 0u32;
        let mut n_p2_pos1_total = 0u32;

        for (&mask, _class) in &by_descent {
            if mask & 1 != 0 { continue; }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 { continue; }

            for &p in &vp {
                if p <= 1 { continue; }

                let sp_a = source_asc(mask, p, n);
                let sp_d = source_desc(mask, p, n);

                // Collect source permutations
                let mut source_perms: Vec<(&Vec<u8>, usize)> = Vec::new(); // (perm, modified_swaps)
                for &sp in &[Some(sp_a), sp_d].iter().filter_map(|x| *x).collect::<Vec<_>>() {
                    if let Some(cls) = prev_by_des.get(&sp) {
                        for pi in cls {
                            let e1 = if epsilon1(pi, p) { 1 } else { 0 };
                            let ms = compute(pi, Stat::Swaps) + e1;
                            source_perms.push((pi, ms));
                        }
                    }
                }

                if source_perms.is_empty() { continue; }

                if p >= 3 {
                    // Refine by pos(n-1)
                    let mut by_q: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
                    for (pi, ms) in &source_perms {
                        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                        by_q.entry(q).or_default().push(*ms);
                    }

                    let mut source_rows: Vec<Vec<i64>> = Vec::new();
                    let mut keys: Vec<u8> = by_q.keys().copied().collect();
                    keys.sort();
                    for q in &keys {
                        source_rows.push(build_poly(by_q.get(q).unwrap()));
                    }

                    if source_rows.len() >= 2 {
                        n_src_total += 1;
                        if check_tp2(&source_rows) {
                            n_src_ok += 1;
                        } else {
                            let s_str = descent_set_to_string(mask, n);
                            println!("  Source TP_2 FAIL: n={} S={} p={}", n, s_str, p);
                        }
                    }
                } else if p == 2 {
                    // Refine by pos(n-1) -- this is the "wrong" refinement for p=2
                    let mut by_q: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
                    for (pi, ms) in &source_perms {
                        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                        by_q.entry(q).or_default().push(*ms);
                    }
                    let mut source_rows: Vec<Vec<i64>> = Vec::new();
                    let mut keys: Vec<u8> = by_q.keys().copied().collect();
                    keys.sort();
                    for q in &keys {
                        source_rows.push(build_poly(by_q.get(q).unwrap()));
                    }
                    if source_rows.len() >= 2 {
                        n_p2_src_total += 1;
                        if check_tp2(&source_rows) {
                            n_p2_src_ok += 1;
                        }
                    }

                    // Refine by pos(1) -- the "right" refinement for p=2
                    let mut by_r: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
                    for (pi, ms) in &source_perms {
                        let r = pi.iter().position(|&v| v == 1).unwrap() as u8 + 1;
                        by_r.entry(r).or_default().push(*ms);
                    }
                    let mut pos1_rows: Vec<Vec<i64>> = Vec::new();
                    let mut keys: Vec<u8> = by_r.keys().copied().collect();
                    keys.sort();
                    for r in &keys {
                        pos1_rows.push(build_poly(by_r.get(r).unwrap()));
                    }
                    if pos1_rows.len() >= 2 {
                        n_p2_pos1_total += 1;
                        if check_tp2(&pos1_rows) {
                            n_p2_pos1_ok += 1;
                        } else {
                            let s_str = descent_set_to_string(mask, n);
                            println!("  p=2 pos(1)-source TP_2 FAIL: n={} S={}", n, s_str);
                        }
                    }
                }
            }
        }

        src_tp2_ok += n_src_ok as u64;
        src_tp2_total += n_src_total as u64;
        p2_src_tp2_ok += n_p2_src_ok as u64;
        p2_src_tp2_total += n_p2_src_total as u64;
        p2_pos1_tp2_ok += n_p2_pos1_ok as u64;
        p2_pos1_tp2_total += n_p2_pos1_total as u64;

        println!("n={}: source(p>=3) TP_2 {}/{}, p=2 pos(n-1) TP_2 {}/{}, p=2 pos(1) TP_2 {}/{}",
            n, n_src_ok, n_src_total, n_p2_src_ok, n_p2_src_total,
            n_p2_pos1_ok, n_p2_pos1_total);
    }

    println!("\n=== SUMMARY ===");
    println!("Source (p>=3, pos(n-1)-refined) TP_2:  {}/{}", src_tp2_ok, src_tp2_total);
    println!("Source (p=2, pos(n-1)-refined) TP_2:   {}/{}", p2_src_tp2_ok, p2_src_tp2_total);
    println!("Source (p=2, pos(1)-refined) TP_2:     {}/{}", p2_pos1_tp2_ok, p2_pos1_tp2_total);
}
