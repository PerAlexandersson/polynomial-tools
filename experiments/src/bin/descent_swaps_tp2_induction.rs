/// Check the inductive step for the TP_2 conjecture.
///
/// At level n, the coefficient matrix of {L_{n,S}^{(p)}} is TP_2.
/// Question: is the coefficient matrix of the SOURCE polynomials also TP_2?
///
/// For each target S with 1 not in S and each valid p in P(S):
///   1. Compute modified source polynomials F~_q^{(p)} (using swaps + eps1,
///      from both S'_p and S''_p sources)
///   2. Form the coefficient matrix of (F~_q^{(p)})_q ordered by q
///   3. Check if this source coefficient matrix is TP_2
///
/// Also check: for p=2 with the pos(1) refinement (G_r),
///   is the coefficient matrix of (G_r)_r TP_2?
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

fn descent_set_to_string(mask: u64, n: u8) -> String {
    let mut s = String::from("{");
    let mut first = true;
    for i in 1..n {
        if (mask >> (i - 1)) & 1 == 1 {
            if !first { s.push(','); }
            s.push_str(&i.to_string());
            first = false;
        }
    }
    s.push('}');
    s
}

fn coeff(poly: &[i64], k: usize) -> i64 {
    if k < poly.len() { poly[k] } else { 0 }
}

fn is_tp2(rows: &[Vec<i64>]) -> bool {
    let nrows = rows.len();
    if nrows < 2 { return true; }
    let max_deg = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    for i1 in 0..nrows {
        for i2 in (i1 + 1)..nrows {
            for j1 in 0..max_deg {
                for j2 in (j1 + 1)..max_deg {
                    let minor = coeff(&rows[i1], j1) * coeff(&rows[i2], j2)
                        - coeff(&rows[i1], j2) * coeff(&rows[i2], j1);
                    if minor < 0 {
                        return false;
                    }
                }
            }
        }
    }
    true
}

/// Find the first negative 2x2 minor (for diagnostics).
fn first_neg_minor(rows: &[Vec<i64>]) -> Option<(usize, usize, usize, usize, i64)> {
    let nrows = rows.len();
    if nrows < 2 { return None; }
    let max_deg = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    for i1 in 0..nrows {
        for i2 in (i1 + 1)..nrows {
            for j1 in 0..max_deg {
                for j2 in (j1 + 1)..max_deg {
                    let minor = coeff(&rows[i1], j1) * coeff(&rows[i2], j2)
                        - coeff(&rows[i1], j2) * coeff(&rows[i2], j1);
                    if minor < 0 {
                        return Some((i1, i2, j1, j2, minor));
                    }
                }
            }
        }
    }
    None
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    println!("=== TP_2 Induction Check: Source Coefficient Matrices ===");
    println!("Checking n = 4 to {}\n", max_n);

    let mut total_src_p_ge3 = 0u64;
    let mut ok_src_p_ge3 = 0u64;
    let mut total_src_p2_posn1 = 0u64;
    let mut ok_src_p2_posn1 = 0u64;
    let mut total_src_p2_pos1 = 0u64;
    let mut ok_src_p2_pos1 = 0u64;

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

        let mut n_src_p_ge3_total = 0u32;
        let mut n_src_p_ge3_ok = 0u32;
        let mut n_p2_posn1_total = 0u32;
        let mut n_p2_posn1_ok = 0u32;
        let mut n_p2_pos1_total = 0u32;
        let mut n_p2_pos1_ok = 0u32;

        for (&mask, _class) in &by_descent {
            // Only consider S with 1 not in S
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            for &p in &vp {
                if p <= 1 {
                    continue;
                }

                // Gather source permutations from both S'_p (ascent source) and S''_p (descent source)
                let sp_a = source_asc(mask, p, n);
                let sp_d = source_desc(mask, p, n);

                let mut source_perms: Vec<(&Vec<u8>, usize)> = Vec::new();
                for &sp in &[Some(sp_a), sp_d]
                    .iter()
                    .filter_map(|x| *x)
                    .collect::<Vec<_>>()
                {
                    if let Some(cls) = prev_by_des.get(&sp) {
                        for pi in cls {
                            let e1 = if epsilon1(pi, p) { 1 } else { 0 };
                            let ms = compute(pi, Stat::Swaps) + e1;
                            source_perms.push((pi, ms));
                        }
                    }
                }

                if source_perms.is_empty() {
                    continue;
                }

                if p >= 3 {
                    // ---- p >= 3: refine by pos(n-1) in the (n-1)-permutation ----
                    let mut by_q: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
                    for (pi, ms) in &source_perms {
                        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                        by_q.entry(q).or_default().push(*ms);
                    }

                    let keys: Vec<u8> = by_q.keys().copied().collect();
                    let source_rows: Vec<Vec<i64>> = keys
                        .iter()
                        .map(|q| build_poly(by_q.get(q).unwrap()))
                        .collect();

                    if source_rows.len() >= 2 {
                        n_src_p_ge3_total += 1;
                        if is_tp2(&source_rows) {
                            n_src_p_ge3_ok += 1;
                        } else {
                            let s_str = descent_set_to_string(mask, n);
                            let (i1, i2, j1, j2, val) =
                                first_neg_minor(&source_rows).unwrap();
                            println!(
                                "  FAIL p>=3 src TP_2: n={} S={} p={} rows({},{}) cols({},{}) minor={}",
                                n, s_str, p, keys[i1], keys[i2], j1, j2, val
                            );
                            for (idx, q) in keys.iter().enumerate() {
                                println!(
                                    "    q={}: {}",
                                    q,
                                    format_poly(&source_rows[idx])
                                );
                            }
                        }
                    }
                } else {
                    // ---- p == 2 ----

                    // Check 1: refine by pos(n-1) (the same refinement as p>=3 uses)
                    {
                        let mut by_q: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
                        for (pi, ms) in &source_perms {
                            let q =
                                pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            by_q.entry(q).or_default().push(*ms);
                        }

                        let keys: Vec<u8> = by_q.keys().copied().collect();
                        let source_rows: Vec<Vec<i64>> = keys
                            .iter()
                            .map(|q| build_poly(by_q.get(q).unwrap()))
                            .collect();

                        if source_rows.len() >= 2 {
                            n_p2_posn1_total += 1;
                            if is_tp2(&source_rows) {
                                n_p2_posn1_ok += 1;
                            } else {
                                let s_str = descent_set_to_string(mask, n);
                                let (i1, i2, j1, j2, val) =
                                    first_neg_minor(&source_rows).unwrap();
                                println!(
                                    "  FAIL p=2 pos(n-1) src TP_2: n={} S={} rows({},{}) cols({},{}) minor={}",
                                    n, s_str, keys[i1], keys[i2], j1, j2, val
                                );
                                for (idx, q) in keys.iter().enumerate() {
                                    println!(
                                        "    q={}: {}",
                                        q,
                                        format_poly(&source_rows[idx])
                                    );
                                }
                            }
                        }
                    }

                    // Check 2: refine by pos(1) -- the natural refinement for p=2
                    {
                        let mut by_r: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
                        for (pi, ms) in &source_perms {
                            let r = pi.iter().position(|&v| v == 1).unwrap() as u8 + 1;
                            by_r.entry(r).or_default().push(*ms);
                        }

                        let keys: Vec<u8> = by_r.keys().copied().collect();
                        let pos1_rows: Vec<Vec<i64>> = keys
                            .iter()
                            .map(|r| build_poly(by_r.get(r).unwrap()))
                            .collect();

                        if pos1_rows.len() >= 2 {
                            n_p2_pos1_total += 1;
                            if is_tp2(&pos1_rows) {
                                n_p2_pos1_ok += 1;
                            } else {
                                let s_str = descent_set_to_string(mask, n);
                                let (i1, i2, j1, j2, val) =
                                    first_neg_minor(&pos1_rows).unwrap();
                                println!(
                                    "  FAIL p=2 pos(1) src TP_2: n={} S={} rows({},{}) cols({},{}) minor={}",
                                    n, s_str, keys[i1], keys[i2], j1, j2, val
                                );
                                for (idx, r) in keys.iter().enumerate() {
                                    println!(
                                        "    r={}: {}",
                                        r,
                                        format_poly(&pos1_rows[idx])
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }

        // Also verify TP_2 of the TARGET coefficient matrix (L^{(p)}) for reference
        let mut n_target_total = 0u32;
        let mut n_target_ok = 0u32;
        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            let mut target_rows: Vec<Vec<i64>> = Vec::new();
            let mut target_labels: Vec<u8> = Vec::new();
            for &p in &vp {
                let vals: Vec<usize> = class
                    .iter()
                    .filter(|pi| {
                        pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p
                    })
                    .map(|pi| compute(pi, Stat::Swaps))
                    .collect();
                if !vals.is_empty() {
                    target_rows.push(build_poly(&vals));
                    target_labels.push(p);
                }
            }

            if target_rows.len() >= 2 {
                n_target_total += 1;
                if is_tp2(&target_rows) {
                    n_target_ok += 1;
                } else {
                    let s_str = descent_set_to_string(mask, n);
                    println!("  FAIL target TP_2: n={} S={}", n, s_str);
                }
            }
        }

        total_src_p_ge3 += n_src_p_ge3_total as u64;
        ok_src_p_ge3 += n_src_p_ge3_ok as u64;
        total_src_p2_posn1 += n_p2_posn1_total as u64;
        ok_src_p2_posn1 += n_p2_posn1_ok as u64;
        total_src_p2_pos1 += n_p2_pos1_total as u64;
        ok_src_p2_pos1 += n_p2_pos1_ok as u64;

        println!(
            "n={}: target TP_2 {}/{} | src(p>=3,pos(n-1)) TP_2 {}/{} | p=2 pos(n-1) TP_2 {}/{} | p=2 pos(1) TP_2 {}/{}",
            n,
            n_target_ok, n_target_total,
            n_src_p_ge3_ok, n_src_p_ge3_total,
            n_p2_posn1_ok, n_p2_posn1_total,
            n_p2_pos1_ok, n_p2_pos1_total,
        );
    }

    println!("\n=== SUMMARY (n=4..{}) ===", max_n);
    println!(
        "Source (p>=3, pos(n-1)-refined) TP_2:  {}/{}",
        ok_src_p_ge3, total_src_p_ge3
    );
    println!(
        "Source (p=2,  pos(n-1)-refined) TP_2:  {}/{}",
        ok_src_p2_posn1, total_src_p2_posn1
    );
    println!(
        "Source (p=2,  pos(1)-refined)   TP_2:  {}/{}",
        ok_src_p2_pos1, total_src_p2_pos1
    );

    let all_ok = ok_src_p_ge3 == total_src_p_ge3
        && ok_src_p2_pos1 == total_src_p2_pos1;
    println!(
        "\nInductive source TP_2 (p>=3 + p=2 pos(1)): {}",
        if all_ok { "ALL PASS" } else { "SOME FAIL" }
    );
}
