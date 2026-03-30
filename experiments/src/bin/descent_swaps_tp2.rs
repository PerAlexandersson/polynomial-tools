/// Approach 6: Total positivity of the coefficient matrix.
///
/// For each (n, S) with 1 not in S, form the matrix M where:
///   Row index: p in P(S), ordered increasingly
///   Column index: k = 0, 1, ..., max_degree
///   M_{p,k} = [t^k] L^{(p)}(t)
///
/// Check: Is M totally positive of order 2 (TP_2)?
/// I.e., are all 2x2 minors non-negative?
///
/// If TP_2 holds, then the polynomials form an interlacing sequence
/// (which is stronger than a compatible family).
///
/// Also check TP_3 (all 3x3 minors non-negative).
///
/// Additionally, check the STRONGER condition:
///   The coefficient matrix with rows ordered by p is TP_2.
///   This would imply interlacing, hence compatibility.
///
/// Separate test: coefficient matrix of the MODIFIED SOURCE vectors F~_q^{(p)}.
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

fn descent_set_to_string(mask: u64, n: u8) -> String {
    let mut s = String::from("{");
    let mut first = true;
    for i in 1..n { if (mask >> (i - 1)) & 1 == 1 { if !first { s.push(','); } s.push_str(&i.to_string()); first = false; } }
    s.push('}'); s
}

fn coeff(poly: &[i64], k: usize) -> i64 {
    if k < poly.len() { poly[k] } else { 0 }
}

/// Check all 2x2 minors of matrix M are non-negative.
/// M has 'rows' rows, each a polynomial (coefficient vector).
fn check_tp2(rows: &[Vec<i64>]) -> (bool, Option<(usize, usize, usize, usize, i64)>) {
    let max_deg = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    for i in 0..rows.len() {
        for j in (i+1)..rows.len() {
            for k1 in 0..max_deg {
                for k2 in (k1+1)..max_deg {
                    let minor = coeff(&rows[i], k1) * coeff(&rows[j], k2)
                              - coeff(&rows[i], k2) * coeff(&rows[j], k1);
                    if minor < 0 {
                        return (false, Some((i, j, k1, k2, minor)));
                    }
                }
            }
        }
    }
    (true, None)
}

/// Check all 3x3 minors are non-negative.
fn check_tp3(rows: &[Vec<i64>]) -> bool {
    let max_deg = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    if rows.len() < 3 || max_deg < 3 { return true; }
    for i in 0..rows.len() {
        for j in (i+1)..rows.len() {
            for l in (j+1)..rows.len() {
                for k1 in 0..max_deg {
                    for k2 in (k1+1)..max_deg {
                        for k3 in (k2+1)..max_deg {
                            let a = coeff(&rows[i], k1);
                            let b = coeff(&rows[i], k2);
                            let c = coeff(&rows[i], k3);
                            let d = coeff(&rows[j], k1);
                            let e = coeff(&rows[j], k2);
                            let f = coeff(&rows[j], k3);
                            let g = coeff(&rows[l], k1);
                            let h = coeff(&rows[l], k2);
                            let k = coeff(&rows[l], k3);
                            let det = a*(e*k - f*h) - b*(d*k - f*g) + c*(d*h - e*g);
                            if det < 0 { return false; }
                        }
                    }
                }
            }
        }
    }
    true
}

fn main() {
    let max_n: u8 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(9);

    println!("=== Approach 6: Total positivity of coefficient matrix ===\n");

    let mut tp2_ok = 0u64;
    let mut tp2_total = 0u64;
    let mut tp3_ok = 0u64;
    let mut tp3_total = 0u64;

    // Also test: just the p>=3 rows (should be TP_2 since they interlace)
    let mut tp2_p3_ok = 0u64;
    let mut tp2_p3_total = 0u64;

    // Test: can we REVERSE the row order and get TP_2?
    let mut tp2_rev_ok = 0u64;
    let mut tp2_rev_total = 0u64;

    for n in 4..=max_n {
        let perms = all_permutations(n);
        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms { by_descent.entry(descent_set_bitmask(pi)).or_default().push(pi); }

        let mut n_tp2_ok = 0u32;
        let mut n_tp2_total = 0u32;
        let mut n_tp3_ok = 0u32;
        let mut n_tp3_total = 0u32;
        let mut n_tp2_p3_ok = 0u32;
        let mut n_tp2_p3_total = 0u32;
        let mut n_tp2_rev_ok = 0u32;
        let mut n_tp2_rev_total = 0u32;

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 { continue; }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 { continue; }

            // Build L^{(p)} for each p in P(S)
            let mut lp_polys: Vec<(u8, Vec<i64>)> = Vec::new();
            for &p in &vp {
                let vals: Vec<usize> = class.iter()
                    .filter(|pi| pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                    .map(|pi| compute(pi, Stat::Swaps))
                    .collect();
                if !vals.is_empty() {
                    lp_polys.push((p, build_poly(&vals)));
                }
            }

            if lp_polys.len() < 2 { continue; }

            // Full coefficient matrix (all p, ordered by p)
            let rows: Vec<Vec<i64>> = lp_polys.iter().map(|(_, poly)| poly.clone()).collect();

            n_tp2_total += 1;
            let (is_tp2, fail_info) = check_tp2(&rows);
            if is_tp2 {
                n_tp2_ok += 1;
            } else {
                if let Some((ri, rj, k1, k2, val)) = fail_info {
                    let s_str = descent_set_to_string(mask, n);
                    let (p_i, _) = &lp_polys[ri];
                    let (p_j, _) = &lp_polys[rj];
                    println!("  TP_2 FAIL: n={} S={} rows p={},p={} cols {},{} minor={}",
                        n, s_str, p_i, p_j, k1, k2, val);
                    println!("    L^({}) = {}", p_i, format_poly(&lp_polys[ri].1));
                    println!("    L^({}) = {}", p_j, format_poly(&lp_polys[rj].1));
                }
            }

            // TP_3
            n_tp3_total += 1;
            if check_tp3(&rows) { n_tp3_ok += 1; }

            // Just p >= 3
            let rows_p3: Vec<Vec<i64>> = lp_polys.iter()
                .filter(|(p, _)| *p >= 3)
                .map(|(_, poly)| poly.clone())
                .collect();
            if rows_p3.len() >= 2 {
                n_tp2_p3_total += 1;
                let (is_tp2_p3, _) = check_tp2(&rows_p3);
                if is_tp2_p3 { n_tp2_p3_ok += 1; }
            }

            // Reversed order
            let mut rows_rev = rows.clone();
            rows_rev.reverse();
            n_tp2_rev_total += 1;
            let (is_tp2_rev, _) = check_tp2(&rows_rev);
            if is_tp2_rev { n_tp2_rev_ok += 1; }
        }

        tp2_ok += n_tp2_ok as u64;
        tp2_total += n_tp2_total as u64;
        tp3_ok += n_tp3_ok as u64;
        tp3_total += n_tp3_total as u64;
        tp2_p3_ok += n_tp2_p3_ok as u64;
        tp2_p3_total += n_tp2_p3_total as u64;
        tp2_rev_ok += n_tp2_rev_ok as u64;
        tp2_rev_total += n_tp2_rev_total as u64;

        println!("n={}: TP_2 {}/{}, TP_3 {}/{}, TP_2(p>=3) {}/{}, TP_2(rev) {}/{}",
            n, n_tp2_ok, n_tp2_total, n_tp3_ok, n_tp3_total,
            n_tp2_p3_ok, n_tp2_p3_total, n_tp2_rev_ok, n_tp2_rev_total);
    }

    println!("\n=== SUMMARY ===");
    println!("TP_2 (all p, natural order):  {}/{}", tp2_ok, tp2_total);
    println!("TP_3 (all p):                 {}/{}", tp3_ok, tp3_total);
    println!("TP_2 (p>=3 only):             {}/{}", tp2_p3_ok, tp2_p3_total);
    println!("TP_2 (reversed order):        {}/{}", tp2_rev_ok, tp2_rev_total);
    println!();
    println!("If TP_2 passes for all, the L^{{(p)}} form an interlacing sequence,");
    println!("which is stronger than a compatible family.");
    println!("If only p>=3 is TP_2, we still need the p=2 gap.");
}
