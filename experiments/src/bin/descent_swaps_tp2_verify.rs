/// Verify the relationship between TP_2 and interlacing for L^{(p)}.
///
/// Check:
/// 1. TP_2 of coefficient matrix
/// 2. Whether consecutive L^{(p)} interlace: L^{(p_i)} << L^{(p_{i+1})}
/// 3. Whether ALL pairs interlace (full interlacing sequence)
/// 4. Whether ALL pairs are compatible (pairwise common interlacing)
///
/// The key question: does TP_2 IMPLY the compatible family property?
///
/// Also test the NORMALIZED version: divide each L^{(p)} by its sum of coefficients
/// (or by its leading coefficient) before forming the matrix.
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

fn check_tp2(rows: &[Vec<i64>]) -> bool {
    let max_deg = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    for i in 0..rows.len() {
        for j in (i + 1)..rows.len() {
            for k1 in 0..max_deg {
                for k2 in (k1 + 1)..max_deg {
                    let minor = coeff(&rows[i], k1) * coeff(&rows[j], k2)
                        - coeff(&rows[i], k2) * coeff(&rows[j], k1);
                    if minor < 0 {
                        return false;
                    }
                }
            }
        }
    }
    true
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
        (1, 10),
        (10, 1),
        (1, 20),
        (20, 1),
        (7, 11),
        (11, 7),
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

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    println!("=== TP_2 vs Interlacing vs Compatibility ===\n");

    let mut tp2_ok = 0u64;
    let mut consec_interlace_ok = 0u64;
    let mut all_interlace_ok = 0u64;
    let mut all_compat_ok = 0u64;
    let mut total = 0u64;

    // Track: cases where tp2 holds but interlacing fails, or vice versa
    let mut tp2_not_interlace = 0u64;
    let mut interlace_not_tp2 = 0u64;

    for n in 4..=max_n {
        let perms = all_permutations(n);
        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            let mut lp_polys: Vec<(u8, Vec<i64>)> = Vec::new();
            for &p in &vp {
                let vals: Vec<usize> = class
                    .iter()
                    .filter(|pi| pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                    .map(|pi| compute(pi, Stat::Swaps))
                    .collect();
                if !vals.is_empty() {
                    lp_polys.push((p, build_poly(&vals)));
                }
            }

            if lp_polys.len() < 2 {
                continue;
            }

            total += 1;

            let rows: Vec<Vec<i64>> = lp_polys.iter().map(|(_, poly)| poly.clone()).collect();
            let is_tp2 = check_tp2(&rows);
            if is_tp2 {
                tp2_ok += 1;
            }

            // Check consecutive interlacing
            let mut consec_ok = true;
            for i in 0..rows.len() - 1 {
                if check_weak_interlacing(&rows[i], &rows[i + 1]) != Some(true) {
                    consec_ok = false;
                    break;
                }
            }
            if consec_ok {
                consec_interlace_ok += 1;
            }

            // Check all pairs interlace
            let mut all_inter = true;
            for i in 0..rows.len() {
                for j in (i + 1)..rows.len() {
                    if check_weak_interlacing(&rows[i], &rows[j]) != Some(true) {
                        all_inter = false;
                        break;
                    }
                }
                if !all_inter {
                    break;
                }
            }
            if all_inter {
                all_interlace_ok += 1;
            }

            // Check all pairs compatible
            let mut all_compat = true;
            for i in 0..rows.len() {
                for j in (i + 1)..rows.len() {
                    if !check_pairwise_compatible(&rows[i], &rows[j]) {
                        all_compat = false;
                        break;
                    }
                }
                if !all_compat {
                    break;
                }
            }
            if all_compat {
                all_compat_ok += 1;
            }

            // Cross-comparison
            if is_tp2 && !all_inter {
                tp2_not_interlace += 1;
                let s_str = descent_set_to_string(mask, n);
                println!("  TP_2 but NOT interlacing: n={} S={}", n, s_str);
                for (p, poly) in &lp_polys {
                    println!("    L^({}) = {}", p, format_poly(poly));
                }
            }
            if all_inter && !is_tp2 {
                interlace_not_tp2 += 1;
            }
        }

        println!(
            "n={}: total={}, TP_2={}, consec_inter={}, all_inter={}, compat={}",
            n, total, tp2_ok, consec_interlace_ok, all_interlace_ok, all_compat_ok
        );
    }

    println!("\n=== SUMMARY over n=4..{} ===", max_n);
    println!("Total families:          {}", total);
    println!("TP_2:                    {}", tp2_ok);
    println!("Consecutive interlacing: {}", consec_interlace_ok);
    println!("All pairs interlace:     {}", all_interlace_ok);
    println!("Compatible family:       {}", all_compat_ok);
    println!("TP_2 but not inter:      {}", tp2_not_interlace);
    println!("Inter but not TP_2:      {}", interlace_not_tp2);
}
