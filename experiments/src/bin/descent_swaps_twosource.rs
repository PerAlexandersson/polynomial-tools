/// Verify the two-source structure of the insertion recurrence.
///
/// For target descent set S in S_n and insertion position p,
/// the source permutations π ∈ S_{n-1} decompose into:
///   - π ∈ D(n-1, S'_p) where position p-1 is an ASCENT in π
///   - π ∈ D(n-1, S''_p) where position p-1 is a DESCENT in π
///     (with S''_p = S'_p ∪ {p-1})
///
/// Key property: the S''_p source has ε₁ = 0 always (pure staircase).
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::format_poly;
use std::collections::BTreeMap;

fn build_poly_from_vals(vals: &[usize]) -> Vec<i64> {
    if vals.is_empty() { return vec![0]; }
    let max_s = *vals.iter().max().unwrap();
    let mut coeffs = vec![0i64; max_s + 1];
    for &s in vals { coeffs[s] += 1; }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 { coeffs.pop(); }
    coeffs
}

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut result = vec![0i64; len];
    for (i, &c) in a.iter().enumerate() { result[i] += c; }
    for (i, &c) in b.iter().enumerate() { result[i] += c; }
    while result.len() > 1 && *result.last().unwrap() == 0 { result.pop(); }
    result
}

fn valid_positions(s_mask: u64, n: u8) -> Vec<u8> {
    let mut positions = Vec::new();
    for p in 1..n {
        let p_in_s = (s_mask >> (p - 1)) & 1 == 1;
        let pm1_in_s = if p >= 2 { (s_mask >> (p - 2)) & 1 == 1 } else { false };
        if p_in_s && !pm1_in_s { positions.push(p); }
    }
    if n >= 2 && (s_mask >> (n - 2)) & 1 == 0 { positions.push(n); }
    positions
}

/// Source descent set S'_p (ascent at p-1)
fn source_asc(s_mask: u64, p: u8, n: u8) -> u64 {
    if n <= 2 { return 0; }
    if p == n { return s_mask; }
    let mut sp = 0u64;
    if p == 1 {
        for j in 2..n {
            if (s_mask >> (j - 1)) & 1 == 1 { sp |= 1 << (j - 2); }
        }
    } else {
        for pos in 1..=(p.saturating_sub(2)) {
            if (s_mask >> (pos - 1)) & 1 == 1 { sp |= 1 << (pos - 1); }
        }
        for j in (p + 1)..n {
            if (s_mask >> (j - 1)) & 1 == 1 { sp |= 1 << (j - 2); }
        }
    }
    sp
}

/// Source descent set S''_p = S'_p ∪ {p-1} (descent at p-1)
fn source_desc(s_mask: u64, p: u8, n: u8) -> Option<u64> {
    if p <= 1 || p >= n { return None; } // Only for 1 < p < n
    let sp = source_asc(s_mask, p, n);
    Some(sp | (1 << (p - 2))) // Add bit for position p-1
}

fn epsilon1(pi: &[u8], p: u8) -> bool {
    let n = pi.len() as u8 + 1;
    if p <= 1 || p >= n { return false; }
    pi[(p - 2) as usize] + 1 == pi[(p - 1) as usize]
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

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    println!("Verifying two-source recurrence structure\n");

    let mut total_checks = 0u64;
    let mut total_ok = 0u64;
    let mut total_eps1_zero_in_desc = 0u64;
    let mut total_eps1_checks_in_desc = 0u64;

    for n in 4..=max_n {
        let perms = all_permutations(n);
        let perms_prev = all_permutations(n - 1);

        // Group current perms
        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent.entry(descent_set_bitmask(pi)).or_default().push(pi);
        }

        // Group prev perms
        let mut by_descent_prev: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms_prev {
            by_descent_prev.entry(descent_set_bitmask(pi)).or_default().push(pi);
        }

        let mut n_ok = 0;
        let mut n_total = 0;
        let mut n_eps1_zero = 0;
        let mut n_eps1_checks = 0;

        for (&mask, target_class) in &by_descent {
            let vp = valid_positions(mask, n);

            for &p in &vp {
                // Compute actual L_{n,S}^{(p)}
                let target_vals: Vec<usize> = target_class.iter()
                    .filter(|pi| pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                    .map(|pi| compute(pi, Stat::Swaps))
                    .collect();
                if target_vals.is_empty() { continue; }
                let target_poly = build_poly_from_vals(&target_vals);

                // Compute recurrence with both sources
                let sp_asc = source_asc(mask, p, n);
                let sp_desc = source_desc(mask, p, n);

                let mut recurrence_poly = vec![0i64; 1];

                // Source 1: S'_p (ascent at p-1)
                if let Some(asc_class) = by_descent_prev.get(&sp_asc) {
                    for pi in asc_class {
                        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                        let eps2 = if q <= p.saturating_sub(2) { 1usize } else { 0 };
                        let e1 = if epsilon1(pi, p) { 1usize } else { 0 };
                        let s = compute(pi, Stat::Swaps) + eps2 + e1;
                        if s >= recurrence_poly.len() {
                            recurrence_poly.resize(s + 1, 0);
                        }
                        recurrence_poly[s] += 1;
                    }
                }

                // Source 2: S''_p (descent at p-1) — SHOULD have ε₁ = 0
                if let Some(sp_d) = sp_desc {
                    if let Some(desc_class) = by_descent_prev.get(&sp_d) {
                        for pi in desc_class {
                            let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                            let eps2 = if q <= p.saturating_sub(2) { 1usize } else { 0 };
                            let e1 = epsilon1(pi, p);
                            n_eps1_checks += 1;
                            if !e1 {
                                n_eps1_zero += 1;
                            } else {
                                // This should NEVER happen
                                let s_str = descent_set_to_string(mask, n);
                                println!("  BUG: eps1=1 in desc source! n={} S={} p={} pi={:?}", n, s_str, p, pi);
                            }
                            let s = compute(pi, Stat::Swaps) + eps2;
                            if s >= recurrence_poly.len() {
                                recurrence_poly.resize(s + 1, 0);
                            }
                            recurrence_poly[s] += 1;
                        }
                    }
                }

                // Trim
                while recurrence_poly.len() > 1 && *recurrence_poly.last().unwrap() == 0 {
                    recurrence_poly.pop();
                }

                n_total += 1;
                if recurrence_poly == target_poly {
                    n_ok += 1;
                } else {
                    let s_str = descent_set_to_string(mask, n);
                    println!("  MISMATCH n={} S={} p={}:", n, s_str, p);
                    println!("    target:     {}", format_poly(&target_poly));
                    println!("    recurrence: {}", format_poly(&recurrence_poly));
                }
            }
        }

        total_checks += n_total;
        total_ok += n_ok;
        total_eps1_zero_in_desc += n_eps1_zero;
        total_eps1_checks_in_desc += n_eps1_checks;
        println!("n={}: recurrence {}/{}, eps1=0 in desc-source: {}/{}",
            n, n_ok, n_total, n_eps1_zero, n_eps1_checks);
    }

    println!("\n=== SUMMARY ===");
    println!("Two-source recurrence correct: {}/{}", total_ok, total_checks);
    println!("eps1 = 0 in descent-source: {}/{}", total_eps1_zero_in_desc, total_eps1_checks_in_desc);
    println!("\nThe S''_p source (descent at p-1) is ALWAYS correction-free.");
    println!("Only the S'_p source (ascent at p-1) has nonzero corrections.");
}
