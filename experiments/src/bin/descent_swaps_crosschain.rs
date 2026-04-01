/// Cross-chain interlacing analysis for the modified-source staircase approach.
///
/// For each target S with 1 not in S and consecutive valid positions p_a < p_b in P(S):
///   Chain A = {f_{p_a, q}}_q   (modified source polys for position p_a)
///   Chain B = {f_{p_b, q}}_q   (modified source polys for position p_b)
///
/// Tests:
/// 1. Do chains A and B individually interlace? (Expected: yes)
/// 2. Does every f_{p_a,q} interlace every f_{p_b,r}? (Full pairwise cross-interlacing)
/// 3. For each shared q, does f_{p_a,q} interlace f_{p_b,q}? (Same-q cross)
/// 4. Does the concatenation (chain A, then chain B) form an interlacing sequence?
/// 5. Ratio analysis: is f_{p_b,q}/f_{p_a,q} non-decreasing in coefficients?
/// 6. Do the staircase outputs L^{(p_a)} and L^{(p_b)} share common roots?
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use num_traits::Signed;
use polynomial_tools::real_rootedness::{
    check_weak_interlacing, format_poly, is_real_rooted, real_roots,
};
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

fn source_asc(s_mask: u64, p: u8, n: u8) -> u64 {
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

fn source_desc(s_mask: u64, p: u8, n: u8) -> Option<u64> {
    if p <= 1 || p >= n {
        return None;
    }
    Some(source_asc(s_mask, p, n) | (1 << (p - 2)))
}

fn epsilon1(pi: &[u8], p: u8) -> bool {
    let n = pi.len() as u8 + 1;
    if p <= 1 || p >= n {
        return false;
    }
    pi[(p - 2) as usize] + 1 == pi[(p - 1) as usize]
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

fn poly_is_zero(p: &[i64]) -> bool {
    p.iter().all(|&c| c == 0)
}

fn poly_degree(p: &[i64]) -> Option<usize> {
    for i in (0..p.len()).rev() {
        if p[i] != 0 {
            return Some(i);
        }
    }
    None
}

/// Check if f and g are pairwise compatible (i.e., all nonnegative combinations are real-rooted).
fn check_pairwise_compatible(f: &[i64], g: &[i64]) -> bool {
    if poly_is_zero(f) || poly_is_zero(g) {
        return true;
    }
    if f.len() <= 1 && g.len() <= 1 {
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

/// Check if a sequence of polynomials forms an interlacing chain.
/// That means each consecutive pair is pairwise compatible (interlacing family).
fn check_chain_interlacing(polys: &[Vec<i64>]) -> bool {
    if polys.len() < 2 {
        return true;
    }
    for i in 0..polys.len() {
        for j in (i + 1)..polys.len() {
            if !check_pairwise_compatible(&polys[i], &polys[j]) {
                return false;
            }
        }
    }
    true
}

/// Compute the staircase output L^{(p)} = sum_q t^{eps2(p,q)} * f_{p,q}
fn staircase_output(f_polys: &[(u8, Vec<i64>)], p: u8) -> Vec<i64> {
    let mut result = vec![0i64; 1];
    for (q, poly) in f_polys {
        let eps2 = if *q <= p.saturating_sub(2) { 1usize } else { 0 };
        for (k, &c) in poly.iter().enumerate() {
            let idx = k + eps2;
            if idx >= result.len() {
                result.resize(idx + 1, 0);
            }
            result[idx] += c;
        }
    }
    while result.len() > 1 && *result.last().unwrap() == 0 {
        result.pop();
    }
    result
}

/// Build the modified source chain for a given target S and position p.
/// Returns (sorted_qs, f_polys) where f_polys[i] = f_{p, sorted_qs[i]}.
fn build_source_chain(
    s_mask: u64,
    p: u8,
    n: u8,
    by_descent_prev: &BTreeMap<u64, Vec<usize>>,
    perm_prev_list: &[&Vec<u8>],
) -> (Vec<u8>, Vec<Vec<i64>>) {
    let sp_a = source_asc(s_mask, p, n);
    let sp_d = source_desc(s_mask, p, n);

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
                let s = compute(pi, Stat::Swaps); // eps1=0 always for descent source
                by_q.entry(q).or_default().push(s);
            }
        }
    }

    let sorted_qs: Vec<u8> = by_q.keys().copied().collect();
    let f_polys: Vec<Vec<i64>> = sorted_qs.iter().map(|q| build_poly(&by_q[q])).collect();

    (sorted_qs, f_polys)
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    println!("=== Cross-chain interlacing analysis ===");
    println!(
        "For consecutive valid positions p_a < p_b, analyze how source chains A and B relate.\n"
    );

    // Global counters
    let mut total_pairs = 0u64;
    let mut chain_a_ok = 0u64;
    let mut chain_b_ok = 0u64;
    let mut full_cross_ok = 0u64;
    let mut same_q_ok = 0u64;
    let mut same_q_tested = 0u64;
    let mut concat_ok = 0u64;
    let mut ratio_nondecr = 0u64;
    let mut ratio_tested = 0u64;
    let mut shared_roots_count = 0u64;
    let mut shared_roots_tested = 0u64;

    for n in 5..=max_n {
        let perms_prev = all_permutations(n - 1);

        let mut by_descent_prev: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
        let mut perm_prev_list: Vec<&Vec<u8>> = Vec::new();
        for pi in &perms_prev {
            let idx = perm_prev_list.len();
            perm_prev_list.push(pi);
            by_descent_prev
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(idx);
        }

        let mut n_pairs = 0u32;
        let mut n_chain_a_ok = 0u32;
        let mut n_chain_b_ok = 0u32;
        let mut n_full_cross_ok = 0u32;
        let mut n_same_q_ok = 0u32;
        let mut n_same_q_tested = 0u32;
        let mut n_concat_ok = 0u32;
        let mut n_ratio_nondecr = 0u32;
        let mut n_ratio_tested = 0u32;
        let mut n_shared_roots = 0u32;
        let mut n_shared_roots_tested = 0u32;

        for s_mask in 0u64..(1 << (n - 1)) {
            if s_mask & 1 != 0 {
                continue; // 1 not in S
            }
            let vp = valid_positions(s_mask, n);
            if vp.len() < 2 {
                continue;
            }

            // Build source chains for all valid positions
            let mut chains: Vec<(u8, Vec<u8>, Vec<Vec<i64>>)> = Vec::new();
            for &p in &vp {
                let (qs, fs) = build_source_chain(s_mask, p, n, &by_descent_prev, &perm_prev_list);
                chains.push((p, qs, fs));
            }

            // Build staircase outputs for all valid positions
            let mut staircase_outputs: Vec<(u8, Vec<i64>)> = Vec::new();
            for (p, qs, fs) in &chains {
                let qf: Vec<(u8, Vec<i64>)> = qs.iter().copied().zip(fs.iter().cloned()).collect();
                let l_p = staircase_output(&qf, *p);
                staircase_outputs.push((*p, l_p));
            }

            // For each consecutive pair of valid positions
            for idx in 0..vp.len() - 1 {
                let p_a = vp[idx];
                let p_b = vp[idx + 1];
                n_pairs += 1;

                let (_, ref qs_a, ref fs_a) = chains[idx];
                let (_, ref qs_b, ref fs_b) = chains[idx + 1];

                let s_str = descent_set_to_string(s_mask, n);

                // ---- TEST 1: Individual chain interlacing ----
                let chain_a_interlaces = check_chain_interlacing(fs_a);
                let chain_b_interlaces = check_chain_interlacing(fs_b);
                if chain_a_interlaces {
                    n_chain_a_ok += 1;
                }
                if chain_b_interlaces {
                    n_chain_b_ok += 1;
                }

                // ---- TEST 2: Full pairwise cross-interlacing ----
                let mut full_cross = true;
                'cross: for fa in fs_a {
                    for fb in fs_b {
                        if !check_pairwise_compatible(fa, fb) {
                            full_cross = false;
                            break 'cross;
                        }
                    }
                }
                if full_cross {
                    n_full_cross_ok += 1;
                } else if n <= 7 {
                    println!("  CROSS FAIL: n={} S={} p_a={} p_b={}", n, s_str, p_a, p_b);
                    for (i, q) in qs_a.iter().enumerate() {
                        println!("    f_{{p_a={},q={}}} = {}", p_a, q, format_poly(&fs_a[i]));
                    }
                    for (i, q) in qs_b.iter().enumerate() {
                        println!("    f_{{p_b={},q={}}} = {}", p_b, q, format_poly(&fs_b[i]));
                    }
                    // Find which pair fails
                    for (ia, fa) in fs_a.iter().enumerate() {
                        for (ib, fb) in fs_b.iter().enumerate() {
                            if !check_pairwise_compatible(fa, fb) {
                                println!(
                                    "    FAIL pair: f_{{p_a,q={}}} vs f_{{p_b,q={}}}",
                                    qs_a[ia], qs_b[ib]
                                );
                            }
                        }
                    }
                }

                // ---- TEST 3: Same-q cross-interlacing ----
                let mut all_same_q_ok = true;
                let mut any_same_q = false;
                for (ia, &qa) in qs_a.iter().enumerate() {
                    if let Some(ib) = qs_b.iter().position(|&qb| qb == qa) {
                        any_same_q = true;
                        n_same_q_tested += 1;
                        if !check_pairwise_compatible(&fs_a[ia], &fs_b[ib]) {
                            all_same_q_ok = false;
                            if n <= 7 {
                                println!(
                                    "  SAME-Q FAIL: n={} S={} p_a={} p_b={} q={}",
                                    n, s_str, p_a, p_b, qa
                                );
                                println!(
                                    "    f_{{p_a,q}} = {}  vs  f_{{p_b,q}} = {}",
                                    format_poly(&fs_a[ia]),
                                    format_poly(&fs_b[ib])
                                );
                            }
                        } else {
                            n_same_q_ok += 1;
                        }
                    }
                }
                if !any_same_q && n <= 7 {
                    // No shared q values at all
                }

                // ---- TEST 4: Concatenation (A then B) interlacing ----
                let mut concat: Vec<Vec<i64>> = fs_a.clone();
                concat.extend(fs_b.iter().cloned());
                let concat_interlaces = check_chain_interlacing(&concat);
                if concat_interlaces {
                    n_concat_ok += 1;
                } else if n <= 7 {
                    println!("  CONCAT FAIL: n={} S={} p_a={} p_b={}", n, s_str, p_a, p_b);
                }

                // ---- TEST 5: Ratio analysis for shared q ----
                for (ia, &qa) in qs_a.iter().enumerate() {
                    if let Some(ib) = qs_b.iter().position(|&qb| qb == qa) {
                        let fa = &fs_a[ia];
                        let fb = &fs_b[ib];
                        if poly_is_zero(fa) || poly_is_zero(fb) {
                            continue;
                        }
                        n_ratio_tested += 1;

                        // Check if the ratio fb/fa is "non-decreasing" in the coefficient sense.
                        // That is: for all k, fb[k]/fa[k] <= fb[k+1]/fa[k+1]
                        // (using cross-multiplication to stay in integers).
                        // This only makes sense if both have same support or we use
                        // the coefficient-wise ratio fb[k]*fa[k+1] >= fb[k+1]*fa[k].
                        // Actually, non-decreasing ratio means:
                        //   fb[k]/fa[k] <= fb[k+1]/fa[k+1]
                        //   i.e., fb[k]*fa[k+1] <= fb[k+1]*fa[k]
                        // Wait, that is the LR ordering condition!
                        // LR ordering (B dominates A): for all k, A[k]*B[k+1] >= A[k+1]*B[k]
                        // where A = f_{p_a,q}, B = f_{p_b,q}.
                        let max_deg = fa.len().max(fb.len());
                        let coeff_a = |k: usize| -> i64 {
                            if k < fa.len() {
                                fa[k]
                            } else {
                                0
                            }
                        };
                        let coeff_b = |k: usize| -> i64 {
                            if k < fb.len() {
                                fb[k]
                            } else {
                                0
                            }
                        };

                        let mut is_lr = true;
                        for k in 0..max_deg.saturating_sub(1) {
                            // LR: A[k]*B[k+1] >= A[k+1]*B[k]
                            let lhs = coeff_a(k) * coeff_b(k + 1);
                            let rhs = coeff_a(k + 1) * coeff_b(k);
                            if lhs < rhs {
                                is_lr = false;
                                break;
                            }
                        }

                        if is_lr {
                            n_ratio_nondecr += 1;
                        } else if n <= 7 {
                            println!(
                                "  RATIO NOT LR: n={} S={} p_a={} p_b={} q={}",
                                n, s_str, p_a, p_b, qa
                            );
                            println!(
                                "    f_{{p_a,q}} = {}  vs  f_{{p_b,q}} = {}",
                                format_poly(fa),
                                format_poly(fb)
                            );
                        }
                    }
                }

                // ---- TEST 6: Shared roots between staircase outputs ----
                let l_a = &staircase_outputs[idx].1;
                let l_b = &staircase_outputs[idx + 1].1;

                if !poly_is_zero(l_a) && !poly_is_zero(l_b) && l_a.len() > 1 && l_b.len() > 1 {
                    n_shared_roots_tested += 1;
                    let roots_a = real_roots(l_a);
                    let roots_b = real_roots(l_b);

                    if let (Some(ra), Some(rb)) = (roots_a, roots_b) {
                        // Check for approximately shared roots.
                        // Use rational comparison: for each root in ra, check if any root in rb
                        // is very close (within a tolerance).
                        let mut shared = 0;
                        let tolerance = num_rational::Ratio::new(
                            num_bigint::BigInt::from(1),
                            num_bigint::BigInt::from(10000),
                        );
                        for root_a in &ra {
                            for root_b in &rb {
                                let diff = root_a - root_b;
                                if diff.abs() < tolerance {
                                    shared += 1;
                                }
                            }
                        }
                        if shared > 0 {
                            n_shared_roots += 1;
                            if n <= 7 {
                                println!(
                                    "  SHARED ROOTS: n={} S={} p_a={} p_b={} ({} shared)",
                                    n, s_str, p_a, p_b, shared
                                );
                                println!("    L^{{(p_a)}} = {}", format_poly(l_a));
                                println!("    L^{{(p_b)}} = {}", format_poly(l_b));
                            }
                        }
                    }
                }

                // Print detailed info for small n
                if n <= 6 {
                    println!("\nn={} S={} p_a={} p_b={}", n, s_str, p_a, p_b);
                    println!("  Chain A (p_a={}):", p_a);
                    for (i, q) in qs_a.iter().enumerate() {
                        println!("    q={}: {}", q, format_poly(&fs_a[i]));
                    }
                    println!("  Chain B (p_b={}):", p_b);
                    for (i, q) in qs_b.iter().enumerate() {
                        println!("    q={}: {}", q, format_poly(&fs_b[i]));
                    }
                    println!("  L^{{(p_a)}} = {}", format_poly(l_a));
                    println!("  L^{{(p_b)}} = {}", format_poly(l_b));
                    println!(
                        "  Chain A interlaces: {}  Chain B interlaces: {}",
                        chain_a_interlaces, chain_b_interlaces
                    );
                    println!("  Full cross-interlacing: {}", full_cross);
                    println!("  Same-q compatible: {}", all_same_q_ok);
                    println!("  Concatenation interlaces: {}", concat_interlaces);

                    // Show interlacing relationship between outputs
                    let interlace_result = if l_a.len() <= 1 || l_b.len() <= 1 {
                        "trivial (degree <= 0)".to_string()
                    } else {
                        let (small, large) = if l_a.len() <= l_b.len() {
                            (l_a, l_b)
                        } else {
                            (l_b, l_a)
                        };
                        match check_weak_interlacing(small, large) {
                            Some(true) => "YES".to_string(),
                            Some(false) => "NO".to_string(),
                            None => "N/A (degree mismatch)".to_string(),
                        }
                    };
                    println!("  L^{{(p_a)}} interlaces L^{{(p_b)}}: {}", interlace_result);
                }
            }
        }

        total_pairs += n_pairs as u64;
        chain_a_ok += n_chain_a_ok as u64;
        chain_b_ok += n_chain_b_ok as u64;
        full_cross_ok += n_full_cross_ok as u64;
        same_q_ok += n_same_q_ok as u64;
        same_q_tested += n_same_q_tested as u64;
        concat_ok += n_concat_ok as u64;
        ratio_nondecr += n_ratio_nondecr as u64;
        ratio_tested += n_ratio_tested as u64;
        shared_roots_count += n_shared_roots as u64;
        shared_roots_tested += n_shared_roots_tested as u64;

        println!("\n--- n={} summary ({} consecutive pairs) ---", n, n_pairs);
        println!(
            "  TEST 1 - Chain A interlaces: {}/{}  Chain B interlaces: {}/{}",
            n_chain_a_ok, n_pairs, n_chain_b_ok, n_pairs
        );
        println!(
            "  TEST 2 - Full cross-interlacing (every f_{{p_a,q}} ~ f_{{p_b,r}}): {}/{}",
            n_full_cross_ok, n_pairs
        );
        println!(
            "  TEST 3 - Same-q cross-interlacing: {}/{}",
            n_same_q_ok, n_same_q_tested
        );
        println!(
            "  TEST 4 - Concatenation (A,B) interlacing: {}/{}",
            n_concat_ok, n_pairs
        );
        println!(
            "  TEST 5 - LR ratio ordering at source level: {}/{}",
            n_ratio_nondecr, n_ratio_tested
        );
        println!(
            "  TEST 6 - Shared roots between L^(p_a) and L^(p_b): {}/{}",
            n_shared_roots, n_shared_roots_tested
        );
        println!();
    }

    println!("========================================");
    println!("=== GRAND SUMMARY (n=5..{}) ===", max_n);
    println!("========================================");
    println!("Total consecutive pairs tested: {}", total_pairs);
    println!();
    println!(
        "TEST 1 - Individual chain interlacing:  A: {}/{}  B: {}/{}",
        chain_a_ok, total_pairs, chain_b_ok, total_pairs
    );
    println!(
        "TEST 2 - Full cross-interlacing:        {}/{}",
        full_cross_ok, total_pairs
    );
    println!(
        "TEST 3 - Same-q cross-interlacing:      {}/{}",
        same_q_ok, same_q_tested
    );
    println!(
        "TEST 4 - Concatenation interlacing:     {}/{}",
        concat_ok, total_pairs
    );
    println!(
        "TEST 5 - LR ratio at source level:      {}/{}",
        ratio_nondecr, ratio_tested
    );
    println!(
        "TEST 6 - Shared roots (L^pa, L^pb):     {}/{}",
        shared_roots_count, shared_roots_tested
    );
    println!();
    println!("INTERPRETATION:");
    println!("  If TEST 2 passes: cross-chain interlacing holds at source level,");
    println!("    so the staircase transform automatically preserves it.");
    println!("  If TEST 2 fails but TEST 4 passes: concatenation structure exists");
    println!("    but not full pairwise compatibility.");
    println!("  If TEST 5 passes: LR ordering at source level directly implies");
    println!("    LR ordering at output level (TN_2 for coefficient matrix).");
    println!("  If TEST 6 has many shared roots: interlacing may be forced by");
    println!("    common root structure from the staircase construction.");
}
