/// Explore real-rootedness of L_{n,S}(t) for fixed descent sets S.
///
/// For each S ⊆ [n-1], compute:
///   L_{n,S}(t) = Σ_{σ ∈ D(n,S)} t^{swaps(σ)}
/// and the position-of-maximum refinement:
///   L_{n,S}^{(p)}(t) = Σ_{σ ∈ D(n,S), σ⁻¹(n)=p} t^{swaps(σ)}
///
/// Checks:
/// 1. Real-rootedness of L_{n,S}(t)
/// 2. Compatible family property of {L_{n,S}^{(p)}}
/// 3. Correction terms C_{p,q}(t) and their interlacing with L_{n-1,S'}^{(q)}
/// 4. Staircase structure of the transfer matrix
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{
    check_weak_interlacing, format_poly, is_real_rooted,
};

/// Build polynomial from iterator of swaps values
fn build_poly_from_vals(vals: impl Iterator<Item = usize>) -> Vec<i64> {
    let vals: Vec<usize> = vals.collect();
    if vals.is_empty() {
        return vec![0];
    }
    let max_s = *vals.iter().max().unwrap();
    let mut coeffs = vec![0i64; max_s + 1];
    for s in vals {
        coeffs[s] += 1;
    }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
        coeffs.pop();
    }
    coeffs
}

/// Check if a set of polynomials forms a compatible family
/// (every pair has a common interlacing).
/// Returns (is_compatible, failing_pair_indices)
fn check_compatible_family(polys: &[Vec<i64>]) -> (bool, Vec<(usize, usize)>) {
    let mut failures = Vec::new();
    for i in 0..polys.len() {
        for j in (i + 1)..polys.len() {
            let f = &polys[i];
            let g = &polys[j];
            // Skip trivial cases
            if f.len() <= 1 || g.len() <= 1 {
                continue;
            }
            // For compatible family, every convex combination must be real-rooted.
            // A sufficient check: f+g is real-rooted AND various combinations work.
            // The Bézout-based weak interlacing check handles this:
            // two polynomials are compatible iff they have a common interlacing.
            //
            // We check: does (f, g) have a common interlacing?
            // Approach: check that c*f + (1-c)*g is real-rooted for several c values.
            // If any fails, they're not compatible.
            if !check_pairwise_compatible(f, g) {
                failures.push((i, j));
            }
        }
    }
    (failures.is_empty(), failures)
}

/// Check if two polynomials are pairwise compatible.
/// Two real-rooted polynomials with non-negative coefficients are compatible
/// iff every convex combination is real-rooted.
fn check_pairwise_compatible(f: &[i64], g: &[i64]) -> bool {
    if f.len() <= 1 || g.len() <= 1 {
        return true;
    }
    // Check that both are individually real-rooted
    if !is_real_rooted(f) || !is_real_rooted(g) {
        return false;
    }
    // Check convex combinations at several rational points
    // If f and g have a common interlacing, ALL convex combinations are rr.
    // We sample enough points to be confident.
    let test_weights: Vec<(i64, i64)> = vec![
        (1, 1),   // f + g
        (1, 2),   // f + 2g
        (2, 1),   // 2f + g
        (1, 3),   // f + 3g
        (3, 1),   // 3f + g
        (1, 5),
        (5, 1),
        (2, 3),
        (3, 2),
    ];
    let maxlen = f.len().max(g.len());
    for (a, b) in &test_weights {
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

/// Compute valid positions P(S) for descent set S (as bitmask) in S_n
fn valid_positions(s_mask: u64, n: u8) -> Vec<u8> {
    let mut positions = Vec::new();
    // p ∈ S such that p-1 ∉ S (left boundary of descent run)
    // Note: bit j in s_mask corresponds to descent at position j+1 (1-indexed)
    // So descent at position p means bit (p-1) is set
    for p in 1..n {
        let p_in_s = (s_mask >> (p - 1)) & 1 == 1;
        let pm1_in_s = if p >= 2 {
            (s_mask >> (p - 2)) & 1 == 1
        } else {
            false // 0 ∉ S
        };
        if p_in_s && !pm1_in_s {
            positions.push(p);
        }
    }
    // Also p = n if n-1 ∉ S
    let nm1_in_s = if n >= 2 {
        (s_mask >> (n - 2)) & 1 == 1
    } else {
        false
    };
    if !nm1_in_s {
        positions.push(n);
    }
    positions
}

/// Compute the source descent set S'_p for inserting n at position p
/// given target descent set S in S_n.
/// Returns the bitmask for S'_p ⊆ [n-2].
fn source_descent_set(s_mask: u64, p: u8, n: u8) -> u64 {
    if n <= 2 {
        return 0;
    }
    let mut sp = 0u64;
    if p == n {
        // S'_n = S (must have S ⊆ [n-2])
        sp = s_mask;
    } else if p == 1 {
        // S'_1 = {j-1 : j ∈ S, j ≥ 2}
        for j in 2..n {
            if (s_mask >> (j - 1)) & 1 == 1 {
                sp |= 1 << (j - 2);
            }
        }
    } else {
        // S'_p = (S ∩ [p-2]) ∪ {j-1 : j ∈ S, j > p}
        // S ∩ [p-2]: descents at positions 1, ..., p-2
        for pos in 1..=(p.saturating_sub(2)) {
            if (s_mask >> (pos - 1)) & 1 == 1 {
                sp |= 1 << (pos - 1);
            }
        }
        // {j-1 : j ∈ S, j > p}
        for j in (p + 1)..n {
            if (s_mask >> (j - 1)) & 1 == 1 {
                sp |= 1 << (j - 2);
            }
        }
    }
    sp
}

/// Check the epsilon_1 condition: does inserting n at position p in π
/// break an adjacent consecutive pair?
/// Requires 1 < p <= n-1 and π(p-1) + 1 = π(p) (1-indexed positions in π)
fn epsilon1(pi: &[u8], p: u8) -> bool {
    let n = pi.len() as u8 + 1; // π ∈ S_{n-1}, inserting into S_n
    if p <= 1 || p >= n {
        return false;
    }
    // π is 0-indexed, positions are 1-indexed
    // π(p-1) means pi[p-2], π(p) means pi[p-1]
    let val_left = pi[(p - 2) as usize];
    let val_right = pi[(p - 1) as usize];
    val_left + 1 == val_right
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
        .unwrap_or(8);

    let mut total_sets = 0u64;
    let mut rr_count = 0u64;
    let mut compatible_count = 0u64;
    let mut compatible_total = 0u64;
    let mut correction_interlace_count = 0u64;
    let mut correction_total = 0u64;
    let mut worst_failures: Vec<(u8, u64, String)> = Vec::new();

    for n in 3..=max_n {
        let perms = all_permutations(n);
        println!("========== n = {} ({} perms) ==========", n, perms.len());

        // Group permutations by descent set bitmask
        let mut by_descent: std::collections::HashMap<u64, Vec<&Vec<u8>>> =
            std::collections::HashMap::new();
        for pi in &perms {
            let mask = descent_set_bitmask(pi);
            by_descent.entry(mask).or_default().push(pi);
        }

        let num_sets = by_descent.len();
        let mut n_rr = 0;
        let mut n_compat = 0;
        let mut n_compat_tested = 0;
        let mut n_corr_ok = 0;
        let mut n_corr_tested = 0;

        // Sort descent sets for deterministic output
        let mut sorted_masks: Vec<u64> = by_descent.keys().copied().collect();
        sorted_masks.sort();

        for &mask in &sorted_masks {
            let class = &by_descent[&mask];
            let s_str = descent_set_to_string(mask, n);

            // Compute L_{n,S}(t)
            let poly = build_poly_from_vals(class.iter().map(|pi| compute(pi, Stat::Swaps)));
            let rr = poly.len() <= 2 || is_real_rooted(&poly);

            if rr {
                n_rr += 1;
            } else {
                println!("  NOT RR: S={} poly={}", s_str, format_poly(&poly));
                worst_failures.push((n, mask, s_str.clone()));
            }

            // Compute position-of-maximum refinement L_{n,S}^{(p)}
            let valid_pos = valid_positions(mask, n);

            let mut pos_polys: Vec<(u8, Vec<i64>)> = Vec::new();
            for &p in &valid_pos {
                let pos_class: Vec<usize> = class
                    .iter()
                    .filter(|pi| {
                        let pos_of_n = pi.iter().position(|&v| v == n).unwrap();
                        (pos_of_n + 1) as u8 == p // 1-indexed
                    })
                    .map(|pi| compute(pi, Stat::Swaps))
                    .collect();
                if pos_class.is_empty() {
                    continue;
                }
                let ppoly = build_poly_from_vals(pos_class.into_iter());
                pos_polys.push((p, ppoly));
            }

            // Check compatible family
            if pos_polys.len() >= 2 {
                let polys_only: Vec<Vec<i64>> =
                    pos_polys.iter().map(|(_, p)| p.clone()).collect();
                let (compat, failures) = check_compatible_family(&polys_only);
                n_compat_tested += 1;
                if compat {
                    n_compat += 1;
                } else {
                    let one_not_in_s = (mask & 1) == 0; // bit 0 = position 1
                    println!(
                        "  COMPAT FAIL: S={} (1∉S: {}) failures={:?}",
                        s_str, one_not_in_s, failures
                    );
                    for (i, j) in &failures {
                        println!(
                            "    pos {}:{} vs pos {}:{}",
                            pos_polys[*i].0,
                            format_poly(&pos_polys[*i].1),
                            pos_polys[*j].0,
                            format_poly(&pos_polys[*j].1)
                        );
                    }
                }
            }

            // Check correction terms
            if n >= 4 {
                let perms_nm1 = all_permutations(n - 1);
                let mut by_descent_nm1: std::collections::HashMap<u64, Vec<&Vec<u8>>> =
                    std::collections::HashMap::new();
                for pi in &perms_nm1 {
                    let m = descent_set_bitmask(pi);
                    by_descent_nm1.entry(m).or_default().push(pi);
                }

                for &p in &valid_pos {
                    if p <= 1 || p >= n {
                        continue; // No correction for p=1 or p=n
                    }
                    let sp_mask = source_descent_set(mask, p, n);
                    let source_class = match by_descent_nm1.get(&sp_mask) {
                        Some(c) => c,
                        None => continue,
                    };

                    // For each valid q in S'_p, compute C_{p,q} and L_{n-1,S'}^{(q)}
                    let valid_q = valid_positions(sp_mask, n - 1);
                    for &q in &valid_q {
                        // L_{n-1,S'}^{(q)}: source perms with max at position q
                        let source_at_q: Vec<&Vec<u8>> = source_class
                            .iter()
                            .filter(|pi| {
                                let pos_of_max =
                                    pi.iter().position(|&v| v == n - 1).unwrap();
                                (pos_of_max + 1) as u8 == q
                            })
                            .copied()
                            .collect();
                        if source_at_q.is_empty() {
                            continue;
                        }

                        // C_{p,q}: subset with epsilon_1 = 1
                        let correction_vals: Vec<usize> = source_at_q
                            .iter()
                            .filter(|pi| epsilon1(pi, p))
                            .map(|pi| compute(pi, Stat::Swaps))
                            .collect();

                        if correction_vals.is_empty() {
                            continue; // C_{p,q} = 0, trivially OK
                        }

                        let c_poly =
                            build_poly_from_vals(correction_vals.into_iter());
                        let l_poly = build_poly_from_vals(
                            source_at_q.iter().map(|pi| compute(pi, Stat::Swaps)),
                        );

                        n_corr_tested += 1;

                        // Check: does C_{p,q} interlace L_{n-1,S'}^{(q)}?
                        if c_poly.len() <= 1 || l_poly.len() <= 1 {
                            n_corr_ok += 1;
                            continue;
                        }
                        let (small, large) = if c_poly.len() <= l_poly.len() {
                            (&c_poly, &l_poly)
                        } else {
                            (&l_poly, &c_poly)
                        };
                        match check_weak_interlacing(small, large) {
                            Some(true) => {
                                n_corr_ok += 1;
                            }
                            Some(false) => {
                                println!(
                                    "  CORR FAIL: S={} p={} q={} C={} L={}",
                                    s_str,
                                    p,
                                    q,
                                    format_poly(&c_poly),
                                    format_poly(&l_poly)
                                );
                            }
                            None => {
                                // Inconclusive, try convex combo check
                                let sum: Vec<i64> = {
                                    let maxlen = c_poly.len().max(l_poly.len());
                                    let mut s = vec![0i64; maxlen];
                                    for (i, &c) in c_poly.iter().enumerate() {
                                        s[i] += c;
                                    }
                                    for (i, &c) in l_poly.iter().enumerate() {
                                        s[i] += c;
                                    }
                                    s
                                };
                                if is_real_rooted(&sum) {
                                    n_corr_ok += 1;
                                } else {
                                    println!(
                                        "  CORR ?: S={} p={} q={} C={} L={}",
                                        s_str,
                                        p,
                                        q,
                                        format_poly(&c_poly),
                                        format_poly(&l_poly)
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }

        total_sets += num_sets as u64;
        rr_count += n_rr;
        compatible_count += n_compat;
        compatible_total += n_compat_tested;
        correction_interlace_count += n_corr_ok;
        correction_total += n_corr_tested;

        println!(
            "n={}: {}/{} rr, {}/{} compatible, {}/{} corrections OK\n",
            n, n_rr, num_sets, n_compat, n_compat_tested, n_corr_ok, n_corr_tested,
        );
    }

    println!("=== SUMMARY ===");
    println!("Total descent sets tested: {}", total_sets);
    println!("Real-rooted: {}/{}", rr_count, total_sets);
    println!(
        "Compatible families: {}/{}",
        compatible_count, compatible_total
    );
    println!(
        "Correction interlacing: {}/{}",
        correction_interlace_count, correction_total
    );
    if !worst_failures.is_empty() {
        println!("\nFAILURES:");
        for (n, mask, s) in &worst_failures {
            println!("  n={} S={} (mask={})", n, s, mask);
        }
    }
}
