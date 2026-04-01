/// Explore induction on |S| for real-rootedness of L_{n,S}(t).
///
/// For a descent set S ⊆ [n-1] with |S| >= 1, consider removing one element
/// i from S to get S' = S \ {i}. Both D(n,S) and D(n,S') are sets of
/// permutations of [n], giving polynomials L_{n,S}(t) and L_{n,S'}(t).
///
/// Questions explored:
/// 1. Does L_{n,S'} interlace L_{n,S} (or vice versa)?
/// 2. Is there a preferred element to remove (min, max, or other)?
/// 3. Can we write L_{n,S} as a "transfer" from L_{n,S'}?
/// 4. Is there a polynomial matrix recurrence S -> S ∪ {i}?
///
/// We restrict to S with 1 ∉ S (the "no leading descent" family).
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly, is_real_rooted};
use std::collections::BTreeMap;

// ============================================================
// Helpers
// ============================================================

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

/// Pretty-print a descent set given as a bitmask.
fn descent_set_str(mask: u64, n: u8) -> String {
    let positions: Vec<String> = (0..n - 1)
        .filter(|&i| mask & (1 << i) != 0)
        .map(|i| (i + 1).to_string())
        .collect();
    if positions.is_empty() {
        "{}".to_string()
    } else {
        format!("{{{}}}", positions.join(","))
    }
}

/// Elements of S (1-indexed positions) from a bitmask.
fn descent_set_elements(mask: u64, n: u8) -> Vec<u8> {
    (1..n).filter(|&i| mask & (1 << (i - 1)) != 0).collect()
}

/// Popcount of a bitmask.
#[allow(dead_code)]
fn popcount(mask: u64) -> u32 {
    mask.count_ones()
}

/// Remove element i (1-indexed) from descent set bitmask.
fn remove_element(mask: u64, i: u8) -> u64 {
    mask & !(1u64 << (i - 1))
}

/// Check if f weakly interlaces g or g weakly interlaces f.
/// Returns: (f_leq_g, g_leq_f) where each is true if that direction holds.
fn check_both_interlacings(f: &[i64], g: &[i64]) -> (bool, bool) {
    if f.iter().all(|&c| c == 0) || g.iter().all(|&c| c == 0) {
        return (true, true);
    }
    let f_leq_g = check_weak_interlacing(f, g) == Some(true);
    let g_leq_f = check_weak_interlacing(g, f) == Some(true);
    (f_leq_g, g_leq_f)
}

/// Polynomial subtraction: f - g.
fn poly_sub(f: &[i64], g: &[i64]) -> Vec<i64> {
    let len = f.len().max(g.len());
    let mut result = vec![0i64; len];
    for (i, &c) in f.iter().enumerate() {
        result[i] += c;
    }
    for (i, &c) in g.iter().enumerate() {
        result[i] -= c;
    }
    while result.len() > 1 && *result.last().unwrap() == 0 {
        result.pop();
    }
    result
}

/// Polynomial addition: f + g.
#[allow(dead_code)]
fn poly_add(f: &[i64], g: &[i64]) -> Vec<i64> {
    let len = f.len().max(g.len());
    let mut result = vec![0i64; len];
    for (i, &c) in f.iter().enumerate() {
        result[i] += c;
    }
    for (i, &c) in g.iter().enumerate() {
        result[i] += c;
    }
    while result.len() > 1 && *result.last().unwrap() == 0 {
        result.pop();
    }
    result
}

/// Multiply polynomial by t (shift coefficients right by 1).
#[allow(dead_code)]
fn poly_mul_t(f: &[i64]) -> Vec<i64> {
    let mut result = vec![0i64; f.len() + 1];
    for (i, &c) in f.iter().enumerate() {
        result[i + 1] = c;
    }
    result
}

/// Degree of a polynomial (ignoring trailing zeros).
fn poly_deg(f: &[i64]) -> Option<usize> {
    for i in (0..f.len()).rev() {
        if f[i] != 0 {
            return Some(i);
        }
    }
    None
}

/// Check if a polynomial is the zero polynomial.
fn is_zero_poly(f: &[i64]) -> bool {
    f.iter().all(|&c| c == 0)
}

// ============================================================
// Part 1: Interlacing between L_{n,S} and L_{n,S\{i}}
// ============================================================

fn part1_interlacing(n: u8, polys: &BTreeMap<u64, Vec<i64>>) {
    println!("--- n = {} ---", n);

    let mut total_pairs = 0u32;
    let mut sp_leq_s = 0u32; // L_{n,S'} ≪ L_{n,S}
    let mut s_leq_sp = 0u32; // L_{n,S} ≪ L_{n,S'}
    let mut both_dir = 0u32; // interlacing in both directions
    let mut neither = 0u32; // no interlacing

    // Track which removal position tends to give interlacing
    let mut by_removal_pos: BTreeMap<String, (u32, u32)> = BTreeMap::new(); // "min"/"max"/"mid" -> (total, success)

    let mut fail_examples: Vec<String> = Vec::new();

    for (&s_mask, s_poly) in polys {
        let elems = descent_set_elements(s_mask, n);
        if elems.is_empty() {
            continue; // |S| = 0, nothing to remove
        }
        // Skip if 1 ∈ S
        if s_mask & 1 != 0 {
            continue;
        }

        let s_rr = s_poly.len() <= 2 || is_real_rooted(s_poly);

        for &i in &elems {
            let sp_mask = remove_element(s_mask, i);
            let sp_poly = match polys.get(&sp_mask) {
                Some(p) => p,
                None => continue,
            };

            let sp_rr = sp_poly.len() <= 2 || is_real_rooted(sp_poly);

            let (fwd, bwd) = check_both_interlacings(sp_poly, s_poly);
            total_pairs += 1;

            if fwd {
                sp_leq_s += 1;
            }
            if bwd {
                s_leq_sp += 1;
            }
            if fwd && bwd {
                both_dir += 1;
            }
            if !fwd && !bwd {
                neither += 1;
            }

            // Classify removal position
            let pos_label = if i == *elems.first().unwrap() && elems.len() > 1 {
                "min"
            } else if i == *elems.last().unwrap() && elems.len() > 1 {
                "max"
            } else if elems.len() <= 1 {
                "only"
            } else {
                "mid"
            };
            let entry = by_removal_pos
                .entry(pos_label.to_string())
                .or_insert((0, 0));
            entry.0 += 1;
            if fwd || bwd {
                entry.1 += 1;
            }

            // Collect failures
            if !fwd && !bwd && fail_examples.len() < 3 {
                fail_examples.push(format!(
                    "S={} rm {} -> S'={}  rr(S)={} rr(S')={}  L_S={} L_S'={}",
                    descent_set_str(s_mask, n),
                    i,
                    descent_set_str(sp_mask, n),
                    s_rr,
                    sp_rr,
                    format_poly(s_poly),
                    format_poly(sp_poly),
                ));
            }
        }
    }

    println!(
        "  pairs: {}  S'≪S: {}  S≪S': {}  both: {}  neither: {}",
        total_pairs, sp_leq_s, s_leq_sp, both_dir, neither
    );
    println!("  by removal position:");
    for (label, (tot, succ)) in &by_removal_pos {
        println!("    {}: {}/{} have some interlacing", label, succ, tot);
    }
    for ex in &fail_examples {
        println!("  FAIL: {}", ex);
    }
}

// ============================================================
// Part 2: Transfer from L_{n,S'} to L_{n,S} where
//         S' = S\{max(S)} or S' = S\{min(S)}
// ============================================================

fn part2_transfer(n: u8, polys: &BTreeMap<u64, Vec<i64>>, _swap_vals: &BTreeMap<u64, Vec<usize>>) {
    println!("--- n = {} ---", n);

    // For each S with |S| >= 1, 1 ∉ S, try to express
    // L_{n,S} = a(t) * L_{n,S'} + b(t) * R(t)
    // where S' = S\{max(S)} or S\{min(S)}.
    // First, just check the difference polynomial D = L_{n,S} - L_{n,S'}.

    let mut remove_max_results: Vec<String> = Vec::new();
    let mut remove_min_results: Vec<String> = Vec::new();

    for (&s_mask, s_poly) in polys {
        let elems = descent_set_elements(s_mask, n);
        if elems.is_empty() || s_mask & 1 != 0 {
            continue;
        }
        if elems.len() < 2 {
            continue; // need at least 2 elements to distinguish min/max
        }

        // Remove max(S)
        let max_elem = *elems.last().unwrap();
        let sp_max_mask = remove_element(s_mask, max_elem);
        if let Some(sp_max_poly) = polys.get(&sp_max_mask) {
            let diff = poly_sub(s_poly, sp_max_poly);
            let diff_rr = diff.len() <= 2 || is_zero_poly(&diff) || is_real_rooted(&diff);

            // Check: is diff a multiple of t? (i.e., diff[0] == 0)
            let diff_div_t = diff.len() > 1 && diff[0] == 0;

            // Check interlacing of diff with s_poly
            let (diff_leq_s, _) = if !is_zero_poly(&diff) {
                check_both_interlacings(&diff, s_poly)
            } else {
                (true, true)
            };

            // Check sign pattern of diff coefficients
            let all_nonneg = diff.iter().all(|&c| c >= 0);
            let all_nonpos = diff.iter().all(|&c| c <= 0);

            if remove_max_results.len() < 5 || !diff_rr || !all_nonneg {
                remove_max_results.push(format!(
                    "  S={} rm max({})  D=L_S-L_S' rr={} nonneg={} nonpos={} div_t={} D≪S={}\n    L_S={}\n    L_S'={}\n    D={}",
                    descent_set_str(s_mask, n),
                    max_elem,
                    diff_rr,
                    all_nonneg,
                    all_nonpos,
                    diff_div_t,
                    diff_leq_s,
                    format_poly(s_poly),
                    format_poly(sp_max_poly),
                    format_poly(&diff),
                ));
            }
        }

        // Remove min(S)
        let min_elem = *elems.first().unwrap();
        let sp_min_mask = remove_element(s_mask, min_elem);
        if let Some(sp_min_poly) = polys.get(&sp_min_mask) {
            let diff = poly_sub(s_poly, sp_min_poly);
            let diff_rr = diff.len() <= 2 || is_zero_poly(&diff) || is_real_rooted(&diff);
            let diff_div_t = diff.len() > 1 && diff[0] == 0;
            let (diff_leq_s, _) = if !is_zero_poly(&diff) {
                check_both_interlacings(&diff, s_poly)
            } else {
                (true, true)
            };
            let all_nonneg = diff.iter().all(|&c| c >= 0);
            let all_nonpos = diff.iter().all(|&c| c <= 0);

            if remove_min_results.len() < 5 || !diff_rr || !all_nonneg {
                remove_min_results.push(format!(
                    "  S={} rm min({})  D=L_S-L_S' rr={} nonneg={} nonpos={} div_t={} D≪S={}\n    L_S={}\n    L_S'={}\n    D={}",
                    descent_set_str(s_mask, n),
                    min_elem,
                    diff_rr,
                    all_nonneg,
                    all_nonpos,
                    diff_div_t,
                    diff_leq_s,
                    format_poly(s_poly),
                    format_poly(sp_min_poly),
                    format_poly(&diff),
                ));
            }
        }
    }

    // Summarize
    println!("  Remove max(S) examples:");
    for r in &remove_max_results {
        println!("{}", r);
    }
    println!("  Remove min(S) examples:");
    for r in &remove_min_results {
        println!("{}", r);
    }
}

// ============================================================
// Part 3: Matrix recurrence S -> S ∪ {i}
// ============================================================

/// For a given S and S' = S ∪ {i}, refine L_{n,S} and L_{n,S'} by position
/// of some value (e.g., position of n, or position of i, etc.) and look
/// for a polynomial matrix connecting the refined vectors.
fn part3_matrix_recurrence(n: u8, all_perms: &[Vec<u8>], polys: &BTreeMap<u64, Vec<i64>>) {
    println!("--- n = {} ---", n);

    // For each S with |S|>=1 and 1 ∉ S, and each i ∈ S,
    // let S' = S\{i}. Refine both L_{n,S} and L_{n,S'} by
    // position of n. Then check if there is a polynomial
    // (column) operation connecting them.

    // Group permutations by (descent_set, pos_of_n)
    let mut by_des_pos: BTreeMap<(u64, usize), Vec<usize>> = BTreeMap::new();
    for pi in all_perms {
        let ds = descent_set_bitmask(pi);
        let pos_n = pi.iter().position(|&v| v == n).unwrap() + 1;
        let sw = compute(pi, Stat::Swaps);
        by_des_pos.entry((ds, pos_n)).or_default().push(sw);
    }

    let mut found_pattern = 0u32;
    let mut total_checked = 0u32;
    let mut examples: Vec<String> = Vec::new();

    for (&s_mask, _s_poly) in polys {
        let elems = descent_set_elements(s_mask, n);
        if elems.is_empty() || s_mask & 1 != 0 {
            continue;
        }

        // Collect position-refined polys for S
        let mut s_pos_polys: BTreeMap<usize, Vec<i64>> = BTreeMap::new();
        for pos in 1..=n as usize {
            if let Some(vals) = by_des_pos.get(&(s_mask, pos)) {
                s_pos_polys.insert(pos, build_poly(vals));
            }
        }

        for &i in &elems {
            let sp_mask = remove_element(s_mask, i);

            // Collect position-refined polys for S'
            let mut sp_pos_polys: BTreeMap<usize, Vec<i64>> = BTreeMap::new();
            for pos in 1..=n as usize {
                if let Some(vals) = by_des_pos.get(&(sp_mask, pos)) {
                    sp_pos_polys.insert(pos, build_poly(vals));
                }
            }

            total_checked += 1;

            // Check: for each position p, is L_{n,S}^{(p)} - L_{n,S'}^{(p)}
            // a simple polynomial (monomial, multiple of t, etc.)?
            let mut all_diffs_nonneg = true;
            let mut all_diffs_rr = true;
            let mut diff_descriptions: Vec<String> = Vec::new();

            let all_positions: Vec<usize> = {
                let mut ps: Vec<usize> = s_pos_polys
                    .keys()
                    .chain(sp_pos_polys.keys())
                    .cloned()
                    .collect();
                ps.sort();
                ps.dedup();
                ps
            };

            for &pos in &all_positions {
                let s_p = s_pos_polys.get(&pos).cloned().unwrap_or(vec![0]);
                let sp_p = sp_pos_polys.get(&pos).cloned().unwrap_or(vec![0]);
                let diff = poly_sub(&s_p, &sp_p);

                if !diff.iter().all(|&c| c >= 0) {
                    all_diffs_nonneg = false;
                }
                if diff.len() > 2 && !is_zero_poly(&diff) && !is_real_rooted(&diff) {
                    all_diffs_rr = false;
                }

                let deg = poly_deg(&diff);
                let nonzero_count = diff.iter().filter(|&&c| c != 0).count();
                diff_descriptions.push(format!(
                    "p={}: deg={:?} #nz={} D={}",
                    pos,
                    deg,
                    nonzero_count,
                    format_poly(&diff),
                ));
            }

            if all_diffs_nonneg && all_diffs_rr {
                found_pattern += 1;
            }

            if examples.len() < 4 {
                examples.push(format!(
                    "  S={} rm {} -> S'={}  all_diffs_nonneg={} all_diffs_rr={}",
                    descent_set_str(s_mask, n),
                    i,
                    descent_set_str(sp_mask, n),
                    all_diffs_nonneg,
                    all_diffs_rr,
                ));
                for d in &diff_descriptions {
                    examples.push(format!("    {}", d));
                }
            }
        }
    }

    println!(
        "  checked {} (S, removal) pairs: {} have all position-diffs nonneg & rr",
        total_checked, found_pattern,
    );
    for ex in &examples {
        println!("{}", ex);
    }
}

// ============================================================
// Part 4: Pairwise compatibility of {L_{n,S'} : i ∈ S}
// ============================================================

fn part4_pairwise_compatibility(n: u8, polys: &BTreeMap<u64, Vec<i64>>) {
    println!("--- n = {} ---", n);

    let mut total_sets = 0u32;
    let mut all_compat = 0u32;
    let mut some_fail = 0u32;
    let mut fail_examples: Vec<String> = Vec::new();

    for (&s_mask, s_poly) in polys {
        let elems = descent_set_elements(s_mask, n);
        if elems.len() < 2 || s_mask & 1 != 0 {
            continue;
        }

        total_sets += 1;

        // Collect all L_{n,S\{i}} for i ∈ S
        let mut child_polys: Vec<(u8, Vec<i64>)> = Vec::new();
        for &i in &elems {
            let sp_mask = remove_element(s_mask, i);
            if let Some(p) = polys.get(&sp_mask) {
                child_polys.push((i, p.clone()));
            }
        }

        // Check pairwise weak interlacing (in either direction)
        let mut compatible = true;
        let mut fail_pair = None;
        for a in 0..child_polys.len() {
            for b in a + 1..child_polys.len() {
                let (fwd, bwd) = check_both_interlacings(&child_polys[a].1, &child_polys[b].1);
                if !fwd && !bwd {
                    compatible = false;
                    fail_pair = Some((child_polys[a].0, child_polys[b].0));
                    break;
                }
            }
            if !compatible {
                break;
            }
        }

        if compatible {
            all_compat += 1;
        } else {
            some_fail += 1;
            if fail_examples.len() < 3 {
                let (i1, i2) = fail_pair.unwrap();
                fail_examples.push(format!(
                    "S={} |S|={}: L_{{S\\{{{}}}}} and L_{{S\\{{{}}}}} not compatible  L_S={}",
                    descent_set_str(s_mask, n),
                    elems.len(),
                    i1,
                    i2,
                    format_poly(s_poly),
                ));
            }
        }
    }

    println!(
        "  sets with |S|>=2: {}  all children compatible: {}  some fail: {}",
        total_sets, all_compat, some_fail,
    );
    for ex in &fail_examples {
        println!("  FAIL: {}", ex);
    }
}

// ============================================================
// Part 5: Degree and coefficient comparison
// ============================================================

fn part5_degree_comparison(n: u8, polys: &BTreeMap<u64, Vec<i64>>) {
    println!("--- n = {} ---", n);

    // For S -> S\{i}, compare degrees and leading coefficients.
    // Also check: does L_{n,S} always have >= degree than L_{n,S'}?

    let mut always_geq_deg = true;
    let mut always_leq_deg = true;
    let mut deg_table: Vec<String> = Vec::new();

    for (&s_mask, s_poly) in polys {
        let elems = descent_set_elements(s_mask, n);
        if elems.is_empty() || s_mask & 1 != 0 {
            continue;
        }

        let s_deg = poly_deg(s_poly);

        for &i in &elems {
            let sp_mask = remove_element(s_mask, i);
            if let Some(sp_poly) = polys.get(&sp_mask) {
                let sp_deg = poly_deg(sp_poly);

                match (s_deg, sp_deg) {
                    (Some(d1), Some(d2)) => {
                        if d1 < d2 {
                            always_geq_deg = false;
                        }
                        if d1 > d2 {
                            always_leq_deg = false;
                        }
                    }
                    _ => {}
                }

                if deg_table.len() < 8 {
                    let diff = poly_sub(s_poly, sp_poly);
                    deg_table.push(format!(
                        "  S={} rm {} -> S'={}  deg(S)={:?} deg(S')={:?} diff_nonneg={} |D(n,S)|={} |D(n,S')|={}",
                        descent_set_str(s_mask, n),
                        i,
                        descent_set_str(sp_mask, n),
                        s_deg,
                        sp_deg,
                        diff.iter().all(|&c| c >= 0),
                        s_poly.iter().sum::<i64>(),
                        sp_poly.iter().sum::<i64>(),
                    ));
                }
            }
        }
    }

    println!(
        "  deg(L_S) >= deg(L_S') always: {}  deg(L_S) <= deg(L_S') always: {}",
        always_geq_deg, always_leq_deg,
    );
    for row in &deg_table {
        println!("{}", row);
    }
}

// ============================================================
// Part 6: Detailed view of small cases
// ============================================================

fn part6_small_cases(n: u8, polys: &BTreeMap<u64, Vec<i64>>) {
    println!("--- n = {} (detailed interlacing table) ---", n);

    for (&s_mask, s_poly) in polys {
        let elems = descent_set_elements(s_mask, n);
        if elems.is_empty() || s_mask & 1 != 0 {
            continue;
        }

        let s_rr = s_poly.len() <= 2 || is_real_rooted(s_poly);

        println!(
            "  S={}  L_S={}  rr={}",
            descent_set_str(s_mask, n),
            format_poly(s_poly),
            s_rr,
        );

        for &i in &elems {
            let sp_mask = remove_element(s_mask, i);
            if let Some(sp_poly) = polys.get(&sp_mask) {
                let sp_rr = sp_poly.len() <= 2 || is_real_rooted(sp_poly);
                let (fwd, bwd) = check_both_interlacings(sp_poly, s_poly);
                let diff = poly_sub(s_poly, sp_poly);
                let diff_rr = diff.len() <= 2 || is_zero_poly(&diff) || is_real_rooted(&diff);

                let interl_str = match (fwd, bwd) {
                    (true, true) => "BOTH",
                    (true, false) => "S'≪S",
                    (false, true) => "S≪S'",
                    (false, false) => "NONE",
                };

                println!(
                    "    rm {} -> S'={}  L_S'={}  rr={}  interlace={}  D={} rr(D)={}",
                    i,
                    descent_set_str(sp_mask, n),
                    format_poly(sp_poly),
                    sp_rr,
                    interl_str,
                    format_poly(&diff),
                    diff_rr,
                );
            }
        }
    }
}

// ============================================================
// Main
// ============================================================

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    let min_n: u8 = 4;

    println!("================================================================");
    println!("  Induction on |S|: interlacing between L_{{n,S}} and L_{{n,S\\{{i}}}}");
    println!("  Restricting to S with 1 not in S");
    println!("  n = {} .. {}", min_n, max_n);
    println!("================================================================\n");

    // Precompute all data
    let mut all_data: BTreeMap<
        u8,
        (
            Vec<Vec<u8>>,
            BTreeMap<u64, Vec<i64>>,
            BTreeMap<u64, Vec<usize>>,
        ),
    > = BTreeMap::new();

    for n in min_n..=max_n {
        let all = all_permutations(n);
        let mut by_des: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
        for pi in &all {
            let ds = descent_set_bitmask(pi);
            let sw = compute(pi, Stat::Swaps);
            by_des.entry(ds).or_default().push(sw);
        }
        let polys: BTreeMap<u64, Vec<i64>> = by_des
            .iter()
            .map(|(&ds, vals)| (ds, build_poly(vals)))
            .collect();
        all_data.insert(n, (all, polys, by_des));
    }

    // ============================================================
    // Part 0: Summary of real-rootedness
    // ============================================================
    println!("════════════════════════════════════════════════════════════════");
    println!("  Part 0: Real-rootedness summary for L_{{n,S}} with 1 not in S");
    println!("════════════════════════════════════════════════════════════════\n");

    for n in min_n..=max_n {
        let (_, polys, _) = all_data.get(&n).unwrap();
        let mut total = 0u32;
        let mut rr = 0u32;
        let mut not_rr_list: Vec<String> = Vec::new();
        for (&ds, poly) in polys {
            if ds & 1 != 0 {
                continue;
            }
            total += 1;
            if poly.len() <= 2 || is_real_rooted(poly) {
                rr += 1;
            } else {
                not_rr_list.push(format!(
                    "S={} L={}",
                    descent_set_str(ds, n),
                    format_poly(poly),
                ));
            }
        }
        println!("n={}: {}/{} real-rooted", n, rr, total);
        for ex in &not_rr_list {
            println!("  NOT RR: {}", ex);
        }
    }

    // ============================================================
    // Part 1: Interlacing between L_{n,S} and L_{n,S\{i}}
    // ============================================================
    println!("\n════════════════════════════════════════════════════════════════");
    println!("  Part 1: Interlacing between L_{{n,S}} and L_{{n,S\\{{i}}}}");
    println!("════════════════════════════════════════════════════════════════\n");

    for n in min_n..=max_n {
        let (_, polys, _) = all_data.get(&n).unwrap();
        part1_interlacing(n, polys);
        println!();
    }

    // ============================================================
    // Part 2: Transfer structure (remove max or min)
    // ============================================================
    println!("════════════════════════════════════════════════════════════════");
    println!("  Part 2: Difference D = L_{{n,S}} - L_{{n,S\\{{max(S)}}}}");
    println!("          and D = L_{{n,S}} - L_{{n,S\\{{min(S)}}}}");
    println!("════════════════════════════════════════════════════════════════\n");

    for n in min_n..=max_n.min(7) {
        let (_, polys, swap_vals) = all_data.get(&n).unwrap();
        part2_transfer(n, polys, swap_vals);
        println!();
    }

    // ============================================================
    // Part 3: Position-refined matrix recurrence
    // ============================================================
    println!("════════════════════════════════════════════════════════════════");
    println!("  Part 3: Position-of-max refined differences S -> S\\{{i}}");
    println!("════════════════════════════════════════════════════════════════\n");

    for n in min_n..=max_n.min(7) {
        let (all_perms, polys, _) = all_data.get(&n).unwrap();
        part3_matrix_recurrence(n, all_perms, polys);
        println!();
    }

    // ============================================================
    // Part 4: Pairwise compatibility of children
    // ============================================================
    println!("════════════════════════════════════════════════════════════════");
    println!("  Part 4: Pairwise compatibility of {{L_{{n,S\\{{i}}}} : i in S}}");
    println!("════════════════════════════════════════════════════════════════\n");

    for n in min_n..=max_n {
        let (_, polys, _) = all_data.get(&n).unwrap();
        part4_pairwise_compatibility(n, polys);
        println!();
    }

    // ============================================================
    // Part 5: Degree and size comparison
    // ============================================================
    println!("════════════════════════════════════════════════════════════════");
    println!("  Part 5: Degree and coefficient comparison S vs S\\{{i}}");
    println!("════════════════════════════════════════════════════════════════\n");

    for n in min_n..=max_n.min(7) {
        let (_, polys, _) = all_data.get(&n).unwrap();
        part5_degree_comparison(n, polys);
        println!();
    }

    // ============================================================
    // Part 6: Detailed view for small n
    // ============================================================
    println!("════════════════════════════════════════════════════════════════");
    println!("  Part 6: Detailed interlacing table for small n");
    println!("════════════════════════════════════════════════════════════════\n");

    for n in min_n..=max_n.min(6) {
        let (_, polys, _) = all_data.get(&n).unwrap();
        part6_small_cases(n, polys);
        println!();
    }

    println!("Done.");
}
