/// Search for palindromic AND real-rooted backtrack polynomial sequences.
///
/// Focuses on des, exc, peak statistics (which tend to be real-rooted).
/// Explores: single length-3 patterns, named families, two-pattern avoidance,
/// length-4 patterns, and Catalan search orders.
use combpoly::permutation::{
    all_permutations, avoiding_permutations, backtrack_stat_poly, backtrack_stat_poly_subset,
    contains_pattern, is_alternating, is_derangement, is_grassmann, is_involution, is_separable,
    is_simsun, is_vexillary, parse_sequence,
};
use combpoly::statistics::Stat;
use polynomial_tools::real_rootedness::{self as polynomial, is_real_rooted};
use std::env;

fn is_palindromic(p: &[i64]) -> bool {
    let n = p.len();
    if n <= 1 {
        return true;
    }
    let end = p.iter().rposition(|&x| x != 0).unwrap_or(0);
    let start = p.iter().position(|&x| x != 0).unwrap_or(0);
    let trimmed = &p[start..=end];
    let m = trimmed.len();
    (0..m / 2).all(|i| trimmed[i] == trimmed[m - 1 - i])
}

fn format_poly(p: &[i64]) -> String {
    p.iter()
        .map(|c| c.to_string())
        .collect::<Vec<_>>()
        .join(", ")
}

/// Check if sequence is palindromic and real-rooted through all n.
/// Returns (all_palindromic, all_real_rooted, first_rr_fail, polys).
fn check_sequence(
    max_n: u8,
    poly_fn: &dyn Fn(u8) -> Vec<i64>,
) -> (bool, bool, Option<u8>, Vec<Vec<i64>>) {
    let mut polys = Vec::new();
    let mut all_palindromic = true;
    let mut all_rr = true;
    let mut first_rr_fail = None;

    for n in 1..=max_n {
        let coeffs = poly_fn(n);
        if !is_palindromic(&coeffs) {
            all_palindromic = false;
            return (false, false, None, polys);
        }
        if all_rr && !is_real_rooted(&coeffs) {
            all_rr = false;
            first_rr_fail = Some(n);
        }
        polys.push(coeffs);
    }
    (all_palindromic, all_rr, first_rr_fail, polys)
}

fn main() {
    let max_n: u8 = env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    let stats: &[(&str, Stat)] = &[
        ("des", Stat::Des),
        ("exc", Stat::Exc),
        ("peak", Stat::Peak),
    ];

    // ============================================================
    // Single length-3 patterns
    // ============================================================
    println!("=== Single length-3 patterns (n=1..{}) ===\n", max_n);

    let patterns3 = ["123", "132", "213", "231", "312", "321"];
    for &(stat_name, stat) in stats {
        for pat_str in &patterns3 {
            let pat = parse_sequence(pat_str);
            let (pal, rr, rr_fail, polys) = check_sequence(max_n, &|n| {
                let s_set = avoiding_permutations(n, &[pat.clone()]);
                backtrack_stat_poly(n, &s_set, stat)
            });
            if pal && polys.len() >= 4 && polys.iter().any(|p| p.len() > 2) {
                let rr_str = if rr {
                    format!("RR ✓")
                } else {
                    format!("RR fails n={}", rr_fail.unwrap())
                };
                println!("Av_{} + {}: palindromic, {}", pat_str, stat_name, rr_str);
                for (i, p) in polys.iter().enumerate() {
                    println!("  n={}: [{}]", i + 1, format_poly(p));
                }
                println!();
            }
        }
    }

    // ============================================================
    // Named families
    // ============================================================
    println!("\n=== Named families (n=1..{}) ===\n", max_n);

    let named_sets: Vec<(&str, Box<dyn Fn(&[u8]) -> bool>)> = vec![
        ("involutions", Box::new(|p: &[u8]| is_involution(p))),
        ("derangements", Box::new(|p: &[u8]| is_derangement(p))),
        ("vexillary", Box::new(|p: &[u8]| is_vexillary(p))),
        ("grassmann", Box::new(|p: &[u8]| is_grassmann(p))),
        ("separable", Box::new(|p: &[u8]| is_separable(p))),
        ("simsun", Box::new(|p: &[u8]| is_simsun(p))),
    ];

    for &(stat_name, stat) in stats {
        for (set_name, filter_fn) in &named_sets {
            let (pal, rr, rr_fail, polys) = check_sequence(max_n, &|n| {
                let all = all_permutations(n);
                let s_set: Vec<Vec<u8>> = all.into_iter().filter(|p| filter_fn(p)).collect();
                if s_set.is_empty() {
                    return vec![0];
                }
                backtrack_stat_poly(n, &s_set, stat)
            });
            if pal && polys.len() >= 4 && polys.iter().any(|p| p.len() > 2) {
                let rr_str = if rr {
                    format!("RR ✓")
                } else {
                    format!("RR fails n={}", rr_fail.unwrap())
                };
                println!("{} + {}: palindromic, {}", set_name, stat_name, rr_str);
                for (i, p) in polys.iter().enumerate() {
                    println!("  n={}: [{}]", i + 1, format_poly(p));
                }
                println!();
            }
        }
    }

    // ============================================================
    // Two-pattern avoidance (length 3)
    // ============================================================
    println!("\n=== Two-pattern avoidance, length 3 (n=1..{}) ===\n", max_n);

    for &(stat_name, stat) in stats {
        for i in 0..patterns3.len() {
            for j in i + 1..patterns3.len() {
                let pat1 = parse_sequence(patterns3[i]);
                let pat2 = parse_sequence(patterns3[j]);
                let (pal, rr, rr_fail, polys) = check_sequence(max_n, &|n| {
                    let s_set = avoiding_permutations(n, &[pat1.clone(), pat2.clone()]);
                    if s_set.is_empty() {
                        return vec![0];
                    }
                    backtrack_stat_poly(n, &s_set, stat)
                });
                if pal && polys.len() >= 4 && polys.iter().any(|p| p.len() > 2) {
                    let rr_str = if rr {
                        format!("RR ✓")
                    } else {
                        format!("RR fails n={}", rr_fail.unwrap())
                    };
                    println!(
                        "Av_{{{},{}}} + {}: palindromic, {}",
                        patterns3[i], patterns3[j], stat_name, rr_str
                    );
                    for (k, p) in polys.iter().enumerate() {
                        println!("  n={}: [{}]", k + 1, format_poly(p));
                    }
                    println!();
                }
            }
        }
    }

    // ============================================================
    // Length-4 single patterns (only check up to max_n-1 since n! grows fast)
    // ============================================================
    let max_n4 = max_n.min(7);
    println!(
        "\n=== Single length-4 patterns (n=1..{}) ===\n",
        max_n4
    );

    // Generate all 24 length-4 patterns
    let mut patterns4: Vec<String> = Vec::new();
    let p4 = all_permutations(4);
    for perm in &p4 {
        patterns4.push(
            perm.iter()
                .map(|&x| x.to_string())
                .collect::<Vec<_>>()
                .join(""),
        );
    }

    for &(stat_name, stat) in stats {
        for pat_str in &patterns4 {
            let pat = parse_sequence(pat_str);
            let (pal, rr, rr_fail, polys) = check_sequence(max_n4, &|n| {
                let s_set = avoiding_permutations(n, &[pat.clone()]);
                backtrack_stat_poly(n, &s_set, stat)
            });
            if pal && polys.len() >= 4 && polys.iter().any(|p| p.len() > 2) {
                let rr_str = if rr {
                    format!("RR ✓")
                } else {
                    format!("RR fails n={}", rr_fail.unwrap())
                };
                println!("Av_{} + {}: palindromic, {}", pat_str, stat_name, rr_str);
                for (k, p) in polys.iter().enumerate() {
                    println!("  n={}: [{}]", k + 1, format_poly(p));
                }
                println!();
            }
        }
    }

    // ============================================================
    // Catalan search orders with des, exc, peak
    // ============================================================
    println!(
        "\n=== Catalan search orders (n=1..{}) ===",
        max_n
    );
    println!("Only showing NON-Narayana palindromic results.\n");

    // Narayana for reference: C(n,k)^2/n = N(n,k)
    // Rows: [1], [1,1], [1,3,1], [1,6,6,1], [1,10,20,10,1], ...

    for &(stat_name, stat) in stats {
        for pi_pat_str in &patterns3 {
            let pi_p = parse_sequence(pi_pat_str);
            for s_pat_str in &patterns3 {
                let s_p = parse_sequence(s_pat_str);
                let (pal, rr, rr_fail, polys) = check_sequence(max_n, &|n| {
                    let all = all_permutations(n);
                    let pi_set: Vec<Vec<u8>> = all
                        .iter()
                        .filter(|p| !contains_pattern(p, &pi_p))
                        .cloned()
                        .collect();
                    let s_set: Vec<Vec<u8>> = all
                        .into_iter()
                        .filter(|p| !contains_pattern(p, &s_p))
                        .collect();
                    backtrack_stat_poly_subset(&pi_set, &s_set, stat)
                });
                if pal && polys.len() >= 4 && polys.iter().any(|p| p.len() > 2) {
                    // Check if it's Narayana or Eulerian (skip those)
                    let is_narayana = polys.len() >= 5
                        && polys[4] == vec![1, 10, 20, 10, 1];
                    let is_eulerian = polys.len() >= 4
                        && polys[3] == vec![1, 11, 11, 1];
                    if is_narayana || is_eulerian {
                        continue;
                    }

                    let rr_str = if rr {
                        format!("RR ✓")
                    } else {
                        format!("RR fails n={}", rr_fail.unwrap())
                    };
                    println!(
                        "Av_{} pi + Av_{} S + {}: palindromic, {}",
                        pi_pat_str, s_pat_str, stat_name, rr_str
                    );
                    for (k, p) in polys.iter().enumerate() {
                        println!("  n={}: [{}]", k + 1, format_poly(p));
                    }
                    println!();
                }
            }
        }
    }
}
