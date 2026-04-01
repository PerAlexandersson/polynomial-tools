/// Long swaps on permutations with a fixed descent set.
///
/// For S ⊆ [n-1], let D(n, S) = {σ ∈ S_n : Des(σ) = S}.
/// Compute the long-swaps polynomial on D(n, S) and check real-rootedness.
///
/// This generalizes alternating (S = {2,4,6,...}) and k-alternating
/// (S = {k,2k,3k,...}).
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{format_poly, is_real_rooted};

fn descent_set_str(s: u64, n: u8) -> String {
    let positions: Vec<String> = (0..n - 1)
        .filter(|&i| s & (1 << i) != 0)
        .map(|i| (i + 1).to_string())
        .collect();
    if positions.is_empty() {
        "∅".to_string()
    } else {
        format!("{{{}}}", positions.join(","))
    }
}

fn build_poly(perms: &[&Vec<u8>]) -> Vec<i64> {
    if perms.is_empty() {
        return vec![0];
    }
    let max_s = perms
        .iter()
        .map(|s| compute(s, Stat::Swaps))
        .max()
        .unwrap_or(0);
    let mut coeffs = vec![0i64; max_s + 1];
    for s in perms {
        coeffs[compute(s, Stat::Swaps)] += 1;
    }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
        coeffs.pop();
    }
    coeffs
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    for n in 3..=max_n {
        let all = all_permutations(n);

        // Group by descent set
        let mut by_des: std::collections::BTreeMap<u64, Vec<&Vec<u8>>> =
            std::collections::BTreeMap::new();
        for s in &all {
            by_des.entry(descent_set_bitmask(s)).or_default().push(s);
        }

        let total_sets = by_des.len();
        let mut rr_count = 0;
        let mut fail_count = 0;
        let mut failures = Vec::new();

        for (&ds, perms) in &by_des {
            let poly = build_poly(perms);
            let rr = if poly.len() <= 2 {
                true
            } else {
                is_real_rooted(&poly)
            };
            if rr {
                rr_count += 1;
            } else {
                fail_count += 1;
                failures.push((ds, perms.len(), poly));
            }
        }

        println!(
            "n={}: {} descent sets, {} real-rooted, {} failures",
            n, total_sets, rr_count, fail_count
        );

        if fail_count > 0 && fail_count <= 10 {
            for (ds, count, poly) in &failures {
                println!(
                    "  FAIL: Des={} |D|={} {}",
                    descent_set_str(*ds, n),
                    count,
                    format_poly(poly)
                );
            }
        } else if fail_count > 10 {
            // Just show first few
            for (ds, count, poly) in failures.iter().take(5) {
                println!(
                    "  FAIL: Des={} |D|={} {}",
                    descent_set_str(*ds, n),
                    count,
                    format_poly(poly)
                );
            }
            println!("  ... and {} more failures", fail_count - 5);
        }
    }
}
