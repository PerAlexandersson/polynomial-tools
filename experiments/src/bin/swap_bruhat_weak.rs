/// Long swaps on Bruhat and weak order ideals.
///
/// For σ ∈ S_n, let T = {π : π ≤ σ} in Bruhat (or weak) order.
/// Compute the long-swaps polynomial on T and check real-rootedness.
use combpoly::order::{bruhat_lower_ideal, weak_lower_ideal};
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, Stat};
use polynomial_tools::real_rootedness::{format_poly, is_real_rooted};

fn build_poly(perms: &[Vec<u8>]) -> Vec<i64> {
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

fn perm_str(p: &[u8]) -> String {
    p.iter().map(|v| v.to_string()).collect()
}

fn main() {
    for n in 3..=6u8 {
        println!("=== n={} ===\n", n);
        let all = all_permutations(n);

        let mut weak_rr = 0;
        let mut weak_fail = 0;
        let mut bruhat_rr = 0;
        let mut bruhat_fail = 0;

        println!(
            "{:<8} {:<6} {:<6} {:<50} {:<5} {:<6} {:<50} {}",
            "σ", "|W|", "|B|", "weak poly", "w-rr", "|B|", "bruhat poly", "b-rr"
        );

        for sigma in &all {
            let weak = weak_lower_ideal(sigma);
            let bruhat = bruhat_lower_ideal(sigma);

            let wp = build_poly(&weak);
            let bp = build_poly(&bruhat);

            let wrr = if wp.len() <= 2 {
                true
            } else {
                is_real_rooted(&wp)
            };
            let brr = if bp.len() <= 2 {
                true
            } else {
                is_real_rooted(&bp)
            };

            if wrr {
                weak_rr += 1;
            } else {
                weak_fail += 1;
            }
            if brr {
                bruhat_rr += 1;
            } else {
                bruhat_fail += 1;
            }

            // Only print failures or interesting cases
            if !wrr || !brr {
                println!(
                    "{:<8} {:<6} {:<6} {:<50} {:<5} {:<6} {:<50} {}",
                    perm_str(sigma),
                    weak.len(),
                    bruhat.len(),
                    format_poly(&wp),
                    if wrr { "✓" } else { "✗" },
                    bruhat.len(),
                    format_poly(&bp),
                    if brr { "✓" } else { "✗" },
                );
            }
        }

        println!(
            "\nWeak:   {} real-rooted, {} failures out of {}",
            weak_rr,
            weak_fail,
            all.len()
        );
        println!(
            "Bruhat: {} real-rooted, {} failures out of {}\n",
            bruhat_rr,
            bruhat_fail,
            all.len()
        );
    }
}
