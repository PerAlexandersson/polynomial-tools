/// Compare h*-vector of the order polytope of the k-alternating poset
/// with the swaps polynomial H^(k)_n(t).
///
/// If they match, real-rootedness of H^(k)_n follows from Stanley's
/// theorem + general results on order polytopes.
use combinatoric_core::poset::Poset;
use combpoly::permutation::{alternating_permutations, k_alternating_permutations};
use combpoly::statistics::{compute, Stat};
use polynomial_tools::real_rootedness::format_poly;

fn build_swaps_poly(perms: &[Vec<u8>]) -> Vec<i64> {
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
    println!("=== Comparing h*-vector of O(P^(k)_n) with H^(k)_n(t) ===\n");

    for k in 2..=4u8 {
        println!("--- k={} ---", k);
        for n in 2..=10usize {
            // Compute H^(k)_n via swaps
            let perms = k_alternating_permutations(n as u8, k);
            let swaps_poly = build_swaps_poly(&perms);

            // Compute h*-vector of order polytope of k-alternating poset
            let poset = Poset::k_alternating(n, k as usize);

            // Verify linear extensions match
            let lin_ext = poset.num_linear_extensions();
            assert_eq!(
                lin_ext,
                perms.len(),
                "k={} n={}: linear extensions ({}) != k-alt perms ({})",
                k,
                n,
                lin_ext,
                perms.len()
            );

            // Compute P-Eulerian (descent poly on linear extensions)
            let p_euler = poset.p_eulerian_polynomial();

            // Compute h* of order polytope
            let hstar = poset.order_polytope_hstar();

            let match_swaps = hstar == swaps_poly;
            let match_peuler = hstar == p_euler;

            println!(
                "  n={:>2}: |A|={:<6}  h*={:<45}  H(t)={:<45}  W(t)={:<45}  h*==H? {}  h*==W? {}",
                n,
                perms.len(),
                format_poly(&hstar),
                format_poly(&swaps_poly),
                format_poly(&p_euler),
                if match_swaps { "✓" } else { "✗" },
                if match_peuler { "✓" } else { "✗" },
            );
        }
        println!();
    }
}
