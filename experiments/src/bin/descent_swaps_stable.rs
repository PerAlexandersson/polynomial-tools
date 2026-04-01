/// Explore multivariate stability of the swaps polynomial:
///   Φ_{n,S}(x_1,...,x_{n-1}) = Σ_{σ ∈ D(n,S)} Π_{i ∈ Swaps(σ)} x_i
///
/// If Φ is real-stable (nonvanishing when all Im(x_i) > 0),
/// then L_{n,S}(t) = Φ(t,...,t) is real-rooted, and the
/// position-of-max pieces form a compatible family.
///
/// For multi-affine polynomials, stability can be checked by verifying
/// that Φ(a_1 + b_1·z, ..., a_{n-1} + b_{n-1}·z) ∈ R[z] is real-rooted
/// for all a_i ∈ R, b_i ≥ 0.
///
/// We also check "half-plane" property and Rayleigh-type conditions.
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute_set, descent_set_bitmask, SetStat};
use polynomial_tools::real_rootedness::{format_poly, is_real_rooted};
use std::collections::HashMap;

/// Multi-affine polynomial represented as map from monomial support (bitmask) to coefficient.
/// Variable x_i corresponds to bit i-1 in the bitmask (i is 1-indexed).
type MultiAffine = HashMap<u64, i64>;

/// Compute the multi-affine swaps polynomial for a descent class.
fn compute_multiaffine(perms: &[&Vec<u8>]) -> MultiAffine {
    let mut poly = MultiAffine::new();
    for pi in perms {
        let swaps_set = compute_set(pi, SetStat::SwapsSet);
        let mut mask = 0u64;
        for &i in &swaps_set {
            mask |= 1 << (i - 1); // i is 1-indexed
        }
        *poly.entry(mask).or_insert(0) += 1;
    }
    poly
}

/// Evaluate the multi-affine polynomial Φ at x_i = a_i + b_i * z,
/// returning the resulting univariate polynomial in z as a coefficient vector.
fn substitute_linear(phi: &MultiAffine, a: &[f64], b: &[f64], num_vars: usize) -> Vec<f64> {
    // Degree is at most num_vars
    let mut result = vec![0.0f64; num_vars + 1];

    for (&monomial_mask, &coeff) in phi {
        // For this monomial Π_{i ∈ S} x_i where S is the support,
        // substitute x_i = a_i + b_i * z.
        // The product Π(a_i + b_i z) is a polynomial in z.
        // We compute it by iterating over subsets of S.
        let support: Vec<usize> = (0..num_vars)
            .filter(|&i| (monomial_mask >> i) & 1 == 1)
            .collect();
        let k = support.len();

        // Expand Π_{i ∈ S} (a_i + b_i z) using binomial-like expansion
        // Coefficient of z^j = Σ over j-subsets T of S: (Π_{i∈T} b_i)(Π_{i∈S\T} a_i)
        let mut mono_poly = vec![0.0f64; k + 1];
        for subset_mask in 0u64..(1u64 << k) {
            let j = subset_mask.count_ones() as usize;
            let mut term = coeff as f64;
            for (idx, &var) in support.iter().enumerate() {
                if (subset_mask >> idx) & 1 == 1 {
                    term *= b[var];
                } else {
                    term *= a[var];
                }
            }
            mono_poly[j] += term;
        }
        for (j, &c) in mono_poly.iter().enumerate() {
            result[j] += c;
        }
    }

    // Trim trailing zeros
    while result.len() > 1 && result.last().unwrap().abs() < 1e-12 {
        result.pop();
    }
    result
}

/// Check if a floating-point polynomial is approximately real-rooted.
/// Uses the companion matrix / eigenvalue approach.
fn is_approx_real_rooted(coeffs: &[f64]) -> bool {
    if coeffs.len() <= 2 {
        return true;
    }
    // Convert to integer approximation and use exact checker
    // Scale up to avoid precision issues
    let max_coeff = coeffs.iter().map(|c| c.abs()).fold(0.0f64, f64::max);
    if max_coeff < 1e-10 {
        return true;
    }
    let scale = 1e6 / max_coeff;
    let int_coeffs: Vec<i64> = coeffs.iter().map(|&c| (c * scale).round() as i64).collect();
    // Trim
    let mut trimmed = int_coeffs;
    while trimmed.len() > 1 && *trimmed.last().unwrap() == 0 {
        trimmed.pop();
    }
    if trimmed.len() <= 2 {
        return true;
    }
    is_real_rooted(&trimmed)
}

/// Check Rayleigh condition: for all i ≠ j,
///   Φ · (∂²Φ/∂x_i∂x_j) ≤ (∂Φ/∂x_i)(∂Φ/∂x_j)
/// evaluated at x = (1,1,...,1).
/// This is a NECESSARY condition for stability of multi-affine polynomials.
fn check_rayleigh(phi: &MultiAffine, num_vars: usize) -> bool {
    // Compute Φ(1,...,1) = Σ coeffs
    let phi_val: i64 = phi.values().sum();
    if phi_val == 0 {
        return true;
    }

    for i in 0..num_vars {
        for j in (i + 1)..num_vars {
            // ∂Φ/∂x_i at (1,...,1) = Σ_{S: i ∈ S} coeff(S)
            let dphi_i: i64 = phi
                .iter()
                .filter(|(&m, _)| (m >> i) & 1 == 1)
                .map(|(_, &c)| c)
                .sum();

            let dphi_j: i64 = phi
                .iter()
                .filter(|(&m, _)| (m >> j) & 1 == 1)
                .map(|(_, &c)| c)
                .sum();

            // ∂²Φ/∂x_i∂x_j at (1,...,1) = Σ_{S: {i,j} ⊆ S} coeff(S)
            let d2phi_ij: i64 = phi
                .iter()
                .filter(|(&m, _)| (m >> i) & 1 == 1 && (m >> j) & 1 == 1)
                .map(|(_, &c)| c)
                .sum();

            // Rayleigh: Φ · d²Φ_ij ≤ dΦ_i · dΦ_j
            if phi_val * d2phi_ij > dphi_i * dphi_j {
                return false;
            }
        }
    }
    true
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

    let num_random_tests = 500; // number of random substitutions per polynomial

    // Deterministic "difficult" test directions
    let special_dirs: Vec<Vec<f64>> = vec![
        // All equal (gives L_{n,S}(t))
        // Already known to be rr, skip
    ];

    use std::collections::BTreeMap;

    for n in 3..=max_n {
        let perms = all_permutations(n);
        let num_vars = (n - 1) as usize;

        println!("========== n = {} ==========", n);

        // Group by descent set
        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            let mask = descent_set_bitmask(pi);
            by_descent.entry(mask).or_default().push(pi);
        }

        let mut n_rayleigh_ok = 0;
        let mut n_rayleigh_total = 0;
        let mut n_stable_ok = 0;
        let mut n_stable_total = 0;

        for (&mask, class) in &by_descent {
            let s_str = descent_set_to_string(mask, n);
            let phi = compute_multiaffine(class);

            if phi.is_empty() {
                continue;
            }

            // Check Rayleigh condition (necessary for stability)
            n_rayleigh_total += 1;
            let rayleigh = check_rayleigh(&phi, num_vars);
            if rayleigh {
                n_rayleigh_ok += 1;
            } else {
                println!("  RAYLEIGH FAIL: S={}", s_str);
            }

            // Check stability via random linear substitutions
            n_stable_total += 1;
            let mut stable = true;

            // Use a simple PRNG for reproducibility
            let mut rng_state: u64 = mask.wrapping_mul(12345).wrapping_add(n as u64 * 67890);
            let next_rand = |state: &mut u64| -> f64 {
                *state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                // Map to [-2, 2] for a values, [0, 2] for b values
                ((*state >> 33) as f64) / (1u64 << 31) as f64
            };

            for _ in 0..num_random_tests {
                let a: Vec<f64> = (0..num_vars)
                    .map(|_| next_rand(&mut rng_state) * 4.0 - 2.0)
                    .collect();
                let b: Vec<f64> = (0..num_vars)
                    .map(|_| next_rand(&mut rng_state).abs() * 2.0)
                    .collect();

                let univar = substitute_linear(&phi, &a, &b, num_vars);
                if univar.len() > 2 && !is_approx_real_rooted(&univar) {
                    stable = false;
                    println!(
                        "  STABILITY FAIL: S={} a={:?} b={:?}",
                        s_str,
                        &a[..a.len().min(5)],
                        &b[..b.len().min(5)]
                    );
                    // Print the failing polynomial
                    let int_coeffs: Vec<i64> = univar.iter().map(|&c| c.round() as i64).collect();
                    println!("    poly ≈ {}", format_poly(&int_coeffs));
                    break;
                }
            }

            if stable {
                n_stable_ok += 1;
            }
        }

        println!(
            "n={}: Rayleigh {}/{}, Stability {}/{} (random, {} tests each)\n",
            n, n_rayleigh_ok, n_rayleigh_total, n_stable_ok, n_stable_total, num_random_tests
        );
    }

    println!("=== SUMMARY ===");
    println!("If all stability checks pass, the multi-affine polynomial");
    println!("Phi_{{n,S}}(x) = sum prod x_i is a strong candidate for a stable lift.");
    println!("A proof of stability would imply real-rootedness of L_{{n,S}}(t)");
    println!("and the compatible family property simultaneously.");
}
