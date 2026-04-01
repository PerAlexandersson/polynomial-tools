/// Test multivariate stability of the long-swaps generating polynomial.
///
/// For S ⊆ [n-1], define:
///   M_{n,S}(x_1,...,x_{n-1}) = Σ_{σ ∈ D(n,S)} Π_{i ∈ SWAPS(σ)} x_i
///
/// A polynomial is real stable if it has no root where all variables
/// have strictly positive imaginary part.
///
/// Equivalent characterization (Borcea-Brändén): f(x₁,...,xₘ) is real stable
/// iff for all a ∈ ℝᵐ and b ∈ ℝ₊ᵐ, the univariate polynomial
/// g(t) = f(a + tb) is either identically zero or has only real roots.
///
/// We test this by:
/// 1. Sampling many (a, b) with b > 0 and checking if g(t) is real-rooted.
/// 2. Testing specific structured substitutions.
/// 3. Checking the "diagonal" substitution x_i = t (which must be real-rooted).
/// 4. Checking mixed substitutions like x_i = s, x_j = t for pairs.
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute_set, descent_set_bitmask, SetStat};
use polynomial_tools::real_rootedness::is_real_rooted;
use std::collections::BTreeMap;

/// Compute the multivariate swaps monomial for a permutation.
/// Returns the exponent vector: e[i] = 1 if i ∈ SWAPS(σ), else 0.
fn swaps_exponent(sigma: &[u8]) -> Vec<u8> {
    let n = sigma.len();
    let swaps = compute_set(sigma, SetStat::SwapsSet);
    let mut exp = vec![0u8; n.saturating_sub(1)];
    for i in swaps {
        // SwapsSet returns 1-indexed values i where i is a long swap
        // index into x_{i} (0-indexed: i-1)
        if i >= 1 && i <= n - 1 {
            exp[i - 1] = 1;
        }
    }
    exp
}

/// Evaluate the multivariate polynomial at a point.
/// M(x_1,...,x_{n-1}) = Σ_{σ ∈ perms} Π_{i ∈ SWAPS(σ)} x_i
fn eval_multivar(perms: &[&Vec<u8>], point: &[f64]) -> f64 {
    let mut sum = 0.0f64;
    for sigma in perms {
        let exp = swaps_exponent(sigma);
        let mut prod = 1.0f64;
        for (i, &e) in exp.iter().enumerate() {
            if e > 0 {
                prod *= point[i];
            }
        }
        sum += prod;
    }
    sum
}

/// Build the univariate polynomial g(t) = M(a + t*b) where a, b ∈ ℝ^m.
/// Since each x_i appears with degree 0 or 1, M is multilinear.
/// g(t) = M(a₁+tb₁, ..., aₘ+tbₘ)
/// For a multilinear polynomial, this is a polynomial in t of degree ≤ m.
fn build_univariate_slice(perms: &[&Vec<u8>], a: &[f64], b: &[f64]) -> Vec<f64> {
    let m = a.len(); // number of variables
                     // degree of g(t) is at most m, but practically at most max_swaps
    let max_deg = m;
    let mut coeffs = vec![0.0f64; max_deg + 1];

    for sigma in perms {
        let exp = swaps_exponent(sigma);
        // This monomial contributes Π_{i: exp[i]=1} (a_i + t b_i)
        // = Σ_{S ⊆ support} (Π_{i∈S} b_i) (Π_{i∉S, in support} a_i) t^|S|
        let support: Vec<usize> = (0..m).filter(|&i| exp[i] > 0).collect();
        let k = support.len();

        // Enumerate subsets of support
        for mask in 0..(1u64 << k) {
            let mut prod = 1.0f64;
            let mut deg = 0usize;
            for (j, &idx) in support.iter().enumerate() {
                if mask & (1 << j) != 0 {
                    prod *= b[idx];
                    deg += 1;
                } else {
                    prod *= a[idx];
                }
            }
            coeffs[deg] += prod;
        }
    }

    // Trim trailing near-zero coefficients
    while coeffs.len() > 1 && coeffs.last().unwrap().abs() < 1e-12 {
        coeffs.pop();
    }

    coeffs
}

/// Check if a floating-point polynomial has all real roots (approximately).
/// Uses the discriminant / companion matrix approach.
fn is_approx_real_rooted(coeffs: &[f64]) -> bool {
    if coeffs.len() <= 2 {
        return true;
    }

    // Convert to integer approximation and use exact check
    let scale = coeffs.iter().map(|c| c.abs()).fold(0.0f64, f64::max);
    if scale < 1e-12 {
        return true; // zero polynomial
    }

    // Scale up and round to integers
    let precision = 1e8;
    let int_coeffs: Vec<i64> = coeffs
        .iter()
        .map(|c| (c / scale * precision).round() as i64)
        .collect();

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

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    println!("═══ Multivariate stability test for M_{{n,S}}(x) ═══\n");

    // Random directions for testing
    let test_directions: Vec<(Vec<f64>, Vec<f64>)> = {
        // Use deterministic "random" directions
        let mut dirs = Vec::new();
        // Axis-aligned: x_i = a_i + t, others fixed
        // Will be generated per-n below
        dirs
    };

    let _ = test_directions; // suppress warning

    for n in 3..=max_n {
        let all = all_permutations(n);
        let m = (n - 1) as usize; // number of variables

        let mut by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for s in &all {
            by_des.entry(descent_set_bitmask(s)).or_default().push(s);
        }

        let mut total_sets = 0u32;
        let mut stable_count = 0u32;
        let mut unstable_count = 0u32;
        let mut unstable_examples: Vec<String> = Vec::new();

        for (&ds, perms) in &by_des {
            total_sets += 1;

            let mut all_slices_rr = true;

            // Test 1: All single-variable substitutions
            // Fix all x_j = c_j for j ≠ i, vary x_i = t
            // Use c_j ∈ {1, 2, 0.5}
            let test_vals = [1.0, 2.0, 0.5, 3.0];

            'outer: for &base_val in &test_vals {
                let base: Vec<f64> = vec![base_val; m];

                for i in 0..m {
                    let mut a = base.clone();
                    let mut b = vec![0.0f64; m];
                    a[i] = 0.0;
                    b[i] = 1.0;

                    let g = build_univariate_slice(&perms, &a, &b);
                    if !is_approx_real_rooted(&g) {
                        all_slices_rr = false;
                        if unstable_examples.len() < 3 {
                            let ds_str = format!("{:b}", ds);
                            unstable_examples.push(format!(
                                "n={} Des={} var x_{} base={} g={:?}",
                                n,
                                ds_str,
                                i + 1,
                                base_val,
                                g
                            ));
                        }
                        break 'outer;
                    }
                }
            }

            if !all_slices_rr {
                unstable_count += 1;
                continue;
            }

            // Test 2: Two-variable substitutions
            // Fix all but x_i, x_j, then set x_i = a_i + t b_i, x_j = a_j + t b_j
            let mut two_var_ok = true;
            let two_var_tests: Vec<(f64, f64, f64, f64)> = vec![
                (0.0, 1.0, 0.0, 1.0), // both from 0, slope 1
                (1.0, 1.0, 0.0, 1.0), // one shifted
                (0.0, 1.0, 0.0, 2.0), // different slopes
                (1.0, 2.0, 2.0, 1.0), // general
                (0.5, 1.0, 1.5, 1.0), // other
            ];

            'outer2: for i in 0..m {
                for j in i + 1..m {
                    for &(ai, bi, aj, bj) in &two_var_tests {
                        let mut a = vec![1.0f64; m]; // other vars = 1
                        let mut b = vec![0.0f64; m];
                        a[i] = ai;
                        b[i] = bi;
                        a[j] = aj;
                        b[j] = bj;

                        let g = build_univariate_slice(&perms, &a, &b);
                        if !is_approx_real_rooted(&g) {
                            two_var_ok = false;
                            if unstable_examples.len() < 3 {
                                unstable_examples.push(format!(
                                    "n={} Des={:b} 2-var x_{},x_{} ({},{},{},{}) FAIL",
                                    n,
                                    ds,
                                    i + 1,
                                    j + 1,
                                    ai,
                                    bi,
                                    aj,
                                    bj
                                ));
                            }
                            break 'outer2;
                        }
                    }
                }
            }

            if !two_var_ok {
                unstable_count += 1;
                continue;
            }

            // Test 3: Random-ish general directions (b > 0)
            let general_tests: Vec<Vec<(f64, f64)>> = vec![
                (0..m).map(|i| (i as f64 + 1.0, 1.0)).collect(),
                (0..m).map(|i| (0.0, (i + 1) as f64)).collect(),
                (0..m).map(|i| ((m - i) as f64, (i + 1) as f64)).collect(),
                (0..m).map(|_| (0.5, 1.5)).collect(),
                (0..m)
                    .map(|i| (if i % 2 == 0 { 0.0 } else { 1.0 }, 1.0))
                    .collect(),
            ];

            let mut general_ok = true;
            for test in &general_tests {
                let a: Vec<f64> = test.iter().map(|(ai, _)| *ai).collect();
                let b: Vec<f64> = test.iter().map(|(_, bi)| *bi).collect();

                let g = build_univariate_slice(&perms, &a, &b);
                if !is_approx_real_rooted(&g) {
                    general_ok = false;
                    break;
                }
            }

            if !general_ok {
                unstable_count += 1;
            } else {
                stable_count += 1;
            }
        }

        println!(
            "n={}: {}/{} pass all stability tests, {} fail",
            n, stable_count, total_sets, unstable_count,
        );
        for ex in &unstable_examples {
            println!("  {}", ex);
        }
    }

    // Part 2: For the full polynomial on S_n (all descent sets combined),
    // check stability
    println!("\n═══ Full S_n multivariate stability ═══\n");

    for n in 3..=max_n {
        let all = all_permutations(n);
        let refs: Vec<&Vec<u8>> = all.iter().collect();
        let m = (n - 1) as usize;

        let mut all_ok = true;
        let base_vals = [1.0, 2.0, 0.5];

        for &base_val in &base_vals {
            for i in 0..m {
                let mut a = vec![base_val; m];
                let mut b = vec![0.0f64; m];
                a[i] = 0.0;
                b[i] = 1.0;

                let g = build_univariate_slice(&refs, &a, &b);
                if !is_approx_real_rooted(&g) {
                    all_ok = false;
                    println!("  n={}: FAIL at var x_{} base={}", n, i + 1, base_val,);
                }
            }
        }

        // General direction tests
        let general_tests: Vec<Vec<(f64, f64)>> = vec![
            (0..m).map(|i| (i as f64 + 1.0, 1.0)).collect(),
            (0..m).map(|i| (0.0, (i + 1) as f64)).collect(),
            (0..m).map(|i| ((m - i) as f64, (i + 1) as f64)).collect(),
        ];
        for test in &general_tests {
            let a: Vec<f64> = test.iter().map(|(ai, _)| *ai).collect();
            let b: Vec<f64> = test.iter().map(|(_, bi)| *bi).collect();
            let g = build_univariate_slice(&refs, &a, &b);
            if !is_approx_real_rooted(&g) {
                all_ok = false;
                println!("  n={}: FAIL at general direction", n);
            }
        }

        if all_ok {
            println!("n={}: full S_n polynomial passes stability tests ✓", n);
        }
    }
}
