/// Test total non-negativity of the coefficient matrix for position-refined
/// descent-swap polynomials.
///
/// For each descent set S with 1 ∉ S in S_n (n = 4..max_n), form the matrix M where:
///   Row i  = coefficients of L_{n,S}^{(p_i)}(t)  (position-refined poly)
///   Col j  = coefficient of t^j
/// with rows ordered by position p_i ∈ P(S).
///
/// Check: is every 2×2 minor of M non-negative (totally non-negative)?
/// Also check 3×3 minors when feasible.
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::format_poly;
use std::collections::HashMap;

/// Build polynomial (coefficient vector) from an iterator of statistic values.
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

/// Compute valid positions P(S) for descent set S (as bitmask) in S_n.
/// These are:
///   - p ∈ S such that p-1 ∉ S  (left boundaries of descent runs), and
///   - p = n if n-1 ∉ S.
fn valid_positions(s_mask: u64, n: u8) -> Vec<u8> {
    let mut positions = Vec::new();
    for p in 1..n {
        let p_in_s = (s_mask >> (p - 1)) & 1 == 1;
        let pm1_in_s = if p >= 2 {
            (s_mask >> (p - 2)) & 1 == 1
        } else {
            false
        };
        if p_in_s && !pm1_in_s {
            positions.push(p);
        }
    }
    // p = n if n-1 ∉ S
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

/// Build the coefficient matrix M for a given descent class.
/// Row i corresponds to position p_i, column j to coefficient of t^j.
/// Returns (positions, matrix).
fn build_coeff_matrix(class: &[&Vec<u8>], positions: &[u8], n: u8) -> (Vec<u8>, Vec<Vec<i64>>) {
    let mut pos_polys: Vec<(u8, Vec<i64>)> = Vec::new();
    for &p in positions {
        let vals: Vec<usize> = class
            .iter()
            .filter(|pi| {
                let pos_of_n = pi.iter().position(|&v| v == n).unwrap();
                (pos_of_n + 1) as u8 == p
            })
            .map(|pi| compute(pi, Stat::Swaps))
            .collect();
        if vals.is_empty() {
            continue;
        }
        let poly = build_poly_from_vals(vals.into_iter());
        pos_polys.push((p, poly));
    }

    if pos_polys.is_empty() {
        return (vec![], vec![]);
    }

    // Find max degree across all position polys
    let max_deg = pos_polys.iter().map(|(_, p)| p.len()).max().unwrap_or(1);

    // Build matrix: pad shorter polys with zeros
    let positions_used: Vec<u8> = pos_polys.iter().map(|(p, _)| *p).collect();
    let matrix: Vec<Vec<i64>> = pos_polys
        .iter()
        .map(|(_, poly)| {
            let mut row = poly.clone();
            row.resize(max_deg, 0);
            row
        })
        .collect();

    (positions_used, matrix)
}

/// Check all 2×2 minors of a matrix for non-negativity.
/// Returns (all_ok, list of failing minors as (i1, i2, j1, j2, value)).
fn check_2x2_minors(matrix: &[Vec<i64>]) -> (bool, Vec<(usize, usize, usize, usize, i64)>) {
    let nrows = matrix.len();
    if nrows < 2 {
        return (true, vec![]);
    }
    let ncols = matrix[0].len();
    if ncols < 2 {
        return (true, vec![]);
    }
    let mut failures = Vec::new();
    for i1 in 0..nrows {
        for i2 in (i1 + 1)..nrows {
            for j1 in 0..ncols {
                for j2 in (j1 + 1)..ncols {
                    let det = matrix[i1][j1] * matrix[i2][j2] - matrix[i1][j2] * matrix[i2][j1];
                    if det < 0 {
                        failures.push((i1, i2, j1, j2, det));
                    }
                }
            }
        }
    }
    (failures.is_empty(), failures)
}

/// Check all 3×3 minors of a matrix for non-negativity.
/// Returns (all_ok, number of failures).
fn check_3x3_minors(matrix: &[Vec<i64>]) -> (bool, usize) {
    let nrows = matrix.len();
    if nrows < 3 {
        return (true, 0);
    }
    let ncols = matrix[0].len();
    if ncols < 3 {
        return (true, 0);
    }
    let mut fail_count = 0usize;
    for i1 in 0..nrows {
        for i2 in (i1 + 1)..nrows {
            for i3 in (i2 + 1)..nrows {
                for j1 in 0..ncols {
                    for j2 in (j1 + 1)..ncols {
                        for j3 in (j2 + 1)..ncols {
                            let det = matrix[i1][j1]
                                * (matrix[i2][j2] * matrix[i3][j3]
                                    - matrix[i2][j3] * matrix[i3][j2])
                                - matrix[i1][j2]
                                    * (matrix[i2][j1] * matrix[i3][j3]
                                        - matrix[i2][j3] * matrix[i3][j1])
                                + matrix[i1][j3]
                                    * (matrix[i2][j1] * matrix[i3][j2]
                                        - matrix[i2][j2] * matrix[i3][j1]);
                            if det < 0 {
                                fail_count += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    (fail_count == 0, fail_count)
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    println!("Testing total non-negativity of position-refined coefficient matrices");
    println!("for L_{{n,S}}^{{(p)}}(t) with 1 ∉ S, n = 4..{}\n", max_n);

    let mut total_tested = 0u64;
    let mut total_tnn_2x2 = 0u64;
    let mut total_tnn_3x3 = 0u64;
    let mut total_checked_3x3 = 0u64;

    for n in 4..=max_n {
        let perms = all_permutations(n);
        println!("========== n = {} ({} perms) ==========", n, perms.len());

        // Group permutations by descent set
        let mut by_descent: HashMap<u64, Vec<&Vec<u8>>> = HashMap::new();
        for pi in &perms {
            let mask = descent_set_bitmask(pi);
            by_descent.entry(mask).or_default().push(pi);
        }

        let mut sorted_masks: Vec<u64> = by_descent.keys().copied().collect();
        sorted_masks.sort();

        let mut n_tested = 0u64;
        let mut n_tnn_2x2 = 0u64;
        let mut n_tnn_3x3 = 0u64;
        let mut n_checked_3x3 = 0u64;

        for &mask in &sorted_masks {
            // Only consider S with 1 ∉ S (bit 0 not set)
            if mask & 1 != 0 {
                continue;
            }

            let class = &by_descent[&mask];
            let s_str = descent_set_to_string(mask, n);
            let positions = valid_positions(mask, n);

            let (pos_used, matrix) = build_coeff_matrix(class, &positions, n);
            if matrix.len() < 2 {
                // Need at least 2 rows for a non-trivial minor check
                continue;
            }

            n_tested += 1;

            // Print matrix
            println!(
                "  S={} positions={:?} matrix {}x{}:",
                s_str,
                pos_used,
                matrix.len(),
                matrix[0].len()
            );
            for (i, row) in matrix.iter().enumerate() {
                let poly_str = format_poly(row);
                println!("    p={}: {} => {:?}", pos_used[i], poly_str, row);
            }

            // Check 2x2 minors
            let (ok_2x2, failures_2x2) = check_2x2_minors(&matrix);
            if ok_2x2 {
                n_tnn_2x2 += 1;
                println!("    2x2 minors: ALL NON-NEGATIVE");
            } else {
                println!("    2x2 minors: {} FAILURES", failures_2x2.len());
                for &(i1, i2, j1, j2, det) in failures_2x2.iter().take(5) {
                    println!(
                        "      rows({},{}) cols({},{}) => det = {} (M[{}][{}]*M[{}][{}] - M[{}][{}]*M[{}][{}] = {}*{} - {}*{} = {})",
                        pos_used[i1], pos_used[i2], j1, j2,
                        det,
                        i1, j1, i2, j2, i1, j2, i2, j1,
                        matrix[i1][j1], matrix[i2][j2],
                        matrix[i1][j2], matrix[i2][j1],
                        det
                    );
                }
            }

            // Check 3x3 minors if matrix is big enough
            if matrix.len() >= 3 && matrix[0].len() >= 3 {
                let (ok_3x3, fail_3x3) = check_3x3_minors(&matrix);
                n_checked_3x3 += 1;
                if ok_3x3 {
                    n_tnn_3x3 += 1;
                    println!("    3x3 minors: ALL NON-NEGATIVE");
                } else {
                    println!("    3x3 minors: {} FAILURES", fail_3x3);
                }
            }

            println!();
        }

        total_tested += n_tested;
        total_tnn_2x2 += n_tnn_2x2;
        total_tnn_3x3 += n_tnn_3x3;
        total_checked_3x3 += n_checked_3x3;

        println!(
            "n={}: {}/{} TNN(2x2), {}/{} TNN(3x3)\n",
            n, n_tnn_2x2, n_tested, n_tnn_3x3, n_checked_3x3,
        );
    }

    println!("=== SUMMARY ===");
    println!(
        "Total descent sets tested (1 ∉ S, ≥2 positions): {}",
        total_tested
    );
    println!(
        "Totally non-negative (2x2): {}/{}",
        total_tnn_2x2, total_tested
    );
    println!(
        "Totally non-negative (3x3): {}/{}",
        total_tnn_3x3, total_checked_3x3
    );
    if total_tnn_2x2 == total_tested {
        println!("\nCONJECTURE HOLDS: All matrices are TNN (2x2 minors)!");
    } else {
        println!(
            "\nCONJECTURE FAILS: {}/{} matrices have negative 2x2 minors.",
            total_tested - total_tnn_2x2,
            total_tested,
        );
    }
    if total_tnn_3x3 == total_checked_3x3 && total_checked_3x3 > 0 {
        println!("CONJECTURE HOLDS: All matrices are TNN (3x3 minors)!");
    } else if total_checked_3x3 > 0 {
        println!(
            "CONJECTURE FAILS: {}/{} matrices have negative 3x3 minors.",
            total_checked_3x3 - total_tnn_3x3,
            total_checked_3x3,
        );
    }
}
