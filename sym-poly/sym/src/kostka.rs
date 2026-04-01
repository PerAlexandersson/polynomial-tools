use std::collections::BTreeMap;
use sym_poly_core::Partition;

/// Compute the Kostka coefficient K(λ, μ) = number of SSYT of shape λ and weight μ.
///
/// Uses the horizontal-strip DP approach.
/// `KostkaCoefficient` in Mathematica.
pub fn kostka_coefficient(lambda: &Partition, mu: &Partition) -> i64 {
    if lambda.size() != mu.size() {
        return 0;
    }
    if lambda.is_empty() && mu.is_empty() {
        return 1;
    }
    let weight = mu.parts();

    let mut current: BTreeMap<Partition, i64> = BTreeMap::new();
    current.insert(Partition::empty(), 1);

    for &w in weight.iter().rev() {
        let mut next: BTreeMap<Partition, i64> = BTreeMap::new();
        for (shape, count) in &current {
            for nu in add_horizontal_strip(shape, w, lambda) {
                *next.entry(nu).or_insert(0) += count;
            }
        }
        current = next;
    }

    current.get(lambda).copied().unwrap_or(0)
}

/// All partitions obtained by adding a horizontal strip of size k to `inner`,
/// that are contained in `outer`.
fn add_horizontal_strip(inner: &Partition, k: u32, outer: &Partition) -> Vec<Partition> {
    if k == 0 {
        return vec![inner.clone()];
    }
    let max_parts = outer.num_parts();
    let mut results = Vec::new();
    add_hstrip_helper(inner, outer, k, 0, max_parts, &mut vec![], &mut results);
    results
}

fn add_hstrip_helper(
    inner: &Partition,
    outer: &Partition,
    remaining: u32,
    row: usize,
    max_rows: usize,
    current: &mut Vec<u32>,
    results: &mut Vec<Partition>,
) {
    if row >= max_rows {
        if remaining == 0 {
            let mut parts = current.clone();
            for i in row..inner.num_parts() {
                parts.push(inner.part(i));
            }
            results.push(Partition::from_sorted(parts));
        }
        return;
    }

    let inner_r = inner.part(row);
    let outer_r = outer.part(row);

    let prev_new = if current.is_empty() {
        u32::MAX
    } else {
        current[current.len() - 1]
    };
    let hstrip_upper = if row == 0 {
        u32::MAX
    } else {
        inner.part(row - 1)
    };

    let lower = inner_r;
    let upper = outer_r.min(prev_new).min(hstrip_upper);

    if lower > upper {
        if inner_r <= prev_new && inner_r <= hstrip_upper && inner_r <= outer_r {
            current.push(inner_r);
            add_hstrip_helper(inner, outer, remaining, row + 1, max_rows, current, results);
            current.pop();
        }
        return;
    }

    let max_add = (upper - lower).min(remaining);

    for add in 0..=max_add {
        let new_r = inner_r + add;
        if new_r > prev_new || new_r > hstrip_upper || new_r > outer_r {
            continue;
        }
        current.push(new_r);
        add_hstrip_helper(
            inner,
            outer,
            remaining - add,
            row + 1,
            max_rows,
            current,
            results,
        );
        current.pop();
    }
}

/// Compute the full Kostka matrix for partitions of degree n.
/// Returns `(partitions, matrix)` where `matrix[i][j] = K(λ_i, μ_j)`.
pub fn kostka_matrix(n: u32) -> (Vec<Partition>, Vec<Vec<i64>>) {
    let partitions = Partition::all_of_size(n);
    let k = partitions.len();
    let mut matrix = vec![vec![0i64; k]; k];
    for i in 0..k {
        for j in 0..k {
            matrix[i][j] = kostka_coefficient(&partitions[i], &partitions[j]);
        }
    }
    (partitions, matrix)
}

/// Inverse Kostka matrix K^{-1}: used for Schur -> monomial conversion.
pub fn inverse_kostka_matrix(n: u32) -> (Vec<Partition>, Vec<Vec<i64>>) {
    let partitions = Partition::all_of_size_dominance_order(n);
    let k = partitions.len();

    let mut kmat = vec![vec![0i64; k]; k];
    for i in 0..k {
        for j in 0..k {
            kmat[i][j] = kostka_coefficient(&partitions[i], &partitions[j]);
        }
    }

    let inv = invert_unitriangular(&kmat);
    (partitions, inv)
}

fn invert_unitriangular(mat: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let n = mat.len();
    let mut inv = vec![vec![0i64; n]; n];

    for i in 0..n {
        inv[i][i] = 1;
    }

    for i in (0..n).rev() {
        for j in (i + 1)..n {
            let mut sum = 0i64;
            for k in (i + 1)..=j {
                sum += mat[i][k] * inv[k][j];
            }
            inv[i][j] = -sum;
        }
    }

    inv
}

/// Character of the symmetric group S_n.
/// χ^λ(μ) = character of irrep λ evaluated at conjugacy class of type μ.
///
/// Uses the Murnaghan-Nakayama rule.
pub fn sn_character(lambda: &Partition, mu: &Partition) -> i64 {
    if lambda.size() != mu.size() {
        return 0;
    }
    if lambda.is_empty() && mu.is_empty() {
        return 1;
    }

    let weight = mu.parts();
    if weight.is_empty() {
        return if lambda.is_empty() { 1 } else { 0 };
    }

    let strip_size = weight[0];
    let remaining_mu = Partition::from_sorted(weight[1..].to_vec());

    let mut result = 0i64;
    for (nu, height) in remove_border_strip(lambda, strip_size) {
        let sign = if height % 2 == 0 { 1 } else { -1 };
        result += sign * sn_character(&nu, &remaining_mu);
    }

    result
}

/// Remove all possible border strips (rim hooks) of size k from partition lambda.
/// Returns (resulting_partition, height) pairs.
fn remove_border_strip(lambda: &Partition, k: u32) -> Vec<(Partition, u32)> {
    if k == 0 {
        return vec![(lambda.clone(), 0)];
    }
    let ell = lambda.num_parts();
    if ell == 0 {
        return vec![];
    }

    let mut results = Vec::new();

    for top in 0..ell {
        for bot in top..ell {
            let mu_bot_signed = lambda.part(top) as i64 + (bot as i64 - top as i64) - k as i64;
            if mu_bot_signed < 0 {
                continue;
            }
            let mu_bot = mu_bot_signed as u32;

            let lambda_bot_next = if bot + 1 < ell {
                lambda.part(bot + 1)
            } else {
                0
            };
            if mu_bot < lambda_bot_next {
                continue;
            }

            if mu_bot >= lambda.part(bot) {
                continue;
            }

            let mut parts = Vec::with_capacity(ell);
            for r in 0..ell {
                if r < top || r > bot {
                    parts.push(lambda.part(r));
                } else if r < bot {
                    parts.push(lambda.part(r + 1) - 1);
                } else {
                    parts.push(mu_bot);
                }
            }
            parts.retain(|&x| x > 0);

            let height = (bot - top) as u32;
            results.push((Partition::from_sorted(parts), height));
        }
    }

    results
}

/// Compute the full S_n character table for partitions of n.
pub fn character_table(n: u32) -> (Vec<Partition>, Vec<Vec<i64>>) {
    let partitions = Partition::all_of_size(n);
    let k = partitions.len();
    let mut table = vec![vec![0i64; k]; k];
    for i in 0..k {
        for j in 0..k {
            table[i][j] = sn_character(&partitions[i], &partitions[j]);
        }
    }
    (partitions, table)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kostka_basic() {
        assert_eq!(
            kostka_coefficient(&Partition::new(vec![2, 1]), &Partition::new(vec![2, 1])),
            1
        );
        assert_eq!(
            kostka_coefficient(&Partition::new(vec![3]), &Partition::new(vec![1, 1, 1])),
            1
        );
        assert_eq!(
            kostka_coefficient(&Partition::new(vec![1, 1, 1]), &Partition::new(vec![3])),
            0
        );
    }

    #[test]
    fn test_kostka_matrix_3() {
        let (_parts, mat) = kostka_matrix(3);
        assert_eq!(mat[0][0], 1);
        assert_eq!(mat[0][2], 1);
        assert_eq!(mat[1][2], 2);
        assert_eq!(mat[2][0], 0);
    }

    #[test]
    fn test_sn_character() {
        assert_eq!(
            sn_character(&Partition::new(vec![2, 1]), &Partition::new(vec![3])),
            -1
        );
        assert_eq!(
            sn_character(&Partition::new(vec![2, 1]), &Partition::new(vec![2, 1])),
            0
        );
        assert_eq!(
            sn_character(&Partition::new(vec![2, 1]), &Partition::new(vec![1, 1, 1])),
            2
        );
    }

    #[test]
    fn test_character_table_orthogonality() {
        for n in 1..=5u32 {
            let (parts, table) = character_table(n);
            let k = parts.len();
            for i in 0..k {
                for j in 0..k {
                    let mut sum = 0f64;
                    for l in 0..k {
                        let z = parts[l].z_coefficient() as f64;
                        sum += (table[i][l] as f64) * (table[j][l] as f64) / z;
                    }
                    let expected = if i == j { 1.0 } else { 0.0 };
                    assert!((sum - expected).abs() < 1e-10, "n={}: ({},{})", n, i, j);
                }
            }
        }
    }
}
