//! Shifted Littlewood--Richardson coefficients.
//!
//! This module implements the interval recurrence for shifted Schur structure
//! constants.  It is intended as a correctness-first library home for the
//! experimental recurrence, not as a replacement for optimized ordinary LR
//! coefficient engines.

use std::collections::HashMap;

use num_bigint::BigInt;
use num_traits::{One, Zero};
use sym_poly_core::Partition;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ShiftedLrError {
    ArithmeticOverflow,
    NonIntegralRecurrence,
    ZeroEvaluationDenominator,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ShiftedLrStats {
    pub coefficient: BigInt,
    pub fixed: Partition,
    pub lower: Partition,
    pub upper: Partition,
    pub interval_partitions: usize,
    pub base_states: usize,
    pub total_pair_states: usize,
    pub peak_layer_states: usize,
    pub layer_states: Vec<usize>,
}

#[derive(Clone, Debug)]
struct IntervalData {
    parts: Vec<Vec<u32>>,
    by_rank: Vec<Vec<usize>>,
    add_neighbors: Vec<Vec<usize>>,
    remove_neighbors: Vec<Vec<usize>>,
    lower_index: usize,
    upper_index: usize,
}

/// Compute the shifted Littlewood--Richardson coefficient
/// `c^upper_{fixed, lower}`.
///
/// The coefficient is for the shifted Schur product.  When
/// `|upper|-|lower| = |fixed|`, this top-degree coefficient agrees with the
/// ordinary Littlewood--Richardson coefficient.
pub fn shifted_littlewood_richardson_coefficient(
    fixed: &Partition,
    lower: &Partition,
    upper: &Partition,
) -> Result<BigInt, ShiftedLrError> {
    Ok(shifted_littlewood_richardson_stats(fixed, lower, upper)?.coefficient)
}

/// Compute a shifted Littlewood--Richardson coefficient and interval-DP stats.
pub fn shifted_littlewood_richardson_stats(
    fixed: &Partition,
    lower: &Partition,
    upper: &Partition,
) -> Result<ShiftedLrStats, ShiftedLrError> {
    if !lower.partition_less_equal(upper) {
        return Ok(zero_stats(fixed, lower, upper));
    }
    let gap = usize::try_from(upper.size() - lower.size())
        .map_err(|_| ShiftedLrError::ArithmeticOverflow)?;
    let fixed_size =
        usize::try_from(fixed.size()).map_err(|_| ShiftedLrError::ArithmeticOverflow)?;
    if gap > fixed_size {
        return Ok(zero_stats(fixed, lower, upper));
    }

    let interval = build_interval(lower, upper)?;
    let base_states = interval.parts.len();
    let mut previous = HashMap::<(usize, usize), BigInt>::with_capacity(base_states);
    for (index, part) in interval.parts.iter().enumerate() {
        let point = Partition::from_sorted(part.clone());
        let value = shifted_schur_evaluation(fixed, &point)?;
        previous.insert((index, index), value);
    }

    let mut layer_states = vec![base_states];
    let mut total_pair_states = base_states;
    let mut peak_layer_states = base_states;

    for distance in 1..=gap {
        let mut current = HashMap::new();
        for upper_rank in distance..interval.by_rank.len() {
            let lower_rank = upper_rank - distance;
            for &upper_index in &interval.by_rank[upper_rank] {
                for &lower_index in &interval.by_rank[lower_rank] {
                    if !partition_leq_vec(
                        &interval.parts[lower_index],
                        &interval.parts[upper_index],
                    ) {
                        continue;
                    }
                    let value = shifted_interval_next_value(
                        &interval,
                        &previous,
                        upper_index,
                        lower_index,
                        distance,
                    )?;
                    current.insert((upper_index, lower_index), value);
                }
            }
        }
        let count = current.len();
        total_pair_states = total_pair_states
            .checked_add(count)
            .ok_or(ShiftedLrError::ArithmeticOverflow)?;
        peak_layer_states = peak_layer_states.max(count);
        layer_states.push(count);
        previous = current;
    }

    let coefficient = previous
        .get(&(interval.upper_index, interval.lower_index))
        .cloned()
        .unwrap_or_else(BigInt::zero);

    Ok(ShiftedLrStats {
        coefficient,
        fixed: fixed.clone(),
        lower: lower.clone(),
        upper: upper.clone(),
        interval_partitions: interval.parts.len(),
        base_states,
        total_pair_states,
        peak_layer_states,
        layer_states,
    })
}

/// Evaluate the shifted Schur function `s^*_shape(point)`.
pub fn shifted_schur_evaluation(
    shape: &Partition,
    point: &Partition,
) -> Result<BigInt, ShiftedLrError> {
    let n = shape.num_parts().max(point.num_parts());
    if n == 0 {
        return Ok(BigInt::one());
    }
    let shape = pad_partition(shape, n);
    let point = pad_partition(point, n);
    let mut numerator = vec![vec![BigInt::zero(); n]; n];
    let mut denominator = vec![vec![BigInt::zero(); n]; n];

    for i in 0..n {
        let x = u64::from(point[i])
            .checked_add(u64::try_from(n - i - 1).map_err(|_| ShiftedLrError::ArithmeticOverflow)?)
            .ok_or(ShiftedLrError::ArithmeticOverflow)?;
        for j in 0..n {
            let degree_shift = n - j - 1;
            let num_degree = usize::try_from(shape[j])
                .map_err(|_| ShiftedLrError::ArithmeticOverflow)?
                .checked_add(degree_shift)
                .ok_or(ShiftedLrError::ArithmeticOverflow)?;
            numerator[i][j] = falling_factorial(x, num_degree);
            denominator[i][j] = falling_factorial(x, degree_shift);
        }
    }

    let denominator = determinant_bareiss(denominator);
    if denominator.is_zero() {
        return Err(ShiftedLrError::ZeroEvaluationDenominator);
    }
    let numerator = determinant_bareiss(numerator);
    if &numerator % &denominator != BigInt::zero() {
        return Err(ShiftedLrError::NonIntegralRecurrence);
    }
    Ok(numerator / denominator)
}

fn shifted_interval_next_value(
    interval: &IntervalData,
    previous: &HashMap<(usize, usize), BigInt>,
    upper_index: usize,
    lower_index: usize,
    distance: usize,
) -> Result<BigInt, ShiftedLrError> {
    let mut numerator = BigInt::zero();
    for &lower_plus in &interval.add_neighbors[lower_index] {
        if partition_leq_vec(&interval.parts[lower_plus], &interval.parts[upper_index]) {
            if let Some(value) = previous.get(&(upper_index, lower_plus)) {
                numerator += value;
            }
        }
    }
    for &upper_minus in &interval.remove_neighbors[upper_index] {
        if partition_leq_vec(&interval.parts[lower_index], &interval.parts[upper_minus]) {
            if let Some(value) = previous.get(&(upper_minus, lower_index)) {
                numerator -= value;
            }
        }
    }

    let divisor = BigInt::from(distance);
    if &numerator % &divisor != BigInt::zero() {
        return Err(ShiftedLrError::NonIntegralRecurrence);
    }
    Ok(numerator / divisor)
}

fn build_interval(lower: &Partition, upper: &Partition) -> Result<IntervalData, ShiftedLrError> {
    let len = lower.num_parts().max(upper.num_parts());
    let lower = pad_partition(lower, len);
    let upper = pad_partition(upper, len);
    let lower_size = partition_sum_usize(&lower)?;
    let upper_size = partition_sum_usize(&upper)?;
    let max_rank = upper_size
        .checked_sub(lower_size)
        .ok_or(ShiftedLrError::ArithmeticOverflow)?;

    let mut parts = Vec::new();
    let mut current = vec![0; len];
    generate_interval_partitions(0, &lower, &upper, &mut current, &mut parts);
    let mut index_by_part = HashMap::with_capacity(parts.len());
    for (index, part) in parts.iter().enumerate() {
        index_by_part.insert(part.clone(), index);
    }

    let mut by_rank = vec![Vec::new(); max_rank + 1];
    for (index, part) in parts.iter().enumerate() {
        let rank = partition_sum_usize(part)?
            .checked_sub(lower_size)
            .ok_or(ShiftedLrError::ArithmeticOverflow)?;
        by_rank[rank].push(index);
    }

    let lower_index = *index_by_part
        .get(&lower)
        .ok_or(ShiftedLrError::ArithmeticOverflow)?;
    let upper_index = *index_by_part
        .get(&upper)
        .ok_or(ShiftedLrError::ArithmeticOverflow)?;

    let mut add_neighbors = vec![Vec::new(); parts.len()];
    let mut remove_neighbors = vec![Vec::new(); parts.len()];
    for (index, part) in parts.iter().enumerate() {
        for row in 0..len {
            let mut next = part.clone();
            next[row] = next[row]
                .checked_add(1)
                .ok_or(ShiftedLrError::ArithmeticOverflow)?;
            if next[row] <= upper[row] && is_partition(&next) {
                if let Some(&next_index) = index_by_part.get(&next) {
                    add_neighbors[index].push(next_index);
                }
            }

            if part[row] > lower[row] {
                let mut previous = part.clone();
                previous[row] -= 1;
                if is_partition(&previous) {
                    if let Some(&previous_index) = index_by_part.get(&previous) {
                        remove_neighbors[index].push(previous_index);
                    }
                }
            }
        }
    }

    Ok(IntervalData {
        parts,
        by_rank,
        add_neighbors,
        remove_neighbors,
        lower_index,
        upper_index,
    })
}

fn generate_interval_partitions(
    row: usize,
    lower: &[u32],
    upper: &[u32],
    current: &mut [u32],
    out: &mut Vec<Vec<u32>>,
) {
    if row == current.len() {
        out.push(current.to_vec());
        return;
    }

    let previous_bound = if row == 0 {
        upper[row]
    } else {
        current[row - 1].min(upper[row])
    };
    let next_lower = lower.get(row + 1).copied().unwrap_or(0);
    let minimum = lower[row].max(next_lower);
    if minimum > previous_bound {
        return;
    }
    for value in (minimum..=previous_bound).rev() {
        current[row] = value;
        generate_interval_partitions(row + 1, lower, upper, current, out);
    }
}

fn determinant_bareiss(mut matrix: Vec<Vec<BigInt>>) -> BigInt {
    let n = matrix.len();
    if n == 0 {
        return BigInt::one();
    }
    if n == 1 {
        return matrix[0][0].clone();
    }

    let mut sign = BigInt::one();
    let mut previous_pivot = BigInt::one();
    for k in 0..n - 1 {
        let Some(pivot_row) = (k..n).find(|&row| !matrix[row][k].is_zero()) else {
            return BigInt::zero();
        };
        if pivot_row != k {
            matrix.swap(k, pivot_row);
            sign = -sign;
        }
        let pivot = matrix[k][k].clone();
        for i in k + 1..n {
            for j in k + 1..n {
                let value =
                    (&matrix[i][j] * &pivot - &matrix[i][k] * &matrix[k][j]) / &previous_pivot;
                matrix[i][j] = value;
            }
        }
        for row in matrix.iter_mut().take(n).skip(k + 1) {
            row[k] = BigInt::zero();
        }
        previous_pivot = pivot;
    }
    sign * matrix[n - 1][n - 1].clone()
}

fn falling_factorial(x: u64, degree: usize) -> BigInt {
    let mut value = BigInt::one();
    let x = BigInt::from(x);
    for offset in 0..degree {
        value *= &x - BigInt::from(offset);
    }
    value
}

fn zero_stats(fixed: &Partition, lower: &Partition, upper: &Partition) -> ShiftedLrStats {
    ShiftedLrStats {
        coefficient: BigInt::zero(),
        fixed: fixed.clone(),
        lower: lower.clone(),
        upper: upper.clone(),
        interval_partitions: 0,
        base_states: 0,
        total_pair_states: 0,
        peak_layer_states: 0,
        layer_states: Vec::new(),
    }
}

fn pad_partition(partition: &Partition, len: usize) -> Vec<u32> {
    let mut result = vec![0; len];
    for (index, &part) in partition.parts().iter().take(len).enumerate() {
        result[index] = part;
    }
    result
}

fn partition_sum_usize(partition: &[u32]) -> Result<usize, ShiftedLrError> {
    let mut sum = 0usize;
    for &part in partition {
        sum = sum
            .checked_add(usize::try_from(part).map_err(|_| ShiftedLrError::ArithmeticOverflow)?)
            .ok_or(ShiftedLrError::ArithmeticOverflow)?;
    }
    Ok(sum)
}

fn partition_leq_vec(left: &[u32], right: &[u32]) -> bool {
    let len = left.len().max(right.len());
    (0..len).all(|index| part_entry(left, index) <= part_entry(right, index))
}

fn part_entry(partition: &[u32], index: usize) -> u32 {
    partition.get(index).copied().unwrap_or(0)
}

fn is_partition(partition: &[u32]) -> bool {
    partition.windows(2).all(|window| window[0] >= window[1])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn part(parts: &[u32]) -> Partition {
        Partition::from_sorted(parts.to_vec())
    }

    #[test]
    fn shifted_schur_base_values_for_one_box() {
        let one = part(&[1]);
        assert_eq!(
            shifted_schur_evaluation(&one, &part(&[1])).unwrap(),
            BigInt::from(1)
        );
        assert_eq!(
            shifted_schur_evaluation(&one, &part(&[2])).unwrap(),
            BigInt::from(2)
        );
        assert_eq!(
            shifted_schur_evaluation(&one, &part(&[1, 1])).unwrap(),
            BigInt::from(2)
        );
    }

    #[test]
    fn shifted_base_layer_coefficients() {
        assert_eq!(
            shifted_littlewood_richardson_coefficient(&part(&[1]), &part(&[2]), &part(&[2]))
                .unwrap(),
            BigInt::from(2)
        );
        assert_eq!(
            shifted_littlewood_richardson_coefficient(&part(&[1]), &part(&[1, 1]), &part(&[1, 1]))
                .unwrap(),
            BigInt::from(2)
        );
    }

    #[test]
    fn top_degree_recovers_known_ordinary_lr_coefficients() {
        assert_eq!(
            shifted_littlewood_richardson_coefficient(&part(&[1]), &part(&[1]), &part(&[2]))
                .unwrap(),
            BigInt::from(1)
        );
        assert_eq!(
            shifted_littlewood_richardson_coefficient(&part(&[1]), &part(&[1]), &part(&[1, 1]))
                .unwrap(),
            BigInt::from(1)
        );
        assert_eq!(
            shifted_littlewood_richardson_coefficient(
                &part(&[2, 1]),
                &part(&[2, 1]),
                &part(&[3, 2, 1])
            )
            .unwrap(),
            BigInt::from(2)
        );
        assert_eq!(
            shifted_littlewood_richardson_coefficient(
                &part(&[2, 1]),
                &part(&[2, 1]),
                &part(&[4, 2])
            )
            .unwrap(),
            BigInt::from(1)
        );
    }
}
