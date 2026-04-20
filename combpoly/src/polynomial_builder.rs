//! Polynomial construction from combinatorial objects.

use crate::statistics::{self, Stat};

/// Build a polynomial from a collection of objects and a statistic.
/// Returns the coefficient vector where `coeffs[k] = #{w : stat(w) = k}`.
pub fn build_generating_polynomial(objects: &[Vec<u8>], stat: Stat) -> Vec<i64> {
    if objects.is_empty() {
        return vec![0];
    }

    let values: Vec<usize> = objects
        .iter()
        .map(|w| statistics::compute(w, stat))
        .collect();
    let max_val = values.iter().copied().max().unwrap_or(0);

    let mut coeffs = vec![0i64; max_val + 1];
    for v in values {
        coeffs[v] += 1;
    }
    coeffs
}
