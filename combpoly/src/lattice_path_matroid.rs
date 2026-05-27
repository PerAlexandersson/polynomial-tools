//! Lattice-path matroids from Dyck area sequences.
//!
//! The main constructor in this module uses the Dyck-path convention where each
//! peak gives an interval in `[n]`.  The lattice-path matroid bases are the
//! complete transversals of these peak intervals.  This is meant as a small,
//! exact computation layer for producing area-sequence-to-`h*` tables.  The
//! default Ehrhart counter uses the alcoved/cyclic-interval inequality
//! description of lattice-path matroid base polytopes; the all-subset
//! rank-inequality code remains available as a small-instance oracle.

use std::collections::BTreeSet;

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, ToPrimitive, Zero};

use crate::catalan::{
    area_sequence_to_dyck_word, for_each_area_sequence, is_area_sequence, peak_cliques,
};

type BigRat = Ratio<BigInt>;

/// A basis-list lattice-path matroid built from peak intervals.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LatticePathMatroid {
    ground_size: usize,
    intervals: Vec<(usize, usize)>,
    basis_masks: Vec<usize>,
}

/// One cyclic-interval inequality
/// `sum_{i in I} x_i <= rank(I)` in the base-polytope description.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CyclicIntervalInequality {
    /// First element of `I`, in 1-based cyclic notation.
    pub start: usize,
    /// Last element of `I`, in 1-based cyclic notation.  If `end < start`,
    /// the interval wraps through `ground_size`.
    pub end: usize,
    /// Bitmask of the elements in `I`, using bit `i-1` for element `i`.
    pub mask: usize,
    /// Matroid rank of `I`.
    pub rank: usize,
}

impl LatticePathMatroid {
    /// Build a lattice-path matroid from 1-based inclusive intervals.
    pub fn from_peak_intervals(
        ground_size: usize,
        intervals: Vec<(usize, usize)>,
    ) -> Result<Self, String> {
        if ground_size >= usize::BITS as usize {
            return Err(format!(
                "ground size {ground_size} is too large for the exact mask backend"
            ));
        }
        for &(start, end) in &intervals {
            if start == 0 || start > end || end > ground_size {
                return Err(format!("invalid peak interval ({start}, {end})"));
            }
        }

        let mut bases = BTreeSet::new();
        let mut used = 0usize;
        complete_transversal_masks(&intervals, 0, &mut used, &mut bases);
        if bases.is_empty() {
            return Err("peak intervals have no complete transversal".to_string());
        }

        Ok(Self {
            ground_size,
            intervals,
            basis_masks: bases.into_iter().collect(),
        })
    }

    /// Build from a Dyck area sequence.
    pub fn from_area_sequence(area: &[u8]) -> Result<Self, String> {
        if !is_area_sequence(area) {
            return Err(format!("not a valid Dyck area sequence: {area:?}"));
        }
        Self::from_peak_intervals(area.len(), peak_cliques(area))
    }

    pub fn ground_size(&self) -> usize {
        self.ground_size
    }

    pub fn rank(&self) -> usize {
        self.intervals.len()
    }

    pub fn intervals(&self) -> &[(usize, usize)] {
        &self.intervals
    }

    pub fn num_bases(&self) -> usize {
        self.basis_masks.len()
    }

    /// Number of bases from the Lindstrom--Gessel--Viennot determinant for
    /// the peak-interval presentation.
    pub fn num_bases_exact(&self) -> BigInt {
        lpm_basis_count_from_peak_intervals(&self.intervals)
    }

    /// Bases as 1-based sorted element lists.
    pub fn bases(&self) -> Vec<Vec<usize>> {
        self.basis_masks
            .iter()
            .map(|&mask| {
                (0..self.ground_size)
                    .filter_map(|i| ((mask & (1usize << i)) != 0).then_some(i + 1))
                    .collect()
            })
            .collect()
    }

    /// Affine dimension of the base polytope.
    pub fn base_polytope_dimension(&self) -> usize {
        affine_dimension_from_masks(&self.basis_masks, self.ground_size)
    }

    /// Cyclic-interval inequalities for the alcoved base-polytope description.
    ///
    /// Together with `sum_i x_i = rank`, `x_i >= 0`, and the singleton upper
    /// bounds included here, these inequalities define the LPM base polytope.
    /// Ranks are computed directly from the peak-interval presentation by
    /// bipartite matching, not by enumerating all bases.
    pub fn cyclic_interval_inequalities(&self) -> Vec<CyclicIntervalInequality> {
        cyclic_interval_inequalities(&self.intervals, self.ground_size)
    }

    /// Exact `h*` vector by Ehrhart inversion over the alcoved cyclic-interval
    /// inequalities.
    pub fn hstar(&self) -> Vec<BigInt> {
        hstar_from_cyclic_interval_inequalities_reciprocity(&self.intervals, self.ground_size)
    }

    /// Exact `h*` vector by generic Ehrhart inversion over all matroid rank
    /// inequalities.  This is exponential in `ground_size`; use only for small
    /// LPMs or as an oracle for the cyclic-interval implementation.
    pub fn hstar_generic_rank_table(&self) -> Vec<BigInt> {
        hstar_from_basis_masks(&self.basis_masks, self.ground_size)
    }
}

/// Compute the exact LPM `h*` vector from a Dyck area sequence using the
/// cyclic-interval inequality description, without enumerating bases.
pub fn lpm_hstar_from_area_sequence(area: &[u8]) -> Result<Vec<BigInt>, String> {
    if !is_area_sequence(area) {
        return Err(format!("not a valid Dyck area sequence: {area:?}"));
    }
    Ok(hstar_from_cyclic_interval_inequalities_reciprocity(
        &peak_cliques(area),
        area.len(),
    ))
}

/// Compute the cyclic-interval inequalities from a Dyck area sequence without
/// enumerating bases.
pub fn lpm_cyclic_interval_inequalities_from_area_sequence(
    area: &[u8],
) -> Result<Vec<CyclicIntervalInequality>, String> {
    if !is_area_sequence(area) {
        return Err(format!("not a valid Dyck area sequence: {area:?}"));
    }
    Ok(cyclic_interval_inequalities(
        &peak_cliques(area),
        area.len(),
    ))
}

/// Count bases of an LPM from its peak intervals using a determinant formula.
///
/// If the bases are the increasing sequences
/// `x_1 < ... < x_r` with `a_i <= x_i <= b_i`, set
/// `alpha_i = a_i - i` and `beta_i = b_i - i`.  The number of such bases is
/// the determinant of
/// `binom(beta_i - alpha_j + 1, j - i + 1)`.
pub fn lpm_basis_count_from_peak_intervals(intervals: &[(usize, usize)]) -> BigInt {
    let rank = intervals.len();
    if rank == 0 {
        return BigInt::one();
    }

    let mut matrix = vec![vec![BigInt::zero(); rank]; rank];
    for i in 0..rank {
        let beta_i = intervals[i].1 as i128 - (i + 1) as i128;
        for j in 0..rank {
            let alpha_j = intervals[j].0 as i128 - (j + 1) as i128;
            matrix[i][j] = binomial_bigint_signed(beta_i - alpha_j + 1, j as i128 - i as i128 + 1);
        }
    }

    determinant_bigint(matrix)
}

/// Count bases from a Dyck area sequence without enumerating bases.
pub fn lpm_basis_count_from_area_sequence(area: &[u8]) -> Result<BigInt, String> {
    if !is_area_sequence(area) {
        return Err(format!("not a valid Dyck area sequence: {area:?}"));
    }
    Ok(lpm_basis_count_from_peak_intervals(&peak_cliques(area)))
}

/// Row in the small area-sequence-to-`h*` table.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LpmHstarRow {
    pub area_sequence: Vec<u8>,
    pub dyck_word: String,
    pub intervals: Vec<(usize, usize)>,
    pub rank: usize,
    pub num_bases: BigInt,
    pub dimension: usize,
    pub hstar: Vec<BigInt>,
}

/// Row comparing `h*(1)` with the snake-contact volume formula.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LpmSnakeContactVolumeRow {
    pub area_sequence: Vec<u8>,
    pub dyck_word: String,
    pub intervals: Vec<(usize, usize)>,
    pub rank: usize,
    pub hstar: Vec<BigInt>,
    pub hstar_volume: BigInt,
    pub snake_volume: BigInt,
}

/// Compute exact `h*` data for all Dyck area sequences of a fixed semilength.
pub fn lpm_hstar_table(semilength: usize) -> Result<Vec<LpmHstarRow>, String> {
    let mut rows = Vec::new();
    let mut first_error = None;

    for_each_area_sequence(semilength, |area| {
        if first_error.is_some() {
            return;
        }
        if !is_area_sequence(area) {
            first_error = Some(format!("not a valid Dyck area sequence: {area:?}"));
            return;
        }
        let intervals = peak_cliques(area);
        if let Some(err) = validate_peak_intervals(semilength, &intervals) {
            first_error = Some(err);
            return;
        }
        let (dimension, hstar) = cyclic_interval_hstar_data_reciprocity(&intervals, semilength);
        let num_bases = lpm_basis_count_from_peak_intervals(&intervals);
        rows.push(LpmHstarRow {
            area_sequence: area.to_vec(),
            dyck_word: area_sequence_to_dyck_word(area),
            rank: intervals.len(),
            intervals,
            num_bases,
            dimension,
            hstar,
        });
    });

    if let Some(err) = first_error {
        Err(err)
    } else {
        Ok(rows)
    }
}

/// Compute the volume predicted by the Dyck snake-contact formula.
pub fn lpm_snake_contact_volume_from_area_sequence(area: &[u8]) -> Result<BigInt, String> {
    if !is_area_sequence(area) {
        return Err(format!("not a valid Dyck area sequence: {area:?}"));
    }
    Ok(lpm_snake_contact_volume_from_peak_intervals(&peak_cliques(
        area,
    )))
}

/// Compute the snake-contact volume from the peak-interval presentation.
pub fn lpm_snake_contact_volume_from_peak_intervals(intervals: &[(usize, usize)]) -> BigInt {
    if intervals.is_empty() {
        return BigInt::one();
    }

    let mut component_sums = Vec::new();
    let mut component_sizes = Vec::new();
    let mut start = 0usize;

    for end in 1..=intervals.len() {
        let split = end == intervals.len() || intervals[end - 1].1 < intervals[end].0;
        if split {
            let block_start = intervals[start].0;
            let block_end = intervals[end - 1].1;
            component_sizes.push(block_end - block_start);
            component_sums.push(snake_contact_component_sum(intervals, start, end));
            start = end;
        }
    }

    let total_size = component_sizes.iter().sum();
    let mut volume = multinomial_bigint(total_size, &component_sizes);
    for component_sum in component_sums {
        volume *= component_sum;
    }
    volume
}

/// Compare the snake-contact volume formula with `h*(1)` for a semilength.
pub fn lpm_snake_contact_volume_table(
    semilength: usize,
) -> Result<Vec<LpmSnakeContactVolumeRow>, String> {
    Ok(lpm_hstar_table(semilength)?
        .into_iter()
        .map(|row| {
            let hstar_volume = row
                .hstar
                .iter()
                .fold(BigInt::zero(), |acc, coeff| acc + coeff);
            let snake_volume = lpm_snake_contact_volume_from_peak_intervals(&row.intervals);
            LpmSnakeContactVolumeRow {
                area_sequence: row.area_sequence,
                dyck_word: row.dyck_word,
                intervals: row.intervals,
                rank: row.rank,
                hstar: row.hstar,
                hstar_volume,
                snake_volume,
            }
        })
        .collect())
}

/// Rank-two uniform matroid formula.
pub fn uniform_rank_two_hstar(n: usize) -> Vec<BigInt> {
    let mut coeffs = vec![BigInt::zero(); n / 2 + 1];
    coeffs[0] = BigInt::one();
    if coeffs.len() > 1 {
        coeffs[1] = BigInt::from(binomial_usize(n, 2) - n as u128);
    }
    for (k, coeff) in coeffs.iter_mut().enumerate().skip(2) {
        *coeff = BigInt::from(binomial_usize(n, 2 * k));
    }
    trim_bigint_trailing_zeros(coeffs)
}

/// Rank-two Schubert matroid formula `S_{n,ell} = h*(M_n[1,ell])`.
pub fn schubert_rank_two_hstar(n: usize, ell: usize) -> Vec<BigInt> {
    let mut coeffs = vec![BigInt::zero(); n / 2 + 1];
    coeffs[0] = BigInt::one();
    if coeffs.len() > 1 {
        let h1 = ell * (n - 2) - ell * (ell - 1) / 2 - 1;
        coeffs[1] = BigInt::from(h1);
    }
    for (p, coeff) in coeffs.iter_mut().enumerate().skip(2) {
        let value = (p - 1..ell)
            .map(|i| binomial_usize(i, p - 1) * binomial_usize(n - 1 - i, p))
            .sum::<u128>();
        *coeff = BigInt::from(value);
    }
    trim_bigint_trailing_zeros(coeffs)
}

/// Rank-two lattice-path formula
/// `h*(M_n[k,ell]) = S_{n,ell} + S_{n,n-k-1} - h*(U_{2,n})`.
pub fn rank_two_lpm_hstar(n: usize, k: usize, ell: usize) -> Result<Vec<BigInt>, String> {
    if !(1 <= k && k < ell && ell < n) {
        return Err(format!(
            "expected rank-two parameters 1 <= k < ell < n, got n={n}, k={k}, ell={ell}"
        ));
    }

    let a = schubert_rank_two_hstar(n, ell);
    let b = schubert_rank_two_hstar(n, n - k - 1);
    let u = uniform_rank_two_hstar(n);
    let length = a.len().max(b.len()).max(u.len());
    let mut coeffs = Vec::with_capacity(length);
    for i in 0..length {
        let mut value = a.get(i).cloned().unwrap_or_else(BigInt::zero);
        value += b.get(i).cloned().unwrap_or_else(BigInt::zero);
        value -= u.get(i).cloned().unwrap_or_else(BigInt::zero);
        coeffs.push(value);
    }
    Ok(trim_bigint_trailing_zeros(coeffs))
}

fn snake_contact_component_sum(intervals: &[(usize, usize)], start: usize, end: usize) -> BigInt {
    let block_start = intervals[start].0;
    let block_end = intervals[end - 1].1;
    let mut contacts = Vec::new();
    let mut sum = BigInt::zero();
    enumerate_component_contacts(
        intervals,
        end,
        start,
        block_start,
        block_start,
        block_end,
        &mut contacts,
        &mut sum,
    );
    sum
}

fn enumerate_component_contacts(
    intervals: &[(usize, usize)],
    end: usize,
    index: usize,
    previous: usize,
    block_start: usize,
    block_end: usize,
    contacts: &mut Vec<usize>,
    sum: &mut BigInt,
) {
    if index + 1 == end {
        *sum += fence_linear_extension_count(block_start, block_end, contacts);
        return;
    }

    let lower = intervals[index + 1].0.max(previous + 1);
    let upper = intervals[index].1.min(block_end - 1);
    for contact in lower..=upper {
        contacts.push(contact);
        enumerate_component_contacts(
            intervals,
            end,
            index + 1,
            contact,
            block_start,
            block_end,
            contacts,
            sum,
        );
        contacts.pop();
    }
}

fn fence_linear_extension_count(
    block_start: usize,
    block_end: usize,
    contacts: &[usize],
) -> BigInt {
    let size = block_end - block_start;
    let contact_positions: BTreeSet<_> = contacts
        .iter()
        .map(|&contact| contact - block_start)
        .collect();
    let descents: BTreeSet<_> = (1..size)
        .filter(|position| !contact_positions.contains(position))
        .collect();
    descent_set_permutation_count(size, &descents)
}

fn descent_set_permutation_count(size: usize, descents: &BTreeSet<usize>) -> BigInt {
    if size == 0 {
        return BigInt::one();
    }

    let mut dp = vec![BigInt::one()];
    for position in 1..size {
        let require_descent = descents.contains(&position);
        let mut next = vec![BigInt::zero(); position + 1];

        if require_descent {
            let mut suffix = BigInt::zero();
            for old_rank in (0..position).rev() {
                suffix += &dp[old_rank];
                next[old_rank] = suffix.clone();
            }
        } else {
            let mut prefix = BigInt::zero();
            for new_rank in 0..=position {
                next[new_rank] = prefix.clone();
                if new_rank < position {
                    prefix += &dp[new_rank];
                }
            }
        }

        dp = next;
    }

    dp.into_iter().sum()
}

fn multinomial_bigint(total: usize, parts: &[usize]) -> BigInt {
    let mut result = BigInt::one();
    let mut remaining = total;
    for &part in parts {
        result *= binomial_bigint(remaining, part);
        remaining -= part;
    }
    result
}

fn complete_transversal_masks(
    intervals: &[(usize, usize)],
    index: usize,
    used: &mut usize,
    bases: &mut BTreeSet<usize>,
) {
    if index == intervals.len() {
        bases.insert(*used);
        return;
    }

    let (start, end) = intervals[index];
    for label in start..=end {
        let bit = 1usize << (label - 1);
        if *used & bit == 0 {
            *used |= bit;
            complete_transversal_masks(intervals, index + 1, used, bases);
            *used &= !bit;
        }
    }
}

fn validate_peak_intervals(ground_size: usize, intervals: &[(usize, usize)]) -> Option<String> {
    if ground_size >= usize::BITS as usize {
        return Some(format!(
            "ground size {ground_size} is too large for the exact mask backend"
        ));
    }
    for &(start, end) in intervals {
        if start == 0 || start > end || end > ground_size {
            return Some(format!("invalid peak interval ({start}, {end})"));
        }
    }
    None
}

fn binomial_bigint_signed(n: i128, k: i128) -> BigInt {
    if n < 0 || k < 0 || k > n {
        return BigInt::zero();
    }
    binomial_bigint(n as usize, k as usize)
}

fn binomial_bigint(n: usize, k: usize) -> BigInt {
    if k > n {
        return BigInt::zero();
    }
    let k = k.min(n - k);
    let mut result = BigInt::one();
    for i in 0..k {
        result *= BigInt::from(n - i);
        result /= BigInt::from(i + 1);
    }
    result
}

fn determinant_bigint(mut matrix: Vec<Vec<BigInt>>) -> BigInt {
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
        let Some(pivot_row) = (k..n).find(|&i| !matrix[i][k].is_zero()) else {
            return BigInt::zero();
        };
        if pivot_row != k {
            matrix.swap(k, pivot_row);
            sign = -sign;
        }

        let pivot = matrix[k][k].clone();
        for i in k + 1..n {
            for j in k + 1..n {
                let numerator = &matrix[i][j] * &pivot - &matrix[i][k] * &matrix[k][j];
                matrix[i][j] = if k == 0 {
                    numerator
                } else {
                    numerator / &previous_pivot
                };
            }
        }

        previous_pivot = pivot;
    }

    sign * matrix[n - 1][n - 1].clone()
}

fn hstar_from_basis_masks(basis_masks: &[usize], ground_size: usize) -> Vec<BigInt> {
    let dimension = affine_dimension_from_masks(basis_masks, ground_size);
    let counts: Vec<BigInt> = (0..=dimension)
        .map(|t| BigInt::from(ehrhart_count(basis_masks, ground_size, t)))
        .collect();

    let mut hstar = Vec::with_capacity(dimension + 1);
    for i in 0..=dimension {
        let mut value = BigInt::zero();
        for k in 0..=i {
            let term = BigInt::from(binomial_usize(dimension + 1, k)) * &counts[i - k];
            if k % 2 == 0 {
                value += term;
            } else {
                value -= term;
            }
        }
        hstar.push(value);
    }
    trim_bigint_trailing_zeros(hstar)
}

fn hstar_from_cyclic_interval_inequalities_reciprocity(
    intervals: &[(usize, usize)],
    ground_size: usize,
) -> Vec<BigInt> {
    cyclic_interval_hstar_data_reciprocity(intervals, ground_size).1
}

#[cfg(test)]
fn cyclic_interval_hstar_data_positive(
    intervals: &[(usize, usize)],
    ground_size: usize,
) -> (usize, Vec<BigInt>) {
    let rank = intervals.len();
    let inequalities = cyclic_interval_inequalities(intervals, ground_size);
    let dimension = affine_dimension_from_cyclic_rank_equalities(intervals, ground_size);
    let counts: Vec<BigInt> = (0..=dimension)
        .map(|t| {
            BigInt::from(ehrhart_count_with_cyclic_intervals(
                ground_size,
                rank,
                t,
                &inequalities,
            ))
        })
        .collect();
    let hstar = hstar_from_ehrhart_counts(&counts, dimension);
    (dimension, hstar)
}

fn cyclic_interval_hstar_data_reciprocity(
    intervals: &[(usize, usize)],
    ground_size: usize,
) -> (usize, Vec<BigInt>) {
    let dimension = affine_dimension_from_cyclic_rank_equalities(intervals, ground_size);
    if dimension == 0 {
        return (0, vec![BigInt::one()]);
    }

    let rank = intervals.len();
    let inequalities = cyclic_interval_inequalities(intervals, ground_size);
    let strict_flags = strict_cyclic_interval_flags(intervals, ground_size, &inequalities);
    let strict_lower_bounds = strict_lower_bound_flags(intervals, ground_size);

    let num_positive = (dimension + 2) / 2;
    let num_negative = (dimension + 1) / 2;
    let sign = if dimension % 2 == 0 {
        BigInt::one()
    } else {
        -BigInt::one()
    };
    let mut points = Vec::with_capacity(dimension + 1);

    for t in 0..num_positive {
        let count = ehrhart_count_with_cyclic_intervals(ground_size, rank, t, &inequalities);
        points.push((t as i64, BigInt::from(count)));
    }

    for k in 1..=num_negative {
        let strict_count = strict_ehrhart_count_with_cyclic_intervals(
            ground_size,
            rank,
            k,
            &inequalities,
            &strict_flags,
            &strict_lower_bounds,
        );
        points.push((-(k as i64), &sign * BigInt::from(strict_count)));
    }

    let ehrhart_coeffs = lagrange_interpolation_big_rational(&points);
    let counts: Vec<BigInt> = (0..=dimension)
        .map(|t| big_rat_to_bigint(eval_big_rat_poly_at_i64(&ehrhart_coeffs, t as i64)))
        .collect();
    let hstar = hstar_from_ehrhart_counts(&counts, dimension);
    (dimension, hstar)
}

fn hstar_from_ehrhart_counts(counts: &[BigInt], dimension: usize) -> Vec<BigInt> {
    let mut hstar = Vec::with_capacity(dimension + 1);
    for i in 0..=dimension {
        let mut value = BigInt::zero();
        for k in 0..=i {
            let term = BigInt::from(binomial_usize(dimension + 1, k)) * &counts[i - k];
            if k % 2 == 0 {
                value += term;
            } else {
                value -= term;
            }
        }
        hstar.push(value);
    }
    trim_bigint_trailing_zeros(hstar)
}

fn lagrange_interpolation_big_rational(points: &[(i64, BigInt)]) -> Vec<BigRat> {
    if points.is_empty() {
        return vec![BigRat::zero()];
    }

    let mut result = vec![BigRat::zero(); points.len()];
    for (i, &(x_i, ref y_i)) in points.iter().enumerate() {
        let mut basis = vec![BigRat::one()];
        for (j, &(x_j, _)) in points.iter().enumerate() {
            if i == j {
                continue;
            }
            let denominator = BigRat::from_integer(BigInt::from(x_i - x_j));
            let constant = BigRat::from_integer(BigInt::from(-x_j)) / denominator.clone();
            let linear = BigRat::one() / denominator;
            let mut next = vec![BigRat::zero(); basis.len() + 1];
            for (degree, coeff) in basis.iter().enumerate() {
                next[degree] += coeff * &constant;
                next[degree + 1] += coeff * &linear;
            }
            basis = next;
        }

        let y = BigRat::from_integer(y_i.clone());
        for (degree, coeff) in basis.iter().enumerate() {
            result[degree] += coeff * &y;
        }
    }

    while result.len() > 1 && result.last().is_some_and(BigRat::is_zero) {
        result.pop();
    }
    result
}

fn eval_big_rat_poly_at_i64(coeffs: &[BigRat], x: i64) -> BigRat {
    let x = BigRat::from_integer(BigInt::from(x));
    coeffs
        .iter()
        .rev()
        .fold(BigRat::zero(), |acc, coeff| acc * &x + coeff)
}

fn big_rat_to_bigint(value: BigRat) -> BigInt {
    assert!(
        value.is_integer(),
        "expected integer Ehrhart value, got {value}"
    );
    value.to_integer()
}

fn cyclic_interval_inequalities(
    intervals: &[(usize, usize)],
    ground_size: usize,
) -> Vec<CyclicIntervalInequality> {
    let mut seen_masks = BTreeSet::new();
    let mut inequalities = Vec::new();

    for start0 in 0..ground_size {
        let mut mask = 0usize;
        for length in 1..ground_size {
            let end0 = (start0 + length - 1) % ground_size;
            mask |= 1usize << end0;
            if seen_masks.insert(mask) {
                inequalities.push(CyclicIntervalInequality {
                    start: start0 + 1,
                    end: end0 + 1,
                    mask,
                    rank: transversal_rank_on_mask(intervals, mask, ground_size),
                });
            }
        }
    }

    inequalities
}

fn affine_dimension_from_cyclic_rank_equalities(
    intervals: &[(usize, usize)],
    ground_size: usize,
) -> usize {
    if ground_size == 0 {
        return 0;
    }

    let rank = intervals.len();
    let full_mask = (1usize << ground_size) - 1;
    let inequalities = cyclic_interval_inequalities(intervals, ground_size);
    let mut rows = vec![vec![1i64; ground_size]];

    for inequality in inequalities {
        let complement = full_mask ^ inequality.mask;
        let complement_rank = transversal_rank_on_mask(intervals, complement, ground_size);
        if inequality.rank + complement_rank == rank {
            rows.push(indicator_row(inequality.mask, ground_size));
        }
    }

    ground_size - rational_rank(&rows)
}

fn strict_cyclic_interval_flags(
    intervals: &[(usize, usize)],
    ground_size: usize,
    inequalities: &[CyclicIntervalInequality],
) -> Vec<bool> {
    let rank = intervals.len();
    let full_mask = if ground_size == 0 {
        0
    } else {
        (1usize << ground_size) - 1
    };

    inequalities
        .iter()
        .map(|inequality| {
            let complement = full_mask ^ inequality.mask;
            let complement_rank = transversal_rank_on_mask(intervals, complement, ground_size);
            inequality.rank + complement_rank != rank
        })
        .collect()
}

fn strict_lower_bound_flags(intervals: &[(usize, usize)], ground_size: usize) -> Vec<bool> {
    (0..ground_size)
        .map(|i| transversal_rank_on_mask(intervals, 1usize << i, ground_size) > 0)
        .collect()
}

fn ehrhart_count_with_cyclic_intervals(
    ground_size: usize,
    rank: usize,
    dilation: usize,
    inequalities: &[CyclicIntervalInequality],
) -> u128 {
    let strict_flags = vec![false; inequalities.len()];
    let strict_lower_bounds = vec![false; ground_size];
    match interval_sum_bounds(
        ground_size,
        rank,
        dilation,
        inequalities,
        &strict_flags,
        &strict_lower_bounds,
    ) {
        Some(bounds) => count_sequences_with_interval_bounds(rank, dilation, &bounds),
        None => 0,
    }
}

fn strict_ehrhart_count_with_cyclic_intervals(
    ground_size: usize,
    rank: usize,
    dilation: usize,
    inequalities: &[CyclicIntervalInequality],
    strict_flags: &[bool],
    strict_lower_bounds: &[bool],
) -> u128 {
    match interval_sum_bounds(
        ground_size,
        rank,
        dilation,
        inequalities,
        strict_flags,
        strict_lower_bounds,
    ) {
        Some(bounds) => count_sequences_with_interval_bounds(rank, dilation, &bounds),
        None => 0,
    }
}

#[derive(Clone, Debug)]
struct IntervalSumBounds {
    lower: Vec<Vec<usize>>,
    upper: Vec<Vec<usize>>,
}

fn interval_sum_bounds(
    ground_size: usize,
    rank: usize,
    dilation: usize,
    inequalities: &[CyclicIntervalInequality],
    strict_flags: &[bool],
    strict_lower_bounds: &[bool],
) -> Option<IntervalSumBounds> {
    let total = rank * dilation;
    let mut lower = vec![vec![0usize; ground_size]; ground_size];
    let mut upper = vec![vec![total; ground_size]; ground_size];

    for (i, &strict) in strict_lower_bounds.iter().enumerate() {
        if strict {
            lower[i][i] = 1;
        }
    }

    for (inequality, &strict) in inequalities.iter().zip(strict_flags) {
        if inequality.start <= inequality.end {
            let start = inequality.start - 1;
            let end = inequality.end - 1;
            let rhs = dilation * inequality.rank;
            let bound = if strict { rhs.checked_sub(1)? } else { rhs };
            upper[start][end] = upper[start][end].min(bound);
            if lower[start][end] > upper[start][end] {
                return None;
            }
        } else {
            let start = inequality.end;
            let end = inequality.start - 2;
            let mut bound = dilation * (rank - inequality.rank);
            if strict {
                bound += 1;
            }
            lower[start][end] = lower[start][end].max(bound);
            if lower[start][end] > upper[start][end] {
                return None;
            }
        }
    }

    Some(IntervalSumBounds { lower, upper })
}

fn count_sequences_with_interval_bounds(
    rank: usize,
    dilation: usize,
    bounds: &IntervalSumBounds,
) -> u128 {
    let ground_size = bounds.lower.len();
    let target = rank * dilation;
    if ground_size == 0 {
        return u128::from(target == 0);
    }

    let mut suffix_min = vec![0usize; ground_size + 1];
    let mut suffix_max = vec![0usize; ground_size + 1];
    for i in (0..ground_size).rev() {
        suffix_min[i] = suffix_min[i + 1] + bounds.lower[i][i];
        suffix_max[i] = suffix_max[i + 1] + bounds.upper[i][i].min(dilation);
    }
    if target < suffix_min[0] || target > suffix_max[0] {
        return 0;
    }

    let mut prefixes = vec![0usize];
    count_sequences_with_interval_bounds_rec(
        0,
        target,
        dilation,
        bounds,
        &suffix_min,
        &suffix_max,
        &mut prefixes,
    )
}

fn count_sequences_with_interval_bounds_rec(
    position: usize,
    target: usize,
    dilation: usize,
    bounds: &IntervalSumBounds,
    suffix_min: &[usize],
    suffix_max: &[usize],
    prefixes: &mut Vec<usize>,
) -> u128 {
    let ground_size = bounds.lower.len();
    let current_sum = *prefixes.last().unwrap();
    if position == ground_size {
        return u128::from(current_sum == target);
    }

    let remaining_min = suffix_min[position + 1];
    if current_sum + remaining_min > target {
        return 0;
    }

    let mut lo = bounds.lower[position][position];
    let hi_from_target = target - current_sum - remaining_min;
    let hi = bounds.upper[position][position]
        .min(dilation)
        .min(hi_from_target);

    let remaining_max = suffix_max[position + 1];
    if target > current_sum + remaining_max {
        lo = lo.max(target - current_sum - remaining_max);
    }
    if lo > hi {
        return 0;
    }

    let mut total = 0u128;
    for value in lo..=hi {
        let next_sum = current_sum + value;
        prefixes.push(next_sum);
        if interval_constraints_ending_at(position, bounds, prefixes) {
            total += count_sequences_with_interval_bounds_rec(
                position + 1,
                target,
                dilation,
                bounds,
                suffix_min,
                suffix_max,
                prefixes,
            );
        }
        prefixes.pop();
    }

    total
}

fn interval_constraints_ending_at(
    end: usize,
    bounds: &IntervalSumBounds,
    prefixes: &[usize],
) -> bool {
    let end_sum = prefixes[end + 1];
    for start in 0..=end {
        let interval_sum = end_sum - prefixes[start];
        if interval_sum < bounds.lower[start][end] || interval_sum > bounds.upper[start][end] {
            return false;
        }
    }
    true
}

#[cfg(test)]
fn ehrhart_count_with_cyclic_intervals_by_compositions(
    ground_size: usize,
    rank: usize,
    dilation: usize,
    inequalities: &[CyclicIntervalInequality],
) -> u128 {
    if dilation == 0 {
        return 1;
    }

    let mut total = 0u128;
    for_each_bounded_composition(ground_size, rank * dilation, dilation, |composition| {
        if satisfies_cyclic_interval_inequalities(composition, dilation, inequalities) {
            total += 1;
        }
    });

    total
}

#[cfg(test)]
fn strict_ehrhart_count_with_cyclic_intervals_by_compositions(
    ground_size: usize,
    rank: usize,
    dilation: usize,
    inequalities: &[CyclicIntervalInequality],
    strict_flags: &[bool],
    strict_lower_bounds: &[bool],
) -> u128 {
    if dilation == 0 {
        return 0;
    }

    let mut total = 0u128;
    for_each_bounded_composition(ground_size, rank * dilation, dilation, |composition| {
        if satisfies_strict_cyclic_interval_inequalities(
            composition,
            dilation,
            inequalities,
            strict_flags,
            strict_lower_bounds,
        ) {
            total += 1;
        }
    });

    total
}

#[cfg(test)]
fn satisfies_strict_cyclic_interval_inequalities(
    composition: &[usize],
    dilation: usize,
    inequalities: &[CyclicIntervalInequality],
    strict_flags: &[bool],
    strict_lower_bounds: &[bool],
) -> bool {
    for (i, &strict) in strict_lower_bounds.iter().enumerate() {
        if strict && composition[i] == 0 {
            return false;
        }
    }

    for (inequality, &strict) in inequalities.iter().zip(strict_flags) {
        let lhs = sum_on_mask(composition, inequality.mask);
        let rhs = dilation * inequality.rank;
        if if strict { lhs >= rhs } else { lhs > rhs } {
            return false;
        }
    }

    true
}

#[cfg(test)]
fn satisfies_cyclic_interval_inequalities(
    composition: &[usize],
    dilation: usize,
    inequalities: &[CyclicIntervalInequality],
) -> bool {
    inequalities
        .iter()
        .all(|ineq| sum_on_mask(composition, ineq.mask) <= dilation * ineq.rank)
}

#[cfg(test)]
fn sum_on_mask(composition: &[usize], mask: usize) -> usize {
    composition
        .iter()
        .enumerate()
        .filter_map(|(i, &value)| ((mask & (1usize << i)) != 0).then_some(value))
        .sum()
}

fn ehrhart_count(basis_masks: &[usize], ground_size: usize, dilation: usize) -> u128 {
    if dilation == 0 {
        return 1;
    }

    let rank = basis_masks[0].count_ones() as usize;
    let rank_table = matroid_rank_table(basis_masks, ground_size);
    let mut total = 0u128;

    for_each_bounded_composition(ground_size, rank * dilation, dilation, |composition| {
        if satisfies_base_polytope_inequalities(composition, dilation, &rank_table) {
            total += 1;
        }
    });

    total
}

fn transversal_rank_on_mask(
    intervals: &[(usize, usize)],
    mask: usize,
    ground_size: usize,
) -> usize {
    let mut matched_interval_for_element = vec![None; ground_size];
    let mut rank = 0usize;

    for interval_index in 0..intervals.len() {
        let mut seen_elements = vec![false; ground_size];
        if augment_transversal_matching(
            interval_index,
            intervals,
            mask,
            &mut matched_interval_for_element,
            &mut seen_elements,
        ) {
            rank += 1;
        }
    }

    rank
}

fn augment_transversal_matching(
    interval_index: usize,
    intervals: &[(usize, usize)],
    mask: usize,
    matched_interval_for_element: &mut [Option<usize>],
    seen_elements: &mut [bool],
) -> bool {
    let (start, end) = intervals[interval_index];
    for label in start..=end {
        let element = label - 1;
        let bit = 1usize << element;
        if (mask & bit) == 0 || seen_elements[element] {
            continue;
        }
        seen_elements[element] = true;
        match matched_interval_for_element[element] {
            None => {
                matched_interval_for_element[element] = Some(interval_index);
                return true;
            }
            Some(previous_interval) => {
                if augment_transversal_matching(
                    previous_interval,
                    intervals,
                    mask,
                    matched_interval_for_element,
                    seen_elements,
                ) {
                    matched_interval_for_element[element] = Some(interval_index);
                    return true;
                }
            }
        }
    }

    false
}

fn satisfies_base_polytope_inequalities(
    composition: &[usize],
    dilation: usize,
    rank_table: &[usize],
) -> bool {
    let num_masks = rank_table.len();
    let mut sums = vec![0usize; num_masks];
    for mask in 1..num_masks {
        let bit = mask & mask.wrapping_neg();
        let i = bit.trailing_zeros() as usize;
        sums[mask] = sums[mask ^ bit] + composition[i];
    }

    for mask in 1..num_masks - 1 {
        if sums[mask] > dilation * rank_table[mask] {
            return false;
        }
    }
    true
}

fn matroid_rank_table(basis_masks: &[usize], ground_size: usize) -> Vec<usize> {
    let num_masks = 1usize << ground_size;
    let mut ranks = vec![0usize; num_masks];
    for (mask, rank) in ranks.iter_mut().enumerate() {
        *rank = basis_masks
            .iter()
            .map(|basis| (mask & basis).count_ones() as usize)
            .max()
            .unwrap_or(0);
    }
    ranks
}

fn for_each_bounded_composition<F>(length: usize, total: usize, upper: usize, mut visitor: F)
where
    F: FnMut(&[usize]),
{
    let mut current = vec![0usize; length];
    bounded_composition_rec(0, total, upper, &mut current, &mut visitor);
}

fn bounded_composition_rec<F>(
    position: usize,
    remaining: usize,
    upper: usize,
    current: &mut [usize],
    visitor: &mut F,
) where
    F: FnMut(&[usize]),
{
    if position + 1 == current.len() {
        if remaining <= upper {
            current[position] = remaining;
            visitor(current);
        }
        return;
    }

    let remaining_positions = current.len() - position - 1;
    let lo = remaining.saturating_sub(upper * remaining_positions);
    let hi = remaining.min(upper);
    for value in lo..=hi {
        current[position] = value;
        bounded_composition_rec(position + 1, remaining - value, upper, current, visitor);
    }
}

fn affine_dimension_from_masks(basis_masks: &[usize], ground_size: usize) -> usize {
    if basis_masks.len() <= 1 {
        return 0;
    }

    let first = indicator_row(basis_masks[0], ground_size);
    let rows: Vec<Vec<i64>> = basis_masks[1..]
        .iter()
        .map(|&mask| {
            indicator_row(mask, ground_size)
                .into_iter()
                .zip(first.iter().copied())
                .map(|(x, y)| x - y)
                .collect()
        })
        .collect();
    rational_rank(&rows)
}

fn indicator_row(mask: usize, ground_size: usize) -> Vec<i64> {
    (0..ground_size)
        .map(|i| if mask & (1usize << i) == 0 { 0 } else { 1 })
        .collect()
}

fn rational_rank(rows: &[Vec<i64>]) -> usize {
    if rows.is_empty() {
        return 0;
    }

    let mut matrix: Vec<Vec<BigRat>> = rows
        .iter()
        .map(|row| {
            row.iter()
                .map(|&x| BigRat::from_integer(BigInt::from(x)))
                .collect()
        })
        .collect();
    let mut rank = 0usize;
    let mut col = 0usize;
    let cols = matrix[0].len();

    while rank < matrix.len() && col < cols {
        let Some(pivot) = (rank..matrix.len()).find(|&i| !matrix[i][col].is_zero()) else {
            col += 1;
            continue;
        };
        matrix.swap(rank, pivot);

        let scale = matrix[rank][col].clone();
        for entry in &mut matrix[rank][col..] {
            *entry /= scale.clone();
        }

        for i in 0..matrix.len() {
            if i == rank || matrix[i][col].is_zero() {
                continue;
            }
            let factor = matrix[i][col].clone();
            for j in col..cols {
                let pivot_entry = matrix[rank][j].clone();
                matrix[i][j] -= factor.clone() * pivot_entry;
            }
        }

        rank += 1;
        col += 1;
    }

    rank
}

fn binomial_usize(n: usize, k: usize) -> u128 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result = 1u128;
    for i in 0..k {
        result = result * (n - i) as u128 / (i + 1) as u128;
    }
    result
}

fn trim_bigint_trailing_zeros(mut coeffs: Vec<BigInt>) -> Vec<BigInt> {
    while coeffs.len() > 1 && coeffs.last().is_some_and(BigInt::is_zero) {
        coeffs.pop();
    }
    coeffs
}

/// Convert a small `h*` vector to `i64` coefficients for use with
/// `polynomial-tools` checks.
pub fn hstar_to_i64(coeffs: &[BigInt]) -> Option<Vec<i64>> {
    coeffs.iter().map(ToPrimitive::to_i64).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ints(values: &[i64]) -> Vec<BigInt> {
        values.iter().copied().map(BigInt::from).collect()
    }

    #[test]
    fn test_peak_interval_lpm_bases() {
        let area = vec![0, 1, 1, 2, 3, 3, 2, 0];
        let matroid = LatticePathMatroid::from_area_sequence(&area).unwrap();
        assert_eq!(
            matroid.intervals(),
            &[(1, 2), (2, 5), (3, 6), (5, 7), (8, 8)]
        );
        assert_eq!(matroid.rank(), 5);
        assert!(matroid.num_bases() > 0);
        assert_eq!(matroid.num_bases_exact(), BigInt::from(matroid.num_bases()));
        assert!(matroid
            .bases()
            .iter()
            .all(|basis| basis.len() == matroid.rank()));
    }

    #[test]
    fn test_tiny_area_hstars() {
        assert_eq!(
            LatticePathMatroid::from_area_sequence(&[0, 0])
                .unwrap()
                .hstar(),
            ints(&[1])
        );
        assert_eq!(
            LatticePathMatroid::from_area_sequence(&[0, 1])
                .unwrap()
                .hstar(),
            ints(&[1])
        );
        assert_eq!(
            LatticePathMatroid::from_area_sequence(&[0, 1, 2])
                .unwrap()
                .hstar(),
            ints(&[1])
        );
    }

    #[test]
    fn test_cyclic_interval_inequalities_example() {
        let matroid = LatticePathMatroid::from_area_sequence(&[0, 1, 1]).unwrap();
        let inequalities = matroid.cyclic_interval_inequalities();
        assert_eq!(inequalities.len(), 6);
        assert!(inequalities
            .iter()
            .any(|ineq| (ineq.start, ineq.end, ineq.rank) == (1, 2, 2)));
        assert!(inequalities
            .iter()
            .any(|ineq| (ineq.start, ineq.end, ineq.rank) == (3, 1, 2)));
    }

    #[test]
    fn test_cyclic_interval_hstar_matches_generic_oracle() {
        for n in 0..=5 {
            for_each_area_sequence(n, |area| {
                let matroid = LatticePathMatroid::from_area_sequence(area).unwrap();
                let (dimension, hstar) =
                    cyclic_interval_hstar_data_positive(matroid.intervals(), matroid.ground_size());
                let (reciprocity_dimension, reciprocity_hstar) =
                    cyclic_interval_hstar_data_reciprocity(
                        matroid.intervals(),
                        matroid.ground_size(),
                    );
                assert_eq!(dimension, matroid.base_polytope_dimension());
                assert_eq!(reciprocity_dimension, dimension);
                assert_eq!(hstar, matroid.hstar_generic_rank_table());
                assert_eq!(reciprocity_hstar, hstar);
                assert_eq!(matroid.hstar(), matroid.hstar_generic_rank_table());
                assert_eq!(
                    lpm_hstar_from_area_sequence(area).unwrap(),
                    reciprocity_hstar
                );
            });
        }
    }

    #[test]
    fn test_interval_bound_dp_matches_composition_counter() {
        for n in 0..=5 {
            for_each_area_sequence(n, |area| {
                let matroid = LatticePathMatroid::from_area_sequence(area).unwrap();
                let ground_size = matroid.ground_size();
                let rank = matroid.rank();
                let dimension =
                    affine_dimension_from_cyclic_rank_equalities(matroid.intervals(), ground_size);
                let inequalities = matroid.cyclic_interval_inequalities();
                let strict_flags =
                    strict_cyclic_interval_flags(matroid.intervals(), ground_size, &inequalities);
                let strict_lower_bounds =
                    strict_lower_bound_flags(matroid.intervals(), ground_size);

                for dilation in 0..=dimension {
                    assert_eq!(
                        ehrhart_count_with_cyclic_intervals(
                            ground_size,
                            rank,
                            dilation,
                            &inequalities
                        ),
                        ehrhart_count_with_cyclic_intervals_by_compositions(
                            ground_size,
                            rank,
                            dilation,
                            &inequalities
                        ),
                        "weak count mismatch for area={area:?}, dilation={dilation}",
                    );
                    if dilation > 0 {
                        assert_eq!(
                            strict_ehrhart_count_with_cyclic_intervals(
                                ground_size,
                                rank,
                                dilation,
                                &inequalities,
                                &strict_flags,
                                &strict_lower_bounds,
                            ),
                            strict_ehrhart_count_with_cyclic_intervals_by_compositions(
                                ground_size,
                                rank,
                                dilation,
                                &inequalities,
                                &strict_flags,
                                &strict_lower_bounds,
                            ),
                            "strict count mismatch for area={area:?}, dilation={dilation}",
                        );
                    }
                }
            });
        }
    }

    #[test]
    fn test_determinant_basis_count_matches_basis_enumeration() {
        for n in 0..=7 {
            for_each_area_sequence(n, |area| {
                let matroid = LatticePathMatroid::from_area_sequence(area).unwrap();
                assert_eq!(
                    lpm_basis_count_from_area_sequence(area).unwrap(),
                    BigInt::from(matroid.num_bases())
                );
                assert_eq!(
                    lpm_basis_count_from_peak_intervals(matroid.intervals()),
                    BigInt::from(matroid.num_bases())
                );
            });
        }
    }

    #[test]
    fn test_snake_contact_volume_examples() {
        assert_eq!(
            lpm_snake_contact_volume_from_area_sequence(&[0, 1, 2, 3, 3]).unwrap(),
            BigInt::from(11)
        );
        assert_eq!(
            lpm_snake_contact_volume_from_area_sequence(&[0, 1, 1, 2, 2]).unwrap(),
            BigInt::from(8)
        );
        assert_eq!(
            lpm_snake_contact_volume_from_area_sequence(&[0, 1, 2, 2, 1, 2, 2]).unwrap(),
            BigInt::from(160)
        );
    }

    #[test]
    fn test_snake_contact_volume_matches_hstar_volume() {
        for n in 0..=7 {
            for row in lpm_snake_contact_volume_table(n).unwrap() {
                assert_eq!(
                    row.snake_volume, row.hstar_volume,
                    "snake-contact volume mismatch for area={:?}",
                    row.area_sequence,
                );
            }
        }
    }

    #[test]
    fn test_uniform_rank_two_formula() {
        assert_eq!(uniform_rank_two_hstar(6), ints(&[1, 9, 15, 1]));
    }

    #[test]
    fn test_rank_two_lpm_formula_example() {
        assert_eq!(rank_two_lpm_hstar(7, 2, 3).unwrap(), ints(&[1, 10, 18, 4]));
    }

    #[test]
    fn test_lpm_hstar_table_counts() {
        assert_eq!(lpm_hstar_table(0).unwrap().len(), 1);
        assert_eq!(lpm_hstar_table(1).unwrap().len(), 1);
        assert_eq!(lpm_hstar_table(2).unwrap().len(), 2);
        assert_eq!(lpm_hstar_table(3).unwrap().len(), 5);
    }
}
