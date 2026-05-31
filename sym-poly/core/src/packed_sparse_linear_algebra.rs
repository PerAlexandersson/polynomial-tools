//! Packed sparse row reduction over small prime fields.
//!
//! This is an experimental backend for primes `P < 256`. Coefficients are
//! stored as `u8`, and each RREF computation builds a small inverse table once.
//! The generic sparse module remains the default; this module exists so we can
//! measure whether byte-sized finite-field rows are worth wiring into larger
//! quotient computations.

use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackedSparseRow {
    pub cols: Vec<u32>,
    pub vals: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackedSparseRrefResult {
    pub num_cols: usize,
    pub rows: Vec<PackedSparseRow>,
    pub pivot_columns: Vec<usize>,
    pub rank: usize,
}

/// A quotient of a byte-packed sparse coordinate space by sparse relations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackedSparseQuotientSpace {
    pub ambient_dimension: usize,
    pub relation_rref: Vec<PackedSparseRow>,
    pub pivot_columns: Vec<usize>,
    pub free_columns: Vec<usize>,
    pub pivot_row_by_column: BTreeMap<usize, usize>,
    pub free_index_by_column: BTreeMap<usize, usize>,
}

impl PackedSparseRow {
    pub fn new<const P: u8, I>(num_cols: usize, entries: I) -> Self
    where
        I: IntoIterator<Item = (usize, u8)>,
    {
        assert!(P > 1, "field modulus must be greater than 1");
        assert!(
            num_cols <= u32::MAX as usize + 1,
            "packed sparse rows support at most u32::MAX + 1 columns"
        );
        let mut terms: BTreeMap<usize, u8> = BTreeMap::new();
        for (col, value) in entries {
            assert!(col < num_cols, "sparse row column out of range");
            let value = value % P;
            if value == 0 {
                continue;
            }
            let entry = terms.entry(col).or_insert(0);
            *entry = add_mod::<P>(*entry, value);
            if *entry == 0 {
                terms.remove(&col);
            }
        }

        let mut cols = Vec::with_capacity(terms.len());
        let mut vals = Vec::with_capacity(terms.len());
        for (col, value) in terms {
            cols.push(col as u32);
            vals.push(value);
        }
        Self { cols, vals }
    }

    pub fn is_empty(&self) -> bool {
        self.cols.is_empty()
    }

    pub fn len(&self) -> usize {
        self.cols.len()
    }

    pub fn leading_column(&self) -> Option<usize> {
        self.cols.first().map(|&col| col as usize)
    }

    pub fn leading_value(&self) -> Option<u8> {
        self.vals.first().copied()
    }

    pub fn coefficient(&self, col: usize) -> u8 {
        let Ok(col) = u32::try_from(col) else {
            return 0;
        };
        self.cols
            .binary_search(&col)
            .map(|index| self.vals[index])
            .unwrap_or(0)
    }

    pub fn to_pairs(&self) -> Vec<(usize, u8)> {
        self.cols
            .iter()
            .copied()
            .zip(self.vals.iter().copied())
            .map(|(col, value)| (col as usize, value))
            .collect()
    }
}

pub fn packed_sparse_rref<const P: u8>(
    num_cols: usize,
    rows: &[PackedSparseRow],
) -> PackedSparseRrefResult {
    packed_sparse_rref_from_rows::<P, _>(num_cols, rows.iter().cloned())
}

pub fn packed_sparse_rref_from_rows<const P: u8, I>(
    num_cols: usize,
    rows: I,
) -> PackedSparseRrefResult
where
    I: IntoIterator<Item = PackedSparseRow>,
{
    assert!(P > 1, "field modulus must be greater than 1");
    let inverses = inverse_table::<P>();
    let mut pivot_rows: BTreeMap<usize, PackedSparseRow> = BTreeMap::new();

    for mut row in rows {
        canonicalize_row::<P>(num_cols, &mut row);
        reduce_by_existing_pivots::<P>(&mut row, &pivot_rows);
        if row.is_empty() {
            continue;
        }

        let pivot_col = row.leading_column().expect("nonempty row has a pivot");
        let pivot_value = row.leading_value().expect("nonempty row has a pivot");
        row = scale_row::<P>(&row, inverses[pivot_value as usize]);

        let existing_pivots = pivot_rows.keys().copied().collect::<Vec<_>>();
        for existing_pivot in existing_pivots {
            let factor = pivot_rows[&existing_pivot].coefficient(pivot_col);
            if factor == 0 {
                continue;
            }
            let updated = subtract_scaled_row::<P>(&pivot_rows[&existing_pivot], factor, &row);
            pivot_rows.insert(existing_pivot, updated);
        }

        pivot_rows.insert(pivot_col, row);
    }

    let pivot_columns = pivot_rows.keys().copied().collect::<Vec<_>>();
    let rows = pivot_columns
        .iter()
        .map(|pivot_col| pivot_rows[pivot_col].clone())
        .collect::<Vec<_>>();
    let rank = pivot_columns.len();

    PackedSparseRrefResult {
        num_cols,
        rows,
        pivot_columns,
        rank,
    }
}

pub fn packed_sparse_rank<const P: u8>(num_cols: usize, rows: &[PackedSparseRow]) -> usize {
    packed_sparse_rref::<P>(num_cols, rows).rank
}

/// Return a packed sparse basis for the right kernel of the matrix with these rows.
pub fn packed_sparse_kernel_basis<const P: u8>(
    num_cols: usize,
    rows: &[PackedSparseRow],
) -> Vec<PackedSparseRow> {
    let reduced = packed_sparse_rref::<P>(num_cols, rows);
    packed_sparse_kernel_basis_from_rref::<P>(&reduced).0
}

/// Return a packed sparse kernel basis together with the corresponding free columns.
pub fn packed_sparse_kernel_basis_with_free_columns<const P: u8>(
    num_cols: usize,
    rows: &[PackedSparseRow],
) -> (Vec<PackedSparseRow>, Vec<usize>) {
    let reduced = packed_sparse_rref::<P>(num_cols, rows);
    packed_sparse_kernel_basis_from_rref::<P>(&reduced)
}

pub fn packed_sparse_kernel_basis_with_free_columns_from_rows<const P: u8, I>(
    num_cols: usize,
    rows: I,
) -> (Vec<PackedSparseRow>, Vec<usize>)
where
    I: IntoIterator<Item = PackedSparseRow>,
{
    let reduced = packed_sparse_rref_from_rows::<P, _>(num_cols, rows);
    packed_sparse_kernel_basis_from_rref::<P>(&reduced)
}

pub fn packed_to_dense<const P: u8>(num_cols: usize, row: &PackedSparseRow) -> Vec<u8> {
    let mut dense = vec![0; num_cols];
    for (col, value) in row.to_pairs() {
        dense[col] = value % P;
    }
    dense
}

pub fn inverse_table<const P: u8>() -> [u8; 256] {
    assert!(P > 1, "field modulus must be greater than 1");
    let mut inverses = [0u8; 256];
    for value in 1..P {
        for candidate in 1..P {
            if mul_mod::<P>(value, candidate) == 1 {
                inverses[value as usize] = candidate;
                break;
            }
        }
        assert!(
            inverses[value as usize] != 0,
            "nonzero element is not invertible; modulus is not prime"
        );
    }
    inverses
}

impl PackedSparseQuotientSpace {
    pub fn from_relations<const P: u8>(
        ambient_dimension: usize,
        relations: &[PackedSparseRow],
    ) -> Self {
        let reduced = packed_sparse_rref::<P>(ambient_dimension, relations);
        let free_columns = complement_columns(ambient_dimension, &reduced.pivot_columns);
        let pivot_row_by_column = reduced
            .pivot_columns
            .iter()
            .copied()
            .enumerate()
            .map(|(row_index, col)| (col, row_index))
            .collect();
        let free_index_by_column = free_columns
            .iter()
            .copied()
            .enumerate()
            .map(|(free_index, col)| (col, free_index))
            .collect();
        Self {
            ambient_dimension,
            relation_rref: reduced.rows,
            pivot_columns: reduced.pivot_columns,
            free_columns,
            pivot_row_by_column,
            free_index_by_column,
        }
    }

    pub fn dimension(&self) -> usize {
        self.free_columns.len()
    }

    pub fn normal_form_sparse<const P: u8>(&self, vector: &PackedSparseRow) -> PackedSparseRow {
        let mut normal = PackedSparseRow::new::<P, _>(self.ambient_dimension, vector.to_pairs());
        while let Some((pivot_position, pivot_row)) = normal
            .cols
            .iter()
            .copied()
            .enumerate()
            .find_map(|(position, col)| {
                self.pivot_row_by_column
                    .get(&(col as usize))
                    .copied()
                    .map(|pivot_row| (position, pivot_row))
            })
        {
            let factor = normal.vals[pivot_position];
            normal = subtract_scaled_row::<P>(&normal, factor, &self.relation_rref[pivot_row]);
        }
        normal
    }

    /// Sparse quotient coordinates indexed by quotient-basis position.
    pub fn quotient_coordinates_sparse<const P: u8>(
        &self,
        vector: &PackedSparseRow,
    ) -> PackedSparseRow {
        let normal = self.normal_form_sparse::<P>(vector);
        let coords = normal
            .to_pairs()
            .into_iter()
            .filter_map(|(ambient_col, coeff)| {
                self.free_index_by_column
                    .get(&ambient_col)
                    .copied()
                    .map(|free_index| (free_index, coeff))
            });
        PackedSparseRow::new::<P, _>(self.dimension(), coords)
    }
}

fn packed_sparse_kernel_basis_from_rref<const P: u8>(
    reduced: &PackedSparseRrefResult,
) -> (Vec<PackedSparseRow>, Vec<usize>) {
    let free_columns = complement_columns(reduced.num_cols, &reduced.pivot_columns);
    let mut basis = Vec::with_capacity(free_columns.len());

    for &free_col in &free_columns {
        let mut entries = Vec::new();
        entries.push((free_col, 1));
        for (pivot_row, &pivot_col) in reduced.pivot_columns.iter().enumerate() {
            let coeff = reduced.rows[pivot_row].coefficient(free_col);
            if coeff != 0 {
                entries.push((pivot_col, neg_mod::<P>(coeff)));
            }
        }
        basis.push(PackedSparseRow::new::<P, _>(reduced.num_cols, entries));
    }

    (basis, free_columns)
}

fn canonicalize_row<const P: u8>(num_cols: usize, row: &mut PackedSparseRow) {
    if row.cols.is_empty() {
        return;
    }
    let canonical = PackedSparseRow::new::<P, _>(num_cols, row.to_pairs());
    *row = canonical;
}

fn reduce_by_existing_pivots<const P: u8>(
    row: &mut PackedSparseRow,
    pivot_rows: &BTreeMap<usize, PackedSparseRow>,
) {
    loop {
        let Some(leading_col) = row.leading_column() else {
            break;
        };
        if let Some(pivot_row) = pivot_rows.get(&leading_col) {
            let factor = row.leading_value().expect("nonempty row has a pivot");
            *row = subtract_scaled_row::<P>(row, factor, pivot_row);
            continue;
        }

        let pivot_to_reduce = row
            .cols
            .iter()
            .copied()
            .zip(row.vals.iter().copied())
            .skip(1)
            .find_map(|(col, value)| {
                pivot_rows
                    .contains_key(&(col as usize))
                    .then_some((value, col as usize))
            });
        let Some((factor, pivot_col)) = pivot_to_reduce else {
            break;
        };
        let pivot_row = &pivot_rows[&pivot_col];
        *row = subtract_scaled_row::<P>(row, factor, pivot_row);
    }
}

fn scale_row<const P: u8>(row: &PackedSparseRow, scale: u8) -> PackedSparseRow {
    let scale = scale % P;
    if scale == 0 {
        return PackedSparseRow {
            cols: Vec::new(),
            vals: Vec::new(),
        };
    }
    if scale == 1 {
        return row.clone();
    }
    PackedSparseRow {
        cols: row.cols.clone(),
        vals: row
            .vals
            .iter()
            .copied()
            .map(|value| mul_mod::<P>(value, scale))
            .filter(|&value| value != 0)
            .collect(),
    }
}

fn subtract_scaled_row<const P: u8>(
    row: &PackedSparseRow,
    scale: u8,
    pivot_row: &PackedSparseRow,
) -> PackedSparseRow {
    let scale = scale % P;
    if scale == 0 {
        return row.clone();
    }

    let mut cols = Vec::with_capacity(row.len() + pivot_row.len());
    let mut vals = Vec::with_capacity(row.len() + pivot_row.len());
    let mut i = 0usize;
    let mut j = 0usize;

    while i < row.len() || j < pivot_row.len() {
        match (row.cols.get(i), pivot_row.cols.get(j)) {
            (Some(&col_a), Some(&col_b)) if col_a == col_b => {
                let value = sub_mod::<P>(row.vals[i], mul_mod::<P>(scale, pivot_row.vals[j]));
                if value != 0 {
                    cols.push(col_a);
                    vals.push(value);
                }
                i += 1;
                j += 1;
            }
            (Some(&col_a), Some(&col_b)) if col_a < col_b => {
                cols.push(col_a);
                vals.push(row.vals[i]);
                i += 1;
            }
            (Some(_), Some(&col_b)) => {
                let value = neg_mod::<P>(mul_mod::<P>(scale, pivot_row.vals[j]));
                if value != 0 {
                    cols.push(col_b);
                    vals.push(value);
                }
                j += 1;
            }
            (Some(&col_a), None) => {
                cols.push(col_a);
                vals.push(row.vals[i]);
                i += 1;
            }
            (None, Some(&col_b)) => {
                let value = neg_mod::<P>(mul_mod::<P>(scale, pivot_row.vals[j]));
                if value != 0 {
                    cols.push(col_b);
                    vals.push(value);
                }
                j += 1;
            }
            (None, None) => break,
        }
    }

    PackedSparseRow { cols, vals }
}

fn add_mod<const P: u8>(a: u8, b: u8) -> u8 {
    let sum = a as u16 + b as u16;
    if sum >= P as u16 {
        (sum - P as u16) as u8
    } else {
        sum as u8
    }
}

fn sub_mod<const P: u8>(a: u8, b: u8) -> u8 {
    if a >= b {
        a - b
    } else {
        P - (b - a)
    }
}

fn neg_mod<const P: u8>(value: u8) -> u8 {
    if value == 0 {
        0
    } else {
        P - value
    }
}

fn mul_mod<const P: u8>(a: u8, b: u8) -> u8 {
    ((a as u16 * b as u16) % P as u16) as u8
}

fn complement_columns(num_cols: usize, pivot_columns: &[usize]) -> Vec<usize> {
    let mut is_pivot = vec![false; num_cols];
    for &col in pivot_columns {
        if col < num_cols {
            is_pivot[col] = true;
        }
    }
    (0..num_cols).filter(|&col| !is_pivot[col]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse_linear_algebra::{
        sparse_kernel_basis_with_free_columns, sparse_rref, sparse_to_dense, sparse_vector,
        SparseQuotientSpace,
    };
    use crate::{PrimeField, Ring};

    type F251 = PrimeField<251>;

    fn generic_row<const P: u8>(num_cols: usize, row: &PackedSparseRow) -> Vec<(usize, F251)> {
        assert_eq!(P, 251);
        sparse_vector(
            num_cols,
            row.to_pairs()
                .into_iter()
                .map(|(col, value)| (col, F251::from_i64(value as i64))),
        )
    }

    fn generic_rows<const P: u8>(
        num_cols: usize,
        rows: &[PackedSparseRow],
    ) -> Vec<Vec<(usize, F251)>> {
        rows.iter()
            .map(|row| generic_row::<P>(num_cols, row))
            .collect()
    }

    #[test]
    fn test_inverse_table() {
        let inverses = inverse_table::<251>();
        assert_eq!(inverses[1], 1);
        assert_eq!(mul_mod::<251>(37, inverses[37]), 1);
        assert_eq!(mul_mod::<251>(250, inverses[250]), 1);
    }

    #[test]
    fn test_packed_sparse_row_canonicalizes_duplicates() {
        let row = PackedSparseRow::new::<251, _>(5, vec![(2, 10), (1, 7), (2, 241)]);

        assert_eq!(row.to_pairs(), vec![(1, 7)]);
    }

    #[test]
    fn test_packed_sparse_row_rejects_columns_that_cannot_be_stored() {
        let Some(too_many_columns) = (u32::MAX as usize).checked_add(2) else {
            return;
        };

        let result = std::panic::catch_unwind(|| {
            PackedSparseRow::new::<251, _>(too_many_columns, Vec::<(usize, u8)>::new());
        });

        assert!(result.is_err());
    }

    #[test]
    fn test_packed_sparse_row_large_lookup_is_zero() {
        let Some(too_large_column) = (u32::MAX as usize).checked_add(1) else {
            return;
        };
        let row = PackedSparseRow::new::<251, _>(3, vec![(1, 7)]);

        assert_eq!(row.coefficient(too_large_column), 0);
    }

    #[test]
    fn test_packed_sparse_rref_matches_generic_prime_field() {
        let rows = vec![
            PackedSparseRow::new::<251, _>(4, vec![(0, 1), (1, 2), (3, 1)]),
            PackedSparseRow::new::<251, _>(4, vec![(0, 2), (1, 4), (3, 2)]),
            PackedSparseRow::new::<251, _>(4, vec![(1, 1), (2, 1)]),
        ];

        let packed = packed_sparse_rref::<251>(4, &rows);
        let generic = sparse_rref(4, &generic_rows::<251>(4, &rows));
        let packed_dense = packed
            .rows
            .iter()
            .map(|row| packed_to_dense::<251>(4, row))
            .collect::<Vec<_>>();
        let generic_dense = generic
            .rows
            .iter()
            .map(|row| {
                sparse_to_dense(4, row)
                    .into_iter()
                    .map(PrimeField::value)
                    .map(|value| value as u8)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        assert_eq!(packed.pivot_columns, generic.pivot_columns);
        assert_eq!(packed_dense, generic_dense);
    }

    #[test]
    fn test_packed_sparse_rref_reduces_later_pivots() {
        let rows = vec![
            PackedSparseRow::new::<251, _>(4, vec![(2, 1), (3, 1)]),
            PackedSparseRow::new::<251, _>(4, vec![(0, 1), (1, 1), (2, 1)]),
        ];

        let packed = packed_sparse_rref::<251>(4, &rows);
        let generic = sparse_rref(4, &generic_rows::<251>(4, &rows));

        assert_eq!(packed.pivot_columns, generic.pivot_columns);
        assert_eq!(packed.rank, generic.rank);
    }

    #[test]
    fn test_packed_sparse_kernel_basis_matches_generic_prime_field() {
        let rows = vec![
            PackedSparseRow::new::<251, _>(4, vec![(0, 1), (1, 1), (3, 1)]),
            PackedSparseRow::new::<251, _>(4, vec![(1, 2), (2, 1)]),
        ];

        let (packed_basis, packed_free_columns) =
            packed_sparse_kernel_basis_with_free_columns::<251>(4, &rows);
        let (generic_basis, generic_free_columns) =
            sparse_kernel_basis_with_free_columns(4, &generic_rows::<251>(4, &rows));
        let packed_dense = packed_basis
            .iter()
            .map(|row| packed_to_dense::<251>(4, row))
            .collect::<Vec<_>>();
        let generic_dense = generic_basis
            .iter()
            .map(|row| {
                sparse_to_dense(4, row)
                    .into_iter()
                    .map(PrimeField::value)
                    .map(|value| value as u8)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        assert_eq!(packed_free_columns, generic_free_columns);
        assert_eq!(packed_dense, generic_dense);
    }

    #[test]
    fn test_packed_sparse_quotient_coordinates_match_generic_prime_field() {
        let relations = vec![
            PackedSparseRow::new::<251, _>(3, vec![(0, 1), (1, 1)]),
            PackedSparseRow::new::<251, _>(3, vec![(1, 1), (2, 1)]),
        ];
        let vector = PackedSparseRow::new::<251, _>(3, vec![(0, 3), (1, 4), (2, 5)]);
        let packed_quotient = PackedSparseQuotientSpace::from_relations::<251>(3, &relations);
        let generic_quotient =
            SparseQuotientSpace::from_relations(3, &generic_rows::<251>(3, &relations));

        let packed_coordinates = packed_quotient.quotient_coordinates_sparse::<251>(&vector);
        let generic_coordinates =
            generic_quotient.quotient_coordinates_sparse(&generic_row::<251>(3, &vector));
        let generic_packed_coordinates = generic_coordinates
            .into_iter()
            .map(|(col, coeff)| (col, coeff.value() as u8));

        assert_eq!(
            packed_coordinates,
            PackedSparseRow::new::<251, _>(packed_quotient.dimension(), generic_packed_coordinates)
        );
    }
}
