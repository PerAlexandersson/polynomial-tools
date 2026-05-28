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

impl PackedSparseRow {
    pub fn new<const P: u8, I>(num_cols: usize, entries: I) -> Self
    where
        I: IntoIterator<Item = (usize, u8)>,
    {
        assert!(P > 1, "field modulus must be greater than 1");
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
        self.cols
            .binary_search(&(col as u32))
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
                    .get(&(col as usize))
                    .map(|pivot_row| (value, pivot_row.clone()))
            });
        let Some((factor, pivot_row)) = pivot_to_reduce else {
            break;
        };
        *row = subtract_scaled_row::<P>(row, factor, &pivot_row);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse_linear_algebra::{sparse_rref, sparse_to_dense, sparse_vector};
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
}
