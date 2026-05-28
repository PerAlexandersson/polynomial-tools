//! Sparse exact linear algebra for modular representation computations.
//!
//! Rows are stored as sorted `(column, coefficient)` pairs with no zero
//! coefficients. The algorithms are generic over [`Field`], but the intended
//! hot path is `PrimeField<P>`: sparse GKM constraints are too large to store
//! as dense rational matrices.

use std::collections::BTreeMap;

use crate::{Field, Ring};

pub type SparseVector<C> = Vec<(usize, C)>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseRrefResult<C> {
    pub num_cols: usize,
    pub rows: Vec<SparseVector<C>>,
    pub pivot_columns: Vec<usize>,
    pub rank: usize,
}

/// A quotient of a sparse coordinate space by a sparse row-span of relations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseQuotientSpace<C> {
    pub ambient_dimension: usize,
    pub relation_rref: Vec<SparseVector<C>>,
    pub pivot_columns: Vec<usize>,
    pub free_columns: Vec<usize>,
}

pub fn sparse_vector<C: Ring, I>(num_cols: usize, entries: I) -> SparseVector<C>
where
    I: IntoIterator<Item = (usize, C)>,
{
    let mut terms: BTreeMap<usize, C> = BTreeMap::new();
    for (col, coeff) in entries {
        assert!(col < num_cols, "sparse vector column out of range");
        if coeff.is_zero() {
            continue;
        }
        let entry = terms.entry(col).or_insert_with(C::zero);
        *entry = entry.clone() + coeff;
    }
    terms
        .into_iter()
        .filter(|(_, coeff)| !coeff.is_zero())
        .collect()
}

pub fn dense_to_sparse<C: Ring>(vector: &[C]) -> SparseVector<C> {
    vector
        .iter()
        .cloned()
        .enumerate()
        .filter(|(_, coeff)| !coeff.is_zero())
        .collect()
}

pub fn sparse_to_dense<C: Ring>(num_cols: usize, vector: &SparseVector<C>) -> Vec<C> {
    let mut dense = vec![C::zero(); num_cols];
    for (col, coeff) in canonical_sparse_vector(num_cols, vector) {
        dense[col] = coeff;
    }
    dense
}

pub fn sparse_rref<C: Field>(num_cols: usize, rows: &[SparseVector<C>]) -> SparseRrefResult<C> {
    sparse_rref_from_rows(num_cols, rows.iter().cloned())
}

pub fn sparse_rref_from_rows<C: Field, I>(num_cols: usize, rows: I) -> SparseRrefResult<C>
where
    I: IntoIterator<Item = SparseVector<C>>,
{
    let mut pivot_rows: BTreeMap<usize, SparseVector<C>> = BTreeMap::new();

    for row in rows {
        let mut reduced = canonical_sparse_vector(num_cols, &row);
        reduce_by_pivots(&mut reduced, &pivot_rows);
        if reduced.is_empty() {
            continue;
        }

        let pivot_col = reduced[0].0;
        let pivot_coeff = reduced[0].1.clone();
        reduced = scale_row(&reduced, C::one() / pivot_coeff);

        let existing_pivots = pivot_rows.keys().copied().collect::<Vec<_>>();
        for existing_pivot in existing_pivots {
            let factor = coefficient_at(&pivot_rows[&existing_pivot], pivot_col);
            if factor.is_zero() {
                continue;
            }
            let updated = subtract_scaled_row(&pivot_rows[&existing_pivot], factor, &reduced);
            pivot_rows.insert(existing_pivot, updated);
        }

        pivot_rows.insert(pivot_col, reduced);
    }

    let pivot_columns = pivot_rows.keys().copied().collect::<Vec<_>>();
    let rows = pivot_columns
        .iter()
        .map(|pivot_col| pivot_rows[pivot_col].clone())
        .collect::<Vec<_>>();
    let rank = pivot_columns.len();

    SparseRrefResult {
        num_cols,
        rows,
        pivot_columns,
        rank,
    }
}

pub fn sparse_rank<C: Field>(num_cols: usize, rows: &[SparseVector<C>]) -> usize {
    sparse_rref(num_cols, rows).rank
}

/// Return a sparse basis for the right kernel of the matrix with these rows.
pub fn sparse_kernel_basis<C: Field>(
    num_cols: usize,
    rows: &[SparseVector<C>],
) -> Vec<SparseVector<C>> {
    let reduced = sparse_rref(num_cols, rows);
    sparse_kernel_basis_from_rref(&reduced).0
}

/// Return a sparse kernel basis together with the corresponding free columns.
pub fn sparse_kernel_basis_with_free_columns<C: Field>(
    num_cols: usize,
    rows: &[SparseVector<C>],
) -> (Vec<SparseVector<C>>, Vec<usize>) {
    let reduced = sparse_rref(num_cols, rows);
    sparse_kernel_basis_from_rref(&reduced)
}

pub fn sparse_kernel_basis_with_free_columns_from_rows<C: Field, I>(
    num_cols: usize,
    rows: I,
) -> (Vec<SparseVector<C>>, Vec<usize>)
where
    I: IntoIterator<Item = SparseVector<C>>,
{
    let reduced = sparse_rref_from_rows(num_cols, rows);
    sparse_kernel_basis_from_rref(&reduced)
}

impl<C: Field> SparseQuotientSpace<C> {
    pub fn from_relations(ambient_dimension: usize, relations: &[SparseVector<C>]) -> Self {
        let reduced = sparse_rref(ambient_dimension, relations);
        let free_columns = complement_columns(ambient_dimension, &reduced.pivot_columns);
        Self {
            ambient_dimension,
            relation_rref: reduced.rows,
            pivot_columns: reduced.pivot_columns,
            free_columns,
        }
    }

    pub fn dimension(&self) -> usize {
        self.free_columns.len()
    }

    pub fn normal_form_sparse(&self, vector: &SparseVector<C>) -> SparseVector<C> {
        let mut normal = canonical_sparse_vector(self.ambient_dimension, vector);
        for (pivot_row, &pivot_col) in self.pivot_columns.iter().enumerate() {
            let factor = coefficient_at(&normal, pivot_col);
            if factor.is_zero() {
                continue;
            }
            normal = subtract_scaled_row(&normal, factor, &self.relation_rref[pivot_row]);
        }
        normal
    }

    /// Sparse quotient coordinates indexed by quotient-basis position.
    pub fn quotient_coordinates_sparse(&self, vector: &SparseVector<C>) -> SparseVector<C> {
        let normal = self.normal_form_sparse(vector);
        let mut coords = Vec::new();
        for (free_index, &ambient_col) in self.free_columns.iter().enumerate() {
            let coeff = coefficient_at(&normal, ambient_col);
            if !coeff.is_zero() {
                coords.push((free_index, coeff));
            }
        }
        coords
    }

    pub fn quotient_coordinates_dense(&self, vector: &SparseVector<C>) -> Vec<C> {
        sparse_to_dense(self.dimension(), &self.quotient_coordinates_sparse(vector))
    }
}

fn sparse_kernel_basis_from_rref<C: Field>(
    reduced: &SparseRrefResult<C>,
) -> (Vec<SparseVector<C>>, Vec<usize>) {
    let free_columns = complement_columns(reduced.num_cols, &reduced.pivot_columns);
    let mut basis = Vec::with_capacity(free_columns.len());

    for &free_col in &free_columns {
        let mut entries = Vec::new();
        entries.push((free_col, C::one()));
        for (pivot_row, &pivot_col) in reduced.pivot_columns.iter().enumerate() {
            let coeff = coefficient_at(&reduced.rows[pivot_row], free_col);
            if !coeff.is_zero() {
                entries.push((pivot_col, -coeff));
            }
        }
        basis.push(sparse_vector(reduced.num_cols, entries));
    }

    (basis, free_columns)
}

fn canonical_sparse_vector<C: Ring>(num_cols: usize, vector: &SparseVector<C>) -> SparseVector<C> {
    sparse_vector(num_cols, vector.iter().cloned())
}

fn reduce_by_pivots<C: Field>(
    row: &mut SparseVector<C>,
    pivot_rows: &BTreeMap<usize, SparseVector<C>>,
) {
    for (&pivot_col, pivot_row) in pivot_rows {
        let factor = coefficient_at(row, pivot_col);
        if factor.is_zero() {
            continue;
        }
        *row = subtract_scaled_row(row, factor, pivot_row);
        if row.is_empty() {
            break;
        }
    }
}

pub fn sparse_coefficient<C: Ring>(row: &SparseVector<C>, col: usize) -> C {
    row.binary_search_by_key(&col, |(entry_col, _)| *entry_col)
        .map(|index| row[index].1.clone())
        .unwrap_or_else(|_| C::zero())
}

fn coefficient_at<C: Ring>(row: &SparseVector<C>, col: usize) -> C {
    sparse_coefficient(row, col)
}

fn scale_row<C: Ring>(row: &SparseVector<C>, scale: C) -> SparseVector<C> {
    if scale.is_zero() {
        return Vec::new();
    }
    row.iter()
        .map(|&(col, ref coeff)| (col, coeff.clone() * scale.clone()))
        .filter(|(_, coeff)| !coeff.is_zero())
        .collect()
}

fn subtract_scaled_row<C: Ring>(
    row: &SparseVector<C>,
    scale: C,
    pivot_row: &SparseVector<C>,
) -> SparseVector<C> {
    if scale.is_zero() {
        return row.clone();
    }

    let mut result = Vec::with_capacity(row.len() + pivot_row.len());
    let mut i = 0usize;
    let mut j = 0usize;

    while i < row.len() || j < pivot_row.len() {
        match (row.get(i), pivot_row.get(j)) {
            (Some(&(col_a, ref value_a)), Some(&(col_b, ref value_b))) if col_a == col_b => {
                let value = value_a.clone() - scale.clone() * value_b.clone();
                if !value.is_zero() {
                    result.push((col_a, value));
                }
                i += 1;
                j += 1;
            }
            (Some(&(col_a, ref value_a)), Some(&(col_b, _))) if col_a < col_b => {
                result.push((col_a, value_a.clone()));
                i += 1;
            }
            (Some(_), Some(&(col_b, ref value_b))) => {
                let value = -(scale.clone() * value_b.clone());
                if !value.is_zero() {
                    result.push((col_b, value));
                }
                j += 1;
            }
            (Some(&(col_a, ref value_a)), None) => {
                result.push((col_a, value_a.clone()));
                i += 1;
            }
            (None, Some(&(col_b, ref value_b))) => {
                let value = -(scale.clone() * value_b.clone());
                if !value.is_zero() {
                    result.push((col_b, value));
                }
                j += 1;
            }
            (None, None) => break,
        }
    }

    result
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
    use num_rational::Ratio;

    use super::*;
    use crate::linear_algebra::{kernel_basis, rref, QuotientSpace};
    use crate::PrimeField;

    type Q = Ratio<i64>;
    type F101 = PrimeField<101>;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn qv(values: &[i64]) -> Vec<Q> {
        values.iter().map(|&value| q(value)).collect()
    }

    fn qm(rows: &[&[i64]]) -> Vec<Vec<Q>> {
        rows.iter().map(|row| qv(row)).collect()
    }

    #[test]
    fn test_sparse_vector_canonicalizes_duplicates_and_zeros() {
        let vector = sparse_vector::<Q, _>(4, vec![(2, q(3)), (0, q(0)), (2, q(-3)), (1, q(5))]);

        assert_eq!(vector, vec![(1, q(5))]);
    }

    #[test]
    fn test_sparse_rref_matches_dense_rref() {
        let dense = qm(&[&[1, 2, 0, 1], &[2, 4, 0, 2], &[0, 1, 1, 0]]);
        let sparse = dense
            .iter()
            .map(|row| dense_to_sparse(row))
            .collect::<Vec<_>>();

        let dense_reduced = rref(&dense);
        let sparse_reduced = sparse_rref(4, &sparse);
        let sparse_as_dense = sparse_reduced
            .rows
            .iter()
            .map(|row| sparse_to_dense(4, row))
            .collect::<Vec<_>>();
        let dense_nonzero = dense_reduced.matrix[..dense_reduced.rank].to_vec();

        assert_eq!(sparse_reduced.pivot_columns, dense_reduced.pivot_columns);
        assert_eq!(sparse_as_dense, dense_nonzero);
    }

    #[test]
    fn test_sparse_kernel_basis_matches_dense_kernel_basis() {
        let dense = qm(&[&[1, 1, 0, 0], &[0, 1, 1, 0]]);
        let sparse = dense
            .iter()
            .map(|row| dense_to_sparse(row))
            .collect::<Vec<_>>();

        let dense_basis = kernel_basis(&dense);
        let sparse_basis = sparse_kernel_basis(4, &sparse)
            .iter()
            .map(|row| sparse_to_dense(4, row))
            .collect::<Vec<_>>();

        assert_eq!(sparse_basis, dense_basis);
    }

    #[test]
    fn test_sparse_quotient_coordinates_match_dense() {
        let relations = qm(&[&[1, 1, 0], &[0, 1, 1]]);
        let sparse_relations = relations
            .iter()
            .map(|row| dense_to_sparse(row))
            .collect::<Vec<_>>();
        let dense_quotient = QuotientSpace::from_relations(3, &relations);
        let sparse_quotient = SparseQuotientSpace::from_relations(3, &sparse_relations);
        let vector = qv(&[3, 4, 5]);
        let sparse_vector = dense_to_sparse(&vector);

        assert_eq!(
            sparse_quotient.quotient_coordinates_dense(&sparse_vector),
            dense_quotient.quotient_coordinates(&vector)
        );
    }

    #[test]
    fn test_sparse_rref_over_prime_field() {
        let rows = vec![
            sparse_vector(3, vec![(0, F101::from_i64(2)), (1, F101::from_i64(4))]),
            sparse_vector(3, vec![(0, F101::from_i64(1)), (2, F101::from_i64(5))]),
        ];
        let reduced = sparse_rref(3, &rows);
        let dense_rows = reduced
            .rows
            .iter()
            .map(|row| sparse_to_dense(3, row))
            .collect::<Vec<_>>();

        assert_eq!(reduced.pivot_columns, vec![0, 1]);
        assert_eq!(
            dense_rows,
            vec![
                vec![F101::one(), F101::zero(), F101::from_i64(5)],
                vec![F101::zero(), F101::one(), F101::from_i64(48)]
            ]
        );
    }

    #[test]
    fn test_sparse_rref_from_rows_matches_slice_entry_point() {
        let rows = vec![
            sparse_vector(4, vec![(0, q(1)), (2, q(3))]),
            sparse_vector(4, vec![(1, q(2)), (2, q(4)), (3, q(6))]),
        ];

        assert_eq!(
            sparse_rref_from_rows(4, rows.clone()),
            sparse_rref(4, &rows)
        );
    }
}
