//! Dense exact linear algebra for representation computations.
//!
//! This module is deliberately small and generic. It supports the row-reduction
//! and quotient-coordinate operations needed by graded `S_n`-module
//! computations, while leaving the older integer transition-matrix helpers in
//! [`crate::matrix`] untouched.

use crate::{Field, Ring};

pub type Matrix<C> = Vec<Vec<C>>;
pub type Vector<C> = Vec<C>;

/// Reduced row-echelon form data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RrefResult<C> {
    pub matrix: Matrix<C>,
    pub pivot_columns: Vec<usize>,
    pub rank: usize,
}

/// A quotient of a coordinate space by a row-span of relations.
///
/// Relations are row vectors in the ambient coordinate basis. Normal forms are
/// obtained by eliminating pivot coordinates using the RREF of the relation
/// matrix, leaving the free coordinates as quotient coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QuotientSpace<C> {
    pub ambient_dimension: usize,
    pub relation_rref: Matrix<C>,
    pub pivot_columns: Vec<usize>,
    pub free_columns: Vec<usize>,
}

pub fn zero_matrix<C: Ring>(rows: usize, cols: usize) -> Matrix<C> {
    vec![vec![C::zero(); cols]; rows]
}

pub fn identity_matrix<C: Ring>(n: usize) -> Matrix<C> {
    let mut matrix = zero_matrix::<C>(n, n);
    for (i, row) in matrix.iter_mut().enumerate() {
        row[i] = C::one();
    }
    matrix
}

pub fn transpose<C: Ring>(matrix: &[Vec<C>]) -> Matrix<C> {
    let cols = rectangular_num_cols(matrix);
    let mut result = zero_matrix::<C>(cols, matrix.len());
    for (i, row) in matrix.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            result[j][i] = value.clone();
        }
    }
    result
}

pub fn matrix_multiply<C: Ring>(a: &[Vec<C>], b: &[Vec<C>]) -> Matrix<C> {
    let a_cols = rectangular_num_cols(a);
    let b_cols = rectangular_num_cols(b);
    assert_eq!(
        a_cols,
        b.len(),
        "matrix dimensions do not match for multiplication"
    );

    let mut result = zero_matrix::<C>(a.len(), b_cols);
    for (i, row_a) in a.iter().enumerate() {
        for (k, value_a) in row_a.iter().enumerate() {
            if value_a.is_zero() {
                continue;
            }
            for j in 0..b_cols {
                result[i][j] = result[i][j].clone() + value_a.clone() * b[k][j].clone();
            }
        }
    }
    result
}

pub fn matrix_vector_multiply<C: Ring>(matrix: &[Vec<C>], vector: &[C]) -> Vector<C> {
    let cols = rectangular_num_cols(matrix);
    assert_eq!(
        cols,
        vector.len(),
        "matrix and vector dimensions do not match"
    );

    matrix
        .iter()
        .map(|row| {
            row.iter()
                .zip(vector.iter())
                .fold(C::zero(), |acc, (a, x)| acc + a.clone() * x.clone())
        })
        .collect()
}

pub fn matrix_trace<C: Ring>(matrix: &[Vec<C>]) -> C {
    let cols = rectangular_num_cols(matrix);
    assert_eq!(matrix.len(), cols, "trace requires a square matrix");
    (0..matrix.len()).fold(C::zero(), |acc, i| acc + matrix[i][i].clone())
}

/// Compute the reduced row-echelon form.
pub fn rref<C: Field>(matrix: &[Vec<C>]) -> RrefResult<C> {
    let num_cols = rectangular_num_cols(matrix);
    let mut a = matrix.to_vec();
    let mut pivot_columns = Vec::new();
    let mut pivot_row = 0usize;

    for col in 0..num_cols {
        let Some(source_row) = (pivot_row..a.len()).find(|&row| !a[row][col].is_zero()) else {
            continue;
        };

        a.swap(pivot_row, source_row);

        let pivot = a[pivot_row][col].clone();
        for j in col..num_cols {
            a[pivot_row][j] = a[pivot_row][j].clone() / pivot.clone();
        }

        let pivot_snapshot = a[pivot_row].clone();
        for row in 0..a.len() {
            if row == pivot_row || a[row][col].is_zero() {
                continue;
            }
            let factor = a[row][col].clone();
            for j in col..num_cols {
                a[row][j] = a[row][j].clone() - factor.clone() * pivot_snapshot[j].clone();
            }
        }

        pivot_columns.push(col);
        pivot_row += 1;
        if pivot_row == a.len() {
            break;
        }
    }

    strip_zero_rows_to_bottom(&mut a);
    let rank = pivot_columns.len();
    RrefResult {
        matrix: a,
        pivot_columns,
        rank,
    }
}

pub fn rank<C: Field>(matrix: &[Vec<C>]) -> usize {
    rref(matrix).rank
}

/// Return a basis for the right kernel of `matrix`.
///
/// Vectors are returned in ambient coordinate order and span all solutions of
/// `matrix * x = 0`.
pub fn kernel_basis<C: Field>(matrix: &[Vec<C>]) -> Vec<Vector<C>> {
    let num_cols = rectangular_num_cols(matrix);
    let reduced = rref(matrix);
    let free_columns = complement_columns(num_cols, &reduced.pivot_columns);
    let mut basis = Vec::with_capacity(free_columns.len());

    for &free_col in &free_columns {
        let mut vector = vec![C::zero(); num_cols];
        vector[free_col] = C::one();
        for (pivot_row, &pivot_col) in reduced.pivot_columns.iter().enumerate() {
            vector[pivot_col] = -reduced.matrix[pivot_row][free_col].clone();
        }
        basis.push(vector);
    }

    basis
}

/// Solve `a * x = b`.
///
/// If there are multiple solutions, free variables are set to zero.
pub fn solve_linear_system<C: Field>(a: &[Vec<C>], b: &[C]) -> Option<Vector<C>> {
    let num_cols = rectangular_num_cols(a);
    assert_eq!(a.len(), b.len(), "right-hand side has the wrong length");

    let augmented: Matrix<C> = a
        .iter()
        .zip(b.iter())
        .map(|(row, rhs)| {
            let mut augmented_row = row.clone();
            augmented_row.push(rhs.clone());
            augmented_row
        })
        .collect();
    let reduced = rref(&augmented);

    for row in &reduced.matrix {
        let all_zero_left = row[..num_cols].iter().all(Ring::is_zero);
        if all_zero_left && !row[num_cols].is_zero() {
            return None;
        }
    }

    let mut solution = vec![C::zero(); num_cols];
    for (pivot_row, &pivot_col) in reduced.pivot_columns.iter().enumerate() {
        if pivot_col < num_cols {
            solution[pivot_col] = reduced.matrix[pivot_row][num_cols].clone();
        }
    }
    Some(solution)
}

impl<C: Field> QuotientSpace<C> {
    pub fn from_relations(ambient_dimension: usize, relations: &[Vec<C>]) -> Self {
        assert!(
            relations.iter().all(|row| row.len() == ambient_dimension),
            "all relations must have ambient dimension length"
        );
        let reduced = rref(relations);
        let free_columns = complement_columns(ambient_dimension, &reduced.pivot_columns);
        Self {
            ambient_dimension,
            relation_rref: reduced.matrix,
            pivot_columns: reduced.pivot_columns,
            free_columns,
        }
    }

    pub fn dimension(&self) -> usize {
        self.free_columns.len()
    }

    pub fn normal_form(&self, vector: &[C]) -> Vector<C> {
        assert_eq!(
            vector.len(),
            self.ambient_dimension,
            "vector has wrong ambient dimension"
        );
        let mut result = vector.to_vec();
        for (pivot_row, &pivot_col) in self.pivot_columns.iter().enumerate() {
            let pivot_value = result[pivot_col].clone();
            if pivot_value.is_zero() {
                continue;
            }
            for col in pivot_col..self.ambient_dimension {
                result[col] = result[col].clone()
                    - pivot_value.clone() * self.relation_rref[pivot_row][col].clone();
            }
        }
        result
    }

    pub fn quotient_coordinates(&self, vector: &[C]) -> Vector<C> {
        let normal = self.normal_form(vector);
        self.free_columns
            .iter()
            .map(|&col| normal[col].clone())
            .collect()
    }
}

fn rectangular_num_cols<C>(matrix: &[Vec<C>]) -> usize {
    let Some(first_row) = matrix.first() else {
        return 0;
    };
    let cols = first_row.len();
    assert!(
        matrix.iter().all(|row| row.len() == cols),
        "matrix rows have inconsistent lengths"
    );
    cols
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

fn strip_zero_rows_to_bottom<C: Ring>(matrix: &mut Matrix<C>) {
    matrix.sort_by_key(|row| row.iter().all(Ring::is_zero));
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Ratio;

    type Q = Ratio<i64>;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn qv(values: &[i64]) -> Vector<Q> {
        values.iter().map(|&value| q(value)).collect()
    }

    fn qm(rows: &[&[i64]]) -> Matrix<Q> {
        rows.iter().map(|row| qv(row)).collect()
    }

    #[test]
    fn test_identity_transpose_multiply_and_trace() {
        let id = identity_matrix::<Q>(2);
        let a = qm(&[&[1, 2], &[3, 4]]);

        assert_eq!(matrix_multiply(&id, &a), a);
        assert_eq!(matrix_multiply(&a, &id), a);
        assert_eq!(transpose(&a), qm(&[&[1, 3], &[2, 4]]));
        assert_eq!(matrix_trace(&a), q(5));
        assert_eq!(matrix_vector_multiply(&a, &qv(&[1, -1])), qv(&[-1, -1]));
    }

    #[test]
    fn test_rref_rank_one_matrix() {
        let reduced = rref(&qm(&[&[1, 2], &[2, 4]]));

        assert_eq!(reduced.rank, 1);
        assert_eq!(reduced.pivot_columns, vec![0]);
        assert_eq!(reduced.matrix, qm(&[&[1, 2], &[0, 0]]));
    }

    #[test]
    fn test_rref_full_rank_matrix() {
        let reduced = rref(&qm(&[&[1, 2], &[3, 4]]));

        assert_eq!(reduced.rank, 2);
        assert_eq!(reduced.pivot_columns, vec![0, 1]);
        assert_eq!(reduced.matrix, identity_matrix::<Q>(2));
    }

    #[test]
    fn test_kernel_basis_for_single_relation() {
        let kernel = kernel_basis(&qm(&[&[1, 1]]));

        assert_eq!(kernel, vec![qv(&[-1, 1])]);
    }

    #[test]
    fn test_kernel_basis_for_rank_one_rectangular_matrix() {
        let kernel = kernel_basis(&qm(&[&[1, 2, 3], &[2, 4, 6]]));

        assert_eq!(kernel, vec![qv(&[-2, 1, 0]), qv(&[-3, 0, 1])]);
    }

    #[test]
    fn test_solve_linear_system_unique_solution() {
        let a = qm(&[&[1, 2], &[3, 4]]);
        let b = qv(&[5, 11]);

        assert_eq!(solve_linear_system(&a, &b), Some(qv(&[1, 2])));
    }

    #[test]
    fn test_solve_linear_system_underdetermined_sets_free_variables_to_zero() {
        let a = qm(&[&[1, 1]]);
        let b = qv(&[3]);

        assert_eq!(solve_linear_system(&a, &b), Some(qv(&[3, 0])));
    }

    #[test]
    fn test_solve_linear_system_inconsistent() {
        let a = qm(&[&[1, 1], &[2, 2]]);
        let b = qv(&[1, 3]);

        assert_eq!(solve_linear_system(&a, &b), None);
    }

    #[test]
    fn test_quotient_space_mod_one_relation() {
        let quotient = QuotientSpace::from_relations(3, &[qv(&[1, 1, 1])]);

        assert_eq!(quotient.dimension(), 2);
        assert_eq!(quotient.pivot_columns, vec![0]);
        assert_eq!(quotient.free_columns, vec![1, 2]);
        assert_eq!(quotient.normal_form(&qv(&[2, 3, 5])), qv(&[0, 1, 3]));
        assert_eq!(quotient.quotient_coordinates(&qv(&[2, 3, 5])), qv(&[1, 3]));
    }

    #[test]
    fn test_transposition_trace_on_one_dimensional_quotient() {
        let quotient = QuotientSpace::from_relations(2, &[qv(&[1, 1])]);
        let basis_lift = qv(&[0, 1]);
        let swapped = qv(&[1, 0]);
        let image_coordinates = quotient.quotient_coordinates(&swapped);
        let action_matrix = vec![image_coordinates];

        assert_eq!(quotient.quotient_coordinates(&basis_lift), qv(&[1]));
        assert_eq!(action_matrix, qm(&[&[-1]]));
        assert_eq!(matrix_trace(&action_matrix), q(-1));
    }
}
