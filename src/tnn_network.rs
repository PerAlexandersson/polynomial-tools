//! Canonical planar-network certificates for lower-unitriangular TNN matrices.
//!
//! This module implements the constructive converse to the LGV lemma in the
//! lower-unitriangular setting: if a finite lower-unitriangular matrix is totally
//! nonnegative, then it is the path matrix of a unique triangular planar network
//! `Gamma_N(lambda)` satisfying the usual zero-propagation condition.
//!
//! The main intended workflow is:
//!
//! 1. start from a monic polynomial sequence,
//! 2. extract its coefficient matrix,
//! 3. check total non-negativity,
//! 4. reconstruct the canonical network,
//! 5. verify the certificate by re-evaluating the path matrix.
//!
//! Note that a network certificate proves TNN of the matrix. Whether that also
//! proves real-rootedness of the original row polynomials depends on the extra
//! structure of the family under study.

use crate::{check_tnn_neville_bigint, Polynomial};
use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Zero};
use std::fmt;

/// Exact rational type used for canonical network weights.
pub type BigRational = Ratio<BigInt>;

/// Canonical triangular planar network `Gamma_N(lambda)`.
///
/// The network has `num_rows = lambda.len() + 1` source rows. Horizontal edges
/// have weight `1`, and `lambda[n][k]` is the weight of the vertical edge
/// `(n + 1, k) -> (n, k)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalPlanarNetwork<T> {
    /// Vertical weights `lambda[n][k]` for `0 <= k <= n < num_rows - 1`.
    pub lambda: Vec<Vec<T>>,
}

impl<T> CanonicalPlanarNetwork<T> {
    /// Number of source rows in the associated triangular network.
    pub fn num_rows(&self) -> usize {
        self.lambda.len() + 1
    }

    /// Access a vertical weight.
    pub fn vertical_weight(&self, n: usize, k: usize) -> Option<&T> {
        self.lambda.get(n).and_then(|row| row.get(k))
    }
}

/// Canonical certificate consisting of a coefficient matrix and its network.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalTnnProof {
    pub coefficient_matrix: Vec<Vec<BigInt>>,
    pub network: CanonicalPlanarNetwork<BigRational>,
}

impl CanonicalTnnProof {
    /// Recompute the path matrix of the stored network.
    pub fn evaluated_path_matrix(&self) -> Vec<Vec<BigRational>> {
        evaluate_path_matrix(&self.network)
    }

    /// Check that the stored network evaluates back to the stored matrix.
    pub fn verify(&self) -> bool {
        verify_path_matrix_certificate(&self.coefficient_matrix, &self.network)
    }
}

/// Errors returned while building or verifying a TNN network certificate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NetworkError {
    EmptyMatrix,
    EmptyPolynomialSequence,
    UnsupportedStartingDegree {
        found: usize,
    },
    UnexpectedPolynomialDegree {
        index: usize,
        expected: usize,
        found: Option<usize>,
    },
    NonMonicPolynomial {
        index: usize,
        leading: Option<BigInt>,
    },
    InvalidRowLength {
        row: usize,
        expected: usize,
        found: usize,
    },
    NegativeEntry {
        row: usize,
        col: usize,
        value: BigInt,
    },
    NonUnitDiagonal {
        row: usize,
        value: BigInt,
    },
    NotTnn {
        reason: String,
    },
    ZeroTailViolation {
        row: usize,
        value: BigInt,
    },
    NegativeReducedEntry {
        row: usize,
        col: usize,
        value: BigRational,
    },
    ReconstructionMismatch,
}

impl fmt::Display for NetworkError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyMatrix => write!(f, "matrix must contain at least one row"),
            Self::EmptyPolynomialSequence => {
                write!(
                    f,
                    "polynomial sequence must contain at least one polynomial"
                )
            }
            Self::UnsupportedStartingDegree { found } => write!(
                f,
                "expected the supplied monic sequence to start at degree 0 or 1, found degree {}",
                found
            ),
            Self::UnexpectedPolynomialDegree {
                index,
                expected,
                found,
            } => write!(
                f,
                "polynomial P_{} should have degree {}, found {:?}",
                index, expected, found
            ),
            Self::NonMonicPolynomial { index, leading } => write!(
                f,
                "polynomial P_{} is not monic; leading coefficient is {:?}",
                index, leading
            ),
            Self::InvalidRowLength {
                row,
                expected,
                found,
            } => write!(
                f,
                "row {} should have length {}, found {}",
                row, expected, found
            ),
            Self::NegativeEntry { row, col, value } => {
                write!(f, "matrix entry ({},{}) = {} is negative", row, col, value)
            }
            Self::NonUnitDiagonal { row, value } => write!(
                f,
                "diagonal entry ({},{}) should be 1, found {}",
                row, row, value
            ),
            Self::NotTnn { reason } => write!(f, "matrix is not totally nonnegative: {}", reason),
            Self::ZeroTailViolation { row, value } => write!(
                f,
                "first-column zero tail violated at row {}: expected 0, found {}",
                row, value
            ),
            Self::NegativeReducedEntry { row, col, value } => write!(
                f,
                "Whitney reduction produced a negative entry at ({},{}): {}",
                row, col, value
            ),
            Self::ReconstructionMismatch => {
                write!(
                    f,
                    "reconstructed network does not evaluate back to the input matrix"
                )
            }
        }
    }
}

impl std::error::Error for NetworkError {}

/// Build the lower-unitriangular coefficient matrix from a monic polynomial sequence.
///
/// The input may start at `P_0 = 1` or at `P_1`. In the latter case this
/// function prepends `P_0 = 1` automatically.
pub fn coefficient_matrix_from_monic_polynomials(
    polys: &[Polynomial<BigInt>],
) -> Result<Vec<Vec<BigInt>>, NetworkError> {
    if polys.is_empty() {
        return Err(NetworkError::EmptyPolynomialSequence);
    }

    let first_degree = polys[0]
        .degree()
        .ok_or(NetworkError::UnexpectedPolynomialDegree {
            index: 0,
            expected: 0,
            found: None,
        })?;

    let offset = match first_degree {
        0 => 0,
        1 => 1,
        d => return Err(NetworkError::UnsupportedStartingDegree { found: d }),
    };

    let mut rows = Vec::with_capacity(polys.len() + offset);
    if offset == 1 {
        rows.push(vec![BigInt::one()]);
    }

    for (supplied_index, poly) in polys.iter().enumerate() {
        let sequence_index = supplied_index + offset;
        let found_degree = poly.degree();
        if found_degree != Some(sequence_index) {
            return Err(NetworkError::UnexpectedPolynomialDegree {
                index: sequence_index,
                expected: sequence_index,
                found: found_degree,
            });
        }

        let leading = poly.leading_coefficient();
        if leading.as_ref() != Some(&BigInt::one()) {
            return Err(NetworkError::NonMonicPolynomial {
                index: sequence_index,
                leading,
            });
        }

        rows.push(poly.coeffs().to_vec());
    }

    Ok(rows)
}

/// Validate the lower-unitriangular / nonnegative shape constraints.
pub fn validate_lower_unitriangular_nonnegative_matrix(
    rows: &[Vec<BigInt>],
) -> Result<(), NetworkError> {
    if rows.is_empty() {
        return Err(NetworkError::EmptyMatrix);
    }

    for (n, row) in rows.iter().enumerate() {
        let expected_len = n + 1;
        if row.len() != expected_len {
            return Err(NetworkError::InvalidRowLength {
                row: n,
                expected: expected_len,
                found: row.len(),
            });
        }

        for (k, value) in row.iter().enumerate() {
            if value < &BigInt::zero() {
                return Err(NetworkError::NegativeEntry {
                    row: n,
                    col: k,
                    value: value.clone(),
                });
            }
        }

        if row[n] != BigInt::one() {
            return Err(NetworkError::NonUnitDiagonal {
                row: n,
                value: row[n].clone(),
            });
        }
    }

    Ok(())
}

/// Reconstruct the canonical network certificate from a lower-unitriangular TNN matrix.
pub fn reconstruct_canonical_tnn_network(
    rows: &[Vec<BigInt>],
) -> Result<CanonicalPlanarNetwork<BigRational>, NetworkError> {
    validate_lower_unitriangular_nonnegative_matrix(rows)?;

    let padded = pad_to_square_matrix(rows);
    check_tnn_neville_bigint(&padded).map_err(|reason| NetworkError::NotTnn { reason })?;

    let first_zero = rows
        .iter()
        .enumerate()
        .skip(1)
        .find(|(_, row)| row[0].is_zero())
        .map(|(idx, _)| idx);
    if let Some(start) = first_zero {
        for (row_idx, row) in rows.iter().enumerate().skip(start) {
            if !row[0].is_zero() {
                return Err(NetworkError::ZeroTailViolation {
                    row: row_idx,
                    value: row[0].clone(),
                });
            }
        }
    }

    let qrows = bigint_triangular_to_q(rows);
    let lambda = reconstruct_lambda_from_q_rows(&qrows)?;
    let network = CanonicalPlanarNetwork { lambda };
    if !verify_path_matrix_certificate(rows, &network) {
        return Err(NetworkError::ReconstructionMismatch);
    }

    Ok(network)
}

/// Build a full canonical TNN certificate from a monic polynomial sequence.
pub fn build_tnn_certificate_from_monic_polynomials(
    polys: &[Polynomial<BigInt>],
) -> Result<CanonicalTnnProof, NetworkError> {
    let coefficient_matrix = coefficient_matrix_from_monic_polynomials(polys)?;
    let network = reconstruct_canonical_tnn_network(&coefficient_matrix)?;
    Ok(CanonicalTnnProof {
        coefficient_matrix,
        network,
    })
}

/// Evaluate the triangular path matrix of `Gamma_N(lambda)`.
pub fn evaluate_path_matrix(
    network: &CanonicalPlanarNetwork<BigRational>,
) -> Vec<Vec<BigRational>> {
    let num_rows = network.num_rows();
    let mut matrix = Vec::with_capacity(num_rows);

    for source in 0..num_rows {
        let mut dp: Vec<Vec<BigRational>> = (0..=source)
            .map(|i| vec![BigRational::zero(); i + 1])
            .collect();
        dp[source][0] = BigRational::one();

        for i in (0..=source).rev() {
            for j in 0..=i {
                let current = dp[i][j].clone();
                if current.is_zero() {
                    continue;
                }

                if j < i {
                    dp[i][j + 1] += current.clone();
                    dp[i - 1][j] += current * network.lambda[i - 1][j].clone();
                }
            }
        }

        matrix.push((0..=source).map(|k| dp[k][k].clone()).collect());
    }

    matrix
}

/// Check whether a network evaluates back to a given integer matrix.
pub fn verify_path_matrix_certificate(
    rows: &[Vec<BigInt>],
    network: &CanonicalPlanarNetwork<BigRational>,
) -> bool {
    evaluate_path_matrix(network) == bigint_triangular_to_q(rows)
}

fn reconstruct_lambda_from_q_rows(
    rows: &[Vec<BigRational>],
) -> Result<Vec<Vec<BigRational>>, NetworkError> {
    let num_rows = rows.len();
    if num_rows <= 1 {
        return Ok(vec![]);
    }

    let first_zero = rows
        .iter()
        .enumerate()
        .skip(1)
        .find(|(_, row)| row[0].is_zero())
        .map(|(idx, _)| idx)
        .unwrap_or(num_rows);

    let mut mu = Vec::with_capacity(num_rows - 1);
    for n in 0..(num_rows - 1) {
        if n < first_zero {
            mu.push(rows[n + 1][0].clone() / rows[n][0].clone());
        } else {
            mu.push(BigRational::zero());
        }
    }

    let mut reduced = Vec::with_capacity(num_rows - 1);
    for n in 0..(num_rows - 1) {
        let mut row = Vec::with_capacity(n + 1);
        for k in 0..=n {
            let entry = rows[n + 1][k + 1].clone() - mu[n].clone() * qcoeff(&rows[n], k + 1);
            if entry < BigRational::zero() {
                return Err(NetworkError::NegativeReducedEntry {
                    row: n,
                    col: k,
                    value: entry,
                });
            }
            row.push(entry);
        }
        reduced.push(row);
    }

    let gamma = reconstruct_lambda_from_q_rows(&reduced)?;
    let mut lambda = Vec::with_capacity(num_rows - 1);
    for n in 0..(num_rows - 1) {
        let mut row = Vec::with_capacity(n + 1);
        row.push(mu[n].clone());
        if n > 0 {
            row.extend(gamma[n - 1].iter().cloned());
        }
        lambda.push(row);
    }

    Ok(lambda)
}

fn pad_to_square_matrix(rows: &[Vec<BigInt>]) -> Vec<Vec<BigInt>> {
    let n = rows.len();
    let mut padded = vec![vec![BigInt::zero(); n]; n];
    for (i, row) in rows.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            padded[i][j] = value.clone();
        }
    }
    padded
}

fn bigint_triangular_to_q(rows: &[Vec<BigInt>]) -> Vec<Vec<BigRational>> {
    rows.iter()
        .map(|row| {
            row.iter()
                .map(|value| BigRational::from_integer(value.clone()))
                .collect()
        })
        .collect()
}

fn qcoeff(row: &[BigRational], index: usize) -> BigRational {
    row.get(index).cloned().unwrap_or_else(BigRational::zero)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bi(v: i64) -> BigInt {
        BigInt::from(v)
    }

    fn br(v: i64) -> BigRational {
        BigRational::from_integer(bi(v))
    }

    fn bi_poly(coeffs: &[i64]) -> Polynomial<BigInt> {
        Polynomial::new(coeffs.iter().map(|&v| bi(v)).collect())
    }

    fn bi_rows(rows: &[&[i64]]) -> Vec<Vec<BigInt>> {
        rows.iter()
            .map(|row| row.iter().map(|&v| bi(v)).collect())
            .collect()
    }

    #[test]
    fn test_coefficient_matrix_from_monic_polynomials_starting_at_zero() {
        let polys = vec![bi_poly(&[1]), bi_poly(&[1, 1]), bi_poly(&[1, 3, 1])];
        let rows = coefficient_matrix_from_monic_polynomials(&polys).unwrap();
        assert_eq!(rows, bi_rows(&[&[1], &[1, 1], &[1, 3, 1]]));
    }

    #[test]
    fn test_coefficient_matrix_from_monic_polynomials_prepends_p0() {
        let polys = vec![
            bi_poly(&[1, 1]),
            bi_poly(&[1, 2, 1]),
            bi_poly(&[1, 3, 3, 1]),
        ];
        let rows = coefficient_matrix_from_monic_polynomials(&polys).unwrap();
        assert_eq!(rows, bi_rows(&[&[1], &[1, 1], &[1, 2, 1], &[1, 3, 3, 1]]));
    }

    #[test]
    fn test_evaluate_pascal_network() {
        let network = CanonicalPlanarNetwork {
            lambda: vec![vec![br(1)], vec![br(1), br(1)], vec![br(1), br(1), br(1)]],
        };
        let matrix = evaluate_path_matrix(&network);
        assert_eq!(
            matrix,
            vec![
                vec![br(1)],
                vec![br(1), br(1)],
                vec![br(1), br(2), br(1)],
                vec![br(1), br(3), br(3), br(1)],
            ]
        );
    }

    #[test]
    fn test_reconstruct_pascal_network() {
        let rows = bi_rows(&[&[1], &[1, 1], &[1, 2, 1], &[1, 3, 3, 1]]);
        let network = reconstruct_canonical_tnn_network(&rows).unwrap();
        assert_eq!(
            network.lambda,
            vec![vec![br(1)], vec![br(1), br(1)], vec![br(1), br(1), br(1)],]
        );
        assert!(verify_path_matrix_certificate(&rows, &network));
    }

    #[test]
    fn test_reconstruct_zero_tail_network() {
        let network = CanonicalPlanarNetwork {
            lambda: vec![vec![br(2)], vec![br(0), br(3)]],
        };
        let rows_q = evaluate_path_matrix(&network);
        let rows: Vec<Vec<BigInt>> = rows_q
            .iter()
            .map(|row| row.iter().map(|q| q.to_integer()).collect())
            .collect();
        assert_eq!(rows, bi_rows(&[&[1], &[2, 1], &[0, 3, 1]]));

        let reconstructed = reconstruct_canonical_tnn_network(&rows).unwrap();
        assert_eq!(reconstructed, network);
    }

    #[test]
    fn test_reconstruct_rejects_non_tnn_matrix() {
        let rows = bi_rows(&[&[1], &[1, 1], &[2, 0, 1]]);
        match reconstruct_canonical_tnn_network(&rows) {
            Err(NetworkError::NotTnn { .. }) => {}
            other => panic!("expected NotTnn, got {:?}", other),
        }
    }

    #[test]
    fn test_build_certificate_from_polynomials() {
        let polys = vec![
            bi_poly(&[1, 1]),
            bi_poly(&[1, 2, 1]),
            bi_poly(&[1, 3, 3, 1]),
        ];
        let proof = build_tnn_certificate_from_monic_polynomials(&polys).unwrap();
        assert!(proof.verify());
        assert_eq!(
            proof.coefficient_matrix,
            bi_rows(&[&[1], &[1, 1], &[1, 2, 1], &[1, 3, 3, 1]])
        );
    }
}
