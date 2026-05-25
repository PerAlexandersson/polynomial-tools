//! Brenti-style real-rootedness certificates from planar weakly `y`-invariant strip digraphs.
//!
//! This module implements a restricted, proof-friendly subclass of Brenti's
//! digraph model:
//!
//! - vertices are `(x, y)` with `x >= 0` and `0 <= y <= H`,
//! - every edge goes from column `x` to column `x + 1`,
//! - the outgoing edges depend only on the height `y`,
//! - planarity is certified by an explicit monotonicity condition on target
//!   heights.
//!
//! The resulting digraph is automatically weakly `y`-invariant. When the
//! monotonicity condition holds and all weights are nonnegative, Brenti's
//! theorem applies: every row of the path matrix is a Pólya frequency sequence,
//! hence every row polynomial has only real nonpositive zeros.

use crate::{tnn_network::BigRational, Polynomial};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::fmt;

/// A single transition in a Brenti strip digraph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BrentiEdge {
    /// Target height in the next column.
    pub to: usize,
    /// Nonnegative edge weight.
    pub weight: BigRational,
}

/// A restricted weakly `y`-invariant digraph on a finite strip of heights.
///
/// From every vertex `(x, y)`, the outgoing edges are exactly `outgoing[y]`,
/// translated one column to the right.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BrentiStripDigraph {
    pub outgoing: Vec<Vec<BrentiEdge>>,
}

impl BrentiStripDigraph {
    /// Number of tracked heights.
    pub fn num_heights(&self) -> usize {
        self.outgoing.len()
    }

    /// Largest tracked height, if any.
    pub fn max_height(&self) -> Option<usize> {
        self.num_heights().checked_sub(1)
    }

    /// Validate nonnegativity, target range, and the monotone-target planarity condition.
    pub fn validate(&self) -> Result<(), BrentiError> {
        if self.outgoing.is_empty() {
            return Err(BrentiError::EmptyDigraph);
        }

        let max_height = self.num_heights() - 1;
        for (from, edges) in self.outgoing.iter().enumerate() {
            for edge in edges {
                if edge.to > max_height {
                    return Err(BrentiError::TargetOutOfRange {
                        from,
                        to: edge.to,
                        max_height,
                    });
                }
                if edge.weight < BigRational::zero() {
                    return Err(BrentiError::NegativeWeight {
                        from,
                        to: edge.to,
                        weight: edge.weight.clone(),
                    });
                }
            }
        }

        for lower in 0..self.num_heights() {
            let lower_max = self.outgoing[lower].iter().map(|edge| edge.to).max();
            let Some(lower_max) = lower_max else {
                continue;
            };

            for upper in (lower + 1)..self.num_heights() {
                let upper_min = self.outgoing[upper].iter().map(|edge| edge.to).min();
                let Some(upper_min) = upper_min else {
                    continue;
                };

                if lower_max > upper_min {
                    return Err(BrentiError::NonPlanar {
                        lower_height: lower,
                        upper_height: upper,
                        lower_max_target: lower_max,
                        upper_min_target: upper_min,
                    });
                }
            }
        }

        Ok(())
    }

    /// Weighted path counts from `(0,0)` to `(n,k)` for all tracked heights `k`.
    pub fn row_counts(&self, n: usize) -> Result<Vec<BigRational>, BrentiError> {
        self.validate()?;

        let mut current = vec![BigRational::zero(); self.num_heights()];
        current[0] = BigRational::one();

        for _ in 0..n {
            let mut next = vec![BigRational::zero(); self.num_heights()];
            for (from, weight_here) in current.iter().enumerate() {
                if weight_here.is_zero() {
                    continue;
                }
                for edge in &self.outgoing[from] {
                    next[edge.to] += weight_here.clone() * edge.weight.clone();
                }
            }
            current = next;
        }

        Ok(current)
    }

    /// The row polynomial `sum_k M_{n,k} t^k`.
    pub fn row_polynomial(&self, n: usize) -> Result<Polynomial<BigRational>, BrentiError> {
        Ok(Polynomial::new(self.row_counts(n)?))
    }

    /// A prefix of row polynomials for `n = 0, ..., last_row`.
    pub fn row_polynomials_prefix(
        &self,
        last_row: usize,
    ) -> Result<Vec<Polynomial<BigRational>>, BrentiError> {
        (0..=last_row).map(|n| self.row_polynomial(n)).collect()
    }
}

/// A verified Brenti-style certificate for a polynomial sequence prefix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BrentiSequenceCertificate {
    pub digraph: BrentiStripDigraph,
    pub row_indices: Vec<usize>,
    pub matched_rows: Vec<Polynomial<BigRational>>,
}

impl BrentiSequenceCertificate {
    /// Re-run the structural checks and row matching.
    pub fn verify_against_bigint_polynomials(
        &self,
        polys: &[Polynomial<BigInt>],
    ) -> Result<(), BrentiError> {
        self.digraph.validate()?;
        if polys.len() != self.row_indices.len() {
            return Err(BrentiError::RowIndexCountMismatch {
                num_polynomials: polys.len(),
                num_row_indices: self.row_indices.len(),
            });
        }

        for (i, (poly, &row_index)) in polys.iter().zip(self.row_indices.iter()).enumerate() {
            let found = self.digraph.row_polynomial(row_index)?;
            if found != bigint_poly_to_q(poly) {
                return Err(BrentiError::SequenceMismatch {
                    sequence_index: i,
                    row_index,
                    expected: bigint_poly_to_q(poly).coeffs().to_vec(),
                    found: found.coeffs().to_vec(),
                });
            }
        }
        Ok(())
    }
}

/// Errors returned by the restricted Brenti certificate machinery.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BrentiError {
    EmptyDigraph,
    TargetOutOfRange {
        from: usize,
        to: usize,
        max_height: usize,
    },
    NegativeWeight {
        from: usize,
        to: usize,
        weight: BigRational,
    },
    NonPlanar {
        lower_height: usize,
        upper_height: usize,
        lower_max_target: usize,
        upper_min_target: usize,
    },
    RowIndexCountMismatch {
        num_polynomials: usize,
        num_row_indices: usize,
    },
    SequenceMismatch {
        sequence_index: usize,
        row_index: usize,
        expected: Vec<BigRational>,
        found: Vec<BigRational>,
    },
}

impl fmt::Display for BrentiError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyDigraph => write!(f, "digraph must contain at least one height"),
            Self::TargetOutOfRange {
                from,
                to,
                max_height,
            } => write!(
                f,
                "edge from height {} targets height {}, but the maximum supported height is {}",
                from, to, max_height
            ),
            Self::NegativeWeight { from, to, weight } => write!(
                f,
                "edge from height {} to height {} has negative weight {}",
                from, to, weight
            ),
            Self::NonPlanar {
                lower_height,
                upper_height,
                lower_max_target,
                upper_min_target,
            } => write!(
                f,
                "nonplanar strip certificate: height {} reaches as high as {}, but height {} reaches as low as {}",
                lower_height, lower_max_target, upper_height, upper_min_target
            ),
            Self::RowIndexCountMismatch {
                num_polynomials,
                num_row_indices,
            } => write!(
                f,
                "got {} polynomials but {} row indices",
                num_polynomials, num_row_indices
            ),
            Self::SequenceMismatch {
                sequence_index,
                row_index,
                expected,
                found,
            } => write!(
                f,
                "polynomial {} does not match Brenti row {}: expected {:?}, found {:?}",
                sequence_index, row_index, expected, found
            ),
        }
    }
}

impl std::error::Error for BrentiError {}

/// Build a verified Brenti-style certificate for a polynomial sequence prefix.
pub fn build_brenti_sequence_certificate(
    digraph: BrentiStripDigraph,
    row_indices: &[usize],
    polys: &[Polynomial<BigInt>],
) -> Result<BrentiSequenceCertificate, BrentiError> {
    digraph.validate()?;

    if polys.len() != row_indices.len() {
        return Err(BrentiError::RowIndexCountMismatch {
            num_polynomials: polys.len(),
            num_row_indices: row_indices.len(),
        });
    }

    let mut matched_rows = Vec::with_capacity(polys.len());
    for (sequence_index, (poly, &row_index)) in polys.iter().zip(row_indices.iter()).enumerate() {
        let found = digraph.row_polynomial(row_index)?;
        let expected = bigint_poly_to_q(poly);
        if found != expected {
            return Err(BrentiError::SequenceMismatch {
                sequence_index,
                row_index,
                expected: expected.coeffs().to_vec(),
                found: found.coeffs().to_vec(),
            });
        }
        matched_rows.push(found);
    }

    Ok(BrentiSequenceCertificate {
        digraph,
        row_indices: row_indices.to_vec(),
        matched_rows,
    })
}

fn bigint_poly_to_q(poly: &Polynomial<BigInt>) -> Polynomial<BigRational> {
    Polynomial::new(
        poly.coeffs()
            .iter()
            .map(|c| BigRational::from_integer(c.clone()))
            .collect(),
    )
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

    fn pascal_strip(max_height: usize) -> BrentiStripDigraph {
        let mut outgoing = Vec::with_capacity(max_height + 1);
        for y in 0..=max_height {
            let mut edges = vec![BrentiEdge {
                to: y,
                weight: br(1),
            }];
            if y < max_height {
                edges.push(BrentiEdge {
                    to: y + 1,
                    weight: br(1),
                });
            }
            outgoing.push(edges);
        }
        BrentiStripDigraph { outgoing }
    }

    #[test]
    fn test_pascal_strip_row_counts() {
        let digraph = pascal_strip(4);
        let row = digraph.row_counts(4).unwrap();
        assert_eq!(row, vec![br(1), br(4), br(6), br(4), br(1)]);
    }

    #[test]
    fn test_pascal_strip_row_polynomial() {
        let digraph = pascal_strip(4);
        let poly = digraph.row_polynomial(4).unwrap();
        assert_eq!(
            poly,
            Polynomial::new(vec![br(1), br(4), br(6), br(4), br(1)])
        );
    }

    #[test]
    fn test_nonplanar_strip_rejected() {
        let digraph = BrentiStripDigraph {
            outgoing: vec![
                vec![BrentiEdge {
                    to: 1,
                    weight: br(1),
                }],
                vec![BrentiEdge {
                    to: 0,
                    weight: br(1),
                }],
            ],
        };
        match digraph.validate() {
            Err(BrentiError::NonPlanar { .. }) => {}
            other => panic!("expected NonPlanar, got {:?}", other),
        }
    }

    #[test]
    fn test_build_brenti_certificate_pascal_rows() {
        let digraph = pascal_strip(4);
        let polys = vec![
            bi_poly(&[1, 1]),
            bi_poly(&[1, 2, 1]),
            bi_poly(&[1, 3, 3, 1]),
            bi_poly(&[1, 4, 6, 4, 1]),
        ];
        let cert = build_brenti_sequence_certificate(digraph, &[1, 2, 3, 4], &polys).unwrap();
        assert_eq!(cert.row_indices, vec![1, 2, 3, 4]);
        assert_eq!(cert.matched_rows.len(), 4);
    }

    #[test]
    fn test_build_brenti_certificate_mismatch() {
        let digraph = pascal_strip(3);
        let polys = vec![bi_poly(&[1, 2])];
        match build_brenti_sequence_certificate(digraph, &[1], &polys) {
            Err(BrentiError::SequenceMismatch { .. }) => {}
            other => panic!("expected SequenceMismatch, got {:?}", other),
        }
    }
}
