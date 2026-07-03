//! Operator-side nonsymmetric Macdonald specializations.
//!
//! The full nonsymmetric Macdonald family depends on an affine `q`-shift
//! operator in addition to the Demazure-Lusztig operators. The current crate
//! does not yet have a dedicated coefficient type for symbolic `q` and
//! Laurent-style affine shifts, so this module exposes the operator recursion
//! for the important `q = 0` specialization:
//!
//! - `E_alpha(x; 0, t)`, the nonsymmetric Hall-Littlewood polynomial
//!
//! This is still built from the correct operator layer and serves as the
//! natural stepping stone toward full `q,t` Macdonald support.

use sym_poly_core::Ring;

use crate::key_polynomial::{dominant_rearrangement, sorting_reduced_word};
use crate::multipoly::MultiPoly;
use crate::operators::ttheta_word;

/// Compute the nonsymmetric Hall-Littlewood polynomial `E_alpha(x; 0, t)`.
pub fn nonsymmetric_hall_littlewood<C: Ring>(alpha: &[u32], t: &C) -> MultiPoly<C> {
    let n = alpha.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    let lambda = dominant_rearrangement(alpha);
    let word = sorting_reduced_word(alpha);
    let init = MultiPoly::x_power(n, lambda);
    ttheta_word(&init, &word.iter().rev().copied().collect::<Vec<_>>(), t)
}

/// Compute `E_alpha(x; 0, t)` using Macdonald-family naming.
pub fn nonsymmetric_macdonald_q0<C: Ring>(alpha: &[u32], t: &C) -> MultiPoly<C> {
    nonsymmetric_hall_littlewood(alpha, t)
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;

    use crate::atom_polynomial::{atom_polynomial, t_atom_polynomial};
    use sym_poly_core::Ssaf;

    fn pow_i64(base: i64, exp: u32) -> i64 {
        (0..exp).fold(1, |acc, _| acc * base)
    }

    fn q0_filling_formula(alpha: &[u32], t: i64) -> MultiPoly<i64> {
        let n = alpha.len();
        let basement: Vec<u32> = (1..=n as u32).collect();
        let mut terms = BTreeMap::new();

        for filling in Ssaf::non_attacking_fillings(alpha, &basement) {
            if filling.major_index() != 0 {
                continue;
            }

            let coeff =
                pow_i64(t, filling.coinversions()) * pow_i64(1 - t, filling.horizontal_descents());
            *terms.entry(filling.weight_vector()).or_insert(0) += coeff;
        }

        MultiPoly::from_terms(n, terms)
    }

    #[test]
    fn test_q0_macdonald_matches_t_atom() {
        let q0: MultiPoly<i64> = nonsymmetric_macdonald_q0(&[1, 0, 2], &2);
        let atom: MultiPoly<i64> = t_atom_polynomial(&[1, 0, 2], &2);
        assert_eq!(q0, atom);
    }

    #[test]
    fn test_q0_macdonald_specializes_to_atom_at_t_zero() {
        let q0: MultiPoly<i64> = nonsymmetric_macdonald_q0(&[0, 2], &0);
        let atom: MultiPoly<i64> = atom_polynomial(&[0, 2]);
        assert_eq!(q0, atom);
    }

    #[test]
    fn test_q0_macdonald_matches_filling_formula() {
        let test_cases: Vec<Vec<u32>> = vec![vec![0, 2], vec![2, 1], vec![1, 0, 2], vec![0, 1, 2]];
        let t_values = [0, 2, -1];

        for alpha in &test_cases {
            for &t in &t_values {
                let operator_side: MultiPoly<i64> = nonsymmetric_macdonald_q0(alpha, &t);
                let filling_side = q0_filling_formula(alpha, t);
                assert_eq!(
                    operator_side, filling_side,
                    "operator/filling mismatch for alpha={alpha:?}, t={t}"
                );
            }
        }
    }

    #[test]
    fn test_q0_macdonald_dominant_is_monomial() {
        let q0: MultiPoly<i64> = nonsymmetric_macdonald_q0(&[3, 2, 1], &5);
        assert_eq!(q0.coefficient(&[3, 2, 1]), 1);
        assert_eq!(q0.terms().len(), 1);
    }
}
