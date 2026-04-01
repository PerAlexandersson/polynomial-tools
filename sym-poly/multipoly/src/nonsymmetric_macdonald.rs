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
    use super::*;

    use crate::atom_polynomial::{atom_polynomial, t_atom_polynomial};

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
    fn test_q0_macdonald_dominant_is_monomial() {
        let q0: MultiPoly<i64> = nonsymmetric_macdonald_q0(&[3, 2, 1], &5);
        assert_eq!(q0.coefficient(&[3, 2, 1]), 1);
        assert_eq!(q0.terms().len(), 1);
    }
}
