//! Demazure atoms and their t-deformations via operator recursion.

use sym_poly_core::Ring;

use crate::key_polynomial::{dominant_rearrangement, sorting_reduced_word};
use crate::multipoly::MultiPoly;
use crate::operators::{theta_word, ttheta_word};

/// Compute the Demazure atom A_α for a weak composition α.
pub fn atom_polynomial<C: Ring>(alpha: &[u32]) -> MultiPoly<C> {
    let n = alpha.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    let lambda = dominant_rearrangement(alpha);
    let word = sorting_reduced_word(alpha);
    let init = MultiPoly::x_power(n, lambda);
    theta_word(&init, &word.iter().rev().copied().collect::<Vec<_>>())
}

/// t-deformed operator atom polynomial.
pub fn t_atom_polynomial<C: Ring>(alpha: &[u32], t: &C) -> MultiPoly<C> {
    let n = alpha.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    let lambda = dominant_rearrangement(alpha);
    let word = sorting_reduced_word(alpha);
    let init = MultiPoly::x_power(n, lambda);
    ttheta_word(&init, &word.iter().rev().copied().collect::<Vec<_>>(), t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom_dominant_is_monomial() {
        let atom: MultiPoly<i64> = atom_polynomial(&[3, 2, 1]);
        assert_eq!(atom.coefficient(&[3, 2, 1]), 1);
        assert_eq!(atom.terms().len(), 1);
    }

    #[test]
    fn test_atom_simple_case() {
        // A_(0,2) = θ_0(x_1^2) = x_1 x_2 + x_2^2
        let atom: MultiPoly<i64> = atom_polynomial(&[0, 2]);
        assert_eq!(atom.coefficient(&[1, 1]), 1);
        assert_eq!(atom.coefficient(&[0, 2]), 1);
        assert_eq!(atom.coefficient(&[2, 0]), 0);
    }

    #[test]
    fn test_t_atom_specialization_at_zero() {
        let atom: MultiPoly<i64> = atom_polynomial(&[0, 2]);
        let t_atom: MultiPoly<i64> = t_atom_polynomial(&[0, 2], &0);
        assert_eq!(atom, t_atom);
    }
}
