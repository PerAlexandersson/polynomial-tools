//! Schubert polynomials via divided differences from the longest permutation.

use combinatoric_core::{
    compose_permutations, inverse_permutation, longest_permutation, reduced_word,
};
use sym_poly_core::Ring;

use crate::multipoly::MultiPoly;
use crate::multipoly_function::MultiPolyFunction;
use crate::operators::partial_word;

/// Compute the Schubert polynomial S_w for a permutation in one-line notation.
///
/// The permutation is given using 1-indexed values, e.g. `[3, 1, 2]`.
pub fn schubert_polynomial<C: Ring>(perm: &[usize]) -> MultiPoly<C> {
    let n = perm.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    let w0 = longest_permutation(n);
    let inv = inverse_permutation(perm);
    let v = compose_permutations(&inv, &w0);
    let word = reduced_word(&v);

    let delta: Vec<u32> = (0..n).map(|i| (n - 1 - i) as u32).collect();
    let top = MultiPoly::x_power(n, delta);
    partial_word(&top, &word)
}

/// Express the Schubert polynomial S_w in the monomial basis as a `MultiPolyFunction`.
pub fn schubert_to_monomial<C: Ring>(perm: &[usize]) -> MultiPolyFunction<C> {
    let poly: MultiPoly<C> = schubert_polynomial(perm);
    MultiPolyFunction::from_multipoly(&poly)
}

/// Express the Schubert polynomial S_w in the key basis.
pub fn schubert_to_key<C: Ring>(perm: &[usize]) -> MultiPolyFunction<C> {
    schubert_to_monomial::<C>(perm).to_key_basis()
}

/// Express the Schubert polynomial S_w in the fundamental slide basis.
pub fn schubert_to_fund_slide<C: Ring>(perm: &[usize]) -> MultiPolyFunction<C> {
    schubert_to_monomial::<C>(perm).to_fund_slide_basis()
}

/// Express the Schubert polynomial S_w in the atom basis.
pub fn schubert_to_atom<C: Ring>(perm: &[usize]) -> MultiPolyFunction<C> {
    schubert_to_monomial::<C>(perm).to_atom_basis()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_schubert_identity() {
        let s: MultiPoly<i64> = schubert_polynomial(&[1, 2, 3]);
        assert_eq!(s.coefficient(&[0, 0, 0]), 1);
        assert_eq!(s.terms().len(), 1);
    }

    #[test]
    fn test_schubert_simple_transposition_s1() {
        let s: MultiPoly<i64> = schubert_polynomial(&[2, 1, 3]);
        assert_eq!(s.coefficient(&[1, 0, 0]), 1);
        assert_eq!(s.terms().len(), 1);
    }

    #[test]
    fn test_schubert_simple_transposition_s2() {
        let s: MultiPoly<i64> = schubert_polynomial(&[1, 3, 2]);
        assert_eq!(s.coefficient(&[1, 0, 0]), 1);
        assert_eq!(s.coefficient(&[0, 1, 0]), 1);
        assert_eq!(s.terms().len(), 2);
    }

    #[test]
    fn test_schubert_longest() {
        let s: MultiPoly<i64> = schubert_polynomial(&[3, 2, 1]);
        assert_eq!(s.coefficient(&[2, 1, 0]), 1);
        assert_eq!(s.terms().len(), 1);
    }
}
