//! Schubert polynomials via divided differences from the longest permutation.

use sym_poly_core::Ring;

use crate::multipoly::MultiPoly;
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
    let v = compose(&inv, &w0);
    let word = reduced_word(&v);

    let delta: Vec<u32> = (0..n).map(|i| (n - 1 - i) as u32).collect();
    let top = MultiPoly::x_power(n, delta);
    partial_word(&top, &word)
}

fn longest_permutation(n: usize) -> Vec<usize> {
    (1..=n).rev().collect()
}

fn inverse_permutation(perm: &[usize]) -> Vec<usize> {
    let n = perm.len();
    let mut inv = vec![0usize; n];
    for (i, &v) in perm.iter().enumerate() {
        inv[v - 1] = i + 1;
    }
    inv
}

fn compose(a: &[usize], b: &[usize]) -> Vec<usize> {
    assert_eq!(a.len(), b.len());
    (0..a.len()).map(|i| a[b[i] - 1]).collect()
}

fn reduced_word(perm: &[usize]) -> Vec<usize> {
    let mut pi = perm.to_vec();
    let mut word = Vec::new();
    loop {
        let mut found = false;
        for i in 0..pi.len().saturating_sub(1) {
            if pi[i] > pi[i + 1] {
                word.push(i);
                pi.swap(i, i + 1);
                found = true;
                break;
            }
        }
        if !found {
            break;
        }
    }
    word
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
