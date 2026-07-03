//! Hivert's quasisymmetric Hall--Littlewood polynomials.
//!
//! We use the normalization from Hivert and from
//! Loehr--Serrano--Warrington: `G_alpha(t)` specializes to the fundamental
//! quasisymmetric function at `t = 0` and to the monomial quasisymmetric
//! function at `t = 1`.

use std::collections::BTreeMap;

use sym_poly_core::{Composition, UnivariatePolynomial};

use crate::{QSymBasis, QSymFunction};

pub type TPolynomial = UnivariatePolynomial<i64>;

/// Expansion of Hivert's `G_alpha(t)` in the fundamental basis.
pub fn hivert_fundamental_expansion(alpha: &Composition) -> QSymFunction<TPolynomial> {
    let mut terms = BTreeMap::new();
    for (beta, bre) in refinements_with_bre(alpha) {
        let sign = if (beta.num_parts() - alpha.num_parts()) % 2 == 0 {
            1
        } else {
            -1
        };
        let exponent = s_stat(&bre);
        terms.insert(
            beta,
            TPolynomial::monomial(usize::try_from(exponent).expect("exponent too large"), sign),
        );
    }
    QSymFunction::from_terms(QSymBasis::Fundamental, terms)
}

/// Expansion of Hivert's `G_alpha(t)` in the monomial basis.
pub fn hivert_monomial_expansion(alpha: &Composition) -> QSymFunction<TPolynomial> {
    hivert_fundamental_expansion(alpha).to_monomial_basis()
}

/// Coefficients of `F_alpha` in the Hivert `G_beta(t)` basis.
///
/// The return value is a map `beta -> coefficient`, where beta runs over
/// refinements of alpha.
pub fn fundamental_in_hivert_expansion(alpha: &Composition) -> BTreeMap<Composition, TPolynomial> {
    refinements_with_bre(alpha)
        .into_iter()
        .map(|(beta, bre)| {
            let exponent = g_stat(&bre);
            (
                beta,
                TPolynomial::monomial(usize::try_from(exponent).expect("exponent too large"), 1),
            )
        })
        .collect()
}

fn refinements_with_bre(alpha: &Composition) -> Vec<(Composition, Vec<u32>)> {
    let mut out = Vec::new();
    let mut current_parts = Vec::new();
    let mut current_bre = Vec::new();
    refinements_rec(
        alpha.parts(),
        0,
        &mut current_parts,
        &mut current_bre,
        &mut out,
    );
    out
}

fn refinements_rec(
    parts: &[u32],
    index: usize,
    current_parts: &mut Vec<u32>,
    current_bre: &mut Vec<u32>,
    out: &mut Vec<(Composition, Vec<u32>)>,
) {
    if index == parts.len() {
        out.push((Composition::new(current_parts.clone()), current_bre.clone()));
        return;
    }

    for refinement in Composition::integer_compositions(parts[index]) {
        current_bre
            .push(u32::try_from(refinement.num_parts()).expect("refinement length too large"));
        let old_len = current_parts.len();
        current_parts.extend_from_slice(refinement.parts());
        refinements_rec(parts, index + 1, current_parts, current_bre, out);
        current_parts.truncate(old_len);
        current_bre.pop();
    }
}

fn s_stat(bre: &[u32]) -> u32 {
    bre.iter()
        .enumerate()
        .map(|(j, &parts_from_block)| {
            let block_index = u32::try_from(j + 1).expect("block index too large");
            block_index * (parts_from_block - 1)
        })
        .sum()
}

fn g_stat(bre: &[u32]) -> u32 {
    let mut total = 0u32;
    let mut position = 0u32;
    for &block_len in bre {
        for offset in 1..block_len {
            total += position + offset;
        }
        position += block_len;
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    fn composition(parts: &[u32]) -> Composition {
        Composition::new(parts.to_vec())
    }

    fn polynomial(coeffs: &[i64]) -> TPolynomial {
        TPolynomial::new(coeffs.to_vec())
    }

    #[test]
    fn test_hivert_g3_fundamental_expansion() {
        let expansion = hivert_fundamental_expansion(&composition(&[3]));
        assert_eq!(expansion.coefficient(&composition(&[3])), polynomial(&[1]));
        assert_eq!(
            expansion.coefficient(&composition(&[2, 1])),
            polynomial(&[0, -1])
        );
        assert_eq!(
            expansion.coefficient(&composition(&[1, 2])),
            polynomial(&[0, -1])
        );
        assert_eq!(
            expansion.coefficient(&composition(&[1, 1, 1])),
            polynomial(&[0, 0, 1])
        );
        assert_eq!(expansion.terms().len(), 4);
    }

    #[test]
    fn test_hivert_g3_monomial_expansion() {
        let expansion = hivert_monomial_expansion(&composition(&[3]));
        assert_eq!(expansion.coefficient(&composition(&[3])), polynomial(&[1]));
        assert_eq!(
            expansion.coefficient(&composition(&[2, 1])),
            polynomial(&[1, -1])
        );
        assert_eq!(
            expansion.coefficient(&composition(&[1, 2])),
            polynomial(&[1, -1])
        );
        assert_eq!(
            expansion.coefficient(&composition(&[1, 1, 1])),
            polynomial(&[1, -2, 1])
        );
        assert_eq!(expansion.terms().len(), 4);
    }

    #[test]
    fn test_fundamental_in_hivert_expansion_degree_three() {
        let expansion = fundamental_in_hivert_expansion(&composition(&[3]));
        assert_eq!(expansion.get(&composition(&[3])), Some(&polynomial(&[1])));
        assert_eq!(
            expansion.get(&composition(&[2, 1])),
            Some(&polynomial(&[0, 1]))
        );
        assert_eq!(
            expansion.get(&composition(&[1, 2])),
            Some(&polynomial(&[0, 1]))
        );
        assert_eq!(
            expansion.get(&composition(&[1, 1, 1])),
            Some(&polynomial(&[0, 0, 0, 1]))
        );
        assert_eq!(expansion.len(), 4);
    }
}
