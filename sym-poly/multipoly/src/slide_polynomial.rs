//! Slide and glide polynomials.

use std::collections::{BTreeMap, BTreeSet};

use combinatoric_core::WeakComposition;
use sym_poly_core::{Composition, Ring};

use crate::basis::MultiPolyBasis;
use crate::multipoly::{checked_total_degree, MultiPoly};
use crate::multipoly_function::MultiPolyFunction;

/// Compute the monomial slide polynomial M_α for a weak composition α.
///
/// M_α = Σ x^b where b dominance-above α and flat(b) = flat(α).
/// See Assaf-Searles (2016), Definition 3.5.
pub fn monomial_slide_polynomial<C: Ring>(alpha: &[u32]) -> MultiPoly<C> {
    let n = alpha.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    let total = checked_total_degree(alpha);
    let flat_alpha: Vec<u32> = alpha.iter().copied().filter(|&x| x > 0).collect();
    let mut result = MultiPoly::zero(n);

    for beta in Composition::weak_integer_compositions(total, n) {
        let beta_parts = beta.parts();
        if !prefix_dominates(beta_parts, alpha) {
            continue;
        }

        let flat_beta: Vec<u32> = beta_parts.iter().copied().filter(|&x| x > 0).collect();
        if flat_beta != flat_alpha {
            continue;
        }

        let monomial = MultiPoly::monomial(n, beta_parts.to_vec(), C::one());
        result = result + monomial;
    }

    result
}

/// Compute the fundamental slide polynomial F_α for a weak composition α.
pub fn fundamental_slide_polynomial<C: Ring>(alpha: &[u32]) -> MultiPoly<C> {
    let n = alpha.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    let total = checked_total_degree(alpha);
    let flat_alpha: Vec<u32> = alpha.iter().copied().filter(|&x| x > 0).collect();
    let mut result = MultiPoly::zero(n);

    for beta in Composition::weak_integer_compositions(total, n) {
        let beta_parts = beta.parts();
        if !prefix_dominates(beta_parts, alpha) {
            continue;
        }

        let flat_beta: Vec<u32> = beta_parts.iter().copied().filter(|&x| x > 0).collect();
        if !is_refinement(&flat_beta, &flat_alpha) {
            continue;
        }

        let monomial = MultiPoly::monomial(n, beta_parts.to_vec(), C::one());
        result = result + monomial;
    }

    result
}

/// Compute the glide polynomial indexed by `alpha`.
///
/// The coefficient of a monomial is a polynomial in `beta`: each red entry in
/// a weak komposition contributes one factor of `beta`.  Setting `beta = 0`
/// recovers the fundamental slide polynomial.
pub fn glide_polynomial<C: Ring>(alpha: &[u32], beta: &C) -> MultiPoly<C> {
    let n = alpha.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    let colored_glides = colored_glides(alpha);
    let mut terms: BTreeMap<Vec<u32>, C> = BTreeMap::new();
    for glide in colored_glides {
        let coeff = ring_power(beta, glide.excess());
        if coeff.is_zero() {
            continue;
        }
        let entry = terms.entry(glide.entries).or_insert_with(C::zero);
        *entry = entry.clone() + coeff;
    }

    MultiPoly::from_terms(n, terms)
}

/// Expand a sparse polynomial in the fundamental slide basis.
///
/// This uses the unitriangularity of fundamental slides: every monomial in
/// `F_alpha` prefix-dominates `alpha`, and `x^alpha` appears with coefficient
/// one.  Processing weak compositions by increasing prefix sums therefore
/// gives a support-driven conversion without building a dense transition
/// matrix.
pub fn fundamental_slide_expansion<C: Ring>(poly: &MultiPoly<C>) -> MultiPolyFunction<C> {
    let num_vars = poly.num_vars();
    if poly.is_zero() {
        return MultiPolyFunction::zero(MultiPolyBasis::FundSlide, num_vars);
    }

    let mut degrees = BTreeSet::new();
    for exp in poly.terms().keys() {
        degrees.insert(checked_total_degree(exp));
    }

    let mut result_terms = BTreeMap::new();
    for degree in degrees {
        let mut current: BTreeMap<Vec<u32>, C> = poly
            .terms()
            .iter()
            .filter(|(exp, _)| checked_total_degree(exp) == degree)
            .map(|(exp, coeff)| (exp.clone(), coeff.clone()))
            .collect();

        let mut weak_compositions = WeakComposition::all_weak_compositions(degree, num_vars);
        weak_compositions.sort_by_key(|alpha| prefix_sums(alpha.parts()));

        for alpha in weak_compositions {
            let coeff = current.get(alpha.parts()).cloned().unwrap_or_else(C::zero);
            if coeff.is_zero() {
                continue;
            }
            result_terms.insert(alpha.clone(), coeff.clone());

            let slide = fundamental_slide_polynomial::<C>(alpha.parts());
            for (beta, beta_coeff) in slide.terms() {
                let new_coeff = current.get(beta).cloned().unwrap_or_else(C::zero)
                    - coeff.clone() * beta_coeff.clone();
                if new_coeff.is_zero() {
                    current.remove(beta);
                } else {
                    current.insert(beta.clone(), new_coeff);
                }
            }
        }

        assert!(
            current.is_empty(),
            "fundamental slide expansion left a nonzero remainder"
        );
    }

    MultiPolyFunction::from_terms(MultiPolyBasis::FundSlide, num_vars, result_terms)
}

fn prefix_dominates(beta: &[u32], alpha: &[u32]) -> bool {
    let mut a_sum = 0u32;
    let mut b_sum = 0u32;
    for i in 0..alpha.len() {
        a_sum = a_sum
            .checked_add(alpha[i])
            .expect("composition weight overflow");
        b_sum = b_sum
            .checked_add(beta[i])
            .expect("composition weight overflow");
        if b_sum < a_sum {
            return false;
        }
    }
    true
}

fn prefix_sums(alpha: &[u32]) -> Vec<u32> {
    let mut running = 0u32;
    alpha
        .iter()
        .map(|&part| {
            running = running
                .checked_add(part)
                .expect("composition weight overflow");
            running
        })
        .collect()
}

fn is_refinement(flat_beta: &[u32], flat_alpha: &[u32]) -> bool {
    if flat_alpha.is_empty() {
        return flat_beta.is_empty();
    }

    let mut i = 0usize;
    let mut j = 0usize;
    let mut running = 0u32;

    while i < flat_beta.len() && j < flat_alpha.len() {
        running = running
            .checked_add(flat_beta[i])
            .expect("composition weight overflow");
        if running == flat_alpha[j] {
            j += 1;
            running = 0;
        } else if running > flat_alpha[j] {
            return false;
        }
        i += 1;
    }

    i == flat_beta.len() && j == flat_alpha.len() && running == 0
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct ColoredGlide {
    entries: Vec<u32>,
    red: Vec<bool>,
}

impl ColoredGlide {
    fn excess(&self) -> usize {
        self.red.iter().filter(|&&is_red| is_red).count()
    }
}

fn colored_glides(alpha: &[u32]) -> BTreeSet<ColoredGlide> {
    let n = alpha.len();
    let nonzero_positions = alpha
        .iter()
        .enumerate()
        .filter_map(|(idx, &part)| (part > 0).then_some(idx))
        .collect::<Vec<_>>();

    if nonzero_positions.is_empty() {
        return BTreeSet::from([ColoredGlide {
            entries: vec![0; n],
            red: vec![false; n],
        }]);
    }

    let mut entries = vec![0u32; n];
    let mut red = vec![false; n];
    let mut result = BTreeSet::new();
    collect_colored_glides(
        alpha,
        &nonzero_positions,
        0,
        0,
        &mut entries,
        &mut red,
        &mut result,
    );
    result
}

fn collect_colored_glides(
    alpha: &[u32],
    nonzero_positions: &[usize],
    block_idx: usize,
    start: usize,
    entries: &mut [u32],
    red: &mut [bool],
    result: &mut BTreeSet<ColoredGlide>,
) {
    if block_idx == nonzero_positions.len() {
        result.insert(ColoredGlide {
            entries: entries.to_vec(),
            red: red.to_vec(),
        });
        return;
    }

    let target_pos = nonzero_positions[block_idx];
    let target_part = alpha[target_pos];
    for end_exclusive in start + 1..=target_pos + 1 {
        let len = end_exclusive - start;
        for (segment_entries, segment_red) in colored_glide_segments(len, target_part) {
            entries[start..end_exclusive].copy_from_slice(&segment_entries);
            red[start..end_exclusive].copy_from_slice(&segment_red);

            collect_colored_glides(
                alpha,
                nonzero_positions,
                block_idx + 1,
                end_exclusive,
                entries,
                red,
                result,
            );

            for idx in start..end_exclusive {
                entries[idx] = 0;
                red[idx] = false;
            }
        }
    }
}

fn colored_glide_segments(len: usize, target_part: u32) -> Vec<(Vec<u32>, Vec<bool>)> {
    let mut segments = Vec::new();
    for excess in 0..=len {
        let Some(total) = target_part.checked_add(excess as u32) else {
            continue;
        };
        for weak in Composition::weak_integer_compositions(total, len) {
            let entries = weak.parts().to_vec();
            let positive_positions = entries
                .iter()
                .enumerate()
                .filter_map(|(idx, &entry)| (entry > 0).then_some(idx))
                .collect::<Vec<_>>();
            if positive_positions.is_empty() || excess > positive_positions.len().saturating_sub(1)
            {
                continue;
            }

            let first_positive = positive_positions[0];
            for red_positions in combinations(&positive_positions[1..], excess) {
                let mut red = vec![false; len];
                for pos in red_positions {
                    red[pos] = true;
                }
                red[first_positive] = false;
                segments.push((entries.clone(), red));
            }
        }
    }
    segments
}

fn combinations(items: &[usize], size: usize) -> Vec<Vec<usize>> {
    if size == 0 {
        return vec![Vec::new()];
    }
    if size > items.len() {
        return Vec::new();
    }

    let mut result = Vec::new();
    let mut current = Vec::new();
    collect_combinations(items, size, 0, &mut current, &mut result);
    result
}

fn collect_combinations(
    items: &[usize],
    size: usize,
    start: usize,
    current: &mut Vec<usize>,
    result: &mut Vec<Vec<usize>>,
) {
    if current.len() == size {
        result.push(current.clone());
        return;
    }
    let remaining = size - current.len();
    for idx in start..=items.len() - remaining {
        current.push(items[idx]);
        collect_combinations(items, size, idx + 1, current, result);
        current.pop();
    }
}

fn ring_power<C: Ring>(base: &C, exponent: usize) -> C {
    let mut result = C::one();
    for _ in 0..exponent {
        result = result * base.clone();
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_monomial_slide_zero() {
        let m: MultiPoly<i64> = monomial_slide_polynomial(&[0, 0, 0]);
        assert_eq!(m.coefficient(&[0, 0, 0]), 1);
        assert_eq!(m.terms().len(), 1);
    }

    #[test]
    fn test_monomial_slide_assaf_example() {
        // From Assaf-Searles 2016, Eq. 3.7:
        // M_{(0,2,0,3)} = x1^2*x2^3 + x1^2*x3^3 + x1^2*x4^3 + x2^2*x3^3 + x2^2*x4^3
        let m: MultiPoly<i64> = monomial_slide_polynomial(&[0, 2, 0, 3]);
        assert_eq!(m.terms().len(), 5);
        assert_eq!(m.coefficient(&[2, 3, 0, 0]), 1);
        assert_eq!(m.coefficient(&[2, 0, 3, 0]), 1);
        assert_eq!(m.coefficient(&[2, 0, 0, 3]), 1);
        assert_eq!(m.coefficient(&[0, 2, 3, 0]), 1);
        assert_eq!(m.coefficient(&[0, 2, 0, 3]), 1);
    }

    #[test]
    fn test_monomial_slide_dominant_is_monomial() {
        // For dominant α (weakly decreasing), M_α = x^α
        let m: MultiPoly<i64> = monomial_slide_polynomial(&[3, 2, 1]);
        assert_eq!(m.terms().len(), 1);
        assert_eq!(m.coefficient(&[3, 2, 1]), 1);
    }

    #[test]
    fn test_fundamental_slide_zero() {
        let f: MultiPoly<i64> = fundamental_slide_polynomial(&[0, 0, 0]);
        assert_eq!(f.coefficient(&[0, 0, 0]), 1);
        assert_eq!(f.terms().len(), 1);
    }

    #[test]
    fn test_fundamental_slide_102() {
        // F_(1,0,2) = x1*x3^2 + x1*x2*x3 + x1*x2^2
        let f: MultiPoly<i64> = fundamental_slide_polynomial(&[1, 0, 2]);
        assert_eq!(f.coefficient(&[1, 0, 2]), 1);
        assert_eq!(f.coefficient(&[1, 1, 1]), 1);
        assert_eq!(f.coefficient(&[1, 2, 0]), 1);
        assert_eq!(f.terms().len(), 3);
    }

    #[test]
    fn test_fundamental_slide_11() {
        // F_(1,1) = x1*x2
        let f: MultiPoly<i64> = fundamental_slide_polynomial(&[1, 1]);
        assert_eq!(f.coefficient(&[1, 1]), 1);
        assert_eq!(f.terms().len(), 1);
    }

    #[test]
    fn test_fundamental_slide_expansion_round_trip() {
        let p = fundamental_slide_polynomial::<i64>(&[1, 0, 2])
            + fundamental_slide_polynomial::<i64>(&[0, 3, 0]).scale(&2);
        let expansion = fundamental_slide_expansion(&p);
        assert_eq!(
            expansion.coefficient(&WeakComposition::from_slice(&[1, 0, 2])),
            1
        );
        assert_eq!(
            expansion.coefficient(&WeakComposition::from_slice(&[0, 3, 0])),
            2
        );
        assert_eq!(expansion.terms().len(), 2);
        assert_eq!(expansion.to_multipoly(), p);
    }

    #[test]
    fn test_fundamental_slide_expansion_matches_dense_transition() {
        let p = MultiPoly::monomial(3, vec![1, 0, 2], 1)
            + MultiPoly::monomial(3, vec![0, 3, 0], 2)
            + MultiPoly::monomial(3, vec![2, 1, 0], 3);
        let fast = fundamental_slide_expansion(&p);
        let dense = MultiPolyFunction::from_multipoly(&p).to_fund_slide_basis();
        assert_eq!(fast, dense);
    }

    #[test]
    #[should_panic(expected = "total degree overflow")]
    fn test_slide_rejects_weight_overflow() {
        let _: MultiPoly<i64> = monomial_slide_polynomial(&[u32::MAX, 1]);
    }

    #[test]
    fn test_glide_beta_zero_is_fundamental_slide() {
        for alpha in &[
            vec![0, 0, 0],
            vec![1, 0, 2],
            vec![0, 2, 0, 3],
            vec![2, 1, 0],
        ] {
            let glide: MultiPoly<i64> = glide_polynomial(alpha, &0);
            let slide: MultiPoly<i64> = fundamental_slide_polynomial(alpha);
            assert_eq!(glide, slide, "glide beta=0 mismatch for {alpha:?}");
        }
    }

    #[test]
    fn test_glide_small_beta_terms() {
        let glide: MultiPoly<i64> = glide_polynomial(&[1, 0, 2], &2);

        assert_eq!(glide.coefficient(&[1, 0, 2]), 1);
        assert_eq!(glide.coefficient(&[1, 1, 1]), 1);
        assert_eq!(glide.coefficient(&[1, 2, 0]), 1);
        assert_eq!(glide.coefficient(&[1, 1, 2]), 2);
        assert_eq!(glide.coefficient(&[1, 2, 1]), 2);
        assert_eq!(glide.terms().len(), 5);
    }

    #[test]
    fn test_colored_glides_are_deduplicated() {
        let glides = colored_glides(&[1, 0, 2]);
        let underlying = glides
            .iter()
            .map(|glide| glide.entries.clone())
            .collect::<BTreeSet<_>>();

        assert!(underlying.contains(&vec![1, 2, 0]));
        assert_eq!(
            glides
                .iter()
                .filter(|glide| glide.entries == vec![1, 2, 0])
                .count(),
            1
        );
    }
}
