//! Stable Kohnert limits and extended Schur functions.
//!
//! For a finite diagram `D`, Assaf--Searles define the Kohnert
//! quasisymmetric function as the stable limit of the Kohnert polynomials of
//! `0^m x D`.  In a finite calculation, the coefficient of `M_alpha` is the
//! coefficient of the monomial with positive row-weight exactly `alpha` once
//! the vertical shift is large enough.
//!
//! Lock diagrams are right-justified diagrams.  Their stable Kohnert limits
//! are the extended Schur functions.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use sym_poly_core::{Composition, Ring};

use crate::basis::QSymBasis;
use crate::qsym_function::QSymFunction;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct KohnertCell {
    pub col: usize,
    pub row: usize,
}

pub type KohnertDiagram = BTreeSet<KohnertCell>;

/// The right-justified lock diagram of a composition.
///
/// Row `i + 1` has `alpha[i]` cells, and all rows are right-justified in the
/// rectangle of width `max(alpha)`.
pub fn right_justified_diagram(alpha: &[u32]) -> KohnertDiagram {
    let width = alpha.iter().copied().max().unwrap_or(0) as usize;
    let mut diagram = KohnertDiagram::new();
    if width == 0 {
        return diagram;
    }

    for (idx, &row_len) in alpha.iter().enumerate() {
        let row_len = row_len as usize;
        let row = idx + 1;
        for col in (width - row_len + 1)..=width {
            diagram.insert(KohnertCell { col, row });
        }
    }

    diagram
}

/// The left-justified key diagram of a weak composition.
///
/// Row `i + 1` has `alpha[i]` cells in columns `1, ..., alpha[i]`.
pub fn left_justified_diagram(alpha: &[u32]) -> KohnertDiagram {
    let mut diagram = KohnertDiagram::new();
    for (idx, &row_len) in alpha.iter().enumerate() {
        let row = idx + 1;
        for col in 1..=row_len as usize {
            diagram.insert(KohnertCell { col, row });
        }
    }
    diagram
}

/// Shift a diagram upward by `shift` rows.
pub fn shift_diagram(diagram: &KohnertDiagram, shift: usize) -> KohnertDiagram {
    diagram
        .iter()
        .map(|cell| KohnertCell {
            col: cell.col,
            row: cell.row + shift,
        })
        .collect()
}

/// One-step Kohnert moves, with rows indexed from bottom to top.
pub fn kohnert_moves(diagram: &KohnertDiagram) -> Vec<KohnertDiagram> {
    let rows = diagram.iter().map(|cell| cell.row).collect::<BTreeSet<_>>();
    let mut moves = Vec::new();

    for row in rows {
        let Some(rightmost) = diagram
            .iter()
            .filter(|cell| cell.row == row)
            .max_by_key(|cell| cell.col)
            .copied()
        else {
            continue;
        };
        if row == 1 {
            continue;
        }

        let Some(target_row) = (1..row).rev().find(|target_row| {
            !diagram.contains(&KohnertCell {
                col: rightmost.col,
                row: *target_row,
            })
        }) else {
            continue;
        };

        let mut next = diagram.clone();
        next.remove(&rightmost);
        next.insert(KohnertCell {
            col: rightmost.col,
            row: target_row,
        });
        moves.push(next);
    }

    moves
}

/// Closure of a diagram under Kohnert moves.
pub fn kohnert_diagrams(
    initial: &KohnertDiagram,
    max_diagrams: usize,
) -> Result<Vec<KohnertDiagram>, String> {
    let mut seen = BTreeSet::new();
    let mut queue = VecDeque::new();
    seen.insert(initial.clone());
    queue.push_back(initial.clone());

    while let Some(diagram) = queue.pop_front() {
        for next in kohnert_moves(&diagram) {
            if seen.insert(next.clone()) {
                if seen.len() > max_diagrams {
                    return Err(format!(
                        "Kohnert diagram cap exceeded: more than {max_diagrams} diagrams"
                    ));
                }
                queue.push_back(next);
            }
        }
    }

    Ok(seen.into_iter().collect())
}

/// Stable Kohnert quasisymmetric function in the monomial QSym basis.
///
/// The calculation shifts the input diagram upward until the positive-prefix
/// row-weight counts stabilize, or until `max_shift` is reached.
pub fn kohnert_quasisymmetric_monomial_with_options<C: Ring>(
    diagram: &KohnertDiagram,
    max_shift: usize,
    max_diagrams: usize,
) -> Result<QSymFunction<C>, String> {
    let mut previous: Option<BTreeMap<Composition, i64>> = None;

    for shift in 0..=max_shift {
        let shifted = shift_diagram(diagram, shift);
        let counts = positive_prefix_counts(&shifted, max_diagrams)?;
        if previous.as_ref() == Some(&counts) {
            return Ok(counts_to_qsym(counts));
        }
        previous = Some(counts);
    }

    Err(format!(
        "stable Kohnert limit did not stabilize by shift {max_shift}"
    ))
}

/// Extended Schur function in the monomial QSym basis.
pub fn extended_schur_monomial_with_options<C: Ring>(
    alpha: &[u32],
    max_shift: usize,
    max_diagrams: usize,
) -> Result<QSymFunction<C>, String> {
    let diagram = right_justified_diagram(alpha);
    kohnert_quasisymmetric_monomial_with_options(&diagram, max_shift, max_diagrams)
}

/// Extended Schur function in the monomial QSym basis, with conservative
/// defaults intended for small checked examples.
pub fn extended_schur_monomial<C: Ring>(alpha: &[u32]) -> QSymFunction<C> {
    extended_schur_monomial_with_options(alpha, default_shift_bound(alpha), 100_000)
        .expect("failed to compute extended Schur stable Kohnert limit")
}

/// Stable Kohnert limit of the left-justified key diagram of `alpha`.
pub fn key_stable_monomial_with_options<C: Ring>(
    alpha: &[u32],
    max_shift: usize,
    max_diagrams: usize,
) -> Result<QSymFunction<C>, String> {
    let diagram = left_justified_diagram(alpha);
    kohnert_quasisymmetric_monomial_with_options(&diagram, max_shift, max_diagrams)
}

/// Stable Kohnert limit of the left-justified key diagram of `alpha`, with
/// conservative defaults intended for small checked examples.
pub fn key_stable_monomial<C: Ring>(alpha: &[u32]) -> QSymFunction<C> {
    key_stable_monomial_with_options(alpha, default_shift_bound(alpha), 100_000)
        .expect("failed to compute key stable Kohnert limit")
}

/// Stable Kohnert limit of the left-justified key diagram in the fundamental
/// QSym basis.
pub fn key_stable<C: Ring>(alpha: &[u32]) -> QSymFunction<C> {
    key_stable_monomial::<C>(alpha).to_fundamental_basis()
}

/// Extended Schur function in the fundamental QSym basis.
pub fn extended_schur<C: Ring>(alpha: &[u32]) -> QSymFunction<C> {
    extended_schur_monomial::<C>(alpha).to_fundamental_basis()
}

/// Row-strict extended Schur function, using `RSE_alpha = psi(E_alpha)`.
pub fn row_strict_extended_schur<C: Ring>(alpha: &[u32]) -> QSymFunction<C> {
    extended_schur::<C>(alpha).psi_involution()
}

/// Flipped extended Schur function, using the rho involution and reversed
/// index convention.
pub fn flipped_extended_schur<C: Ring>(alpha: &[u32]) -> QSymFunction<C> {
    let rev_alpha = reverse_alpha(alpha);
    extended_schur::<C>(&rev_alpha).rho_involution()
}

/// Backward extended Schur function, using the omega involution and reversed
/// index convention.
pub fn backward_extended_schur<C: Ring>(alpha: &[u32]) -> QSymFunction<C> {
    let rev_alpha = reverse_alpha(alpha);
    extended_schur::<C>(&rev_alpha).omega_involution()
}

fn default_shift_bound(alpha: &[u32]) -> usize {
    alpha.iter().copied().sum::<u32>() as usize + alpha.len() + 2
}

fn reverse_alpha(alpha: &[u32]) -> Vec<u32> {
    alpha.iter().rev().copied().collect()
}

fn positive_prefix_counts(
    initial: &KohnertDiagram,
    max_diagrams: usize,
) -> Result<BTreeMap<Composition, i64>, String> {
    let mut counts = BTreeMap::new();
    for diagram in kohnert_diagrams(initial, max_diagrams)? {
        if let Some(weight) = positive_prefix_weight(&diagram) {
            let entry = counts.entry(weight).or_insert(0);
            *entry += 1;
        }
    }
    Ok(counts)
}

fn positive_prefix_weight(diagram: &KohnertDiagram) -> Option<Composition> {
    if diagram.is_empty() {
        return Some(Composition::empty());
    }

    let max_row = diagram.iter().map(|cell| cell.row).max().unwrap_or(0);
    let mut weight = vec![0u32; max_row];
    for cell in diagram {
        weight[cell.row - 1] += 1;
    }
    while weight.last() == Some(&0) {
        weight.pop();
    }
    if weight.iter().any(|&part| part == 0) {
        return None;
    }
    Some(Composition::new(weight))
}

fn counts_to_qsym<C: Ring>(counts: BTreeMap<Composition, i64>) -> QSymFunction<C> {
    let terms = counts
        .into_iter()
        .filter(|(_, count)| *count != 0)
        .map(|(alpha, count)| (alpha, C::from_i64(count)))
        .collect();
    QSymFunction::from_terms(QSymBasis::Monomial, terms)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::schur_qsym::qsym_schur;

    fn comp(parts: &[u32]) -> Composition {
        Composition::new(parts.to_vec())
    }

    fn all_coefficients_nonnegative(f: &QSymFunction<i64>) -> bool {
        f.terms().values().all(|&coeff| coeff >= 0)
    }

    #[test]
    fn test_right_justified_diagram() {
        let diagram = right_justified_diagram(&[3, 1]);
        assert!(diagram.contains(&KohnertCell { col: 1, row: 1 }));
        assert!(diagram.contains(&KohnertCell { col: 2, row: 1 }));
        assert!(diagram.contains(&KohnertCell { col: 3, row: 1 }));
        assert!(diagram.contains(&KohnertCell { col: 3, row: 2 }));
        assert_eq!(diagram.len(), 4);
    }

    #[test]
    fn test_left_justified_diagram() {
        let diagram = left_justified_diagram(&[3, 1]);
        assert!(diagram.contains(&KohnertCell { col: 1, row: 1 }));
        assert!(diagram.contains(&KohnertCell { col: 2, row: 1 }));
        assert!(diagram.contains(&KohnertCell { col: 3, row: 1 }));
        assert!(diagram.contains(&KohnertCell { col: 1, row: 2 }));
        assert_eq!(diagram.len(), 4);
    }

    #[test]
    fn test_extended_schur_monomial_small_examples() {
        let e21 = extended_schur_monomial::<i64>(&[2, 1]);
        assert_eq!(e21.coefficient(&comp(&[2, 1])), 1);
        assert_eq!(e21.coefficient(&comp(&[1, 2])), 1);
        assert_eq!(e21.coefficient(&comp(&[1, 1, 1])), 2);
        assert_eq!(e21.terms().len(), 3);

        let e12 = extended_schur_monomial::<i64>(&[1, 2]);
        assert_eq!(e12.coefficient(&comp(&[1, 2])), 1);
        assert_eq!(e12.coefficient(&comp(&[1, 1, 1])), 1);
        assert_eq!(e12.terms().len(), 2);
    }

    #[test]
    fn test_key_stable_partition_shapes_are_schur_small_degrees() {
        // Assaf--Searles recall that Demazure characters stabilize to Schur
        // functions.  Left-justified diagrams are the Kohnert diagrams for
        // Demazure characters.
        for n in 1..=4 {
            for alpha in Composition::integer_compositions(n) {
                let parts = alpha.parts();
                if parts.windows(2).all(|window| window[0] >= window[1]) {
                    let stable_key = key_stable::<i64>(parts);
                    let schur = schur_from_partition(parts);
                    assert_eq!(
                        stable_key.terms(),
                        schur.terms(),
                        "partition shape {parts:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_key_stable_ignores_zero_rows_before_sorting_to_partition() {
        // The stable limit of a Demazure character depends only on the sorted
        // positive weight.  This checks a non-dominant key diagram with an
        // internal zero row.
        let stable_key = key_stable::<i64>(&[0, 2, 1]);
        let schur = schur_from_partition(&[2, 1]);
        assert_eq!(stable_key, schur);
    }

    #[test]
    fn test_assaf_searles_example_e_212_fundamental_expansion() {
        // Assaf--Searles, Definition 6.13, give this as the worked stable
        // lock-polynomial example immediately before Corollary 6.14.
        let e212 = extended_schur::<i64>(&[2, 1, 2]);
        assert_eq!(e212.coefficient(&comp(&[2, 1, 2])), 1);
        assert_eq!(e212.coefficient(&comp(&[1, 2, 2])), 1);
        assert_eq!(e212.coefficient(&comp(&[1, 1, 2, 1])), 1);
        assert_eq!(e212.terms().len(), 3);
    }

    #[test]
    fn test_assaf_searles_adjacent_swap_f_positive_example() {
        // Assaf--Searles, Proposition 6.21:
        // E_(2,2,1) - E_(2,1,2) = F_(2,2,1) + F_(1,2,1,1).
        let diff = extended_schur::<i64>(&[2, 2, 1]) - extended_schur::<i64>(&[2, 1, 2]);
        assert_eq!(diff.basis(), QSymBasis::Fundamental);
        assert_eq!(diff.coefficient(&comp(&[2, 2, 1])), 1);
        assert_eq!(diff.coefficient(&comp(&[1, 2, 1, 1])), 1);
        assert_eq!(diff.terms().len(), 2);
        assert!(all_coefficients_nonnegative(&diff));
    }

    #[test]
    fn test_assaf_searles_reverse_hook_single_fundamental_examples() {
        // Assaf--Searles, Proposition 6.22, says E_alpha is a single
        // fundamental F_alpha exactly for reverse hooks alpha = (1^k, ell).
        let reverse_hook = extended_schur::<i64>(&[1, 1, 3]);
        assert_eq!(reverse_hook.coefficient(&comp(&[1, 1, 3])), 1);
        assert_eq!(reverse_hook.terms().len(), 1);

        let not_reverse_hook = extended_schur::<i64>(&[2, 1, 2]);
        assert_eq!(not_reverse_hook.terms().len(), 3);
        assert_eq!(not_reverse_hook.coefficient(&comp(&[2, 1, 2])), 1);
    }

    #[test]
    fn test_extended_schur_partition_shapes_are_schur_small_degrees() {
        for n in 1..=4 {
            for alpha in Composition::integer_compositions(n) {
                let parts = alpha.parts();
                if parts.windows(2).all(|window| window[0] >= window[1]) {
                    let ext = extended_schur::<i64>(parts);
                    let schur = schur_from_partition(parts);
                    assert_eq!(ext.terms(), schur.terms(), "partition shape {parts:?}");
                }
            }
        }
    }

    #[test]
    fn test_extended_schur_positive_in_monomial_and_fundamental() {
        for n in 1..=4 {
            for alpha in Composition::integer_compositions(n) {
                let mon = extended_schur_monomial::<i64>(alpha.parts());
                assert!(
                    all_coefficients_nonnegative(&mon),
                    "monomial positivity failed for {:?}",
                    alpha.parts()
                );
                let fund = mon.to_fundamental_basis();
                assert!(
                    all_coefficients_nonnegative(&fund),
                    "fundamental positivity failed for {:?}: {fund}",
                    alpha.parts()
                );
            }
        }
    }

    #[test]
    fn test_extended_schur_family_fundamental_positive() {
        for n in 1..=4 {
            for alpha in Composition::integer_compositions(n) {
                let families = [
                    extended_schur::<i64>(alpha.parts()),
                    row_strict_extended_schur::<i64>(alpha.parts()),
                    flipped_extended_schur::<i64>(alpha.parts()),
                    backward_extended_schur::<i64>(alpha.parts()),
                ];
                for function in families {
                    let fund = function.to_fundamental_basis();
                    assert!(
                        all_coefficients_nonnegative(&fund),
                        "fundamental positivity failed for {:?}: {fund}",
                        alpha.parts()
                    );
                }
            }
        }
    }

    #[test]
    fn test_extended_schur_involution_family_relations() {
        let alpha = [3, 1, 2];
        let rev_alpha = [2, 1, 3];
        assert_eq!(
            row_strict_extended_schur::<i64>(&alpha).terms(),
            extended_schur::<i64>(&alpha).psi_involution().terms()
        );
        assert_eq!(
            flipped_extended_schur::<i64>(&alpha).terms(),
            extended_schur::<i64>(&rev_alpha).rho_involution().terms()
        );
        assert_eq!(
            backward_extended_schur::<i64>(&alpha).terms(),
            extended_schur::<i64>(&rev_alpha).omega_involution().terms()
        );
    }

    #[test]
    fn test_daugherty_partition_specializations_under_involutions() {
        // Daugherty records the extended-Schur family as a system of four
        // bases under psi, rho, and omega.  On symmetric functions, rho is
        // the identity while psi and omega restrict to classical omega.
        let lambda = [3, 2];
        let lambda_rev = [2, 3];
        let lambda_conj = conjugate_partition(&lambda);

        let s_lambda = schur_from_partition(&lambda);
        let s_lambda_conj = schur_from_partition(&lambda_conj);

        assert_eq!(row_strict_extended_schur::<i64>(&lambda), s_lambda_conj);
        assert_eq!(flipped_extended_schur::<i64>(&lambda_rev), s_lambda);
        assert_eq!(backward_extended_schur::<i64>(&lambda_rev), s_lambda_conj);
    }

    #[test]
    fn test_stable_kohnert_limit_ignores_added_empty_bottom_rows() {
        // Definition 3.10 stabilizes K_{0^m x D}; shifting D first merely
        // increases the number of initially empty rows.
        let diagram = right_justified_diagram(&[2, 1]);
        let shifted = shift_diagram(&diagram, 2);
        let stable = kohnert_quasisymmetric_monomial_with_options::<i64>(&diagram, 8, 10_000)
            .expect("unshifted stable Kohnert limit");
        let stable_shifted =
            kohnert_quasisymmetric_monomial_with_options::<i64>(&shifted, 8, 10_000)
                .expect("shifted stable Kohnert limit");
        assert_eq!(stable, stable_shifted);
    }

    fn schur_from_partition(lambda: &[u32]) -> QSymFunction<i64> {
        let n = lambda.iter().sum();
        let mut schur = QSymFunction::zero(QSymBasis::Fundamental);
        for beta in Composition::integer_compositions(n) {
            if same_partition(lambda, beta.parts()) {
                schur = schur + qsym_schur::<i64>(beta.parts(), n);
            }
        }
        schur
    }

    fn conjugate_partition(lambda: &[u32]) -> Vec<u32> {
        let width = lambda.iter().copied().max().unwrap_or(0);
        (1..=width)
            .map(|col| lambda.iter().filter(|&&part| part >= col).count() as u32)
            .collect()
    }

    fn same_partition(alpha: &[u32], beta: &[u32]) -> bool {
        let mut alpha_sorted = alpha.to_vec();
        let mut beta_sorted = beta.to_vec();
        alpha_sorted.sort_by(|left, right| right.cmp(left));
        beta_sorted.sort_by(|left, right| right.cmp(left));
        alpha_sorted == beta_sorted
    }
}
