//! Stembridge peak quasisymmetric functions.
//!
//! The peak function indexed by a peak set `Lambda` in degree `n` is returned
//! in the fundamental basis:
//!
//! `K_Lambda = 2^(|Lambda| + 1) sum F_alpha`,
//!
//! where the sum is over compositions `alpha` of `n` such that
//! `Lambda ⊆ Des(alpha) △ (Des(alpha) + 1)`.

use std::collections::{BTreeMap, BTreeSet};

use sym_poly_core::{Composition, Ring};

use crate::basis::QSymBasis;
use crate::qsym_function::QSymFunction;

/// Return whether `peak_set` is a valid peak set in degree `degree`.
///
/// A peak set is a subset of `{2, ..., degree - 1}` with no adjacent entries.
pub fn is_peak_set(peak_set: &[u32], degree: u32) -> bool {
    let set: BTreeSet<u32> = peak_set.iter().copied().collect();
    if set.len() != peak_set.len() {
        return false;
    }
    set.iter().all(|&p| 2 <= p && p < degree) && set.iter().all(|&p| !set.contains(&(p + 1)))
}

/// Stembridge peak quasisymmetric function in the fundamental basis.
pub fn peak_quasisymmetric<C: Ring>(peak_set: &[u32], degree: u32) -> QSymFunction<C> {
    assert!(
        is_peak_set(peak_set, degree),
        "invalid peak set {:?} in degree {}",
        peak_set,
        degree
    );

    let peak_set: BTreeSet<u32> = peak_set.iter().copied().collect();
    let coeff = peak_coefficient::<C>(peak_set.len());
    let mut terms = BTreeMap::new();

    for alpha in Composition::integer_compositions(degree) {
        if peak_set.is_subset(&descent_symmetric_difference(&alpha)) {
            terms.insert(alpha, coeff.clone());
        }
    }

    QSymFunction::from_terms(QSymBasis::Fundamental, terms)
}

fn peak_coefficient<C: Ring>(peak_set_size: usize) -> C {
    let exponent = u32::try_from(peak_set_size + 1).expect("peak set size overflow");
    let value = 2_i64
        .checked_pow(exponent)
        .expect("peak coefficient overflow");
    C::from_i64(value)
}

fn descent_symmetric_difference(alpha: &Composition) -> BTreeSet<u32> {
    let descents = alpha.composition_to_descent_set();
    let shifted: BTreeSet<u32> = descents.iter().map(|&d| d + 1).collect();
    descents.symmetric_difference(&shifted).copied().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_peak_set_validation() {
        assert!(is_peak_set(&[], 4));
        assert!(is_peak_set(&[2], 4));
        assert!(is_peak_set(&[3], 4));
        assert!(!is_peak_set(&[1], 4));
        assert!(!is_peak_set(&[4], 4));
        assert!(!is_peak_set(&[2, 3], 5));
        assert!(!is_peak_set(&[2, 2], 5));
    }

    #[test]
    fn test_peak_quasisymmetric_degree_three_empty() {
        let k_empty: QSymFunction<i64> = peak_quasisymmetric(&[], 3);

        assert_eq!(k_empty.coefficient(&Composition::new(vec![1, 1, 1])), 2);
        assert_eq!(k_empty.coefficient(&Composition::new(vec![1, 2])), 2);
        assert_eq!(k_empty.coefficient(&Composition::new(vec![2, 1])), 2);
        assert_eq!(k_empty.coefficient(&Composition::new(vec![3])), 2);
        assert_eq!(k_empty.terms().len(), 4);
    }

    #[test]
    fn test_peak_quasisymmetric_degree_three_peak_two() {
        let k_two: QSymFunction<i64> = peak_quasisymmetric(&[2], 3);

        assert_eq!(k_two.coefficient(&Composition::new(vec![1, 2])), 4);
        assert_eq!(k_two.coefficient(&Composition::new(vec![2, 1])), 4);
        assert_eq!(k_two.terms().len(), 2);
    }
}
