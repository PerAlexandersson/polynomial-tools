//! Maps between Sym and QSym.
//!
//! - `sym_to_qsym`: Sym -> QSym (inclusion). In monomial bases:
//!   m_λ ↦ Σ_{α: sort(α)=λ} M_α
//!
//! - `qsym_to_sym`: QSym -> Sym (forgetful map). In monomial bases:
//!   M_α ↦ m_{sort(α)}

use std::collections::BTreeMap;

use sym_poly_core::{Composition, Partition, Ring};
use sym_poly_sym::Basis;
use sym_poly_sym::SymmetricFunction;

use crate::basis::QSymBasis;
use crate::qsym_function::QSymFunction;

/// Include a symmetric function into QSym (monomial basis).
///
/// m_λ ↦ Σ_{α: sort(α)=λ} M_α
///
/// This is a ring homomorphism (Sym is a subalgebra of QSym).
pub fn sym_to_qsym<C: Ring>(f: &SymmetricFunction<C>) -> QSymFunction<C> {
    let in_m = f.to_monomial_basis();
    let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();

    for (lambda, coeff) in in_m.terms() {
        // Find all compositions α with sort(α) = λ
        for alpha in compositions_sorting_to(lambda) {
            let entry = result_terms.entry(alpha).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone();
        }
    }

    QSymFunction::from_terms(QSymBasis::Monomial, result_terms)
}

/// Forgetful map QSym -> Sym (monomial basis).
///
/// M_α ↦ m_{sort(α)}
///
/// This is a ring homomorphism.
pub fn qsym_to_sym<C: Ring>(f: &QSymFunction<C>) -> SymmetricFunction<C> {
    let in_m = f.to_monomial_basis();
    let mut result_terms: BTreeMap<Partition, C> = BTreeMap::new();

    for (alpha, coeff) in in_m.terms() {
        let lambda = alpha.to_partition();
        let entry = result_terms.entry(lambda).or_insert_with(C::zero);
        *entry = entry.clone() + coeff.clone();
    }

    SymmetricFunction::from_terms(Basis::Monomial, result_terms)
}

/// All compositions whose parts, when sorted in decreasing order, give λ.
/// These are exactly the permutations of the parts of λ (as a multiset).
fn compositions_sorting_to(lambda: &Partition) -> Vec<Composition> {
    let parts = lambda.parts().to_vec();
    let mut perms = Vec::new();
    multiset_permutations(&parts, &mut Vec::new(), &mut perms);
    perms
}

/// Generate all distinct permutations of a multiset given as a sorted (desc) slice.
fn multiset_permutations(
    remaining: &[u32],
    current: &mut Vec<u32>,
    results: &mut Vec<Composition>,
) {
    if remaining.is_empty() {
        results.push(Composition::new(current.clone()));
        return;
    }
    let mut prev = None;
    for i in 0..remaining.len() {
        if Some(remaining[i]) == prev {
            continue; // skip duplicates
        }
        prev = Some(remaining[i]);
        current.push(remaining[i]);
        let mut rest: Vec<u32> = remaining.to_vec();
        rest.remove(i);
        multiset_permutations(&rest, current, results);
        current.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compositions_sorting_to_21() {
        let lambda = Partition::new(vec![2, 1]);
        let comps = compositions_sorting_to(&lambda);
        // Permutations of {2, 1}: (2,1) and (1,2)
        assert_eq!(comps.len(), 2);
        assert!(comps.contains(&Composition::new(vec![2, 1])));
        assert!(comps.contains(&Composition::new(vec![1, 2])));
    }

    #[test]
    fn test_compositions_sorting_to_211() {
        let lambda = Partition::new(vec![2, 1, 1]);
        let comps = compositions_sorting_to(&lambda);
        // Permutations of {2, 1, 1}: 3!/2! = 3
        assert_eq!(comps.len(), 3);
    }

    #[test]
    fn test_sym_to_qsym_m21() {
        // m_(2,1) -> M_(2,1) + M_(1,2)
        let m: SymmetricFunction<i64> =
            SymmetricFunction::monomial_symmetric(Partition::new(vec![2, 1]));
        let q = sym_to_qsym(&m);
        assert_eq!(q.coefficient(&Composition::new(vec![2, 1])), 1);
        assert_eq!(q.coefficient(&Composition::new(vec![1, 2])), 1);
        assert_eq!(q.terms().len(), 2);
    }

    #[test]
    fn test_sym_to_qsym_m3() {
        // m_(3) -> M_(3)
        let m: SymmetricFunction<i64> =
            SymmetricFunction::monomial_symmetric(Partition::new(vec![3]));
        let q = sym_to_qsym(&m);
        assert_eq!(q.coefficient(&Composition::new(vec![3])), 1);
        assert_eq!(q.terms().len(), 1);
    }

    #[test]
    fn test_qsym_to_sym_m21() {
        // M_(2,1) -> m_(2,1)
        let q: QSymFunction<i64> = QSymFunction::monomial_qsym(Composition::new(vec![2, 1]));
        let s = qsym_to_sym(&q);
        assert_eq!(s.coefficient(&Partition::new(vec![2, 1])), 1);
        assert_eq!(s.terms().len(), 1);
    }

    #[test]
    fn test_roundtrip_sym_qsym_sym() {
        // qsym_to_sym(sym_to_qsym(m_λ)) = |{α : sort(α) = λ}| * m_λ
        // For λ = (2,1): 2 compositions sort to it, so coefficient doubles
        let f: SymmetricFunction<i64> =
            SymmetricFunction::monomial_symmetric(Partition::new(vec![2, 1]));
        let q = sym_to_qsym(&f);
        let back = qsym_to_sym(&q);
        assert_eq!(back.coefficient(&Partition::new(vec![2, 1])), 2);

        // For λ = (3): only 1 composition sorts to it
        let f3: SymmetricFunction<i64> =
            SymmetricFunction::monomial_symmetric(Partition::new(vec![3]));
        let back3 = qsym_to_sym(&sym_to_qsym(&f3));
        assert_eq!(back3.coefficient(&Partition::new(vec![3])), 1);

        // For λ = (1,1,1): only 1 composition sorts to it
        let f111: SymmetricFunction<i64> =
            SymmetricFunction::monomial_symmetric(Partition::new(vec![1, 1, 1]));
        let back111 = qsym_to_sym(&sym_to_qsym(&f111));
        assert_eq!(back111.coefficient(&Partition::new(vec![1, 1, 1])), 1);
    }

    #[test]
    fn test_qsym_to_sym_of_symmetric() {
        // If we start with a symmetric function, include in QSym, and apply
        // qsym_to_sym, each m_λ coefficient gets multiplied by the number
        // of distinct permutations of λ's parts.
        for lambda in Partition::all_of_size(4) {
            let f: SymmetricFunction<i64> = SymmetricFunction::monomial_symmetric(lambda.clone());
            let q = sym_to_qsym(&f);
            let back = qsym_to_sym(&q);
            let num_perms = compositions_sorting_to(&lambda).len() as i64;
            assert_eq!(
                back.coefficient(&lambda),
                num_perms,
                "scaling for m_{}",
                lambda
            );
        }
    }
}
