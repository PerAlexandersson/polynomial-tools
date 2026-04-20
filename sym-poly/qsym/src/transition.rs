//! Basis conversions for quasisymmetric functions.
//!
//! The fundamental relationship is:
//!   F_α = Σ_{β refines α} M_β
//! and by Mobius inversion on the refinement poset:
//!   M_α = Σ_{β refines α} (-1)^{ℓ(β)-ℓ(α)} F_β

use std::collections::BTreeMap;

use sym_poly_core::{Composition, Ring};

use crate::basis::QSymBasis;
use crate::qsym_function::QSymFunction;

/// Convert a QSym function to a different basis.
pub fn convert<C: Ring>(f: &QSymFunction<C>, target: QSymBasis) -> QSymFunction<C> {
    if f.basis() == target {
        return f.clone();
    }
    if f.is_zero() {
        return QSymFunction::zero(target);
    }

    match (f.basis(), target) {
        (QSymBasis::Fundamental, QSymBasis::Monomial) => fundamental_to_monomial(f),
        (QSymBasis::Monomial, QSymBasis::Fundamental) => monomial_to_fundamental(f),

        // PowerSumPsi conversions (via monomial as bridge)
        (QSymBasis::PowerSumPsi, QSymBasis::Monomial) => crate::power_sum::psi_to_monomial(f),
        (QSymBasis::Monomial, QSymBasis::PowerSumPsi) => crate::power_sum::monomial_to_psi(f),
        (QSymBasis::PowerSumPsi, QSymBasis::Fundamental) => {
            let in_m = crate::power_sum::psi_to_monomial(f);
            monomial_to_fundamental(&in_m)
        }
        (QSymBasis::Fundamental, QSymBasis::PowerSumPsi) => {
            let in_m = fundamental_to_monomial(f);
            crate::power_sum::monomial_to_psi(&in_m)
        }

        // PowerSumPhi conversions (via monomial as bridge)
        (QSymBasis::PowerSumPhi, QSymBasis::Monomial) => crate::power_sum::phi_to_monomial(f),
        (QSymBasis::Monomial, QSymBasis::PowerSumPhi) => crate::power_sum::monomial_to_phi(f),
        (QSymBasis::PowerSumPhi, QSymBasis::Fundamental) => {
            let in_m = crate::power_sum::phi_to_monomial(f);
            monomial_to_fundamental(&in_m)
        }
        (QSymBasis::Fundamental, QSymBasis::PowerSumPhi) => {
            let in_m = fundamental_to_monomial(f);
            crate::power_sum::monomial_to_phi(&in_m)
        }

        // Cross power-sum conversions (Ψ ↔ Φ via monomial)
        (QSymBasis::PowerSumPsi, QSymBasis::PowerSumPhi) => {
            let in_m = crate::power_sum::psi_to_monomial(f);
            crate::power_sum::monomial_to_phi(&in_m)
        }
        (QSymBasis::PowerSumPhi, QSymBasis::PowerSumPsi) => {
            let in_m = crate::power_sum::phi_to_monomial(f);
            crate::power_sum::monomial_to_psi(&in_m)
        }

        _ => unreachable!("all QSym basis pairs covered"),
    }
}

/// F_α = Σ_{β refines α} M_β
fn fundamental_to_monomial<C: Ring>(f: &QSymFunction<C>) -> QSymFunction<C> {
    let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();

    for (alpha, coeff) in f.terms() {
        // Sum over all refinements β of α
        for beta in alpha.composition_refinements() {
            let entry = result_terms.entry(beta).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone();
        }
    }

    QSymFunction::from_terms(QSymBasis::Monomial, result_terms)
}

/// M_α = Σ_{β refines α} (-1)^{ℓ(β)-ℓ(α)} F_β
fn monomial_to_fundamental<C: Ring>(f: &QSymFunction<C>) -> QSymFunction<C> {
    let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();

    for (alpha, coeff) in f.terms() {
        let ell_alpha = alpha.num_parts();
        for beta in alpha.composition_refinements() {
            let ell_beta = beta.num_parts();
            let sign_exp = ell_beta - ell_alpha;
            let sign = if sign_exp % 2 == 0 {
                C::one()
            } else {
                C::minus_one()
            };
            let entry = result_terms.entry(beta).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone() * sign;
        }
    }

    QSymFunction::from_terms(QSymBasis::Fundamental, result_terms)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_f_to_m_single_part() {
        // F_(n) = Σ_{α ⊨ n} M_α (all compositions refine (n))
        let f3: QSymFunction<i64> = QSymFunction::fundamental_qsym(Composition::new(vec![3]));
        let in_m = f3.to_monomial_basis();
        // All compositions of 3: (3), (2,1), (1,2), (1,1,1)
        assert_eq!(in_m.coefficient(&Composition::new(vec![3])), 1);
        assert_eq!(in_m.coefficient(&Composition::new(vec![2, 1])), 1);
        assert_eq!(in_m.coefficient(&Composition::new(vec![1, 2])), 1);
        assert_eq!(in_m.coefficient(&Composition::new(vec![1, 1, 1])), 1);
    }

    #[test]
    fn test_f_to_m_21() {
        // F_(2,1): refinements of (2,1) are (2,1) and (1,1,1)
        let f21: QSymFunction<i64> = QSymFunction::fundamental_qsym(Composition::new(vec![2, 1]));
        let in_m = f21.to_monomial_basis();
        assert_eq!(in_m.coefficient(&Composition::new(vec![2, 1])), 1);
        assert_eq!(in_m.coefficient(&Composition::new(vec![1, 1, 1])), 1);
        assert_eq!(in_m.coefficient(&Composition::new(vec![1, 2])), 0);
        assert_eq!(in_m.coefficient(&Composition::new(vec![3])), 0);
    }

    #[test]
    fn test_m_to_f_roundtrip() {
        // M -> F -> M should be identity
        let m: QSymFunction<i64> = QSymFunction::monomial_qsym(Composition::new(vec![2, 1]));
        let in_f = m.to_fundamental_basis();
        let back = in_f.to_monomial_basis();
        assert_eq!(back.coefficient(&Composition::new(vec![2, 1])), 1);
        assert_eq!(back.terms().len(), 1);
    }

    #[test]
    fn test_f_to_m_roundtrip() {
        // F -> M -> F should be identity
        let f: QSymFunction<i64> = QSymFunction::fundamental_qsym(Composition::new(vec![1, 2]));
        let in_m = f.to_monomial_basis();
        let back = in_m.to_fundamental_basis();
        assert_eq!(back.coefficient(&Composition::new(vec![1, 2])), 1);
        assert_eq!(back.terms().len(), 1);
    }

    #[test]
    fn test_m_to_f_111() {
        // M_(1,1,1) = F_(1,1,1) since (1,1,1) has no proper refinements
        // (it is already the finest composition)
        let m: QSymFunction<i64> = QSymFunction::monomial_qsym(Composition::new(vec![1, 1, 1]));
        let in_f = m.to_fundamental_basis();
        assert_eq!(in_f.coefficient(&Composition::new(vec![1, 1, 1])), 1);
        assert_eq!(in_f.terms().len(), 1);
    }

    #[test]
    fn test_m_to_f_3() {
        // M_(3) in fundamental basis:
        // Refinements of (3): (3), (2,1), (1,2), (1,1,1)
        // M_(3) = F_(3) - F_(2,1) - F_(1,2) + F_(1,1,1)
        let m: QSymFunction<i64> = QSymFunction::monomial_qsym(Composition::new(vec![3]));
        let in_f = m.to_fundamental_basis();
        assert_eq!(in_f.coefficient(&Composition::new(vec![3])), 1);
        assert_eq!(in_f.coefficient(&Composition::new(vec![2, 1])), -1);
        assert_eq!(in_f.coefficient(&Composition::new(vec![1, 2])), -1);
        assert_eq!(in_f.coefficient(&Composition::new(vec![1, 1, 1])), 1);
    }

    #[test]
    fn test_roundtrip_all_compositions_of_4() {
        // Roundtrip for every composition of 4
        for alpha in Composition::integer_compositions(4) {
            let m: QSymFunction<i64> = QSymFunction::monomial_qsym(alpha.clone());
            let back = m.to_fundamental_basis().to_monomial_basis();
            assert_eq!(back, m, "roundtrip failed for {}", alpha);
        }
    }
}
