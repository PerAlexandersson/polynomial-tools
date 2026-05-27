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
use crate::schur_qsym::CompositionTableau;

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

        (QSymBasis::QuasisymmetricSchur, QSymBasis::Fundamental) => {
            quasisymmetric_schur_to_fundamental(f)
        }
        (QSymBasis::DualImmaculate, QSymBasis::Fundamental) => dual_immaculate_to_fundamental(f),
        (QSymBasis::Fundamental, QSymBasis::QuasisymmetricSchur) => {
            fundamental_to_tableau_basis(f, QSymBasis::QuasisymmetricSchur)
        }
        (QSymBasis::Fundamental, QSymBasis::DualImmaculate) => {
            fundamental_to_tableau_basis(f, QSymBasis::DualImmaculate)
        }

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

        (QSymBasis::QuasisymmetricSchur, _) => {
            let in_f = quasisymmetric_schur_to_fundamental(f);
            convert(&in_f, target)
        }
        (QSymBasis::DualImmaculate, _) => {
            let in_f = dual_immaculate_to_fundamental(f);
            convert(&in_f, target)
        }
        (_, QSymBasis::QuasisymmetricSchur) => {
            let in_f = convert(f, QSymBasis::Fundamental);
            fundamental_to_tableau_basis(&in_f, QSymBasis::QuasisymmetricSchur)
        }
        (_, QSymBasis::DualImmaculate) => {
            let in_f = convert(f, QSymBasis::Fundamental);
            fundamental_to_tableau_basis(&in_f, QSymBasis::DualImmaculate)
        }

        _ => unreachable!("same-basis QSym conversion should have returned early"),
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

fn quasisymmetric_schur_to_fundamental<C: Ring>(f: &QSymFunction<C>) -> QSymFunction<C> {
    let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();

    for (alpha, coeff) in f.terms() {
        let expansion = quasisymmetric_schur_basis_element::<C>(alpha);
        for (beta, beta_coeff) in expansion.terms() {
            let entry = result_terms.entry(beta.clone()).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone() * beta_coeff.clone();
        }
    }

    QSymFunction::from_terms(QSymBasis::Fundamental, result_terms)
}

fn dual_immaculate_to_fundamental<C: Ring>(f: &QSymFunction<C>) -> QSymFunction<C> {
    let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();

    for (alpha, coeff) in f.terms() {
        let expansion = crate::schur_qsym::dual_immaculate::<C>(alpha.parts());
        for (beta, beta_coeff) in expansion.terms() {
            let entry = result_terms.entry(beta.clone()).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone() * beta_coeff.clone();
        }
    }

    QSymFunction::from_terms(QSymBasis::Fundamental, result_terms)
}

fn quasisymmetric_schur_basis_element<C: Ring>(alpha: &Composition) -> QSymFunction<C> {
    let n = alpha.size();
    if n == 0 {
        return QSymFunction::basis_element(QSymBasis::Fundamental, Composition::empty());
    }

    let tableaux = CompositionTableau::enumerate(alpha.parts(), n);
    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();

    for tableau in &tableaux {
        let entries: Vec<u32> = tableau
            .rows
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect();
        if entries.len() != n as usize {
            continue;
        }
        let mut seen = vec![false; n as usize + 1];
        let mut is_standard = true;
        for value in entries {
            if value == 0 || value > n || seen[value as usize] {
                is_standard = false;
                break;
            }
            seen[value as usize] = true;
        }
        if !is_standard {
            continue;
        }

        let descent = tableau.descent_composition();
        let entry = terms.entry(descent).or_insert_with(C::zero);
        *entry = entry.clone() + C::one();
    }

    QSymFunction::from_terms(QSymBasis::Fundamental, terms)
}

fn fundamental_to_tableau_basis<C: Ring>(
    f: &QSymFunction<C>,
    target: QSymBasis,
) -> QSymFunction<C> {
    assert_eq!(
        f.basis(),
        QSymBasis::Fundamental,
        "tableau-basis conversion expects fundamental input"
    );
    assert!(
        matches!(
            target,
            QSymBasis::QuasisymmetricSchur | QSymBasis::DualImmaculate
        ),
        "unsupported tableau basis target"
    );

    let mut by_degree: BTreeMap<u32, BTreeMap<Composition, C>> = BTreeMap::new();
    for (alpha, coeff) in f.terms() {
        by_degree
            .entry(alpha.size())
            .or_default()
            .insert(alpha.clone(), coeff.clone());
    }

    let mut result_terms = BTreeMap::new();
    for (degree, terms) in by_degree {
        let homogeneous = QSymFunction::from_terms(QSymBasis::Fundamental, terms);
        let expanded = expand_fundamental_homogeneous(&homogeneous, degree, target)
            .expect("failed to expand in requested QSym basis");
        for (alpha, coeff) in expanded {
            let entry = result_terms.entry(alpha).or_insert_with(C::zero);
            *entry = entry.clone() + coeff;
        }
    }

    QSymFunction::from_terms(target, result_terms)
}

fn expand_fundamental_homogeneous<C: Ring>(
    f: &QSymFunction<C>,
    degree: u32,
    target: QSymBasis,
) -> Option<BTreeMap<Composition, C>> {
    let comps = Composition::integer_compositions(degree);
    let mut basis_in_f: BTreeMap<Composition, QSymFunction<i64>> = BTreeMap::new();

    for alpha in comps {
        let expansion = match target {
            QSymBasis::QuasisymmetricSchur => quasisymmetric_schur_basis_element::<i64>(&alpha),
            QSymBasis::DualImmaculate => crate::schur_qsym::dual_immaculate::<i64>(alpha.parts()),
            _ => unreachable!("unsupported tableau basis target"),
        };
        basis_in_f.insert(alpha, expansion);
    }

    gaussian_expand_generic(f, &basis_in_f)
}

fn gaussian_expand_generic<C: Ring>(
    f_fund: &QSymFunction<C>,
    basis_in_f: &BTreeMap<Composition, QSymFunction<i64>>,
) -> Option<BTreeMap<Composition, C>> {
    let mut all_f_comps: BTreeMap<Composition, usize> = BTreeMap::new();
    let mut idx = 0;
    for basis_element in basis_in_f.values() {
        for comp in basis_element.terms().keys() {
            if !all_f_comps.contains_key(comp) {
                all_f_comps.insert(comp.clone(), idx);
                idx += 1;
            }
        }
    }
    for comp in f_fund.terms().keys() {
        if !all_f_comps.contains_key(comp) {
            all_f_comps.insert(comp.clone(), idx);
            idx += 1;
        }
    }

    let rows = all_f_comps.len();
    let basis_list: Vec<Composition> = basis_in_f.keys().cloned().collect();
    let cols = basis_list.len();

    let mut mat = vec![vec![0i128; cols]; rows];
    let mut rhs = vec![C::zero(); rows];

    for (col, basis_comp) in basis_list.iter().enumerate() {
        let basis_element = basis_in_f.get(basis_comp).unwrap();
        for (fund_comp, &coeff) in basis_element.terms() {
            let row = all_f_comps[fund_comp];
            mat[row][col] = coeff as i128;
        }
    }
    for (fund_comp, coeff) in f_fund.terms() {
        let row = all_f_comps[fund_comp];
        rhs[row] = coeff.clone();
    }

    let mut pivot_row = 0usize;
    let mut pivot_cols = Vec::new();
    for col in 0..cols {
        let Some(found) = (pivot_row..rows).find(|&row| mat[row][col] != 0) else {
            continue;
        };
        mat.swap(pivot_row, found);
        rhs.swap(pivot_row, found);
        pivot_cols.push((pivot_row, col));

        let pivot = mat[pivot_row][col];
        for row in 0..rows {
            if row == pivot_row {
                continue;
            }
            let factor = mat[row][col];
            if factor == 0 {
                continue;
            }

            for c in 0..cols {
                mat[row][c] = mat[row][c] * pivot - factor * mat[pivot_row][c];
            }
            rhs[row] = scale_ring_i128(&rhs[row], pivot) - scale_ring_i128(&rhs[pivot_row], factor);

            let divisor = row_gcd(&mat[row]);
            if divisor > 1 {
                for c in 0..cols {
                    mat[row][c] /= divisor;
                }
                rhs[row] = rhs[row].exact_div_i64(i128_to_i64(divisor));
            }
        }

        pivot_row += 1;
    }

    for row in 0..rows {
        if mat[row].iter().all(|&entry| entry == 0) && !rhs[row].is_zero() {
            return None;
        }
    }

    if pivot_cols.len() != cols {
        return None;
    }

    let mut coeffs = BTreeMap::new();
    for (row, col) in pivot_cols {
        let pivot = mat[row][col];
        let coeff = rhs[row].exact_div_i64(i128_to_i64(pivot));
        if !coeff.is_zero() {
            coeffs.insert(basis_list[col].clone(), coeff);
        }
    }

    Some(coeffs)
}

fn scale_ring_i128<C: Ring>(value: &C, scalar: i128) -> C {
    if scalar == 0 || value.is_zero() {
        C::zero()
    } else {
        value.clone() * C::from_i64(i128_to_i64(scalar))
    }
}

fn i128_to_i64(value: i128) -> i64 {
    assert!(
        value >= i64::MIN as i128 && value <= i64::MAX as i128,
        "integer matrix coefficient out of i64 range"
    );
    value as i64
}

fn gcd_i128(mut a: i128, mut b: i128) -> i128 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

fn row_gcd(row: &[i128]) -> i128 {
    row.iter().copied().fold(0, gcd_i128)
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

    #[test]
    fn test_quasisymmetric_schur_roundtrip_degree_4() {
        for alpha in Composition::integer_compositions(4) {
            let qs: QSymFunction<i64> =
                QSymFunction::basis_element(QSymBasis::QuasisymmetricSchur, alpha.clone());
            let back = qs.to_fundamental_basis().to_quasisymmetric_schur_basis();
            assert_eq!(back, qs, "roundtrip failed for {}", alpha);
        }
    }

    #[test]
    fn test_dual_immaculate_roundtrip_degree_4() {
        for alpha in Composition::integer_compositions(4) {
            let di: QSymFunction<i64> =
                QSymFunction::basis_element(QSymBasis::DualImmaculate, alpha.clone());
            let back = di.to_fundamental_basis().to_dual_immaculate_basis();
            assert_eq!(back, di, "roundtrip failed for {}", alpha);
        }
    }

    #[test]
    fn test_polynomial_coefficients_expand_to_new_bases() {
        use sym_poly_core::UnivariatePolynomial;

        let mut terms = BTreeMap::new();
        terms.insert(
            Composition::new(vec![2, 1]),
            UnivariatePolynomial::<i64>::new(vec![1, 1]),
        );
        let f = QSymFunction::from_terms(QSymBasis::Fundamental, terms);
        let _qs = f.to_quasisymmetric_schur_basis();
        let _di = f.to_dual_immaculate_basis();
    }
}
