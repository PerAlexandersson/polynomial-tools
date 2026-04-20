//! QSym power sum bases Ψ (type 1) and Φ (type 2).
//!
//! Two quasisymmetric power sum bases indexed by compositions, as defined in:
//!
//! > Ballantine, Daugherty, Hicks, Mason, Niese.
//! > *Quasisymmetric Power Sums.* J. Combin. Theory Ser. A (2020).
//! > <https://doi.org/10.1016/j.jcta.2020.105273>
//!
//! Both are defined via the same structural formula over coarsenings of α,
//! differing only in the block-weight function:
//!
//! ```text
//! Ψ_α = z_α  Σ_{β coarsens α}  1 / Π_B π(B)   · M_{sum(β)}
//! Φ_α = z_α  Σ_{β coarsens α}  1 / Π_B sp(B)  · M_{sum(β)}
//! ```
//!
//! where:
//! - β ranges over all coarsenings of α (merging adjacent parts)
//! - π(B) = Π_{i=1}^m (b_1 + ... + b_i) — product of partial sums
//! - sp(B) = m! · Π_{i=1}^m b_i — factorial of block length times product of parts
//! - z_α = Π i^{m_i} · m_i! — standard z-coefficient
//!
//! Both Ψ ↔ M and Φ ↔ M conversions involve rational coefficients.

use std::collections::BTreeMap;

use sym_poly_core::{Composition, Ring};

use crate::basis::QSymBasis;
use crate::qsym_function::QSymFunction;

// =========================================================================
// Public API
// =========================================================================

/// Compute Ψ_α (type 1 power sum) expanded in the monomial basis.
///
/// Requires rational coefficients since the expansion involves 1/π(B) factors.
pub fn psi_in_monomial_basis<C: Ring>(alpha: &Composition) -> QSymFunction<C> {
    power_sum_in_monomial_basis::<C>(alpha, pi_value)
}

/// Compute Ψ̃_α = Ψ_α / z_α (normalized type 1 power sum) in the monomial basis.
///
/// Ψ̃_α = Σ_{β coarsens α} (1 / Π_B π(B)) · M_{sum(β)}
///
/// The normalized basis has the property that for naturally labeled posets P,
/// ω(Γ(P)) is Ψ̃-positive with non-negative integer coefficients.
pub fn psi_normalized_in_monomial_basis<C: Ring>(alpha: &Composition) -> QSymFunction<C> {
    normalized_power_sum_in_monomial_basis::<C>(alpha, pi_value)
}

/// Compute Φ̃_α = Φ_α / z_α (normalized type 2 power sum) in the monomial basis.
pub fn phi_normalized_in_monomial_basis<C: Ring>(alpha: &Composition) -> QSymFunction<C> {
    normalized_power_sum_in_monomial_basis::<C>(alpha, sp_value)
}

/// Compute Φ_α (type 2 power sum) expanded in the monomial basis.
///
/// Requires rational coefficients since the expansion involves 1/sp(B) factors.
pub fn phi_in_monomial_basis<C: Ring>(alpha: &Composition) -> QSymFunction<C> {
    power_sum_in_monomial_basis::<C>(alpha, sp_value)
}

/// Convert a QSymFunction from Ψ basis to M basis.
pub fn psi_to_monomial<C: Ring>(f: &QSymFunction<C>) -> QSymFunction<C> {
    expand_to_monomial(f, pi_value)
}

/// Convert a QSymFunction from Φ basis to M basis.
pub fn phi_to_monomial<C: Ring>(f: &QSymFunction<C>) -> QSymFunction<C> {
    expand_to_monomial(f, sp_value)
}

/// Convert a QSymFunction from M basis to Ψ basis.
///
/// Uses per-degree rational matrix inversion.
/// Requires rational coefficients (Ratio<BigInt>, Ratio<i64>).
pub fn monomial_to_psi<C: Ring>(f: &QSymFunction<C>) -> QSymFunction<C> {
    monomial_to_power_sum(f, QSymBasis::PowerSumPsi, pi_value)
}

/// Convert a QSymFunction from M basis to Φ basis.
///
/// Uses per-degree rational matrix inversion.
/// Requires rational coefficients (Ratio<BigInt>, Ratio<i64>).
pub fn monomial_to_phi<C: Ring>(f: &QSymFunction<C>) -> QSymFunction<C> {
    monomial_to_power_sum(f, QSymBasis::PowerSumPhi, sp_value)
}

// =========================================================================
// Block-weight functions
// =========================================================================

/// π(B) = product of partial sums of block B.
///
/// For B = (b_1, ..., b_m): π(B) = b_1 · (b_1+b_2) · ... · (b_1+...+b_m).
fn pi_value(block: &[u32]) -> i64 {
    let mut sum = 0i64;
    let mut product = 1i64;
    for &b in block {
        sum += b as i64;
        product *= sum;
    }
    product
}

/// sp(B) = |B|! · Π b_i — factorial of block length times product of parts.
///
/// For B = (b_1, ..., b_m): sp(B) = m! · b_1 · b_2 · ... · b_m.
fn sp_value(block: &[u32]) -> i64 {
    let m = block.len() as i64;
    let mut result = 1i64;
    for k in 1..=m {
        result *= k;
    }
    for &b in block {
        result *= b as i64;
    }
    result
}

// =========================================================================
// Generic implementation (shared by Ψ and Φ)
// =========================================================================

/// Normalized power sum: Σ_{β coarsens α} 1/Π_B w(B) · M_{sum(β)} (no z_α factor).
fn normalized_power_sum_in_monomial_basis<C: Ring>(
    alpha: &Composition,
    block_weight: fn(&[u32]) -> i64,
) -> QSymFunction<C> {
    if alpha.is_empty() {
        return QSymFunction::basis_element(QSymBasis::Monomial, Composition::empty());
    }

    let parts = alpha.parts();
    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();

    for coarsening in composition_coarsenings(parts) {
        let mut result_parts = Vec::new();
        let mut denom: i64 = 1;

        for block in &coarsening {
            result_parts.push(block.iter().sum::<u32>());
            denom *= block_weight(block);
        }

        let comp = Composition::new(result_parts);
        let coeff = C::one().exact_div_i64(denom);

        let entry = terms.entry(comp).or_insert_with(C::zero);
        *entry = entry.clone() + coeff;
    }

    QSymFunction::from_terms(QSymBasis::Monomial, terms)
}

/// Generic power sum expansion: z_α · Σ_{β coarsens α} 1/Π_B w(B) · M_{sum(β)}.
fn power_sum_in_monomial_basis<C: Ring>(
    alpha: &Composition,
    block_weight: fn(&[u32]) -> i64,
) -> QSymFunction<C> {
    if alpha.is_empty() {
        return QSymFunction::basis_element(QSymBasis::Monomial, Composition::empty());
    }

    let parts = alpha.parts();
    let z = z_coefficient(parts);

    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();

    for coarsening in composition_coarsenings(parts) {
        let mut result_parts = Vec::new();
        let mut denom: i64 = 1;

        for block in &coarsening {
            result_parts.push(block.iter().sum::<u32>());
            denom *= block_weight(block);
        }

        let comp = Composition::new(result_parts);
        let coeff = C::from_i64(z).exact_div_i64(denom);

        let entry = terms.entry(comp).or_insert_with(C::zero);
        *entry = entry.clone() + coeff;
    }

    QSymFunction::from_terms(QSymBasis::Monomial, terms)
}

/// Expand from a power sum basis to monomial basis.
fn expand_to_monomial<C: Ring>(
    f: &QSymFunction<C>,
    block_weight: fn(&[u32]) -> i64,
) -> QSymFunction<C> {
    let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();

    for (alpha, coeff) in f.terms() {
        let expansion = power_sum_in_monomial_basis::<C>(alpha, block_weight);
        for (beta, c) in expansion.terms() {
            let entry = result_terms.entry(beta.clone()).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone() * c.clone();
        }
    }

    QSymFunction::from_terms(QSymBasis::Monomial, result_terms)
}

/// Convert from monomial to a power sum basis via per-degree matrix inversion.
fn monomial_to_power_sum<C: Ring>(
    f: &QSymFunction<C>,
    target_basis: QSymBasis,
    block_weight: fn(&[u32]) -> i64,
) -> QSymFunction<C> {
    if f.is_zero() {
        return QSymFunction::zero(target_basis);
    }

    let mut by_degree: BTreeMap<u32, Vec<(Composition, C)>> = BTreeMap::new();
    for (comp, coeff) in f.terms() {
        by_degree
            .entry(comp.size())
            .or_default()
            .push((comp.clone(), coeff.clone()));
    }

    let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();

    for (deg, terms) in by_degree {
        let compositions = Composition::integer_compositions(deg);
        let k = compositions.len();
        let comp_index: BTreeMap<&Composition, usize> = compositions
            .iter()
            .enumerate()
            .map(|(i, c)| (c, i))
            .collect();

        // Build power-sum → M transition matrix as (num, den) pairs
        let ps_to_m = build_transition_matrix(&compositions, &comp_index, block_weight);

        // Invert to get M → power-sum matrix
        let m_to_ps = invert_rational_matrix(&ps_to_m);

        // Build input vector
        let mut input = vec![C::zero(); k];
        for (comp, coeff) in &terms {
            if let Some(&idx) = comp_index.get(comp) {
                input[idx] = coeff.clone();
            }
        }

        // Matrix-vector multiply
        for i in 0..k {
            let mut sum = C::zero();
            for j in 0..k {
                if !input[j].is_zero() {
                    let (num, den) = m_to_ps[j][i];
                    if num != 0 {
                        let term = input[j].clone() * C::from_i64(num);
                        sum = sum + term.exact_div_i64(den);
                    }
                }
            }
            if !sum.is_zero() {
                result_terms.insert(compositions[i].clone(), sum);
            }
        }
    }

    QSymFunction::from_terms(target_basis, result_terms)
}

// =========================================================================
// Combinatorial helpers
// =========================================================================

/// z_α = Π i^{m_i} · m_i! where m_i = #{j : α_j = i}.
fn z_coefficient(parts: &[u32]) -> i64 {
    let mut counts: BTreeMap<u32, u32> = BTreeMap::new();
    for &p in parts {
        *counts.entry(p).or_insert(0) += 1;
    }
    let mut z: i64 = 1;
    for (&val, &mult) in &counts {
        for _ in 0..mult {
            z *= val as i64;
        }
        for k in 1..=mult as i64 {
            z *= k;
        }
    }
    z
}

/// All coarsenings of a composition (merging adjacent parts).
///
/// For k parts there are 2^{k-1} coarsenings, one for each subset of
/// the k-1 gaps between consecutive parts.
fn composition_coarsenings(parts: &[u32]) -> Vec<Vec<Vec<u32>>> {
    let k = parts.len();
    if k == 0 {
        return vec![vec![]];
    }
    let num_gaps = k - 1;
    let mut result = Vec::with_capacity(1 << num_gaps);

    for mask in 0..(1u32 << num_gaps) {
        let mut coarsening = Vec::new();
        let mut block = vec![parts[0]];
        for i in 0..num_gaps {
            if mask & (1 << i) != 0 {
                block.push(parts[i + 1]);
            } else {
                coarsening.push(block);
                block = vec![parts[i + 1]];
            }
        }
        coarsening.push(block);
        result.push(coarsening);
    }

    result
}

// =========================================================================
// Rational matrix utilities
// =========================================================================

/// Build the power-sum → M transition matrix as (num, den) pairs.
fn build_transition_matrix(
    compositions: &[Composition],
    comp_index: &BTreeMap<&Composition, usize>,
    block_weight: fn(&[u32]) -> i64,
) -> Vec<Vec<(i64, i64)>> {
    let k = compositions.len();
    let mut matrix = vec![vec![(0i64, 1i64); k]; k];

    for (j, alpha) in compositions.iter().enumerate() {
        let parts = alpha.parts();
        let z = z_coefficient(parts);

        for coarsening in composition_coarsenings(parts) {
            let mut result_parts = Vec::new();
            let mut denom: i64 = 1;

            for block in &coarsening {
                result_parts.push(block.iter().sum::<u32>());
                denom *= block_weight(block);
            }

            let comp = Composition::new(result_parts);
            if let Some(&i) = comp_index.get(&comp) {
                let (ref mut num, ref mut den) = matrix[j][i];
                *num = (*num) * denom + z * (*den);
                *den = (*den) * denom;
                let g = gcd(num.unsigned_abs(), den.unsigned_abs()) as i64;
                if g > 1 {
                    *num /= g;
                    *den /= g;
                }
                if *den < 0 {
                    *num = -(*num);
                    *den = -(*den);
                }
            }
        }
    }

    matrix
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Invert a rational matrix (entries as (num, den) pairs) via Gaussian elimination.
fn invert_rational_matrix(m: &[Vec<(i64, i64)>]) -> Vec<Vec<(i64, i64)>> {
    use num_rational::Ratio;
    type Q = Ratio<i64>;

    let n = m.len();
    if n == 0 {
        return vec![];
    }

    let mut aug: Vec<Vec<Q>> = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = Vec::with_capacity(2 * n);
        for j in 0..n {
            let (num, den) = m[i][j];
            row.push(Q::new(num, den));
        }
        for j in 0..n {
            row.push(if i == j {
                Q::from_integer(1)
            } else {
                Q::from_integer(0)
            });
        }
        aug.push(row);
    }

    for col in 0..n {
        let mut pivot = None;
        for row in col..n {
            if aug[row][col] != Q::from_integer(0) {
                pivot = Some(row);
                break;
            }
        }
        let pivot = pivot.expect("power-sum → M matrix is singular");
        aug.swap(col, pivot);

        let diag = aug[col][col].clone();
        for j in 0..2 * n {
            aug[col][j] = aug[col][j].clone() / diag.clone();
        }

        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = aug[row][col].clone();
            for j in 0..2 * n {
                let val = aug[col][j].clone() * factor.clone();
                aug[row][j] = aug[row][j].clone() - val;
            }
        }
    }

    let mut inv = vec![vec![(0i64, 1i64); n]; n];
    for i in 0..n {
        for j in 0..n {
            let val = &aug[i][n + j];
            inv[i][j] = (*val.numer(), *val.denom());
        }
    }
    inv
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Ratio;
    type Q = Ratio<i64>;

    // -- Helper tests --

    #[test]
    fn test_pi_value() {
        assert_eq!(pi_value(&[2, 1]), 6); // 2 * 3
        assert_eq!(pi_value(&[3]), 3);
        assert_eq!(pi_value(&[1, 1, 1]), 6); // 1 * 2 * 3
    }

    #[test]
    fn test_sp_value() {
        assert_eq!(sp_value(&[2, 1]), 2 * 2 * 1); // 2! * 2 * 1 = 4
        assert_eq!(sp_value(&[3]), 1 * 3); // 1! * 3 = 3
        assert_eq!(sp_value(&[1, 1, 1]), 6 * 1 * 1 * 1); // 3! * 1 * 1 * 1 = 6
    }

    #[test]
    fn test_z_coefficient() {
        assert_eq!(z_coefficient(&[2, 1]), 2);
        assert_eq!(z_coefficient(&[1, 1, 1]), 6);
        assert_eq!(z_coefficient(&[3]), 3);
        assert_eq!(z_coefficient(&[2, 2]), 8);
    }

    #[test]
    fn test_coarsenings_count() {
        assert_eq!(composition_coarsenings(&[2, 1]).len(), 2);
        assert_eq!(composition_coarsenings(&[1, 1, 1]).len(), 4);
        assert_eq!(composition_coarsenings(&[3]).len(), 1);
    }

    // -- Ψ tests --

    #[test]
    fn test_psi_single_part() {
        // Ψ_(n) = M_(n)
        let psi3: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![3]));
        assert_eq!(
            psi3.coefficient(&Composition::new(vec![3])),
            Q::from_integer(1)
        );
        assert_eq!(psi3.terms().len(), 1);
    }

    #[test]
    fn test_psi_21_ne_psi_12() {
        let psi21: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![2, 1]));
        let psi12: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![1, 2]));
        assert_ne!(psi21, psi12);

        // Ψ_(2,1) = M_(2,1) + 1/3 M_(3)
        assert_eq!(
            psi21.coefficient(&Composition::new(vec![2, 1])),
            Q::from_integer(1)
        );
        assert_eq!(psi21.coefficient(&Composition::new(vec![3])), Q::new(1, 3));

        // Ψ_(1,2) = M_(1,2) + 2/3 M_(3)
        assert_eq!(
            psi12.coefficient(&Composition::new(vec![1, 2])),
            Q::from_integer(1)
        );
        assert_eq!(psi12.coefficient(&Composition::new(vec![3])), Q::new(2, 3));
    }

    #[test]
    fn test_psi_11() {
        // Ψ_(1,1) = 2*M_(1,1) + M_(2)
        let psi11: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![1, 1]));
        assert_eq!(
            psi11.coefficient(&Composition::new(vec![1, 1])),
            Q::from_integer(2)
        );
        assert_eq!(
            psi11.coefficient(&Composition::new(vec![2])),
            Q::from_integer(1)
        );
    }

    #[test]
    fn test_psi_111() {
        // Ψ_(1,1,1) = 6*M_(1,1,1) + 3*M_(2,1) + 3*M_(1,2) + M_(3)
        let psi: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![1, 1, 1]));
        assert_eq!(
            psi.coefficient(&Composition::new(vec![1, 1, 1])),
            Q::from_integer(6)
        );
        assert_eq!(
            psi.coefficient(&Composition::new(vec![2, 1])),
            Q::from_integer(3)
        );
        assert_eq!(
            psi.coefficient(&Composition::new(vec![1, 2])),
            Q::from_integer(3)
        );
        assert_eq!(
            psi.coefficient(&Composition::new(vec![3])),
            Q::from_integer(1)
        );
    }

    #[test]
    fn test_psi_roundtrip_degree3() {
        for alpha in Composition::integer_compositions(3) {
            let m: QSymFunction<Q> = QSymFunction::monomial_qsym(alpha.clone());
            let in_psi = monomial_to_psi(&m);
            let back = psi_to_monomial(&in_psi);
            assert_eq!(
                back.coefficient(&alpha),
                Q::from_integer(1),
                "Ψ roundtrip failed for M_{}",
                alpha
            );
            assert_eq!(back.terms().len(), 1, "Ψ extra terms for M_{}", alpha);
        }
    }

    // -- Φ tests --

    #[test]
    fn test_phi_single_part() {
        // Φ_(n) should also equal M_(n) since the only coarsening of (n) is itself,
        // and both π({n}) = n and sp({n}) = 1!*n = n.
        let phi3: QSymFunction<Q> = phi_in_monomial_basis(&Composition::new(vec![3]));
        assert_eq!(
            phi3.coefficient(&Composition::new(vec![3])),
            Q::from_integer(1)
        );
        assert_eq!(phi3.terms().len(), 1);
    }

    #[test]
    fn test_phi_21_ne_phi_12() {
        let phi21: QSymFunction<Q> = phi_in_monomial_basis(&Composition::new(vec![2, 1]));
        let phi12: QSymFunction<Q> = phi_in_monomial_basis(&Composition::new(vec![1, 2]));
        assert_ne!(phi21, phi12);
    }

    #[test]
    fn test_phi_21() {
        // Φ_(2,1): coarsenings of (2,1):
        //   {(2),(1)}: sp({2})=2, sp({1})=1. denom=2. coeff = z/2 = 2/2 = 1. comp=(2,1).
        //   {(2,1)}:   sp({2,1})=2!*2*1=4. denom=4. coeff = z/4 = 2/4 = 1/2. comp=(3).
        // Φ_(2,1) = M_(2,1) + 1/2 M_(3)
        let phi21: QSymFunction<Q> = phi_in_monomial_basis(&Composition::new(vec![2, 1]));
        assert_eq!(
            phi21.coefficient(&Composition::new(vec![2, 1])),
            Q::from_integer(1)
        );
        assert_eq!(phi21.coefficient(&Composition::new(vec![3])), Q::new(1, 2));
    }

    #[test]
    fn test_phi_12() {
        // Φ_(1,2): coarsenings of (1,2):
        //   {(1),(2)}: sp({1})=1, sp({2})=2. denom=2. coeff = 2/2 = 1. comp=(1,2).
        //   {(1,2)}:   sp({1,2})=2!*1*2=4. denom=4. coeff = 2/4 = 1/2. comp=(3).
        // Φ_(1,2) = M_(1,2) + 1/2 M_(3)
        let phi12: QSymFunction<Q> = phi_in_monomial_basis(&Composition::new(vec![1, 2]));
        assert_eq!(
            phi12.coefficient(&Composition::new(vec![1, 2])),
            Q::from_integer(1)
        );
        assert_eq!(phi12.coefficient(&Composition::new(vec![3])), Q::new(1, 2));
    }

    #[test]
    fn test_phi_roundtrip_degree3() {
        for alpha in Composition::integer_compositions(3) {
            let m: QSymFunction<Q> = QSymFunction::monomial_qsym(alpha.clone());
            let in_phi = monomial_to_phi(&m);
            let back = phi_to_monomial(&in_phi);
            assert_eq!(
                back.coefficient(&alpha),
                Q::from_integer(1),
                "Φ roundtrip failed for M_{}",
                alpha
            );
            assert_eq!(back.terms().len(), 1, "Φ extra terms for M_{}", alpha);
        }
    }

    // -- Cross-check: Ψ ≠ Φ --

    #[test]
    fn test_psi_ne_phi() {
        let psi: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![2, 1]));
        let phi: QSymFunction<Q> = phi_in_monomial_basis(&Composition::new(vec![2, 1]));
        assert_ne!(psi, phi);
    }

    // -- Omega involution --

    #[test]
    fn test_omega_involution() {
        // ω² = id on a fundamental basis element
        let f: QSymFunction<Q> = QSymFunction::fundamental_qsym(Composition::new(vec![2, 1]));
        let omega_f = f.omega_involution();
        let back = omega_f.omega_involution();
        assert_eq!(back, f, "ω² should be identity");
    }

    #[test]
    fn test_omega_on_psi() {
        // ω(Ψ_α) = (-1)^{n-ℓ(α)} Ψ_{α^r}
        // For α = (2,1), n=3, ℓ=2: sign = (-1)^{3-2} = -1, α^r = (1,2)
        let psi21: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![2, 1]));
        let omega_psi = psi21.omega_involution();
        let psi12: QSymFunction<Q> = psi_in_monomial_basis(&Composition::new(vec![1, 2]));
        let expected = psi12.scale(&Q::from_integer(-1));
        assert_eq!(
            omega_psi.to_monomial_basis(),
            expected.to_monomial_basis(),
            "ω(Ψ_(2,1)) should be -Ψ_(1,2)"
        );
    }

    #[test]
    fn test_omega_psi_identity_degree3() {
        // Verify ω(Ψ_α) = (-1)^{n-ℓ(α)} Ψ_{α^r} for all compositions of 3
        for alpha in Composition::integer_compositions(3) {
            let psi_alpha: QSymFunction<Q> = psi_in_monomial_basis(&alpha);
            let omega_psi = psi_alpha.omega_involution().to_monomial_basis();

            let n = alpha.size();
            let ell = alpha.num_parts() as u32;
            let sign = if (n - ell) % 2 == 0 {
                Q::from_integer(1)
            } else {
                Q::from_integer(-1)
            };
            let reversed_parts: Vec<u32> = alpha.parts().iter().rev().copied().collect();
            let alpha_r = Composition::new(reversed_parts);
            let psi_rev: QSymFunction<Q> = psi_in_monomial_basis(&alpha_r);
            let expected = psi_rev.scale(&sign).to_monomial_basis();

            assert_eq!(
                omega_psi, expected,
                "ω(Ψ_{}) should be ({})Ψ_{}",
                alpha, sign, alpha_r
            );
        }
    }

    // -- Normalized Ψ̃ = Ψ/z --

    #[test]
    fn test_psi_normalized_21() {
        // Ψ̃_(2,1) = Ψ_(2,1)/z_(2,1) = (M_(2,1) + 1/3 M_(3)) / 2
        //          = 1/2 M_(2,1) + 1/6 M_(3)
        let psi_n: QSymFunction<Q> =
            psi_normalized_in_monomial_basis(&Composition::new(vec![2, 1]));
        assert_eq!(
            psi_n.coefficient(&Composition::new(vec![2, 1])),
            Q::new(1, 2)
        );
        assert_eq!(psi_n.coefficient(&Composition::new(vec![3])), Q::new(1, 6));
    }

    #[test]
    fn test_p_partition_psi_normalized_positive() {
        // Γ(P) is Ψ̃-positive for naturally labeled posets:
        //   Γ(P) = Σ c_α · Ψ̃_α  with c_α ∈ ℤ_{≥0}
        // i.e., the Ψ_α coefficient of Γ(P) times z_α is a non-negative integer.
        let posets: Vec<(usize, Vec<(usize, usize)>)> = vec![
            (3, vec![(0, 1), (1, 2)]),         // chain
            (3, vec![]),                       // antichain
            (3, vec![(0, 2), (1, 2)]),         // V-poset
            (4, vec![(0, 1), (1, 2), (2, 3)]), // chain 4
            (4, vec![(0, 2), (1, 3)]),         // two disjoint edges
            (4, vec![(0, 1), (0, 2), (0, 3)]), // star
        ];
        for (n, covers) in &posets {
            let gamma: QSymFunction<Q> =
                crate::p_partition::p_partition_generating_function(*n, covers);
            let in_psi = gamma.to_monomial_basis().to_basis(QSymBasis::PowerSumPsi);
            for (alpha, coeff) in in_psi.terms() {
                let z = z_coefficient(alpha.parts());
                let scaled = coeff.clone() * Q::from_integer(z);
                assert!(
                    *scaled.denom() == 1,
                    "poset {:?}: z·coeff of Ψ_{} not integer: {}",
                    covers,
                    alpha,
                    scaled
                );
                assert!(
                    *scaled.numer() >= 0,
                    "poset {:?}: z·coeff of Ψ_{} negative: {}",
                    covers,
                    alpha,
                    scaled
                );
            }
        }
    }
}
