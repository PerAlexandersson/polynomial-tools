//! Macdonald eigenoperators on the modified Macdonald basis.
//!
//! This module implements the diagonal data for the modified Macdonald basis
//! `Htilde_mu`: `nabla`, `Delta_f`, and `Delta'_f`. Conversion between Schur
//! and modified Macdonald bases is a later layer; the functions here are the
//! exact operator core used once such expansions are available.

use std::collections::BTreeMap;

use num_rational::Ratio;
use sym_poly_core::{Partition, Ring, UnivariatePolynomial};

use crate::SymmetricFunction;

pub type Rational = Ratio<i64>;

/// Polynomials in `q,t`, represented as polynomials in `t` with coefficients in `q`.
pub type QtPolynomial = UnivariatePolynomial<UnivariatePolynomial<Rational>>;

/// A finite expansion in the modified Macdonald basis `Htilde_mu`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModifiedMacdonaldExpansion {
    terms: BTreeMap<Partition, QtPolynomial>,
}

impl ModifiedMacdonaldExpansion {
    pub fn zero() -> Self {
        Self {
            terms: BTreeMap::new(),
        }
    }

    pub fn from_terms(terms: BTreeMap<Partition, QtPolynomial>) -> Self {
        let mut result = Self { terms };
        result.strip_zeros();
        result
    }

    pub fn basis_element(partition: Partition) -> Self {
        Self::scaled_basis_element(partition, <QtPolynomial as Ring>::one())
    }

    pub fn scaled_basis_element(partition: Partition, coefficient: QtPolynomial) -> Self {
        if coefficient.is_zero() {
            return Self::zero();
        }
        Self::from_terms(BTreeMap::from([(partition, coefficient)]))
    }

    pub fn terms(&self) -> &BTreeMap<Partition, QtPolynomial> {
        &self.terms
    }

    pub fn coefficient(&self, partition: &Partition) -> QtPolynomial {
        self.terms
            .get(partition)
            .cloned()
            .unwrap_or_else(<QtPolynomial as Ring>::zero)
    }

    /// Apply `nabla(Htilde_mu) = T_mu Htilde_mu`.
    pub fn nabla(&self) -> Self {
        scale_modified_macdonald_expansion(self, |partition| nabla_eigenvalue(partition))
    }

    /// Apply `Delta_f(Htilde_mu) = f[B_mu] Htilde_mu`.
    pub fn delta(&self, f: &SymmetricFunction<QtPolynomial>) -> Self {
        scale_modified_macdonald_expansion(self, |partition| delta_eigenvalue(f, partition))
    }

    /// Apply `Delta'_f(Htilde_mu) = f[B_mu - 1] Htilde_mu`.
    pub fn delta_prime(&self, f: &SymmetricFunction<QtPolynomial>) -> Self {
        scale_modified_macdonald_expansion(self, |partition| delta_prime_eigenvalue(f, partition))
    }

    fn strip_zeros(&mut self) {
        self.terms.retain(|_, coeff| !coeff.is_zero());
    }
}

pub fn qt_constant(value: i64) -> QtPolynomial {
    <QtPolynomial as Ring>::from_i64(value)
}

pub fn qt_monomial(q_degree: usize, t_degree: usize) -> QtPolynomial {
    let q_part = UnivariatePolynomial::monomial(q_degree, Rational::from_integer(1));
    UnivariatePolynomial::monomial(t_degree, q_part)
}

pub fn qt_coefficient(polynomial: &QtPolynomial, q_degree: usize, t_degree: usize) -> Rational {
    polynomial.coeff(t_degree).coeff(q_degree)
}

/// The alphabet `B_mu = sum_{c in mu} q^{a'(c)} t^{l'(c)}` as monomial values.
pub fn macdonald_b_alphabet(partition: &Partition) -> Vec<QtPolynomial> {
    partition
        .diagram_boxes()
        .into_iter()
        .map(|(row, col)| qt_monomial(col, row))
        .collect()
}

/// The scalar `B_mu(q,t) = sum_{c in mu} q^{a'(c)} t^{l'(c)}`.
pub fn macdonald_b_eigenvalue(partition: &Partition) -> QtPolynomial {
    macdonald_b_power_sum(partition, 1)
}

/// The scalar `T_mu = q^{n(mu')} t^{n(mu)}` for `nabla`.
pub fn nabla_eigenvalue(partition: &Partition) -> QtPolynomial {
    let q_degree = partition.conjugate_partition().partition_n() as usize;
    let t_degree = partition.partition_n() as usize;
    qt_monomial(q_degree, t_degree)
}

/// The scalar by which `Delta_f` acts on `Htilde_mu`.
pub fn delta_eigenvalue(
    f: &SymmetricFunction<QtPolynomial>,
    partition: &Partition,
) -> QtPolynomial {
    plethystic_evaluate_from_power_sums(f, |power| macdonald_b_power_sum(partition, power))
}

/// The scalar by which `Delta'_f` acts on `Htilde_mu`.
pub fn delta_prime_eigenvalue(
    f: &SymmetricFunction<QtPolynomial>,
    partition: &Partition,
) -> QtPolynomial {
    plethystic_evaluate_from_power_sums(f, |power| {
        macdonald_b_power_sum(partition, power) - <QtPolynomial as Ring>::one()
    })
}

pub fn nabla_modified_macdonald(
    expansion: &ModifiedMacdonaldExpansion,
) -> ModifiedMacdonaldExpansion {
    expansion.nabla()
}

pub fn delta_modified_macdonald(
    f: &SymmetricFunction<QtPolynomial>,
    expansion: &ModifiedMacdonaldExpansion,
) -> ModifiedMacdonaldExpansion {
    expansion.delta(f)
}

pub fn delta_prime_modified_macdonald(
    f: &SymmetricFunction<QtPolynomial>,
    expansion: &ModifiedMacdonaldExpansion,
) -> ModifiedMacdonaldExpansion {
    expansion.delta_prime(f)
}

fn scale_modified_macdonald_expansion<F>(
    expansion: &ModifiedMacdonaldExpansion,
    mut eigenvalue: F,
) -> ModifiedMacdonaldExpansion
where
    F: FnMut(&Partition) -> QtPolynomial,
{
    let terms = expansion
        .terms()
        .iter()
        .map(|(partition, coefficient)| {
            (
                partition.clone(),
                coefficient.clone() * eigenvalue(partition),
            )
        })
        .collect();
    ModifiedMacdonaldExpansion::from_terms(terms)
}

fn macdonald_b_power_sum(partition: &Partition, power: u32) -> QtPolynomial {
    let mut result = <QtPolynomial as Ring>::zero();
    for (row, col) in partition.diagram_boxes() {
        result = result + qt_monomial(col * power as usize, row * power as usize);
    }
    result
}

fn plethystic_evaluate_from_power_sums<F>(
    f: &SymmetricFunction<QtPolynomial>,
    mut power_sum_value: F,
) -> QtPolynomial
where
    F: FnMut(u32) -> QtPolynomial,
{
    let in_power_sum_basis = f.to_power_sum_basis();
    let mut result = <QtPolynomial as Ring>::zero();
    for (partition, coefficient) in in_power_sum_basis.terms() {
        let mut term = coefficient.clone();
        for &part in partition.parts() {
            term = term * power_sum_value(part);
        }
        result = result + term;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Basis;

    fn p(parts: &[u32]) -> Partition {
        Partition::new(parts.to_vec())
    }

    fn r(value: i64) -> Rational {
        Rational::from_integer(value)
    }

    #[test]
    fn test_qt_monomial_coefficients() {
        let monomial = qt_monomial(2, 3);

        assert_eq!(qt_coefficient(&monomial, 2, 3), r(1));
        assert_eq!(qt_coefficient(&monomial, 1, 3), r(0));
        assert_eq!(qt_coefficient(&monomial, 2, 2), r(0));
    }

    #[test]
    fn test_macdonald_b_and_nabla_eigenvalues() {
        let lambda = p(&[2, 1]);
        let b = macdonald_b_eigenvalue(&lambda);

        assert_eq!(qt_coefficient(&b, 0, 0), r(1));
        assert_eq!(qt_coefficient(&b, 1, 0), r(1));
        assert_eq!(qt_coefficient(&b, 0, 1), r(1));
        assert_eq!(b, qt_constant(1) + qt_monomial(1, 0) + qt_monomial(0, 1));
        assert_eq!(nabla_eigenvalue(&lambda), qt_monomial(1, 1));
    }

    #[test]
    fn test_delta_eigenvalues() {
        let lambda = p(&[2, 1]);
        let e1 = SymmetricFunction::<QtPolynomial>::elementary_symmetric(p(&[1]));

        assert_eq!(
            delta_eigenvalue(&e1, &lambda),
            macdonald_b_eigenvalue(&lambda)
        );
        assert_eq!(
            delta_prime_eigenvalue(&e1, &lambda),
            qt_monomial(1, 0) + qt_monomial(0, 1)
        );
    }

    #[test]
    fn test_delta_eigenvalue_uses_power_sum_conversion() {
        let lambda = p(&[2]);
        let e2 = SymmetricFunction::<QtPolynomial>::elementary_symmetric(p(&[2]));

        assert_eq!(delta_eigenvalue(&e2, &lambda), qt_monomial(1, 0));
    }

    #[test]
    fn test_modified_macdonald_expansion_operators() {
        let lambda = p(&[2, 1]);
        let expansion = ModifiedMacdonaldExpansion::basis_element(lambda.clone());
        let e1 = SymmetricFunction::<QtPolynomial>::basis_element(Basis::Elementary, p(&[1]));

        assert_eq!(
            expansion.nabla().coefficient(&lambda),
            nabla_eigenvalue(&lambda)
        );
        assert_eq!(
            expansion.delta(&e1).coefficient(&lambda),
            macdonald_b_eigenvalue(&lambda)
        );
        assert_eq!(
            expansion.delta_prime(&e1).coefficient(&lambda),
            qt_monomial(1, 0) + qt_monomial(0, 1)
        );
    }
}
