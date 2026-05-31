//! A polynomial expressed in one of the nonsymmetric polynomial bases.
//!
//! `MultiPolyFunction<C>` is the nonsymmetric analogue of `SymmetricFunction<C>`:
//! it stores a linear combination of basis elements (key, atom, slide, or monomial)
//! indexed by weak compositions, and supports conversion between bases.

use std::collections::BTreeMap;
use std::fmt;
use std::ops::{Add, Neg, Sub};

use combinatoric_core::WeakComposition;
use sym_poly_core::Ring;

use crate::basis::MultiPolyBasis;
use crate::multipoly::MultiPoly;
use crate::transition;

/// A polynomial expressed as a formal linear combination of basis elements
/// indexed by weak compositions.
///
/// The `num_vars` field determines the ambient polynomial ring (number of variables).
#[derive(Debug, Clone)]
pub struct MultiPolyFunction<C: Ring> {
    basis: MultiPolyBasis,
    num_vars: usize,
    terms: BTreeMap<WeakComposition, C>,
}

impl<C: Ring> MultiPolyFunction<C> {
    // -------------------------------------------------------------------
    // Constructors
    // -------------------------------------------------------------------

    /// Create from explicit terms. Strips zero coefficients.
    pub fn from_terms(
        basis: MultiPolyBasis,
        num_vars: usize,
        terms: BTreeMap<WeakComposition, C>,
    ) -> Self {
        assert!(
            terms.keys().all(|alpha| alpha.num_vars() == num_vars),
            "all weak compositions must have length num_vars"
        );
        let mut f = MultiPolyFunction {
            basis,
            num_vars,
            terms,
        };
        f.strip_zeros();
        f
    }

    /// The zero element.
    pub fn zero(basis: MultiPolyBasis, num_vars: usize) -> Self {
        MultiPolyFunction {
            basis,
            num_vars,
            terms: BTreeMap::new(),
        }
    }

    /// A single basis element with coefficient 1.
    pub fn basis_element(basis: MultiPolyBasis, alpha: &[u32]) -> Self {
        let wc = WeakComposition::from_slice(alpha);
        let mut terms = BTreeMap::new();
        terms.insert(wc, C::one());
        MultiPolyFunction {
            basis,
            num_vars: alpha.len(),
            terms,
        }
    }

    /// A single basis element scaled by a coefficient.
    pub fn scaled_basis_element(basis: MultiPolyBasis, alpha: &[u32], coeff: C) -> Self {
        if coeff.is_zero() {
            return Self::zero(basis, alpha.len());
        }
        let wc = WeakComposition::from_slice(alpha);
        let mut terms = BTreeMap::new();
        terms.insert(wc, coeff);
        MultiPolyFunction {
            basis,
            num_vars: alpha.len(),
            terms,
        }
    }

    /// Key polynomial κ_α as a single basis element.
    pub fn key(alpha: &[u32]) -> Self {
        Self::basis_element(MultiPolyBasis::Key, alpha)
    }

    /// Demazure atom A_α as a single basis element.
    pub fn atom(alpha: &[u32]) -> Self {
        Self::basis_element(MultiPolyBasis::Atom, alpha)
    }

    /// Fundamental slide polynomial F_α as a single basis element.
    pub fn fund_slide(alpha: &[u32]) -> Self {
        Self::basis_element(MultiPolyBasis::FundSlide, alpha)
    }

    /// Monomial slide polynomial M_α as a single basis element.
    pub fn mon_slide(alpha: &[u32]) -> Self {
        Self::basis_element(MultiPolyBasis::MonSlide, alpha)
    }

    /// Monomial x^α as a single basis element.
    pub fn monomial(alpha: &[u32]) -> Self {
        Self::basis_element(MultiPolyBasis::Monomial, alpha)
    }

    // -------------------------------------------------------------------
    // Bridge: MultiPoly <-> MultiPolyFunction
    // -------------------------------------------------------------------

    /// Convert a `MultiPoly<C>` (monomial expansion) into a `MultiPolyFunction`
    /// in the Monomial basis.
    pub fn from_multipoly(p: &MultiPoly<C>) -> Self {
        let num_vars = p.num_vars();
        let terms: BTreeMap<WeakComposition, C> = p
            .terms()
            .iter()
            .map(|(exp, c)| (WeakComposition::from_slice(exp), c.clone()))
            .collect();
        MultiPolyFunction {
            basis: MultiPolyBasis::Monomial,
            num_vars,
            terms,
        }
    }

    /// Expand to the monomial representation as a `MultiPoly<C>`.
    pub fn to_multipoly(&self) -> MultiPoly<C> {
        let mon = self.to_monomial_basis();
        let terms: BTreeMap<Vec<u32>, C> = mon
            .terms
            .into_iter()
            .map(|(wc, c)| (wc.parts().to_vec(), c))
            .collect();
        MultiPoly::from_terms(self.num_vars, terms)
    }

    // -------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------

    /// Which basis these terms are expressed in.
    pub fn basis(&self) -> MultiPolyBasis {
        self.basis
    }

    /// Number of variables.
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// The terms map (weak composition → coefficient).
    pub fn terms(&self) -> &BTreeMap<WeakComposition, C> {
        &self.terms
    }

    /// Consume self and return the terms map.
    pub fn into_terms(self) -> BTreeMap<WeakComposition, C> {
        self.terms
    }

    /// Coefficient of a given weak composition.
    pub fn coefficient(&self, alpha: &WeakComposition) -> C {
        self.terms.get(alpha).cloned().unwrap_or_else(C::zero)
    }

    /// Whether this is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Number of nonzero terms.
    pub fn num_terms(&self) -> usize {
        self.terms.len()
    }

    /// The maximum degree among nonzero terms, or `None` if zero.
    pub fn degree(&self) -> Option<u32> {
        self.terms.keys().map(|wc| wc.degree()).max()
    }

    /// Whether all coefficients are positive.
    pub fn positive_coefficients(&self) -> bool
    where
        C: Ord,
    {
        self.terms.values().all(|c| *c > C::zero())
    }

    /// Multiply all coefficients by a scalar.
    pub fn scale(&self, scalar: &C) -> Self {
        if scalar.is_zero() {
            return Self::zero(self.basis, self.num_vars);
        }
        let terms = self
            .terms
            .iter()
            .map(|(wc, c)| (wc.clone(), c.clone() * scalar.clone()))
            .filter(|(_, c)| !c.is_zero())
            .collect();
        MultiPolyFunction {
            basis: self.basis,
            num_vars: self.num_vars,
            terms,
        }
    }

    // -------------------------------------------------------------------
    // Basis conversion
    // -------------------------------------------------------------------

    /// Convert to a different basis.
    pub fn to_basis(&self, target: MultiPolyBasis) -> Self {
        transition::convert(self, target)
    }

    /// Convenience: convert to monomial basis.
    pub fn to_monomial_basis(&self) -> Self {
        self.to_basis(MultiPolyBasis::Monomial)
    }

    /// Convenience: convert to key basis.
    pub fn to_key_basis(&self) -> Self {
        self.to_basis(MultiPolyBasis::Key)
    }

    /// Convenience: convert to atom basis.
    pub fn to_atom_basis(&self) -> Self {
        self.to_basis(MultiPolyBasis::Atom)
    }

    /// Convenience: convert to fundamental slide basis.
    pub fn to_fund_slide_basis(&self) -> Self {
        self.to_basis(MultiPolyBasis::FundSlide)
    }

    /// Convenience: convert to monomial slide basis.
    pub fn to_mon_slide_basis(&self) -> Self {
        self.to_basis(MultiPolyBasis::MonSlide)
    }

    // -------------------------------------------------------------------
    // Internal
    // -------------------------------------------------------------------

    fn strip_zeros(&mut self) {
        self.terms.retain(|_, c| !c.is_zero());
    }
}

// =========================================================================
// Arithmetic
// =========================================================================

impl<C: Ring> Add for MultiPolyFunction<C> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        assert_eq!(self.basis, other.basis, "cannot add different bases");
        assert_eq!(
            self.num_vars, other.num_vars,
            "cannot add different num_vars"
        );
        let mut result = self;
        for (wc, c) in other.terms {
            let entry = result.terms.entry(wc).or_insert_with(C::zero);
            *entry = entry.clone() + c;
        }
        result.strip_zeros();
        result
    }
}

impl<C: Ring> Sub for MultiPolyFunction<C> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        assert_eq!(self.basis, other.basis, "cannot subtract different bases");
        assert_eq!(
            self.num_vars, other.num_vars,
            "cannot subtract different num_vars"
        );
        let mut result = self;
        for (wc, c) in other.terms {
            let entry = result.terms.entry(wc).or_insert_with(C::zero);
            *entry = entry.clone() - c;
        }
        result.strip_zeros();
        result
    }
}

impl<C: Ring> Neg for MultiPolyFunction<C> {
    type Output = Self;

    fn neg(self) -> Self {
        let terms = self.terms.into_iter().map(|(wc, c)| (wc, -c)).collect();
        MultiPolyFunction {
            basis: self.basis,
            num_vars: self.num_vars,
            terms,
        }
    }
}

impl<C: Ring> PartialEq for MultiPolyFunction<C> {
    fn eq(&self, other: &Self) -> bool {
        self.basis == other.basis && self.num_vars == other.num_vars && self.terms == other.terms
    }
}

impl<C: Ring> Eq for MultiPolyFunction<C> {}

// =========================================================================
// Display
// =========================================================================

impl<C: Ring> fmt::Display for MultiPolyFunction<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let sym = self.basis.symbol();
        let mut first = true;
        for (wc, coeff) in &self.terms {
            if first {
                write!(f, "({}) {}_{}", coeff, sym, wc)?;
                first = false;
            } else {
                write!(f, " + ({}) {}_{}", coeff, sym, wc)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "all weak compositions must have length num_vars")]
    fn test_from_terms_rejects_wrong_weak_composition_length() {
        let mut terms = BTreeMap::new();
        terms.insert(WeakComposition::from_slice(&[1]), 1i64);

        let _ = MultiPolyFunction::from_terms(MultiPolyBasis::Monomial, 2, terms);
    }
}
