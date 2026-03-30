//! Generic formal linear combination of basis elements.
//!
//! `FormalSum<I, C>` represents Σ c_α B_α where α ranges over indices
//! of type `I` and coefficients live in ring `C`. This is the foundation
//! for `SymmetricFunction`, `QSymFunction`, etc.

use std::collections::BTreeMap;
use std::fmt;
use std::ops::{Add, Neg, Sub};

use crate::index::BasisIndex;
use crate::ring::Ring;

/// A formal linear combination of basis elements indexed by `I` with
/// coefficients in ring `C`.
#[derive(Debug, Clone)]
pub struct FormalSum<I: BasisIndex, C: Ring> {
    terms: BTreeMap<I, C>,
}

impl<I: BasisIndex, C: Ring> FormalSum<I, C> {
    /// Create an empty formal sum (the zero element).
    pub fn zero() -> Self {
        FormalSum {
            terms: BTreeMap::new(),
        }
    }

    /// Create from a map of terms. Strips zero coefficients.
    pub fn from_terms(terms: BTreeMap<I, C>) -> Self {
        let mut fs = FormalSum { terms };
        fs.strip_zeros();
        fs
    }

    /// Create a single basis element with coefficient 1.
    pub fn basis_element(index: I) -> Self {
        let mut terms = BTreeMap::new();
        terms.insert(index, C::one());
        FormalSum { terms }
    }

    /// Create a single basis element scaled by a coefficient.
    pub fn scaled_basis_element(index: I, coeff: C) -> Self {
        if coeff.is_zero() {
            return Self::zero();
        }
        let mut terms = BTreeMap::new();
        terms.insert(index, coeff);
        FormalSum { terms }
    }

    /// The underlying terms map.
    pub fn terms(&self) -> &BTreeMap<I, C> {
        &self.terms
    }

    /// Mutable access to the terms map.
    pub fn terms_mut(&mut self) -> &mut BTreeMap<I, C> {
        &mut self.terms
    }

    /// Consume self and return the terms map.
    pub fn into_terms(self) -> BTreeMap<I, C> {
        self.terms
    }

    /// Coefficient of a given index. Returns zero if absent.
    pub fn coefficient(&self, index: &I) -> C {
        self.terms.get(index).cloned().unwrap_or_else(C::zero)
    }

    /// Whether this is the zero element.
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Number of nonzero terms.
    pub fn num_terms(&self) -> usize {
        self.terms.len()
    }

    /// Remove all terms with zero coefficient.
    pub fn strip_zeros(&mut self) {
        self.terms.retain(|_, c| !c.is_zero());
    }

    /// Multiply all coefficients by a scalar.
    pub fn scale(&self, scalar: &C) -> Self {
        if scalar.is_zero() {
            return Self::zero();
        }
        let terms = self
            .terms
            .iter()
            .map(|(i, c)| (i.clone(), c.clone() * scalar.clone()))
            .filter(|(_, c)| !c.is_zero())
            .collect();
        FormalSum { terms }
    }

    /// Whether all coefficients are positive (requires `C: Ord`).
    pub fn positive_coefficients(&self) -> bool
    where
        C: Ord,
    {
        self.terms.values().all(|c| *c > C::zero())
    }

    /// The maximum degree among nonzero terms, or `None` if zero.
    pub fn degree(&self) -> Option<u32> {
        self.terms.keys().map(|i| i.degree()).max()
    }

    /// Whether all nonzero terms have the same degree.
    pub fn is_homogeneous(&self) -> bool {
        if self.terms.len() <= 1 {
            return true;
        }
        let mut degs = self.terms.keys().map(|i| i.degree());
        let first = degs.next().unwrap();
        degs.all(|d| d == first)
    }

    /// Add a single term, combining with any existing coefficient.
    pub fn add_term(&mut self, index: I, coeff: C) {
        if coeff.is_zero() {
            return;
        }
        let entry = self.terms.entry(index).or_insert_with(C::zero);
        *entry = entry.clone() + coeff;
        // Don't strip here for performance; caller can strip_zeros() when done
    }
}

// ---------------------------------------------------------------------------
// Arithmetic: Add, Sub, Neg
// ---------------------------------------------------------------------------

impl<I: BasisIndex, C: Ring> Add for FormalSum<I, C> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut result = self;
        for (i, c) in other.terms {
            let entry = result.terms.entry(i).or_insert_with(C::zero);
            *entry = entry.clone() + c;
        }
        result.strip_zeros();
        result
    }
}

impl<I: BasisIndex, C: Ring> Sub for FormalSum<I, C> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let mut result = self;
        for (i, c) in other.terms {
            let entry = result.terms.entry(i).or_insert_with(C::zero);
            *entry = entry.clone() - c;
        }
        result.strip_zeros();
        result
    }
}

impl<I: BasisIndex, C: Ring> Neg for FormalSum<I, C> {
    type Output = Self;

    fn neg(self) -> Self {
        let terms = self
            .terms
            .into_iter()
            .map(|(i, c)| (i, -c))
            .collect();
        FormalSum { terms }
    }
}

impl<I: BasisIndex, C: Ring> PartialEq for FormalSum<I, C> {
    fn eq(&self, other: &Self) -> bool {
        self.terms == other.terms
    }
}

impl<I: BasisIndex, C: Ring> Eq for FormalSum<I, C> {}

// ---------------------------------------------------------------------------
// Display
// ---------------------------------------------------------------------------

impl<I: BasisIndex, C: Ring> fmt::Display for FormalSum<I, C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let mut first = true;
        for (index, coeff) in &self.terms {
            if first {
                write!(f, "({}) [{}]", coeff, index)?;
                first = false;
            } else {
                write!(f, " + ({}) [{}]", coeff, index)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition::Partition;

    #[test]
    fn test_zero() {
        let z: FormalSum<Partition, i64> = FormalSum::zero();
        assert!(z.is_zero());
        assert_eq!(z.num_terms(), 0);
    }

    #[test]
    fn test_basis_element() {
        let p = Partition::new(vec![2, 1]);
        let f: FormalSum<Partition, i64> = FormalSum::basis_element(p.clone());
        assert_eq!(f.coefficient(&p), 1);
        assert_eq!(f.num_terms(), 1);
    }

    #[test]
    fn test_add_sub() {
        let p1 = Partition::new(vec![2, 1]);
        let p2 = Partition::new(vec![3]);

        let f1: FormalSum<Partition, i64> = FormalSum::scaled_basis_element(p1.clone(), 3);
        let f2: FormalSum<Partition, i64> = FormalSum::scaled_basis_element(p2.clone(), 5);

        let sum = f1.clone() + f2.clone();
        assert_eq!(sum.coefficient(&p1), 3);
        assert_eq!(sum.coefficient(&p2), 5);

        let diff = sum - f2;
        assert_eq!(diff.coefficient(&p1), 3);
        assert_eq!(diff.coefficient(&p2), 0);
        assert_eq!(diff.num_terms(), 1);
    }

    #[test]
    fn test_neg() {
        let p = Partition::new(vec![2, 1]);
        let f: FormalSum<Partition, i64> = FormalSum::scaled_basis_element(p.clone(), 7);
        let neg = -f;
        assert_eq!(neg.coefficient(&p), -7);
    }

    #[test]
    fn test_scale() {
        let p = Partition::new(vec![2, 1]);
        let f: FormalSum<Partition, i64> = FormalSum::scaled_basis_element(p.clone(), 3);
        let scaled = f.scale(&4);
        assert_eq!(scaled.coefficient(&p), 12);
    }

    #[test]
    fn test_cancel() {
        let p = Partition::new(vec![2, 1]);
        let f1: FormalSum<Partition, i64> = FormalSum::scaled_basis_element(p.clone(), 5);
        let f2: FormalSum<Partition, i64> = FormalSum::scaled_basis_element(p.clone(), -5);
        let sum = f1 + f2;
        assert!(sum.is_zero());
    }

    #[test]
    fn test_display() {
        let z: FormalSum<Partition, i64> = FormalSum::zero();
        assert_eq!(format!("{}", z), "0");
    }
}
