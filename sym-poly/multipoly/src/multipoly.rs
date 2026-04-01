//! Sparse multivariate polynomial type.
//!
//! `MultiPoly<C>` represents polynomials in x_1, ..., x_n with coefficients
//! in a ring C. Monomials are stored as exponent vectors in a BTreeMap for
//! efficient sparse representation with lexicographic ordering.

use std::collections::BTreeMap;
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

use sym_poly_core::Ring;

/// A sparse multivariate polynomial in x_1, ..., x_{num_vars}.
///
/// Monomials are represented as `Vec<u32>` exponent vectors, stored
/// in a `BTreeMap` for lexicographic ordering and O(log n) access.
#[derive(Debug, Clone)]
pub struct MultiPoly<C: Ring> {
    num_vars: usize,
    pub(crate) terms: BTreeMap<Vec<u32>, C>,
}

impl<C: Ring> MultiPoly<C> {
    /// Create the zero polynomial in n variables.
    pub fn zero(num_vars: usize) -> Self {
        MultiPoly {
            num_vars,
            terms: BTreeMap::new(),
        }
    }

    /// Create from a map of exponent vectors to coefficients. Strips zeros.
    pub fn from_terms(num_vars: usize, terms: BTreeMap<Vec<u32>, C>) -> Self {
        let mut p = MultiPoly { num_vars, terms };
        p.strip_zeros();
        p
    }

    /// Create a monomial c * x_1^{e_1} * ... * x_n^{e_n}.
    pub fn monomial(num_vars: usize, exponents: Vec<u32>, coeff: C) -> Self {
        if coeff.is_zero() {
            return Self::zero(num_vars);
        }
        assert_eq!(exponents.len(), num_vars);
        let mut terms = BTreeMap::new();
        terms.insert(exponents, coeff);
        MultiPoly { num_vars, terms }
    }

    /// Create the variable x_i (0-indexed).
    pub fn var(num_vars: usize, i: usize) -> Self {
        assert!(i < num_vars);
        let mut exp = vec![0u32; num_vars];
        exp[i] = 1;
        Self::monomial(num_vars, exp, C::one())
    }

    /// Create the monomial x^α = x_1^{α_1} * ... * x_n^{α_n}.
    pub fn x_power(num_vars: usize, exponents: Vec<u32>) -> Self {
        assert_eq!(exponents.len(), num_vars);
        Self::monomial(num_vars, exponents, C::one())
    }

    /// Create a constant polynomial.
    pub fn constant(num_vars: usize, c: C) -> Self {
        if c.is_zero() {
            return Self::zero(num_vars);
        }
        Self::monomial(num_vars, vec![0; num_vars], c)
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    pub fn num_vars(&self) -> usize {
        self.num_vars
    }
    pub fn terms(&self) -> &BTreeMap<Vec<u32>, C> {
        &self.terms
    }
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn coefficient(&self, exponents: &[u32]) -> C {
        self.terms.get(exponents).cloned().unwrap_or_else(C::zero)
    }

    /// Total degree of the polynomial.
    pub fn total_degree(&self) -> Option<u32> {
        self.terms.keys().map(|e| e.iter().sum::<u32>()).max()
    }

    pub(crate) fn strip_zeros(&mut self) {
        self.terms.retain(|_, c| !c.is_zero());
    }

    // -----------------------------------------------------------------------
    // Variable swap: s_i swaps x_i and x_{i+1}
    // -----------------------------------------------------------------------

    /// Apply the simple transposition s_i: swap x_i and x_{i+1} (0-indexed).
    pub fn swap_vars(&self, i: usize) -> Self {
        assert!(i + 1 < self.num_vars, "swap index out of bounds");
        let terms = self
            .terms
            .iter()
            .map(|(exp, c)| {
                let mut new_exp = exp.clone();
                new_exp.swap(i, i + 1);
                (new_exp, c.clone())
            })
            .collect();
        MultiPoly {
            num_vars: self.num_vars,
            terms,
        }
    }

    /// Multiply by x_i (0-indexed).
    pub fn mul_var(&self, i: usize) -> Self {
        assert!(i < self.num_vars);
        let terms = self
            .terms
            .iter()
            .map(|(exp, c)| {
                let mut new_exp = exp.clone();
                new_exp[i] += 1;
                (new_exp, c.clone())
            })
            .collect();
        MultiPoly {
            num_vars: self.num_vars,
            terms,
        }
    }

    pub fn scale(&self, scalar: &C) -> Self {
        if scalar.is_zero() {
            return Self::zero(self.num_vars);
        }
        let terms = self
            .terms
            .iter()
            .map(|(e, c)| (e.clone(), c.clone() * scalar.clone()))
            .collect();
        Self::from_terms(self.num_vars, terms)
    }
}

// ---------------------------------------------------------------------------
// Arithmetic
// ---------------------------------------------------------------------------

impl<C: Ring> Add for MultiPoly<C> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        assert_eq!(self.num_vars, rhs.num_vars);
        let mut terms = self.terms;
        for (e, c) in rhs.terms {
            let entry = terms.entry(e).or_insert_with(C::zero);
            *entry = entry.clone() + c;
        }
        Self::from_terms(self.num_vars, terms)
    }
}

impl<C: Ring> Sub for MultiPoly<C> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}

impl<C: Ring> Neg for MultiPoly<C> {
    type Output = Self;
    fn neg(self) -> Self {
        let terms = self.terms.into_iter().map(|(e, c)| (e, -c)).collect();
        MultiPoly {
            num_vars: self.num_vars,
            terms,
        }
    }
}

impl<C: Ring> Mul for MultiPoly<C> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        assert_eq!(self.num_vars, rhs.num_vars);
        let n = self.num_vars;
        let mut terms: BTreeMap<Vec<u32>, C> = BTreeMap::new();
        for (e1, c1) in &self.terms {
            if c1.is_zero() {
                continue;
            }
            for (e2, c2) in &rhs.terms {
                let exp: Vec<u32> = (0..n).map(|i| e1[i] + e2[i]).collect();
                let coeff = c1.clone() * c2.clone();
                let entry = terms.entry(exp).or_insert_with(C::zero);
                *entry = entry.clone() + coeff;
            }
        }
        Self::from_terms(n, terms)
    }
}

impl<C: Ring> PartialEq for MultiPoly<C> {
    fn eq(&self, other: &Self) -> bool {
        self.num_vars == other.num_vars && self.terms == other.terms
    }
}
impl<C: Ring> Eq for MultiPoly<C> {}

impl<C: Ring> fmt::Display for MultiPoly<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let mut first = true;
        for (exp, coeff) in &self.terms {
            if !first {
                write!(f, " + ")?;
            }
            first = false;
            // Format monomial
            let mut mono = String::new();
            for (i, &e) in exp.iter().enumerate() {
                if e == 0 {
                    continue;
                }
                if !mono.is_empty() {
                    mono.push('*');
                }
                if e == 1 {
                    mono.push_str(&format!("x{}", i + 1));
                } else {
                    mono.push_str(&format!("x{}^{}", i + 1, e));
                }
            }
            if mono.is_empty() {
                write!(f, "{}", coeff)?;
            } else if *coeff == C::one() {
                write!(f, "{}", mono)?;
            } else if *coeff == C::minus_one() {
                write!(f, "-{}", mono)?;
            } else {
                write!(f, "{}*{}", coeff, mono)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero() {
        let z: MultiPoly<i64> = MultiPoly::zero(3);
        assert!(z.is_zero());
    }

    #[test]
    fn test_var() {
        let x1: MultiPoly<i64> = MultiPoly::var(3, 0);
        assert_eq!(x1.coefficient(&[1, 0, 0]), 1);
        assert_eq!(x1.coefficient(&[0, 1, 0]), 0);
    }

    #[test]
    fn test_add() {
        let x1: MultiPoly<i64> = MultiPoly::var(2, 0);
        let x2: MultiPoly<i64> = MultiPoly::var(2, 1);
        let sum = x1 + x2;
        assert_eq!(sum.coefficient(&[1, 0]), 1);
        assert_eq!(sum.coefficient(&[0, 1]), 1);
    }

    #[test]
    fn test_mul() {
        let x1: MultiPoly<i64> = MultiPoly::var(2, 0);
        let x2: MultiPoly<i64> = MultiPoly::var(2, 1);
        let prod = x1 * x2;
        assert_eq!(prod.coefficient(&[1, 1]), 1);
        assert_eq!(prod.terms().len(), 1);
    }

    #[test]
    fn test_swap_vars() {
        // x1^2 * x2 under s_0 becomes x1 * x2^2
        let f: MultiPoly<i64> = MultiPoly::x_power(3, vec![2, 1, 0]);
        let swapped = f.swap_vars(0);
        assert_eq!(swapped.coefficient(&[1, 2, 0]), 1);
    }

    #[test]
    fn test_mul_var() {
        let f: MultiPoly<i64> = MultiPoly::x_power(2, vec![1, 0]);
        let g = f.mul_var(1);
        assert_eq!(g.coefficient(&[1, 1]), 1);
    }

    #[test]
    fn test_constant() {
        let c: MultiPoly<i64> = MultiPoly::constant(2, 5);
        assert_eq!(c.coefficient(&[0, 0]), 5);
    }
}
