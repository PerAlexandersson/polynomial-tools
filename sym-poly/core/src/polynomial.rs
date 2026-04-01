use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

use crate::Ring;

/// A sparse-by-trailing-zeros univariate polynomial in one indeterminate.
///
/// Coefficients are stored in ascending degree order.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct UnivariatePolynomial<C: Ring> {
    coeffs: Vec<C>,
}

impl<C: Ring> UnivariatePolynomial<C> {
    /// Construct a polynomial from coefficients in ascending degree order.
    pub fn new(coeffs: Vec<C>) -> Self {
        let mut poly = UnivariatePolynomial { coeffs };
        poly.strip_zeros();
        poly
    }

    /// The zero polynomial.
    pub fn zero() -> Self {
        UnivariatePolynomial { coeffs: Vec::new() }
    }

    /// A constant polynomial.
    pub fn constant(c: C) -> Self {
        if c.is_zero() {
            Self::zero()
        } else {
            Self::new(vec![c])
        }
    }

    /// The monomial `coeff * q^degree`.
    pub fn monomial(degree: usize, coeff: C) -> Self {
        if coeff.is_zero() {
            return Self::zero();
        }
        let mut coeffs = vec![C::zero(); degree + 1];
        coeffs[degree] = coeff;
        Self::new(coeffs)
    }

    /// The indeterminate `q`.
    pub fn variable() -> Self {
        Self::monomial(1, C::one())
    }

    /// The coefficients in ascending degree order.
    pub fn coeffs(&self) -> &[C] {
        &self.coeffs
    }

    /// Coefficient of `q^degree`.
    pub fn coeff(&self, degree: usize) -> C {
        self.coeffs.get(degree).cloned().unwrap_or_else(C::zero)
    }

    /// Degree of the polynomial, or `None` for zero.
    pub fn degree(&self) -> Option<usize> {
        if self.coeffs.is_empty() {
            None
        } else {
            Some(self.coeffs.len() - 1)
        }
    }

    /// Evaluate at a ring element.
    pub fn evaluate(&self, x: &C) -> C {
        let mut result = C::zero();
        for coeff in self.coeffs.iter().rev() {
            result = result * x.clone() + coeff.clone();
        }
        result
    }

    fn strip_zeros(&mut self) {
        while self.coeffs.last().is_some_and(|c| c.is_zero()) {
            self.coeffs.pop();
        }
    }
}

impl<C: Ring> Add for UnivariatePolynomial<C> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let n = self.coeffs.len().max(rhs.coeffs.len());
        let coeffs = (0..n).map(|i| self.coeff(i) + rhs.coeff(i)).collect();
        Self::new(coeffs)
    }
}

impl<C: Ring> Sub for UnivariatePolynomial<C> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let n = self.coeffs.len().max(rhs.coeffs.len());
        let coeffs = (0..n).map(|i| self.coeff(i) - rhs.coeff(i)).collect();
        Self::new(coeffs)
    }
}

impl<C: Ring> Neg for UnivariatePolynomial<C> {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(self.coeffs.into_iter().map(|c| -c).collect())
    }
}

impl<C: Ring> Mul for UnivariatePolynomial<C> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        if self.coeffs.is_empty() || rhs.coeffs.is_empty() {
            return Self::zero();
        }
        let mut coeffs = vec![C::zero(); self.coeffs.len() + rhs.coeffs.len() - 1];
        for (i, a) in self.coeffs.into_iter().enumerate() {
            for (j, b) in rhs.coeffs.iter().cloned().enumerate() {
                coeffs[i + j] = coeffs[i + j].clone() + a.clone() * b;
            }
        }
        Self::new(coeffs)
    }
}

impl<C: Ring> Ring for UnivariatePolynomial<C> {
    fn zero() -> Self {
        Self::zero()
    }

    fn one() -> Self {
        Self::constant(C::one())
    }

    fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    fn from_i64(n: i64) -> Self {
        Self::constant(C::from_i64(n))
    }

    fn exact_div_i64(&self, divisor: i64) -> Self {
        Self::new(
            self.coeffs
                .iter()
                .map(|c| c.exact_div_i64(divisor))
                .collect(),
        )
    }
}

impl<C: Ring> fmt::Display for UnivariatePolynomial<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.coeffs.is_empty() {
            return write!(f, "0");
        }

        let mut first = true;
        for degree in (0..self.coeffs.len()).rev() {
            let coeff = &self.coeffs[degree];
            if coeff.is_zero() {
                continue;
            }

            if !first {
                write!(f, " + ")?;
            }
            first = false;

            match degree {
                0 => write!(f, "{coeff}")?,
                1 => {
                    if *coeff == C::one() {
                        write!(f, "q")?
                    } else if *coeff == C::minus_one() {
                        write!(f, "-q")?
                    } else {
                        write!(f, "{coeff}q")?
                    }
                }
                _ => {
                    if *coeff == C::one() {
                        write!(f, "q^{degree}")?
                    } else if *coeff == C::minus_one() {
                        write!(f, "-q^{degree}")?
                    } else {
                        write!(f, "{coeff}q^{degree}")?
                    }
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_mul_and_display() {
        let a = UnivariatePolynomial::<i64>::new(vec![1, 2]);
        let b = UnivariatePolynomial::<i64>::new(vec![0, 1, 1]);
        assert_eq!((a.clone() + b.clone()).coeffs(), &[1, 3, 1]);
        assert_eq!((a * b).coeffs(), &[0, 1, 3, 2]);
        assert_eq!(
            format!("{}", UnivariatePolynomial::<i64>::new(vec![1, 1, 2])),
            "2q^2 + q + 1"
        );
    }

    #[test]
    fn test_ring_exact_div() {
        let p = UnivariatePolynomial::<i64>::new(vec![2, 4, 6]);
        assert_eq!(p.exact_div_i64(2), UnivariatePolynomial::new(vec![1, 2, 3]));
    }

    #[test]
    fn test_evaluate() {
        let p = UnivariatePolynomial::<i64>::new(vec![1, 2, 3]);
        assert_eq!(p.evaluate(&2), 17);
    }
}
