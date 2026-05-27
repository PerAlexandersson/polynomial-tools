use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::{Field, Ring, UnivariatePolynomial};

/// Formal fractions over a coefficient ring.
///
/// This is intended as a lightweight fraction-field layer for coefficient
/// rings such as `Q[q,t]`. Fractions are compared by cross multiplication; no
/// polynomial gcd reduction is attempted.
#[derive(Debug, Clone)]
pub struct RationalFunction<P: Ring> {
    numerator: P,
    denominator: P,
}

impl<P: Ring> RationalFunction<P> {
    pub fn new(numerator: P, denominator: P) -> Self {
        assert!(!denominator.is_zero(), "zero denominator");
        if numerator.is_zero() {
            return Self {
                numerator,
                denominator: P::one(),
            };
        }
        Self {
            numerator,
            denominator,
        }
    }

    pub fn from_polynomial(polynomial: P) -> Self {
        Self::new(polynomial, P::one())
    }

    pub fn numerator(&self) -> &P {
        &self.numerator
    }

    pub fn denominator(&self) -> &P {
        &self.denominator
    }

    pub fn inverse(self) -> Self {
        assert!(!self.numerator.is_zero(), "division by zero");
        Self::new(self.denominator, self.numerator)
    }
}

impl<P: Ring> Add for RationalFunction<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let numerator = self.numerator.clone() * rhs.denominator.clone()
            + rhs.numerator.clone() * self.denominator.clone();
        let denominator = self.denominator * rhs.denominator;
        Self::new(numerator, denominator)
    }
}

impl<P: Ring> Sub for RationalFunction<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}

impl<P: Ring> Neg for RationalFunction<P> {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(-self.numerator, self.denominator)
    }
}

impl<P: Ring> Mul for RationalFunction<P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self::new(
            self.numerator * rhs.numerator,
            self.denominator * rhs.denominator,
        )
    }
}

impl<P: Ring> Div for RationalFunction<P> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        assert!(!rhs.numerator.is_zero(), "division by zero");
        Self::new(
            self.numerator * rhs.denominator,
            self.denominator * rhs.numerator,
        )
    }
}

impl<P: Ring> PartialEq for RationalFunction<P> {
    fn eq(&self, other: &Self) -> bool {
        self.numerator.clone() * other.denominator.clone()
            == other.numerator.clone() * self.denominator.clone()
    }
}

impl<P: Ring> Eq for RationalFunction<P> {}

impl<P: Ring> Ring for RationalFunction<P> {
    fn zero() -> Self {
        Self::from_polynomial(P::zero())
    }

    fn one() -> Self {
        Self::from_polynomial(P::one())
    }

    fn is_zero(&self) -> bool {
        self.numerator.is_zero()
    }

    fn from_i64(n: i64) -> Self {
        Self::from_polynomial(P::from_i64(n))
    }

    fn exact_div_i64(&self, divisor: i64) -> Self {
        assert!(divisor != 0, "division by zero");
        Self::new(
            self.numerator.clone(),
            self.denominator.clone() * P::from_i64(divisor),
        )
    }
}

impl<P: Ring> Field for RationalFunction<P> {}

impl<P: Ring> fmt::Display for RationalFunction<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.denominator == P::one() {
            write!(f, "{}", self.numerator)
        } else {
            write!(f, "({})/({})", self.numerator, self.denominator)
        }
    }
}

/// Polynomials in `q,t`, represented as polynomials in `t` whose coefficients
/// are polynomials in `q`.
pub type QtPolynomial<C> = UnivariatePolynomial<UnivariatePolynomial<C>>;

/// Rational functions in `q,t`.
pub type QtRationalFunction<C> = RationalFunction<QtPolynomial<C>>;

pub fn qt_polynomial_constant<C: Ring>(coefficient: C) -> QtPolynomial<C> {
    UnivariatePolynomial::constant(UnivariatePolynomial::constant(coefficient))
}

pub fn qt_constant<C: Ring>(value: i64) -> QtPolynomial<C> {
    qt_polynomial_constant(C::from_i64(value))
}

pub fn qt_monomial<C: Ring>(q_degree: usize, t_degree: usize, coefficient: C) -> QtPolynomial<C> {
    let q_part = UnivariatePolynomial::monomial(q_degree, coefficient);
    UnivariatePolynomial::monomial(t_degree, q_part)
}

pub fn qt_unit_monomial<C: Ring>(q_degree: usize, t_degree: usize) -> QtPolynomial<C> {
    qt_monomial(q_degree, t_degree, C::one())
}

pub fn qt_coefficient<C: Ring>(
    polynomial: &QtPolynomial<C>,
    q_degree: usize,
    t_degree: usize,
) -> C {
    polynomial.coeff(t_degree).coeff(q_degree)
}

pub fn qt_rational_monomial<C: Ring>(q_degree: usize, t_degree: usize) -> QtRationalFunction<C> {
    QtRationalFunction::from_polynomial(qt_unit_monomial(q_degree, t_degree))
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Ratio;

    type Q = Ratio<i64>;
    type Qt = QtPolynomial<Q>;
    type K = QtRationalFunction<Q>;

    fn q(value: i64) -> Q {
        Q::from_integer(value)
    }

    #[test]
    fn test_qt_polynomial_coefficients() {
        let polynomial = qt_monomial(2, 3, q(5)) + qt_unit_monomial(1, 0);

        assert_eq!(qt_coefficient(&polynomial, 2, 3), q(5));
        assert_eq!(qt_coefficient(&polynomial, 1, 0), q(1));
        assert_eq!(qt_coefficient(&polynomial, 3, 2), q(0));
    }

    #[test]
    fn test_rational_function_cross_multiply_equality() {
        let q_var = qt_unit_monomial::<Q>(1, 0);
        let quotient = RationalFunction::new(q_var.clone(), q_var);

        assert_eq!(quotient, K::one());
    }

    #[test]
    fn test_qt_rational_function_arithmetic() {
        let q_var = K::from_polynomial(qt_unit_monomial::<Q>(1, 0));
        let t_var = K::from_polynomial(qt_unit_monomial::<Q>(0, 1));
        let fraction = q_var.clone() / t_var.clone();

        assert_eq!(fraction * t_var, q_var);
        assert_eq!(
            K::from_polynomial(Qt::from_i64(2)).exact_div_i64(2),
            K::one()
        );
    }
}
