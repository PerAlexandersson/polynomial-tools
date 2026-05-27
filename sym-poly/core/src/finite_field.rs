use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::{Field, Ring};

/// An element of the prime field `F_P`.
///
/// The modulus is a const generic so the existing `Ring`/`Field` traits can
/// construct `zero()` and `one()` without carrying runtime coefficient-ring
/// state. The implementation assumes `P` is prime; division by a nonzero
/// element panics only if this assumption is violated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PrimeField<const P: u64> {
    value: u64,
}

impl<const P: u64> PrimeField<P> {
    pub fn new<T: Into<i128>>(value: T) -> Self {
        assert!(P > 1, "field modulus must be greater than 1");
        let modulus = P as i128;
        let mut reduced = value.into() % modulus;
        if reduced < 0 {
            reduced += modulus;
        }
        Self {
            value: reduced as u64,
        }
    }

    pub fn modulus() -> u64 {
        P
    }

    pub fn value(self) -> u64 {
        self.value
    }

    pub fn balanced_value(self) -> i64 {
        if self.value <= P / 2 {
            self.value as i64
        } else {
            -((P - self.value) as i64)
        }
    }

    pub fn inverse(self) -> Self {
        assert!(!self.is_zero(), "division by zero");
        let (gcd, inverse, _) = extended_gcd(self.value as i128, P as i128);
        assert!(
            gcd == 1,
            "nonzero element is not invertible; modulus is not prime"
        );
        Self::new(inverse)
    }
}

impl<const P: u64> Add for PrimeField<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let value = ((self.value as u128 + rhs.value as u128) % P as u128) as u64;
        Self { value }
    }
}

impl<const P: u64> Sub for PrimeField<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let value = ((P as u128 + self.value as u128 - rhs.value as u128) % P as u128) as u64;
        Self { value }
    }
}

impl<const P: u64> Neg for PrimeField<P> {
    type Output = Self;

    fn neg(self) -> Self {
        if self.is_zero() {
            self
        } else {
            Self {
                value: P - self.value,
            }
        }
    }
}

impl<const P: u64> Mul for PrimeField<P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let value = ((self.value as u128 * rhs.value as u128) % P as u128) as u64;
        Self { value }
    }
}

impl<const P: u64> Div for PrimeField<P> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        self * rhs.inverse()
    }
}

impl<const P: u64> Ring for PrimeField<P> {
    fn zero() -> Self {
        Self { value: 0 }
    }

    fn one() -> Self {
        Self::new(1)
    }

    fn is_zero(&self) -> bool {
        self.value == 0
    }

    fn from_i64(n: i64) -> Self {
        Self::new(n as i128)
    }

    fn exact_div_i64(&self, divisor: i64) -> Self {
        assert!(divisor != 0, "division by zero");
        *self / Self::from_i64(divisor)
    }
}

impl<const P: u64> Field for PrimeField<P> {}

impl<const P: u64> fmt::Display for PrimeField<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.balanced_value())
    }
}

fn extended_gcd(a: i128, b: i128) -> (i128, i128, i128) {
    let (mut old_r, mut r) = (a, b);
    let (mut old_s, mut s) = (1, 0);
    let (mut old_t, mut t) = (0, 1);

    while r != 0 {
        let quotient = old_r / r;
        (old_r, r) = (r, old_r - quotient * r);
        (old_s, s) = (s, old_s - quotient * s);
        (old_t, t) = (t, old_t - quotient * t);
    }

    if old_r < 0 {
        (-old_r, -old_s, -old_t)
    } else {
        (old_r, old_s, old_t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type F101 = PrimeField<101>;

    #[test]
    fn test_prime_field_arithmetic() {
        assert_eq!((F101::from_i64(100) + F101::from_i64(2)).value(), 1);
        assert_eq!((F101::from_i64(3) - F101::from_i64(5)).value(), 99);
        assert_eq!((-F101::from_i64(7)).value(), 94);
        assert_eq!((F101::from_i64(12) * F101::from_i64(9)).value(), 7);
    }

    #[test]
    fn test_prime_field_division() {
        let x = F101::from_i64(37);

        assert_eq!(x / x, F101::one());
        assert_eq!(F101::from_i64(6).exact_div_i64(3), F101::from_i64(2));
    }

    #[test]
    fn test_balanced_value() {
        assert_eq!(F101::from_i64(50).balanced_value(), 50);
        assert_eq!(F101::from_i64(51).balanced_value(), -50);
        assert_eq!(F101::from_i64(-3).balanced_value(), -3);
    }
}
