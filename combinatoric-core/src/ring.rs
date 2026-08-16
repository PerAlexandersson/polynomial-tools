use num_bigint::BigInt;
use num_traits::ToPrimitive;
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

/// A commutative ring with identity.
///
/// Types implementing `Ring` can serve as coefficients for symmetric functions.
/// Standard choices: `i64` for exploration, `BigInt` for exact computation,
/// `Polynomial<BigInt>` for Hall-Littlewood / Jack (future).
///
/// TODO: For combinatorial counting APIs in `combinatoric-core`, prefer the
/// uniform pattern "exact implementation first, narrow wrapper second":
/// `foo_bigint(...) -> Vec<BigInt>` as the source of truth, and
/// `foo(...) -> Vec<i64>` only as a checked convenience projection.
pub trait Ring:
    Clone
    + Eq
    + fmt::Debug
    + fmt::Display
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + Sized
{
    /// The additive identity (0).
    fn zero() -> Self;
    /// The multiplicative identity (1).
    fn one() -> Self;
    /// Whether this element is the additive identity.
    fn is_zero(&self) -> bool;
    /// Embed a 64-bit integer into this ring.
    fn from_i64(n: i64) -> Self;

    /// Embed an arbitrary-precision integer into this ring.
    ///
    /// Bounded coefficient rings use the checked default conversion. Exact
    /// arbitrary-precision rings override this method.
    fn from_bigint(n: &BigInt) -> Self {
        Self::from_i64(
            n.to_i64()
                .expect("integer coefficient does not fit in this coefficient ring"),
        )
    }

    /// Additive inverse exists by `Neg`, but this is convenient.
    fn minus_one() -> Self {
        Self::zero() - Self::one()
    }

    /// Exact division by an integer. Required for power-sum basis conversions
    /// which involve 1/z_μ factors.
    ///
    /// For integer types (i64, BigInt): panics if division is not exact.
    /// For rational types: always works.
    fn exact_div_i64(&self, divisor: i64) -> Self;

    /// Exact division by an arbitrary-precision integer.
    fn exact_div_bigint(&self, divisor: &BigInt) -> Self {
        self.exact_div_i64(
            divisor
                .to_i64()
                .expect("integer divisor does not fit in this coefficient ring"),
        )
    }
}

// ---------------------------------------------------------------------------
// i64
// ---------------------------------------------------------------------------

impl Ring for i64 {
    fn zero() -> Self {
        0
    }
    fn one() -> Self {
        1
    }
    fn is_zero(&self) -> bool {
        *self == 0
    }
    fn from_i64(n: i64) -> Self {
        n
    }
    fn exact_div_i64(&self, divisor: i64) -> Self {
        assert!(divisor != 0, "division by zero");
        assert!(
            self % divisor == 0,
            "inexact division: {} / {} (power-sum conversions require rational coefficients)",
            self,
            divisor
        );
        self / divisor
    }
}

// ---------------------------------------------------------------------------
// BigInt
// ---------------------------------------------------------------------------

impl Ring for BigInt {
    fn zero() -> Self {
        BigInt::from(0)
    }
    fn one() -> Self {
        BigInt::from(1)
    }
    fn is_zero(&self) -> bool {
        *self == BigInt::from(0)
    }
    fn from_i64(n: i64) -> Self {
        BigInt::from(n)
    }
    fn from_bigint(n: &BigInt) -> Self {
        n.clone()
    }
    fn exact_div_i64(&self, divisor: i64) -> Self {
        let d = BigInt::from(divisor);
        let (q, r) = num_integer::Integer::div_rem(self, &d);
        assert!(
            r == BigInt::from(0),
            "inexact division: {} / {} (power-sum conversions require rational coefficients)",
            self,
            divisor
        );
        q
    }
    fn exact_div_bigint(&self, divisor: &BigInt) -> Self {
        let (q, r) = num_integer::Integer::div_rem(self, divisor);
        assert!(r == BigInt::from(0), "inexact arbitrary-precision division");
        q
    }
}

// ---------------------------------------------------------------------------
// Rational<BigInt>
// ---------------------------------------------------------------------------

use num_rational::Ratio;

impl Ring for Ratio<BigInt> {
    fn zero() -> Self {
        Ratio::new(BigInt::from(0), BigInt::from(1))
    }
    fn one() -> Self {
        Ratio::new(BigInt::from(1), BigInt::from(1))
    }
    fn is_zero(&self) -> bool {
        self.numer().is_zero()
    }
    fn from_i64(n: i64) -> Self {
        Ratio::from_integer(BigInt::from(n))
    }
    fn from_bigint(n: &BigInt) -> Self {
        Ratio::from_integer(n.clone())
    }
    fn exact_div_i64(&self, divisor: i64) -> Self {
        self.clone() / Ratio::from_integer(BigInt::from(divisor))
    }
    fn exact_div_bigint(&self, divisor: &BigInt) -> Self {
        self.clone() / Ratio::from_integer(divisor.clone())
    }
}

impl Ring for Ratio<i64> {
    fn zero() -> Self {
        Ratio::new(0, 1)
    }
    fn one() -> Self {
        Ratio::new(1, 1)
    }
    fn is_zero(&self) -> bool {
        *self.numer() == 0
    }
    fn from_i64(n: i64) -> Self {
        Ratio::from_integer(n)
    }
    fn exact_div_i64(&self, divisor: i64) -> Self {
        self.clone() / Ratio::from_integer(divisor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_i64_ring() {
        assert_eq!(i64::zero() + i64::one(), 1);
        assert_eq!(i64::from_i64(-3) * i64::from_i64(7), -21);
        assert!(i64::zero().is_zero());
        assert!(!i64::one().is_zero());
    }

    #[test]
    fn test_bigint_ring() {
        let a = BigInt::from_i64(100);
        let b = BigInt::from_i64(200);
        assert_eq!(a + b, BigInt::from(300));
    }
}
