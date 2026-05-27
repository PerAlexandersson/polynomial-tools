use std::ops::Div;

use num_bigint::BigInt;
use num_rational::Ratio;

use crate::Ring;

/// A coefficient type that forms a field.
///
/// This is intentionally separate from [`Ring`]: many existing combinatorial
/// algorithms use integer coefficients, but row reduction and quotient-space
/// coordinates need division by arbitrary nonzero pivots.
pub trait Field: Ring + Div<Output = Self> {}

impl Field for Ratio<i64> {}

impl Field for Ratio<BigInt> {}
