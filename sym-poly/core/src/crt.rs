use num_integer::Integer;
use num_rational::Ratio;

/// Combine two congruences with coprime moduli.
///
/// Returns the unique residue modulo `m1 * m2`, normalized to
/// `0 <= residue < modulus`.
pub fn chinese_remainder_pair(
    residue1: i128,
    modulus1: i128,
    residue2: i128,
    modulus2: i128,
) -> Option<(i128, i128)> {
    assert!(modulus1 > 0, "first modulus must be positive");
    assert!(modulus2 > 0, "second modulus must be positive");
    if modulus1.gcd(&modulus2) != 1 {
        return None;
    }

    let r1 = residue_mod(residue1, modulus1);
    let r2 = residue_mod(residue2, modulus2);
    let inverse = modular_inverse(modulus1, modulus2)?;
    let adjustment = residue_mod((r2 - r1) * inverse, modulus2);
    let modulus = modulus1 * modulus2;
    Some((residue_mod(r1 + modulus1 * adjustment, modulus), modulus))
}

/// Combine a nonempty list of pairwise coprime congruences.
pub fn chinese_remainder(congruences: &[(i128, i128)]) -> Option<(i128, i128)> {
    let (&(first_residue, first_modulus), rest) = congruences.split_first()?;
    let mut combined = (residue_mod(first_residue, first_modulus), first_modulus);
    for &(residue, modulus) in rest {
        combined = chinese_remainder_pair(combined.0, combined.1, residue, modulus)?;
    }
    Some(combined)
}

/// Return the representative in `[-modulus/2, modulus/2]`.
pub fn symmetric_residue(residue: i128, modulus: i128) -> i128 {
    assert!(modulus > 0, "modulus must be positive");
    let normalized = residue_mod(residue, modulus);
    if 2 * normalized > modulus {
        normalized - modulus
    } else {
        normalized
    }
}

/// Rational reconstruction from a residue modulo `modulus`.
///
/// If `residue = a / b (mod modulus)` and `|a|, b <= sqrt(modulus / 2)`,
/// this returns `a / b`. Otherwise it returns `None`.
pub fn rational_reconstruction(residue: i128, modulus: i128) -> Option<Ratio<i128>> {
    assert!(modulus > 1, "modulus must be greater than 1");
    let bound = integer_sqrt((modulus - 1) / 2);
    let target = symmetric_residue(residue, modulus);

    let (mut r0, mut r1) = (modulus, target);
    let (mut s0, mut s1) = (0i128, 1i128);

    while r1.abs() > bound {
        if r1 == 0 {
            return None;
        }
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (s0, s1) = (s1, s0 - q * s1);
    }

    let mut numerator = r1;
    let mut denominator = s1;
    if denominator < 0 {
        numerator = -numerator;
        denominator = -denominator;
    }

    if denominator == 0 || numerator.abs() > bound || denominator > bound {
        return None;
    }
    if numerator.gcd(&denominator) != 1 {
        return None;
    }
    if residue_mod(denominator * residue - numerator, modulus) != 0 {
        return None;
    }

    Some(Ratio::new(numerator, denominator))
}

fn modular_inverse(value: i128, modulus: i128) -> Option<i128> {
    let (gcd, inverse, _) = extended_gcd(value, modulus);
    (gcd == 1).then(|| residue_mod(inverse, modulus))
}

fn residue_mod(value: i128, modulus: i128) -> i128 {
    let mut residue = value % modulus;
    if residue < 0 {
        residue += modulus;
    }
    residue
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

fn integer_sqrt(n: i128) -> i128 {
    assert!(n >= 0, "square root input must be nonnegative");
    if n < 2 {
        return n;
    }

    let mut lo = 1i128;
    let mut hi = n.min(1 << 63);
    while lo <= hi {
        let mid = lo + (hi - lo) / 2;
        match mid.checked_mul(mid).map(|square| square.cmp(&n)) {
            Some(std::cmp::Ordering::Equal) => return mid,
            Some(std::cmp::Ordering::Less) => lo = mid + 1,
            _ => hi = mid - 1,
        }
    }
    hi
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chinese_remainder_pair() {
        let (residue, modulus) = chinese_remainder_pair(2, 5, 3, 7).unwrap();

        assert_eq!(modulus, 35);
        assert_eq!(residue % 5, 2);
        assert_eq!(residue % 7, 3);
        assert_eq!(residue, 17);
    }

    #[test]
    fn test_chinese_remainder_sequence() {
        let (residue, modulus) = chinese_remainder(&[(1, 5), (2, 7), (3, 11)]).unwrap();

        assert_eq!(modulus, 385);
        assert_eq!(residue % 5, 1);
        assert_eq!(residue % 7, 2);
        assert_eq!(residue % 11, 3);
    }

    #[test]
    fn test_symmetric_residue() {
        assert_eq!(symmetric_residue(98, 101), -3);
        assert_eq!(symmetric_residue(49, 101), 49);
    }

    #[test]
    fn test_rational_reconstruction() {
        let modulus = 1009;
        let inverse_7 = modular_inverse(7, modulus).unwrap();
        let residue = (5 * inverse_7) % modulus;

        assert_eq!(
            rational_reconstruction(residue, modulus),
            Some(Ratio::new(5, 7))
        );
    }
}
