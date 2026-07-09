//! Exact root-of-unity checks for cyclic sieving.
//!
//! The implementation avoids floating point arithmetic.  For a primitive
//! `r`-th root of unity, `P(zeta) = m` if and only if `P(q) - m` is divisible by
//! the cyclotomic polynomial `Phi_r(q)`.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};
use std::collections::BTreeMap;

/// Exact evaluation data for one root of unity.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RootOfUnityEvaluation {
    pub group_order: usize,
    pub power: usize,
    pub root_order: usize,
    pub integer_value: Option<BigInt>,
    pub remainder: Vec<BigInt>,
}

/// One CSP equality check against a supplied fixed-point count.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CyclicSievingPowerCheck {
    pub power: usize,
    pub root_order: usize,
    pub expected_fixed_points: BigInt,
    pub evaluation: RootOfUnityEvaluation,
    pub holds: bool,
}

/// CSP report for one polynomial and one cyclic group order.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CyclicSievingReport {
    pub order: usize,
    pub coefficients: Vec<BigInt>,
    pub fixed_counts: Option<Vec<BigInt>>,
    pub evaluations: Vec<RootOfUnityEvaluation>,
    pub checks: Vec<CyclicSievingPowerCheck>,
    pub holds: Option<bool>,
}

/// Sequence-level CSP report for one row and several candidate group orders.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CyclicSievingSequenceItem {
    pub row: usize,
    pub index: isize,
    pub coefficients: Vec<BigInt>,
    pub candidate_orders: Vec<CyclicSievingReport>,
}

fn trim(mut p: Vec<BigInt>) -> Vec<BigInt> {
    while p.last().is_some_and(|c| c.is_zero()) {
        p.pop();
    }
    p
}

fn degree(p: &[BigInt]) -> Option<usize> {
    p.iter().rposition(|c| !c.is_zero())
}

fn divisors(n: usize) -> Vec<usize> {
    let mut out = Vec::new();
    for d in 1..=n {
        if n % d == 0 {
            out.push(d);
        }
    }
    out
}

#[cfg(test)]
fn poly_mul(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![BigInt::zero(); a.len() + b.len() - 1];
    for (i, ca) in a.iter().enumerate() {
        if ca.is_zero() {
            continue;
        }
        for (j, cb) in b.iter().enumerate() {
            if !cb.is_zero() {
                out[i + j] += ca * cb;
            }
        }
    }
    trim(out)
}

fn exact_div_monic(mut numerator: Vec<BigInt>, divisor: &[BigInt]) -> Option<Vec<BigInt>> {
    numerator = trim(numerator);
    let divisor = trim(divisor.to_vec());
    let dd = degree(&divisor)?;
    if divisor[dd] != BigInt::one() {
        return None;
    }
    let nd = match degree(&numerator) {
        Some(nd) => nd,
        None => return Some(Vec::new()),
    };
    if nd < dd {
        return Some(Vec::new());
    }

    let mut quotient = vec![BigInt::zero(); nd - dd + 1];
    for shift in (0..=nd - dd).rev() {
        let coeff = numerator[dd + shift].clone();
        if coeff.is_zero() {
            continue;
        }
        quotient[shift] = coeff.clone();
        for (j, dcoeff) in divisor.iter().enumerate().take(dd + 1) {
            numerator[j + shift] -= &coeff * dcoeff;
        }
    }
    if trim(numerator).is_empty() {
        Some(trim(quotient))
    } else {
        None
    }
}

fn remainder_monic(mut numerator: Vec<BigInt>, divisor: &[BigInt]) -> Vec<BigInt> {
    numerator = trim(numerator);
    let divisor = trim(divisor.to_vec());
    let Some(dd) = degree(&divisor) else {
        panic!("division by zero polynomial");
    };
    assert_eq!(divisor[dd], BigInt::one(), "divisor must be monic");
    while let Some(nd) = degree(&numerator) {
        if nd < dd {
            break;
        }
        let shift = nd - dd;
        let coeff = numerator[nd].clone();
        if coeff.is_zero() {
            numerator.pop();
            continue;
        }
        for (j, dcoeff) in divisor.iter().enumerate().take(dd + 1) {
            numerator[j + shift] -= &coeff * dcoeff;
        }
        numerator = trim(numerator);
    }
    trim(numerator)
}

fn cyclotomic_cached(n: usize, cache: &mut BTreeMap<usize, Vec<BigInt>>) -> Vec<BigInt> {
    if let Some(poly) = cache.get(&n) {
        return poly.clone();
    }
    let poly = if n == 1 {
        vec![BigInt::from(-1), BigInt::one()]
    } else {
        let mut p = vec![BigInt::zero(); n + 1];
        p[0] = BigInt::from(-1);
        p[n] = BigInt::one();
        for d in divisors(n).into_iter().filter(|&d| d < n) {
            let phi_d = cyclotomic_cached(d, cache);
            p = exact_div_monic(p, &phi_d).expect("cyclotomic division is exact");
        }
        p
    };
    cache.insert(n, poly.clone());
    poly
}

/// Return the `n`-th cyclotomic polynomial in ascending coefficient order.
pub fn cyclotomic_polynomial_bigint(n: usize) -> Vec<BigInt> {
    assert!(n >= 1, "cyclotomic index must be positive");
    cyclotomic_cached(n, &mut BTreeMap::new())
}

fn evaluate_at_one(coeffs: &[BigInt]) -> BigInt {
    coeffs.iter().fold(BigInt::zero(), |acc, coeff| acc + coeff)
}

/// Exactly reduce `P(q)` at `q = zeta_order^power`.
pub fn root_of_unity_evaluation_bigint(
    coeffs: &[BigInt],
    group_order: usize,
    power: usize,
) -> RootOfUnityEvaluation {
    assert!(group_order >= 1, "group order must be positive");
    let power = power % group_order;
    let root_order = if power == 0 {
        1
    } else {
        group_order / group_order.gcd(&power)
    };
    if root_order == 1 {
        let value = evaluate_at_one(coeffs);
        return RootOfUnityEvaluation {
            group_order,
            power,
            root_order,
            integer_value: Some(value.clone()),
            remainder: vec![value],
        };
    }

    let phi = cyclotomic_polynomial_bigint(root_order);
    let remainder = remainder_monic(coeffs.to_vec(), &phi);
    let integer_value = match remainder.as_slice() {
        [] => Some(BigInt::zero()),
        [constant] => Some(constant.clone()),
        _ => None,
    };
    RootOfUnityEvaluation {
        group_order,
        power,
        root_order,
        integer_value,
        remainder,
    }
}

/// Check one CSP candidate.  If `fixed_counts` is `None`, this only profiles
/// root-of-unity evaluations.
pub fn cyclic_sieving_report_bigint(
    coeffs: &[BigInt],
    order: usize,
    fixed_counts: Option<&[BigInt]>,
) -> Result<CyclicSievingReport, String> {
    if order == 0 {
        return Err("cyclic group order must be positive".to_string());
    }
    if let Some(counts) = fixed_counts {
        if counts.len() != order {
            return Err(format!(
                "expected {order} fixed-point counts for order {order}, got {}",
                counts.len()
            ));
        }
    }

    let coefficients = trim(coeffs.to_vec());
    let mut evaluations = Vec::with_capacity(order);
    let mut checks = Vec::new();
    for power in 0..order {
        let evaluation = root_of_unity_evaluation_bigint(&coefficients, order, power);
        if let Some(counts) = fixed_counts {
            let expected = counts[power].clone();
            let holds = evaluation.integer_value.as_ref() == Some(&expected);
            checks.push(CyclicSievingPowerCheck {
                power,
                root_order: evaluation.root_order,
                expected_fixed_points: expected,
                evaluation: evaluation.clone(),
                holds,
            });
        }
        evaluations.push(evaluation);
    }
    let holds = fixed_counts.map(|_| checks.iter().all(|check| check.holds));
    Ok(CyclicSievingReport {
        order,
        coefficients,
        fixed_counts: fixed_counts.map(|counts| counts.to_vec()),
        evaluations,
        checks,
        holds,
    })
}

/// Build sequence CSP profiles for candidate orders `index + offset`.
pub fn cyclic_sieving_sequence_reports_bigint(
    polynomials: &[Vec<BigInt>],
    first_index: isize,
    offsets: &[isize],
    fixed_counts: &BTreeMap<(isize, usize), Vec<BigInt>>,
) -> Vec<CyclicSievingSequenceItem> {
    polynomials
        .iter()
        .enumerate()
        .map(|(row, coeffs)| {
            let index = first_index + row as isize;
            let mut candidate_orders = Vec::new();
            for &offset in offsets {
                let order_index = index + offset;
                if order_index <= 0 {
                    continue;
                }
                let order = order_index as usize;
                let counts = fixed_counts.get(&(index, order)).map(Vec::as_slice);
                if let Ok(report) = cyclic_sieving_report_bigint(coeffs, order, counts) {
                    candidate_orders.push(report);
                }
            }
            CyclicSievingSequenceItem {
                row,
                index,
                coefficients: trim(coeffs.clone()),
                candidate_orders,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b(values: &[i64]) -> Vec<BigInt> {
        values.iter().map(|&v| BigInt::from(v)).collect()
    }

    #[test]
    fn cyclotomic_polynomials_are_exact() {
        assert_eq!(cyclotomic_polynomial_bigint(1), b(&[-1, 1]));
        assert_eq!(cyclotomic_polynomial_bigint(2), b(&[1, 1]));
        assert_eq!(cyclotomic_polynomial_bigint(3), b(&[1, 1, 1]));
        assert_eq!(cyclotomic_polynomial_bigint(4), b(&[1, 0, 1]));
        assert_eq!(
            poly_mul(&cyclotomic_polynomial_bigint(3), &b(&[-1, 1])),
            b(&[-1, 0, 0, 1])
        );
    }

    #[test]
    fn root_of_unity_integer_profiles() {
        let p = b(&[1, 1, 1, 1]);
        let at_one = root_of_unity_evaluation_bigint(&p, 4, 0);
        assert_eq!(at_one.integer_value, Some(BigInt::from(4)));

        let at_minus_one = root_of_unity_evaluation_bigint(&p, 4, 2);
        assert_eq!(at_minus_one.integer_value, Some(BigInt::from(0)));

        let at_i = root_of_unity_evaluation_bigint(&p, 4, 1);
        assert_eq!(at_i.integer_value, Some(BigInt::from(0)));
    }

    #[test]
    fn checks_basic_csp() {
        let p = b(&[1, 1]);
        let counts = b(&[2, 0]);
        let report = cyclic_sieving_report_bigint(&p, 2, Some(&counts)).unwrap();
        assert_eq!(report.holds, Some(true));
    }
}
