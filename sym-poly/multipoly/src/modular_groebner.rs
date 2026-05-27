//! Modular helpers for Groebner-basis experiments.
//!
//! These utilities support the standard modular workflow: reduce rational or
//! integer input modulo one or more good primes, compare leading-monomial
//! profiles, and lift monomial-aligned modular outputs by CRT or rational
//! reconstruction. They deliberately do not try to certify a full modular
//! Groebner reconstruction yet.

use std::collections::BTreeMap;
use std::fmt;

use num_bigint::BigInt;
use num_rational::Ratio;
use sym_poly_core::{
    chinese_remainder_pair, rational_reconstruction, symmetric_residue, PrimeField, Ring,
};

use crate::groebner::{GroebnerBasis, GroebnerOptions};
use crate::{MonomialOrder, MultiPoly};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModularGroebnerError {
    DenominatorDivisibleByPrime { denominator: i64, prime: u64 },
    IncompatibleVariableCounts { first: usize, second: usize },
    IncompatibleSupports,
    IncompatibleBasisLengths { first: usize, second: usize },
    NonCoprimePrimes { first: u64, second: u64 },
    RationalReconstructionFailed { monomial: Vec<u32>, modulus: i128 },
}

impl fmt::Display for ModularGroebnerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DenominatorDivisibleByPrime { denominator, prime } => write!(
                f,
                "denominator {denominator} is zero modulo the prime {prime}"
            ),
            Self::IncompatibleVariableCounts { first, second } => write!(
                f,
                "polynomials have incompatible variable counts {first} and {second}"
            ),
            Self::IncompatibleSupports => write!(f, "modular polynomials have different supports"),
            Self::IncompatibleBasisLengths { first, second } => {
                write!(f, "bases have incompatible lengths {first} and {second}")
            }
            Self::NonCoprimePrimes { first, second } => {
                write!(f, "moduli {first} and {second} are not coprime")
            }
            Self::RationalReconstructionFailed { monomial, modulus } => write!(
                f,
                "could not rationally reconstruct coefficient at {:?} modulo {}",
                monomial, modulus
            ),
        }
    }
}

pub fn reduce_i64_polynomial_mod_prime<const P: u64>(
    polynomial: &MultiPoly<i64>,
) -> MultiPoly<PrimeField<P>> {
    let terms = polynomial
        .terms()
        .iter()
        .filter_map(|(monomial, coefficient)| {
            let reduced = PrimeField::<P>::from_i64(*coefficient);
            (!reduced.is_zero()).then_some((monomial.clone(), reduced))
        })
        .collect();
    MultiPoly::from_terms(polynomial.num_vars(), terms)
}

pub fn reduce_i64_polynomials_mod_prime<const P: u64>(
    polynomials: &[MultiPoly<i64>],
) -> Vec<MultiPoly<PrimeField<P>>> {
    polynomials
        .iter()
        .map(reduce_i64_polynomial_mod_prime::<P>)
        .collect()
}

pub fn reduce_rational_i64_mod_prime<const P: u64>(
    value: &Ratio<i64>,
) -> Result<PrimeField<P>, ModularGroebnerError> {
    let denominator = *value.denom();
    let denominator_mod = PrimeField::<P>::from_i64(denominator);
    if denominator_mod.is_zero() {
        return Err(ModularGroebnerError::DenominatorDivisibleByPrime {
            denominator,
            prime: P,
        });
    }
    Ok(PrimeField::<P>::from_i64(*value.numer()) / denominator_mod)
}

pub fn reduce_rational_i64_polynomial_mod_prime<const P: u64>(
    polynomial: &MultiPoly<Ratio<i64>>,
) -> Result<MultiPoly<PrimeField<P>>, ModularGroebnerError> {
    let mut terms = BTreeMap::new();
    for (monomial, coefficient) in polynomial.terms() {
        let reduced = reduce_rational_i64_mod_prime::<P>(coefficient)?;
        if !reduced.is_zero() {
            terms.insert(monomial.clone(), reduced);
        }
    }
    Ok(MultiPoly::from_terms(polynomial.num_vars(), terms))
}

pub fn reduce_rational_i64_polynomials_mod_prime<const P: u64>(
    polynomials: &[MultiPoly<Ratio<i64>>],
) -> Result<Vec<MultiPoly<PrimeField<P>>>, ModularGroebnerError> {
    polynomials
        .iter()
        .map(reduce_rational_i64_polynomial_mod_prime::<P>)
        .collect()
}

pub fn modular_groebner_basis_i64_mod_prime<const P: u64>(
    generators: &[MultiPoly<i64>],
    order: MonomialOrder,
    options: GroebnerOptions,
) -> GroebnerBasis<PrimeField<P>> {
    GroebnerBasis::with_options(
        reduce_i64_polynomials_mod_prime::<P>(generators),
        order,
        options,
    )
}

pub fn modular_groebner_basis_rational_i64_mod_prime<const P: u64>(
    generators: &[MultiPoly<Ratio<i64>>],
    order: MonomialOrder,
    options: GroebnerOptions,
) -> Result<GroebnerBasis<PrimeField<P>>, ModularGroebnerError> {
    Ok(GroebnerBasis::with_options(
        reduce_rational_i64_polynomials_mod_prime::<P>(generators)?,
        order,
        options,
    ))
}

pub fn groebner_leading_monomials<C: Ring>(basis: &GroebnerBasis<C>) -> Vec<Vec<u32>> {
    basis
        .leading_terms
        .iter()
        .map(|leading_term| leading_term.exponents.clone())
        .collect()
}

pub fn modular_leading_monomials_match_i64<const P: u64, const Q: u64>(
    generators: &[MultiPoly<i64>],
    order: MonomialOrder,
    options: GroebnerOptions,
) -> bool {
    let first = modular_groebner_basis_i64_mod_prime::<P>(generators, order, options);
    let second = modular_groebner_basis_i64_mod_prime::<Q>(generators, order, options);
    groebner_leading_monomials(&first) == groebner_leading_monomials(&second)
}

pub fn modular_leading_monomials_match_rational_i64<const P: u64, const Q: u64>(
    generators: &[MultiPoly<Ratio<i64>>],
    order: MonomialOrder,
    options: GroebnerOptions,
) -> Result<bool, ModularGroebnerError> {
    let first = modular_groebner_basis_rational_i64_mod_prime::<P>(generators, order, options)?;
    let second = modular_groebner_basis_rational_i64_mod_prime::<Q>(generators, order, options)?;
    Ok(groebner_leading_monomials(&first) == groebner_leading_monomials(&second))
}

pub fn crt_lift_prime_field_polynomial_pair<const P: u64, const Q: u64>(
    first: &MultiPoly<PrimeField<P>>,
    second: &MultiPoly<PrimeField<Q>>,
) -> Result<MultiPoly<BigInt>, ModularGroebnerError> {
    assert_same_ring(first, second)?;
    assert_same_support(first, second)?;

    let mut terms = BTreeMap::new();
    for monomial in first.terms().keys() {
        let first_value = first.coefficient(monomial).value() as i128;
        let second_value = second.coefficient(monomial).value() as i128;
        let (residue, modulus) =
            chinese_remainder_pair(first_value, P as i128, second_value, Q as i128).ok_or(
                ModularGroebnerError::NonCoprimePrimes {
                    first: P,
                    second: Q,
                },
            )?;
        let coefficient = BigInt::from(symmetric_residue(residue, modulus));
        if coefficient != BigInt::from(0) {
            terms.insert(monomial.clone(), coefficient);
        }
    }

    Ok(MultiPoly::from_terms(first.num_vars(), terms))
}

pub fn rational_reconstruct_prime_field_polynomial_pair<const P: u64, const Q: u64>(
    first: &MultiPoly<PrimeField<P>>,
    second: &MultiPoly<PrimeField<Q>>,
) -> Result<MultiPoly<Ratio<BigInt>>, ModularGroebnerError> {
    assert_same_ring(first, second)?;
    assert_same_support(first, second)?;

    let mut terms = BTreeMap::new();
    for monomial in first.terms().keys() {
        let first_value = first.coefficient(monomial).value() as i128;
        let second_value = second.coefficient(monomial).value() as i128;
        let (residue, modulus) =
            chinese_remainder_pair(first_value, P as i128, second_value, Q as i128).ok_or(
                ModularGroebnerError::NonCoprimePrimes {
                    first: P,
                    second: Q,
                },
            )?;
        let reconstructed = rational_reconstruction(residue, modulus).ok_or_else(|| {
            ModularGroebnerError::RationalReconstructionFailed {
                monomial: monomial.clone(),
                modulus,
            }
        })?;
        let coefficient = Ratio::new(
            BigInt::from(*reconstructed.numer()),
            BigInt::from(*reconstructed.denom()),
        );
        if !coefficient.is_zero() {
            terms.insert(monomial.clone(), coefficient);
        }
    }

    Ok(MultiPoly::from_terms(first.num_vars(), terms))
}

pub fn crt_lift_prime_field_basis_pair<const P: u64, const Q: u64>(
    first: &[MultiPoly<PrimeField<P>>],
    second: &[MultiPoly<PrimeField<Q>>],
) -> Result<Vec<MultiPoly<BigInt>>, ModularGroebnerError> {
    if first.len() != second.len() {
        return Err(ModularGroebnerError::IncompatibleBasisLengths {
            first: first.len(),
            second: second.len(),
        });
    }
    first
        .iter()
        .zip(second.iter())
        .map(|(first_polynomial, second_polynomial)| {
            crt_lift_prime_field_polynomial_pair(first_polynomial, second_polynomial)
        })
        .collect()
}

pub fn rational_reconstruct_prime_field_basis_pair<const P: u64, const Q: u64>(
    first: &[MultiPoly<PrimeField<P>>],
    second: &[MultiPoly<PrimeField<Q>>],
) -> Result<Vec<MultiPoly<Ratio<BigInt>>>, ModularGroebnerError> {
    if first.len() != second.len() {
        return Err(ModularGroebnerError::IncompatibleBasisLengths {
            first: first.len(),
            second: second.len(),
        });
    }
    first
        .iter()
        .zip(second.iter())
        .map(|(first_polynomial, second_polynomial)| {
            rational_reconstruct_prime_field_polynomial_pair(first_polynomial, second_polynomial)
        })
        .collect()
}

fn assert_same_ring<const P: u64, const Q: u64>(
    first: &MultiPoly<PrimeField<P>>,
    second: &MultiPoly<PrimeField<Q>>,
) -> Result<(), ModularGroebnerError> {
    if first.num_vars() != second.num_vars() {
        return Err(ModularGroebnerError::IncompatibleVariableCounts {
            first: first.num_vars(),
            second: second.num_vars(),
        });
    }
    Ok(())
}

fn assert_same_support<const P: u64, const Q: u64>(
    first: &MultiPoly<PrimeField<P>>,
    second: &MultiPoly<PrimeField<Q>>,
) -> Result<(), ModularGroebnerError> {
    if first.terms().keys().eq(second.terms().keys()) {
        Ok(())
    } else {
        Err(ModularGroebnerError::IncompatibleSupports)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type Q = Ratio<i64>;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn mono(num_vars: usize, exponents: &[u32], coefficient: i64) -> MultiPoly<Q> {
        MultiPoly::monomial(num_vars, exponents.to_vec(), q(coefficient))
    }

    fn artin_s3_generators() -> Vec<MultiPoly<Q>> {
        vec![
            mono(3, &[1, 0, 0], 1) + mono(3, &[0, 1, 0], 1) + mono(3, &[0, 0, 1], 1),
            mono(3, &[1, 1, 0], 1) + mono(3, &[1, 0, 1], 1) + mono(3, &[0, 1, 1], 1),
            mono(3, &[1, 1, 1], 1),
        ]
    }

    #[test]
    fn test_reduce_rational_rejects_bad_prime() {
        let coefficient = Ratio::new(1, 5);

        assert_eq!(
            reduce_rational_i64_mod_prime::<5>(&coefficient),
            Err(ModularGroebnerError::DenominatorDivisibleByPrime {
                denominator: 5,
                prime: 5
            })
        );
        assert_eq!(
            reduce_rational_i64_mod_prime::<7>(&coefficient).unwrap(),
            PrimeField::<7>::from_i64(3)
        );
    }

    #[test]
    fn test_modular_groebner_leading_monomials_match_artin_s3() {
        let generators = artin_s3_generators();
        let basis_101 = modular_groebner_basis_rational_i64_mod_prime::<101>(
            &generators,
            MonomialOrder::Lex,
            GroebnerOptions::default(),
        )
        .unwrap();
        let basis_103 = modular_groebner_basis_rational_i64_mod_prime::<103>(
            &generators,
            MonomialOrder::Lex,
            GroebnerOptions::default(),
        )
        .unwrap();
        let leading_monomials = groebner_leading_monomials(&basis_101);

        assert_eq!(leading_monomials, groebner_leading_monomials(&basis_103));
        assert_eq!(
            leading_monomials,
            vec![vec![1, 0, 0], vec![0, 2, 0], vec![0, 0, 3]]
        );
        assert_eq!(
            modular_leading_monomials_match_rational_i64::<101, 103>(
                &generators,
                MonomialOrder::Lex,
                GroebnerOptions::default()
            ),
            Ok(true)
        );
    }

    #[test]
    fn test_crt_lift_prime_field_polynomial_pair() {
        let first = MultiPoly::monomial(2, vec![1, 0], PrimeField::<101>::from_i64(-3))
            + MultiPoly::monomial(2, vec![0, 1], PrimeField::<101>::from_i64(5));
        let second = MultiPoly::monomial(2, vec![1, 0], PrimeField::<103>::from_i64(-3))
            + MultiPoly::monomial(2, vec![0, 1], PrimeField::<103>::from_i64(5));

        let lifted = crt_lift_prime_field_polynomial_pair(&first, &second).unwrap();

        assert_eq!(
            lifted,
            MultiPoly::monomial(2, vec![1, 0], BigInt::from(-3))
                + MultiPoly::monomial(2, vec![0, 1], BigInt::from(5))
        );
    }

    #[test]
    fn test_rational_reconstruct_prime_field_polynomial_pair() {
        let first = MultiPoly::monomial(
            1,
            vec![2],
            PrimeField::<101>::from_i64(2) / PrimeField::<101>::from_i64(3),
        );
        let second = MultiPoly::monomial(
            1,
            vec![2],
            PrimeField::<103>::from_i64(2) / PrimeField::<103>::from_i64(3),
        );

        let lifted = rational_reconstruct_prime_field_polynomial_pair(&first, &second).unwrap();

        assert_eq!(
            lifted,
            MultiPoly::monomial(1, vec![2], Ratio::new(BigInt::from(2), BigInt::from(3)))
        );
    }
}
