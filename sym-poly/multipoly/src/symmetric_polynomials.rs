//! Small symmetric-polynomial constructors in ordinary variables.

use std::collections::BTreeMap;

use sym_poly_core::Ring;

use crate::MultiPoly;

/// The elementary symmetric polynomial `e_k(x_1, ..., x_n)`.
///
/// By convention `e_0 = 1`; for `k > n`, this returns zero.
pub fn elementary_symmetric_polynomial<C: Ring>(num_vars: usize, degree: usize) -> MultiPoly<C> {
    if degree == 0 {
        return MultiPoly::constant(num_vars, C::one());
    }
    if degree > num_vars {
        return MultiPoly::zero(num_vars);
    }

    let mut terms = BTreeMap::new();
    let mut chosen = Vec::with_capacity(degree);
    elementary_subsets(num_vars, degree, 0, &mut chosen, &mut terms);
    MultiPoly::from_terms(num_vars, terms)
}

/// The Artin coinvariant ideal generators `e_1, ..., e_n`.
pub fn elementary_symmetric_generators<C: Ring>(num_vars: usize) -> Vec<MultiPoly<C>> {
    (1..=num_vars)
        .map(|degree| elementary_symmetric_polynomial(num_vars, degree))
        .collect()
}

fn elementary_subsets<C: Ring>(
    num_vars: usize,
    degree: usize,
    start: usize,
    chosen: &mut Vec<usize>,
    terms: &mut BTreeMap<Vec<u32>, C>,
) {
    if chosen.len() == degree {
        let mut exponents = vec![0u32; num_vars];
        for &i in chosen.iter() {
            exponents[i] = 1;
        }
        terms.insert(exponents, C::one());
        return;
    }

    let remaining = degree - chosen.len();
    for i in start..=num_vars - remaining {
        chosen.push(i);
        elementary_subsets(num_vars, degree, i + 1, chosen, terms);
        chosen.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Ratio;

    type Q = Ratio<i64>;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    #[test]
    fn test_elementary_symmetric_polynomial() {
        let e2 = elementary_symmetric_polynomial::<Q>(3, 2);
        let expected = MultiPoly::monomial(3, vec![1, 1, 0], q(1))
            + MultiPoly::monomial(3, vec![1, 0, 1], q(1))
            + MultiPoly::monomial(3, vec![0, 1, 1], q(1));

        assert_eq!(e2, expected);
        assert_eq!(
            elementary_symmetric_polynomial::<Q>(3, 0),
            MultiPoly::constant(3, q(1))
        );
        assert_eq!(
            elementary_symmetric_polynomial::<Q>(3, 4),
            MultiPoly::zero(3)
        );
    }

    #[test]
    fn test_elementary_symmetric_generators() {
        let generators = elementary_symmetric_generators::<Q>(3);

        assert_eq!(generators.len(), 3);
        assert_eq!(generators[0].total_degree(), Some(1));
        assert_eq!(generators[1].total_degree(), Some(2));
        assert_eq!(generators[2], MultiPoly::monomial(3, vec![1, 1, 1], q(1)));
    }
}
