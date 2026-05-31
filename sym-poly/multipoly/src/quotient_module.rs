//! Finite polynomial quotient modules with indexed `S_n` actions.

use std::collections::BTreeMap;
use std::fmt;

use sym_poly_core::linear_algebra::Matrix;
use sym_poly_core::{Field, Partition};

use crate::groebner::GroebnerBasis;
use crate::indexed_variables::{
    ideal_generators_are_sn_invariant,
    quotient_action_matrices_by_index_permutation_and_multidegree,
    quotient_action_matrices_by_multidegree_and_cycle_type, quotient_basis_multidegrees,
    IndexedVariables,
};
use crate::quotient::{quotient_basis, quotient_coordinates, QuotientBasis};
use crate::{MonomialOrder, MultiPoly};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PolynomialQuotientSnModuleError {
    NonFiniteStandardMonomialBasis,
    IdealNotSnInvariant,
}

impl fmt::Display for PolynomialQuotientSnModuleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFiniteStandardMonomialBasis => {
                write!(f, "could not certify a finite standard-monomial basis")
            }
            Self::IdealNotSnInvariant => write!(f, "ideal generators are not S_n-invariant"),
        }
    }
}

/// A finite quotient `R/I` with variables grouped as `x_{a,i}`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolynomialQuotientSnModule<C: Field> {
    variables: IndexedVariables,
    ideal_generators: Vec<MultiPoly<C>>,
    groebner_basis: GroebnerBasis<C>,
    quotient_basis: QuotientBasis,
}

impl<C: Field> PolynomialQuotientSnModule<C> {
    /// Build a quotient module and require `S_n`-invariance of the generators.
    pub fn new(
        variables: IndexedVariables,
        ideal_generators: Vec<MultiPoly<C>>,
        order: MonomialOrder,
    ) -> Result<Self, PolynomialQuotientSnModuleError> {
        Self::build(variables, ideal_generators, order, true)
    }

    /// Build a quotient module without checking `S_n`-invariance.
    pub fn new_without_invariance_check(
        variables: IndexedVariables,
        ideal_generators: Vec<MultiPoly<C>>,
        order: MonomialOrder,
    ) -> Result<Self, PolynomialQuotientSnModuleError> {
        Self::build(variables, ideal_generators, order, false)
    }

    pub fn variables(&self) -> &IndexedVariables {
        &self.variables
    }

    pub fn ideal_generators(&self) -> &[MultiPoly<C>] {
        &self.ideal_generators
    }

    pub fn groebner_basis(&self) -> &GroebnerBasis<C> {
        &self.groebner_basis
    }

    pub fn quotient_basis(&self) -> &QuotientBasis {
        &self.quotient_basis
    }

    pub fn dimension(&self) -> usize {
        self.quotient_basis.dimension()
    }

    pub fn normal_form(&self, polynomial: &MultiPoly<C>) -> MultiPoly<C> {
        self.groebner_basis.normal_form(polynomial)
    }

    pub fn coordinates(&self, polynomial: &MultiPoly<C>) -> Option<Vec<C>> {
        quotient_coordinates(polynomial, &self.groebner_basis, &self.quotient_basis)
    }

    pub fn multidegrees(&self) -> BTreeMap<Vec<u32>, Vec<usize>> {
        quotient_basis_multidegrees(&self.variables, &self.quotient_basis)
    }

    pub fn hilbert_series_by_multidegree(&self) -> BTreeMap<Vec<u32>, usize> {
        self.multidegrees()
            .into_iter()
            .map(|(degree, indices)| (degree, indices.len()))
            .collect()
    }

    pub fn action_matrices_by_index_permutation_and_multidegree(
        &self,
        index_permutation: &[usize],
    ) -> Option<BTreeMap<Vec<u32>, Matrix<C>>> {
        quotient_action_matrices_by_index_permutation_and_multidegree(
            &self.variables,
            &self.groebner_basis,
            &self.quotient_basis,
            index_permutation,
        )
    }

    pub fn action_matrices_by_multidegree_and_cycle_type(
        &self,
    ) -> Option<BTreeMap<Vec<u32>, BTreeMap<Partition, Matrix<C>>>> {
        quotient_action_matrices_by_multidegree_and_cycle_type(
            &self.variables,
            &self.groebner_basis,
            &self.quotient_basis,
        )
    }

    pub fn ideal_generators_are_sn_invariant(&self) -> bool {
        ideal_generators_are_sn_invariant(
            &self.variables,
            &self.ideal_generators,
            &self.groebner_basis,
        )
    }

    fn build(
        variables: IndexedVariables,
        ideal_generators: Vec<MultiPoly<C>>,
        order: MonomialOrder,
        check_invariance: bool,
    ) -> Result<Self, PolynomialQuotientSnModuleError> {
        assert!(
            ideal_generators
                .iter()
                .all(|polynomial| polynomial.num_vars() == variables.num_vars()),
            "all generators must use the indexed variable count"
        );
        if ideal_generators.is_empty() && variables.num_vars() > 0 {
            return Err(PolynomialQuotientSnModuleError::NonFiniteStandardMonomialBasis);
        }

        let groebner_basis = GroebnerBasis::new(ideal_generators.clone(), order);
        let quotient_basis = quotient_basis(&groebner_basis)
            .ok_or(PolynomialQuotientSnModuleError::NonFiniteStandardMonomialBasis)?;

        if check_invariance
            && !ideal_generators_are_sn_invariant(&variables, &ideal_generators, &groebner_basis)
        {
            return Err(PolynomialQuotientSnModuleError::IdealNotSnInvariant);
        }

        Ok(Self {
            variables,
            ideal_generators,
            groebner_basis,
            quotient_basis,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elementary_symmetric_generators;
    use num_rational::Ratio;
    use sym_poly_core::linear_algebra::matrix_trace;

    type Q = Ratio<i64>;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn mono(exponents: &[u32], coefficient: i64) -> MultiPoly<Q> {
        MultiPoly::monomial(exponents.len(), exponents.to_vec(), q(coefficient))
    }

    #[test]
    fn test_artin_s2_quotient_module() {
        let module = PolynomialQuotientSnModule::new(
            IndexedVariables::new(1, 2),
            elementary_symmetric_generators::<Q>(2),
            MonomialOrder::Lex,
        )
        .unwrap();

        assert_eq!(module.dimension(), 2);
        assert!(module.ideal_generators_are_sn_invariant());
        assert_eq!(
            module.multidegrees(),
            BTreeMap::from([(vec![0], vec![0]), (vec![1], vec![1])])
        );
        assert_eq!(
            module.hilbert_series_by_multidegree(),
            BTreeMap::from([(vec![0], 1), (vec![1], 1)])
        );

        let blocks = module
            .action_matrices_by_index_permutation_and_multidegree(&[1, 0])
            .unwrap();
        assert_eq!(matrix_trace(&blocks[&vec![0]]), q(1));
        assert_eq!(matrix_trace(&blocks[&vec![1]]), q(-1));
    }

    #[test]
    fn test_artin_s2_quotient_module_supports_standard_orders() {
        for order in MonomialOrder::STANDARD_ORDERS {
            let module = PolynomialQuotientSnModule::new(
                IndexedVariables::new(1, 2),
                elementary_symmetric_generators::<Q>(2),
                order,
            )
            .unwrap();

            assert_eq!(module.dimension(), 2, "{order} gave the wrong dimension");
            assert!(module.ideal_generators_are_sn_invariant());
            assert_eq!(
                module.hilbert_series_by_multidegree(),
                BTreeMap::from([(vec![0], 1), (vec![1], 1)])
            );
        }
    }

    #[test]
    fn test_rejects_non_invariant_ideal() {
        let result = PolynomialQuotientSnModule::new(
            IndexedVariables::new(1, 2),
            vec![mono(&[1, 0], 1), mono(&[0, 2], 1)],
            MonomialOrder::Lex,
        );

        assert_eq!(
            result.unwrap_err(),
            PolynomialQuotientSnModuleError::IdealNotSnInvariant
        );
    }

    #[test]
    fn test_rejects_empty_generators_in_positive_variables() {
        let result = PolynomialQuotientSnModule::<Q>::new(
            IndexedVariables::new(1, 2),
            Vec::new(),
            MonomialOrder::Lex,
        );

        assert_eq!(
            result.unwrap_err(),
            PolynomialQuotientSnModuleError::NonFiniteStandardMonomialBasis
        );
    }
}
