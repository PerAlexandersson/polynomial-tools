use num_rational::Ratio;
use std::collections::BTreeMap;
use sym_poly_core::Partition;
use sym_poly_multipoly::{
    quotient_action_matrices_by_multidegree_and_cycle_type, quotient_basis, GroebnerBasis,
    IndexedVariables, MonomialOrder, MultiPoly, PolynomialQuotientSnModule,
};
use sym_poly_sym::multigraded_frobenius_from_trace_matrices;

type Q = Ratio<i64>;

fn q(n: i64) -> Q {
    Q::from_integer(n)
}

fn p(parts: &[u32]) -> Partition {
    Partition::new(parts.to_vec())
}

fn mono(exponents: &[u32], coefficient: i64) -> MultiPoly<Q> {
    MultiPoly::monomial(4, exponents.to_vec(), q(coefficient))
}

fn diagonal_coinvariant_s2_generators() -> Vec<MultiPoly<Q>> {
    vec![
        mono(&[1, 0, 0, 0], 1) + mono(&[0, 1, 0, 0], 1),
        mono(&[0, 0, 1, 0], 1) + mono(&[0, 0, 0, 1], 1),
        mono(&[2, 0, 0, 0], 1) + mono(&[0, 2, 0, 0], 1),
        mono(&[1, 0, 1, 0], 1) + mono(&[0, 1, 0, 1], 1),
        mono(&[0, 0, 2, 0], 1) + mono(&[0, 0, 0, 2], 1),
    ]
}

fn diagonal_power_sum_generator(
    variables: &IndexedVariables,
    x_degree: u32,
    y_degree: u32,
) -> MultiPoly<Q> {
    let mut polynomial = MultiPoly::zero(variables.num_vars());
    for index in 0..variables.num_indices() {
        let mut exponents = vec![0u32; variables.num_vars()];
        exponents[variables.variable_index(0, index)] = x_degree;
        exponents[variables.variable_index(1, index)] = y_degree;
        polynomial = polynomial + MultiPoly::monomial(variables.num_vars(), exponents, q(1));
    }
    polynomial
}

fn diagonal_power_sum_generators(num_indices: usize) -> Vec<MultiPoly<Q>> {
    let variables = IndexedVariables::new(2, num_indices);
    let mut generators = Vec::new();
    for total_degree in 1..=num_indices as u32 {
        for x_degree in (0..=total_degree).rev() {
            generators.push(diagonal_power_sum_generator(
                &variables,
                x_degree,
                total_degree - x_degree,
            ));
        }
    }
    generators
}

#[test]
fn diagonal_coinvariant_s2_multigraded_frobenius() {
    let variables = IndexedVariables::new(2, 2);
    let gb = GroebnerBasis::new(diagonal_coinvariant_s2_generators(), MonomialOrder::Lex);
    let basis = quotient_basis(&gb).expect("diagonal coinvariant quotient should be finite");

    assert_eq!(basis.dimension(), 3);

    let matrices =
        quotient_action_matrices_by_multidegree_and_cycle_type(&variables, &gb, &basis).unwrap();
    let frobenius = multigraded_frobenius_from_trace_matrices(&matrices);

    let degree_00 = frobenius[&vec![0, 0]].to_schur_basis();
    let degree_10 = frobenius[&vec![1, 0]].to_schur_basis();
    let degree_01 = frobenius[&vec![0, 1]].to_schur_basis();

    assert_eq!(frobenius.len(), 3);
    assert_eq!(degree_00.coefficient(&p(&[2])), q(1));
    assert_eq!(degree_00.terms().len(), 1);
    assert_eq!(degree_10.coefficient(&p(&[1, 1])), q(1));
    assert_eq!(degree_10.terms().len(), 1);
    assert_eq!(degree_01.coefficient(&p(&[1, 1])), q(1));
    assert_eq!(degree_01.terms().len(), 1);
}

#[test]
#[ignore = "benchmark: run explicitly when checking the diagonal coinvariant S_3 Groebner path"]
fn diagonal_coinvariant_s3_dimension_benchmark() {
    let module = PolynomialQuotientSnModule::new(
        IndexedVariables::new(2, 3),
        diagonal_power_sum_generators(3),
        MonomialOrder::Lex,
    )
    .expect("diagonal coinvariant quotient should be finite and S_3-invariant");

    assert_eq!(module.dimension(), 16);
    assert_eq!(
        module.hilbert_series_by_multidegree(),
        BTreeMap::from([
            (vec![0, 0], 1),
            (vec![0, 1], 2),
            (vec![0, 2], 2),
            (vec![0, 3], 1),
            (vec![1, 0], 2),
            (vec![1, 1], 3),
            (vec![1, 2], 1),
            (vec![2, 0], 2),
            (vec![2, 1], 1),
            (vec![3, 0], 1),
        ])
    );

    let matrices = module
        .action_matrices_by_multidegree_and_cycle_type()
        .expect("diagonal quotient should carry the diagonal S_3 action");
    let frobenius = multigraded_frobenius_from_trace_matrices(&matrices);
    let trivial_shape = p(&[3]);
    for (degree, frobenius_degree) in &frobenius {
        let schur = frobenius_degree.to_schur_basis();
        let expected = if degree == &vec![0, 0] { q(1) } else { q(0) };
        assert_eq!(schur.coefficient(&trivial_shape), expected);
    }
    let degree_00 = frobenius[&vec![0, 0]].to_schur_basis();

    assert_eq!(degree_00.coefficient(&p(&[3])), q(1));
    assert_eq!(degree_00.terms().len(), 1);
}
