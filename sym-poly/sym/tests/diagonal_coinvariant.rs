use num_rational::Ratio;
use std::collections::BTreeMap;
use sym_poly_core::linear_algebra::matrix_trace;
use sym_poly_core::sn_action::conjugacy_class_representatives;
use sym_poly_core::{chinese_remainder_pair, symmetric_residue, Partition, PrimeField, Ring};
use sym_poly_multipoly::{
    graded_quotient_component, quotient_action_matrices_by_multidegree_and_cycle_type,
    quotient_basis, GradedQuotientComponent, GroebnerBasis, IndexedVariables, MonomialOrder,
    MultiPoly, PolynomialQuotientSnModule,
};
use sym_poly_sym::multigraded_frobenius_from_trace_matrices;

type Q = Ratio<i64>;

fn q(n: i64) -> Q {
    Q::from_integer(n)
}

fn fp<const P: u64>(n: i64) -> PrimeField<P> {
    PrimeField::<P>::from_i64(n)
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

fn diagonal_power_sum_generator<C: Ring>(
    variables: &IndexedVariables,
    x_degree: u32,
    y_degree: u32,
) -> MultiPoly<C> {
    let mut polynomial = MultiPoly::zero(variables.num_vars());
    for index in 0..variables.num_indices() {
        let mut exponents = vec![0u32; variables.num_vars()];
        exponents[variables.variable_index(0, index)] = x_degree;
        exponents[variables.variable_index(1, index)] = y_degree;
        polynomial = polynomial + MultiPoly::monomial(variables.num_vars(), exponents, C::one());
    }
    polynomial
}

fn diagonal_power_sum_generators<C: Ring>(num_indices: usize) -> Vec<MultiPoly<C>> {
    let variables = IndexedVariables::new(2, num_indices);
    let mut generators = Vec::new();
    for total_degree in 1..=num_indices as u32 {
        for x_degree in (0..=total_degree).rev() {
            generators.push(diagonal_power_sum_generator::<C>(
                &variables,
                x_degree,
                total_degree - x_degree,
            ));
        }
    }
    generators
}

fn diagonal_coinvariant_s3_hilbert() -> BTreeMap<Vec<u32>, usize> {
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
}

fn factorial(n: usize) -> i64 {
    (1..=n as i64).product()
}

fn sign_character_value(cycle_type: &Partition) -> i64 {
    let n = cycle_type.size() as usize;
    if (n - cycle_type.parts().len()) % 2 == 0 {
        1
    } else {
        -1
    }
}

fn sign_multiplicity_mod_prime<const P: u64>(
    variables: &IndexedVariables,
    component: &GradedQuotientComponent<PrimeField<P>>,
) -> PrimeField<P> {
    let mut signed_trace_sum = PrimeField::<P>::zero();
    for (cycle_type, representative) in conjugacy_class_representatives(variables.num_indices()) {
        let class_size = factorial(variables.num_indices()) / cycle_type.z_coefficient() as i64;
        let sign_value = sign_character_value(&cycle_type);
        let action = component.action_matrix_by_index_permutation(variables, &representative);
        signed_trace_sum =
            signed_trace_sum + matrix_trace(&action) * fp::<P>(class_size * sign_value);
    }
    signed_trace_sum / fp::<P>(factorial(variables.num_indices()))
}

fn diagonal_coinvariant_s4_sign_hilbert_mod_prime<const P: u64>(
) -> BTreeMap<Vec<u32>, PrimeField<P>> {
    let variables = IndexedVariables::new(2, 4);
    let generators = diagonal_power_sum_generators::<PrimeField<P>>(4);
    let mut sign_hilbert = BTreeMap::new();

    for total_degree in 0..=6 {
        for x_degree in 0..=total_degree {
            let degree = vec![x_degree, total_degree - x_degree];
            let component = graded_quotient_component(&variables, &generators, &degree);
            let multiplicity = sign_multiplicity_mod_prime(&variables, &component);
            if multiplicity != PrimeField::<P>::zero() {
                sign_hilbert.insert(degree, multiplicity);
            }
        }
    }

    sign_hilbert
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
        diagonal_power_sum_generators::<Q>(3),
        MonomialOrder::Lex,
    )
    .expect("diagonal coinvariant quotient should be finite and S_3-invariant");

    assert_eq!(module.dimension(), 16);
    assert_eq!(
        module.hilbert_series_by_multidegree(),
        diagonal_coinvariant_s3_hilbert()
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

#[test]
#[ignore = "benchmark: run explicitly when comparing monomial orders on the diagonal S_3 path"]
fn diagonal_coinvariant_s3_dimension_all_standard_orders() {
    for order in MonomialOrder::STANDARD_ORDERS {
        let module = PolynomialQuotientSnModule::new(
            IndexedVariables::new(2, 3),
            diagonal_power_sum_generators::<Q>(3),
            order,
        )
        .expect("diagonal coinvariant quotient should be finite and S_3-invariant");

        assert_eq!(module.dimension(), 16, "{order} gave the wrong dimension");
        assert_eq!(
            module.hilbert_series_by_multidegree(),
            diagonal_coinvariant_s3_hilbert(),
            "{order} gave the wrong Hilbert series"
        );
    }
}

#[test]
#[ignore = "benchmark: modular S_4 diagonal coinvariant degree-by-degree path with CRT"]
fn diagonal_coinvariant_s4_sign_modular_crt_benchmark() {
    let sign_mod_32003 = diagonal_coinvariant_s4_sign_hilbert_mod_prime::<32003>();
    let sign_mod_32009 = diagonal_coinvariant_s4_sign_hilbert_mod_prime::<32009>();
    assert_eq!(
        sign_mod_32003.keys().collect::<Vec<_>>(),
        sign_mod_32009.keys().collect::<Vec<_>>()
    );

    let sign_hilbert: BTreeMap<_, _> = sign_mod_32003
        .iter()
        .map(|(degree, value_32003)| {
            let value_32009 = sign_mod_32009[degree];
            let (residue, modulus) = chinese_remainder_pair(
                value_32003.value() as i128,
                32003,
                value_32009.value() as i128,
                32009,
            )
            .unwrap();
            (degree.clone(), symmetric_residue(residue, modulus))
        })
        .collect();

    assert_eq!(
        sign_hilbert,
        BTreeMap::from([
            (vec![0, 6], 1),
            (vec![1, 3], 1),
            (vec![1, 4], 1),
            (vec![1, 5], 1),
            (vec![2, 2], 1),
            (vec![2, 3], 1),
            (vec![2, 4], 1),
            (vec![3, 1], 1),
            (vec![3, 2], 1),
            (vec![3, 3], 1),
            (vec![4, 1], 1),
            (vec![4, 2], 1),
            (vec![5, 1], 1),
            (vec![6, 0], 1),
        ])
    );
}
