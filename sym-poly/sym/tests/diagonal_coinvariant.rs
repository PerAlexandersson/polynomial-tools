use num_rational::Ratio;
use sym_poly_core::Partition;
use sym_poly_multipoly::{
    quotient_action_matrices_by_multidegree_and_cycle_type, quotient_basis, GroebnerBasis,
    IndexedVariables, MonomialOrder, MultiPoly,
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
