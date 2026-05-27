use std::collections::BTreeMap;

use num_rational::Ratio;
use sym_poly_core::sn_action::{identity_permutation, simple_transposition};
use sym_poly_core::Partition;
use sym_poly_multipoly::{
    quotient_action_matrix_by_permutation, quotient_basis, quotient_basis_degrees, GroebnerBasis,
    MonomialOrder, MultiPoly,
};
use sym_poly_sym::frobenius_from_trace_matrices;

type Q = Ratio<i64>;

fn q(n: i64) -> Q {
    Q::from_integer(n)
}

fn p(parts: &[u32]) -> Partition {
    Partition::new(parts.to_vec())
}

fn mono(exponents: &[u32], coefficient: i64) -> MultiPoly<Q> {
    MultiPoly::monomial(2, exponents.to_vec(), q(coefficient))
}

fn trace_on_indices(matrix: &[Vec<Q>], indices: &[usize]) -> Q {
    indices
        .iter()
        .fold(q(0), |acc, &index| acc + matrix[index][index].clone())
}

#[test]
fn artin_coinvariant_s2_graded_frobenius() {
    let e1 = mono(&[1, 0], 1) + mono(&[0, 1], 1);
    let e2 = mono(&[1, 1], 1);
    let gb = GroebnerBasis::new(vec![e1, e2], MonomialOrder::Lex);
    let basis = quotient_basis(&gb).expect("Artin quotient should be finite");

    assert_eq!(basis.monomials, vec![vec![0, 0], vec![0, 1]]);
    assert_eq!(
        quotient_basis_degrees(&basis),
        BTreeMap::from([(0, vec![0]), (1, vec![1])])
    );

    let id_action =
        quotient_action_matrix_by_permutation(&gb, &basis, &identity_permutation(2)).unwrap();
    let s_action =
        quotient_action_matrix_by_permutation(&gb, &basis, &simple_transposition(2, 0)).unwrap();
    let degrees = quotient_basis_degrees(&basis);

    let degree_zero_chars = BTreeMap::from([
        (p(&[1, 1]), trace_on_indices(&id_action, &degrees[&0])),
        (p(&[2]), trace_on_indices(&s_action, &degrees[&0])),
    ]);
    let degree_one_chars = BTreeMap::from([
        (p(&[1, 1]), trace_on_indices(&id_action, &degrees[&1])),
        (p(&[2]), trace_on_indices(&s_action, &degrees[&1])),
    ]);

    let degree_zero = frobenius_from_trace_matrices(
        &degree_zero_chars
            .into_iter()
            .map(|(cycle_type, trace)| (cycle_type, vec![vec![trace]]))
            .collect(),
    )
    .to_schur_basis();
    let degree_one = frobenius_from_trace_matrices(
        &degree_one_chars
            .into_iter()
            .map(|(cycle_type, trace)| (cycle_type, vec![vec![trace]]))
            .collect(),
    )
    .to_schur_basis();

    assert_eq!(degree_zero.coefficient(&p(&[2])), q(1));
    assert_eq!(degree_zero.terms().len(), 1);
    assert_eq!(degree_one.coefficient(&p(&[1, 1])), q(1));
    assert_eq!(degree_one.terms().len(), 1);
}
