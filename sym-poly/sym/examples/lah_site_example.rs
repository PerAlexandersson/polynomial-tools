use sym_poly_core::Partition;
use sym_poly_sym::{lah_symmetric_elementary, lah_symmetric_monomial};

fn partition(parts: &[u32]) -> Partition {
    Partition::new(parts.to_vec())
}

fn main() {
    let expected_elementary = [
        vec![(vec![3], 6), (vec![2, 1], 8), (vec![1, 1, 1], 1)],
        vec![(vec![2], 8), (vec![1, 1], 7)],
        vec![(vec![1], 6)],
        vec![(vec![], 1)],
    ];

    let expected_monomial = [
        vec![(vec![3], 1), (vec![2, 1], 11), (vec![1, 1, 1], 36)],
        vec![(vec![2], 7), (vec![1, 1], 22)],
        vec![(vec![1], 6)],
        vec![(vec![], 1)],
    ];

    for k in 1..=4 {
        let elementary = lah_symmetric_elementary(4, k);
        let monomial = lah_symmetric_monomial(4, k);

        for (parts, coeff) in &expected_elementary[k - 1] {
            assert_eq!(elementary.coefficient(&partition(parts)), *coeff);
        }
        assert_eq!(elementary.terms().len(), expected_elementary[k - 1].len());

        for (parts, coeff) in &expected_monomial[k - 1] {
            assert_eq!(monomial.coefficient(&partition(parts)), *coeff);
        }
        assert_eq!(monomial.terms().len(), expected_monomial[k - 1].len());

        println!("L_{{4,{k}}} in e-basis: {elementary}");
        println!("L_{{4,{k}}} in m-basis: {monomial}");
    }
}
