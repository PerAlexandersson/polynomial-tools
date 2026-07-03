use combinatoric_core::Partition;
use sym_poly_sym::{shifted_multiset_tableau_distribution, ShiftedMultisetTableau};

fn partition(parts: &[u32]) -> Partition {
    Partition::new(parts.to_vec())
}

fn main() {
    let shape = partition(&[2, 1]);
    let tableaux = ShiftedMultisetTableau::enumerate_smt0(&shape, 2, 4);
    let distribution = shifted_multiset_tableau_distribution(&shape, 2, 4);

    assert_eq!(tableaux.len(), 8);
    assert_eq!(distribution[&(vec![3, 1], vec![1, 0])], 1);
    assert_eq!(distribution[&(vec![3, 1], vec![0, 1])], 1);
    assert_eq!(distribution[&(vec![2, 2], vec![1, 0])], 2);
    assert_eq!(distribution[&(vec![2, 2], vec![0, 1])], 2);
    assert_eq!(distribution[&(vec![1, 3], vec![1, 0])], 1);
    assert_eq!(distribution[&(vec![1, 3], vec![0, 1])], 1);

    println!(
        "SMT_0(2,1), max entry 2, degree 4: {} tableaux",
        tableaux.len()
    );
    for ((weight, diagonal_weight), coefficient) in &distribution {
        println!(
            "  coeff {coefficient}: x^({},{}) t^({},{})",
            weight[0], weight[1], diagonal_weight[0], diagonal_weight[1]
        );
    }
}
