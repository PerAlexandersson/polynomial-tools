use combinatoric_core::Partition;
use sym_poly_sym::{shifted_multiset_tableau_distribution, ShiftedMultisetTableau};

fn partition(parts: &[u32]) -> Partition {
    Partition::new(parts.to_vec())
}

fn main() {
    let shape = partition(&[2, 1]);
    let degree_three_tableaux = ShiftedMultisetTableau::enumerate_smt0(&shape, 2, 3);
    let degree_three = shifted_multiset_tableau_distribution(&shape, 2, 3);
    let degree_four_tableaux = ShiftedMultisetTableau::enumerate_smt0(&shape, 2, 4);
    let degree_four = shifted_multiset_tableau_distribution(&shape, 2, 4);
    let degree_five_tableaux = ShiftedMultisetTableau::enumerate_smt0(&shape, 2, 5);
    let degree_five = shifted_multiset_tableau_distribution(&shape, 2, 5);

    assert_eq!(degree_three_tableaux.len(), 2);
    assert_eq!(degree_three[&(vec![2, 1], vec![0, 0])], 1);
    assert_eq!(degree_three[&(vec![1, 2], vec![0, 0])], 1);

    assert_eq!(degree_four_tableaux.len(), 8);
    assert_eq!(degree_four[&(vec![3, 1], vec![1, 0])], 1);
    assert_eq!(degree_four[&(vec![3, 1], vec![0, 1])], 1);
    assert_eq!(degree_four[&(vec![2, 2], vec![1, 0])], 2);
    assert_eq!(degree_four[&(vec![2, 2], vec![0, 1])], 2);
    assert_eq!(degree_four[&(vec![1, 3], vec![1, 0])], 1);
    assert_eq!(degree_four[&(vec![1, 3], vec![0, 1])], 1);

    assert_eq!(degree_five_tableaux.len(), 20);
    assert_eq!(degree_five[&(vec![4, 1], vec![2, 0])], 1);
    assert_eq!(degree_five[&(vec![4, 1], vec![1, 1])], 1);
    assert_eq!(degree_five[&(vec![4, 1], vec![0, 2])], 1);
    assert_eq!(degree_five[&(vec![3, 2], vec![2, 0])], 2);
    assert_eq!(degree_five[&(vec![3, 2], vec![1, 1])], 3);
    assert_eq!(degree_five[&(vec![3, 2], vec![0, 2])], 2);
    assert_eq!(degree_five[&(vec![2, 3], vec![2, 0])], 2);
    assert_eq!(degree_five[&(vec![2, 3], vec![1, 1])], 3);
    assert_eq!(degree_five[&(vec![2, 3], vec![0, 2])], 2);
    assert_eq!(degree_five[&(vec![1, 4], vec![2, 0])], 1);
    assert_eq!(degree_five[&(vec![1, 4], vec![1, 1])], 1);
    assert_eq!(degree_five[&(vec![1, 4], vec![0, 2])], 1);

    println!(
        "SMT_0(2,1), max entry 2: degrees 3, 4, 5 have {}, {}, {} tableaux",
        degree_three_tableaux.len(),
        degree_four_tableaux.len(),
        degree_five_tableaux.len()
    );
    for ((weight, diagonal_weight), coefficient) in &degree_four {
        println!(
            "  degree 4 coeff {coefficient}: x^({},{}) t^({},{})",
            weight[0], weight[1], diagonal_weight[0], diagonal_weight[1]
        );
    }
}
