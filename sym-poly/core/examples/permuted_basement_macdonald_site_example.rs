use sym_poly_core::ssaf::Ssaf;

fn factor_string(arm: u32, leg: u32) -> String {
    format!("(1-t)/(1-q^{} t^{})", leg + 1, arm + 1)
}

fn factor_product_string(data: &[(u32, u32)]) -> String {
    if data.is_empty() {
        return "1".to_string();
    }

    data.iter()
        .map(|&(arm, leg)| factor_string(arm, leg))
        .collect::<Vec<_>>()
        .join(" ")
}

fn row_entries(filling: &Ssaf) -> Vec<u32> {
    filling
        .rows()
        .iter()
        .flat_map(|row| row.iter().skip(1).copied())
        .collect()
}

fn main() {
    // Moura--Mandelshtam, Example 2.13: alpha = 1101 and sigma = [2,4,1,3].
    let shape = [1, 1, 0, 1];
    let basement = [2, 4, 1, 3];

    let mut fillings = Ssaf::non_attacking_fillings(&shape, &basement);
    fillings.sort_by_key(row_entries);

    assert_eq!(fillings.len(), 4);

    let expected = [
        (vec![1, 4, 3], vec![1, 0, 1, 1], 0, 0, vec![(2, 0)]),
        (vec![2, 1, 3], vec![1, 1, 1, 0], 0, 1, vec![(1, 0)]),
        (vec![2, 4, 3], vec![0, 1, 1, 1], 0, 0, vec![]),
        (vec![4, 1, 3], vec![1, 0, 1, 1], 1, 2, vec![(2, 0), (1, 0)]),
    ];

    for (filling, (entries, weight, maj, coinv, factors)) in fillings.iter().zip(expected) {
        let factor_data = filling.arm_leg_data();

        assert_eq!(row_entries(filling), entries);
        assert_eq!(filling.weight_vector(), weight);
        assert_eq!(filling.major_index(), maj);
        assert_eq!(filling.coinversions(), coinv);
        assert_eq!(factor_data, factors);

        println!("{filling}");
        println!("  entries      = {:?}", row_entries(filling));
        println!("  weight       = {:?}", filling.weight_vector());
        println!(
            "  maj, coinv   = {}, {}",
            filling.major_index(),
            filling.coinversions()
        );
        println!("  factors      = {}", factor_product_string(&factor_data));
        println!();
    }
}
