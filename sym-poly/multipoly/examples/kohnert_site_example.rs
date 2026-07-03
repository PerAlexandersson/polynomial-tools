use std::collections::BTreeMap;

use sym_poly_multipoly::{kohnert_diagrams_for_composition, kohnert_weight_counts};

fn format_monomial(weight: &[u32]) -> String {
    let factors = weight
        .iter()
        .enumerate()
        .filter_map(|(idx, &exp)| match exp {
            0 => None,
            1 => Some(format!("x_{}", idx + 1)),
            _ => Some(format!("x_{}^{}", idx + 1, exp)),
        })
        .collect::<Vec<_>>();
    if factors.is_empty() {
        "1".to_string()
    } else {
        factors.join(" ")
    }
}

fn main() {
    let alpha = [0, 2];
    let diagrams =
        kohnert_diagrams_for_composition(&alpha, 10).expect("small Kohnert example should fit cap");
    assert_eq!(diagrams.len(), 3);

    let counts = kohnert_weight_counts(&alpha, 10).expect("small Kohnert example should fit cap");
    let expected = BTreeMap::from([(vec![0, 2], 1), (vec![1, 1], 1), (vec![2], 1)]);
    assert_eq!(counts, expected);

    println!("Kohnert weights for alpha=(0,2):");
    for (weight, count) in counts {
        println!("  {} * {}", count, format_monomial(&weight));
    }
}
