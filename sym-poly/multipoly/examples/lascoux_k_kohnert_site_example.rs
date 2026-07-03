use std::collections::BTreeMap;

use sym_poly_multipoly::k_kohnert_weight_counts;

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
    let alpha = [0, 2, 1];
    let counts = k_kohnert_weight_counts(&alpha, 100).expect("small example should fit cap");
    let expected = BTreeMap::from([
        ((0, vec![0, 2, 1]), 1),
        ((0, vec![1, 1, 1]), 1),
        ((0, vec![1, 2]), 1),
        ((0, vec![2, 0, 1]), 1),
        ((0, vec![2, 1]), 1),
        ((1, vec![1, 2, 1]), 2),
        ((1, vec![2, 1, 1]), 2),
        ((1, vec![2, 2]), 1),
        ((2, vec![2, 2, 1]), 1),
    ]);
    assert_eq!(counts, expected);

    println!("K-Kohnert weights for alpha=(0,2,1):");
    for ((ghosts, weight), count) in counts {
        println!(
            "  beta^{} * {} * {}",
            ghosts,
            count,
            format_monomial(&weight)
        );
    }
}
