use std::collections::{BTreeMap, BTreeSet};

use sym_poly_multipoly::{ghost_diagram_weight, k_kohnert_weight_counts, Cell, GhostDiagram};

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

fn cells(cells: &[(usize, usize)]) -> BTreeSet<Cell> {
    cells.iter().map(|&(col, row)| Cell { col, row }).collect()
}

fn diagram_pair_weight(real: &BTreeSet<Cell>, extra: &BTreeSet<Cell>) -> Vec<u32> {
    ghost_diagram_weight(&GhostDiagram::new(real.clone(), extra.clone()))
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

    // Pan--Yu, Examples 3.2 and 3.5, for the same alpha.
    let kohnert_cells = cells(&[(1, 1), (1, 2), (2, 1)]);
    let ghost_cells = cells(&[(1, 3), (2, 2)]);
    let top_k_kohnert = GhostDiagram::new(kohnert_cells, ghost_cells);
    assert_eq!(top_k_kohnert.ghost_count(), 2);
    assert_eq!(ghost_diagram_weight(&top_k_kohnert), vec![2, 2, 1]);

    let leading_cells = cells(&[(1, 2), (1, 3), (2, 2)]);
    let extra_cells = cells(&[(1, 1), (2, 1)]);
    assert_eq!(extra_cells.len(), 2);
    assert_eq!(
        diagram_pair_weight(&leading_cells, &extra_cells),
        vec![2, 2, 1]
    );

    println!(
        "Pan--Yu bijection example: beta^2 * {} on both sides",
        format_monomial(&[2, 2, 1])
    );
}
