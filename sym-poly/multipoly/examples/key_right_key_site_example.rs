use std::collections::BTreeMap;

use sym_poly_core::{Partition, Tableau};
use sym_poly_multipoly::key_polynomial;

fn format_weight(weight: &[u32]) -> String {
    format!(
        "({})",
        weight
            .iter()
            .map(u32::to_string)
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn format_tableau(tableau: &Tableau) -> String {
    tableau
        .rows()
        .iter()
        .map(|row| {
            format!(
                "[{}]",
                row.iter().map(u32::to_string).collect::<Vec<_>>().join(",")
            )
        })
        .collect::<Vec<_>>()
        .join(" / ")
}

fn tableau_weight(tableau: &Tableau, n: usize) -> Vec<u32> {
    let mut weight = vec![0; n];
    for &entry in tableau.rows().iter().flatten() {
        if entry as usize <= n {
            weight[entry as usize - 1] += 1;
        }
    }
    weight
}

fn tableau_entrywise_leq(left: &Tableau, right: &Tableau) -> bool {
    left.rows().len() == right.rows().len()
        && left.rows().iter().zip(right.rows()).all(|(lrow, rrow)| {
            lrow.len() == rrow.len() && lrow.iter().zip(rrow).all(|(l, r)| l <= r)
        })
}

fn main() {
    let alpha = vec![1, 0, 2];
    let shape = Partition::from_sorted(vec![2, 1]);
    let key_bound = Tableau::key_tableau_from_weight(&alpha);
    let key_poly = key_polynomial::<i64>(&alpha);

    let tableaux = Tableau::semistandard_tableaux(&shape, alpha.len() as u32);
    let mut accepted_weight_counts = BTreeMap::new();

    println!("alpha = {}", format_weight(&alpha));
    println!("key(alpha) = {}", format_tableau(&key_bound));
    println!("SSYT of shape (2,1), entries at most 3:");

    for tableau in tableaux {
        let right_key_weight = tableau.right_key_weight_via_ssaf(alpha.len());
        let right_key = tableau.right_key_via_ssaf(alpha.len());
        let accepted = tableau_entrywise_leq(&right_key, &key_bound);
        if accepted {
            *accepted_weight_counts
                .entry(tableau_weight(&tableau, alpha.len()))
                .or_insert(0) += 1;
        }
        println!(
            "  {} has K_+ weight {} ({})",
            format_tableau(&tableau),
            format_weight(right_key_weight.parts()),
            if accepted { "included" } else { "excluded" }
        );
    }

    assert_eq!(accepted_weight_counts, *key_poly.terms());
}
