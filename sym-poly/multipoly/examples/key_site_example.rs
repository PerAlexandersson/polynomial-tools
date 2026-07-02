use std::collections::BTreeMap;

use sym_poly_core::{Ssaf, Tableau};
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
        String::from("1")
    } else {
        factors.join(" ")
    }
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

fn format_ssaf(ssaf: &Ssaf) -> String {
    ssaf.rows()
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

fn main() {
    let alpha = vec![1, 0, 2];

    let key_tableau = Tableau::key_tableau_from_weight(&alpha);
    assert!(key_tableau.is_key_tableau());
    println!("alpha = {}", format_weight(&alpha));
    println!("key tableau rows: {}", format_tableau(&key_tableau));

    let polynomial = key_polynomial::<i64>(&alpha);
    let fillings = Ssaf::key_fillings(&alpha);
    let mut ssaf_weight_counts = BTreeMap::new();
    for filling in &fillings {
        *ssaf_weight_counts
            .entry(filling.weight_vector())
            .or_insert(0) += 1;
    }
    assert_eq!(polynomial.terms(), &ssaf_weight_counts);

    println!("key fillings:");
    for filling in &fillings {
        println!(
            "  {} has weight {}",
            format_ssaf(filling),
            format_weight(&filling.weight_vector())
        );
    }

    println!("key polynomial:");
    for (weight, coefficient) in polynomial.terms() {
        println!("  {} {}", coefficient, format_monomial(weight));
    }
}
