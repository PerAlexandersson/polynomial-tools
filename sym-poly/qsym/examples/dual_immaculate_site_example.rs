use sym_poly_core::Composition;
use sym_poly_qsym::{dual_immaculate, ImmaculateTableau};

fn format_composition(composition: &Composition) -> String {
    composition
        .parts()
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join("")
}

fn format_tableau(tableau: &ImmaculateTableau) -> String {
    tableau
        .rows
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
    let alpha = vec![2, 1];
    let tableaux = ImmaculateTableau::enumerate(&alpha);
    let expansion = dual_immaculate::<i64>(&alpha);

    assert_eq!(tableaux.len(), 2);
    assert_eq!(expansion.coefficient(&Composition::new(vec![1, 2])), 1);
    assert_eq!(expansion.coefficient(&Composition::new(vec![2, 1])), 1);
    assert_eq!(expansion.terms().len(), 2);

    println!("shape alpha = (2,1)");
    println!("standard immaculate tableaux:");
    for tableau in &tableaux {
        println!(
            "  {}: Des = {:?}, composition = {}",
            format_tableau(tableau),
            tableau.descent_set(),
            format_composition(&tableau.descent_composition())
        );
    }

    println!("dual immaculate expansion:");
    for (composition, coefficient) in expansion.terms() {
        println!("  {} F_{}", coefficient, format_composition(composition));
    }
}
