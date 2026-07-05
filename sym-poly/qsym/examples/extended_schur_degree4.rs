use sym_poly_core::Composition;
use sym_poly_qsym::{extended_schur, extended_schur_monomial};

fn format_composition(composition: &Composition) -> String {
    composition
        .parts()
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join("")
}

fn format_term(prefix: &str, composition: &Composition, coefficient: i64) -> String {
    let basis = format!("{prefix}_{{{}}}", format_composition(composition));
    match coefficient {
        1 => basis,
        -1 => format!("-{basis}"),
        c => format!("{c}{basis}"),
    }
}

fn main() {
    let degree = 4;
    let compositions = Composition::integer_compositions(degree);
    assert_eq!(compositions.len(), 8);

    println!("Monomial expansions:");
    for alpha in &compositions {
        let expansion = extended_schur_monomial::<i64>(alpha.parts());
        let terms = expansion
            .terms()
            .iter()
            .map(|(composition, coefficient)| format_term("\\qmonom", composition, *coefficient))
            .collect::<Vec<_>>()
            .join(" + ");

        println!("${}$ & ${}$ \\\\", format_composition(alpha), terms);
    }

    println!();
    println!("Fundamental expansions:");
    for alpha in compositions {
        let expansion = extended_schur::<i64>(alpha.parts());
        let terms = expansion
            .terms()
            .iter()
            .map(|(composition, coefficient)| format_term("\\gessel", composition, *coefficient))
            .collect::<Vec<_>>()
            .join(" + ");

        println!("${}$ & ${}$ \\\\", format_composition(&alpha), terms);
    }
}
