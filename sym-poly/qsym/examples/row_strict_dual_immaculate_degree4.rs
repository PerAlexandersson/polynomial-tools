use sym_poly_core::Composition;
use sym_poly_qsym::row_strict_dual_immaculate;

fn format_composition(composition: &Composition) -> String {
    composition
        .parts()
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join("")
}

fn format_qmonom_term(composition: &Composition, coefficient: i64) -> String {
    let basis = format!("\\qmonom_{{{}}}", format_composition(composition));
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

    for alpha in compositions {
        let expansion = row_strict_dual_immaculate::<i64>(alpha.parts(), degree);
        let terms = expansion
            .terms()
            .iter()
            .map(|(composition, coefficient)| format_qmonom_term(composition, *coefficient))
            .collect::<Vec<_>>()
            .join(" + ");

        println!("${}$ & ${}$ \\\\", format_composition(&alpha), terms);
    }
}
