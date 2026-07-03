use sym_poly_core::Composition;
use sym_poly_qsym::young_qsym_schur;

fn format_composition(composition: &Composition) -> String {
    composition
        .parts()
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join("")
}

fn format_gessel_term(composition: &Composition, coefficient: i64) -> String {
    let basis = format!("\\gessel_{{{}}}", format_composition(composition));
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
        let expansion = young_qsym_schur::<i64>(alpha.parts());
        let terms = expansion
            .terms()
            .iter()
            .map(|(composition, coefficient)| format_gessel_term(composition, *coefficient))
            .collect::<Vec<_>>()
            .join(" + ");

        println!("${}$ & ${}$ \\\\", format_composition(&alpha), terms);
    }
}
