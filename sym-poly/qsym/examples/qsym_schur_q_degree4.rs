use sym_poly_core::Composition;
use sym_poly_qsym::qsym_schur_q_peak_expansion;

fn format_composition(composition: &Composition) -> String {
    composition
        .parts()
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join("")
}

fn format_peak_term(composition: &Composition, coefficient: i64) -> String {
    let basis = format!("K_{{{}}}", format_composition(composition));
    match coefficient {
        1 => basis,
        -1 => format!("-{basis}"),
        c => format!("{c}{basis}"),
    }
}

fn main() {
    let peak_compositions = [
        Composition::new(vec![4]),
        Composition::new(vec![2, 2]),
        Composition::new(vec![3, 1]),
    ];

    for alpha in peak_compositions {
        let expansion = qsym_schur_q_peak_expansion::<i64>(alpha.parts());
        let terms = expansion
            .iter()
            .map(|(composition, coefficient)| format_peak_term(composition, *coefficient))
            .collect::<Vec<_>>()
            .join(" + ");

        println!("${}$ & ${}$ \\\\", format_composition(&alpha), terms);
    }
}
