use sym_poly_core::Composition;
use sym_poly_qsym::peak_quasisymmetric;

fn format_composition(composition: &Composition) -> String {
    composition
        .parts()
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join("")
}

fn format_peak_set(peak_set: &[u32]) -> String {
    if peak_set.is_empty() {
        return "\\emptyset".to_string();
    }
    format!(
        "\\{{{}\\}}",
        peak_set
            .iter()
            .map(u32::to_string)
            .collect::<Vec<_>>()
            .join(",")
    )
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
    let peak_sets = [vec![], vec![2], vec![3]];

    for peak_set in peak_sets {
        let expansion = peak_quasisymmetric::<i64>(&peak_set, degree);
        let terms = expansion
            .terms()
            .iter()
            .map(|(composition, coefficient)| format_gessel_term(composition, *coefficient))
            .collect::<Vec<_>>()
            .join(" + ");

        println!("${}$ & ${}$ \\\\", format_peak_set(&peak_set), terms);
    }
}
