use sym_poly_core::Composition;
use sym_poly_qsym::{
    combinatorial_power_sum_in_monomial_basis, reverse_combinatorial_power_sum_in_monomial_basis,
};

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

fn print_table(title: &str, reverse: bool) {
    println!("{title}:");
    for alpha in Composition::integer_compositions(4) {
        let expansion = if reverse {
            reverse_combinatorial_power_sum_in_monomial_basis::<i64>(&alpha)
        } else {
            combinatorial_power_sum_in_monomial_basis::<i64>(&alpha)
        };
        let terms = expansion
            .terms()
            .iter()
            .map(|(composition, coefficient)| format_term("\\qmonom", composition, *coefficient))
            .collect::<Vec<_>>()
            .join(" + ");

        println!("${}$ & ${}$ \\\\", format_composition(&alpha), terms);
    }
}

fn main() {
    print_table("Combinatorial power sums p_alpha", false);
    println!();
    print_table("Reverse combinatorial power sums p^r_alpha", true);
}
