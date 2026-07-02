use sym_poly_core::Composition;
use sym_poly_qsym::{p_partition_generating_function, p_partition_linear_extensions};

fn format_composition(composition: &Composition) -> String {
    composition
        .parts()
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join("")
}

fn format_extension(extension: &[usize]) -> String {
    extension
        .iter()
        .map(|v| (v + 1).to_string())
        .collect::<Vec<_>>()
        .join("")
}

fn main() {
    // The V-poset with 1 < 3 and 2 < 3.
    let covers = [(0, 2), (1, 2)];
    let data = p_partition_linear_extensions(3, &covers);

    println!("covers: 1 < 3, 2 < 3");
    println!("linear extensions:");
    for row in &data {
        let descent_set = if row.descent_set.is_empty() {
            String::from("{}")
        } else {
            format!(
                "{{{}}}",
                row.descent_set
                    .iter()
                    .map(u32::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            )
        };
        println!(
            "  {}: Des = {}, composition = {}",
            format_extension(&row.extension),
            descent_set,
            format_composition(&row.descent_composition)
        );
    }

    let gamma = p_partition_generating_function::<i64>(3, &covers);
    println!("fundamental expansion:");
    for (composition, coefficient) in gamma.terms() {
        println!("  {} F_{}", coefficient, format_composition(composition));
    }
}
