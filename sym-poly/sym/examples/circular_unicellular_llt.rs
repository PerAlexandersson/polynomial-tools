use std::env;

use combinatoric_core::Graph;
use sym_poly_core::UnivariatePolynomial;
use sym_poly_sym::{
    circular_unicellular_llt, circular_unicellular_llt_character_values_by_degree,
    circular_unicellular_llt_frobenius_target, SymmetricFunction,
};

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let areas = parse_input(&args).unwrap_or_else(|message| {
        eprintln!("{message}");
        eprintln!("usage:");
        eprintln!("  cargo run -p sym-poly-sym --example circular_unicellular_llt");
        eprintln!("  cargo run -p sym-poly-sym --example circular_unicellular_llt -- 3");
        eprintln!("  cargo run -p sym-poly-sym --example circular_unicellular_llt -- 1,1,1");
        std::process::exit(2);
    });

    for area in areas {
        print_area_sequence(&area);
    }
}

fn parse_input(args: &[String]) -> Result<Vec<Vec<u8>>, String> {
    if args.is_empty() {
        return Ok(vec![vec![1, 1, 1]]);
    }

    if args.len() == 1 && !args[0].contains(',') {
        let n = args[0]
            .parse::<usize>()
            .map_err(|_| format!("expected a rank or an area sequence, got {:?}", args[0]))?;
        if n > u8::MAX as usize {
            return Err("rank must be at most 255 for u8 area sequences".to_string());
        }
        return Ok(all_circular_area_sequences(n));
    }

    let area = if args.len() == 1 {
        args[0]
            .split(',')
            .map(parse_area_entry)
            .collect::<Result<Vec<_>, _>>()?
    } else {
        args.iter()
            .map(|arg| parse_area_entry(arg))
            .collect::<Result<Vec<_>, _>>()?
    };

    Ok(vec![area])
}

fn parse_area_entry(text: &str) -> Result<u8, String> {
    text.parse::<u8>()
        .map_err(|_| format!("expected a nonnegative area entry, got {text:?}"))
}

fn print_area_sequence(area: &[u8]) {
    let Some(llt) = circular_unicellular_llt(area) else {
        println!(
            "area {:?}: invalid or non-symmetric circular LLT data",
            area
        );
        return;
    };
    let frobenius = circular_unicellular_llt_frobenius_target(area).unwrap();
    let characters = circular_unicellular_llt_character_values_by_degree(area).unwrap();

    println!("area {:?}", area);
    println!("Schur-positive: {}", is_schur_positive(&frobenius));
    println!("Circular unicellular LLT in m-basis:");
    println!("  {}", format_monomial_function(&llt));
    println!("Schur Frobenius:");
    for (degree, schur) in frobenius {
        println!("  q^{degree}: {schur}");
    }
    println!("Characters by cycle type:");
    for (degree, values) in characters {
        let parts = values
            .into_iter()
            .map(|(cycle_type, value)| format!("chi({cycle_type})={value}"))
            .collect::<Vec<_>>()
            .join(", ");
        println!("  q^{degree}: {parts}");
    }
    println!();
}

fn format_monomial_function(f: &SymmetricFunction<UnivariatePolynomial<i64>>) -> String {
    let monomial = f.to_monomial_basis();
    if monomial.is_zero() {
        return "0".to_string();
    }

    monomial
        .terms()
        .iter()
        .map(|(partition, coefficient)| {
            let coefficient = coefficient.to_string();
            let coefficient = if coefficient.contains(" + ") {
                format!("({coefficient})")
            } else {
                coefficient
            };
            if coefficient == "1" {
                format!("m[{partition}]")
            } else {
                format!("{coefficient}*m[{partition}]")
            }
        })
        .collect::<Vec<_>>()
        .join(" + ")
}

fn is_schur_positive(
    frobenius: &std::collections::BTreeMap<u32, sym_poly_sym::SymmetricFunction<i64>>,
) -> bool {
    frobenius
        .values()
        .all(|degree| degree.terms().values().all(|&coefficient| coefficient >= 0))
}

fn all_circular_area_sequences(n: usize) -> Vec<Vec<u8>> {
    if n == 0 {
        return vec![Vec::new()];
    }

    let mut result = Vec::new();
    let mut current = vec![0u8; n];
    circular_area_sequences_rec(0, &mut current, &mut result);
    result
}

fn circular_area_sequences_rec(index: usize, current: &mut [u8], result: &mut Vec<Vec<u8>>) {
    if index == current.len() {
        if Graph::is_circular_unit_interval_area_sequence(current) {
            result.push(current.to_vec());
        }
        return;
    }

    for value in 0..current.len() {
        current[index] = value as u8;
        circular_area_sequences_rec(index + 1, current, result);
    }
}
