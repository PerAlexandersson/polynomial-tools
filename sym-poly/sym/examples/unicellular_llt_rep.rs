use std::env;

use sym_poly_sym::{unicellular_llt_character_values_by_degree, unicellular_llt_frobenius_target};

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let areas = parse_input(&args).unwrap_or_else(|message| {
        eprintln!("{message}");
        eprintln!("usage:");
        eprintln!("  cargo run -p sym-poly-sym --example unicellular_llt_rep");
        eprintln!("  cargo run -p sym-poly-sym --example unicellular_llt_rep -- 4");
        eprintln!("  cargo run -p sym-poly-sym --example unicellular_llt_rep -- 0,1,1");
        std::process::exit(2);
    });

    for area in areas {
        print_area_sequence(&area);
    }
}

fn parse_input(args: &[String]) -> Result<Vec<Vec<u8>>, String> {
    if args.is_empty() {
        return Ok(all_area_sequences(3));
    }

    if args.len() == 1 && !args[0].contains(',') {
        let n = args[0]
            .parse::<usize>()
            .map_err(|_| format!("expected a rank or an area sequence, got {:?}", args[0]))?;
        if n > u8::MAX as usize {
            return Err("rank must be at most 255 for u8 area sequences".to_string());
        }
        return Ok(all_area_sequences(n));
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
    let Some(frobenius) = unicellular_llt_frobenius_target(area) else {
        println!("area {:?}: invalid area sequence", area);
        return;
    };
    let characters = unicellular_llt_character_values_by_degree(area).unwrap();

    println!("area {:?}", area);
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

fn all_area_sequences(n: usize) -> Vec<Vec<u8>> {
    if n == 0 {
        return vec![Vec::new()];
    }

    let mut result = Vec::new();
    let mut current = vec![0u8; n];
    area_sequences_rec(1, &mut current, &mut result);
    result
}

fn area_sequences_rec(index: usize, current: &mut [u8], result: &mut Vec<Vec<u8>>) {
    if index == current.len() {
        result.push(current.to_vec());
        return;
    }

    let upper = index.min(current[index - 1] as usize + 1);
    for value in 0..=upper {
        current[index] = value as u8;
        area_sequences_rec(index + 1, current, result);
    }
}
