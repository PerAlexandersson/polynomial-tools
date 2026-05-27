use std::collections::BTreeMap;
use std::env;

use num_rational::Ratio;
use sym_poly_sym::SymmetricFunction;
use sym_poly_sym::{hessenberg_area_dot_frobenius_target, hessenberg_gkm_dot_frobenius};

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let areas = parse_input(&args).unwrap_or_else(|message| {
        eprintln!("{message}");
        eprintln!("usage:");
        eprintln!("  cargo run -p sym-poly-sym --example hessenberg_dot");
        eprintln!("  cargo run -p sym-poly-sym --example hessenberg_dot -- 4");
        eprintln!("  cargo run -p sym-poly-sym --example hessenberg_dot -- 0,1,1");
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
    let Some(frobenius) = hessenberg_gkm_dot_frobenius(area) else {
        println!("area {:?}: invalid area sequence", area);
        return;
    };
    let target = hessenberg_area_dot_frobenius_target(area);
    let matches_target = target
        .as_ref()
        .is_some_and(|target| schur_strings_i64(target) == schur_strings_ratio(&frobenius));

    println!("area {:?}", area);
    println!("matches omega X_G(q): {matches_target}");
    for (degree, schur) in schur_strings_ratio(&frobenius) {
        println!("  q^{degree}: {schur}");
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

fn schur_strings_ratio(
    frobenius: &BTreeMap<u32, SymmetricFunction<Ratio<i64>>>,
) -> BTreeMap<u32, String> {
    frobenius
        .iter()
        .map(|(&degree, function)| (degree, function.to_schur_basis().to_string()))
        .collect()
}

fn schur_strings_i64(frobenius: &BTreeMap<u32, SymmetricFunction<i64>>) -> BTreeMap<u32, String> {
    frobenius
        .iter()
        .map(|(&degree, function)| (degree, function.to_schur_basis().to_string()))
        .collect()
}
