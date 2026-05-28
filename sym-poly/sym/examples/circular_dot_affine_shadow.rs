use std::collections::BTreeMap;
use std::env;

use combinatoric_core::Graph;
use num_traits::ToPrimitive;
use sym_poly_sym::{affine_shadow_circular_gkm_hilbert, circular_area_dot_frobenius_target};

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let areas = parse_input(&args).unwrap_or_else(|message| {
        eprintln!("{message}");
        eprintln!("usage:");
        eprintln!("  cargo run -p sym-poly-sym --example circular_dot_affine_shadow -- 1,1,1");
        eprintln!("  cargo run -p sym-poly-sym --example circular_dot_affine_shadow -- 3");
        std::process::exit(2);
    });

    let mut matches = 0usize;
    for area in &areas {
        if print_area_sequence(area) {
            matches += 1;
        }
    }
    if areas.len() > 1 {
        println!("summary: {matches}/{} Hilbert series matched", areas.len());
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

fn print_area_sequence(area: &[u8]) -> bool {
    let Some(hilbert) = affine_shadow_circular_gkm_hilbert(area, 2) else {
        println!("area {:?}: invalid circular area sequence", area);
        return false;
    };
    let Some(target) = circular_area_dot_frobenius_target(area) else {
        println!("area {:?}: could not compute target", area);
        return false;
    };
    let target_hilbert = frobenius_hilbert(&target);
    let matches_target = hilbert == target_hilbert;

    println!("area {:?}", area);
    println!("affine-shadow Hilbert matches target dimensions: {matches_target}");
    println!("affine shadow: {}", format_hilbert(&hilbert));
    if !matches_target {
        println!("target:        {}", format_hilbert(&target_hilbert));
    }
    println!();

    matches_target
}

fn frobenius_hilbert(
    frobenius: &BTreeMap<u32, sym_poly_sym::SymmetricFunction<i64>>,
) -> BTreeMap<u32, usize> {
    frobenius
        .iter()
        .filter_map(|(&degree, function)| {
            let schur = function.to_schur_basis();
            let dimension = schur
                .terms()
                .iter()
                .map(|(partition, &multiplicity)| {
                    let specht_dimension = partition
                        .count_syt()
                        .to_usize()
                        .expect("small-rank Specht dimension fits in usize");
                    usize::try_from(multiplicity).expect("multiplicity is nonnegative")
                        * specht_dimension
                })
                .sum::<usize>();
            (dimension != 0).then_some((degree, dimension))
        })
        .collect()
}

fn format_hilbert(hilbert: &BTreeMap<u32, usize>) -> String {
    if hilbert.is_empty() {
        return "0".to_string();
    }

    hilbert
        .iter()
        .map(|(&degree, &coefficient)| match (coefficient, degree) {
            (1, 0) => "1".to_string(),
            (c, 0) => c.to_string(),
            (1, 1) => "q".to_string(),
            (c, 1) => format!("{c}q"),
            (1, d) => format!("q^{d}"),
            (c, d) => format!("{c}q^{d}"),
        })
        .collect::<Vec<_>>()
        .join(" + ")
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
