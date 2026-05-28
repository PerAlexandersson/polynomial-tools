use std::collections::BTreeMap;
use std::env;

use combinatoric_core::Graph;
use num_rational::Ratio;
use sym_poly_sym::{
    circular_area_dot_frobenius_target, naive_circular_gkm_dot_frobenius, SymmetricFunction,
};

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let (summary_only, input_args) = parse_mode(&args);
    let areas = parse_input(input_args).unwrap_or_else(|message| {
        eprintln!("{message}");
        eprintln!("usage:");
        eprintln!("  cargo run -p sym-poly-sym --example circular_dot_naive");
        eprintln!("  cargo run -p sym-poly-sym --example circular_dot_naive -- 3");
        eprintln!("  cargo run -p sym-poly-sym --example circular_dot_naive -- --summary 4");
        eprintln!("  cargo run -p sym-poly-sym --example circular_dot_naive -- 1,1,1");
        std::process::exit(2);
    });

    let mut matches = 0usize;
    let mut mismatches = Vec::new();
    let mut naive_gkm_cache: BTreeMap<Vec<(usize, usize)>, Option<BTreeMap<u32, String>>> =
        BTreeMap::new();
    for area in &areas {
        let matches_target = if summary_only {
            area_matches_target_cached(area, &mut naive_gkm_cache)
        } else {
            print_area_sequence(area)
        };
        if matches_target {
            matches += 1;
        } else {
            mismatches.push(area.clone());
        }
    }

    if summary_only || areas.len() > 1 {
        println!("summary: {matches}/{} matched", areas.len());
        if summary_only && !mismatches.is_empty() {
            println!("mismatches:");
            for area in mismatches {
                println!("  {:?}", area);
            }
        }
    }
}

fn parse_mode(args: &[String]) -> (bool, &[String]) {
    if args.first().is_some_and(|arg| arg == "--summary") {
        (true, &args[1..])
    } else {
        (false, args)
    }
}

fn parse_input(args: &[String]) -> Result<Vec<Vec<u8>>, String> {
    if args.is_empty() {
        return Ok(all_circular_area_sequences(3));
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
    if !Graph::is_circular_unit_interval_area_sequence(area) {
        println!("area {:?}: invalid circular area sequence", area);
        return false;
    }

    let Some((computed_schur, target_schur, matches_target)) = area_comparison(area) else {
        println!("area {:?}: naive GKM computation failed", area);
        return false;
    };

    println!("area {:?}", area);
    println!("matches omega X_Gamma(q): {matches_target}");
    println!("naive GKM:");
    for (degree, schur) in &computed_schur {
        println!("  q^{degree}: {schur}");
    }
    if !matches_target {
        println!("target omega X_Gamma(q):");
        for (degree, schur) in &target_schur {
            println!("  q^{degree}: {schur}");
        }
    }
    println!();

    matches_target
}

fn area_matches_target_cached(
    area: &[u8],
    naive_gkm_cache: &mut BTreeMap<Vec<(usize, usize)>, Option<BTreeMap<u32, String>>>,
) -> bool {
    if !Graph::is_circular_unit_interval_area_sequence(area) {
        return false;
    }

    let Some(graph) = Graph::circular_unit_interval(area) else {
        return false;
    };
    let computed_schur = naive_gkm_cache
        .entry(graph.edges().to_vec())
        .or_insert_with(|| {
            naive_circular_gkm_dot_frobenius(area).map(|computed| schur_strings_ratio(&computed))
        })
        .clone();
    let Some(computed_schur) = computed_schur else {
        return false;
    };
    let Some(target) = circular_area_dot_frobenius_target(area) else {
        return false;
    };
    computed_schur == schur_strings_i64(&target)
}

fn area_comparison(area: &[u8]) -> Option<(BTreeMap<u32, String>, BTreeMap<u32, String>, bool)> {
    if !Graph::is_circular_unit_interval_area_sequence(area) {
        return None;
    }
    let computed = naive_circular_gkm_dot_frobenius(area)?;
    let target = circular_area_dot_frobenius_target(area)?;
    let computed_schur = schur_strings_ratio(&computed);
    let target_schur = schur_strings_i64(&target);
    let matches_target = computed_schur == target_schur;
    Some((computed_schur, target_schur, matches_target))
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
