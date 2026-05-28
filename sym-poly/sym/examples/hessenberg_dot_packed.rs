use std::collections::BTreeMap;
use std::env;
use std::time::Instant;

use num_rational::Ratio;
use sym_poly_sym::{
    hessenberg_area_dot_frobenius_target, hessenberg_gkm_dot_character_values_packed_mod_prime,
    hessenberg_gkm_dot_frobenius, hessenberg_gkm_dot_frobenius_packed, SymmetricFunction,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    Compare,
    GenericOnly,
    PackedModPrimeOnly,
    PackedOnly,
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let (mode, area_args) = parse_mode(&args).unwrap_or_else(|message| {
        eprintln!("{message}");
        print_usage_and_exit();
    });
    let area = parse_area(area_args).unwrap_or_else(|message| {
        eprintln!("{message}");
        print_usage_and_exit();
    });

    println!("area {area:?}");

    if mode == Mode::PackedModPrimeOnly {
        let started = Instant::now();
        let Some(values) = hessenberg_gkm_dot_character_values_packed_mod_prime::<251>(&area)
        else {
            println!("invalid area sequence");
            return;
        };
        println!(
            "packed mod 251 dot-action traces: {:.3?}",
            started.elapsed()
        );
        for (degree, degree_values) in values {
            println!("  q^{degree}: {} nonzero class traces", degree_values.len());
        }
        return;
    }

    let target = hessenberg_area_dot_frobenius_target(&area);

    let packed_result = if mode != Mode::GenericOnly {
        let packed_started = Instant::now();
        let Some(packed) = hessenberg_gkm_dot_frobenius_packed(&area) else {
            println!("invalid area sequence");
            return;
        };
        let packed_elapsed = packed_started.elapsed();
        let packed_schur = schur_strings_ratio(&packed);
        let matches_target = target
            .as_ref()
            .is_some_and(|target| schur_strings_i64(target) == packed_schur);

        println!("packed small-prime dot action: {:.3?}", packed_elapsed);
        println!("packed matches omega X_G(q): {matches_target}");
        Some((packed_schur, packed_elapsed))
    } else {
        None
    };

    let generic_result = if mode != Mode::PackedOnly {
        let generic_started = Instant::now();
        let Some(generic) = hessenberg_gkm_dot_frobenius(&area) else {
            println!("invalid area sequence");
            return;
        };
        let generic_elapsed = generic_started.elapsed();
        let generic_schur = schur_strings_ratio(&generic);
        println!("generic modular dot action: {:.3?}", generic_elapsed);
        Some((generic_schur, generic_elapsed))
    } else {
        None
    };

    if let (Some((packed_schur, _)), Some((generic_schur, _))) = (&packed_result, &generic_result) {
        println!("packed matches generic: {}", packed_schur == generic_schur);
    }

    let display_schur = packed_result
        .as_ref()
        .map(|(schur, _)| schur)
        .or_else(|| generic_result.as_ref().map(|(schur, _)| schur))
        .expect("at least one computation was requested");
    for (degree, schur) in display_schur {
        println!("  q^{degree}: {schur}");
    }
}

fn parse_mode(args: &[String]) -> Result<(Mode, &[String]), String> {
    match args.first().map(String::as_str) {
        Some("--packed-only") => Ok((Mode::PackedOnly, &args[1..])),
        Some("--packed-mod-prime-only") => Ok((Mode::PackedModPrimeOnly, &args[1..])),
        Some("--generic-only") => Ok((Mode::GenericOnly, &args[1..])),
        Some("--compare") => Ok((Mode::Compare, &args[1..])),
        Some(flag) if flag.starts_with("--") => Err(format!("unknown flag {flag:?}")),
        _ => Ok((Mode::Compare, args)),
    }
}

fn print_usage_and_exit() -> ! {
    eprintln!("usage:");
    eprintln!("  cargo run -p sym-poly-sym --example hessenberg_dot_packed");
    eprintln!("  cargo run -p sym-poly-sym --example hessenberg_dot_packed -- --packed-only");
    eprintln!(
        "  cargo run -p sym-poly-sym --example hessenberg_dot_packed -- --packed-mod-prime-only"
    );
    eprintln!(
        "  cargo run -p sym-poly-sym --example hessenberg_dot_packed -- --packed-only 0,1,1,2,3,3"
    );
    std::process::exit(2);
}

fn parse_area(args: &[String]) -> Result<Vec<u8>, String> {
    if args.is_empty() {
        return Ok(vec![0, 1, 1, 2, 3, 3]);
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

    Ok(area)
}

fn parse_area_entry(text: &str) -> Result<u8, String> {
    text.parse::<u8>()
        .map_err(|_| format!("expected a nonnegative area entry, got {text:?}"))
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
