use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use flagged_lorentzian::{
    crystal_e, crystal_f, enumerate_tableaux, families::skew_shapes_of_size, RowFlaggedSkewShape,
    SkewShape, TableauRecord,
};

type Content = Vec<u32>;

#[derive(Debug, Clone)]
struct Args {
    lambda: Option<Vec<u32>>,
    mu: Vec<u32>,
    alphabet: usize,
    lower_label: u32,
    max_skew_size: u32,
    max_outer_extra: u32,
    connected_only: bool,
    tableau_limit: Option<usize>,
    stop_at_first_partial: bool,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            lambda: None,
            mu: Vec::new(),
            alphabet: 5,
            lower_label: 1,
            max_skew_size: 7,
            max_outer_extra: 8,
            connected_only: false,
            tableau_limit: None,
            stop_at_first_partial: false,
        }
    }
}

fn main() {
    let args = match Args::parse_from_env() {
        Ok(args) => args,
        Err(message) => {
            eprintln!("{message}");
            eprintln!();
            eprintln!("{}", help_text());
            std::process::exit(2);
        }
    };
    let started = Instant::now();

    if let Some(lambda) = args.lambda.clone() {
        let shape = SkewShape::from_parts(lambda, args.mu.clone());
        let flagged_shape = RowFlaggedSkewShape::ordinary(shape, args.alphabet);
        let outcome = scan_one_shape(&flagged_shape, &args);
        print_single_outcome(&flagged_shape, &outcome, started.elapsed().as_secs_f64());
        std::process::exit(if outcome.limit_exceeded { 2 } else { 0 });
    }

    let mut shapes_checked = 0usize;
    let mut fibers_checked = 0usize;
    let mut total_negative_pairs = 0u128;
    let mut total_mapped_pairs = 0u128;
    let mut skipped_by_limit = 0usize;

    for skew_size in 0..=args.max_skew_size {
        let max_outer_size = skew_size + args.max_outer_extra;
        for shape in skew_shapes_of_size(skew_size, max_outer_size) {
            if args.connected_only && !shape.is_connected() {
                continue;
            }

            let flagged_shape = RowFlaggedSkewShape::ordinary(shape, args.alphabet);
            let outcome = scan_one_shape(&flagged_shape, &args);
            if outcome.limit_exceeded {
                skipped_by_limit += 1;
                continue;
            }

            shapes_checked += 1;
            fibers_checked += outcome.fibers_checked;
            total_negative_pairs += outcome.negative_pairs;
            total_mapped_pairs += outcome.mapped_pairs;

            if args.stop_at_first_partial && outcome.first_partial.is_some() {
                println!("PARTIAL");
                println!("outer={:?}", flagged_shape.shape().outer().parts());
                println!("inner={:?}", flagged_shape.shape().inner().parts());
                println!("{}", outcome.first_partial.unwrap());
                println!(
                    "shapes_checked={shapes_checked}, fibers_checked={fibers_checked}, total_negative_pairs={total_negative_pairs}, total_mapped_pairs={total_mapped_pairs}, skipped_by_limit={skipped_by_limit}, elapsed={:.3}s",
                    started.elapsed().as_secs_f64()
                );
                std::process::exit(0);
            }
        }

        println!(
            "skew_size={skew_size}: shapes_checked={shapes_checked}, fibers_checked={fibers_checked}, total_negative_pairs={total_negative_pairs}, total_mapped_pairs={total_mapped_pairs}, skipped_by_limit={skipped_by_limit}, elapsed={:.3}s",
            started.elapsed().as_secs_f64()
        );
    }

    println!("DONE");
    println!(
        "shapes_checked={shapes_checked}, fibers_checked={fibers_checked}, total_negative_pairs={total_negative_pairs}, total_mapped_pairs={total_mapped_pairs}, skipped_by_limit={skipped_by_limit}, elapsed={:.3}s",
        started.elapsed().as_secs_f64()
    );
}

impl Args {
    fn parse_from_env() -> Result<Self, String> {
        let mut args = Args::default();
        let mut iter = std::env::args().skip(1);
        while let Some(flag) = iter.next() {
            match flag.as_str() {
                "--help" | "-h" => {
                    println!("{}", help_text());
                    std::process::exit(0);
                }
                "--lambda" => args.lambda = Some(parse_u32_vec(&take_value(&mut iter, &flag)?)?),
                "--mu" => args.mu = parse_u32_vec(&take_value(&mut iter, &flag)?)?,
                "--alphabet" => {
                    args.alphabet = take_value(&mut iter, &flag)?
                        .parse()
                        .map_err(|err| format!("invalid --alphabet: {err}"))?
                }
                "--lower-label" => {
                    args.lower_label = take_value(&mut iter, &flag)?
                        .parse()
                        .map_err(|err| format!("invalid --lower-label: {err}"))?
                }
                "--max-skew-size" => {
                    args.max_skew_size = take_value(&mut iter, &flag)?
                        .parse()
                        .map_err(|err| format!("invalid --max-skew-size: {err}"))?
                }
                "--max-outer-extra" => {
                    args.max_outer_extra = take_value(&mut iter, &flag)?
                        .parse()
                        .map_err(|err| format!("invalid --max-outer-extra: {err}"))?
                }
                "--connected-only" => args.connected_only = true,
                "--tableau-limit" => {
                    args.tableau_limit = Some(
                        take_value(&mut iter, &flag)?
                            .parse()
                            .map_err(|err| format!("invalid --tableau-limit: {err}"))?,
                    )
                }
                "--stop-at-first-partial" => args.stop_at_first_partial = true,
                other => return Err(format!("unknown argument `{other}`")),
            }
        }
        if args.lower_label as usize >= args.alphabet {
            return Err("--lower-label must be smaller than --alphabet".to_string());
        }
        Ok(args)
    }
}

#[derive(Debug)]
struct ShapeScanOutcome {
    tableaux: usize,
    fibers_checked: usize,
    negative_pairs: u128,
    positive_pairs: u128,
    mapped_pairs: u128,
    collision_pairs: u128,
    first_partial: Option<PartialFiber>,
    limit_exceeded: bool,
}

#[derive(Debug, Clone)]
struct PartialFiber {
    beta: Content,
    negative_pairs: u128,
    positive_pairs: u128,
    mapped_pairs: u128,
    blocked_first: u128,
    blocked_second: u128,
    blocked_both: u128,
}

impl std::fmt::Display for PartialFiber {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "beta={:?}, negative_pairs={}, positive_pairs={}, mapped_pairs={}, blocked_first={}, blocked_second={}, blocked_both={}",
            self.beta,
            self.negative_pairs,
            self.positive_pairs,
            self.mapped_pairs,
            self.blocked_first,
            self.blocked_second,
            self.blocked_both
        )
    }
}

fn scan_one_shape(flagged_shape: &RowFlaggedSkewShape, args: &Args) -> ShapeScanOutcome {
    let tableaux = match enumerate_tableaux(flagged_shape, args.tableau_limit) {
        Ok(tableaux) => tableaux,
        Err(_) => {
            return ShapeScanOutcome {
                tableaux: 0,
                fibers_checked: 0,
                negative_pairs: 0,
                positive_pairs: 0,
                mapped_pairs: 0,
                collision_pairs: 0,
                first_partial: None,
                limit_exceeded: true,
            }
        }
    };

    let mut by_content = BTreeMap::<Content, Vec<usize>>::new();
    for (idx, tableau) in tableaux.iter().enumerate() {
        by_content
            .entry(tableau.content.clone())
            .or_default()
            .push(idx);
    }

    let mut seen_beta = BTreeSet::new();
    let mut fibers_checked = 0usize;
    let mut negative_pairs = 0u128;
    let mut positive_pairs = 0u128;
    let mut mapped_pairs = 0u128;
    let mut collision_pairs = 0u128;
    let mut first_partial = None;

    let lower = args.lower_label as usize - 1;
    let upper = lower + 1;

    for content in by_content.keys() {
        let Some(beta) = subtract_units(content, lower, 2) else {
            continue;
        };
        let ii = add_units(&beta, lower, 2);
        let mixed = add_units(&add_units(&beta, lower, 1), upper, 1);
        let jj = add_units(&beta, upper, 2);
        if content != &ii || !seen_beta.insert(beta.clone()) {
            continue;
        }

        let (Some(ii_indices), Some(mixed_indices), Some(jj_indices)) = (
            by_content.get(&ii),
            by_content.get(&mixed),
            by_content.get(&jj),
        ) else {
            continue;
        };

        fibers_checked += 1;
        let partial = scan_fiber(
            flagged_shape.shape(),
            &tableaux,
            &beta,
            ii_indices,
            mixed_indices,
            jj_indices,
            args.lower_label,
        );
        negative_pairs += partial.negative_pairs;
        positive_pairs += partial.positive_pairs;
        mapped_pairs += partial.mapped_pairs;
        collision_pairs += partial.mapped_pairs.saturating_sub(unique_image_count(
            flagged_shape.shape(),
            &tableaux,
            ii_indices,
            jj_indices,
            args.lower_label,
        ));

        if partial.mapped_pairs < partial.negative_pairs && first_partial.is_none() {
            first_partial = Some(partial);
        }
    }

    ShapeScanOutcome {
        tableaux: tableaux.len(),
        fibers_checked,
        negative_pairs,
        positive_pairs,
        mapped_pairs,
        collision_pairs,
        first_partial,
        limit_exceeded: false,
    }
}

fn scan_fiber(
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    beta: &[u32],
    ii_indices: &[usize],
    mixed_indices: &[usize],
    jj_indices: &[usize],
    lower_label: u32,
) -> PartialFiber {
    let mixed_values: BTreeSet<_> = mixed_indices
        .iter()
        .map(|&idx| tableaux[idx].values.clone())
        .collect();
    let mut mapped_pairs = 0u128;
    let mut blocked_first = 0u128;
    let mut blocked_second = 0u128;
    let mut blocked_both = 0u128;

    for &left_idx in ii_indices {
        let left_image = crystal_f(shape, &tableaux[left_idx].values, lower_label);
        for &right_idx in jj_indices {
            let right_image = crystal_e(shape, &tableaux[right_idx].values, lower_label);
            match (&left_image, &right_image) {
                (Some(left_values), Some(right_values))
                    if mixed_values.contains(left_values)
                        && mixed_values.contains(right_values) =>
                {
                    mapped_pairs += 1;
                }
                (None, None) => blocked_both += 1,
                (None, _) => blocked_first += 1,
                (_, None) => blocked_second += 1,
                _ => {}
            }
        }
    }

    PartialFiber {
        beta: beta.to_vec(),
        negative_pairs: ii_indices.len() as u128 * jj_indices.len() as u128,
        positive_pairs: mixed_indices.len() as u128 * mixed_indices.len() as u128,
        mapped_pairs,
        blocked_first,
        blocked_second,
        blocked_both,
    }
}

fn unique_image_count(
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    ii_indices: &[usize],
    jj_indices: &[usize],
    lower_label: u32,
) -> u128 {
    let mut images = BTreeSet::new();
    for &left_idx in ii_indices {
        let left_image = crystal_f(shape, &tableaux[left_idx].values, lower_label);
        for &right_idx in jj_indices {
            let right_image = crystal_e(shape, &tableaux[right_idx].values, lower_label);
            if let (Some(left_values), Some(right_values)) = (&left_image, &right_image) {
                images.insert((left_values.clone(), right_values.clone()));
            }
        }
    }
    images.len() as u128
}

fn subtract_units(content: &[u32], index: usize, amount: u32) -> Option<Content> {
    let mut out = content.to_vec();
    *out.get_mut(index)? = out[index].checked_sub(amount)?;
    Some(out)
}

fn add_units(content: &[u32], index: usize, amount: u32) -> Content {
    let mut out = content.to_vec();
    out[index] += amount;
    out
}

fn print_single_outcome(
    flagged_shape: &RowFlaggedSkewShape,
    outcome: &ShapeScanOutcome,
    elapsed_seconds: f64,
) {
    if outcome.limit_exceeded {
        println!("SKIPPED: tableau limit exceeded");
        return;
    }

    println!("outer={:?}", flagged_shape.shape().outer().parts());
    println!("inner={:?}", flagged_shape.shape().inner().parts());
    println!("alphabet={}", flagged_shape.alphabet());
    println!("tableaux={}", outcome.tableaux);
    println!("fibers_checked={}", outcome.fibers_checked);
    println!("negative_pairs={}", outcome.negative_pairs);
    println!("positive_pairs={}", outcome.positive_pairs);
    println!("mapped_pairs={}", outcome.mapped_pairs);
    println!("collision_pairs={}", outcome.collision_pairs);
    println!("elapsed={elapsed_seconds:.3}s");
    if let Some(partial) = &outcome.first_partial {
        println!("first_partial={partial}");
    } else {
        println!("all_negative_pairs_mapped_by_f_tensor_e");
    }
}

fn take_value(iter: &mut impl Iterator<Item = String>, flag: &str) -> Result<String, String> {
    iter.next()
        .ok_or_else(|| format!("missing value after {flag}"))
}

fn parse_u32_vec(input: &str) -> Result<Vec<u32>, String> {
    if input.trim().is_empty() {
        return Ok(Vec::new());
    }

    input
        .split(',')
        .map(|part| {
            part.trim()
                .parse::<u32>()
                .map_err(|err| format!("invalid integer `{part}`: {err}"))
        })
        .collect()
}

fn help_text() -> &'static str {
    "Scan how much of the 2x2 injection is covered by (A,C) -> (f_i A, e_i C).

USAGE:
  crystal_transfer_scan [OPTIONS]

OPTIONS:
  --lambda PARTS             Outer partition, e.g. 4,3,1. If omitted, scan a family.
  --mu PARTS                 Inner partition, e.g. 3,1. Default: empty.
  --alphabet N               Alphabet size. Default: 5.
  --lower-label I            Use crystal operators on labels I,I+1. Default: 1.
  --max-skew-size N          Maximum skew size for family scans. Default: 7.
  --max-outer-extra N        Family outer sizes go up to skew_size + this. Default: 8.
  --connected-only           Restrict family scans to connected skew shapes.
  --tableau-limit N          Skip any shape whose tableau enumeration exceeds N.
  --stop-at-first-partial    Stop as soon as f_i x e_i is not total on a fiber.
  --help                     Print this help.

EXAMPLE:
  cargo run -q -p flagged-lorentzian --bin crystal_transfer_scan -- \\
    --lambda 3,2 --alphabet 5"
}
