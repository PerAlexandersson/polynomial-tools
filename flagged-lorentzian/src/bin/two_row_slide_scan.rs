use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use flagged_lorentzian::{
    active_subword_descent_data_for_values, descent_data_for_values, enumerate_tableaux,
    families::skew_shapes_of_size, is_gt_array, DescentData, DescentStatistic, RowFlaggedSkewShape,
    SkewGtPattern, SkewShape,
};

type Content = Vec<u32>;
type GtRows = Vec<Vec<u32>>;
type DescentPair = (DescentData, DescentData);

#[derive(Debug, Clone)]
struct Args {
    lambda: Option<Vec<u32>>,
    mu: Vec<u32>,
    row_flags: Option<Vec<u32>>,
    alphabet: usize,
    lower_label: u32,
    max_skew_size: u32,
    max_outer_extra: u32,
    exactly_two_rows: bool,
    max_beta_active: Option<u32>,
    check_injective: bool,
    check_matching: bool,
    check_greedy: bool,
    check_descent_output: bool,
    check_target_intervals: bool,
    check_target_swap_closed: bool,
    check_unordered_target_intervals: bool,
    check_augment_target_increase: bool,
    check_additive_rank: bool,
    check_layered_additive_rank: bool,
    check_carrier_word_additive_rank: bool,
    source_order: SourceOrder,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
    descent_mode: DescentMode,
    augment_potential: TargetPotential,
    tableau_limit: Option<usize>,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            lambda: None,
            mu: Vec::new(),
            row_flags: None,
            alphabet: 5,
            lower_label: 4,
            max_skew_size: 5,
            max_outer_extra: 8,
            exactly_two_rows: false,
            max_beta_active: None,
            check_injective: false,
            check_matching: false,
            check_greedy: false,
            check_descent_output: false,
            check_target_intervals: false,
            check_target_swap_closed: false,
            check_unordered_target_intervals: false,
            check_augment_target_increase: false,
            check_additive_rank: false,
            check_layered_additive_rank: false,
            check_carrier_word_additive_rank: false,
            source_order: SourceOrder::Lex,
            carrier_order: CarrierOrder::ShortFirst,
            terminal_order: TerminalOrder::BottomFirst,
            descent_mode: DescentMode::None,
            augment_potential: TargetPotential::Lex,
            tableau_limit: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CarrierOrder {
    ShortFirst,
    LongFirst,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TerminalOrder {
    BottomFirst,
    TopFirst,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SourceOrder {
    Lex,
    ReverseLex,
    FewestTargets,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DescentMode {
    None,
    Global,
    Componentwise,
    ActiveGlobal,
    ActiveComponentwise,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TargetPotential {
    Lex,
    Colex,
    SumRight,
    SumLeft,
    SumReverseLeft,
    SumReverseRight,
    RightReverseLeft,
    MaxMin,
    MinMax,
}

#[derive(Debug, Clone)]
struct ShapeOutcome {
    fibers_checked: usize,
    negative_pairs: usize,
    failures: Vec<SlideFailure>,
    limit_exceeded: bool,
}

#[derive(Debug, Clone)]
struct SlideFailure {
    beta: Content,
    message: String,
    left_pos: usize,
    right_pos: usize,
    left_values: Vec<u32>,
    right_values: Vec<u32>,
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
        let flagged_shape = flagged_shape_for_args(shape, &args);
        let outcome = scan_one_shape(&flagged_shape, &args);
        print_single_outcome(&flagged_shape, &outcome, started.elapsed().as_secs_f64());
        std::process::exit(if outcome.failures.is_empty() { 0 } else { 1 });
    }

    let mut shapes_checked = 0usize;
    let mut skipped_by_limit = 0usize;
    let mut fibers_checked = 0usize;
    let mut negative_pairs = 0usize;

    for skew_size in 0..=args.max_skew_size {
        let max_outer_size = skew_size + args.max_outer_extra;
        for shape in skew_shapes_of_size(skew_size, max_outer_size) {
            if shape.row_count() > 2 {
                continue;
            }
            if args.exactly_two_rows && shape.row_count() != 2 {
                continue;
            }
            let flagged_shape = flagged_shape_for_args(shape, &args);
            let outcome = scan_one_shape(&flagged_shape, &args);
            if outcome.limit_exceeded {
                skipped_by_limit += 1;
                continue;
            }
            shapes_checked += 1;
            fibers_checked += outcome.fibers_checked;
            negative_pairs += outcome.negative_pairs;
            if let Some(failure) = outcome.failures.first() {
                println!("FAIL");
                println!("outer={:?}", flagged_shape.shape().outer().parts());
                println!("inner={:?}", flagged_shape.shape().inner().parts());
                println!("flags={:?}", flagged_shape.row_flags());
                println!(
                    "beta={:?}, {}, first_pair=({}, {})",
                    failure.beta, failure.message, failure.left_pos, failure.right_pos
                );
                println!(
                    "left={}, right={}",
                    format_tableau(flagged_shape.shape(), &failure.left_values),
                    format_tableau(flagged_shape.shape(), &failure.right_values)
                );
                println!(
                    "shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, fibers_checked={fibers_checked}, negative_pairs={negative_pairs}, elapsed={:.3}s",
                    started.elapsed().as_secs_f64()
                );
                std::process::exit(1);
            }
        }
        println!(
            "skew_size={skew_size}: shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, fibers_checked={fibers_checked}, negative_pairs={negative_pairs}, elapsed={:.3}s",
            started.elapsed().as_secs_f64()
        );
    }

    println!("PASS");
    println!(
        "shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, fibers_checked={fibers_checked}, negative_pairs={negative_pairs}, elapsed={:.3}s",
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
                "--row-flags" => {
                    args.row_flags = Some(parse_u32_vec(&take_value(&mut iter, &flag)?)?)
                }
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
                "--exactly-two-rows" => args.exactly_two_rows = true,
                "--check-injective" => args.check_injective = true,
                "--check-matching" => args.check_matching = true,
                "--check-greedy" => args.check_greedy = true,
                "--check-descent-output" => args.check_descent_output = true,
                "--check-target-intervals" => args.check_target_intervals = true,
                "--check-target-swap-closed" => args.check_target_swap_closed = true,
                "--check-unordered-target-intervals" => {
                    args.check_unordered_target_intervals = true
                }
                "--check-augment-target-increase" => args.check_augment_target_increase = true,
                "--check-additive-rank" => args.check_additive_rank = true,
                "--check-layered-additive-rank" => args.check_layered_additive_rank = true,
                "--check-carrier-word-additive-rank" => {
                    args.check_carrier_word_additive_rank = true
                }
                "--source-order" => {
                    args.source_order = SourceOrder::parse(&take_value(&mut iter, &flag)?)?
                }
                "--carrier-order" => {
                    args.carrier_order = CarrierOrder::parse(&take_value(&mut iter, &flag)?)?
                }
                "--terminal-order" => {
                    args.terminal_order = TerminalOrder::parse(&take_value(&mut iter, &flag)?)?
                }
                "--descent-mode" => {
                    args.descent_mode = DescentMode::parse(&take_value(&mut iter, &flag)?)?
                }
                "--augment-potential" => {
                    args.augment_potential = TargetPotential::parse(&take_value(&mut iter, &flag)?)?
                }
                "--max-beta-active" => {
                    args.max_beta_active = Some(
                        take_value(&mut iter, &flag)?
                            .parse()
                            .map_err(|err| format!("invalid --max-beta-active: {err}"))?,
                    )
                }
                "--tableau-limit" => {
                    args.tableau_limit = Some(
                        take_value(&mut iter, &flag)?
                            .parse()
                            .map_err(|err| format!("invalid --tableau-limit: {err}"))?,
                    )
                }
                other => return Err(format!("unknown argument `{other}`")),
            }
        }
        if args.lower_label as usize >= args.alphabet {
            return Err("--lower-label must be smaller than --alphabet".to_string());
        }
        Ok(args)
    }
}

impl CarrierOrder {
    fn parse(input: &str) -> Result<Self, String> {
        match input {
            "short-first" => Ok(Self::ShortFirst),
            "long-first" => Ok(Self::LongFirst),
            other => Err(format!(
                "invalid --carrier-order `{other}`; expected `short-first` or `long-first`"
            )),
        }
    }
}

impl TerminalOrder {
    fn parse(input: &str) -> Result<Self, String> {
        match input {
            "bottom-first" => Ok(Self::BottomFirst),
            "top-first" => Ok(Self::TopFirst),
            other => Err(format!(
                "invalid --terminal-order `{other}`; expected `bottom-first` or `top-first`"
            )),
        }
    }
}

impl SourceOrder {
    fn parse(input: &str) -> Result<Self, String> {
        match input {
            "lex" => Ok(Self::Lex),
            "reverse-lex" => Ok(Self::ReverseLex),
            "fewest-targets" => Ok(Self::FewestTargets),
            other => Err(format!(
                "invalid --source-order `{other}`; expected `lex`, `reverse-lex`, or `fewest-targets`"
            )),
        }
    }
}

impl DescentMode {
    fn parse(input: &str) -> Result<Self, String> {
        match input {
            "none" => Ok(Self::None),
            "global" => Ok(Self::Global),
            "componentwise" => Ok(Self::Componentwise),
            "active-global" => Ok(Self::ActiveGlobal),
            "active-componentwise" => Ok(Self::ActiveComponentwise),
            other => Err(format!(
                "invalid --descent-mode `{other}`; expected `none`, `global`, `componentwise`, `active-global`, or `active-componentwise`"
            )),
        }
    }

    fn reading_orders(self, shape: &SkewShape) -> Option<Vec<Vec<usize>>> {
        match self {
            Self::None => None,
            Self::Global | Self::ActiveGlobal => {
                Some(DescentStatistic::Global.reading_orders(shape))
            }
            Self::Componentwise | Self::ActiveComponentwise => {
                Some(DescentStatistic::Componentwise.reading_orders(shape))
            }
        }
    }

    fn descent_data(
        self,
        values: &[u32],
        reading_orders: &[Vec<usize>],
        lower_label: u32,
    ) -> Option<DescentData> {
        match self {
            Self::None => None,
            Self::Global | Self::Componentwise => {
                Some(descent_data_for_values(values, reading_orders))
            }
            Self::ActiveGlobal | Self::ActiveComponentwise => Some(
                active_subword_descent_data_for_values(values, reading_orders, lower_label),
            ),
        }
    }
}

impl TargetPotential {
    fn parse(input: &str) -> Result<Self, String> {
        match input {
            "lex" => Ok(Self::Lex),
            "colex" => Ok(Self::Colex),
            "sum-right" => Ok(Self::SumRight),
            "sum-left" => Ok(Self::SumLeft),
            "sum-reverse-left" => Ok(Self::SumReverseLeft),
            "sum-reverse-right" => Ok(Self::SumReverseRight),
            "right-reverse-left" => Ok(Self::RightReverseLeft),
            "max-min" => Ok(Self::MaxMin),
            "min-max" => Ok(Self::MinMax),
            other => Err(format!(
                "invalid --augment-potential `{other}`; expected `lex`, `colex`, `sum-right`, `sum-left`, `sum-reverse-left`, `sum-reverse-right`, `right-reverse-left`, `max-min`, or `min-max`"
            )),
        }
    }

    fn key(self, target: usize, mixed_count: usize) -> [usize; 3] {
        let left = target / mixed_count;
        let right = target % mixed_count;
        match self {
            Self::Lex => [left, right, 0],
            Self::Colex => [right, left, 0],
            Self::SumRight => [left + right, right, 0],
            Self::SumLeft => [left + right, left, 0],
            Self::SumReverseLeft => [left + right, mixed_count - 1 - left, 0],
            Self::SumReverseRight => [left + right, mixed_count - 1 - right, 0],
            Self::RightReverseLeft => [right, mixed_count - 1 - left, 0],
            Self::MaxMin => [left.max(right), left.min(right), 0],
            Self::MinMax => [left.min(right), left.max(right), 0],
        }
    }
}

fn scan_one_shape(flagged_shape: &RowFlaggedSkewShape, args: &Args) -> ShapeOutcome {
    let tableaux = match enumerate_tableaux(flagged_shape, args.tableau_limit) {
        Ok(tableaux) => tableaux,
        Err(_) => {
            return ShapeOutcome {
                fibers_checked: 0,
                negative_pairs: 0,
                failures: Vec::new(),
                limit_exceeded: true,
            }
        }
    };
    let gt_rows: Vec<_> = tableaux
        .iter()
        .map(|tableau| {
            SkewGtPattern::from_tableau(flagged_shape.shape(), &tableau.values, args.alphabet)
                .rows()
                .to_vec()
        })
        .collect();
    let reading_orders = args.descent_mode.reading_orders(flagged_shape.shape());
    let descents: Vec<_> = tableaux
        .iter()
        .map(|tableau| {
            reading_orders.as_ref().and_then(|reading_orders| {
                args.descent_mode
                    .descent_data(&tableau.values, reading_orders, args.lower_label)
            })
        })
        .collect();
    let mut by_content = BTreeMap::<Content, Vec<usize>>::new();
    let mut by_content_and_gt = BTreeMap::<Content, BTreeSet<GtRows>>::new();
    for (idx, tableau) in tableaux.iter().enumerate() {
        by_content
            .entry(tableau.content.clone())
            .or_default()
            .push(idx);
        by_content_and_gt
            .entry(tableau.content.clone())
            .or_default()
            .insert(gt_rows[idx].clone());
    }

    let lower = args.lower_label as usize - 1;
    let upper = lower + 1;
    let mut seen_beta = BTreeSet::new();
    let mut fibers_checked = 0usize;
    let mut negative_pairs = 0usize;
    let mut failures = Vec::new();

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
        if let Some(max_beta_active) = args.max_beta_active {
            if beta[lower] + beta[upper] > max_beta_active {
                continue;
            }
        }

        let (Some(ii_indices), Some(jj_indices), Some(mixed_indices), Some(_mixed_gts)) = (
            by_content.get(&ii),
            by_content.get(&jj),
            by_content.get(&mixed),
            by_content_and_gt.get(&mixed),
        ) else {
            continue;
        };

        let mixed_gt_to_pos: BTreeMap<_, _> = mixed_indices
            .iter()
            .enumerate()
            .map(|(pos, &idx)| (gt_rows[idx].clone(), pos))
            .collect();
        let mut used_targets = BTreeMap::<(usize, usize), (usize, usize)>::new();

        fibers_checked += 1;
        if args.check_matching
            || args.check_greedy
            || args.check_descent_output
            || args.check_target_intervals
            || args.check_target_swap_closed
            || args.check_unordered_target_intervals
            || args.check_augment_target_increase
            || args.check_additive_rank
            || args.check_layered_additive_rank
            || args.check_carrier_word_additive_rank
        {
            let source_count = ii_indices.len() * jj_indices.len();
            negative_pairs += source_count;
            let failure = if args.check_greedy {
                greedy_carrier_matching_failure(
                    ii_indices,
                    jj_indices,
                    mixed_indices.len(),
                    mixed_indices,
                    &gt_rows,
                    &descents,
                    &mixed_gt_to_pos,
                    args.lower_label as usize,
                    args.carrier_order,
                    args.terminal_order,
                    args.source_order,
                )
            } else if args.check_target_intervals {
                carrier_target_interval_failure(
                    ii_indices,
                    jj_indices,
                    mixed_indices.len(),
                    mixed_indices,
                    &gt_rows,
                    &descents,
                    &mixed_gt_to_pos,
                    args.lower_label as usize,
                    args.carrier_order,
                    args.terminal_order,
                )
            } else if args.check_target_swap_closed {
                carrier_target_swap_failure(
                    ii_indices,
                    jj_indices,
                    mixed_indices.len(),
                    mixed_indices,
                    &gt_rows,
                    &descents,
                    &mixed_gt_to_pos,
                    args.lower_label as usize,
                    args.carrier_order,
                    args.terminal_order,
                )
            } else if args.check_unordered_target_intervals {
                carrier_unordered_target_interval_failure(
                    ii_indices,
                    jj_indices,
                    mixed_indices.len(),
                    mixed_indices,
                    &gt_rows,
                    &descents,
                    &mixed_gt_to_pos,
                    args.lower_label as usize,
                    args.carrier_order,
                    args.terminal_order,
                )
            } else if args.check_augment_target_increase {
                carrier_augment_target_increase_failure(
                    ii_indices,
                    jj_indices,
                    mixed_indices.len(),
                    mixed_indices,
                    &gt_rows,
                    &descents,
                    &mixed_gt_to_pos,
                    args.lower_label as usize,
                    args.carrier_order,
                    args.terminal_order,
                    args.augment_potential,
                )
            } else if args.check_additive_rank {
                carrier_additive_rank_failure(
                    ii_indices,
                    jj_indices,
                    mixed_indices.len(),
                    mixed_indices,
                    &gt_rows,
                    &descents,
                    &mixed_gt_to_pos,
                    args.lower_label as usize,
                    args.carrier_order,
                    args.terminal_order,
                )
            } else if args.check_layered_additive_rank {
                carrier_layered_additive_rank_failure(
                    ii_indices,
                    jj_indices,
                    mixed_indices.len(),
                    mixed_indices,
                    &gt_rows,
                    &descents,
                    &mixed_gt_to_pos,
                    args.lower_label as usize,
                    args.carrier_order,
                    args.terminal_order,
                )
            } else if args.check_carrier_word_additive_rank {
                carrier_word_additive_rank_failure(
                    ii_indices,
                    jj_indices,
                    mixed_indices.len(),
                    mixed_indices,
                    &gt_rows,
                    &descents,
                    &mixed_gt_to_pos,
                    args.lower_label as usize,
                    args.carrier_order,
                    args.terminal_order,
                )
            } else if args.check_descent_output {
                canonical_matching_descent_output_failure(
                    ii_indices,
                    jj_indices,
                    mixed_indices.len(),
                    mixed_indices,
                    &gt_rows,
                    &descents,
                    &mixed_gt_to_pos,
                    args.lower_label as usize,
                    args.carrier_order,
                    args.terminal_order,
                )
            } else {
                carrier_matching_failure(
                    ii_indices,
                    jj_indices,
                    mixed_indices.len(),
                    mixed_indices,
                    &gt_rows,
                    &descents,
                    &mixed_gt_to_pos,
                    args.lower_label as usize,
                    args.carrier_order,
                    args.terminal_order,
                )
            };
            if let Some((left_pos, right_pos, message)) = failure {
                let left_idx = ii_indices[left_pos];
                let right_idx = jj_indices[right_pos];
                failures.push(SlideFailure {
                    beta: beta.clone(),
                    message,
                    left_pos,
                    right_pos,
                    left_values: tableaux[left_idx].values.clone(),
                    right_values: tableaux[right_idx].values.clone(),
                });
                return ShapeOutcome {
                    fibers_checked,
                    negative_pairs,
                    failures,
                    limit_exceeded: false,
                };
            }
            continue;
        }

        for (left_pos, &left_idx) in ii_indices.iter().enumerate() {
            for (right_pos, &right_idx) in jj_indices.iter().enumerate() {
                negative_pairs += 1;
                let source_descent_pair = unordered_descent_pair(
                    descents[left_idx].as_ref(),
                    descents[right_idx].as_ref(),
                );
                let target = first_carrier_target(
                    &gt_rows[left_idx],
                    &gt_rows[right_idx],
                    &mixed_gt_to_pos,
                    mixed_indices,
                    &descents,
                    source_descent_pair.as_ref(),
                    args.lower_label as usize,
                    args.carrier_order,
                    args.terminal_order,
                );
                let message = match target {
                    Some(target) if args.check_injective => {
                        if let Some(previous) = used_targets.insert(target, (left_pos, right_pos)) {
                            Some(format!(
                                "carrier collision at target {target:?}; previous source {previous:?}"
                            ))
                        } else {
                            None
                        }
                    }
                    Some(_) => None,
                    None => Some("no carrier target".to_string()),
                };
                let Some(message) = message else { continue };
                failures.push(SlideFailure {
                    beta: beta.clone(),
                    message,
                    left_pos,
                    right_pos,
                    left_values: tableaux[left_idx].values.clone(),
                    right_values: tableaux[right_idx].values.clone(),
                });
                return ShapeOutcome {
                    fibers_checked,
                    negative_pairs,
                    failures,
                    limit_exceeded: false,
                };
            }
        }
    }

    ShapeOutcome {
        fibers_checked,
        negative_pairs,
        failures,
        limit_exceeded: false,
    }
}

fn first_carrier_target(
    left: &[Vec<u32>],
    right: &[Vec<u32>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    mixed_indices: &[usize],
    descents: &[Option<DescentData>],
    source_descent_pair: Option<&DescentPair>,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
) -> Option<(usize, usize)> {
    if left.first().map_or(0, Vec::len) != 2
        || right.first().map_or(0, Vec::len) != 2
        || lower_label >= left.len()
        || lower_label >= right.len()
    {
        return None;
    }

    for carrier in carrier_differences(left.len(), lower_label, carrier_order, terminal_order) {
        let Some(left_image) = apply_carrier(left, &carrier, -1) else {
            continue;
        };
        let Some(&left_pos) = mixed_gt_to_pos.get(&left_image) else {
            continue;
        };
        let Some(right_image) = apply_carrier(right, &carrier, 1) else {
            continue;
        };
        let Some(&right_pos) = mixed_gt_to_pos.get(&right_image) else {
            continue;
        };
        if !descent_pair_allowed(
            left_pos,
            right_pos,
            mixed_indices,
            descents,
            source_descent_pair,
        ) {
            continue;
        }
        return Some((left_pos, right_pos));
    }
    None
}

fn carrier_matching_failure(
    ii_indices: &[usize],
    jj_indices: &[usize],
    mixed_count: usize,
    mixed_indices: &[usize],
    gt_rows: &[GtRows],
    descents: &[Option<DescentData>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
) -> Option<(usize, usize, String)> {
    let source_pairs: Vec<_> = (0..ii_indices.len())
        .flat_map(|left_pos| (0..jj_indices.len()).map(move |right_pos| (left_pos, right_pos)))
        .collect();
    let mut edges_by_source = Vec::with_capacity(source_pairs.len());
    for &(left_pos, right_pos) in &source_pairs {
        let left_idx = ii_indices[left_pos];
        let right_idx = jj_indices[right_pos];
        let edges = carrier_targets(
            &gt_rows[left_idx],
            &gt_rows[right_idx],
            mixed_gt_to_pos,
            mixed_indices,
            descents,
            unordered_descent_pair(descents[left_idx].as_ref(), descents[right_idx].as_ref())
                .as_ref(),
            mixed_count,
            lower_label,
            carrier_order,
            terminal_order,
        );
        if edges.is_empty() {
            return Some((left_pos, right_pos, "no carrier target".to_string()));
        }
        edges_by_source.push(edges);
    }

    let mut target_match = vec![None; mixed_count * mixed_count];
    for source in 0..source_pairs.len() {
        let mut seen = vec![false; target_match.len()];
        if !augment_carrier_matching(source, &edges_by_source, &mut seen, &mut target_match) {
            let (left_pos, right_pos) = source_pairs[source];
            return Some((
                left_pos,
                right_pos,
                format!("carrier graph has no augmenting match at source {source}"),
            ));
        }
    }
    None
}

fn carrier_target_interval_failure(
    ii_indices: &[usize],
    jj_indices: &[usize],
    mixed_count: usize,
    mixed_indices: &[usize],
    gt_rows: &[GtRows],
    descents: &[Option<DescentData>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
) -> Option<(usize, usize, String)> {
    for left_pos in 0..ii_indices.len() {
        for right_pos in 0..jj_indices.len() {
            let left_idx = ii_indices[left_pos];
            let right_idx = jj_indices[right_pos];
            let edges = carrier_targets(
                &gt_rows[left_idx],
                &gt_rows[right_idx],
                mixed_gt_to_pos,
                mixed_indices,
                descents,
                unordered_descent_pair(descents[left_idx].as_ref(), descents[right_idx].as_ref())
                    .as_ref(),
                mixed_count,
                lower_label,
                carrier_order,
                terminal_order,
            );
            if edges.is_empty() {
                return Some((left_pos, right_pos, "no carrier target".to_string()));
            }
            if !is_contiguous_interval(&edges) {
                let first = edges[0];
                let last = edges[edges.len() - 1];
                return Some((
                    left_pos,
                    right_pos,
                    format!(
                        "carrier targets are not a lex interval: len={}, span={}..={}",
                        edges.len(),
                        target_pair_string(first, mixed_count),
                        target_pair_string(last, mixed_count)
                    ),
                ));
            }
        }
    }
    None
}

fn carrier_target_swap_failure(
    ii_indices: &[usize],
    jj_indices: &[usize],
    mixed_count: usize,
    mixed_indices: &[usize],
    gt_rows: &[GtRows],
    descents: &[Option<DescentData>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
) -> Option<(usize, usize, String)> {
    for left_pos in 0..ii_indices.len() {
        for right_pos in 0..jj_indices.len() {
            let left_idx = ii_indices[left_pos];
            let right_idx = jj_indices[right_pos];
            let edges = carrier_targets(
                &gt_rows[left_idx],
                &gt_rows[right_idx],
                mixed_gt_to_pos,
                mixed_indices,
                descents,
                unordered_descent_pair(descents[left_idx].as_ref(), descents[right_idx].as_ref())
                    .as_ref(),
                mixed_count,
                lower_label,
                carrier_order,
                terminal_order,
            );
            if edges.is_empty() {
                return Some((left_pos, right_pos, "no carrier target".to_string()));
            }
            let edge_set: BTreeSet<_> = edges.iter().copied().collect();
            for target in edges {
                let swapped = swap_target(target, mixed_count);
                if !edge_set.contains(&swapped) {
                    return Some((
                        left_pos,
                        right_pos,
                        format!(
                            "carrier targets are not swap-closed: has {}, missing {}",
                            target_pair_string(target, mixed_count),
                            target_pair_string(swapped, mixed_count)
                        ),
                    ));
                }
            }
        }
    }
    None
}

fn carrier_unordered_target_interval_failure(
    ii_indices: &[usize],
    jj_indices: &[usize],
    mixed_count: usize,
    mixed_indices: &[usize],
    gt_rows: &[GtRows],
    descents: &[Option<DescentData>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
) -> Option<(usize, usize, String)> {
    for left_pos in 0..ii_indices.len() {
        for right_pos in 0..jj_indices.len() {
            let left_idx = ii_indices[left_pos];
            let right_idx = jj_indices[right_pos];
            let edges = carrier_targets(
                &gt_rows[left_idx],
                &gt_rows[right_idx],
                mixed_gt_to_pos,
                mixed_indices,
                descents,
                unordered_descent_pair(descents[left_idx].as_ref(), descents[right_idx].as_ref())
                    .as_ref(),
                mixed_count,
                lower_label,
                carrier_order,
                terminal_order,
            );
            if edges.is_empty() {
                return Some((left_pos, right_pos, "no carrier target".to_string()));
            }
            let unordered_edges: Vec<_> = edges
                .iter()
                .map(|&target| unordered_target_index(target, mixed_count))
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            if !is_contiguous_interval(&unordered_edges) {
                let first = unordered_edges[0];
                let last = unordered_edges[unordered_edges.len() - 1];
                return Some((
                    left_pos,
                    right_pos,
                    format!(
                        "unordered carrier targets are not an interval: len={}, span={}..={}",
                        unordered_edges.len(),
                        first,
                        last
                    ),
                ));
            }
        }
    }
    None
}

fn carrier_augment_target_increase_failure(
    ii_indices: &[usize],
    jj_indices: &[usize],
    mixed_count: usize,
    mixed_indices: &[usize],
    gt_rows: &[GtRows],
    descents: &[Option<DescentData>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
    augment_potential: TargetPotential,
) -> Option<(usize, usize, String)> {
    let source_pairs: Vec<_> = (0..ii_indices.len())
        .flat_map(|left_pos| (0..jj_indices.len()).map(move |right_pos| (left_pos, right_pos)))
        .collect();
    let mut edges_by_source = Vec::with_capacity(source_pairs.len());
    for &(left_pos, right_pos) in &source_pairs {
        let left_idx = ii_indices[left_pos];
        let right_idx = jj_indices[right_pos];
        let edges = carrier_targets(
            &gt_rows[left_idx],
            &gt_rows[right_idx],
            mixed_gt_to_pos,
            mixed_indices,
            descents,
            unordered_descent_pair(descents[left_idx].as_ref(), descents[right_idx].as_ref())
                .as_ref(),
            mixed_count,
            lower_label,
            carrier_order,
            terminal_order,
        );
        if edges.is_empty() {
            return Some((left_pos, right_pos, "no carrier target".to_string()));
        }
        edges_by_source.push(edges);
    }

    let mut target_match = vec![None; mixed_count * mixed_count];
    let mut source_match = vec![None; source_pairs.len()];
    for source in 0..source_pairs.len() {
        let mut seen = vec![false; target_match.len()];
        let mut violation = None;
        if !augment_carrier_matching_with_increase_check(
            source,
            &edges_by_source,
            &mut seen,
            &mut target_match,
            &mut source_match,
            mixed_count,
            augment_potential,
            &mut violation,
        ) {
            let (left_pos, right_pos) = source_pairs[source];
            let message = violation.unwrap_or_else(|| {
                format!("carrier graph has no augmenting match at source {source}")
            });
            return Some((left_pos, right_pos, message));
        }
    }
    None
}

fn carrier_additive_rank_failure(
    ii_indices: &[usize],
    jj_indices: &[usize],
    mixed_count: usize,
    mixed_indices: &[usize],
    gt_rows: &[GtRows],
    descents: &[Option<DescentData>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
) -> Option<(usize, usize, String)> {
    if ii_indices.is_empty() || jj_indices.is_empty() {
        return None;
    }

    let mut diagonal_sums = vec![vec![0i32; jj_indices.len()]; ii_indices.len()];
    for left_pos in 0..ii_indices.len() {
        for right_pos in 0..jj_indices.len() {
            let left_idx = ii_indices[left_pos];
            let right_idx = jj_indices[right_pos];
            let edges = carrier_targets(
                &gt_rows[left_idx],
                &gt_rows[right_idx],
                mixed_gt_to_pos,
                mixed_indices,
                descents,
                unordered_descent_pair(descents[left_idx].as_ref(), descents[right_idx].as_ref())
                    .as_ref(),
                mixed_count,
                lower_label,
                carrier_order,
                terminal_order,
            );
            if edges.is_empty() {
                return Some((left_pos, right_pos, "no carrier target".to_string()));
            }
            let target_sums: BTreeSet<_> = edges
                .iter()
                .map(|&target| target / mixed_count + target % mixed_count)
                .collect();
            if target_sums.len() != 1 {
                return Some((
                    left_pos,
                    right_pos,
                    format!("carrier targets lie on multiple target diagonals: {target_sums:?}"),
                ));
            }
            diagonal_sums[left_pos][right_pos] = *target_sums.first().unwrap() as i32;
        }
    }

    let base = diagonal_sums[0][0];
    for left_pos in 0..ii_indices.len() {
        for right_pos in 0..jj_indices.len() {
            let expected = diagonal_sums[left_pos][0] + diagonal_sums[0][right_pos] - base;
            if diagonal_sums[left_pos][right_pos] != expected {
                return Some((
                    left_pos,
                    right_pos,
                    format!(
                        "target diagonal is not additive: got {}, expected {}",
                        diagonal_sums[left_pos][right_pos], expected
                    ),
                ));
            }
        }
    }

    None
}

fn carrier_layered_additive_rank_failure(
    ii_indices: &[usize],
    jj_indices: &[usize],
    mixed_count: usize,
    mixed_indices: &[usize],
    gt_rows: &[GtRows],
    descents: &[Option<DescentData>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
) -> Option<(usize, usize, String)> {
    for start_level in 1..=lower_label {
        let mut entries = Vec::new();
        for left_pos in 0..ii_indices.len() {
            for right_pos in 0..jj_indices.len() {
                let left_idx = ii_indices[left_pos];
                let right_idx = jj_indices[right_pos];
                let targets = carrier_targets_with_start(
                    &gt_rows[left_idx],
                    &gt_rows[right_idx],
                    mixed_gt_to_pos,
                    mixed_indices,
                    descents,
                    unordered_descent_pair(
                        descents[left_idx].as_ref(),
                        descents[right_idx].as_ref(),
                    )
                    .as_ref(),
                    mixed_count,
                    lower_label,
                    carrier_order,
                    terminal_order,
                );
                if targets.is_empty() {
                    return Some((left_pos, right_pos, "no carrier target".to_string()));
                }
                let layer_diagonals: BTreeSet<_> = targets
                    .iter()
                    .filter_map(|&(start, target)| {
                        (start == start_level)
                            .then_some((target / mixed_count + target % mixed_count) as i32)
                    })
                    .collect();
                match layer_diagonals.len() {
                    0 => {}
                    1 => entries.push((
                        left_pos,
                        right_pos,
                        *layer_diagonals.first().unwrap(),
                    )),
                    _ => {
                        return Some((
                            left_pos,
                            right_pos,
                            format!(
                                "start level {start_level} has multiple target diagonals: {layer_diagonals:?}"
                            ),
                        ))
                    }
                }
            }
        }

        if let Some((left_pos, right_pos, message)) =
            layered_additive_entries_failure(&entries, ii_indices.len(), jj_indices.len())
        {
            return Some((
                left_pos,
                right_pos,
                format!("start level {start_level}: {message}"),
            ));
        }
    }

    None
}

fn carrier_word_additive_rank_failure(
    ii_indices: &[usize],
    jj_indices: &[usize],
    mixed_count: usize,
    mixed_indices: &[usize],
    gt_rows: &[GtRows],
    descents: &[Option<DescentData>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
) -> Option<(usize, usize, String)> {
    let mut entries_by_word = BTreeMap::<usize, Vec<(usize, usize, i32)>>::new();
    for left_pos in 0..ii_indices.len() {
        for right_pos in 0..jj_indices.len() {
            let left_idx = ii_indices[left_pos];
            let right_idx = jj_indices[right_pos];
            let targets = carrier_targets_with_word(
                &gt_rows[left_idx],
                &gt_rows[right_idx],
                mixed_gt_to_pos,
                mixed_indices,
                descents,
                unordered_descent_pair(descents[left_idx].as_ref(), descents[right_idx].as_ref())
                    .as_ref(),
                mixed_count,
                lower_label,
                carrier_order,
                terminal_order,
            );
            if targets.is_empty() {
                return Some((left_pos, right_pos, "no carrier target".to_string()));
            }
            for (word, target) in targets {
                entries_by_word.entry(word).or_default().push((
                    left_pos,
                    right_pos,
                    (target / mixed_count + target % mixed_count) as i32,
                ));
            }
        }
    }

    for (word, entries) in entries_by_word {
        if let Some((left_pos, right_pos, message)) =
            layered_additive_entries_failure(&entries, ii_indices.len(), jj_indices.len())
        {
            return Some((
                left_pos,
                right_pos,
                format!("carrier word {word}: {message}"),
            ));
        }
    }

    None
}

fn layered_additive_entries_failure(
    entries: &[(usize, usize, i32)],
    left_count: usize,
    right_count: usize,
) -> Option<(usize, usize, String)> {
    let mut left_adj = vec![Vec::<(usize, i32)>::new(); left_count];
    let mut right_adj = vec![Vec::<(usize, i32)>::new(); right_count];
    for &(left_pos, right_pos, diagonal) in entries {
        left_adj[left_pos].push((right_pos, diagonal));
        right_adj[right_pos].push((left_pos, diagonal));
    }

    let mut left_rank = vec![None::<i32>; left_count];
    let mut right_rank = vec![None::<i32>; right_count];
    for start_left in 0..left_count {
        if left_rank[start_left].is_some() || left_adj[start_left].is_empty() {
            continue;
        }
        left_rank[start_left] = Some(0);
        let mut stack = vec![RankNode::Left(start_left)];
        while let Some(node) = stack.pop() {
            match node {
                RankNode::Left(left_pos) => {
                    let left_value = left_rank[left_pos].unwrap();
                    for &(right_pos, diagonal) in &left_adj[left_pos] {
                        let expected = diagonal - left_value;
                        match right_rank[right_pos] {
                            Some(value) if value != expected => {
                                return Some((
                                    left_pos,
                                    right_pos,
                                    format!(
                                        "non-additive layer: right rank conflict, got {value}, expected {expected}"
                                    ),
                                ))
                            }
                            Some(_) => {}
                            None => {
                                right_rank[right_pos] = Some(expected);
                                stack.push(RankNode::Right(right_pos));
                            }
                        }
                    }
                }
                RankNode::Right(right_pos) => {
                    let right_value = right_rank[right_pos].unwrap();
                    for &(left_pos, diagonal) in &right_adj[right_pos] {
                        let expected = diagonal - right_value;
                        match left_rank[left_pos] {
                            Some(value) if value != expected => {
                                return Some((
                                    left_pos,
                                    right_pos,
                                    format!(
                                        "non-additive layer: left rank conflict, got {value}, expected {expected}"
                                    ),
                                ))
                            }
                            Some(_) => {}
                            None => {
                                left_rank[left_pos] = Some(expected);
                                stack.push(RankNode::Left(left_pos));
                            }
                        }
                    }
                }
            }
        }
    }

    for start_right in 0..right_count {
        if right_rank[start_right].is_some() || right_adj[start_right].is_empty() {
            continue;
        }
        right_rank[start_right] = Some(0);
        let mut stack = vec![RankNode::Right(start_right)];
        while let Some(node) = stack.pop() {
            match node {
                RankNode::Left(left_pos) => {
                    let left_value = left_rank[left_pos].unwrap();
                    for &(right_pos, diagonal) in &left_adj[left_pos] {
                        let expected = diagonal - left_value;
                        match right_rank[right_pos] {
                            Some(value) if value != expected => {
                                return Some((
                                    left_pos,
                                    right_pos,
                                    format!(
                                        "non-additive right component: got {value}, expected {expected}"
                                    ),
                                ))
                            }
                            Some(_) => {}
                            None => {
                                right_rank[right_pos] = Some(expected);
                                stack.push(RankNode::Right(right_pos));
                            }
                        }
                    }
                }
                RankNode::Right(right_pos) => {
                    let right_value = right_rank[right_pos].unwrap();
                    for &(left_pos, diagonal) in &right_adj[right_pos] {
                        let expected = diagonal - right_value;
                        match left_rank[left_pos] {
                            Some(value) if value != expected => {
                                return Some((
                                    left_pos,
                                    right_pos,
                                    format!(
                                        "non-additive right component: got {value}, expected {expected}"
                                    ),
                                ))
                            }
                            Some(_) => {}
                            None => {
                                left_rank[left_pos] = Some(expected);
                                stack.push(RankNode::Left(left_pos));
                            }
                        }
                    }
                }
            }
        }
    }

    None
}

#[derive(Debug, Clone, Copy)]
enum RankNode {
    Left(usize),
    Right(usize),
}

fn canonical_matching_descent_output_failure(
    ii_indices: &[usize],
    jj_indices: &[usize],
    mixed_count: usize,
    mixed_indices: &[usize],
    gt_rows: &[GtRows],
    descents: &[Option<DescentData>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
) -> Option<(usize, usize, String)> {
    let source_pairs: Vec<_> = (0..ii_indices.len())
        .flat_map(|left_pos| (0..jj_indices.len()).map(move |right_pos| (left_pos, right_pos)))
        .collect();
    let mut edges_by_source = Vec::with_capacity(source_pairs.len());
    for &(left_pos, right_pos) in &source_pairs {
        let left_idx = ii_indices[left_pos];
        let right_idx = jj_indices[right_pos];
        let edges = carrier_targets(
            &gt_rows[left_idx],
            &gt_rows[right_idx],
            mixed_gt_to_pos,
            mixed_indices,
            descents,
            None,
            mixed_count,
            lower_label,
            carrier_order,
            terminal_order,
        );
        if edges.is_empty() {
            return Some((left_pos, right_pos, "no carrier target".to_string()));
        }
        edges_by_source.push(edges);
    }

    let mut target_match = vec![None; mixed_count * mixed_count];
    for source in 0..source_pairs.len() {
        let mut seen = vec![false; target_match.len()];
        if !augment_carrier_matching(source, &edges_by_source, &mut seen, &mut target_match) {
            let (left_pos, right_pos) = source_pairs[source];
            return Some((
                left_pos,
                right_pos,
                format!("carrier graph has no augmenting match at source {source}"),
            ));
        }
    }

    let mut source_to_target = vec![None; source_pairs.len()];
    for (target, &source) in target_match.iter().enumerate() {
        if let Some(source) = source {
            source_to_target[source] = Some(target);
        }
    }

    for (source, &(left_pos, right_pos)) in source_pairs.iter().enumerate() {
        let Some(target) = source_to_target[source] else {
            return Some((
                left_pos,
                right_pos,
                format!("canonical matching left source {source} unmatched"),
            ));
        };
        let target_left_pos = target / mixed_count;
        let target_right_pos = target % mixed_count;
        let source_left_idx = ii_indices[left_pos];
        let source_right_idx = jj_indices[right_pos];
        let target_left_idx = mixed_indices[target_left_pos];
        let target_right_idx = mixed_indices[target_right_pos];
        let source_descent_pair = unordered_descent_pair(
            descents[source_left_idx].as_ref(),
            descents[source_right_idx].as_ref(),
        );
        let target_descent_pair = unordered_descent_pair(
            descents[target_left_idx].as_ref(),
            descents[target_right_idx].as_ref(),
        );
        if source_descent_pair != target_descent_pair {
            return Some((
                left_pos,
                right_pos,
                format!(
                    "canonical unrestricted matching changes descents at source {source}: target=({target_left_pos},{target_right_pos})"
                ),
            ));
        }
    }

    None
}

fn is_contiguous_interval(values: &[usize]) -> bool {
    values.windows(2).all(|window| window[1] == window[0] + 1)
}

fn target_pair_string(target: usize, mixed_count: usize) -> String {
    format!("({}, {})", target / mixed_count, target % mixed_count)
}

fn swap_target(target: usize, mixed_count: usize) -> usize {
    let left = target / mixed_count;
    let right = target % mixed_count;
    right * mixed_count + left
}

fn unordered_target_index(target: usize, mixed_count: usize) -> usize {
    let mut left = target / mixed_count;
    let mut right = target % mixed_count;
    if left > right {
        std::mem::swap(&mut left, &mut right);
    }
    let skipped: usize = (0..left).map(|row| mixed_count - row).sum();
    skipped + (right - left)
}

fn greedy_carrier_matching_failure(
    ii_indices: &[usize],
    jj_indices: &[usize],
    mixed_count: usize,
    mixed_indices: &[usize],
    gt_rows: &[GtRows],
    descents: &[Option<DescentData>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
    source_order: SourceOrder,
) -> Option<(usize, usize, String)> {
    let source_pairs: Vec<_> = (0..ii_indices.len())
        .flat_map(|left_pos| (0..jj_indices.len()).map(move |right_pos| (left_pos, right_pos)))
        .collect();
    let mut edges_by_source = Vec::with_capacity(source_pairs.len());
    for &(left_pos, right_pos) in &source_pairs {
        let left_idx = ii_indices[left_pos];
        let right_idx = jj_indices[right_pos];
        let edges = carrier_targets_in_order(
            &gt_rows[left_idx],
            &gt_rows[right_idx],
            mixed_gt_to_pos,
            mixed_indices,
            descents,
            unordered_descent_pair(descents[left_idx].as_ref(), descents[right_idx].as_ref())
                .as_ref(),
            mixed_count,
            lower_label,
            carrier_order,
            terminal_order,
        );
        if edges.is_empty() {
            return Some((left_pos, right_pos, "no carrier target".to_string()));
        }
        edges_by_source.push(edges);
    }

    let mut source_ordering: Vec<_> = (0..source_pairs.len()).collect();
    match source_order {
        SourceOrder::Lex => {}
        SourceOrder::ReverseLex => source_ordering.reverse(),
        SourceOrder::FewestTargets => source_ordering
            .sort_by_key(|&source| (edges_by_source[source].len(), source_pairs[source])),
    }

    let mut used_targets = vec![false; mixed_count * mixed_count];
    for source in source_ordering {
        let Some(&target) = edges_by_source[source]
            .iter()
            .find(|&&target| !used_targets[target])
        else {
            let (left_pos, right_pos) = source_pairs[source];
            return Some((
                left_pos,
                right_pos,
                format!("greedy carrier matching stuck at source {source}"),
            ));
        };
        used_targets[target] = true;
    }
    None
}

fn augment_carrier_matching(
    source: usize,
    edges_by_source: &[Vec<usize>],
    seen: &mut [bool],
    target_match: &mut [Option<usize>],
) -> bool {
    for &target in &edges_by_source[source] {
        if seen[target] {
            continue;
        }
        seen[target] = true;
        let can_use = match target_match[target] {
            None => true,
            Some(previous) => {
                augment_carrier_matching(previous, edges_by_source, seen, target_match)
            }
        };
        if can_use {
            target_match[target] = Some(source);
            return true;
        }
    }
    false
}

fn augment_carrier_matching_with_increase_check(
    source: usize,
    edges_by_source: &[Vec<usize>],
    seen: &mut [bool],
    target_match: &mut [Option<usize>],
    source_match: &mut [Option<usize>],
    mixed_count: usize,
    augment_potential: TargetPotential,
    violation: &mut Option<String>,
) -> bool {
    for &target in &edges_by_source[source] {
        if seen[target] {
            continue;
        }
        seen[target] = true;
        let can_use = match target_match[target] {
            None => true,
            Some(previous) => augment_carrier_matching_with_increase_check(
                previous,
                edges_by_source,
                seen,
                target_match,
                source_match,
                mixed_count,
                augment_potential,
                violation,
            ),
        };
        if !can_use {
            if violation.is_some() {
                return false;
            }
            continue;
        }

        if let Some(old_target) = source_match[source] {
            let old_key = augment_potential.key(old_target, mixed_count);
            let new_key = augment_potential.key(target, mixed_count);
            if new_key <= old_key {
                *violation = Some(format!(
                    "augmenting move does not increase {augment_potential:?} potential for source {source}: {} -> {}",
                    target_pair_string(old_target, mixed_count),
                    target_pair_string(target, mixed_count)
                ));
                return false;
            }
        }

        target_match[target] = Some(source);
        source_match[source] = Some(target);
        return true;
    }
    false
}

fn carrier_targets(
    left: &[Vec<u32>],
    right: &[Vec<u32>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    mixed_indices: &[usize],
    descents: &[Option<DescentData>],
    source_descent_pair: Option<&DescentPair>,
    mixed_count: usize,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
) -> Vec<usize> {
    let mut out = carrier_targets_in_order(
        left,
        right,
        mixed_gt_to_pos,
        mixed_indices,
        descents,
        source_descent_pair,
        mixed_count,
        lower_label,
        carrier_order,
        terminal_order,
    );
    out.sort_unstable();
    out
}

fn carrier_targets_with_start(
    left: &[Vec<u32>],
    right: &[Vec<u32>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    mixed_indices: &[usize],
    descents: &[Option<DescentData>],
    source_descent_pair: Option<&DescentPair>,
    mixed_count: usize,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut seen = BTreeSet::new();
    for carrier in carrier_differences(left.len(), lower_label, carrier_order, terminal_order) {
        let Some(start_level) = carrier_start_level(&carrier) else {
            continue;
        };
        let Some(left_image) = apply_carrier(left, &carrier, -1) else {
            continue;
        };
        let Some(&left_pos) = mixed_gt_to_pos.get(&left_image) else {
            continue;
        };
        let Some(right_image) = apply_carrier(right, &carrier, 1) else {
            continue;
        };
        let Some(&right_pos) = mixed_gt_to_pos.get(&right_image) else {
            continue;
        };
        if !descent_pair_allowed(
            left_pos,
            right_pos,
            mixed_indices,
            descents,
            source_descent_pair,
        ) {
            continue;
        }
        let target = left_pos * mixed_count + right_pos;
        if seen.insert((start_level, target)) {
            out.push((start_level, target));
        }
    }
    out
}

fn carrier_targets_with_word(
    left: &[Vec<u32>],
    right: &[Vec<u32>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    mixed_indices: &[usize],
    descents: &[Option<DescentData>],
    source_descent_pair: Option<&DescentPair>,
    mixed_count: usize,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut seen = BTreeSet::new();
    for (word, carrier) in
        carrier_differences(left.len(), lower_label, carrier_order, terminal_order)
            .into_iter()
            .enumerate()
    {
        let Some(left_image) = apply_carrier(left, &carrier, -1) else {
            continue;
        };
        let Some(&left_pos) = mixed_gt_to_pos.get(&left_image) else {
            continue;
        };
        let Some(right_image) = apply_carrier(right, &carrier, 1) else {
            continue;
        };
        let Some(&right_pos) = mixed_gt_to_pos.get(&right_image) else {
            continue;
        };
        if !descent_pair_allowed(
            left_pos,
            right_pos,
            mixed_indices,
            descents,
            source_descent_pair,
        ) {
            continue;
        }
        let target = left_pos * mixed_count + right_pos;
        if seen.insert((word, target)) {
            out.push((word, target));
        }
    }
    out
}

fn carrier_start_level(carrier: &[[i32; 2]]) -> Option<usize> {
    carrier.iter().position(|&delta| delta != [0, 0])
}

fn carrier_targets_in_order(
    left: &[Vec<u32>],
    right: &[Vec<u32>],
    mixed_gt_to_pos: &BTreeMap<GtRows, usize>,
    mixed_indices: &[usize],
    descents: &[Option<DescentData>],
    source_descent_pair: Option<&DescentPair>,
    mixed_count: usize,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
) -> Vec<usize> {
    let mut out = Vec::new();
    let mut seen = BTreeSet::new();
    for carrier in carrier_differences(left.len(), lower_label, carrier_order, terminal_order) {
        let Some(left_image) = apply_carrier(left, &carrier, -1) else {
            continue;
        };
        let Some(&left_pos) = mixed_gt_to_pos.get(&left_image) else {
            continue;
        };
        let Some(right_image) = apply_carrier(right, &carrier, 1) else {
            continue;
        };
        let Some(&right_pos) = mixed_gt_to_pos.get(&right_image) else {
            continue;
        };
        if !descent_pair_allowed(
            left_pos,
            right_pos,
            mixed_indices,
            descents,
            source_descent_pair,
        ) {
            continue;
        }
        let target = left_pos * mixed_count + right_pos;
        if seen.insert(target) {
            out.push(target);
        }
    }
    out
}

fn descent_pair_allowed(
    left_pos: usize,
    right_pos: usize,
    mixed_indices: &[usize],
    descents: &[Option<DescentData>],
    source_descent_pair: Option<&DescentPair>,
) -> bool {
    let Some(source_descent_pair) = source_descent_pair else {
        return true;
    };
    let left_idx = mixed_indices[left_pos];
    let right_idx = mixed_indices[right_pos];
    unordered_descent_pair(descents[left_idx].as_ref(), descents[right_idx].as_ref()).as_ref()
        == Some(source_descent_pair)
}

fn unordered_descent_pair(
    left: Option<&DescentData>,
    right: Option<&DescentData>,
) -> Option<DescentPair> {
    let (left, right) = (left?.clone(), right?.clone());
    if left <= right {
        Some((left, right))
    } else {
        Some((right, left))
    }
}

fn carrier_differences(
    row_count: usize,
    lower_label: usize,
    carrier_order: CarrierOrder,
    terminal_order: TerminalOrder,
) -> Vec<Vec<[i32; 2]>> {
    let mut carriers = Vec::new();
    let start_levels: Vec<_> = match carrier_order {
        CarrierOrder::ShortFirst => (1..=lower_label).rev().collect(),
        CarrierOrder::LongFirst => (1..=lower_label).collect(),
    };
    for start_level in start_levels {
        let mut carrier = vec![[0, 0]; row_count];
        extend_carrier(
            &mut carriers,
            &mut carrier,
            start_level,
            lower_label,
            terminal_order,
        );
    }
    carriers
}

fn extend_carrier(
    carriers: &mut Vec<Vec<[i32; 2]>>,
    carrier: &mut Vec<[i32; 2]>,
    level: usize,
    lower_label: usize,
    terminal_order: TerminalOrder,
) {
    if level == lower_label {
        for terminal in terminal_deltas(terminal_order) {
            carrier[level] = terminal;
            carriers.push(carrier.clone());
        }
        carrier[level] = [0, 0];
        return;
    }

    carrier[level] = [1, -1];
    extend_carrier(carriers, carrier, level + 1, lower_label, terminal_order);
    carrier[level] = [-1, 1];
    extend_carrier(carriers, carrier, level + 1, lower_label, terminal_order);
    carrier[level] = [0, 0];
}

fn terminal_deltas(terminal_order: TerminalOrder) -> [[i32; 2]; 2] {
    match terminal_order {
        TerminalOrder::BottomFirst => [[0, 1], [1, 0]],
        TerminalOrder::TopFirst => [[1, 0], [0, 1]],
    }
}

fn apply_carrier(rows: &[Vec<u32>], carrier: &[[i32; 2]], sign: i32) -> Option<GtRows> {
    let mut out = rows.to_vec();
    for (row, delta) in out.iter_mut().zip(carrier) {
        for col in 0..2 {
            let value = row[col] as i32 + sign * delta[col];
            if value < 0 {
                return None;
            }
            row[col] = value as u32;
        }
    }
    is_gt_array(&out).then_some(out)
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

fn flagged_shape_for_args(shape: SkewShape, args: &Args) -> RowFlaggedSkewShape {
    match &args.row_flags {
        Some(row_flags) => RowFlaggedSkewShape::new(shape, row_flags.clone(), args.alphabet),
        None => RowFlaggedSkewShape::ordinary(shape, args.alphabet),
    }
}

fn print_single_outcome(
    flagged_shape: &RowFlaggedSkewShape,
    outcome: &ShapeOutcome,
    elapsed_seconds: f64,
) {
    if outcome.limit_exceeded {
        println!("SKIPPED: tableau limit exceeded");
        return;
    }
    println!("outer={:?}", flagged_shape.shape().outer().parts());
    println!("inner={:?}", flagged_shape.shape().inner().parts());
    println!("flags={:?}", flagged_shape.row_flags());
    println!("alphabet={}", flagged_shape.alphabet());
    println!("fibers_checked={}", outcome.fibers_checked);
    println!("negative_pairs={}", outcome.negative_pairs);
    println!("elapsed={elapsed_seconds:.3}s");
    if let Some(failure) = outcome.failures.first() {
        println!("FAIL");
        println!(
            "beta={:?}, {}, first_pair=({}, {})",
            failure.beta, failure.message, failure.left_pos, failure.right_pos
        );
        println!(
            "left={}, right={}",
            format_tableau(flagged_shape.shape(), &failure.left_values),
            format_tableau(flagged_shape.shape(), &failure.right_values)
        );
    } else {
        println!("PASS");
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

fn format_tableau(shape: &SkewShape, values: &[u32]) -> String {
    let rows = shape
        .cells()
        .iter()
        .map(|cell| cell.row)
        .max()
        .map_or(0, |row| row + 1);
    let cols = shape
        .cells()
        .iter()
        .map(|cell| cell.col)
        .max()
        .map_or(0, |col| col + 1);

    let mut out = Vec::new();
    for row in 0..rows {
        let mut row_values = Vec::new();
        for col in 0..cols {
            let value = shape
                .cell_index(flagged_lorentzian::Cell { row, col })
                .map(|idx| values[idx].to_string())
                .unwrap_or_else(|| ".".to_string());
            row_values.push(value);
        }
        out.push(row_values.join(" "));
    }
    format!("[{}]", out.join(" / "))
}

fn help_text() -> &'static str {
    "Check carrier-path models in two-row shapes.

USAGE:
  two_row_slide_scan [OPTIONS]

OPTIONS:
  --lambda PARTS          Outer partition, e.g. 4,2. If omitted, scan a family.
  --mu PARTS              Inner partition, e.g. 2. Default: empty.
  --row-flags PARTS       Row upper flags. Omit for ordinary skew tableaux.
  --alphabet N            Alphabet size. Default: 5.
  --lower-label I         Active labels are I,I+1. Default: 4.
  --max-skew-size N       Maximum skew size for family scans. Default: 5.
  --max-outer-extra N     Family outer sizes go up to skew_size + this. Default: 8.
  --exactly-two-rows      Skip one-row shapes in family scans.
  --check-injective       Check that the first carrier target gives an injection.
  --check-matching        Check that the carrier graph has a matching.
  --check-greedy          Check greedy carrier matching.
  --check-descent-output  Check whether unrestricted canonical matching
                          preserves the selected descent mode.
  --check-target-intervals
                          Check whether each source has a lex-interval
                          carrier neighborhood.
  --check-target-swap-closed
                          Check whether each source neighborhood is closed
                          under swapping the two mixed tableaux.
  --check-unordered-target-intervals
                          Check interval neighborhoods after identifying
                          swapped mixed pairs.
  --check-augment-target-increase
                          Check whether every displaced source in the
                          canonical augmenting matching moves to a later
                          target.
  --check-additive-rank   Check whether each source neighborhood lies on one
                          target anti-diagonal, with anti-diagonal index
                          additive in the two source tableaux.
  --check-layered-additive-rank
                          Check the same additive anti-diagonal condition
                          separately for each carrier start level.
  --check-carrier-word-additive-rank
                          Check additive anti-diagonal behavior separately
                          for each full carrier word.
  --augment-potential P   Potential for --check-augment-target-increase:
                          `lex`, `colex`, `sum-right`, `sum-left`,
                          `sum-reverse-left`, `sum-reverse-right`,
                          `right-reverse-left`, `max-min`, or `min-max`.
                          Default: lex.
  --source-order ORDER    Use `lex`, `reverse-lex`, or `fewest-targets`.
                          Default: lex.
  --carrier-order ORDER   Use `short-first` or `long-first`. Default: short-first.
  --terminal-order ORDER  Use `bottom-first` or `top-first`. Default: bottom-first.
  --descent-mode MODE     Preserve no descents, `global`, `componentwise`,
                          `active-global`, or `active-componentwise`. Default: none.
  --max-beta-active N     Restrict to fibers with beta_i + beta_{i+1} <= N.
  --tableau-limit N       Skip shapes whose tableau enumeration exceeds N.
  --help                  Print this help."
}
