use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use flagged_lorentzian::{
    active_subword_descent_data_for_values, add_patterns, descent_data_for_values,
    elementary_row_exchange_neighbors, enumerate_tableaux, families::skew_shapes_of_size,
    pair_envelope, sharp_flag, DescentData, DescentStatistic, RowFlaggedSkewShape, SkewGtPattern,
    SkewShape,
};

type Content = Vec<u32>;
type GtSum = Vec<Vec<u32>>;
type ActiveRow = Vec<u32>;
type Envelope = Vec<u32>;
type DescentPair = (DescentData, DescentData);

#[derive(Debug, Clone)]
struct Args {
    lambda: Option<Vec<u32>>,
    mu: Vec<u32>,
    row_flags: Option<Vec<u32>>,
    alphabet: usize,
    lower_label: u32,
    invariant_level: Option<u32>,
    upper_invariant: bool,
    max_exchange_depth: usize,
    max_skew_size: u32,
    max_outer_extra: u32,
    max_rows: Option<usize>,
    min_beta_active: Option<u32>,
    max_beta_active: Option<u32>,
    connected_only: bool,
    matching_mode: MatchingMode,
    envelope_mode: EnvelopeMode,
    descent_mode: DescentMode,
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
            invariant_level: None,
            upper_invariant: false,
            max_exchange_depth: 1,
            max_skew_size: 5,
            max_outer_extra: 7,
            max_rows: None,
            min_beta_active: None,
            max_beta_active: None,
            connected_only: false,
            matching_mode: MatchingMode::Flow,
            envelope_mode: EnvelopeMode::Exact,
            descent_mode: DescentMode::None,
            tableau_limit: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MatchingMode {
    Flow,
    Greedy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EnvelopeMode {
    Exact,
    Nonincrease,
    Ignore,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DescentMode {
    None,
    Global,
    Componentwise,
    ActiveGlobal,
    ActiveComponentwise,
}

#[derive(Debug, Clone)]
struct SingleData {
    gt: SkewGtPattern,
    sharp_flag: Envelope,
    descent: Option<DescentData>,
}

#[derive(Debug, Clone)]
struct PairData {
    active_row: ActiveRow,
    envelope: Envelope,
    descent_pair: Option<DescentPair>,
    gt_sum: GtSum,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct ExchangeKey {
    active_row: ActiveRow,
    envelope: Envelope,
    descent_pair: Option<DescentPair>,
    gt_sum: GtSum,
}

impl PairData {
    fn exchange_key(&self) -> ExchangeKey {
        ExchangeKey {
            active_row: self.active_row.clone(),
            envelope: self.envelope.clone(),
            descent_pair: self.descent_pair.clone(),
            gt_sum: self.gt_sum.clone(),
        }
    }
}

#[derive(Debug, Clone)]
struct FiberFailure {
    beta: Content,
    negative_pairs: usize,
    positive_pairs: usize,
    matched_pairs: usize,
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
        std::process::exit(if outcome.failure.is_some() { 1 } else { 0 });
    }

    let mut shapes_checked = 0usize;
    let mut skipped_by_limit = 0usize;
    let mut fibers_checked = 0usize;
    let mut negative_pairs = 0usize;
    let mut positive_pairs = 0usize;

    for skew_size in 0..=args.max_skew_size {
        let max_outer_size = skew_size + args.max_outer_extra;
        for shape in skew_shapes_of_size(skew_size, max_outer_size) {
            if let Some(max_rows) = args.max_rows {
                if shape.row_count() > max_rows {
                    continue;
                }
            }
            if args.connected_only && !shape.is_connected() {
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
            positive_pairs += outcome.positive_pairs;

            if let Some(failure) = outcome.failure {
                println!("FAIL");
                println!("outer={:?}", flagged_shape.shape().outer().parts());
                println!("inner={:?}", flagged_shape.shape().inner().parts());
                println!("flags={:?}", flagged_shape.row_flags());
                println!(
                    "beta={:?}, negative_pairs={}, positive_pairs={}, matched_pairs={}",
                    failure.beta,
                    failure.negative_pairs,
                    failure.positive_pairs,
                    failure.matched_pairs
                );
                println!(
                    "shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, fibers_checked={fibers_checked}, negative_pairs={negative_pairs}, positive_pairs={positive_pairs}, elapsed={:.3}s",
                    started.elapsed().as_secs_f64()
                );
                std::process::exit(1);
            }
        }
        println!(
            "skew_size={skew_size}: shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, fibers_checked={fibers_checked}, negative_pairs={negative_pairs}, positive_pairs={positive_pairs}, elapsed={:.3}s",
            started.elapsed().as_secs_f64()
        );
    }

    println!("PASS");
    println!(
        "shapes_checked={shapes_checked}, skipped_by_limit={skipped_by_limit}, fibers_checked={fibers_checked}, negative_pairs={negative_pairs}, positive_pairs={positive_pairs}, elapsed={:.3}s",
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
                "--invariant-level" => {
                    args.invariant_level = Some(
                        take_value(&mut iter, &flag)?
                            .parse()
                            .map_err(|err| format!("invalid --invariant-level: {err}"))?,
                    )
                }
                "--upper-invariant" => args.upper_invariant = true,
                "--max-exchange-depth" => {
                    args.max_exchange_depth = take_value(&mut iter, &flag)?
                        .parse()
                        .map_err(|err| format!("invalid --max-exchange-depth: {err}"))?
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
                "--max-rows" => {
                    args.max_rows = Some(
                        take_value(&mut iter, &flag)?
                            .parse()
                            .map_err(|err| format!("invalid --max-rows: {err}"))?,
                    )
                }
                "--min-beta-active" => {
                    args.min_beta_active = Some(
                        take_value(&mut iter, &flag)?
                            .parse()
                            .map_err(|err| format!("invalid --min-beta-active: {err}"))?,
                    )
                }
                "--max-beta-active" => {
                    args.max_beta_active = Some(
                        take_value(&mut iter, &flag)?
                            .parse()
                            .map_err(|err| format!("invalid --max-beta-active: {err}"))?,
                    )
                }
                "--connected-only" => args.connected_only = true,
                "--matching-mode" => {
                    args.matching_mode = MatchingMode::parse(&take_value(&mut iter, &flag)?)?
                }
                "--envelope-mode" => {
                    args.envelope_mode = EnvelopeMode::parse(&take_value(&mut iter, &flag)?)?
                }
                "--active-row-only" => args.envelope_mode = EnvelopeMode::Ignore,
                "--descent-mode" => {
                    args.descent_mode = DescentMode::parse(&take_value(&mut iter, &flag)?)?
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
        if args.invariant_level().is_none() {
            return Err("--invariant-level must be at most --alphabet".to_string());
        }
        Ok(args)
    }

    fn invariant_level(&self) -> Option<u32> {
        let level = if self.upper_invariant {
            self.lower_label + 1
        } else {
            self.invariant_level.unwrap_or(self.lower_label)
        };
        (level as usize <= self.alphabet).then_some(level)
    }
}

impl MatchingMode {
    fn parse(input: &str) -> Result<Self, String> {
        match input {
            "flow" => Ok(Self::Flow),
            "greedy" => Ok(Self::Greedy),
            other => Err(format!(
                "invalid --matching-mode `{other}`; expected `flow` or `greedy`"
            )),
        }
    }
}

impl EnvelopeMode {
    fn parse(input: &str) -> Result<Self, String> {
        match input {
            "exact" => Ok(Self::Exact),
            "nonincrease" => Ok(Self::Nonincrease),
            "ignore" | "none" => Ok(Self::Ignore),
            other => Err(format!(
                "invalid --envelope-mode `{other}`; expected `exact`, `nonincrease`, or `ignore`"
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

#[derive(Debug, Clone)]
struct ShapeOutcome {
    fibers_checked: usize,
    negative_pairs: usize,
    positive_pairs: usize,
    failure: Option<FiberFailure>,
    limit_exceeded: bool,
}

fn scan_one_shape(flagged_shape: &RowFlaggedSkewShape, args: &Args) -> ShapeOutcome {
    let tableaux = match enumerate_tableaux(flagged_shape, args.tableau_limit) {
        Ok(tableaux) => tableaux,
        Err(_) => {
            return ShapeOutcome {
                fibers_checked: 0,
                negative_pairs: 0,
                positive_pairs: 0,
                failure: None,
                limit_exceeded: true,
            }
        }
    };
    let reading_orders = args.descent_mode.reading_orders(flagged_shape.shape());
    let single_data: Vec<_> = tableaux
        .iter()
        .map(|tableau| SingleData {
            gt: SkewGtPattern::from_tableau(flagged_shape.shape(), &tableau.values, args.alphabet),
            sharp_flag: sharp_flag(flagged_shape.shape(), &tableau.values),
            descent: reading_orders.as_ref().and_then(|reading_orders| {
                args.descent_mode
                    .descent_data(&tableau.values, reading_orders, args.lower_label)
            }),
        })
        .collect();

    let mut by_content = BTreeMap::<Content, Vec<usize>>::new();
    for (idx, tableau) in tableaux.iter().enumerate() {
        by_content
            .entry(tableau.content.clone())
            .or_default()
            .push(idx);
    }

    let lower = args.lower_label as usize - 1;
    let upper = lower + 1;
    let mut seen_beta = BTreeSet::new();
    let mut fibers_checked = 0usize;
    let mut negative_pairs = 0usize;
    let mut positive_pairs = 0usize;

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
        if let Some(min_beta_active) = args.min_beta_active {
            if beta[lower] + beta[upper] < min_beta_active {
                continue;
            }
        }

        let (Some(ii_indices), Some(mixed_indices), Some(jj_indices)) = (
            by_content.get(&ii),
            by_content.get(&mixed),
            by_content.get(&jj),
        ) else {
            continue;
        };

        fibers_checked += 1;
        let invariant_level = args.invariant_level().expect("validated invariant level");
        let negative = pair_data(&single_data, ii_indices, jj_indices, invariant_level);
        let positive = pair_data(&single_data, mixed_indices, mixed_indices, invariant_level);
        negative_pairs += negative.len();
        positive_pairs += positive.len();
        let matched = exchange_matching_size(
            &negative,
            &positive,
            invariant_level as usize,
            args.max_exchange_depth,
            args.matching_mode,
            args.envelope_mode,
        );
        if matched < negative.len() {
            return ShapeOutcome {
                fibers_checked,
                negative_pairs,
                positive_pairs,
                failure: Some(FiberFailure {
                    beta,
                    negative_pairs: negative.len(),
                    positive_pairs: positive.len(),
                    matched_pairs: matched,
                }),
                limit_exceeded: false,
            };
        }
    }

    ShapeOutcome {
        fibers_checked,
        negative_pairs,
        positive_pairs,
        failure: None,
        limit_exceeded: false,
    }
}

fn pair_data(
    single_data: &[SingleData],
    left_indices: &[usize],
    right_indices: &[usize],
    lower_label: u32,
) -> Vec<PairData> {
    let mut out = Vec::new();
    for &left_idx in left_indices {
        for &right_idx in right_indices {
            let left = &single_data[left_idx];
            let right = &single_data[right_idx];
            let gt_sum = add_patterns(&left.gt, &right.gt);
            out.push(PairData {
                active_row: gt_sum[lower_label as usize].clone(),
                envelope: pair_envelope(&left.sharp_flag, &right.sharp_flag),
                descent_pair: unordered_descent_pair(left.descent.as_ref(), right.descent.as_ref()),
                gt_sum,
            });
        }
    }
    out
}

fn exchange_matching_size(
    negative: &[PairData],
    positive: &[PairData],
    fixed_level: usize,
    max_depth: usize,
    matching_mode: MatchingMode,
    envelope_mode: EnvelopeMode,
) -> usize {
    match matching_mode {
        MatchingMode::Flow => {
            flow_exchange_matching_size(negative, positive, fixed_level, max_depth, envelope_mode)
        }
        MatchingMode::Greedy => {
            greedy_exchange_matching_size(negative, positive, fixed_level, max_depth, envelope_mode)
        }
    }
}

fn flow_exchange_matching_size(
    negative: &[PairData],
    positive: &[PairData],
    fixed_level: usize,
    max_depth: usize,
    envelope_mode: EnvelopeMode,
) -> usize {
    let negative_counts = key_counts(negative);
    let positive_counts = key_counts(positive);
    let negative_keys: Vec<_> = negative_counts.keys().cloned().collect();
    let positive_keys: Vec<_> = positive_counts.keys().cloned().collect();
    let mut positive_groups = BTreeMap::<(ActiveRow, Option<DescentPair>), Vec<usize>>::new();
    for (idx, key) in positive_keys.iter().enumerate() {
        positive_groups
            .entry((key.active_row.clone(), key.descent_pair.clone()))
            .or_default()
            .push(idx);
    }

    let source = 0usize;
    let negative_offset = 1usize;
    let positive_offset = negative_offset + negative_keys.len();
    let sink = positive_offset + positive_keys.len();
    let mut flow = Dinic::new(sink + 1);
    let total_negative = negative.len();

    for (idx, key) in negative_keys.iter().enumerate() {
        flow.add_edge(source, negative_offset + idx, negative_counts[key]);
    }
    for (idx, key) in positive_keys.iter().enumerate() {
        flow.add_edge(positive_offset + idx, sink, positive_counts[key]);
    }

    let mut memo = BTreeMap::<(GtSum, GtSum), bool>::new();
    for (neg_idx, neg) in negative_keys.iter().enumerate() {
        let Some(candidate_indices) =
            positive_groups.get(&(neg.active_row.clone(), neg.descent_pair.clone()))
        else {
            continue;
        };
        for &pos_idx in candidate_indices {
            let pos = &positive_keys[pos_idx];
            if !envelope_allowed(&neg.envelope, &pos.envelope, envelope_mode) {
                continue;
            }
            if gt_l1_distance(&neg.gt_sum, &pos.gt_sum) > 2 * max_depth as u32 {
                continue;
            }
            let reachability_key = (neg.gt_sum.clone(), pos.gt_sum.clone());
            let reachable = match memo.get(&reachability_key) {
                Some(&reachable) => reachable,
                None => {
                    let reachable =
                        exchange_reachable(&neg.gt_sum, &pos.gt_sum, fixed_level, max_depth);
                    memo.insert(reachability_key, reachable);
                    reachable
                }
            };
            if reachable {
                flow.add_edge(
                    negative_offset + neg_idx,
                    positive_offset + pos_idx,
                    total_negative,
                );
            }
        }
    }

    flow.max_flow(source, sink)
}

fn greedy_exchange_matching_size(
    negative: &[PairData],
    positive: &[PairData],
    fixed_level: usize,
    max_depth: usize,
    envelope_mode: EnvelopeMode,
) -> usize {
    let negative_counts = key_counts(negative);
    let positive_counts = key_counts(positive);
    let negative_keys: Vec<_> = negative_counts.keys().cloned().collect();
    let positive_keys: Vec<_> = positive_counts.keys().cloned().collect();
    let mut positive_remaining: Vec<_> = positive_keys
        .iter()
        .map(|key| positive_counts[key])
        .collect();
    let mut positive_groups = BTreeMap::<(ActiveRow, Option<DescentPair>), Vec<usize>>::new();
    for (idx, key) in positive_keys.iter().enumerate() {
        positive_groups
            .entry((key.active_row.clone(), key.descent_pair.clone()))
            .or_default()
            .push(idx);
    }

    let mut memo = BTreeMap::<(GtSum, GtSum), bool>::new();
    let mut matched = 0usize;
    for neg in negative_keys {
        let mut needed = negative_counts[&neg];
        let Some(candidate_indices) =
            positive_groups.get(&(neg.active_row.clone(), neg.descent_pair.clone()))
        else {
            continue;
        };
        for &pos_idx in candidate_indices {
            if needed == 0 {
                break;
            }
            if positive_remaining[pos_idx] == 0 {
                continue;
            }
            let pos = &positive_keys[pos_idx];
            if !envelope_allowed(&neg.envelope, &pos.envelope, envelope_mode) {
                continue;
            }
            if gt_l1_distance(&neg.gt_sum, &pos.gt_sum) > 2 * max_depth as u32 {
                continue;
            }
            let reachability_key = (neg.gt_sum.clone(), pos.gt_sum.clone());
            let reachable = match memo.get(&reachability_key) {
                Some(&reachable) => reachable,
                None => {
                    let reachable =
                        exchange_reachable(&neg.gt_sum, &pos.gt_sum, fixed_level, max_depth);
                    memo.insert(reachability_key, reachable);
                    reachable
                }
            };
            if reachable {
                let take = needed.min(positive_remaining[pos_idx]);
                positive_remaining[pos_idx] -= take;
                needed -= take;
                matched += take;
            }
        }
    }
    matched
}

fn envelope_allowed(negative: &[u32], positive: &[u32], mode: EnvelopeMode) -> bool {
    match mode {
        EnvelopeMode::Exact => negative == positive,
        EnvelopeMode::Nonincrease => positive.iter().zip(negative).all(|(&p, &n)| p <= n),
        EnvelopeMode::Ignore => true,
    }
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

fn gt_l1_distance(left: &[Vec<u32>], right: &[Vec<u32>]) -> u32 {
    left.iter()
        .zip(right)
        .flat_map(|(left, right)| left.iter().zip(right))
        .map(|(&left, &right)| left.abs_diff(right))
        .sum()
}

fn key_counts(pairs: &[PairData]) -> BTreeMap<ExchangeKey, usize> {
    let mut counts = BTreeMap::new();
    for pair in pairs {
        *counts.entry(pair.exchange_key()).or_insert(0) += 1;
    }
    counts
}

fn exchange_reachable(
    start: &[Vec<u32>],
    target: &[Vec<u32>],
    fixed_level: usize,
    max_depth: usize,
) -> bool {
    if start == target {
        return true;
    }
    let mut seen = BTreeSet::from([start.to_vec()]);
    let mut frontier = vec![start.to_vec()];
    for _ in 1..=max_depth {
        let mut next = Vec::new();
        for array in frontier {
            for neighbor in elementary_row_exchange_neighbors(&array, fixed_level) {
                if !seen.insert(neighbor.clone()) {
                    continue;
                }
                if neighbor == target {
                    return true;
                }
                next.push(neighbor);
            }
        }
        frontier = next;
    }
    false
}

#[derive(Debug, Clone)]
struct FlowEdge {
    to: usize,
    rev: usize,
    cap: usize,
}

#[derive(Debug, Clone)]
struct Dinic {
    graph: Vec<Vec<FlowEdge>>,
}

impl Dinic {
    fn new(vertex_count: usize) -> Self {
        Self {
            graph: vec![Vec::new(); vertex_count],
        }
    }

    fn add_edge(&mut self, from: usize, to: usize, cap: usize) {
        if cap == 0 {
            return;
        }
        let from_rev = self.graph[to].len();
        let to_rev = self.graph[from].len();
        self.graph[from].push(FlowEdge {
            to,
            rev: from_rev,
            cap,
        });
        self.graph[to].push(FlowEdge {
            to: from,
            rev: to_rev,
            cap: 0,
        });
    }

    fn max_flow(&mut self, source: usize, sink: usize) -> usize {
        let mut flow = 0usize;
        loop {
            let levels = self.levels(source);
            if levels[sink].is_none() {
                return flow;
            }
            let mut next_edge = vec![0usize; self.graph.len()];
            loop {
                let pushed = self.push(source, sink, usize::MAX, &levels, &mut next_edge);
                if pushed == 0 {
                    break;
                }
                flow += pushed;
            }
        }
    }

    fn levels(&self, source: usize) -> Vec<Option<usize>> {
        let mut levels = vec![None; self.graph.len()];
        let mut queue = std::collections::VecDeque::from([source]);
        levels[source] = Some(0);
        while let Some(vertex) = queue.pop_front() {
            let next_level = levels[vertex].unwrap() + 1;
            for edge in &self.graph[vertex] {
                if edge.cap > 0 && levels[edge.to].is_none() {
                    levels[edge.to] = Some(next_level);
                    queue.push_back(edge.to);
                }
            }
        }
        levels
    }

    fn push(
        &mut self,
        vertex: usize,
        sink: usize,
        flow: usize,
        levels: &[Option<usize>],
        next_edge: &mut [usize],
    ) -> usize {
        if vertex == sink {
            return flow;
        }
        while next_edge[vertex] < self.graph[vertex].len() {
            let edge_idx = next_edge[vertex];
            let edge = self.graph[vertex][edge_idx].clone();
            if edge.cap > 0 && levels[edge.to] == levels[vertex].map(|level| level + 1) {
                let pushed = self.push(edge.to, sink, flow.min(edge.cap), levels, next_edge);
                if pushed > 0 {
                    self.graph[vertex][edge_idx].cap -= pushed;
                    let rev = edge.rev;
                    self.graph[edge.to][rev].cap += pushed;
                    return pushed;
                }
            }
            next_edge[vertex] += 1;
        }
        0
    }
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
    println!("positive_pairs={}", outcome.positive_pairs);
    println!("elapsed={elapsed_seconds:.3}s");
    if let Some(failure) = &outcome.failure {
        println!("FAIL");
        println!(
            "beta={:?}, negative_pairs={}, positive_pairs={}, matched_pairs={}",
            failure.beta, failure.negative_pairs, failure.positive_pairs, failure.matched_pairs
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

fn help_text() -> &'static str {
    "Check matching using bounded elementary GT row-exchange reachability.

USAGE:
  gt_exchange_matching_scan [OPTIONS]

OPTIONS:
  --lambda PARTS          Outer partition, e.g. 4,2,1. If omitted, scan a family.
  --mu PARTS              Inner partition, e.g. 2. Default: empty.
  --row-flags PARTS       Row upper flags. Omit for ordinary skew tableaux.
  --alphabet N            Alphabet size. Default: 5.
  --lower-label I         Use active row for labels I,I+1. Default: 4.
  --invariant-level R     Preserve pair-sum GT row R. Default: lower label.
  --upper-invariant       Alias for `--invariant-level I+1`.
  --max-exchange-depth N  Maximum elementary GT exchange depth. Default: 1.
  --max-skew-size N       Maximum skew size for family scans. Default: 5.
  --max-outer-extra N     Family outer sizes go up to skew_size + this. Default: 7.
  --max-rows N            Restrict family scans to shapes with at most N rows.
  --min-beta-active N     Restrict to fibers with beta_i + beta_{i+1} >= N.
  --max-beta-active N     Restrict to fibers with beta_i + beta_{i+1} <= N.
  --connected-only        Restrict family scans to connected skew shapes.
  --matching-mode MODE    Use `flow` or deterministic `greedy`. Default: flow.
  --envelope-mode MODE    Use `exact`, `nonincrease`, or `ignore`. Default: exact.
  --active-row-only       Alias for `--envelope-mode ignore`.
  --descent-mode MODE     Preserve no descents, `global`, `componentwise`,
                          `active-global`, or `active-componentwise`. Default: none.
  --tableau-limit N       Skip shapes whose tableau enumeration exceeds N.
  --help                  Print this help."
}
