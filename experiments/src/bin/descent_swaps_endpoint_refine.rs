//! Explore endpoint-relative refinements for the fixed-descent swaps problem.
//!
//! For a consecutive valid pair `p_a < p_b` in `P(S)`, the draft now proves
//! that the varying-source relation is local: the active descent run starting
//! at `p_a` ends at some position `r`, and the source descent sets differ by
//! moving that endpoint.
//!
//! This binary asks whether the natural "peak-like" local refinement is:
//!   - source type: base vs augmented (`S'` vs `S''`),
//!   - previous-maximum position `q = pos(n-1)`,
//!   - coarse endpoint-relative bin of `q`:
//!       FL = q <= p_a-2
//!       RW = p_a-1 <= q <= r
//!       GW = r+1 <= q <= p_b-2
//!       FR = q >= p_b-1
//!
//! The main outputs are:
//!   1. how often the q-refined families are interlacing,
//!   2. how often the endpoint-bin sequence FL<RW<GW<FR is interlacing within
//!      each source family,
//!   3. pairwise interlacing statistics between endpoint-bin state classes.

use combpoly::fixed_descent::{
    augmented_source_descent_set_for_target_descent_set_and_insertion_position,
    base_source_descent_set_for_target_descent_set_and_insertion_position,
    consecutive_valid_insertion_position_source_descent_set_data, descent_positions_from_bitmask,
    insertion_breaks_consecutive_ascending_pair_at_boundary,
    valid_insertion_positions_for_target_descent_set,
};
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::check_weak_interlacing;
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum Family {
    LeftBase,
    LeftAugmented,
    RightBase,
    RightAugmented,
}

impl Family {
    fn all() -> [Family; 4] {
        [
            Family::LeftBase,
            Family::LeftAugmented,
            Family::RightBase,
            Family::RightAugmented,
        ]
    }

    fn short_name(self) -> &'static str {
        match self {
            Family::LeftBase => "LB",
            Family::LeftAugmented => "LA",
            Family::RightBase => "RB",
            Family::RightAugmented => "RA",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum Bin {
    FarLeft,
    RunWindow,
    GapWindow,
    FarRight,
}

impl Bin {
    fn all() -> [Bin; 4] {
        [Bin::FarLeft, Bin::RunWindow, Bin::GapWindow, Bin::FarRight]
    }

    fn short_name(self) -> &'static str {
        match self {
            Bin::FarLeft => "FL",
            Bin::RunWindow => "RW",
            Bin::GapWindow => "GW",
            Bin::FarRight => "FR",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum AuxKind {
    Prefix,
    SuffixLr,
}

impl AuxKind {
    fn all() -> [AuxKind; 2] {
        [AuxKind::Prefix, AuxKind::SuffixLr]
    }

    fn short_name(self) -> &'static str {
        match self {
            AuxKind::Prefix => "P",
            AuxKind::SuffixLr => "S",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum EndpointPairPattern {
    TopAdjacentAscending,
    TopAdjacentDescending,
    PreviousMaxLeftFar,
    PreviousMaxRightFar,
    ConsecutiveAscending,
    ConsecutiveDescending,
    NonconsecutiveAscending,
    NonconsecutiveDescending,
}

impl EndpointPairPattern {
    fn all() -> [EndpointPairPattern; 8] {
        [
            EndpointPairPattern::TopAdjacentAscending,
            EndpointPairPattern::TopAdjacentDescending,
            EndpointPairPattern::PreviousMaxLeftFar,
            EndpointPairPattern::PreviousMaxRightFar,
            EndpointPairPattern::ConsecutiveAscending,
            EndpointPairPattern::ConsecutiveDescending,
            EndpointPairPattern::NonconsecutiveAscending,
            EndpointPairPattern::NonconsecutiveDescending,
        ]
    }

    fn short_name(self) -> &'static str {
        match self {
            EndpointPairPattern::TopAdjacentAscending => "top_up",
            EndpointPairPattern::TopAdjacentDescending => "top_down",
            EndpointPairPattern::PreviousMaxLeftFar => "maxL_far",
            EndpointPairPattern::PreviousMaxRightFar => "maxR_far",
            EndpointPairPattern::ConsecutiveAscending => "cons_up",
            EndpointPairPattern::ConsecutiveDescending => "cons_down",
            EndpointPairPattern::NonconsecutiveAscending => "far_up",
            EndpointPairPattern::NonconsecutiveDescending => "far_down",
        }
    }
}

fn trim(mut p: Vec<i64>) -> Vec<i64> {
    while p.len() > 1 && p.last() == Some(&0) {
        p.pop();
    }
    if p.is_empty() {
        vec![0]
    } else {
        p
    }
}

fn poly_degree(p: &[i64]) -> usize {
    let mut d = p.len();
    while d > 1 && p[d - 1] == 0 {
        d -= 1;
    }
    d - 1
}

fn poly_mul_t(p: &[i64]) -> Vec<i64> {
    let mut out = vec![0; p.len() + 1];
    for (i, &c) in p.iter().enumerate() {
        out[i + 1] = c;
    }
    trim(out)
}

/// Weak interlacing check with the same-degree Wagner reduction.
fn interlaces_weak(f: &[i64], g: &[i64]) -> bool {
    let f = trim(f.to_vec());
    let g = trim(g.to_vec());
    if f == [0] {
        return g == [0] || check_weak_interlacing(&[], &g) == Some(true);
    }
    if g == [0] {
        return true;
    }

    let df = poly_degree(&f);
    let dg = poly_degree(&g);
    if dg == df + 1 {
        check_weak_interlacing(&f, &g) == Some(true)
    } else if dg == df {
        check_weak_interlacing(&g, &poly_mul_t(&f)) == Some(true)
    } else {
        false
    }
}

fn build_poly(values: &[usize]) -> Vec<i64> {
    if values.is_empty() {
        return vec![0];
    }
    let mut p = vec![0; values.iter().copied().max().unwrap() + 1];
    for &v in values {
        p[v] += 1;
    }
    trim(p)
}

fn is_pairwise_interlacing_sequence(polys: &[(String, Vec<i64>)]) -> bool {
    for i in 0..polys.len() {
        for j in i + 1..polys.len() {
            if !interlaces_weak(&polys[i].1, &polys[j].1) {
                return false;
            }
        }
    }
    true
}

fn family_total_poly(
    permutations_by_descent_set: &BTreeMap<u64, Vec<Vec<u8>>>,
    family_descent_set_mask: Option<u64>,
    insertion_position: u8,
) -> Vec<i64> {
    let Some(mask) = family_descent_set_mask else {
        return vec![0];
    };
    let Some(permutations) = permutations_by_descent_set.get(&mask) else {
        return vec![0];
    };
    let values: Vec<usize> = permutations
        .iter()
        .map(|pi| {
            compute(pi, Stat::Swaps)
                + usize::from(insertion_breaks_consecutive_ascending_pair_at_boundary(
                    pi,
                    insertion_position,
                ))
        })
        .collect();
    build_poly(&values)
}

fn family_q_sequence(
    permutations_by_descent_set: &BTreeMap<u64, Vec<Vec<u8>>>,
    family_descent_set_mask: Option<u64>,
    insertion_position: u8,
    n: u8,
) -> Vec<(String, Vec<i64>)> {
    let Some(mask) = family_descent_set_mask else {
        return Vec::new();
    };
    let Some(permutations) = permutations_by_descent_set.get(&mask) else {
        return Vec::new();
    };
    let mut by_q: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
    for pi in permutations {
        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
        let modified_swaps = compute(pi, Stat::Swaps)
            + usize::from(insertion_breaks_consecutive_ascending_pair_at_boundary(
                pi,
                insertion_position,
            ));
        by_q.entry(q).or_default().push(modified_swaps);
    }
    by_q.into_iter()
        .map(|(q, values)| (format!("q={}", q), build_poly(&values)))
        .collect()
}

fn classify_endpoint_relative_bin(q: u8, p_a: u8, p_b: u8, r: u8) -> Bin {
    if p_a >= 2 && q <= p_a - 2 {
        Bin::FarLeft
    } else if q <= r {
        Bin::RunWindow
    } else if p_b >= 2 && q <= p_b - 2 {
        Bin::GapWindow
    } else {
        Bin::FarRight
    }
}

fn family_bin_poly(
    permutations_by_descent_set: &BTreeMap<u64, Vec<Vec<u8>>>,
    family_descent_set_mask: Option<u64>,
    insertion_position: u8,
    p_a: u8,
    p_b: u8,
    r: u8,
    target_bin: Bin,
    n: u8,
) -> Vec<i64> {
    let Some(mask) = family_descent_set_mask else {
        return vec![0];
    };
    let Some(permutations) = permutations_by_descent_set.get(&mask) else {
        return vec![0];
    };
    let mut values = Vec::new();
    for pi in permutations {
        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
        if classify_endpoint_relative_bin(q, p_a, p_b, r) != target_bin {
            continue;
        }
        values.push(
            compute(pi, Stat::Swaps)
                + usize::from(insertion_breaks_consecutive_ascending_pair_at_boundary(
                    pi,
                    insertion_position,
                )),
        );
    }
    build_poly(&values)
}

fn family_bin_poly_with_epsilon1_filter(
    permutations_by_descent_set: &BTreeMap<u64, Vec<Vec<u8>>>,
    family_descent_set_mask: Option<u64>,
    insertion_position: u8,
    p_a: u8,
    p_b: u8,
    r: u8,
    target_bin: Bin,
    required_epsilon1: bool,
    n: u8,
) -> Vec<i64> {
    let Some(mask) = family_descent_set_mask else {
        return vec![0];
    };
    let Some(permutations) = permutations_by_descent_set.get(&mask) else {
        return vec![0];
    };
    let mut values = Vec::new();
    for pi in permutations {
        let epsilon1 =
            insertion_breaks_consecutive_ascending_pair_at_boundary(pi, insertion_position);
        if epsilon1 != required_epsilon1 {
            continue;
        }
        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
        if classify_endpoint_relative_bin(q, p_a, p_b, r) != target_bin {
            continue;
        }
        values.push(compute(pi, Stat::Swaps) + usize::from(epsilon1));
    }
    build_poly(&values)
}

fn endpoint_adjacent_values_are_consecutive(
    source_permutation: &[u8],
    endpoint_position: u8,
) -> bool {
    if endpoint_position == 0 || endpoint_position as usize >= source_permutation.len() {
        return false;
    }
    let left = source_permutation[(endpoint_position - 1) as usize];
    let right = source_permutation[endpoint_position as usize];
    left.abs_diff(right) == 1
}

fn endpoint_adjacent_pair_contains_previous_max(
    source_permutation: &[u8],
    endpoint_position: u8,
    n: u8,
) -> bool {
    if endpoint_position == 0 || endpoint_position as usize >= source_permutation.len() {
        return false;
    }
    let left = source_permutation[(endpoint_position - 1) as usize];
    let right = source_permutation[endpoint_position as usize];
    left == n - 1 || right == n - 1
}

fn endpoint_previous_max_is_on_left(
    source_permutation: &[u8],
    endpoint_position: u8,
    n: u8,
) -> Option<bool> {
    if endpoint_position == 0 || endpoint_position as usize >= source_permutation.len() {
        return None;
    }
    let left = source_permutation[(endpoint_position - 1) as usize];
    let right = source_permutation[endpoint_position as usize];
    if left == n - 1 {
        Some(true)
    } else if right == n - 1 {
        Some(false)
    } else {
        None
    }
}

fn classify_endpoint_pair_pattern(
    source_permutation: &[u8],
    endpoint_position: u8,
    n: u8,
) -> Option<EndpointPairPattern> {
    if endpoint_position == 0 || endpoint_position as usize >= source_permutation.len() {
        return None;
    }
    let left = source_permutation[(endpoint_position - 1) as usize];
    let right = source_permutation[endpoint_position as usize];

    if left == n - 2 && right == n - 1 {
        return Some(EndpointPairPattern::TopAdjacentAscending);
    }
    if left == n - 1 && right == n - 2 {
        return Some(EndpointPairPattern::TopAdjacentDescending);
    }
    if left == n - 1 {
        return Some(EndpointPairPattern::PreviousMaxLeftFar);
    }
    if right == n - 1 {
        return Some(EndpointPairPattern::PreviousMaxRightFar);
    }
    if left + 1 == right {
        return Some(EndpointPairPattern::ConsecutiveAscending);
    }
    if right + 1 == left {
        return Some(EndpointPairPattern::ConsecutiveDescending);
    }
    if left < right {
        Some(EndpointPairPattern::NonconsecutiveAscending)
    } else {
        Some(EndpointPairPattern::NonconsecutiveDescending)
    }
}

fn family_bin_poly_with_endpoint_consecutive_filter(
    permutations_by_descent_set: &BTreeMap<u64, Vec<Vec<u8>>>,
    family_descent_set_mask: Option<u64>,
    insertion_position: u8,
    p_a: u8,
    p_b: u8,
    r: u8,
    target_bin: Bin,
    required_endpoint_consecutive: bool,
    n: u8,
) -> Vec<i64> {
    let Some(mask) = family_descent_set_mask else {
        return vec![0];
    };
    let Some(permutations) = permutations_by_descent_set.get(&mask) else {
        return vec![0];
    };
    let mut values = Vec::new();
    for pi in permutations {
        let endpoint_consecutive = endpoint_adjacent_values_are_consecutive(pi, r);
        if endpoint_consecutive != required_endpoint_consecutive {
            continue;
        }
        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
        if classify_endpoint_relative_bin(q, p_a, p_b, r) != target_bin {
            continue;
        }
        values.push(
            compute(pi, Stat::Swaps)
                + usize::from(insertion_breaks_consecutive_ascending_pair_at_boundary(
                    pi,
                    insertion_position,
                )),
        );
    }
    build_poly(&values)
}

fn family_bin_poly_with_endpoint_local_filters(
    permutations_by_descent_set: &BTreeMap<u64, Vec<Vec<u8>>>,
    family_descent_set_mask: Option<u64>,
    insertion_position: u8,
    p_a: u8,
    p_b: u8,
    r: u8,
    target_bin: Bin,
    required_endpoint_consecutive: bool,
    required_contains_previous_max: bool,
    n: u8,
) -> Vec<i64> {
    let Some(mask) = family_descent_set_mask else {
        return vec![0];
    };
    let Some(permutations) = permutations_by_descent_set.get(&mask) else {
        return vec![0];
    };
    let mut values = Vec::new();
    for pi in permutations {
        let endpoint_consecutive = endpoint_adjacent_values_are_consecutive(pi, r);
        if endpoint_consecutive != required_endpoint_consecutive {
            continue;
        }
        let contains_previous_max = endpoint_adjacent_pair_contains_previous_max(pi, r, n);
        if contains_previous_max != required_contains_previous_max {
            continue;
        }
        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
        if classify_endpoint_relative_bin(q, p_a, p_b, r) != target_bin {
            continue;
        }
        values.push(
            compute(pi, Stat::Swaps)
                + usize::from(insertion_breaks_consecutive_ascending_pair_at_boundary(
                    pi,
                    insertion_position,
                )),
        );
    }
    build_poly(&values)
}

fn family_bin_poly_with_endpoint_side_filters(
    permutations_by_descent_set: &BTreeMap<u64, Vec<Vec<u8>>>,
    family_descent_set_mask: Option<u64>,
    insertion_position: u8,
    p_a: u8,
    p_b: u8,
    r: u8,
    target_bin: Bin,
    required_endpoint_consecutive: bool,
    required_previous_max_on_left: bool,
    n: u8,
) -> Vec<i64> {
    let Some(mask) = family_descent_set_mask else {
        return vec![0];
    };
    let Some(permutations) = permutations_by_descent_set.get(&mask) else {
        return vec![0];
    };
    let mut values = Vec::new();
    for pi in permutations {
        let endpoint_consecutive = endpoint_adjacent_values_are_consecutive(pi, r);
        if endpoint_consecutive != required_endpoint_consecutive {
            continue;
        }
        let Some(previous_max_on_left) = endpoint_previous_max_is_on_left(pi, r, n) else {
            continue;
        };
        if previous_max_on_left != required_previous_max_on_left {
            continue;
        }
        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
        if classify_endpoint_relative_bin(q, p_a, p_b, r) != target_bin {
            continue;
        }
        values.push(
            compute(pi, Stat::Swaps)
                + usize::from(insertion_breaks_consecutive_ascending_pair_at_boundary(
                    pi,
                    insertion_position,
                )),
        );
    }
    build_poly(&values)
}

fn family_bin_poly_with_endpoint_pattern_filter(
    permutations_by_descent_set: &BTreeMap<u64, Vec<Vec<u8>>>,
    family_descent_set_mask: Option<u64>,
    insertion_position: u8,
    p_a: u8,
    p_b: u8,
    r: u8,
    target_bin: Bin,
    required_pattern: EndpointPairPattern,
    n: u8,
) -> Vec<i64> {
    let Some(mask) = family_descent_set_mask else {
        return vec![0];
    };
    let Some(permutations) = permutations_by_descent_set.get(&mask) else {
        return vec![0];
    };
    let mut values = Vec::new();
    for pi in permutations {
        let Some(pattern) = classify_endpoint_pair_pattern(pi, r, n) else {
            continue;
        };
        if pattern != required_pattern {
            continue;
        }
        let q = pi.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
        if classify_endpoint_relative_bin(q, p_a, p_b, r) != target_bin {
            continue;
        }
        values.push(
            compute(pi, Stat::Swaps)
                + usize::from(insertion_breaks_consecutive_ascending_pair_at_boundary(
                    pi,
                    insertion_position,
                )),
        );
    }
    build_poly(&values)
}

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let mut out = vec![0; a.len().max(b.len())];
    for (i, &c) in a.iter().enumerate() {
        out[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        out[i] += c;
    }
    trim(out)
}

fn prefix_sums(polys: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let mut out = Vec::new();
    let mut running = vec![0];
    for poly in polys {
        running = poly_add(&running, poly);
        out.push(running.clone());
    }
    out
}

fn suffix_sums(polys: &[Vec<i64>]) -> Vec<Vec<i64>> {
    let mut out_rev = Vec::new();
    let mut running = vec![0];
    for poly in polys.iter().rev() {
        running = poly_add(&running, poly);
        out_rev.push(running.clone());
    }
    out_rev.reverse();
    out_rev
}

fn permutations_of_families() -> Vec<[Family; 4]> {
    fn rec(pos: usize, arr: &mut [Family; 4], out: &mut Vec<[Family; 4]>) {
        if pos == arr.len() {
            out.push(*arr);
            return;
        }
        for i in pos..arr.len() {
            arr.swap(pos, i);
            rec(pos + 1, arr, out);
            arr.swap(pos, i);
        }
    }
    let mut arr = Family::all();
    let mut out = Vec::new();
    rec(0, &mut arr, &mut out);
    out
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    let family_orders = permutations_of_families();
    let mut order_counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut total_cases = 0usize;
    let mut q_interlacing_counts: BTreeMap<&'static str, (usize, usize)> = BTreeMap::new();
    let mut pair_counts: BTreeMap<&'static str, (usize, usize)> = BTreeMap::new();
    let mut fixed_p_combined_orders: BTreeMap<&'static str, (usize, usize)> = BTreeMap::new();
    let mut bin_pair_counts: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let mut family_bin_sequence_counts: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let mut family_aux_sequence_counts: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let mut aux_pair_counts: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let mut aux_matched_index_counts: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let mut eps_split_aux_counts: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let mut endpoint_split_aux_counts: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let mut endpoint_hybrid_split_aux_counts: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let mut endpoint_side_split_aux_counts: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let mut endpoint_pattern_aux_counts: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let mut sample_failures = Vec::new();

    for n in 4..=max_n {
        let perms_prev = all_permutations(n - 1);
        let mut permutations_by_descent_set: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in perms_prev {
            permutations_by_descent_set
                .entry(descent_set_bitmask(&pi))
                .or_default()
                .push(pi);
        }

        for target_descent_set_mask in 0u64..(1 << (n - 1)) {
            if target_descent_set_mask & 1 != 0 {
                continue;
            }
            let valid_positions =
                valid_insertion_positions_for_target_descent_set(target_descent_set_mask, n);

            for &p in &valid_positions {
                let base_mask = Some(
                    base_source_descent_set_for_target_descent_set_and_insertion_position(
                        target_descent_set_mask,
                        p,
                        n,
                    ),
                );
                let augmented_mask =
                    augmented_source_descent_set_for_target_descent_set_and_insertion_position(
                        target_descent_set_mask,
                        p,
                        n,
                    );
                if augmented_mask.is_none() {
                    continue;
                }

                let base_seq = family_q_sequence(&permutations_by_descent_set, base_mask, p, n);
                let aug_seq = family_q_sequence(&permutations_by_descent_set, augmented_mask, p, n);
                if base_seq.is_empty() || aug_seq.is_empty() {
                    continue;
                }

                let candidate_orders: [(&str, Vec<(String, Vec<i64>)>); 4] = [
                    (
                        "B_up_A_up",
                        base_seq
                            .iter()
                            .cloned()
                            .chain(aug_seq.iter().cloned())
                            .collect(),
                    ),
                    (
                        "A_up_B_up",
                        aug_seq
                            .iter()
                            .cloned()
                            .chain(base_seq.iter().cloned())
                            .collect(),
                    ),
                    (
                        "B_down_A_up",
                        base_seq
                            .iter()
                            .rev()
                            .cloned()
                            .chain(aug_seq.iter().cloned())
                            .collect(),
                    ),
                    (
                        "A_up_B_down",
                        aug_seq
                            .iter()
                            .cloned()
                            .chain(base_seq.iter().rev().cloned())
                            .collect(),
                    ),
                ];

                for (label, seq) in candidate_orders {
                    let entry = fixed_p_combined_orders.entry(label).or_insert((0, 0));
                    entry.0 += 1;
                    if is_pairwise_interlacing_sequence(&seq) {
                        entry.1 += 1;
                    }
                }
            }

            for consecutive_pair in valid_positions.windows(2) {
                let p_a = consecutive_pair[0];
                let p_b = consecutive_pair[1];
                let local_data = consecutive_valid_insertion_position_source_descent_set_data(
                    target_descent_set_mask,
                    p_a,
                    p_b,
                    n,
                )
                .unwrap();

                let family_masks = [
                    (
                        Family::LeftBase,
                        Some(local_data.left_base_source_descent_set_mask),
                        p_a,
                    ),
                    (
                        Family::LeftAugmented,
                        local_data.left_augmented_source_descent_set_mask,
                        p_a,
                    ),
                    (
                        Family::RightBase,
                        Some(local_data.right_base_source_descent_set_mask),
                        p_b,
                    ),
                    (
                        Family::RightAugmented,
                        local_data.right_augmented_source_descent_set_mask,
                        p_b,
                    ),
                ];

                let mut family_totals = BTreeMap::new();
                let mut family_bins = BTreeMap::new();
                let mut family_aux_states: BTreeMap<(Family, AuxKind, usize), Vec<i64>> =
                    BTreeMap::new();
                for &(family, mask, insertion_position) in &family_masks {
                    family_totals.insert(
                        family,
                        family_total_poly(&permutations_by_descent_set, mask, insertion_position),
                    );

                    let q_sequence = family_q_sequence(
                        &permutations_by_descent_set,
                        mask,
                        insertion_position,
                        n,
                    );
                    if q_sequence.len() >= 2 {
                        let entry = q_interlacing_counts
                            .entry(family.short_name())
                            .or_insert((0, 0));
                        entry.0 += 1;
                        if is_pairwise_interlacing_sequence(&q_sequence) {
                            entry.1 += 1;
                        }
                    }

                    let mut bin_sequence = Vec::new();
                    let mut ordered_bin_polys = Vec::new();
                    for bin in Bin::all() {
                        let poly = family_bin_poly(
                            &permutations_by_descent_set,
                            mask,
                            insertion_position,
                            p_a,
                            p_b,
                            local_data.descent_run_end_position,
                            bin,
                            n,
                        );
                        family_bins.insert((family, bin), poly.clone());
                        ordered_bin_polys.push(poly.clone());
                        if poly != [0] {
                            bin_sequence.push((bin.short_name().to_string(), poly));
                        }
                    }
                    if bin_sequence.len() >= 2 {
                        let key = format!("{}:FL<RW<GW<FR", family.short_name());
                        let entry = family_bin_sequence_counts.entry(key).or_insert((0, 0));
                        entry.0 += 1;
                        if is_pairwise_interlacing_sequence(&bin_sequence) {
                            entry.1 += 1;
                        }
                    }

                    let prefix = prefix_sums(&ordered_bin_polys);
                    let prefix_seq: Vec<(String, Vec<i64>)> = prefix
                        .iter()
                        .enumerate()
                        .map(|(i, poly)| (format!("P{}", i + 1), poly.clone()))
                        .collect();
                    let key = format!("{}:prefix", family.short_name());
                    let entry = family_aux_sequence_counts.entry(key).or_insert((0, 0));
                    entry.0 += 1;
                    if is_pairwise_interlacing_sequence(&prefix_seq) {
                        entry.1 += 1;
                    }
                    for (i, poly) in prefix.iter().enumerate() {
                        family_aux_states.insert((family, AuxKind::Prefix, i), poly.clone());
                    }

                    let suffix = suffix_sums(&ordered_bin_polys);
                    let suffix_seq_lr: Vec<(String, Vec<i64>)> = suffix
                        .iter()
                        .enumerate()
                        .map(|(i, poly)| (format!("S{}", i + 1), poly.clone()))
                        .collect();
                    let key = format!("{}:suffix_lr", family.short_name());
                    let entry = family_aux_sequence_counts.entry(key).or_insert((0, 0));
                    entry.0 += 1;
                    if is_pairwise_interlacing_sequence(&suffix_seq_lr) {
                        entry.1 += 1;
                    }
                    for (i, poly) in suffix.iter().enumerate() {
                        family_aux_states.insert((family, AuxKind::SuffixLr, i), poly.clone());
                    }

                    let suffix_seq_rl: Vec<(String, Vec<i64>)> = suffix
                        .into_iter()
                        .enumerate()
                        .rev()
                        .map(|(i, poly)| (format!("S{}", i + 1), poly))
                        .collect();
                    let key = format!("{}:suffix_rl", family.short_name());
                    let entry = family_aux_sequence_counts.entry(key).or_insert((0, 0));
                    entry.0 += 1;
                    if is_pairwise_interlacing_sequence(&suffix_seq_rl) {
                        entry.1 += 1;
                    }
                }

                for &(family, mask, insertion_position) in &family_masks {
                    if !matches!(family, Family::LeftBase | Family::RightBase) {
                        continue;
                    }
                    for epsilon1_value in [false, true] {
                        let mut ordered_bin_polys = Vec::new();
                        for bin in Bin::all() {
                            ordered_bin_polys.push(family_bin_poly_with_epsilon1_filter(
                                &permutations_by_descent_set,
                                mask,
                                insertion_position,
                                p_a,
                                p_b,
                                local_data.descent_run_end_position,
                                bin,
                                epsilon1_value,
                                n,
                            ));
                        }

                        let prefix = prefix_sums(&ordered_bin_polys);
                        let prefix_seq: Vec<(String, Vec<i64>)> = prefix
                            .iter()
                            .enumerate()
                            .map(|(i, poly)| (format!("P{}", i + 1), poly.clone()))
                            .collect();
                        let key = format!(
                            "{}:eps{}:prefix",
                            family.short_name(),
                            if epsilon1_value { 1 } else { 0 }
                        );
                        let entry = eps_split_aux_counts.entry(key).or_insert((0, 0));
                        entry.0 += 1;
                        if is_pairwise_interlacing_sequence(&prefix_seq) {
                            entry.1 += 1;
                        }

                        let suffix = suffix_sums(&ordered_bin_polys);
                        let suffix_seq: Vec<(String, Vec<i64>)> = suffix
                            .iter()
                            .enumerate()
                            .map(|(i, poly)| (format!("S{}", i + 1), poly.clone()))
                            .collect();
                        let key = format!(
                            "{}:eps{}:suffix_lr",
                            family.short_name(),
                            if epsilon1_value { 1 } else { 0 }
                        );
                        let entry = eps_split_aux_counts.entry(key).or_insert((0, 0));
                        entry.0 += 1;
                        if is_pairwise_interlacing_sequence(&suffix_seq) {
                            entry.1 += 1;
                        }

                        if family == Family::RightBase {
                            let la_prefixes: Vec<Vec<i64>> = (0..4)
                                .map(|i| {
                                    family_aux_states
                                        .get(&(Family::LeftAugmented, AuxKind::Prefix, i))
                                        .unwrap()
                                        .clone()
                                })
                                .collect();
                            for (i, rb_poly) in prefix.iter().enumerate() {
                                let la_poly = &la_prefixes[i];
                                if *la_poly == [0] || *rb_poly == [0] {
                                    continue;
                                }
                                let key = format!(
                                    "LA:P{}<<RB_eps{}:P{}",
                                    i + 1,
                                    if epsilon1_value { 1 } else { 0 },
                                    i + 1
                                );
                                let entry = eps_split_aux_counts.entry(key).or_insert((0, 0));
                                entry.0 += 1;
                                if interlaces_weak(la_poly, rb_poly) {
                                    entry.1 += 1;
                                }
                            }

                            let la_suffixes: Vec<Vec<i64>> = (0..4)
                                .map(|i| {
                                    family_aux_states
                                        .get(&(Family::LeftAugmented, AuxKind::SuffixLr, i))
                                        .unwrap()
                                        .clone()
                                })
                                .collect();
                            for (i, rb_poly) in suffix.iter().enumerate() {
                                let la_poly = &la_suffixes[i];
                                if *la_poly == [0] || *rb_poly == [0] {
                                    continue;
                                }
                                let key = format!(
                                    "LA:S{}<<RB_eps{}:S{}",
                                    i + 1,
                                    if epsilon1_value { 1 } else { 0 },
                                    i + 1
                                );
                                let entry = eps_split_aux_counts.entry(key).or_insert((0, 0));
                                entry.0 += 1;
                                if interlaces_weak(la_poly, rb_poly) {
                                    entry.1 += 1;
                                }
                            }
                        }
                    }
                }

                for &(family, mask, insertion_position) in &family_masks {
                    if !matches!(family, Family::LeftAugmented | Family::RightBase) {
                        continue;
                    }
                    for endpoint_consecutive in [false, true] {
                        let mut ordered_bin_polys = Vec::new();
                        for bin in Bin::all() {
                            ordered_bin_polys.push(
                                family_bin_poly_with_endpoint_consecutive_filter(
                                    &permutations_by_descent_set,
                                    mask,
                                    insertion_position,
                                    p_a,
                                    p_b,
                                    local_data.descent_run_end_position,
                                    bin,
                                    endpoint_consecutive,
                                    n,
                                ),
                            );
                        }

                        let prefix = prefix_sums(&ordered_bin_polys);
                        let prefix_seq: Vec<(String, Vec<i64>)> = prefix
                            .iter()
                            .enumerate()
                            .map(|(i, poly)| (format!("P{}", i + 1), poly.clone()))
                            .collect();
                        let key = format!(
                            "{}:end{}:prefix",
                            family.short_name(),
                            if endpoint_consecutive { 1 } else { 0 }
                        );
                        let entry = endpoint_split_aux_counts.entry(key).or_insert((0, 0));
                        entry.0 += 1;
                        if is_pairwise_interlacing_sequence(&prefix_seq) {
                            entry.1 += 1;
                        }

                        let suffix = suffix_sums(&ordered_bin_polys);
                        let suffix_seq: Vec<(String, Vec<i64>)> = suffix
                            .iter()
                            .enumerate()
                            .map(|(i, poly)| (format!("S{}", i + 1), poly.clone()))
                            .collect();
                        let key = format!(
                            "{}:end{}:suffix_lr",
                            family.short_name(),
                            if endpoint_consecutive { 1 } else { 0 }
                        );
                        let entry = endpoint_split_aux_counts.entry(key).or_insert((0, 0));
                        entry.0 += 1;
                        if is_pairwise_interlacing_sequence(&suffix_seq) {
                            entry.1 += 1;
                        }

                        if family == Family::RightBase {
                            let la_prefixes: Vec<Vec<i64>> = (0..4)
                                .map(|i| {
                                    family_aux_states
                                        .get(&(Family::LeftAugmented, AuxKind::Prefix, i))
                                        .unwrap()
                                        .clone()
                                })
                                .collect();
                            for (i, rb_poly) in prefix.iter().enumerate() {
                                let la_poly = &la_prefixes[i];
                                if *la_poly == [0] || *rb_poly == [0] {
                                    continue;
                                }
                                let key = format!(
                                    "LA:P{}<<RB_end{}:P{}",
                                    i + 1,
                                    if endpoint_consecutive { 1 } else { 0 },
                                    i + 1
                                );
                                let entry = endpoint_split_aux_counts.entry(key).or_insert((0, 0));
                                entry.0 += 1;
                                if interlaces_weak(la_poly, rb_poly) {
                                    entry.1 += 1;
                                }
                            }

                            let la_suffixes: Vec<Vec<i64>> = (0..4)
                                .map(|i| {
                                    family_aux_states
                                        .get(&(Family::LeftAugmented, AuxKind::SuffixLr, i))
                                        .unwrap()
                                        .clone()
                                })
                                .collect();
                            for (i, rb_poly) in suffix.iter().enumerate() {
                                let la_poly = &la_suffixes[i];
                                if *la_poly == [0] || *rb_poly == [0] {
                                    continue;
                                }
                                let key = format!(
                                    "LA:S{}<<RB_end{}:S{}",
                                    i + 1,
                                    if endpoint_consecutive { 1 } else { 0 },
                                    i + 1
                                );
                                let entry = endpoint_split_aux_counts.entry(key).or_insert((0, 0));
                                entry.0 += 1;
                                if interlaces_weak(la_poly, rb_poly) {
                                    entry.1 += 1;
                                }
                            }
                        }
                    }
                }

                for &(family, mask, insertion_position) in &family_masks {
                    if !matches!(family, Family::LeftAugmented | Family::RightBase) {
                        continue;
                    }
                    for pattern in EndpointPairPattern::all() {
                        let mut ordered_bin_polys = Vec::new();
                        for bin in Bin::all() {
                            ordered_bin_polys.push(family_bin_poly_with_endpoint_pattern_filter(
                                &permutations_by_descent_set,
                                mask,
                                insertion_position,
                                p_a,
                                p_b,
                                local_data.descent_run_end_position,
                                bin,
                                pattern,
                                n,
                            ));
                        }

                        let prefix = prefix_sums(&ordered_bin_polys);
                        let prefix_seq: Vec<(String, Vec<i64>)> = prefix
                            .iter()
                            .enumerate()
                            .map(|(i, poly)| (format!("P{}", i + 1), poly.clone()))
                            .collect();
                        let key =
                            format!("{}:{}:prefix", family.short_name(), pattern.short_name());
                        let entry = endpoint_pattern_aux_counts.entry(key).or_insert((0, 0));
                        entry.0 += 1;
                        if is_pairwise_interlacing_sequence(&prefix_seq) {
                            entry.1 += 1;
                        }

                        let suffix = suffix_sums(&ordered_bin_polys);
                        let suffix_seq: Vec<(String, Vec<i64>)> = suffix
                            .iter()
                            .enumerate()
                            .map(|(i, poly)| (format!("S{}", i + 1), poly.clone()))
                            .collect();
                        let key =
                            format!("{}:{}:suffix_lr", family.short_name(), pattern.short_name());
                        let entry = endpoint_pattern_aux_counts.entry(key).or_insert((0, 0));
                        entry.0 += 1;
                        if is_pairwise_interlacing_sequence(&suffix_seq) {
                            entry.1 += 1;
                        }

                        if family == Family::RightBase {
                            let la_prefixes: Vec<Vec<i64>> = (0..4)
                                .map(|i| {
                                    family_aux_states
                                        .get(&(Family::LeftAugmented, AuxKind::Prefix, i))
                                        .unwrap()
                                        .clone()
                                })
                                .collect();
                            for (i, rb_poly) in prefix.iter().enumerate() {
                                let la_poly = &la_prefixes[i];
                                if *la_poly == [0] || *rb_poly == [0] {
                                    continue;
                                }
                                let key = format!(
                                    "LA:P{}<<RB_{}:P{}",
                                    i + 1,
                                    pattern.short_name(),
                                    i + 1
                                );
                                let entry =
                                    endpoint_pattern_aux_counts.entry(key).or_insert((0, 0));
                                entry.0 += 1;
                                if interlaces_weak(la_poly, rb_poly) {
                                    entry.1 += 1;
                                }
                            }

                            let la_suffixes: Vec<Vec<i64>> = (0..4)
                                .map(|i| {
                                    family_aux_states
                                        .get(&(Family::LeftAugmented, AuxKind::SuffixLr, i))
                                        .unwrap()
                                        .clone()
                                })
                                .collect();
                            for (i, rb_poly) in suffix.iter().enumerate() {
                                let la_poly = &la_suffixes[i];
                                if *la_poly == [0] || *rb_poly == [0] {
                                    continue;
                                }
                                let key = format!(
                                    "LA:S{}<<RB_{}:S{}",
                                    i + 1,
                                    pattern.short_name(),
                                    i + 1
                                );
                                let entry =
                                    endpoint_pattern_aux_counts.entry(key).or_insert((0, 0));
                                entry.0 += 1;
                                if interlaces_weak(la_poly, rb_poly) {
                                    entry.1 += 1;
                                }
                            }
                        }
                    }
                }

                for &(family, mask, insertion_position) in &family_masks {
                    if !matches!(family, Family::LeftAugmented | Family::RightBase) {
                        continue;
                    }
                    for endpoint_consecutive in [false, true] {
                        for previous_max_on_left in [false, true] {
                            let mut ordered_bin_polys = Vec::new();
                            for bin in Bin::all() {
                                ordered_bin_polys.push(family_bin_poly_with_endpoint_side_filters(
                                    &permutations_by_descent_set,
                                    mask,
                                    insertion_position,
                                    p_a,
                                    p_b,
                                    local_data.descent_run_end_position,
                                    bin,
                                    endpoint_consecutive,
                                    previous_max_on_left,
                                    n,
                                ));
                            }

                            let prefix = prefix_sums(&ordered_bin_polys);
                            let prefix_seq: Vec<(String, Vec<i64>)> = prefix
                                .iter()
                                .enumerate()
                                .map(|(i, poly)| (format!("P{}", i + 1), poly.clone()))
                                .collect();
                            let key = format!(
                                "{}:end{}:side{}:prefix",
                                family.short_name(),
                                if endpoint_consecutive { 1 } else { 0 },
                                if previous_max_on_left { "L" } else { "R" }
                            );
                            let entry = endpoint_side_split_aux_counts.entry(key).or_insert((0, 0));
                            entry.0 += 1;
                            if is_pairwise_interlacing_sequence(&prefix_seq) {
                                entry.1 += 1;
                            }

                            let suffix = suffix_sums(&ordered_bin_polys);
                            let suffix_seq: Vec<(String, Vec<i64>)> = suffix
                                .iter()
                                .enumerate()
                                .map(|(i, poly)| (format!("S{}", i + 1), poly.clone()))
                                .collect();
                            let key = format!(
                                "{}:end{}:side{}:suffix_lr",
                                family.short_name(),
                                if endpoint_consecutive { 1 } else { 0 },
                                if previous_max_on_left { "L" } else { "R" }
                            );
                            let entry = endpoint_side_split_aux_counts.entry(key).or_insert((0, 0));
                            entry.0 += 1;
                            if is_pairwise_interlacing_sequence(&suffix_seq) {
                                entry.1 += 1;
                            }

                            if family == Family::RightBase {
                                let la_prefixes: Vec<Vec<i64>> = (0..4)
                                    .map(|i| {
                                        family_aux_states
                                            .get(&(Family::LeftAugmented, AuxKind::Prefix, i))
                                            .unwrap()
                                            .clone()
                                    })
                                    .collect();
                                for (i, rb_poly) in prefix.iter().enumerate() {
                                    let la_poly = &la_prefixes[i];
                                    if *la_poly == [0] || *rb_poly == [0] {
                                        continue;
                                    }
                                    let key = format!(
                                        "LA:P{}<<RB_end{}:side{}:P{}",
                                        i + 1,
                                        if endpoint_consecutive { 1 } else { 0 },
                                        if previous_max_on_left { "L" } else { "R" },
                                        i + 1
                                    );
                                    let entry =
                                        endpoint_side_split_aux_counts.entry(key).or_insert((0, 0));
                                    entry.0 += 1;
                                    if interlaces_weak(la_poly, rb_poly) {
                                        entry.1 += 1;
                                    }
                                }

                                let la_suffixes: Vec<Vec<i64>> = (0..4)
                                    .map(|i| {
                                        family_aux_states
                                            .get(&(Family::LeftAugmented, AuxKind::SuffixLr, i))
                                            .unwrap()
                                            .clone()
                                    })
                                    .collect();
                                for (i, rb_poly) in suffix.iter().enumerate() {
                                    let la_poly = &la_suffixes[i];
                                    if *la_poly == [0] || *rb_poly == [0] {
                                        continue;
                                    }
                                    let key = format!(
                                        "LA:S{}<<RB_end{}:side{}:S{}",
                                        i + 1,
                                        if endpoint_consecutive { 1 } else { 0 },
                                        if previous_max_on_left { "L" } else { "R" },
                                        i + 1
                                    );
                                    let entry =
                                        endpoint_side_split_aux_counts.entry(key).or_insert((0, 0));
                                    entry.0 += 1;
                                    if interlaces_weak(la_poly, rb_poly) {
                                        entry.1 += 1;
                                    }
                                }
                            }
                        }
                    }
                }

                for &(family, mask, insertion_position) in &family_masks {
                    if !matches!(family, Family::LeftAugmented | Family::RightBase) {
                        continue;
                    }
                    for endpoint_consecutive in [false, true] {
                        for contains_previous_max in [false, true] {
                            let mut ordered_bin_polys = Vec::new();
                            for bin in Bin::all() {
                                ordered_bin_polys.push(
                                    family_bin_poly_with_endpoint_local_filters(
                                        &permutations_by_descent_set,
                                        mask,
                                        insertion_position,
                                        p_a,
                                        p_b,
                                        local_data.descent_run_end_position,
                                        bin,
                                        endpoint_consecutive,
                                        contains_previous_max,
                                        n,
                                    ),
                                );
                            }

                            let prefix = prefix_sums(&ordered_bin_polys);
                            let prefix_seq: Vec<(String, Vec<i64>)> = prefix
                                .iter()
                                .enumerate()
                                .map(|(i, poly)| (format!("P{}", i + 1), poly.clone()))
                                .collect();
                            let key = format!(
                                "{}:end{}:max{}:prefix",
                                family.short_name(),
                                if endpoint_consecutive { 1 } else { 0 },
                                if contains_previous_max { 1 } else { 0 }
                            );
                            let entry = endpoint_hybrid_split_aux_counts
                                .entry(key)
                                .or_insert((0, 0));
                            entry.0 += 1;
                            if is_pairwise_interlacing_sequence(&prefix_seq) {
                                entry.1 += 1;
                            }

                            let suffix = suffix_sums(&ordered_bin_polys);
                            let suffix_seq: Vec<(String, Vec<i64>)> = suffix
                                .iter()
                                .enumerate()
                                .map(|(i, poly)| (format!("S{}", i + 1), poly.clone()))
                                .collect();
                            let key = format!(
                                "{}:end{}:max{}:suffix_lr",
                                family.short_name(),
                                if endpoint_consecutive { 1 } else { 0 },
                                if contains_previous_max { 1 } else { 0 }
                            );
                            let entry = endpoint_hybrid_split_aux_counts
                                .entry(key)
                                .or_insert((0, 0));
                            entry.0 += 1;
                            if is_pairwise_interlacing_sequence(&suffix_seq) {
                                entry.1 += 1;
                            }

                            if family == Family::RightBase {
                                let la_prefixes: Vec<Vec<i64>> = (0..4)
                                    .map(|i| {
                                        family_aux_states
                                            .get(&(Family::LeftAugmented, AuxKind::Prefix, i))
                                            .unwrap()
                                            .clone()
                                    })
                                    .collect();
                                for (i, rb_poly) in prefix.iter().enumerate() {
                                    let la_poly = &la_prefixes[i];
                                    if *la_poly == [0] || *rb_poly == [0] {
                                        continue;
                                    }
                                    let key = format!(
                                        "LA:P{}<<RB_end{}:max{}:P{}",
                                        i + 1,
                                        if endpoint_consecutive { 1 } else { 0 },
                                        if contains_previous_max { 1 } else { 0 },
                                        i + 1
                                    );
                                    let entry = endpoint_hybrid_split_aux_counts
                                        .entry(key)
                                        .or_insert((0, 0));
                                    entry.0 += 1;
                                    if interlaces_weak(la_poly, rb_poly) {
                                        entry.1 += 1;
                                    }
                                }

                                let la_suffixes: Vec<Vec<i64>> = (0..4)
                                    .map(|i| {
                                        family_aux_states
                                            .get(&(Family::LeftAugmented, AuxKind::SuffixLr, i))
                                            .unwrap()
                                            .clone()
                                    })
                                    .collect();
                                for (i, rb_poly) in suffix.iter().enumerate() {
                                    let la_poly = &la_suffixes[i];
                                    if *la_poly == [0] || *rb_poly == [0] {
                                        continue;
                                    }
                                    let key = format!(
                                        "LA:S{}<<RB_end{}:max{}:S{}",
                                        i + 1,
                                        if endpoint_consecutive { 1 } else { 0 },
                                        if contains_previous_max { 1 } else { 0 },
                                        i + 1
                                    );
                                    let entry = endpoint_hybrid_split_aux_counts
                                        .entry(key)
                                        .or_insert((0, 0));
                                    entry.0 += 1;
                                    if interlaces_weak(la_poly, rb_poly) {
                                        entry.1 += 1;
                                    }
                                }
                            }
                        }
                    }
                }

                let mut state_labels: Vec<(String, Vec<i64>)> = Vec::new();
                for family in Family::all() {
                    for bin in Bin::all() {
                        state_labels.push((
                            format!("{}-{}", family.short_name(), bin.short_name()),
                            family_bins.get(&(family, bin)).unwrap().clone(),
                        ));
                    }
                }

                for i in 0..state_labels.len() {
                    for j in i + 1..state_labels.len() {
                        if state_labels[i].1 == [0] || state_labels[j].1 == [0] {
                            continue;
                        }
                        let label = format!("{}<<{}", state_labels[i].0, state_labels[j].0);
                        let entry = bin_pair_counts.entry(label).or_insert((0, 0));
                        entry.0 += 1;
                        if interlaces_weak(&state_labels[i].1, &state_labels[j].1) {
                            entry.1 += 1;
                        }
                    }
                }

                let mut aux_labels: Vec<(String, Vec<i64>)> = Vec::new();
                for family in Family::all() {
                    for aux_kind in AuxKind::all() {
                        for i in 0..4 {
                            aux_labels.push((
                                format!(
                                    "{}-{}{}",
                                    family.short_name(),
                                    aux_kind.short_name(),
                                    i + 1
                                ),
                                family_aux_states
                                    .get(&(family, aux_kind, i))
                                    .unwrap()
                                    .clone(),
                            ));
                        }
                    }
                }
                for i in 0..aux_labels.len() {
                    for j in i + 1..aux_labels.len() {
                        if aux_labels[i].1 == [0] || aux_labels[j].1 == [0] {
                            continue;
                        }
                        let label = format!("{}<<{}", aux_labels[i].0, aux_labels[j].0);
                        let entry = aux_pair_counts.entry(label).or_insert((0, 0));
                        entry.0 += 1;
                        if interlaces_weak(&aux_labels[i].1, &aux_labels[j].1) {
                            entry.1 += 1;
                        }
                    }
                }

                for aux_kind in AuxKind::all() {
                    for i in 0..4 {
                        for family_i in 0..Family::all().len() {
                            for family_j in family_i + 1..Family::all().len() {
                                let f1 = Family::all()[family_i];
                                let f2 = Family::all()[family_j];
                                let p1 = family_aux_states.get(&(f1, aux_kind, i)).unwrap();
                                let p2 = family_aux_states.get(&(f2, aux_kind, i)).unwrap();
                                if *p1 == [0] || *p2 == [0] {
                                    continue;
                                }
                                let label = format!(
                                    "{}:{}{}<<{}:{}{}",
                                    f1.short_name(),
                                    aux_kind.short_name(),
                                    i + 1,
                                    f2.short_name(),
                                    aux_kind.short_name(),
                                    i + 1
                                );
                                let entry = aux_matched_index_counts.entry(label).or_insert((0, 0));
                                entry.0 += 1;
                                if interlaces_weak(p1, p2) {
                                    entry.1 += 1;
                                }
                            }
                        }
                    }
                }

                for &(f1, f2, label) in &[
                    (Family::LeftBase, Family::LeftAugmented, "LB<<LA"),
                    (Family::RightBase, Family::RightAugmented, "RB<<RA"),
                    (Family::LeftBase, Family::RightBase, "LB<<RB"),
                    (Family::LeftAugmented, Family::RightAugmented, "LA<<RA"),
                    (Family::LeftBase, Family::RightAugmented, "LB<<RA"),
                    (Family::LeftAugmented, Family::RightBase, "LA<<RB"),
                ] {
                    let entry = pair_counts.entry(label).or_insert((0, 0));
                    entry.0 += 1;
                    if interlaces_weak(&family_totals[&f1], &family_totals[&f2]) {
                        entry.1 += 1;
                    }
                }

                total_cases += 1;
                for order in &family_orders {
                    let labelled_polys: Vec<(String, Vec<i64>)> = order
                        .iter()
                        .map(|family| {
                            (
                                family.short_name().to_string(),
                                family_totals.get(family).unwrap().clone(),
                            )
                        })
                        .collect();
                    if is_pairwise_interlacing_sequence(&labelled_polys) {
                        let key = order
                            .iter()
                            .map(|family| family.short_name())
                            .collect::<Vec<_>>()
                            .join("<");
                        *order_counts.entry(key).or_insert(0) += 1;
                    }
                }

                if sample_failures.len() < 12 {
                    let lb_la = interlaces_weak(
                        &family_totals[&Family::LeftBase],
                        &family_totals[&Family::LeftAugmented],
                    );
                    let rb_ra = interlaces_weak(
                        &family_totals[&Family::RightBase],
                        &family_totals[&Family::RightAugmented],
                    );
                    let lb_rb = interlaces_weak(
                        &family_totals[&Family::LeftBase],
                        &family_totals[&Family::RightBase],
                    );
                    if !(lb_la && rb_ra && lb_rb) {
                        sample_failures.push(format!(
                            "n={} S={:?} pa={} pb={} r={} | LB={} LA={} RB={} RA={} | order bits: lb_la={} rb_ra={} lb_rb={}",
                            n,
                            descent_positions_from_bitmask(target_descent_set_mask, n),
                            p_a,
                            p_b,
                            local_data.descent_run_end_position,
                            family_totals[&Family::LeftBase]
                                .iter()
                                .enumerate()
                                .filter(|(_, c)| **c != 0)
                                .map(|(k, c)| format!("{}:{}", k, c))
                                .collect::<Vec<_>>()
                                .join(","),
                            family_totals[&Family::LeftAugmented]
                                .iter()
                                .enumerate()
                                .filter(|(_, c)| **c != 0)
                                .map(|(k, c)| format!("{}:{}", k, c))
                                .collect::<Vec<_>>()
                                .join(","),
                            family_totals[&Family::RightBase]
                                .iter()
                                .enumerate()
                                .filter(|(_, c)| **c != 0)
                                .map(|(k, c)| format!("{}:{}", k, c))
                                .collect::<Vec<_>>()
                                .join(","),
                            family_totals[&Family::RightAugmented]
                                .iter()
                                .enumerate()
                                .filter(|(_, c)| **c != 0)
                                .map(|(k, c)| format!("{}:{}", k, c))
                                .collect::<Vec<_>>()
                                .join(","),
                            lb_la,
                            rb_ra,
                            lb_rb
                        ));
                    }
                }
            }
        }
    }

    println!("=== Fixed-descent peak-like local source test ===");
    println!("total consecutive local cases: {}", total_cases);
    println!();

    println!("Pairwise total-polynomial interlacing:");
    for (label, (checked, passed)) in pair_counts {
        println!("  {:>7}: {}/{}", label, passed, checked);
    }
    println!();

    println!("Within-family q-sequence interlacing:");
    for family in Family::all() {
        if let Some((checked, passed)) = q_interlacing_counts.get(family.short_name()) {
            println!("  {:>2}: {}/{}", family.short_name(), passed, checked);
        }
    }
    println!();

    println!("Fixed-p combined base/augmented q-sequence orders:");
    for label in ["B_up_A_up", "A_up_B_up", "B_down_A_up", "A_up_B_down"] {
        if let Some((checked, passed)) = fixed_p_combined_orders.get(label) {
            println!("  {:>11}: {}/{}", label, passed, checked);
        }
    }
    println!();

    println!("Within-family endpoint-bin sequence interlacing:");
    for family in Family::all() {
        let key = format!("{}:FL<RW<GW<FR", family.short_name());
        if let Some((checked, passed)) = family_bin_sequence_counts.get(&key) {
            println!("  {:>14}: {}/{}", key, passed, checked);
        }
    }
    println!();

    println!("Within-family auxiliary prefix/postfix sequences:");
    for family in Family::all() {
        for suffix in ["prefix", "suffix_lr", "suffix_rl"] {
            let key = format!("{}:{}", family.short_name(), suffix);
            if let Some((checked, passed)) = family_aux_sequence_counts.get(&key) {
                println!("  {:>14}: {}/{}", key, passed, checked);
            }
        }
    }
    println!();

    println!("Epsilon1-split auxiliary sequences and LA->RB bridge:");
    let mut eps_split_vec: Vec<_> = eps_split_aux_counts.into_iter().collect();
    eps_split_vec.sort_by(|a, b| {
        let ar = a.1 .1 as f64 / a.1 .0 as f64;
        let br = b.1 .1 as f64 / b.1 .0 as f64;
        br.partial_cmp(&ar)
            .unwrap()
            .then_with(|| b.1 .0.cmp(&a.1 .0))
            .then_with(|| a.0.cmp(&b.0))
    });
    for (label, (checked, passed)) in eps_split_vec.into_iter().take(24) {
        println!("  {:<24} {}/{}", label, passed, checked);
    }
    println!();

    println!("Endpoint-consecutive split auxiliary sequences and LA->RB bridge:");
    let mut endpoint_split_vec: Vec<_> = endpoint_split_aux_counts.into_iter().collect();
    endpoint_split_vec.sort_by(|a, b| {
        let ar = a.1 .1 as f64 / a.1 .0 as f64;
        let br = b.1 .1 as f64 / b.1 .0 as f64;
        br.partial_cmp(&ar)
            .unwrap()
            .then_with(|| b.1 .0.cmp(&a.1 .0))
            .then_with(|| a.0.cmp(&b.0))
    });
    for (label, (checked, passed)) in endpoint_split_vec.into_iter().take(24) {
        println!("  {:<24} {}/{}", label, passed, checked);
    }
    println!();

    println!("Endpoint hybrid split (consecutivity + previous-max touch):");
    let mut endpoint_hybrid_split_vec: Vec<_> =
        endpoint_hybrid_split_aux_counts.into_iter().collect();
    endpoint_hybrid_split_vec.sort_by(|a, b| {
        let ar = a.1 .1 as f64 / a.1 .0 as f64;
        let br = b.1 .1 as f64 / b.1 .0 as f64;
        br.partial_cmp(&ar)
            .unwrap()
            .then_with(|| b.1 .0.cmp(&a.1 .0))
            .then_with(|| a.0.cmp(&b.0))
    });
    for (label, (checked, passed)) in endpoint_hybrid_split_vec.into_iter().take(32) {
        println!("  {:<30} {}/{}", label, passed, checked);
    }
    println!();

    println!("Endpoint side split (consecutivity + side of previous max):");
    let mut endpoint_side_split_vec: Vec<_> = endpoint_side_split_aux_counts.into_iter().collect();
    endpoint_side_split_vec.sort_by(|a, b| {
        let ar = a.1 .1 as f64 / a.1 .0 as f64;
        let br = b.1 .1 as f64 / b.1 .0 as f64;
        br.partial_cmp(&ar)
            .unwrap()
            .then_with(|| b.1 .0.cmp(&a.1 .0))
            .then_with(|| a.0.cmp(&b.0))
    });
    for (label, (checked, passed)) in endpoint_side_split_vec.into_iter().take(32) {
        println!("  {:<30} {}/{}", label, passed, checked);
    }
    println!();

    println!("Endpoint pattern split (exact pair type):");
    let mut endpoint_pattern_vec: Vec<_> = endpoint_pattern_aux_counts.into_iter().collect();
    endpoint_pattern_vec.sort_by(|a, b| {
        let ar = a.1 .1 as f64 / a.1 .0 as f64;
        let br = b.1 .1 as f64 / b.1 .0 as f64;
        br.partial_cmp(&ar)
            .unwrap()
            .then_with(|| b.1 .0.cmp(&a.1 .0))
            .then_with(|| a.0.cmp(&b.0))
    });
    for (label, (checked, passed)) in endpoint_pattern_vec.into_iter().take(40) {
        println!("  {:<30} {}/{}", label, passed, checked);
    }
    println!();

    let mut bin_pair_counts_vec: Vec<_> = bin_pair_counts
        .into_iter()
        .filter(|(_, (checked, _))| *checked >= 20)
        .collect();
    bin_pair_counts_vec.sort_by(|a, b| {
        let ar = a.1 .1 as f64 / a.1 .0 as f64;
        let br = b.1 .1 as f64 / b.1 .0 as f64;
        br.partial_cmp(&ar)
            .unwrap()
            .then_with(|| b.1 .0.cmp(&a.1 .0))
            .then_with(|| a.0.cmp(&b.0))
    });
    println!("Best endpoint-bin pairwise interlacing relations (checked >= 20):");
    for (label, (checked, passed)) in bin_pair_counts_vec.into_iter().take(20) {
        println!("  {:<18} {}/{}", label, passed, checked);
    }
    println!();

    let mut aux_pair_counts_vec: Vec<_> = aux_pair_counts
        .into_iter()
        .filter(|(_, (checked, _))| *checked >= 20)
        .collect();
    aux_pair_counts_vec.sort_by(|a, b| {
        let ar = a.1 .1 as f64 / a.1 .0 as f64;
        let br = b.1 .1 as f64 / b.1 .0 as f64;
        br.partial_cmp(&ar)
            .unwrap()
            .then_with(|| b.1 .0.cmp(&a.1 .0))
            .then_with(|| a.0.cmp(&b.0))
    });
    println!("Best auxiliary-state pairwise interlacing relations (checked >= 20):");
    for (label, (checked, passed)) in aux_pair_counts_vec.into_iter().take(25) {
        println!("  {:<18} {}/{}", label, passed, checked);
    }
    println!();

    let mut aux_matched_index_counts_vec: Vec<_> = aux_matched_index_counts
        .into_iter()
        .filter(|(_, (checked, _))| *checked >= 10)
        .collect();
    aux_matched_index_counts_vec.sort_by(|a, b| {
        let ar = a.1 .1 as f64 / a.1 .0 as f64;
        let br = b.1 .1 as f64 / b.1 .0 as f64;
        br.partial_cmp(&ar)
            .unwrap()
            .then_with(|| b.1 .0.cmp(&a.1 .0))
            .then_with(|| a.0.cmp(&b.0))
    });
    println!("Best matched-index cross-family auxiliary relations (checked >= 10):");
    for (label, (checked, passed)) in aux_matched_index_counts_vec.into_iter().take(25) {
        println!("  {:<20} {}/{}", label, passed, checked);
    }
    println!();

    let mut order_counts_vec: Vec<_> = order_counts.into_iter().collect();
    order_counts_vec.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    println!("Top total-polynomial interlacing orders:");
    for (order, count) in order_counts_vec.into_iter().take(12) {
        println!("  {:<20} {}", order, count);
    }
    println!();

    println!("Sample failures of the naive two-type ordering:");
    for failure in sample_failures {
        println!("  {}", failure);
    }
}
