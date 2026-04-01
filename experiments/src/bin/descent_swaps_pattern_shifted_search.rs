//! Search for shifted combinations A + t B inside the promising downward
//! exact-pattern sector on the right-base side.

use combpoly::fixed_descent::{
    consecutive_valid_insertion_position_source_descent_set_data,
    valid_insertion_positions_for_target_descent_set,
};
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::check_weak_interlacing;
use std::collections::BTreeMap;

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
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum DownwardPattern {
    TopDown,
    MaxLeftFar,
    ConsecutiveDown,
    FarDown,
}

impl DownwardPattern {
    fn all() -> [DownwardPattern; 4] {
        [
            DownwardPattern::TopDown,
            DownwardPattern::MaxLeftFar,
            DownwardPattern::ConsecutiveDown,
            DownwardPattern::FarDown,
        ]
    }

    fn short_name(self) -> &'static str {
        match self {
            DownwardPattern::TopDown => "top_down",
            DownwardPattern::MaxLeftFar => "maxL_far",
            DownwardPattern::ConsecutiveDown => "cons_down",
            DownwardPattern::FarDown => "far_down",
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

fn insertion_breaks_consecutive_ascending_pair_at_boundary(
    source_permutation: &[u8],
    insertion_position: u8,
) -> bool {
    if insertion_position == 1 || insertion_position as usize > source_permutation.len() {
        return false;
    }
    let left = source_permutation[(insertion_position - 2) as usize];
    let right = source_permutation[(insertion_position - 1) as usize];
    left + 1 == right
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

fn classify_downward_pattern(
    source_permutation: &[u8],
    endpoint_position: u8,
    n: u8,
) -> Option<DownwardPattern> {
    if endpoint_position == 0 || endpoint_position as usize >= source_permutation.len() {
        return None;
    }
    let left = source_permutation[(endpoint_position - 1) as usize];
    let right = source_permutation[endpoint_position as usize];

    if left == n - 1 && right == n - 2 {
        Some(DownwardPattern::TopDown)
    } else if left == n - 1 {
        Some(DownwardPattern::MaxLeftFar)
    } else if right + 1 == left {
        Some(DownwardPattern::ConsecutiveDown)
    } else if left > right {
        Some(DownwardPattern::FarDown)
    } else {
        None
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

fn family_pattern_poly(
    permutations_by_descent_set: &BTreeMap<u64, Vec<Vec<u8>>>,
    family_descent_set_mask: Option<u64>,
    insertion_position: u8,
    r: u8,
    target_pattern: DownwardPattern,
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
        let Some(pattern) = classify_downward_pattern(pi, r, n) else {
            continue;
        };
        if pattern != target_pattern {
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

fn mask_name(mask: usize) -> String {
    let mut pieces = Vec::new();
    for (i, pattern) in DownwardPattern::all().into_iter().enumerate() {
        if mask & (1 << i) != 0 {
            pieces.push(pattern.short_name());
        }
    }
    if pieces.is_empty() {
        "0".to_string()
    } else {
        pieces.join("+")
    }
}

fn shifted_name(mask_a: usize, mask_b: usize) -> String {
    format!("{} + t*{}", mask_name(mask_a), mask_name(mask_b))
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    let mut p2_counts: BTreeMap<(usize, usize), (usize, usize)> = BTreeMap::new();
    let mut p3_counts: BTreeMap<(usize, usize), (usize, usize)> = BTreeMap::new();
    let mut p4_counts: BTreeMap<(usize, usize), (usize, usize)> = BTreeMap::new();
    let mut s1_counts: BTreeMap<(usize, usize), (usize, usize)> = BTreeMap::new();
    let mut both_p2_p3_counts: BTreeMap<(usize, usize), (usize, usize)> = BTreeMap::new();
    let mut total_cases = 0usize;

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

                let left_augmented_mask = local_data.left_augmented_source_descent_set_mask;
                let right_base_mask = Some(local_data.right_base_source_descent_set_mask);
                let r = local_data.descent_run_end_position;

                let mut left_bin_polys = Vec::new();
                for bin in Bin::all() {
                    left_bin_polys.push(family_bin_poly(
                        &permutations_by_descent_set,
                        left_augmented_mask,
                        p_a,
                        p_a,
                        p_b,
                        r,
                        bin,
                        n,
                    ));
                }
                let left_prefixes = prefix_sums(&left_bin_polys);
                let left_suffixes = suffix_sums(&left_bin_polys);
                let la_p2 = left_prefixes[1].clone();
                let la_p3 = left_prefixes[2].clone();
                let la_p4 = left_prefixes[3].clone();
                let la_s1 = left_suffixes[0].clone();

                let mut downward_polys = Vec::new();
                for pattern in DownwardPattern::all() {
                    downward_polys.push(family_pattern_poly(
                        &permutations_by_descent_set,
                        right_base_mask,
                        p_b,
                        r,
                        pattern,
                        n,
                    ));
                }

                let mut subset_sums = vec![vec![0]; 1 << 4];
                for mask in 1usize..(1 << 4) {
                    let lsb = mask & (!mask + 1);
                    let bit_index = lsb.trailing_zeros() as usize;
                    let prev = mask ^ lsb;
                    subset_sums[mask] = poly_add(&subset_sums[prev], &downward_polys[bit_index]);
                }

                for mask_a in 0usize..(1 << 4) {
                    for mask_b in 0usize..(1 << 4) {
                        if mask_a == 0 && mask_b == 0 {
                            continue;
                        }
                        if mask_a & mask_b != 0 {
                            continue;
                        }
                        let combo =
                            poly_add(&subset_sums[mask_a], &poly_mul_t(&subset_sums[mask_b]));
                        if combo == [0] {
                            continue;
                        }

                        let key = (mask_a, mask_b);

                        if la_p2 != [0] {
                            let entry = p2_counts.entry(key).or_insert((0, 0));
                            entry.0 += 1;
                            if interlaces_weak(&la_p2, &combo) {
                                entry.1 += 1;
                            }
                        }
                        if la_p3 != [0] {
                            let entry = p3_counts.entry(key).or_insert((0, 0));
                            entry.0 += 1;
                            if interlaces_weak(&la_p3, &combo) {
                                entry.1 += 1;
                            }
                        }
                        if la_p4 != [0] {
                            let entry = p4_counts.entry(key).or_insert((0, 0));
                            entry.0 += 1;
                            if interlaces_weak(&la_p4, &combo) {
                                entry.1 += 1;
                            }
                        }
                        if la_s1 != [0] {
                            let entry = s1_counts.entry(key).or_insert((0, 0));
                            entry.0 += 1;
                            if interlaces_weak(&la_s1, &combo) {
                                entry.1 += 1;
                            }
                        }
                        if la_p2 != [0] && la_p3 != [0] {
                            let entry = both_p2_p3_counts.entry(key).or_insert((0, 0));
                            entry.0 += 1;
                            if interlaces_weak(&la_p2, &combo) && interlaces_weak(&la_p3, &combo) {
                                entry.1 += 1;
                            }
                        }
                    }
                }

                total_cases += 1;
            }
        }
    }

    println!("=== Shifted downward-sector search ===");
    println!("total consecutive local cases: {}", total_cases);
    println!();

    let show_top = |title: &str, counts: BTreeMap<(usize, usize), (usize, usize)>| {
        let mut vec: Vec<_> = counts
            .into_iter()
            .filter(|(_, (checked, _))| *checked >= 20)
            .collect();
        vec.sort_by(|a, b| {
            let ar = a.1 .1 as f64 / a.1 .0 as f64;
            let br = b.1 .1 as f64 / b.1 .0 as f64;
            br.partial_cmp(&ar)
                .unwrap()
                .then_with(|| b.1 .0.cmp(&a.1 .0))
                .then_with(|| a.0.cmp(&b.0))
        });
        println!("{}", title);
        for ((mask_a, mask_b), (checked, passed)) in vec.into_iter().take(20) {
            println!(
                "  {:<60} {}/{}",
                shifted_name(mask_a, mask_b),
                passed,
                checked
            );
        }
        println!();
    };

    show_top("Best shifted objects for LA:P2:", p2_counts);
    show_top("Best shifted objects for LA:P3:", p3_counts);
    show_top("Best shifted objects for LA:P4:", p4_counts);
    show_top("Best shifted objects for LA:S1:", s1_counts);
    show_top(
        "Best shifted objects passing both LA:P2 and LA:P3:",
        both_p2_p3_counts,
    );
}
