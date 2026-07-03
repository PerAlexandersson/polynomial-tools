use std::collections::{BTreeMap, BTreeSet};

type PeakExpansion = BTreeMap<Vec<u32>, i64>;

fn permutations(n: u32) -> Vec<Vec<u32>> {
    fn backtrack(prefix: &mut Vec<u32>, used: &mut [bool], results: &mut Vec<Vec<u32>>) {
        if prefix.len() + 1 == used.len() {
            results.push(prefix.clone());
            return;
        }

        for value in 1..used.len() {
            if used[value] {
                continue;
            }
            used[value] = true;
            prefix.push(value as u32);
            backtrack(prefix, used, results);
            prefix.pop();
            used[value] = false;
        }
    }

    let mut results = Vec::new();
    let mut used = vec![false; n as usize + 1];
    backtrack(&mut Vec::new(), &mut used, &mut results);
    results
}

fn standardize(values: &[u32]) -> Vec<u32> {
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    values
        .iter()
        .map(|value| {
            sorted
                .iter()
                .position(|candidate| candidate == value)
                .expect("subsequence value should occur") as u32
                + 1
        })
        .collect()
}

fn contains_pattern(permutation: &[u32], pattern: &[u32]) -> bool {
    assert_eq!(pattern.len(), 3, "site examples only use patterns in S_3");

    for i in 0..permutation.len() {
        for j in i + 1..permutation.len() {
            for k in j + 1..permutation.len() {
                if standardize(&[permutation[i], permutation[j], permutation[k]]) == pattern {
                    return true;
                }
            }
        }
    }

    false
}

fn avoids_patterns(permutation: &[u32], patterns: &[Vec<u32>]) -> bool {
    patterns
        .iter()
        .all(|pattern| !contains_pattern(permutation, pattern))
}

fn descent_set(permutation: &[u32]) -> BTreeSet<u32> {
    permutation
        .windows(2)
        .enumerate()
        .filter_map(|(index, pair)| {
            if pair[0] > pair[1] {
                Some(index as u32 + 1)
            } else {
                None
            }
        })
        .collect()
}

fn peak_set(descents: &BTreeSet<u32>) -> Vec<u32> {
    descents
        .iter()
        .filter(|&&descent| descent != 1 && !descents.contains(&(descent - 1)))
        .copied()
        .collect()
}

fn pattern_avoiding_peak_expansion(n: u32, patterns: &[Vec<u32>]) -> PeakExpansion {
    let mut expansion = PeakExpansion::new();
    for permutation in permutations(n) {
        if avoids_patterns(&permutation, patterns) {
            let peaks = peak_set(&descent_set(&permutation));
            *expansion.entry(peaks).or_insert(0) += 1;
        }
    }
    expansion
}

fn shifted_cells(lambda: &[u32]) -> Vec<(usize, usize)> {
    let mut cells = Vec::new();
    for (row, &length) in lambda.iter().enumerate() {
        for offset in 0..length as usize {
            cells.push((row, row + offset));
        }
    }
    cells
}

fn is_standard_shifted_tableau(
    cells: &[(usize, usize)],
    assignment: &BTreeMap<(usize, usize), u32>,
) -> bool {
    for &(row, col) in cells {
        if let Some(right) = assignment.get(&(row, col + 1)) {
            if assignment[&(row, col)] >= *right {
                return false;
            }
        }
        if let Some(above) = assignment.get(&(row + 1, col)) {
            if assignment[&(row, col)] >= *above {
                return false;
            }
        }
    }
    true
}

fn shifted_tableau_peak_set(assignment: &BTreeMap<(usize, usize), u32>, n: u32) -> Vec<u32> {
    let mut rows_by_entry = vec![0usize; n as usize + 1];
    for (&(row, _), &entry) in assignment {
        rows_by_entry[entry as usize] = row;
    }

    let descents: BTreeSet<u32> = (1..n)
        .filter(|&entry| rows_by_entry[entry as usize + 1] > rows_by_entry[entry as usize])
        .collect();
    peak_set(&descents)
}

fn schur_q_peak_expansion(lambda: &[u32]) -> PeakExpansion {
    let cells = shifted_cells(lambda);
    let n = cells.len() as u32;
    let mut expansion = PeakExpansion::new();

    for values in permutations(n) {
        let assignment: BTreeMap<(usize, usize), u32> =
            cells.iter().copied().zip(values.into_iter()).collect();
        if is_standard_shifted_tableau(&cells, &assignment) {
            let peaks = shifted_tableau_peak_set(&assignment, n);
            *expansion.entry(peaks).or_insert(0) += 1;
        }
    }

    expansion
}

fn scale(expansion: &PeakExpansion, coefficient: i64) -> PeakExpansion {
    expansion
        .iter()
        .map(|(peaks, count)| (peaks.clone(), coefficient * count))
        .collect()
}

fn add(lhs: &PeakExpansion, rhs: &PeakExpansion) -> PeakExpansion {
    let mut result = lhs.clone();
    for (peaks, count) in rhs {
        *result.entry(peaks.clone()).or_insert(0) += count;
    }
    result.retain(|_, count| *count != 0);
    result
}

fn schur_q_formula(terms: &[(i64, Vec<u32>)]) -> PeakExpansion {
    terms
        .iter()
        .fold(PeakExpansion::new(), |acc, (coeff, lambda)| {
            add(&acc, &scale(&schur_q_peak_expansion(lambda), *coeff))
        })
}

fn format_pattern(pattern: &[u32]) -> String {
    pattern.iter().map(u32::to_string).collect::<String>()
}

fn format_patterns(patterns: &[Vec<u32>]) -> String {
    format!(
        "{{{}}}",
        patterns
            .iter()
            .map(|pattern| format_pattern(pattern))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn format_peak_set(peaks: &[u32]) -> String {
    if peaks.is_empty() {
        "\\emptyset".to_string()
    } else {
        format!(
            "{{{}}}",
            peaks
                .iter()
                .map(u32::to_string)
                .collect::<Vec<_>>()
                .join(",")
        )
    }
}

fn format_peak_expansion(expansion: &PeakExpansion) -> String {
    expansion
        .iter()
        .map(|(peaks, coefficient)| format!("{coefficient}K_{}", format_peak_set(peaks)))
        .collect::<Vec<_>>()
        .join(" + ")
}

fn format_q_terms(terms: &[(i64, Vec<u32>)]) -> String {
    terms
        .iter()
        .map(|(coefficient, lambda)| {
            format!(
                "{coefficient}Q_{{({})}}",
                lambda
                    .iter()
                    .map(u32::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            )
        })
        .collect::<Vec<_>>()
        .join(" + ")
}

fn main() {
    let examples = [
        (vec![vec![1, 2, 3]], vec![(4, vec![4]), (5, vec![3, 1])]),
        (vec![vec![1, 3, 2], vec![2, 3, 1]], vec![(8, vec![4])]),
        (
            vec![vec![1, 2, 3], vec![1, 3, 2], vec![3, 1, 2]],
            vec![(2, vec![4]), (1, vec![3, 1])],
        ),
    ];

    for (patterns, q_terms) in examples {
        let direct = pattern_avoiding_peak_expansion(4, &patterns);
        let formula = schur_q_formula(&q_terms);
        assert_eq!(
            direct,
            formula,
            "pattern-avoiding peak formula failed for Pi={}",
            format_patterns(&patterns)
        );
        println!(
            "Pi={} : R_4 = {} = {}",
            format_patterns(&patterns),
            format_q_terms(&q_terms),
            format_peak_expansion(&direct)
        );
    }
}
