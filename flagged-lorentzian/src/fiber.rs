use std::collections::BTreeSet;
use std::fmt;

use crate::descent::DescentData;
use crate::enumeration::{Content, ContentStatisticCounts, Count, StatisticCounts};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiberFailure {
    pub beta: Content,
    /// Zero-indexed first defect label.
    pub i: usize,
    /// Zero-indexed second defect label.
    pub j: usize,
    pub first_statistic: DescentData,
    pub second_statistic: DescentData,
    pub negative_count: Count,
    pub positive_count: Count,
}

impl fmt::Display for FiberFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "beta={:?}, labels=({},{}), stats=({},{}), negative={}, positive={}",
            self.beta,
            self.i + 1,
            self.j + 1,
            self.first_statistic,
            self.second_statistic,
            self.negative_count,
            self.positive_count
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiberScanReport {
    pub tests_checked: usize,
    pub failure: Option<FiberFailure>,
}

impl FiberScanReport {
    pub fn passed(&self) -> bool {
        self.failure.is_none()
    }
}

pub fn check_two_by_two_fiber_inequalities(
    counts: &ContentStatisticCounts,
    alphabet: usize,
) -> FiberScanReport {
    let mut tests_checked = 0usize;
    let mut seen = BTreeSet::new();

    for mixed_content in counts.keys() {
        for i in 0..alphabet {
            for j in (i + 1)..alphabet {
                let Some(beta) = subtract_unit_pair(mixed_content, i, j) else {
                    continue;
                };
                if !seen.insert((beta.clone(), i, j)) {
                    continue;
                }

                let ii = add_unit_pair(&beta, i, i);
                let ij = add_unit_pair(&beta, i, j);
                let jj = add_unit_pair(&beta, j, j);

                if !(counts.contains_key(&ii)
                    && counts.contains_key(&ij)
                    && counts.contains_key(&jj))
                {
                    continue;
                }

                tests_checked += 1;
                if let Some(failure) = check_two_by_two_fiber_at(counts, &beta, i, j) {
                    return FiberScanReport {
                        tests_checked,
                        failure: Some(failure),
                    };
                }
            }
        }
    }

    FiberScanReport {
        tests_checked,
        failure: None,
    }
}

/// Check one fixed \(2\times2\) fiber inequality.
///
/// The indices `i` and `j` are zero-indexed defect labels.
pub fn check_two_by_two_fiber_at(
    counts: &ContentStatisticCounts,
    beta: &[u32],
    i: usize,
    j: usize,
) -> Option<FiberFailure> {
    let ii = add_unit_pair(beta, i, i);
    let ij = add_unit_pair(beta, i, j);
    let jj = add_unit_pair(beta, j, j);
    let (Some(ii_counts), Some(ij_counts), Some(jj_counts)) =
        (counts.get(&ii), counts.get(&ij), counts.get(&jj))
    else {
        return None;
    };
    check_single_two_by_two_fiber(beta, i, j, ii_counts, ij_counts, jj_counts)
}

fn check_single_two_by_two_fiber(
    beta: &[u32],
    i: usize,
    j: usize,
    ii_counts: &StatisticCounts,
    ij_counts: &StatisticCounts,
    jj_counts: &StatisticCounts,
) -> Option<FiberFailure> {
    let mut statistics = BTreeSet::new();
    statistics.extend(ii_counts.keys().cloned());
    statistics.extend(ij_counts.keys().cloned());
    statistics.extend(jj_counts.keys().cloned());
    let statistics: Vec<_> = statistics.into_iter().collect();

    for (left_idx, left_statistic) in statistics.iter().enumerate() {
        for right_statistic in statistics.iter().skip(left_idx) {
            let negative_count =
                unordered_pair_count(ii_counts, jj_counts, left_statistic, right_statistic);
            let positive_count =
                unordered_pair_count(ij_counts, ij_counts, left_statistic, right_statistic);
            if negative_count > positive_count {
                return Some(FiberFailure {
                    beta: beta.to_vec(),
                    i,
                    j,
                    first_statistic: left_statistic.clone(),
                    second_statistic: right_statistic.clone(),
                    negative_count,
                    positive_count,
                });
            }
        }
    }

    None
}

fn unordered_pair_count(
    left_counts: &StatisticCounts,
    right_counts: &StatisticCounts,
    left_statistic: &DescentData,
    right_statistic: &DescentData,
) -> Count {
    let left_right =
        count_for(left_counts, left_statistic) * count_for(right_counts, right_statistic);
    if left_statistic == right_statistic {
        left_right
    } else {
        left_right
            + count_for(left_counts, right_statistic) * count_for(right_counts, left_statistic)
    }
}

fn count_for(counts: &StatisticCounts, statistic: &DescentData) -> Count {
    counts.get(statistic).copied().unwrap_or(0)
}

fn subtract_unit_pair(content: &[u32], i: usize, j: usize) -> Option<Content> {
    let mut beta = content.to_vec();
    *beta.get_mut(i)? = beta[i].checked_sub(1)?;
    *beta.get_mut(j)? = beta[j].checked_sub(1)?;
    Some(beta)
}

fn add_unit_pair(beta: &[u32], i: usize, j: usize) -> Content {
    let mut content = beta.to_vec();
    content[i] += 1;
    content[j] += 1;
    content
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::descent::DescentStatistic;
    use crate::enumeration::{enumerate_content_statistic_counts, EnumerationOptions};
    use crate::shape::{RowFlaggedSkewShape, SkewShape};

    #[test]
    fn global_descents_fail_for_disconnected_counterexample() {
        let shape = SkewShape::from_parts(vec![4, 3, 1], vec![3, 1]);
        let flagged = RowFlaggedSkewShape::ordinary(shape, 5);
        let counts = enumerate_content_statistic_counts(
            &flagged,
            EnumerationOptions::new(DescentStatistic::Global),
        )
        .unwrap();
        let report = check_two_by_two_fiber_inequalities(&counts, 5);
        assert!(report.failure.is_some());
        let failure = check_two_by_two_fiber_at(&counts, &[2, 0, 0, 0, 0], 0, 1)
            .expect("expected specific global descent failure");
        assert_eq!(failure.beta, vec![2, 0, 0, 0, 0]);
        assert_eq!((failure.i, failure.j), (0, 1));
        assert_eq!(failure.negative_count, 1);
        assert_eq!(failure.positive_count, 0);
    }

    #[test]
    fn componentwise_descents_repair_disconnected_counterexample() {
        let shape = SkewShape::from_parts(vec![4, 3, 1], vec![3, 1]);
        let flagged = RowFlaggedSkewShape::ordinary(shape, 5);
        let counts = enumerate_content_statistic_counts(
            &flagged,
            EnumerationOptions::new(DescentStatistic::Componentwise),
        )
        .unwrap();
        let report = check_two_by_two_fiber_inequalities(&counts, 5);
        assert!(report.passed(), "{:?}", report.failure);
    }
}
