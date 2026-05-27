use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::enumeration::{Content, Count, TableauRecord};
use crate::shape::SkewShape;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DefectColumnData(Vec<usize>);

impl DefectColumnData {
    pub fn columns(&self) -> &[usize] {
        &self.0
    }
}

impl fmt::Display for DefectColumnData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let columns: Vec<String> = self.0.iter().map(|col| (col + 1).to_string()).collect();
        write!(f, "{{{}}}", columns.join(","))
    }
}

pub type DefectColumnCounts = BTreeMap<DefectColumnData, Count>;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DefectColumnCountData(usize);

impl DefectColumnCountData {
    pub fn count(&self) -> usize {
        self.0
    }
}

impl fmt::Display for DefectColumnCountData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

pub type DefectColumnNumberCounts = BTreeMap<DefectColumnCountData, Count>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DefectColumnFailure {
    pub beta: Content,
    /// Zero-indexed first defect label.
    pub i: usize,
    /// Zero-indexed second defect label.
    pub j: usize,
    pub first_statistic: DefectColumnData,
    pub second_statistic: DefectColumnData,
    pub negative_count: Count,
    pub positive_count: Count,
}

impl fmt::Display for DefectColumnFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "beta={:?}, labels=({},{}), column-stats=({},{}), negative={}, positive={}",
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
pub struct DefectColumnCountFailure {
    pub beta: Content,
    /// Zero-indexed first defect label.
    pub i: usize,
    /// Zero-indexed second defect label.
    pub j: usize,
    pub first_statistic: DefectColumnCountData,
    pub second_statistic: DefectColumnCountData,
    pub negative_count: Count,
    pub positive_count: Count,
}

impl fmt::Display for DefectColumnCountFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "beta={:?}, labels=({},{}), column-counts=({},{}), negative={}, positive={}",
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
pub struct PairDefectColumnFailure {
    pub beta: Content,
    /// Zero-indexed first defect label.
    pub i: usize,
    /// Zero-indexed second defect label.
    pub j: usize,
    pub statistic: DefectColumnData,
    pub negative_count: Count,
    pub positive_count: Count,
}

impl fmt::Display for PairDefectColumnFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "beta={:?}, labels=({},{}), pair-column-stat={}, negative={}, positive={}",
            self.beta,
            self.i + 1,
            self.j + 1,
            self.statistic,
            self.negative_count,
            self.positive_count
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairDefectColumnCountFailure {
    pub beta: Content,
    /// Zero-indexed first defect label.
    pub i: usize,
    /// Zero-indexed second defect label.
    pub j: usize,
    pub statistic: DefectColumnCountData,
    pub negative_count: Count,
    pub positive_count: Count,
}

impl fmt::Display for PairDefectColumnCountFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "beta={:?}, labels=({},{}), pair-column-count={}, negative={}, positive={}",
            self.beta,
            self.i + 1,
            self.j + 1,
            self.statistic,
            self.negative_count,
            self.positive_count
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DefectColumnScanReport {
    pub tests_checked: usize,
    pub failure: Option<DefectColumnFailure>,
}

impl DefectColumnScanReport {
    pub fn passed(&self) -> bool {
        self.failure.is_none()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DefectColumnCountScanReport {
    pub tests_checked: usize,
    pub failure: Option<DefectColumnCountFailure>,
}

impl DefectColumnCountScanReport {
    pub fn passed(&self) -> bool {
        self.failure.is_none()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairDefectColumnScanReport {
    pub tests_checked: usize,
    pub failure: Option<PairDefectColumnFailure>,
}

impl PairDefectColumnScanReport {
    pub fn passed(&self) -> bool {
        self.failure.is_none()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairDefectColumnCountScanReport {
    pub tests_checked: usize,
    pub failure: Option<PairDefectColumnCountFailure>,
}

impl PairDefectColumnCountScanReport {
    pub fn passed(&self) -> bool {
        self.failure.is_none()
    }
}

pub fn columns_containing_both_defects(
    shape: &SkewShape,
    values: &[u32],
    i: usize,
    j: usize,
) -> DefectColumnData {
    let i_value = i as u32 + 1;
    let j_value = j as u32 + 1;
    let mut column_flags = BTreeMap::<usize, u8>::new();

    for (cell, &value) in shape.cells().iter().zip(values) {
        let flag = match value {
            value if value == i_value => 1,
            value if value == j_value => 2,
            _ => 0,
        };
        if flag != 0 {
            *column_flags.entry(cell.col).or_insert(0) |= flag;
        }
    }

    DefectColumnData(
        column_flags
            .into_iter()
            .filter_map(|(col, flags)| (flags == 3).then_some(col))
            .collect(),
    )
}

pub fn number_of_columns_containing_both_defects(
    shape: &SkewShape,
    values: &[u32],
    i: usize,
    j: usize,
) -> DefectColumnCountData {
    DefectColumnCountData(columns_containing_both_defects(shape, values, i, j).0.len())
}

pub fn pair_columns_containing_both_defects(
    shape: &SkewShape,
    left_values: &[u32],
    right_values: &[u32],
    i: usize,
    j: usize,
) -> DefectColumnData {
    let i_value = i as u32 + 1;
    let j_value = j as u32 + 1;
    let mut column_flags = BTreeMap::<usize, u8>::new();

    for (idx, cell) in shape.cells().iter().enumerate() {
        for value in [left_values[idx], right_values[idx]] {
            let flag = match value {
                value if value == i_value => 1,
                value if value == j_value => 2,
                _ => 0,
            };
            if flag != 0 {
                *column_flags.entry(cell.col).or_insert(0) |= flag;
            }
        }
    }

    DefectColumnData(
        column_flags
            .into_iter()
            .filter_map(|(col, flags)| (flags == 3).then_some(col))
            .collect(),
    )
}

pub fn pair_number_of_columns_containing_both_defects(
    shape: &SkewShape,
    left_values: &[u32],
    right_values: &[u32],
    i: usize,
    j: usize,
) -> DefectColumnCountData {
    DefectColumnCountData(
        pair_columns_containing_both_defects(shape, left_values, right_values, i, j)
            .0
            .len(),
    )
}

pub fn check_two_by_two_defect_column_fibers(
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    alphabet: usize,
) -> DefectColumnScanReport {
    let mut tests_checked = 0usize;
    let mut seen = BTreeSet::new();

    for tableau in tableaux {
        for i in 0..alphabet {
            for j in (i + 1)..alphabet {
                let Some(beta) = subtract_unit_pair(&tableau.content, i, j) else {
                    continue;
                };
                if !seen.insert((beta.clone(), i, j)) {
                    continue;
                }

                let ii = add_unit_pair(&beta, i, i);
                let ij = add_unit_pair(&beta, i, j);
                let jj = add_unit_pair(&beta, j, j);

                let ii_counts = defect_column_counts_for_content(shape, tableaux, &ii, i, j);
                let ij_counts = defect_column_counts_for_content(shape, tableaux, &ij, i, j);
                let jj_counts = defect_column_counts_for_content(shape, tableaux, &jj, i, j);

                if ii_counts.is_empty() || ij_counts.is_empty() || jj_counts.is_empty() {
                    continue;
                }

                tests_checked += 1;
                if let Some(failure) =
                    check_single_two_by_two_fiber(&beta, i, j, &ii_counts, &ij_counts, &jj_counts)
                {
                    return DefectColumnScanReport {
                        tests_checked,
                        failure: Some(failure),
                    };
                }
            }
        }
    }

    DefectColumnScanReport {
        tests_checked,
        failure: None,
    }
}

pub fn check_two_by_two_pair_defect_column_fibers(
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    alphabet: usize,
) -> PairDefectColumnScanReport {
    let mut tests_checked = 0usize;
    let mut seen = BTreeSet::new();

    for tableau in tableaux {
        for i in 0..alphabet {
            for j in (i + 1)..alphabet {
                let Some(beta) = subtract_unit_pair(&tableau.content, i, j) else {
                    continue;
                };
                if !seen.insert((beta.clone(), i, j)) {
                    continue;
                }

                let ii = add_unit_pair(&beta, i, i);
                let ij = add_unit_pair(&beta, i, j);
                let jj = add_unit_pair(&beta, j, j);
                let ii_tableaux = tableaux_for_content(tableaux, &ii);
                let ij_tableaux = tableaux_for_content(tableaux, &ij);
                let jj_tableaux = tableaux_for_content(tableaux, &jj);

                if ii_tableaux.is_empty() || ij_tableaux.is_empty() || jj_tableaux.is_empty() {
                    continue;
                }

                tests_checked += 1;
                let negative_counts =
                    pair_defect_column_counts(shape, &ii_tableaux, &jj_tableaux, i, j);
                let positive_counts =
                    pair_defect_column_counts(shape, &ij_tableaux, &ij_tableaux, i, j);
                if let Some(failure) =
                    check_pair_fiber_counts(&beta, i, j, &negative_counts, &positive_counts)
                {
                    return PairDefectColumnScanReport {
                        tests_checked,
                        failure: Some(failure),
                    };
                }
            }
        }
    }

    PairDefectColumnScanReport {
        tests_checked,
        failure: None,
    }
}

pub fn check_two_by_two_defect_column_count_fibers(
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    alphabet: usize,
) -> DefectColumnCountScanReport {
    let mut tests_checked = 0usize;
    let mut seen = BTreeSet::new();

    for tableau in tableaux {
        for i in 0..alphabet {
            for j in (i + 1)..alphabet {
                let Some(beta) = subtract_unit_pair(&tableau.content, i, j) else {
                    continue;
                };
                if !seen.insert((beta.clone(), i, j)) {
                    continue;
                }

                let ii = add_unit_pair(&beta, i, i);
                let ij = add_unit_pair(&beta, i, j);
                let jj = add_unit_pair(&beta, j, j);

                let ii_counts = defect_column_count_counts_for_content(shape, tableaux, &ii, i, j);
                let ij_counts = defect_column_count_counts_for_content(shape, tableaux, &ij, i, j);
                let jj_counts = defect_column_count_counts_for_content(shape, tableaux, &jj, i, j);

                if ii_counts.is_empty() || ij_counts.is_empty() || jj_counts.is_empty() {
                    continue;
                }

                tests_checked += 1;
                if let Some(failure) = check_single_two_by_two_count_fiber(
                    &beta, i, j, &ii_counts, &ij_counts, &jj_counts,
                ) {
                    return DefectColumnCountScanReport {
                        tests_checked,
                        failure: Some(failure),
                    };
                }
            }
        }
    }

    DefectColumnCountScanReport {
        tests_checked,
        failure: None,
    }
}

pub fn check_two_by_two_pair_defect_column_count_fibers(
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    alphabet: usize,
) -> PairDefectColumnCountScanReport {
    let mut tests_checked = 0usize;
    let mut seen = BTreeSet::new();

    for tableau in tableaux {
        for i in 0..alphabet {
            for j in (i + 1)..alphabet {
                let Some(beta) = subtract_unit_pair(&tableau.content, i, j) else {
                    continue;
                };
                if !seen.insert((beta.clone(), i, j)) {
                    continue;
                }

                let ii = add_unit_pair(&beta, i, i);
                let ij = add_unit_pair(&beta, i, j);
                let jj = add_unit_pair(&beta, j, j);
                let ii_tableaux = tableaux_for_content(tableaux, &ii);
                let ij_tableaux = tableaux_for_content(tableaux, &ij);
                let jj_tableaux = tableaux_for_content(tableaux, &jj);

                if ii_tableaux.is_empty() || ij_tableaux.is_empty() || jj_tableaux.is_empty() {
                    continue;
                }

                tests_checked += 1;
                let negative_counts =
                    pair_defect_column_count_counts(shape, &ii_tableaux, &jj_tableaux, i, j);
                let positive_counts =
                    pair_defect_column_count_counts(shape, &ij_tableaux, &ij_tableaux, i, j);
                if let Some(failure) =
                    check_pair_count_fiber_counts(&beta, i, j, &negative_counts, &positive_counts)
                {
                    return PairDefectColumnCountScanReport {
                        tests_checked,
                        failure: Some(failure),
                    };
                }
            }
        }
    }

    PairDefectColumnCountScanReport {
        tests_checked,
        failure: None,
    }
}

pub fn defect_column_counts_for_content(
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    content: &[u32],
    i: usize,
    j: usize,
) -> DefectColumnCounts {
    let mut counts = DefectColumnCounts::new();
    for tableau in tableaux
        .iter()
        .filter(|tableau| tableau.content.as_slice() == content)
    {
        let statistic = columns_containing_both_defects(shape, &tableau.values, i, j);
        *counts.entry(statistic).or_insert(0) += 1;
    }
    counts
}

fn tableaux_for_content<'a>(
    tableaux: &'a [TableauRecord],
    content: &[u32],
) -> Vec<&'a TableauRecord> {
    tableaux
        .iter()
        .filter(|tableau| tableau.content.as_slice() == content)
        .collect()
}

fn pair_defect_column_counts(
    shape: &SkewShape,
    left_tableaux: &[&TableauRecord],
    right_tableaux: &[&TableauRecord],
    i: usize,
    j: usize,
) -> DefectColumnCounts {
    let mut counts = DefectColumnCounts::new();
    for left in left_tableaux {
        for right in right_tableaux {
            let statistic =
                pair_columns_containing_both_defects(shape, &left.values, &right.values, i, j);
            *counts.entry(statistic).or_insert(0) += 1;
        }
    }
    counts
}

fn pair_defect_column_count_counts(
    shape: &SkewShape,
    left_tableaux: &[&TableauRecord],
    right_tableaux: &[&TableauRecord],
    i: usize,
    j: usize,
) -> DefectColumnNumberCounts {
    let mut counts = DefectColumnNumberCounts::new();
    for left in left_tableaux {
        for right in right_tableaux {
            let statistic = pair_number_of_columns_containing_both_defects(
                shape,
                &left.values,
                &right.values,
                i,
                j,
            );
            *counts.entry(statistic).or_insert(0) += 1;
        }
    }
    counts
}

pub fn defect_column_count_counts_for_content(
    shape: &SkewShape,
    tableaux: &[TableauRecord],
    content: &[u32],
    i: usize,
    j: usize,
) -> DefectColumnNumberCounts {
    let mut counts = DefectColumnNumberCounts::new();
    for tableau in tableaux
        .iter()
        .filter(|tableau| tableau.content.as_slice() == content)
    {
        let statistic = number_of_columns_containing_both_defects(shape, &tableau.values, i, j);
        *counts.entry(statistic).or_insert(0) += 1;
    }
    counts
}

fn check_single_two_by_two_fiber(
    beta: &[u32],
    i: usize,
    j: usize,
    ii_counts: &DefectColumnCounts,
    ij_counts: &DefectColumnCounts,
    jj_counts: &DefectColumnCounts,
) -> Option<DefectColumnFailure> {
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
                return Some(DefectColumnFailure {
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

fn check_pair_fiber_counts(
    beta: &[u32],
    i: usize,
    j: usize,
    negative_counts: &DefectColumnCounts,
    positive_counts: &DefectColumnCounts,
) -> Option<PairDefectColumnFailure> {
    let mut statistics = BTreeSet::new();
    statistics.extend(negative_counts.keys().cloned());
    statistics.extend(positive_counts.keys().cloned());

    for statistic in statistics {
        let negative_count = count_for(negative_counts, &statistic);
        let positive_count = count_for(positive_counts, &statistic);
        if negative_count > positive_count {
            return Some(PairDefectColumnFailure {
                beta: beta.to_vec(),
                i,
                j,
                statistic,
                negative_count,
                positive_count,
            });
        }
    }

    None
}

fn check_single_two_by_two_count_fiber(
    beta: &[u32],
    i: usize,
    j: usize,
    ii_counts: &DefectColumnNumberCounts,
    ij_counts: &DefectColumnNumberCounts,
    jj_counts: &DefectColumnNumberCounts,
) -> Option<DefectColumnCountFailure> {
    let mut statistics = BTreeSet::new();
    statistics.extend(ii_counts.keys().cloned());
    statistics.extend(ij_counts.keys().cloned());
    statistics.extend(jj_counts.keys().cloned());
    let statistics: Vec<_> = statistics.into_iter().collect();

    for (left_idx, left_statistic) in statistics.iter().enumerate() {
        for right_statistic in statistics.iter().skip(left_idx) {
            let negative_count =
                unordered_count_pair_count(ii_counts, jj_counts, left_statistic, right_statistic);
            let positive_count =
                unordered_count_pair_count(ij_counts, ij_counts, left_statistic, right_statistic);
            if negative_count > positive_count {
                return Some(DefectColumnCountFailure {
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

fn check_pair_count_fiber_counts(
    beta: &[u32],
    i: usize,
    j: usize,
    negative_counts: &DefectColumnNumberCounts,
    positive_counts: &DefectColumnNumberCounts,
) -> Option<PairDefectColumnCountFailure> {
    let mut statistics = BTreeSet::new();
    statistics.extend(negative_counts.keys().cloned());
    statistics.extend(positive_counts.keys().cloned());

    for statistic in statistics {
        let negative_count = count_for_count(negative_counts, &statistic);
        let positive_count = count_for_count(positive_counts, &statistic);
        if negative_count > positive_count {
            return Some(PairDefectColumnCountFailure {
                beta: beta.to_vec(),
                i,
                j,
                statistic,
                negative_count,
                positive_count,
            });
        }
    }

    None
}

fn unordered_pair_count(
    left_counts: &DefectColumnCounts,
    right_counts: &DefectColumnCounts,
    left_statistic: &DefectColumnData,
    right_statistic: &DefectColumnData,
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

fn count_for(counts: &DefectColumnCounts, statistic: &DefectColumnData) -> Count {
    counts.get(statistic).copied().unwrap_or(0)
}

fn unordered_count_pair_count(
    left_counts: &DefectColumnNumberCounts,
    right_counts: &DefectColumnNumberCounts,
    left_statistic: &DefectColumnCountData,
    right_statistic: &DefectColumnCountData,
) -> Count {
    let left_right = count_for_count(left_counts, left_statistic)
        * count_for_count(right_counts, right_statistic);
    if left_statistic == right_statistic {
        left_right
    } else {
        left_right
            + count_for_count(left_counts, right_statistic)
                * count_for_count(right_counts, left_statistic)
    }
}

fn count_for_count(counts: &DefectColumnNumberCounts, statistic: &DefectColumnCountData) -> Count {
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
    use crate::enumeration::enumerate_tableaux;
    use crate::shape::{RowFlaggedSkewShape, SkewShape};

    #[test]
    fn records_columns_containing_both_labels() {
        let shape = SkewShape::from_parts(vec![2, 2], vec![]);
        let values = vec![1, 2, 2, 3];
        let statistic = columns_containing_both_defects(&shape, &values, 0, 1);
        assert_eq!(statistic.columns(), &[0]);
    }

    #[test]
    fn counts_columns_containing_both_labels() {
        let shape = SkewShape::from_parts(vec![2, 2], vec![]);
        let values = vec![1, 2, 2, 3];
        let statistic = number_of_columns_containing_both_defects(&shape, &values, 0, 1);
        assert_eq!(statistic.count(), 1);
    }

    #[test]
    fn records_pair_columns_containing_both_labels() {
        let shape = SkewShape::from_parts(vec![2, 2], vec![]);
        let left_values = vec![1, 2, 2, 3];
        let right_values = vec![2, 1, 3, 2];
        let statistic =
            pair_columns_containing_both_defects(&shape, &left_values, &right_values, 0, 1);
        assert_eq!(statistic.columns(), &[0, 1]);
    }

    #[test]
    fn tiny_shape_has_defect_column_fiber_check() {
        let shape = SkewShape::from_parts(vec![2, 1], vec![]);
        let flagged = RowFlaggedSkewShape::ordinary(shape.clone(), 3);
        let tableaux = enumerate_tableaux(&flagged, None).unwrap();
        let report = check_two_by_two_defect_column_fibers(&shape, &tableaux, 3);
        assert!(report.passed(), "{:?}", report.failure);
    }
}
