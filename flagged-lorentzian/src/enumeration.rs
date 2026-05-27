use std::collections::BTreeMap;
use std::error::Error;
use std::fmt;

use crate::descent::{descent_data_for_values, DescentData, DescentStatistic};
use crate::shape::{Cell, RowFlaggedSkewShape};

pub type Content = Vec<u32>;
pub type Count = u128;
pub type StatisticCounts = BTreeMap<DescentData, Count>;
pub type ContentStatisticCounts = BTreeMap<Content, StatisticCounts>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TableauRecord {
    pub values: Vec<u32>,
    pub content: Content,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EnumerationOptions {
    pub statistic: DescentStatistic,
    pub tableau_limit: Option<usize>,
}

impl EnumerationOptions {
    pub fn new(statistic: DescentStatistic) -> Self {
        Self {
            statistic,
            tableau_limit: None,
        }
    }

    pub fn with_tableau_limit(mut self, tableau_limit: Option<usize>) -> Self {
        self.tableau_limit = tableau_limit;
        self
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnumerationLimitExceeded {
    pub limit: usize,
}

impl fmt::Display for EnumerationLimitExceeded {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "tableau enumeration exceeded limit {}", self.limit)
    }
}

impl Error for EnumerationLimitExceeded {}

pub fn enumerate_content_statistic_counts(
    flagged_shape: &RowFlaggedSkewShape,
    options: EnumerationOptions,
) -> Result<ContentStatisticCounts, EnumerationLimitExceeded> {
    let shape = flagged_shape.shape();
    let reading_orders = options.statistic.reading_orders(shape);
    let mut counts = ContentStatisticCounts::new();
    let mut values = vec![0u32; shape.size()];
    let mut content = vec![0u32; flagged_shape.alphabet()];
    let mut tableau_count = 0usize;

    enumerate_from_cell(
        flagged_shape,
        &reading_orders,
        &mut values,
        &mut content,
        &mut counts,
        &mut tableau_count,
        options.tableau_limit,
        0,
    )?;

    Ok(counts)
}

pub fn enumerate_tableaux(
    flagged_shape: &RowFlaggedSkewShape,
    tableau_limit: Option<usize>,
) -> Result<Vec<TableauRecord>, EnumerationLimitExceeded> {
    let shape = flagged_shape.shape();
    let mut tableaux = Vec::new();
    let mut values = vec![0u32; shape.size()];
    let mut content = vec![0u32; flagged_shape.alphabet()];

    enumerate_tableaux_from_cell(
        flagged_shape,
        &mut values,
        &mut content,
        &mut tableaux,
        tableau_limit,
        0,
    )?;

    Ok(tableaux)
}

#[allow(clippy::too_many_arguments)]
fn enumerate_from_cell(
    flagged_shape: &RowFlaggedSkewShape,
    reading_orders: &[Vec<usize>],
    values: &mut [u32],
    content: &mut [u32],
    counts: &mut ContentStatisticCounts,
    tableau_count: &mut usize,
    tableau_limit: Option<usize>,
    cell_index: usize,
) -> Result<(), EnumerationLimitExceeded> {
    let shape = flagged_shape.shape();
    if cell_index == shape.size() {
        *tableau_count += 1;
        if let Some(limit) = tableau_limit {
            if *tableau_count > limit {
                return Err(EnumerationLimitExceeded { limit });
            }
        }

        let descent_data = descent_data_for_values(values, reading_orders);
        *counts
            .entry(content.to_vec())
            .or_default()
            .entry(descent_data)
            .or_insert(0) += 1;
        return Ok(());
    }

    let cell = shape.cells()[cell_index];
    let lower_bound = lower_bound(flagged_shape, values, cell);
    let upper_bound = flagged_shape.max_entry_for_row(cell.row);

    for value in lower_bound..=upper_bound {
        values[cell_index] = value;
        content[value as usize - 1] += 1;
        enumerate_from_cell(
            flagged_shape,
            reading_orders,
            values,
            content,
            counts,
            tableau_count,
            tableau_limit,
            cell_index + 1,
        )?;
        content[value as usize - 1] -= 1;
        values[cell_index] = 0;
    }

    Ok(())
}

fn enumerate_tableaux_from_cell(
    flagged_shape: &RowFlaggedSkewShape,
    values: &mut [u32],
    content: &mut [u32],
    tableaux: &mut Vec<TableauRecord>,
    tableau_limit: Option<usize>,
    cell_index: usize,
) -> Result<(), EnumerationLimitExceeded> {
    let shape = flagged_shape.shape();
    if cell_index == shape.size() {
        if let Some(limit) = tableau_limit {
            if tableaux.len() + 1 > limit {
                return Err(EnumerationLimitExceeded { limit });
            }
        }

        tableaux.push(TableauRecord {
            values: values.to_vec(),
            content: content.to_vec(),
        });
        return Ok(());
    }

    let cell = shape.cells()[cell_index];
    let lower_bound = lower_bound(flagged_shape, values, cell);
    let upper_bound = flagged_shape.max_entry_for_row(cell.row);

    for value in lower_bound..=upper_bound {
        values[cell_index] = value;
        content[value as usize - 1] += 1;
        enumerate_tableaux_from_cell(
            flagged_shape,
            values,
            content,
            tableaux,
            tableau_limit,
            cell_index + 1,
        )?;
        content[value as usize - 1] -= 1;
        values[cell_index] = 0;
    }

    Ok(())
}

fn lower_bound(flagged_shape: &RowFlaggedSkewShape, values: &[u32], cell: Cell) -> u32 {
    let shape = flagged_shape.shape();
    let mut lower_bound = 1u32;

    if cell.col > 0 {
        let left = Cell {
            row: cell.row,
            col: cell.col - 1,
        };
        if let Some(left_index) = shape.cell_index(left) {
            lower_bound = lower_bound.max(values[left_index]);
        }
    }

    if cell.row > 0 {
        let above = Cell {
            row: cell.row - 1,
            col: cell.col,
        };
        if let Some(above_index) = shape.cell_index(above) {
            lower_bound = lower_bound.max(values[above_index] + 1);
        }
    }

    lower_bound
}

pub fn total_tableau_count(counts: &ContentStatisticCounts) -> u128 {
    counts
        .values()
        .flat_map(|stat_counts| stat_counts.values())
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shape::{RowFlaggedSkewShape, SkewShape};

    #[test]
    fn enumerates_tiny_flagged_rectangle() {
        let shape = SkewShape::from_parts(vec![2, 2], vec![]);
        let flagged = RowFlaggedSkewShape::new(shape, vec![2, 4], 4);
        let counts = enumerate_content_statistic_counts(
            &flagged,
            EnumerationOptions::new(DescentStatistic::Componentwise),
        )
        .unwrap();
        assert_eq!(total_tableau_count(&counts), 14);
    }
}
