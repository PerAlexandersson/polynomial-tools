use std::fmt;

use crate::shape::SkewShape;

/// Descent data used for fiber refinements.
///
/// A value is a tuple of bitmasks, one per word.  A bit at position `p-1`
/// records a descent at one-indexed position `p`.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DescentData(Vec<u64>);

impl DescentData {
    pub fn new(masks: Vec<u64>) -> Self {
        Self(masks)
    }

    pub fn masks(&self) -> &[u64] {
        &self.0
    }

    pub fn major_index(&self) -> u32 {
        self.0
            .iter()
            .map(|&mask| {
                (0..64)
                    .filter(|bit| (mask & (1u64 << bit)) != 0)
                    .map(|bit| bit as u32 + 1)
                    .sum::<u32>()
            })
            .sum()
    }

    pub fn descent_count(&self) -> u32 {
        self.0.iter().map(|mask| mask.count_ones()).sum()
    }
}

impl fmt::Display for DescentData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let parts: Vec<String> = self
            .0
            .iter()
            .map(|&mask| format_descent_mask(mask))
            .collect();
        write!(f, "({})", parts.join(","))
    }
}

/// Which descent statistic to record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DescentStatistic {
    /// One row-reading word for the whole skew shape.
    Global,
    /// One row-reading word for each connected component of the skew shape.
    Componentwise,
}

impl DescentStatistic {
    pub fn reading_orders(self, shape: &SkewShape) -> Vec<Vec<usize>> {
        match self {
            DescentStatistic::Global => vec![shape.reading_order()],
            DescentStatistic::Componentwise => shape.component_reading_orders(),
        }
    }
}

pub fn descent_data_for_values(values: &[u32], reading_orders: &[Vec<usize>]) -> DescentData {
    DescentData::new(
        reading_orders
            .iter()
            .map(|order| descent_mask_for_order(values, order))
            .collect(),
    )
}

pub fn active_subword_descent_data_for_values(
    values: &[u32],
    reading_orders: &[Vec<usize>],
    lower_label: u32,
) -> DescentData {
    let upper_label = lower_label + 1;
    DescentData::new(
        reading_orders
            .iter()
            .map(|order| {
                active_subword_descent_mask_for_order(values, order, lower_label, upper_label)
            })
            .collect(),
    )
}

pub fn descent_mask_for_order(values: &[u32], order: &[usize]) -> u64 {
    assert!(
        order.len() <= 65,
        "descent masks support words of length at most 65"
    );
    let mut mask = 0u64;
    for position in 0..order.len().saturating_sub(1) {
        if values[order[position]] > values[order[position + 1]] {
            mask |= 1u64 << position;
        }
    }
    mask
}

pub fn active_subword_descent_mask_for_order(
    values: &[u32],
    order: &[usize],
    lower_label: u32,
    upper_label: u32,
) -> u64 {
    let mut mask = 0u64;
    let mut previous = None;
    let mut active_position = 0usize;

    for &idx in order {
        let value = values[idx];
        if value != lower_label && value != upper_label {
            continue;
        }

        assert!(
            active_position < 65,
            "active-subword descent masks support active words of length at most 65"
        );

        if let Some(previous_value) = previous {
            if previous_value > value {
                mask |= 1u64 << (active_position - 1);
            }
        }
        previous = Some(value);
        active_position += 1;
    }

    mask
}

pub fn format_descent_mask(mask: u64) -> String {
    let positions: Vec<String> = (0..64)
        .filter(|bit| (mask & (1u64 << bit)) != 0)
        .map(|bit| (bit + 1).to_string())
        .collect();
    format!("{{{}}}", positions.join(","))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn descent_mask_uses_one_indexed_positions() {
        let values = vec![2, 1, 2, 1];
        let order = vec![0, 1, 2, 3];
        let mask = descent_mask_for_order(&values, &order);
        assert_eq!(format_descent_mask(mask), "{1,3}");
    }

    #[test]
    fn active_subword_mask_ignores_inactive_letters() {
        let values = vec![4, 2, 5, 1, 4, 5, 3, 4];
        let order = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let mask = active_subword_descent_mask_for_order(&values, &order, 4, 5);

        // The active subword is 4,5,4,5,4, with descents at filtered
        // positions 2 and 4.
        assert_eq!(format_descent_mask(mask), "{2,4}");
    }
}
