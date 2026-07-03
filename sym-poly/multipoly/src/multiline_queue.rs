//! Multiline queues and the Ferrari--Martin labeling at `t = 0`.
//!
//! Rows are stored from bottom to top and columns are indexed from `1` to `n`,
//! following the conventions in Mandelshtam--Valencia-Porras.

use std::collections::{BTreeMap, BTreeSet};

/// A multiline queue on a ring with `n` columns.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MultilineQueue {
    n: usize,
    rows: Vec<BTreeSet<usize>>,
}

/// A pairing step in the Ferrari--Martin algorithm.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct MultilineQueuePairing {
    /// The row where the paired labeled ball starts, indexed from bottom to top.
    pub origin_row: usize,
    /// The label carried by the paired ball.
    pub label: usize,
    /// Starting column of the paired ball.
    pub from_col: usize,
    /// Target column in the row immediately below.
    pub to_col: usize,
    /// Whether the pairing wraps around the ring.
    pub wraps: bool,
}

impl MultilineQueuePairing {
    /// Indicator used in pairing displays.
    pub fn wrap_indicator(self) -> usize {
        usize::from(self.wraps)
    }

    /// Contribution to the Ferrari--Martin major index.
    pub fn major_contribution(self) -> usize {
        if self.wraps {
            self.label
                .checked_sub(self.origin_row)
                .expect("Ferrari--Martin labels are at least their origin row")
                + 1
        } else {
            0
        }
    }
}

/// The labeled multiline queue produced by the Ferrari--Martin algorithm.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FerrariMartinLabeling {
    n: usize,
    labels_by_row: Vec<BTreeMap<usize, usize>>,
    pairings: Vec<MultilineQueuePairing>,
}

impl FerrariMartinLabeling {
    /// Number of columns in the ambient ring.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Labels in each row, stored from bottom to top as maps `column -> label`.
    pub fn labels_by_row(&self) -> &[BTreeMap<usize, usize>] {
        &self.labels_by_row
    }

    /// Pairings performed by the algorithm, in algorithm order.
    pub fn pairings(&self) -> &[MultilineQueuePairing] {
        &self.pairings
    }

    /// The projected multiline-queue word in the bottom row.
    ///
    /// Empty bottom-row columns are recorded as `0`.
    pub fn projection_word(&self) -> Vec<usize> {
        if self.labels_by_row.is_empty() {
            return Vec::new();
        }
        let bottom = &self.labels_by_row[0];
        (1..=self.n)
            .map(|col| bottom.get(&col).copied().unwrap_or(0))
            .collect()
    }

    /// The Ferrari--Martin major index.
    pub fn major_index(&self) -> usize {
        self.pairings
            .iter()
            .copied()
            .map(MultilineQueuePairing::major_contribution)
            .sum()
    }

    /// Pairing data as `(origin row, label, wrap indicator)` triples.
    pub fn pairing_summary(&self) -> Vec<(usize, usize, usize)> {
        self.pairings
            .iter()
            .copied()
            .map(|pairing| (pairing.origin_row, pairing.label, pairing.wrap_indicator()))
            .collect()
    }
}

impl MultilineQueue {
    /// Create a multiline queue from row sets.
    ///
    /// Rows are stored from bottom to top.  Row sizes must be weakly decreasing,
    /// as for the row sizes `lambda'` of a partition shape.
    pub fn new(n: usize, rows: Vec<BTreeSet<usize>>) -> Self {
        assert!(n > 0, "a multiline queue must have at least one column");
        for row in &rows {
            for &col in row {
                assert!(
                    (1..=n).contains(&col),
                    "multiline-queue columns must lie in 1..=n"
                );
            }
        }
        for pair in rows.windows(2) {
            assert!(
                pair[0].len() >= pair[1].len(),
                "multiline-queue row sizes must be weakly decreasing"
            );
        }
        Self { n, rows }
    }

    /// Create a multiline queue from rows of column numbers.
    pub fn from_rows(n: usize, rows: &[&[usize]]) -> Self {
        let row_sets = rows
            .iter()
            .map(|row| {
                let row_set = row.iter().copied().collect::<BTreeSet<_>>();
                assert_eq!(
                    row_set.len(),
                    row.len(),
                    "a multiline-queue row cannot contain duplicate columns"
                );
                row_set
            })
            .collect();
        Self::new(n, row_sets)
    }

    /// Number of columns in the ambient ring.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Rows of the queue, stored from bottom to top.
    pub fn rows(&self) -> &[BTreeSet<usize>] {
        &self.rows
    }

    /// The row word, displayed with one vector per row.
    pub fn row_word_rows(&self) -> Vec<Vec<usize>> {
        self.rows
            .iter()
            .map(|row| row.iter().copied().collect())
            .collect()
    }

    /// The row word as a single flattened word.
    pub fn row_word(&self) -> Vec<usize> {
        self.row_word_rows().into_iter().flatten().collect()
    }

    /// The column word, displayed with one vector per column.
    ///
    /// Entries are row numbers, and columns are read left to right while each
    /// column is read from top to bottom.
    pub fn column_word_columns(&self) -> Vec<Vec<usize>> {
        (1..=self.n)
            .map(|col| {
                (1..=self.rows.len())
                    .rev()
                    .filter(|row| self.rows[row - 1].contains(&col))
                    .collect()
            })
            .collect()
    }

    /// The column word as a single flattened word.
    pub fn column_word(&self) -> Vec<usize> {
        self.column_word_columns().into_iter().flatten().collect()
    }

    /// The monomial weight from counting balls in each column.
    pub fn content_weight(&self) -> Vec<u32> {
        let mut weight = vec![0; self.n];
        for row in &self.rows {
            for &col in row {
                weight[col - 1] += 1;
            }
        }
        weight
    }

    /// Run the Ferrari--Martin labeling algorithm.
    pub fn ferrari_martin_labeling(&self) -> FerrariMartinLabeling {
        let row_count = self.rows.len();
        let mut labels_by_row = vec![BTreeMap::new(); row_count];
        let mut pairings = Vec::new();

        for row in (2..=row_count).rev() {
            let row_idx = row - 1;
            for &col in &self.rows[row_idx] {
                labels_by_row[row_idx].entry(col).or_insert(row);
            }

            let mut labeled_balls = labels_by_row[row_idx]
                .iter()
                .map(|(&col, &label)| (label, col))
                .collect::<Vec<_>>();
            labeled_balls.sort_by(|(label_a, col_a), (label_b, col_b)| {
                label_b.cmp(label_a).then_with(|| col_a.cmp(col_b))
            });

            for (label, from_col) in labeled_balls {
                let below_idx = row_idx - 1;
                let to_col = (0..self.n)
                    .map(|offset| ((from_col - 1 + offset) % self.n) + 1)
                    .find(|candidate| {
                        self.rows[below_idx].contains(candidate)
                            && !labels_by_row[below_idx].contains_key(candidate)
                    })
                    .expect("a lower-row target must exist for every labeled ball");

                labels_by_row[below_idx].insert(to_col, label);
                pairings.push(MultilineQueuePairing {
                    origin_row: row,
                    label,
                    from_col,
                    to_col,
                    wraps: to_col < from_col,
                });
            }
        }

        if let Some(bottom_row) = self.rows.first() {
            for &col in bottom_row {
                labels_by_row[0].entry(col).or_insert(1);
            }
        }

        FerrariMartinLabeling {
            n: self.n,
            labels_by_row,
            pairings,
        }
    }

    /// The pair `(column weight, Ferrari--Martin major index)`.
    pub fn weight_data(&self) -> (Vec<u32>, usize) {
        (
            self.content_weight(),
            self.ferrari_martin_labeling().major_index(),
        )
    }
}

/// Generate all multiline queues with the prescribed row sizes.
///
/// The row sizes are listed from bottom to top.  The result is capped to avoid
/// accidental generation of a very large product of binomial coefficients.
pub fn multiline_queues_with_row_sizes(
    n: usize,
    row_sizes: &[usize],
    max_queues: usize,
) -> Result<Vec<MultilineQueue>, String> {
    if n == 0 {
        return Err("a multiline queue must have at least one column".to_string());
    }
    if max_queues == 0 {
        return Err("max_queues must be positive".to_string());
    }
    if row_sizes.iter().any(|&size| size > n) {
        return Err("row sizes must be at most the number of columns".to_string());
    }
    if row_sizes.windows(2).any(|pair| pair[0] < pair[1]) {
        return Err("row sizes must be weakly decreasing".to_string());
    }

    let row_choices = row_sizes
        .iter()
        .map(|&size| subsets_of_size(n, size, max_queues))
        .collect::<Result<Vec<_>, _>>()?;
    let mut queues = Vec::new();
    let mut rows = Vec::new();
    generate_queue_rows(n, &row_choices, 0, &mut rows, &mut queues, max_queues)?;
    Ok(queues)
}

/// Count multiline queues by `(column weight, major index)`.
pub fn multiline_queue_weight_counts(
    n: usize,
    row_sizes: &[usize],
    max_queues: usize,
) -> Result<BTreeMap<(Vec<u32>, usize), usize>, String> {
    let mut counts = BTreeMap::new();
    for queue in multiline_queues_with_row_sizes(n, row_sizes, max_queues)? {
        *counts.entry(queue.weight_data()).or_insert(0) += 1;
    }
    Ok(counts)
}

fn generate_queue_rows(
    n: usize,
    row_choices: &[Vec<BTreeSet<usize>>],
    row_idx: usize,
    rows: &mut Vec<BTreeSet<usize>>,
    queues: &mut Vec<MultilineQueue>,
    max_queues: usize,
) -> Result<(), String> {
    if row_idx == row_choices.len() {
        queues.push(MultilineQueue::new(n, rows.clone()));
        if queues.len() > max_queues {
            return Err(format!(
                "multiline queue cap exceeded: more than {max_queues} queues"
            ));
        }
        return Ok(());
    }

    for row in &row_choices[row_idx] {
        rows.push(row.clone());
        generate_queue_rows(n, row_choices, row_idx + 1, rows, queues, max_queues)?;
        rows.pop();
    }
    Ok(())
}

fn subsets_of_size(
    n: usize,
    size: usize,
    max_subsets: usize,
) -> Result<Vec<BTreeSet<usize>>, String> {
    let mut subsets = Vec::new();
    let mut current = Vec::new();
    build_subsets(n, size, 1, &mut current, &mut subsets, max_subsets)?;
    Ok(subsets)
}

fn build_subsets(
    n: usize,
    size: usize,
    start_col: usize,
    current: &mut Vec<usize>,
    subsets: &mut Vec<BTreeSet<usize>>,
    max_subsets: usize,
) -> Result<(), String> {
    if current.len() == size {
        subsets.push(current.iter().copied().collect());
        if subsets.len() > max_subsets {
            return Err(format!(
                "row subset cap exceeded: more than {max_subsets} subsets"
            ));
        }
        return Ok(());
    }

    let remaining = size - current.len();
    let last_start = n + 1 - remaining;
    for col in start_col..=last_start {
        current.push(col);
        build_subsets(n, size, col + 1, current, subsets, max_subsets)?;
        current.pop();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn source_example_queue() -> MultilineQueue {
        MultilineQueue::from_rows(6, &[&[1, 2, 3, 4], &[1, 3, 5, 6], &[2, 3], &[3, 5]])
    }

    #[test]
    fn row_and_column_words_match_source_example() {
        let queue = source_example_queue();
        assert_eq!(
            queue.row_word_rows(),
            vec![vec![1, 2, 3, 4], vec![1, 3, 5, 6], vec![2, 3], vec![3, 5]]
        );
        assert_eq!(
            queue.column_word_columns(),
            vec![
                vec![2, 1],
                vec![3, 1],
                vec![4, 3, 2, 1],
                vec![1],
                vec![4, 2],
                vec![2]
            ]
        );
    }

    #[test]
    fn ferrari_martin_weight_matches_source_example() {
        let queue = source_example_queue();
        let labeling = queue.ferrari_martin_labeling();
        assert_eq!(
            labeling.pairing_summary(),
            vec![
                (4, 4, 0),
                (4, 4, 1),
                (3, 4, 0),
                (3, 4, 0),
                (2, 4, 0),
                (2, 4, 1),
                (2, 2, 0),
                (2, 2, 1)
            ]
        );
        assert_eq!(labeling.major_index(), 5);
        assert_eq!(labeling.projection_word(), vec![4, 2, 4, 2, 0, 0]);
        assert_eq!(queue.content_weight(), vec![2, 2, 4, 1, 2, 1]);
        assert_eq!(queue.weight_data(), (vec![2, 2, 4, 1, 2, 1], 5));
    }

    #[test]
    fn generate_small_multiline_queues() {
        let queues = multiline_queues_with_row_sizes(3, &[2, 1], 20).unwrap();
        assert_eq!(queues.len(), 9);

        let counts = multiline_queue_weight_counts(3, &[2, 1], 20).unwrap();
        assert_eq!(counts.values().sum::<usize>(), 9);
        assert_eq!(counts.get(&(vec![2, 1, 0], 0)), Some(&1));
        assert_eq!(counts.get(&(vec![1, 1, 1], 1)), Some(&1));
    }
}
