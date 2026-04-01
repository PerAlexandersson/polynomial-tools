use std::collections::BTreeSet;
use std::fmt;

use crate::{Composition, Partition};

/// A left-justified Young tableau with positive integer entries.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Tableau {
    rows: Vec<Vec<u32>>,
}

impl Tableau {
    /// Create a tableau from rows.
    ///
    /// Rows must have weakly decreasing lengths and use positive entries.
    pub fn new(mut rows: Vec<Vec<u32>>) -> Self {
        while rows.last().is_some_and(|row| row.is_empty()) {
            rows.pop();
        }

        let mut prev_len = usize::MAX;
        for row in &rows {
            assert!(
                row.iter().all(|&x| x > 0),
                "tableau entries must be positive"
            );
            assert!(
                row.len() <= prev_len,
                "tableau rows must have weakly decreasing lengths"
            );
            prev_len = row.len();
        }

        Tableau { rows }
    }

    /// The empty tableau.
    pub fn empty() -> Self {
        Tableau { rows: Vec::new() }
    }

    /// Access the tableau rows.
    pub fn rows(&self) -> &[Vec<u32>] {
        &self.rows
    }

    /// Return the shape of the tableau.
    pub fn shape(&self) -> Partition {
        Partition::from_sorted(self.rows.iter().map(|row| row.len() as u32).collect())
    }

    /// Number of boxes.
    pub fn size(&self) -> u32 {
        self.rows.iter().map(|row| row.len() as u32).sum()
    }

    /// Weight composition, where the i-th part counts entries equal to i + 1.
    pub fn weight(&self) -> Composition {
        let max_entry = self
            .rows
            .iter()
            .flat_map(|row| row.iter())
            .copied()
            .max()
            .unwrap_or(0) as usize;

        let mut counts = vec![0u32; max_entry];
        for row in &self.rows {
            for &entry in row {
                counts[entry as usize - 1] += 1;
            }
        }
        Composition::new(counts)
    }

    /// Row reading word, read left-to-right in each row from top to bottom.
    pub fn row_reading_word(&self) -> Vec<u32> {
        self.rows
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect()
    }

    /// Reverse row reading word, read left-to-right in each row from bottom to top.
    pub fn reverse_row_reading_word(&self) -> Vec<u32> {
        self.rows
            .iter()
            .rev()
            .flat_map(|row| row.iter().copied())
            .collect()
    }

    /// Whether rows are weakly increasing.
    pub fn is_row_weakly_increasing(&self) -> bool {
        self.rows
            .iter()
            .all(|row| row.windows(2).all(|pair| pair[0] <= pair[1]))
    }

    /// Whether columns are strictly increasing.
    pub fn is_column_strictly_increasing(&self) -> bool {
        for r in 1..self.rows.len() {
            for c in 0..self.rows[r].len() {
                if self.rows[r - 1].len() > c && self.rows[r - 1][c] >= self.rows[r][c] {
                    return false;
                }
            }
        }
        true
    }

    /// Whether this is a semistandard Young tableau.
    pub fn is_semistandard(&self) -> bool {
        self.is_row_weakly_increasing() && self.is_column_strictly_increasing()
    }

    /// Whether this is a standard Young tableau.
    pub fn is_standard(&self) -> bool {
        if !self.is_semistandard() {
            return false;
        }

        let n = self.size() as usize;
        let mut seen = vec![false; n + 1];
        for row in &self.rows {
            for &entry in row {
                let idx = entry as usize;
                if idx == 0 || idx > n || seen[idx] {
                    return false;
                }
                seen[idx] = true;
            }
        }
        seen.into_iter().skip(1).all(|b| b)
    }

    /// Descent set of a standard tableau.
    pub fn descent_set(&self) -> BTreeSet<u32> {
        assert!(
            self.is_standard(),
            "descent_set requires a standard tableau"
        );

        let mut positions = vec![(0usize, 0usize); self.size() as usize + 1];
        for (r, row) in self.rows.iter().enumerate() {
            for (c, &entry) in row.iter().enumerate() {
                positions[entry as usize] = (r, c);
            }
        }

        let mut descents = BTreeSet::new();
        for i in 1..self.size() {
            if positions[(i + 1) as usize].0 > positions[i as usize].0 {
                descents.insert(i);
            }
        }
        descents
    }

    /// Enumerate all standard Young tableaux of the given shape.
    pub fn standard_tableaux(shape: &Partition) -> Vec<Tableau> {
        standard_tableaux_from_shape(shape.parts())
    }

    /// Enumerate semistandard tableaux of the given shape with entries in
    /// `{1, ..., max_entry}`.
    pub fn semistandard_tableaux(shape: &Partition, max_entry: u32) -> Vec<Tableau> {
        if max_entry == 0 {
            return if shape.is_empty() {
                vec![Tableau::empty()]
            } else {
                Vec::new()
            };
        }

        let mut grid: Vec<Vec<u32>> = shape
            .parts()
            .iter()
            .map(|&part| vec![0; part as usize])
            .collect();
        let cells: Vec<(usize, usize)> = shape
            .parts()
            .iter()
            .enumerate()
            .flat_map(|(r, &part)| (0..part as usize).map(move |c| (r, c)))
            .collect();
        let mut out = Vec::new();
        semistandard_helper(&mut grid, &cells, 0, max_entry, &mut out);
        out.into_iter().map(Tableau::new).collect()
    }
}

fn standard_tableaux_from_shape(shape: &[u32]) -> Vec<Tableau> {
    let size: u32 = shape.iter().sum();
    if size == 0 {
        return vec![Tableau::empty()];
    }

    let mut result = Vec::new();
    for (row, _) in corner_cells(shape) {
        let mut smaller = shape.to_vec();
        smaller[row] -= 1;
        while smaller.last() == Some(&0) {
            smaller.pop();
        }

        for tableau in standard_tableaux_from_shape(&smaller) {
            let mut rows = tableau.rows().to_vec();
            while rows.len() <= row {
                rows.push(Vec::new());
            }
            rows[row].push(size);
            result.push(Tableau::new(rows));
        }
    }
    result
}

fn corner_cells(shape: &[u32]) -> Vec<(usize, usize)> {
    let mut corners = Vec::new();
    for row in 0..shape.len() {
        let next = if row + 1 < shape.len() {
            shape[row + 1]
        } else {
            0
        };
        if shape[row] > next {
            corners.push((row, shape[row] as usize - 1));
        }
    }
    corners
}

fn semistandard_helper(
    grid: &mut [Vec<u32>],
    cells: &[(usize, usize)],
    idx: usize,
    max_entry: u32,
    out: &mut Vec<Vec<Vec<u32>>>,
) {
    if idx == cells.len() {
        out.push(grid.to_vec());
        return;
    }

    let (row, col) = cells[idx];
    let mut min_value = 1u32;
    if col > 0 {
        min_value = min_value.max(grid[row][col - 1]);
    }
    if row > 0 && grid[row - 1].len() > col {
        min_value = min_value.max(grid[row - 1][col] + 1);
    }

    for value in min_value..=max_entry {
        grid[row][col] = value;
        semistandard_helper(grid, cells, idx + 1, max_entry, out);
    }
}

impl fmt::Display for Tableau {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, row) in self.rows.iter().enumerate() {
            if i > 0 {
                writeln!(f)?;
            }
            for (j, entry) in row.iter().enumerate() {
                if j > 0 {
                    write!(f, " ")?;
                }
                write!(f, "{entry}")?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tableau_shape_weight_and_words() {
        let tab = Tableau::new(vec![vec![1, 1, 3], vec![2, 4]]);
        assert_eq!(tab.shape(), Partition::from_sorted(vec![3, 2]));
        assert_eq!(tab.weight(), Composition::new(vec![2, 1, 1, 1]));
        assert_eq!(tab.row_reading_word(), vec![1, 1, 3, 2, 4]);
        assert_eq!(tab.reverse_row_reading_word(), vec![2, 4, 1, 1, 3]);
    }

    #[test]
    fn test_semistandard_and_standard_checks() {
        let ssyt = Tableau::new(vec![vec![1, 1, 3], vec![2, 4]]);
        assert!(ssyt.is_semistandard());
        assert!(!ssyt.is_standard());

        let syt = Tableau::new(vec![vec![1, 2], vec![3]]);
        assert!(syt.is_standard());
        assert_eq!(syt.descent_set(), BTreeSet::from([2]));
    }

    #[test]
    fn test_standard_tableaux_shape_21() {
        let tabs = Tableau::standard_tableaux(&Partition::from_sorted(vec![2, 1]));
        assert_eq!(tabs.len(), 2);
        assert!(tabs.iter().all(Tableau::is_standard));
    }

    #[test]
    fn test_semistandard_tableaux_shape_21_entries_up_to_2() {
        let tabs = Tableau::semistandard_tableaux(&Partition::from_sorted(vec![2, 1]), 2);
        assert_eq!(tabs.len(), 2);
        assert!(tabs.iter().all(Tableau::is_semistandard));
    }
}
