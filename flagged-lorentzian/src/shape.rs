use std::collections::{BTreeMap, VecDeque};

use sym_poly_core::Partition;

/// A cell coordinate in English notation, indexed from zero.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Cell {
    pub row: usize,
    pub col: usize,
}

/// A skew Young diagram \(\lambda/\mu\).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SkewShape {
    outer: Partition,
    inner: Partition,
    cells: Vec<Cell>,
    index_by_cell: BTreeMap<Cell, usize>,
}

impl SkewShape {
    /// Construct a skew shape from outer and inner partitions.
    pub fn new(outer: Partition, inner: Partition) -> Self {
        assert!(
            inner.partition_less_equal(&outer),
            "inner shape must fit inside outer shape"
        );

        let cells: Vec<Cell> = outer
            .skew_diagram_boxes(&inner)
            .into_iter()
            .map(|(row, col)| Cell { row, col })
            .collect();
        let index_by_cell = cells
            .iter()
            .copied()
            .enumerate()
            .map(|(idx, cell)| (cell, idx))
            .collect();

        Self {
            outer,
            inner,
            cells,
            index_by_cell,
        }
    }

    /// Construct a skew shape from already sorted partition parts.
    pub fn from_parts(outer: Vec<u32>, inner: Vec<u32>) -> Self {
        Self::new(Partition::from_sorted(outer), Partition::from_sorted(inner))
    }

    pub fn outer(&self) -> &Partition {
        &self.outer
    }

    pub fn inner(&self) -> &Partition {
        &self.inner
    }

    pub fn cells(&self) -> &[Cell] {
        &self.cells
    }

    pub fn size(&self) -> usize {
        self.cells.len()
    }

    pub fn row_count(&self) -> usize {
        self.outer.num_parts()
    }

    pub fn cell_index(&self, cell: Cell) -> Option<usize> {
        self.index_by_cell.get(&cell).copied()
    }

    pub fn contains_cell(&self, cell: Cell) -> bool {
        self.index_by_cell.contains_key(&cell)
    }

    /// Cells in row-major order, top-to-bottom and left-to-right.
    pub fn filling_order(&self) -> Vec<usize> {
        (0..self.cells.len()).collect()
    }

    /// Cells in tableau reading order: bottom-to-top, left-to-right in each row.
    pub fn reading_order(&self) -> Vec<usize> {
        let mut order: Vec<usize> = (0..self.cells.len()).collect();
        order.sort_by_key(|&idx| {
            let cell = self.cells[idx];
            (usize::MAX - cell.row, cell.col)
        });
        order
    }

    /// Connected components of the skew diagram using edge adjacency.
    ///
    /// Components are ordered by their topmost row and then leftmost column.
    pub fn connected_components(&self) -> Vec<Vec<usize>> {
        let mut unvisited: BTreeMap<Cell, usize> = self
            .cells
            .iter()
            .copied()
            .enumerate()
            .map(|(idx, cell)| (cell, idx))
            .collect();
        let mut components = Vec::new();

        while let Some((&seed, &seed_idx)) = unvisited.iter().next() {
            unvisited.remove(&seed);
            let mut queue = VecDeque::from([seed]);
            let mut component = vec![seed_idx];

            while let Some(cell) = queue.pop_front() {
                for neighbor in neighbors(cell) {
                    if let Some(idx) = unvisited.remove(&neighbor) {
                        component.push(idx);
                        queue.push_back(neighbor);
                    }
                }
            }

            component.sort_by_key(|&idx| {
                let cell = self.cells[idx];
                (cell.row, cell.col)
            });
            components.push(component);
        }

        components.sort_by_key(|component| {
            component
                .iter()
                .map(|&idx| self.cells[idx])
                .map(|cell| (cell.row, cell.col))
                .min()
                .unwrap_or((usize::MAX, usize::MAX))
        });
        components
    }

    /// Reading orders for each connected component.
    pub fn component_reading_orders(&self) -> Vec<Vec<usize>> {
        self.connected_components()
            .into_iter()
            .map(|mut component| {
                component.sort_by_key(|&idx| {
                    let cell = self.cells[idx];
                    (usize::MAX - cell.row, cell.col)
                });
                component
            })
            .collect()
    }

    pub fn is_connected(&self) -> bool {
        self.size() == 0 || self.connected_components().len() == 1
    }
}

/// A skew shape with row upper flags and an alphabet size.
///
/// The row flag in row `r` is the largest allowed entry in that row.  Entries
/// are positive and bounded above by `alphabet`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RowFlaggedSkewShape {
    shape: SkewShape,
    row_flags: Vec<u32>,
    alphabet: usize,
}

impl RowFlaggedSkewShape {
    pub fn new(shape: SkewShape, mut row_flags: Vec<u32>, alphabet: usize) -> Self {
        assert!(alphabet > 0, "alphabet must be positive");
        while row_flags.len() < shape.row_count() {
            row_flags.push(alphabet as u32);
        }
        row_flags.truncate(shape.row_count());
        assert!(
            row_flags.iter().all(|&flag| flag > 0),
            "row flags must be positive"
        );
        assert!(
            row_flags.iter().all(|&flag| flag <= alphabet as u32),
            "row flags cannot exceed the alphabet size"
        );

        Self {
            shape,
            row_flags,
            alphabet,
        }
    }

    pub fn ordinary(shape: SkewShape, alphabet: usize) -> Self {
        let row_flags = vec![alphabet as u32; shape.row_count()];
        Self::new(shape, row_flags, alphabet)
    }

    pub fn from_parts(
        outer: Vec<u32>,
        inner: Vec<u32>,
        row_flags: Vec<u32>,
        alphabet: usize,
    ) -> Self {
        Self::new(SkewShape::from_parts(outer, inner), row_flags, alphabet)
    }

    pub fn shape(&self) -> &SkewShape {
        &self.shape
    }

    pub fn row_flags(&self) -> &[u32] {
        &self.row_flags
    }

    pub fn alphabet(&self) -> usize {
        self.alphabet
    }

    pub fn max_entry_for_row(&self, row: usize) -> u32 {
        self.row_flags
            .get(row)
            .copied()
            .unwrap_or(self.alphabet as u32)
    }
}

fn neighbors(cell: Cell) -> impl Iterator<Item = Cell> {
    let mut out = Vec::with_capacity(4);
    if cell.row > 0 {
        out.push(Cell {
            row: cell.row - 1,
            col: cell.col,
        });
    }
    out.push(Cell {
        row: cell.row + 1,
        col: cell.col,
    });
    if cell.col > 0 {
        out.push(Cell {
            row: cell.row,
            col: cell.col - 1,
        });
    }
    out.push(Cell {
        row: cell.row,
        col: cell.col + 1,
    });
    out.into_iter()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_disconnected_skew_shape() {
        let shape = SkewShape::from_parts(vec![4, 3, 1], vec![3, 1]);
        assert_eq!(shape.size(), 4);
        assert!(!shape.is_connected());
        assert_eq!(shape.connected_components().len(), 3);
    }
}
