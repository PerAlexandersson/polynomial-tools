use std::collections::{BTreeSet, VecDeque};

use crate::shape::{Cell, SkewShape};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrystalDirection {
    E,
    F,
}

/// Apply the type-A crystal operator \(e_i\) to a tableau.
///
/// The parameter `lower_label` is the one-indexed label \(i\).  The operator
/// acts on the pair of labels \(i,i+1\).
pub fn crystal_e(shape: &SkewShape, values: &[u32], lower_label: u32) -> Option<Vec<u32>> {
    apply_crystal_operator(shape, values, lower_label, CrystalDirection::E)
}

/// Apply the type-A crystal operator \(f_i\) to a tableau.
///
/// The parameter `lower_label` is the one-indexed label \(i\).  The operator
/// acts on the pair of labels \(i,i+1\).
pub fn crystal_f(shape: &SkewShape, values: &[u32], lower_label: u32) -> Option<Vec<u32>> {
    apply_crystal_operator(shape, values, lower_label, CrystalDirection::F)
}

pub fn apply_crystal_operator(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
    direction: CrystalDirection,
) -> Option<Vec<u32>> {
    let upper_label = lower_label + 1;
    let mut unpaired_lower = Vec::new();
    let mut unpaired_upper = Vec::new();

    for cell_index in shape.reading_order() {
        match values[cell_index] {
            value if value == lower_label => {
                if unpaired_upper.pop().is_none() {
                    unpaired_lower.push(cell_index);
                }
            }
            value if value == upper_label => unpaired_upper.push(cell_index),
            _ => {}
        }
    }

    let mut next_values = values.to_vec();
    match direction {
        CrystalDirection::E => {
            let cell_index = *unpaired_upper.first()?;
            next_values[cell_index] = lower_label;
        }
        CrystalDirection::F => {
            let cell_index = unpaired_lower.pop()?;
            next_values[cell_index] = upper_label;
        }
    }

    Some(next_values)
}

pub fn active_component_crystal_e_images(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
) -> Vec<Vec<u32>> {
    active_component_crystal_images(shape, values, lower_label, CrystalDirection::E)
}

pub fn active_component_crystal_f_images(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
) -> Vec<Vec<u32>> {
    active_component_crystal_images(shape, values, lower_label, CrystalDirection::F)
}

pub fn bender_knuth_e_images(shape: &SkewShape, values: &[u32], lower_label: u32) -> Vec<Vec<u32>> {
    bender_knuth_unit_images(shape, values, lower_label, CrystalDirection::E)
}

pub fn bender_knuth_f_images(shape: &SkewShape, values: &[u32], lower_label: u32) -> Vec<Vec<u32>> {
    bender_knuth_unit_images(shape, values, lower_label, CrystalDirection::F)
}

pub fn bender_knuth_unit_images(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
    direction: CrystalDirection,
) -> Vec<Vec<u32>> {
    let upper_label = lower_label + 1;
    let mut images = BTreeSet::new();

    for row in 0..shape.row_count() {
        let mut row_cells: Vec<_> = shape
            .cells()
            .iter()
            .enumerate()
            .filter(|(_, cell)| cell.row == row)
            .map(|(idx, cell)| (idx, cell.col))
            .collect();
        row_cells.sort_by_key(|(_, col)| *col);

        let candidate = match direction {
            CrystalDirection::E => row_cells
                .iter()
                .filter(|(idx, _)| {
                    values[*idx] == upper_label
                        && !is_vertically_paired_upper(shape, values, *idx, lower_label)
                })
                .map(|(idx, _)| *idx)
                .next(),
            CrystalDirection::F => row_cells
                .iter()
                .rev()
                .filter(|(idx, _)| {
                    values[*idx] == lower_label
                        && !is_vertically_paired_lower(shape, values, *idx, lower_label)
                })
                .map(|(idx, _)| *idx)
                .next(),
        };

        if let Some(cell_index) = candidate {
            let mut image = values.to_vec();
            image[cell_index] = match direction {
                CrystalDirection::E => lower_label,
                CrystalDirection::F => upper_label,
            };
            if is_semistandard(shape, &image) {
                images.insert(image);
            }
        }
    }

    images.into_iter().collect()
}

pub fn active_component_crystal_images(
    shape: &SkewShape,
    values: &[u32],
    lower_label: u32,
    direction: CrystalDirection,
) -> Vec<Vec<u32>> {
    let global_reading_order = shape.reading_order();
    let mut images = BTreeSet::new();

    for component in active_components(shape, values, lower_label) {
        let component_set: BTreeSet<_> = component.into_iter().collect();
        let component_order: Vec<_> = global_reading_order
            .iter()
            .copied()
            .filter(|idx| component_set.contains(idx))
            .collect();
        if let Some(image) =
            apply_crystal_operator_on_order(values, &component_order, lower_label, direction)
        {
            if is_semistandard(shape, &image) {
                images.insert(image);
            }
        }
    }

    images.into_iter().collect()
}

fn is_vertically_paired_lower(
    shape: &SkewShape,
    values: &[u32],
    cell_index: usize,
    lower_label: u32,
) -> bool {
    let cell = shape.cells()[cell_index];
    shape
        .cell_index(Cell {
            row: cell.row + 1,
            col: cell.col,
        })
        .is_some_and(|below_idx| values[below_idx] == lower_label + 1)
}

fn is_vertically_paired_upper(
    shape: &SkewShape,
    values: &[u32],
    cell_index: usize,
    lower_label: u32,
) -> bool {
    let cell = shape.cells()[cell_index];
    cell.row > 0
        && shape
            .cell_index(Cell {
                row: cell.row - 1,
                col: cell.col,
            })
            .is_some_and(|above_idx| values[above_idx] == lower_label)
}

pub fn apply_crystal_operator_on_order(
    values: &[u32],
    order: &[usize],
    lower_label: u32,
    direction: CrystalDirection,
) -> Option<Vec<u32>> {
    let upper_label = lower_label + 1;
    let mut unpaired_lower = Vec::new();
    let mut unpaired_upper = Vec::new();

    for &cell_index in order {
        match values[cell_index] {
            value if value == lower_label => {
                if unpaired_upper.pop().is_none() {
                    unpaired_lower.push(cell_index);
                }
            }
            value if value == upper_label => unpaired_upper.push(cell_index),
            _ => {}
        }
    }

    let mut next_values = values.to_vec();
    match direction {
        CrystalDirection::E => {
            let cell_index = *unpaired_upper.first()?;
            next_values[cell_index] = lower_label;
        }
        CrystalDirection::F => {
            let cell_index = unpaired_lower.pop()?;
            next_values[cell_index] = upper_label;
        }
    }

    Some(next_values)
}

pub fn active_components(shape: &SkewShape, values: &[u32], lower_label: u32) -> Vec<Vec<usize>> {
    let upper_label = lower_label + 1;
    let mut active = BTreeSet::new();
    for (idx, &value) in values.iter().enumerate() {
        if value == lower_label || value == upper_label {
            active.insert(idx);
        }
    }

    let mut components = Vec::new();
    while let Some(&seed) = active.iter().next() {
        active.remove(&seed);
        let mut component = vec![seed];
        let mut queue = VecDeque::from([seed]);

        while let Some(cell_index) = queue.pop_front() {
            for neighbor_index in neighbor_indices(shape, cell_index) {
                if active.remove(&neighbor_index) {
                    component.push(neighbor_index);
                    queue.push_back(neighbor_index);
                }
            }
        }

        component.sort_by_key(|&idx| {
            let cell = shape.cells()[idx];
            (cell.row, cell.col)
        });
        components.push(component);
    }

    components.sort_by_key(|component| {
        component
            .iter()
            .map(|&idx| shape.cells()[idx])
            .map(|cell| (cell.row, cell.col))
            .min()
            .unwrap_or((usize::MAX, usize::MAX))
    });
    components
}

pub fn is_semistandard(shape: &SkewShape, values: &[u32]) -> bool {
    shape.cells().iter().enumerate().all(|(idx, &cell)| {
        let value = values[idx];
        let right_ok = shape
            .cell_index(Cell {
                row: cell.row,
                col: cell.col + 1,
            })
            .is_none_or(|right_idx| value <= values[right_idx]);
        let below_ok = shape
            .cell_index(Cell {
                row: cell.row + 1,
                col: cell.col,
            })
            .is_none_or(|below_idx| value < values[below_idx]);
        right_ok && below_ok
    })
}

fn neighbor_indices(shape: &SkewShape, cell_index: usize) -> Vec<usize> {
    let cell = shape.cells()[cell_index];
    let mut neighbors = Vec::with_capacity(4);
    if cell.row > 0 {
        if let Some(idx) = shape.cell_index(Cell {
            row: cell.row - 1,
            col: cell.col,
        }) {
            neighbors.push(idx);
        }
    }
    if let Some(idx) = shape.cell_index(Cell {
        row: cell.row + 1,
        col: cell.col,
    }) {
        neighbors.push(idx);
    }
    if cell.col > 0 {
        if let Some(idx) = shape.cell_index(Cell {
            row: cell.row,
            col: cell.col - 1,
        }) {
            neighbors.push(idx);
        }
    }
    if let Some(idx) = shape.cell_index(Cell {
        row: cell.row,
        col: cell.col + 1,
    }) {
        neighbors.push(idx);
    }
    neighbors
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f_moves_single_row_along_string() {
        let shape = SkewShape::from_parts(vec![2], vec![]);
        assert_eq!(crystal_f(&shape, &[1, 1], 1), Some(vec![1, 2]));
        assert_eq!(crystal_f(&shape, &[1, 2], 1), Some(vec![2, 2]));
        assert_eq!(crystal_f(&shape, &[2, 2], 1), None);
    }

    #[test]
    fn e_moves_single_row_along_string() {
        let shape = SkewShape::from_parts(vec![2], vec![]);
        assert_eq!(crystal_e(&shape, &[2, 2], 1), Some(vec![1, 2]));
        assert_eq!(crystal_e(&shape, &[1, 2], 1), Some(vec![1, 1]));
        assert_eq!(crystal_e(&shape, &[1, 1], 1), None);
    }

    #[test]
    fn both_operators_are_undefined_on_frozen_column() {
        let shape = SkewShape::from_parts(vec![1, 1], vec![]);
        assert_eq!(crystal_f(&shape, &[1, 2], 1), None);
        assert_eq!(crystal_e(&shape, &[1, 2], 1), None);
    }

    #[test]
    fn active_components_ignore_frozen_gap() {
        let shape = SkewShape::from_parts(vec![4, 2], vec![2]);
        let values = vec![1, 1, 2, 2];
        let components = active_components(&shape, &values, 1);
        assert_eq!(components.len(), 2);
        assert_eq!(crystal_f(&shape, &values, 1), None);
        assert!(active_component_crystal_f_images(&shape, &values, 1).contains(&vec![1, 2, 2, 2]));
    }

    #[test]
    fn semistandard_check_detects_column_equality() {
        let shape = SkewShape::from_parts(vec![1, 1], vec![]);
        assert!(is_semistandard(&shape, &[1, 2]));
        assert!(!is_semistandard(&shape, &[1, 1]));
    }

    #[test]
    fn bender_knuth_free_cell_move_sees_connected_obstruction() {
        let shape = SkewShape::from_parts(vec![4, 2], vec![1]);
        let values = vec![1, 1, 5, 2, 2];
        assert_eq!(crystal_f(&shape, &values, 1), None);
        assert!(bender_knuth_f_images(&shape, &values, 1).contains(&vec![1, 2, 5, 2, 2]));
    }
}
