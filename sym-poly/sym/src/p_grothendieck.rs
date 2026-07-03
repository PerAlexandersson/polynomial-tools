//! Shifted multiset tableaux for weak symmetric P-Grothendieck examples.
//!
//! This implements the small tableau model from Hawkes's definition of
//! shifted multiset tableaux.  It is intended for exact small examples rather
//! than large-scale enumeration.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt;

use sym_poly_core::Partition;

/// An entry in the ordered alphabet `1' < 1 < 2' < 2 < ...`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ShiftedMultisetEntry {
    value: u32,
    primed: bool,
}

impl ShiftedMultisetEntry {
    pub fn primed(value: u32) -> Self {
        assert!(value > 0, "entry value must be positive");
        Self {
            value,
            primed: true,
        }
    }

    pub fn unprimed(value: u32) -> Self {
        assert!(value > 0, "entry value must be positive");
        Self {
            value,
            primed: false,
        }
    }

    pub fn value(&self) -> u32 {
        self.value
    }

    pub fn is_primed(&self) -> bool {
        self.primed
    }

    fn order_key(&self) -> (u32, u8) {
        (self.value, if self.primed { 0 } else { 1 })
    }

    fn less_u(&self, other: &Self) -> bool {
        self < other || (self.value == other.value && !self.primed && !other.primed)
    }

    fn less_p(&self, other: &Self) -> bool {
        self < other || (self.value == other.value && self.primed && other.primed)
    }
}

impl Ord for ShiftedMultisetEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.order_key().cmp(&other.order_key())
    }
}

impl PartialOrd for ShiftedMultisetEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for ShiftedMultisetEntry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.primed {
            write!(f, "{}'", self.value)
        } else {
            write!(f, "{}", self.value)
        }
    }
}

/// A shifted multiset tableau of strict shape.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ShiftedMultisetTableau {
    shape: Partition,
    boxes: BTreeMap<(usize, usize), Vec<ShiftedMultisetEntry>>,
}

impl ShiftedMultisetTableau {
    /// Enumerate tableaux in Hawkes's `SMT_0(shape)` with entries bounded by
    /// `max_entry` and with exactly `total_entries` entries.
    pub fn enumerate_smt0(shape: &Partition, max_entry: u32, total_entries: usize) -> Vec<Self> {
        assert!(is_strict_partition(shape), "shape must be strict");
        assert!(max_entry > 0, "max_entry must be positive");

        let cells = shifted_shape_cells(shape);
        assert!(
            total_entries >= cells.len(),
            "total_entries must be at least the number of boxes"
        );

        let max_box_size = total_entries - cells.len() + 1;
        let candidate_boxes = shifted_multiset_boxes(max_entry, max_box_size);
        let mut assigned = BTreeMap::new();
        let mut results = Vec::new();
        enumerate_tableaux_rec(
            shape,
            &cells,
            &candidate_boxes,
            total_entries,
            0,
            0,
            &mut assigned,
            &mut results,
        );
        results
    }

    pub fn shape(&self) -> &Partition {
        &self.shape
    }

    pub fn boxes(&self) -> &BTreeMap<(usize, usize), Vec<ShiftedMultisetEntry>> {
        &self.boxes
    }

    /// The vector `(w_1, ..., w_max_entry)` counting occurrences of each value.
    pub fn weight(&self, max_entry: u32) -> Vec<u32> {
        let mut weight = vec![0; max_entry as usize];
        for entry in self.boxes.values().flatten() {
            if entry.value <= max_entry {
                weight[entry.value as usize - 1] += 1;
            }
        }
        weight
    }

    /// Hawkes's diagonal weight `(T_1, ..., T_l)`, where `l` is the longest part.
    pub fn diagonal_weight(&self) -> Vec<u32> {
        let longest = self.shape.part(0) as usize;
        let mut entries_by_diagonal = vec![0u32; longest];
        let mut boxes_by_diagonal = vec![0u32; longest];
        for (&(row, col), entries) in &self.boxes {
            let offset = col
                .checked_sub(row)
                .expect("shifted cell columns should be at least their row index");
            let label = longest - offset;
            let index = label - 1;
            entries_by_diagonal[index] += entries.len() as u32;
            boxes_by_diagonal[index] += 1;
        }
        entries_by_diagonal
            .into_iter()
            .zip(boxes_by_diagonal)
            .map(|(entries, boxes)| entries - boxes)
            .collect()
    }
}

/// Distribution of tableaux by `(x-weight, diagonal t-weight)`.
pub fn shifted_multiset_tableau_distribution(
    shape: &Partition,
    max_entry: u32,
    total_entries: usize,
) -> BTreeMap<(Vec<u32>, Vec<u32>), usize> {
    let mut distribution = BTreeMap::new();
    for tableau in ShiftedMultisetTableau::enumerate_smt0(shape, max_entry, total_entries) {
        let key = (tableau.weight(max_entry), tableau.diagonal_weight());
        *distribution.entry(key).or_insert(0) += 1;
    }
    distribution
}

fn enumerate_tableaux_rec(
    shape: &Partition,
    cells: &[(usize, usize)],
    candidate_boxes: &[Vec<ShiftedMultisetEntry>],
    total_entries: usize,
    index: usize,
    entries_so_far: usize,
    assigned: &mut BTreeMap<(usize, usize), Vec<ShiftedMultisetEntry>>,
    results: &mut Vec<ShiftedMultisetTableau>,
) {
    if index == cells.len() {
        if entries_so_far == total_entries && rows_have_unprimed_minima(shape, assigned) {
            results.push(ShiftedMultisetTableau {
                shape: shape.clone(),
                boxes: assigned.clone(),
            });
        }
        return;
    }

    let remaining_after = cells.len() - index - 1;
    let cell = cells[index];
    for box_entries in candidate_boxes {
        let next_count = entries_so_far + box_entries.len();
        if next_count + remaining_after > total_entries {
            continue;
        }
        if next_count > total_entries {
            continue;
        }
        if !is_compatible_with_assigned(cell, box_entries, assigned) {
            continue;
        }

        assigned.insert(cell, box_entries.clone());
        enumerate_tableaux_rec(
            shape,
            cells,
            candidate_boxes,
            total_entries,
            index + 1,
            next_count,
            assigned,
            results,
        );
        assigned.remove(&cell);
    }
}

fn shifted_multiset_boxes(max_entry: u32, max_box_size: usize) -> Vec<Vec<ShiftedMultisetEntry>> {
    let mut alphabet = Vec::with_capacity(2 * max_entry as usize);
    for value in 1..=max_entry {
        alphabet.push(ShiftedMultisetEntry::primed(value));
        alphabet.push(ShiftedMultisetEntry::unprimed(value));
    }

    let mut boxes = Vec::new();
    let mut current = Vec::new();
    enumerate_boxes_rec(&alphabet, 0, max_box_size, &mut current, &mut boxes);
    boxes
}

fn enumerate_boxes_rec(
    alphabet: &[ShiftedMultisetEntry],
    start: usize,
    max_box_size: usize,
    current: &mut Vec<ShiftedMultisetEntry>,
    boxes: &mut Vec<Vec<ShiftedMultisetEntry>>,
) {
    if !current.is_empty() {
        boxes.push(current.clone());
    }
    if current.len() == max_box_size {
        return;
    }

    for index in start..alphabet.len() {
        let entry = alphabet[index];
        if entry.is_primed() && current.contains(&entry) {
            continue;
        }
        current.push(entry);
        let next_start = if entry.is_primed() { index + 1 } else { index };
        enumerate_boxes_rec(alphabet, next_start, max_box_size, current, boxes);
        current.pop();
    }
}

fn shifted_shape_cells(shape: &Partition) -> Vec<(usize, usize)> {
    let mut cells = Vec::new();
    for (row, &length) in shape.parts().iter().enumerate() {
        for offset in 0..length as usize {
            cells.push((row, row + offset));
        }
    }
    cells
}

fn is_strict_partition(shape: &Partition) -> bool {
    shape.parts().windows(2).all(|window| window[0] > window[1])
}

fn is_compatible_with_assigned(
    cell: (usize, usize),
    box_entries: &[ShiftedMultisetEntry],
    assigned: &BTreeMap<(usize, usize), Vec<ShiftedMultisetEntry>>,
) -> bool {
    let (row, col) = cell;
    if col > row {
        if let Some(left) = assigned.get(&(row, col - 1)) {
            if !row_relation(left, box_entries) {
                return false;
            }
        }
    }
    if row > 0 {
        if let Some(above) = assigned.get(&(row - 1, col)) {
            if !column_relation(above, box_entries) {
                return false;
            }
        }
    }
    true
}

fn row_relation(left: &[ShiftedMultisetEntry], right: &[ShiftedMultisetEntry]) -> bool {
    right
        .iter()
        .all(|z| left.iter().all(|entry| entry.less_u(z)))
}

fn column_relation(above: &[ShiftedMultisetEntry], below: &[ShiftedMultisetEntry]) -> bool {
    below
        .iter()
        .all(|z| above.iter().any(|entry| entry.less_p(z)))
}

fn rows_have_unprimed_minima(
    shape: &Partition,
    assigned: &BTreeMap<(usize, usize), Vec<ShiftedMultisetEntry>>,
) -> bool {
    for row in 0..shape.num_parts() {
        let Some(row_minimum) = assigned
            .iter()
            .filter(|((cell_row, _), _)| *cell_row == row)
            .flat_map(|(_, entries)| entries.iter())
            .min()
        else {
            return false;
        };
        if row_minimum.is_primed() {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn partition(parts: &[u32]) -> Partition {
        Partition::new(parts.to_vec())
    }

    #[test]
    fn test_hawkes_degree_four_example_count() {
        let shape = partition(&[2, 1]);
        let tableaux = ShiftedMultisetTableau::enumerate_smt0(&shape, 2, 4);
        assert_eq!(tableaux.len(), 8);
    }

    #[test]
    fn test_hawkes_degree_four_example_distribution() {
        let shape = partition(&[2, 1]);
        let distribution = shifted_multiset_tableau_distribution(&shape, 2, 4);

        assert_eq!(distribution[&(vec![3, 1], vec![1, 0])], 1, "x_1^3 x_2 t_1");
        assert_eq!(distribution[&(vec![3, 1], vec![0, 1])], 1, "x_1^3 x_2 t_2");
        assert_eq!(
            distribution[&(vec![2, 2], vec![1, 0])],
            2,
            "x_1^2 x_2^2 t_1"
        );
        assert_eq!(
            distribution[&(vec![2, 2], vec![0, 1])],
            2,
            "x_1^2 x_2^2 t_2"
        );
        assert_eq!(distribution[&(vec![1, 3], vec![1, 0])], 1, "x_1 x_2^3 t_1");
        assert_eq!(distribution[&(vec![1, 3], vec![0, 1])], 1, "x_1 x_2^3 t_2");
        assert_eq!(distribution.len(), 6);
    }
}
