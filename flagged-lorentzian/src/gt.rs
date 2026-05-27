use crate::shape::SkewShape;

pub type GtRow = Vec<u32>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SkewGtPattern {
    rows: Vec<GtRow>,
}

impl SkewGtPattern {
    pub fn from_tableau(shape: &SkewShape, values: &[u32], alphabet: usize) -> Self {
        assert_eq!(
            shape.size(),
            values.len(),
            "tableau values must match the skew shape size"
        );

        let mut current = inner_row(shape);
        let mut rows = Vec::with_capacity(alphabet + 1);
        rows.push(current.clone());

        for level in 1..=alphabet as u32 {
            for (cell, &value) in shape.cells().iter().zip(values) {
                if value == level {
                    current[cell.row] += 1;
                }
            }
            rows.push(current.clone());
        }

        Self { rows }
    }

    pub fn rows(&self) -> &[GtRow] {
        &self.rows
    }

    pub fn row(&self, level: u32) -> &[u32] {
        &self.rows[level as usize]
    }

    pub fn alphabet(&self) -> usize {
        self.rows.len().saturating_sub(1)
    }
}

pub fn inner_row(shape: &SkewShape) -> GtRow {
    (0..shape.row_count())
        .map(|idx| shape.inner().part(idx))
        .collect()
}

pub fn active_gt_row(shape: &SkewShape, values: &[u32], lower_label: u32) -> GtRow {
    SkewGtPattern::from_tableau(shape, values, lower_label as usize)
        .row(lower_label)
        .to_vec()
}

pub fn sharp_flag(shape: &SkewShape, values: &[u32]) -> Vec<u32> {
    let mut raw = vec![0u32; shape.row_count()];
    for (cell, &value) in shape.cells().iter().zip(values) {
        raw[cell.row] = raw[cell.row].max(value);
    }

    let mut sharp = raw;
    for row in 1..sharp.len() {
        sharp[row] = sharp[row].max(sharp[row - 1]);
    }
    sharp
}

pub fn pair_envelope(left: &[u32], right: &[u32]) -> Vec<u32> {
    left.iter()
        .zip(right)
        .map(|(&left, &right)| left.max(right))
        .collect()
}

pub fn add_rows(left: &[u32], right: &[u32]) -> GtRow {
    left.iter()
        .zip(right)
        .map(|(&left, &right)| left + right)
        .collect()
}

pub fn subtract_rows(left: &[u32], right: &[u32]) -> Vec<i32> {
    left.iter()
        .zip(right)
        .map(|(&left, &right)| left as i32 - right as i32)
        .collect()
}

pub fn add_patterns(left: &SkewGtPattern, right: &SkewGtPattern) -> Vec<GtRow> {
    assert_eq!(
        left.rows.len(),
        right.rows.len(),
        "GT patterns must use the same alphabet"
    );
    left.rows
        .iter()
        .zip(right.rows())
        .map(|(left, right)| add_rows(left, right))
        .collect()
}

pub fn subtract_pattern_sums(left: &[GtRow], right: &[GtRow]) -> Vec<Vec<i32>> {
    assert_eq!(
        left.len(),
        right.len(),
        "GT pattern sums must have the same number of rows"
    );
    left.iter()
        .zip(right)
        .map(|(left, right)| subtract_rows(left, right))
        .collect()
}

pub fn is_gt_array(rows: &[GtRow]) -> bool {
    if rows.is_empty() {
        return true;
    }
    let width = rows[0].len();
    if rows.iter().any(|row| row.len() != width) {
        return false;
    }
    if rows.iter().any(|row| !is_weakly_decreasing(row)) {
        return false;
    }

    for level in 1..rows.len() {
        let upper = &rows[level];
        let lower = &rows[level - 1];
        for col in 0..width {
            if upper[col] < lower[col] {
                return false;
            }
            if col + 1 < width && lower[col] < upper[col + 1] {
                return false;
            }
        }
    }
    true
}

pub fn elementary_row_exchange_neighbors(rows: &[GtRow], fixed_level: usize) -> Vec<Vec<GtRow>> {
    let mut neighbors = Vec::new();
    if rows.is_empty() {
        return neighbors;
    }
    let width = rows[0].len();
    if width < 2 {
        return neighbors;
    }

    for level in 1..rows.len().saturating_sub(1) {
        if level == fixed_level {
            continue;
        }
        for source in 0..width {
            if rows[level][source] == 0 {
                continue;
            }
            for target in 0..width {
                if source == target {
                    continue;
                }
                let mut candidate = rows.to_vec();
                candidate[level][source] -= 1;
                candidate[level][target] += 1;
                if is_gt_array(&candidate) {
                    neighbors.push(candidate);
                }
            }
        }
    }

    neighbors.sort();
    neighbors.dedup();
    neighbors
}

fn is_weakly_decreasing(row: &[u32]) -> bool {
    row.windows(2).all(|window| window[0] >= window[1])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SkewShape;

    #[test]
    fn gt_rows_record_cumulative_shapes() {
        let shape = SkewShape::from_parts(vec![2, 1], vec![]);
        let values = vec![1, 3, 2];
        let gt = SkewGtPattern::from_tableau(&shape, &values, 3);

        assert_eq!(gt.rows(), &[vec![0, 0], vec![1, 0], vec![1, 1], vec![2, 1]]);
    }

    #[test]
    fn sharp_flag_is_prefix_maximum_of_row_maxima() {
        let shape = SkewShape::from_parts(vec![2, 1], vec![]);
        let values = vec![3, 3, 2];

        assert_eq!(sharp_flag(&shape, &values), vec![3, 3]);
    }

    #[test]
    fn recognizes_gt_arrays() {
        let valid = vec![vec![4, 0, 0], vec![5, 1, 0], vec![6, 1, 1], vec![6, 3, 1]];
        let invalid = vec![vec![4, 0, 0], vec![4, 2, 0], vec![6, 1, 1], vec![6, 3, 1]];

        assert!(is_gt_array(&valid));
        assert!(!is_gt_array(&invalid));
    }

    #[test]
    fn elementary_exchange_preserves_gt_validity() {
        let rows = vec![
            vec![4, 0, 0],
            vec![5, 1, 0],
            vec![6, 1, 1],
            vec![6, 3, 1],
            vec![6, 4, 2],
            vec![8, 4, 2],
        ];

        let neighbors = elementary_row_exchange_neighbors(&rows, 4);
        assert!(neighbors.iter().all(|neighbor| is_gt_array(neighbor)));
        assert!(neighbors.iter().any(|neighbor| {
            neighbor[2] == vec![6, 2, 0] && neighbor[0] == rows[0] && neighbor[4] == rows[4]
        }));
    }
}
