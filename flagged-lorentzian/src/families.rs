use sym_poly_core::Partition;

use crate::shape::SkewShape;

/// Generate all partitions of `n`.
pub fn partitions_of_size(n: u32) -> Vec<Partition> {
    let mut results = Vec::new();
    partitions_with_max_part(n, n, &mut Vec::new(), &mut results);
    results
}

/// Generate all subpartitions of `outer`.
pub fn subpartitions(outer: &Partition) -> Vec<Partition> {
    let mut results = Vec::new();
    subpartitions_from_row(outer, 0, u32::MAX, &mut Vec::new(), &mut results);
    results
}

/// Generate skew shapes with fixed skew size and bounded outer size.
pub fn skew_shapes_of_size(skew_size: u32, max_outer_size: u32) -> Vec<SkewShape> {
    let mut shapes = Vec::new();
    for outer_size in skew_size..=max_outer_size {
        for outer in partitions_of_size(outer_size) {
            for inner in subpartitions(&outer) {
                if outer.size() - inner.size() == skew_size {
                    shapes.push(SkewShape::new(outer.clone(), inner));
                }
            }
        }
    }
    shapes
}

fn partitions_with_max_part(
    remaining: u32,
    max_part: u32,
    current: &mut Vec<u32>,
    results: &mut Vec<Partition>,
) {
    if remaining == 0 {
        results.push(Partition::from_sorted(current.clone()));
        return;
    }

    for part in (1..=remaining.min(max_part)).rev() {
        current.push(part);
        partitions_with_max_part(remaining - part, part, current, results);
        current.pop();
    }
}

fn subpartitions_from_row(
    outer: &Partition,
    row: usize,
    previous_part: u32,
    current: &mut Vec<u32>,
    results: &mut Vec<Partition>,
) {
    if row == outer.num_parts() {
        results.push(Partition::from_sorted(current.clone()));
        return;
    }

    let max_part = previous_part.min(outer.part(row));
    for part in (0..=max_part).rev() {
        current.push(part);
        subpartitions_from_row(outer, row + 1, part, current, results);
        current.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partitions_of_four() {
        assert_eq!(partitions_of_size(4).len(), 5);
    }

    #[test]
    fn skew_shapes_include_disconnected_counterexample() {
        let found = skew_shapes_of_size(4, 8)
            .into_iter()
            .any(|shape| shape.outer().parts() == [4, 3, 1] && shape.inner().parts() == [3, 1]);
        assert!(found);
    }
}
