use std::fmt;

/// A content interval `[start, end]` on a Hamel--Goulden cutting strip.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ContentInterval {
    start: i32,
    end: i32,
}

impl ContentInterval {
    /// Create a non-empty content interval.
    pub fn new(start: i32, end: i32) -> Self {
        assert!(start <= end, "content interval must be non-empty");
        Self { start, end }
    }

    /// The initial content.
    pub fn start(&self) -> i32 {
        self.start
    }

    /// The terminal content.
    pub fn end(&self) -> i32 {
        self.end
    }
}

impl fmt::Display for ContentInterval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{},{}]", self.start, self.end)
    }
}

/// A segment of a Hamel--Goulden cutting strip.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CuttingStripSegment {
    /// The segment `theta[start, end]`.
    Segment(ContentInterval),
    /// The empty ribbon `theta[q + 1, q]`, whose Schur function is `1`.
    Empty,
    /// An undefined segment `theta[p, q]` with `p > q + 1`, whose Schur
    /// function is treated as `0` in the determinant.
    Undefined,
}

impl CuttingStripSegment {
    /// Construct the segment `theta[start, end]` using the standard empty and
    /// undefined conventions.
    pub fn from_bounds(start: i32, end: i32) -> Self {
        if start <= end {
            CuttingStripSegment::Segment(ContentInterval::new(start, end))
        } else if start == end + 1 {
            CuttingStripSegment::Empty
        } else {
            CuttingStripSegment::Undefined
        }
    }
}

impl fmt::Display for CuttingStripSegment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CuttingStripSegment::Segment(interval) => write!(f, "theta{interval}"),
            CuttingStripSegment::Empty => write!(f, "empty"),
            CuttingStripSegment::Undefined => write!(f, "undefined"),
        }
    }
}

/// Content-interval data for an outside decomposition.
///
/// If the decomposition has ribbons `theta_i` with intervals
/// `[p(theta_i), q(theta_i)]`, the Hamel--Goulden matrix has `(i,j)` entry
/// `theta_i # theta_j = theta[p(theta_j), q(theta_i)]`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutsideDecomposition {
    ribbons: Vec<ContentInterval>,
}

impl OutsideDecomposition {
    /// Create an outside decomposition from ribbon content intervals.
    pub fn from_intervals<I>(intervals: I) -> Self
    where
        I: IntoIterator<Item = (i32, i32)>,
    {
        let ribbons = intervals
            .into_iter()
            .map(|(start, end)| ContentInterval::new(start, end))
            .collect();
        Self { ribbons }
    }

    /// Create an outside decomposition from ribbon cells.
    ///
    /// Cells are `(row, column)` pairs in any consistent integer coordinate
    /// system, with content `column - row`.  This helper records only content
    /// intervals; it does not verify that the cells form valid outside ribbons.
    pub fn from_ribbon_cells(ribbons: &[Vec<(i32, i32)>]) -> Self {
        let intervals = ribbons.iter().map(|cells| {
            assert!(!cells.is_empty(), "ribbon cell list must be non-empty");
            let mut contents = cells.iter().map(|&(row, col)| col - row);
            let first = contents.next().expect("checked non-empty");
            let (start, end) = contents.fold((first, first), |(lo, hi), c| (lo.min(c), hi.max(c)));
            (start, end)
        });
        Self::from_intervals(intervals)
    }

    /// The ribbon intervals, in determinant order.
    pub fn ribbons(&self) -> &[ContentInterval] {
        &self.ribbons
    }

    /// The Hamel--Goulden matrix of cutting-strip segments.
    pub fn determinant_matrix(&self) -> Vec<Vec<CuttingStripSegment>> {
        self.ribbons
            .iter()
            .map(|row_ribbon| {
                self.ribbons
                    .iter()
                    .map(|col_ribbon| {
                        CuttingStripSegment::from_bounds(col_ribbon.start, row_ribbon.end)
                    })
                    .collect()
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::{ContentInterval, CuttingStripSegment, OutsideDecomposition};

    #[test]
    fn hamel_goulden_matrix_matches_site_example() {
        let decomposition = OutsideDecomposition::from_intervals([(-3, 2), (-1, 1)]);
        assert_eq!(
            decomposition.ribbons(),
            &[ContentInterval::new(-3, 2), ContentInterval::new(-1, 1)]
        );
        assert_eq!(
            decomposition.determinant_matrix(),
            vec![
                vec![
                    CuttingStripSegment::Segment(ContentInterval::new(-3, 2)),
                    CuttingStripSegment::Segment(ContentInterval::new(-1, 2)),
                ],
                vec![
                    CuttingStripSegment::Segment(ContentInterval::new(-3, 1)),
                    CuttingStripSegment::Segment(ContentInterval::new(-1, 1)),
                ],
            ]
        );
    }

    #[test]
    fn cutting_strip_segment_conventions() {
        assert_eq!(
            CuttingStripSegment::from_bounds(2, 3),
            CuttingStripSegment::Segment(ContentInterval::new(2, 3))
        );
        assert_eq!(
            CuttingStripSegment::from_bounds(4, 3),
            CuttingStripSegment::Empty
        );
        assert_eq!(
            CuttingStripSegment::from_bounds(5, 3),
            CuttingStripSegment::Undefined
        );
    }

    #[test]
    fn intervals_can_be_read_from_cells() {
        let decomposition = OutsideDecomposition::from_ribbon_cells(&[
            vec![(3, 0), (2, 0), (1, 1), (0, 2)],
            vec![(2, 1), (1, 2)],
        ]);
        assert_eq!(
            decomposition.ribbons(),
            &[ContentInterval::new(-3, 2), ContentInterval::new(-1, 1)]
        );
    }
}
