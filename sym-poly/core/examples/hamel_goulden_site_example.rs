use sym_poly_core::{ContentInterval, CuttingStripSegment, OutsideDecomposition};

fn format_row(row: &[CuttingStripSegment]) -> String {
    row.iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join("  ")
}

fn main() {
    let decomposition = OutsideDecomposition::from_intervals([(-3, 2), (-1, 1)]);
    let matrix = decomposition.determinant_matrix();

    println!("ribbons:");
    for (idx, interval) in decomposition.ribbons().iter().enumerate() {
        println!("  theta_{} = theta{}", idx + 1, interval);
    }

    println!("Hamel--Goulden matrix:");
    for row in &matrix {
        println!("  {}", format_row(row));
    }

    assert_eq!(
        matrix,
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
