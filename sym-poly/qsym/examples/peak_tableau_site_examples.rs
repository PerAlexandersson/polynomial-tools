use sym_poly_qsym::{StandardPeakCompositionTableau, StandardPeakYoungCompositionTableau};

fn main() {
    let spct = StandardPeakCompositionTableau {
        rows: vec![vec![1, 2, 3, 4], vec![5, 6, 7], vec![8]],
    };
    let spyct = StandardPeakYoungCompositionTableau {
        rows: spct.rows.clone(),
    };

    assert!(StandardPeakCompositionTableau::enumerate(&[4, 3, 1]).contains(&spct));
    assert!(StandardPeakYoungCompositionTableau::enumerate(&[4, 3, 1]).contains(&spyct));

    println!("\\begin{{ytableau}}");
    for row in spct.rows.iter().rev() {
        println!(
            "{} \\\\",
            row.iter().map(u32::to_string).collect::<Vec<_>>().join("&")
        );
    }
    println!("\\end{{ytableau}}");
    println!("SPCT upward descents: {:?}", spct.upward_descent_set());
    println!("SPCT peak composition: {}", spct.upward_peak_composition());
    println!("SPYCT left descents: {:?}", spyct.left_descent_set());
    println!("SPYCT peak composition: {}", spyct.left_peak_composition());
}
