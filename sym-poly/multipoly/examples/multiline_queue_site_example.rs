use sym_poly_multipoly::MultilineQueue;

fn format_weight(weight: &[u32]) -> String {
    weight
        .iter()
        .enumerate()
        .filter_map(|(idx, &exp)| match exp {
            0 => None,
            1 => Some(format!("x_{}", idx + 1)),
            _ => Some(format!("x_{}^{}", idx + 1, exp)),
        })
        .collect::<Vec<_>>()
        .join(" ")
}

fn main() {
    // Mandelshtam--Valencia-Porras, Example 3.16.
    let queue = MultilineQueue::from_rows(6, &[&[1, 2, 3, 4], &[1, 3, 5, 6], &[2, 3], &[3, 5]]);
    let labeling = queue.ferrari_martin_labeling();

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
    assert_eq!(queue.content_weight(), vec![2, 2, 4, 1, 2, 1]);

    println!("rows, bottom to top: {:?}", queue.row_word_rows());
    println!("column word by columns: {:?}", queue.column_word_columns());
    println!("Ferrari--Martin pairings: {:?}", labeling.pairing_summary());
    println!(
        "weight: {} q^{}",
        format_weight(&queue.content_weight()),
        labeling.major_index()
    );
}
