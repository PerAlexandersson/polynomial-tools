use sym_poly_multipoly::borodin_wheeler_l_weight;

fn main() {
    // Borodin--Wheeler, Examples 3.2 and 3.3.
    let example_3_2 = borodin_wheeler_l_weight(&[1, 1, 1], 1, &[2, 0, 1], 2)
        .expect("Example 3.2 has nonzero weight");
    let example_3_3 = borodin_wheeler_l_weight(&[2, 1, 2], 0, &[1, 1, 2], 1)
        .expect("Example 3.3 has nonzero weight");

    assert_eq!(example_3_2.factored_string(), "x (1 - t) t");
    assert_eq!(example_3_3.factored_string(), "x (1 - t^2) t^3");

    println!(
        "L_x((1,1,1), 1; (2,0,1), 2) = {}",
        example_3_2.factored_string()
    );
    println!(
        "L_x((2,1,2), 0; (1,1,2), 1) = {}",
        example_3_3.factored_string()
    );
}
