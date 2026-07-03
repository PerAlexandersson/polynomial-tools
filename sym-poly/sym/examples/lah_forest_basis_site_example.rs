use combinatoric_core::{Graph, Partition};
use sym_poly_sym::{
    chromatic_symmetric, lah_forest_basis_elementary, lah_symmetric_forest_basis_expansion,
    SymmetricFunction,
};

fn partition(parts: &[u32]) -> Partition {
    Partition::new(parts.to_vec())
}

fn main() {
    let ell2 = lah_forest_basis_elementary(&partition(&[2]));
    let ell11 = lah_forest_basis_elementary(&partition(&[1, 1]));
    let k2_chromatic = chromatic_symmetric::<i64>(&Graph::complete(2)).to_elementary_basis();
    let lah_expansion = ell2.clone() - ell11.clone();

    assert_eq!(
        k2_chromatic, lah_expansion,
        "X_K2 should equal ell_2 - ell_11"
    );

    let two_e2 = SymmetricFunction::elementary_symmetric(partition(&[2])).scale(&2);
    assert_eq!(k2_chromatic, two_e2);

    let l42 = lah_symmetric_forest_basis_expansion(4, 2);
    assert_eq!(l42[&partition(&[2])], 4);
    assert_eq!(l42[&partition(&[1, 1])], 3);

    println!("ell_2 = {ell2}");
    println!("ell_11 = {ell11}");
    println!("X_K2 = ell_2 - ell_11 = {k2_chromatic}");
    println!("L_{{4,2}} = 4 ell_2 + 3 ell_11");
}
