use sym_poly_multipoly::{
    flagged_schur, flagged_skew_schur, flagged_skew_tableaux, flagged_tableaux, MultiPoly,
    MultiPolyFunction,
};

fn format_monomial(exp: &[u32]) -> String {
    let factors = exp
        .iter()
        .enumerate()
        .filter_map(|(idx, &power)| match power {
            0 => None,
            1 => Some(format!("x_{}", idx + 1)),
            _ => Some(format!("x_{}^{}", idx + 1, power)),
        })
        .collect::<Vec<_>>();
    if factors.is_empty() {
        String::from("1")
    } else {
        factors.join(" ")
    }
}

fn print_poly(label: &str, poly: &MultiPoly<i64>) {
    println!("{label}");
    for (exp, coeff) in poly.terms() {
        println!("  {coeff} {}", format_monomial(exp));
    }
}

fn main() {
    let flagged = flagged_schur::<i64>(&[2, 1], &[2, 3], 3);
    assert_eq!(flagged_tableaux(&[2, 1], &[2, 3]).len(), 5);
    assert_eq!(flagged.terms().len(), 5);
    print_poly("s_(2,1), b=(2,3):", &flagged);

    let skew = flagged_skew_schur::<i64>(&[3, 2], &[1], &[2, 3], 3);
    assert_eq!(flagged_skew_tableaux(&[3, 2], &[1], &[2, 3]).len(), 13);
    print_poly("s_(3,2)/(1), b=(2,3):", &skew);

    let ordinary_skew = flagged_skew_schur::<i64>(&[3, 2], &[1], &[3, 3], 3);
    print_poly("ordinary skew specialization b=(3,3):", &ordinary_skew);

    let key_expansion = MultiPolyFunction::from_multipoly(&skew).to_key_basis();
    assert!(key_expansion.positive_coefficients());
    println!("key expansion:");
    for (alpha, coeff) in key_expansion.terms() {
        println!("  {coeff} kappa_{}", alpha);
    }
}
