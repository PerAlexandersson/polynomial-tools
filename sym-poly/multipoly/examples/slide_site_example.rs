use sym_poly_multipoly::{
    fundamental_slide_polynomial, glide_polynomial, monomial_slide_polynomial, MultiPoly,
};

fn main() {
    let alpha = [1, 0, 2];
    let monomial: MultiPoly<i64> = monomial_slide_polynomial(&alpha);
    let fundamental: MultiPoly<i64> = fundamental_slide_polynomial(&alpha);

    assert_eq!(monomial.coefficient(&[1, 0, 2]), 1);
    assert_eq!(monomial.coefficient(&[1, 2, 0]), 1);
    assert_eq!(monomial.terms().len(), 2);

    assert_eq!(fundamental.coefficient(&[1, 0, 2]), 1);
    assert_eq!(fundamental.coefficient(&[1, 1, 1]), 1);
    assert_eq!(fundamental.coefficient(&[1, 2, 0]), 1);
    assert_eq!(fundamental.terms().len(), 3);

    let glide_at_zero: MultiPoly<i64> = glide_polynomial(&alpha, &0);
    assert_eq!(glide_at_zero, fundamental);

    let glide_at_two: MultiPoly<i64> = glide_polynomial(&alpha, &2);
    assert_eq!(glide_at_two.coefficient(&[1, 0, 2]), 1);
    assert_eq!(glide_at_two.coefficient(&[1, 1, 1]), 1);
    assert_eq!(glide_at_two.coefficient(&[1, 2, 0]), 1);
    assert_eq!(glide_at_two.coefficient(&[1, 1, 2]), 2);
    assert_eq!(glide_at_two.coefficient(&[1, 2, 1]), 2);

    println!("M_(1,0,2) = {monomial}");
    println!("F_(1,0,2) = {fundamental}");
    println!("G_(1,0,2)(beta=2) = {glide_at_two}");
}
