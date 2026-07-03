use sym_poly_multipoly::{fundamental_slide_polynomial, monomial_slide_polynomial, MultiPoly};

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

    println!("M_(1,0,2) = {monomial}");
    println!("F_(1,0,2) = {fundamental}");
}
