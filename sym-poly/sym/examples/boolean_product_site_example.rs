use sym_poly_core::Partition;
use sym_poly_multipoly::MultiPoly;
use sym_poly_sym::SymmetricFunction;

fn main() {
    let n = 3;
    let product = (MultiPoly::<i64>::var(n, 0) + MultiPoly::var(n, 1))
        * (MultiPoly::var(n, 0) + MultiPoly::var(n, 2))
        * (MultiPoly::var(n, 1) + MultiPoly::var(n, 2));

    for exponent in [
        [2, 1, 0],
        [2, 0, 1],
        [1, 2, 0],
        [0, 2, 1],
        [1, 0, 2],
        [0, 1, 2],
    ] {
        assert_eq!(product.coefficient(&exponent), 1);
    }
    assert_eq!(product.coefficient(&[1, 1, 1]), 2);
    assert_eq!(product.terms().len(), 7);

    let schur_21 = SymmetricFunction::<i64>::schur_symmetric(Partition::new(vec![2, 1]));
    let monomial = schur_21.to_monomial_basis();
    assert_eq!(monomial.coefficient(&Partition::new(vec![2, 1])), 1);
    assert_eq!(monomial.coefficient(&Partition::new(vec![1, 1, 1])), 2);
    assert_eq!(monomial.terms().len(), 2);

    println!("B_{{3,2}} = {product}");
    println!("s_21 in the monomial basis: {monomial}");
}
