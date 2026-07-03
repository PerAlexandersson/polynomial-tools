use sym_poly_core::{Partition, UnivariatePolynomial};
use sym_poly_sym::unicellular_llt;

fn p(parts: &[u32]) -> Partition {
    Partition::new(parts.to_vec())
}

fn q(coeffs: &[i64]) -> UnivariatePolynomial<i64> {
    UnivariatePolynomial::new(coeffs.to_vec())
}

fn main() {
    let llt = unicellular_llt(&[0, 1, 2]);
    assert_eq!(llt.coefficient(&p(&[3])), q(&[1]));
    assert_eq!(llt.coefficient(&p(&[2, 1])), q(&[1, 1, 1]));
    assert_eq!(llt.coefficient(&p(&[1, 1, 1])), q(&[1, 2, 2, 1]));

    let schur = llt.to_schur_basis();
    assert_eq!(schur.coefficient(&p(&[3])), q(&[1]));
    assert_eq!(schur.coefficient(&p(&[2, 1])), q(&[0, 1, 1]));
    assert_eq!(schur.coefficient(&p(&[1, 1, 1])), q(&[0, 0, 0, 1]));

    println!("monomial expansion:");
    println!("{llt}");
    println!("Schur expansion:");
    println!("{schur}");
}
