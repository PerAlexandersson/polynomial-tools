use sym_poly_core::Composition;
use sym_poly_qsym::{
    fundamental_in_hivert_expansion, hivert_fundamental_expansion, hivert_monomial_expansion,
};

fn composition(parts: &[u32]) -> Composition {
    Composition::new(parts.to_vec())
}

fn main() {
    let alpha = composition(&[3]);

    let g_in_f = hivert_fundamental_expansion(&alpha);
    let g_in_m = hivert_monomial_expansion(&alpha);
    let f_in_g = fundamental_in_hivert_expansion(&alpha);

    assert_eq!(g_in_f.coefficient(&composition(&[3])).coeffs(), &[1]);
    assert_eq!(g_in_f.coefficient(&composition(&[2, 1])).coeffs(), &[0, -1]);
    assert_eq!(g_in_f.coefficient(&composition(&[1, 2])).coeffs(), &[0, -1]);
    assert_eq!(
        g_in_f.coefficient(&composition(&[1, 1, 1])).coeffs(),
        &[0, 0, 1]
    );
    assert_eq!(g_in_m.coefficient(&composition(&[3])).coeffs(), &[1]);
    assert_eq!(g_in_m.coefficient(&composition(&[2, 1])).coeffs(), &[1, -1]);
    assert_eq!(g_in_m.coefficient(&composition(&[1, 2])).coeffs(), &[1, -1]);
    assert_eq!(
        g_in_m.coefficient(&composition(&[1, 1, 1])).coeffs(),
        &[1, -2, 1]
    );
    assert_eq!(f_in_g[&composition(&[3])].coeffs(), &[1]);
    assert_eq!(f_in_g[&composition(&[2, 1])].coeffs(), &[0, 1]);
    assert_eq!(f_in_g[&composition(&[1, 2])].coeffs(), &[0, 1]);
    assert_eq!(f_in_g[&composition(&[1, 1, 1])].coeffs(), &[0, 0, 0, 1]);

    println!("G_3 in the F-basis: {g_in_f}");
    println!("G_3 in the M-basis: {g_in_m}");
    println!("F_3 in the G-basis:");
    for (beta, coeff) in f_in_g {
        println!("  {coeff} G_{beta}");
    }
}
