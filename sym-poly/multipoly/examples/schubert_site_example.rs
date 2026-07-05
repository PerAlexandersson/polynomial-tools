use sym_poly_multipoly::{schubert_polynomial, MultiPoly};

fn main() {
    let s1 = schubert_polynomial::<i64>(&[2, 1, 3]);
    assert_eq!(s1, monomial(3, &[1, 0, 0]));

    let s2 = schubert_polynomial::<i64>(&[1, 3, 2]);
    assert_eq!(s2, monomial(3, &[1, 0, 0]) + monomial(3, &[0, 1, 0]));

    let s231 = schubert_polynomial::<i64>(&[2, 3, 1]);
    assert_eq!(s231, monomial(3, &[1, 1, 0]));

    let s312 = schubert_polynomial::<i64>(&[3, 1, 2]);
    assert_eq!(s312, monomial(3, &[2, 0, 0]));

    let involution_321 = s231.clone() + s312.clone();
    assert_eq!(
        involution_321,
        monomial(3, &[1, 1, 0]) + monomial(3, &[2, 0, 0])
    );

    let word = [5, 2, 1, 3, 4, 5];
    let perm = apply_simple_word_left_to_right(6, &word);
    assert_eq!(perm, vec![3, 1, 4, 6, 5, 2]);
    assert_eq!(word.len(), inversion_count(&perm));

    println!("S_213 = {}", s1);
    println!("S_132 = {}", s2);
    println!("S_231 + S_312 = {}", involution_321);
    println!("521345 is a reduced word for {:?}", perm);
}

fn monomial(num_vars: usize, exponents: &[u32]) -> MultiPoly<i64> {
    MultiPoly::x_power(num_vars, exponents.to_vec())
}

fn apply_simple_word_left_to_right(n: usize, word: &[usize]) -> Vec<usize> {
    let mut perm: Vec<usize> = (1..=n).collect();
    for &generator in word {
        assert!(
            (1..n).contains(&generator),
            "simple transposition s_{generator} is out of range for S_{n}"
        );
        perm.swap(generator - 1, generator);
    }
    perm
}

fn inversion_count(perm: &[usize]) -> usize {
    let mut count = 0usize;
    for i in 0..perm.len() {
        for j in i + 1..perm.len() {
            if perm[i] > perm[j] {
                count += 1;
            }
        }
    }
    count
}
