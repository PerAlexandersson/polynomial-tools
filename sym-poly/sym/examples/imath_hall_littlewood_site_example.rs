use std::collections::BTreeMap;

type Shape = Vec<u32>;
type Coefficient = BTreeMap<(usize, usize), i64>; // (theta degree, t degree)
type Expansion = BTreeMap<Shape, Coefficient>;

fn monomial(theta_degree: usize, t_degree: usize, coefficient: i64) -> Coefficient {
    BTreeMap::from([((theta_degree, t_degree), coefficient)])
}

fn add_term(coefficient: &mut Coefficient, theta_degree: usize, t_degree: usize, scalar: i64) {
    *coefficient.entry((theta_degree, t_degree)).or_insert(0) += scalar;
    coefficient.retain(|_, value| *value != 0);
}

fn imath_v_two_row(a: u32, b: u32) -> Expansion {
    assert!(a >= b);
    let mut expansion = BTreeMap::from([(vec![a, b], monomial(0, 0, 1))]);

    for s in 1..=b {
        let mut coefficient = Coefficient::new();
        add_term(&mut coefficient, s as usize, s as usize, 1);
        add_term(&mut coefficient, s as usize, s as usize - 1, -1);
        expansion.insert(vec![a - s, b - s], coefficient);
    }

    expansion
}

fn specialize_theta_zero(expansion: &Expansion) -> Expansion {
    expansion
        .iter()
        .filter_map(|(shape, coefficient)| {
            let coefficient = coefficient
                .iter()
                .filter_map(|(&(theta_degree, t_degree), &scalar)| {
                    (theta_degree == 0).then_some(((theta_degree, t_degree), scalar))
                })
                .collect::<Coefficient>();
            (!coefficient.is_empty()).then_some((shape.clone(), coefficient))
        })
        .collect()
}

fn specialize_t_zero(expansion: &Expansion) -> Expansion {
    expansion
        .iter()
        .filter_map(|(shape, coefficient)| {
            let coefficient = coefficient
                .iter()
                .filter_map(|(&(theta_degree, t_degree), &scalar)| {
                    (t_degree == 0).then_some(((theta_degree, t_degree), scalar))
                })
                .collect::<Coefficient>();
            (!coefficient.is_empty()).then_some((shape.clone(), coefficient))
        })
        .collect()
}

fn main() {
    let expansion = imath_v_two_row(2, 1);
    assert_eq!(
        expansion,
        BTreeMap::from([
            (vec![1, 0], BTreeMap::from([((1, 0), -1), ((1, 1), 1)])),
            (vec![2, 1], monomial(0, 0, 1)),
        ])
    );

    assert_eq!(
        specialize_theta_zero(&expansion),
        BTreeMap::from([(vec![2, 1], monomial(0, 0, 1))])
    );
    assert_eq!(
        specialize_t_zero(&expansion),
        BTreeMap::from([
            (vec![1, 0], monomial(1, 0, -1)),
            (vec![2, 1], monomial(0, 0, 1)),
        ])
    );

    println!("V^i_(2,1) = V_(2,1) + (t - 1) theta V_(1)");
    println!("theta=0 leaves V_(2,1)");
    println!("t=0 gives the i-Schur specialization V_(2,1) - theta V_(1)");
}
