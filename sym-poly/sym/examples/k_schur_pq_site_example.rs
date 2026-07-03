use std::collections::{BTreeMap, BTreeSet};

type Partition = Vec<u32>;
type BetaPolynomial = BTreeMap<usize, i64>;
type Expansion = BTreeMap<Partition, BetaPolynomial>;

fn is_strict(partition: &[u32]) -> bool {
    partition.windows(2).all(|pair| pair[0] > pair[1])
}

fn shifted_diagram(partition: &[u32]) -> BTreeSet<(u32, u32)> {
    let mut cells = BTreeSet::new();
    for (row, &part) in partition.iter().enumerate() {
        let row = row as u32 + 1;
        for col in row..row + part {
            cells.insert((row, col));
        }
    }
    cells
}

fn skew_size_and_columns(lambda: &[u32], mu: &[u32]) -> (usize, usize) {
    let lambda_cells = shifted_diagram(lambda);
    let mu_cells = shifted_diagram(mu);
    let skew_cells = lambda_cells
        .difference(&mu_cells)
        .copied()
        .collect::<Vec<_>>();
    let columns = skew_cells
        .iter()
        .map(|&(_, col)| col)
        .collect::<BTreeSet<_>>();

    (skew_cells.len(), columns.len())
}

fn beta_monomial(degree: usize, coefficient: i64) -> BetaPolynomial {
    BTreeMap::from([(degree, coefficient)])
}

fn q_to_p_coefficients(mu: &[u32]) -> Expansion {
    assert!(is_strict(mu));
    let len = mu.len();
    let base_power = 1_i64 << len;
    let mut expansion = Expansion::new();

    for mask in 0..(1_usize << len) {
        let mut lambda = mu.to_vec();
        for (row, part) in lambda.iter_mut().enumerate() {
            if (mask >> row) & 1 == 1 {
                *part += 1;
            }
        }

        if !is_strict(&lambda) {
            continue;
        }

        let (size, columns) = skew_size_and_columns(&lambda, mu);
        let sign = if (size + columns) % 2 == 0 { 1 } else { -1 };
        let coefficient = sign * (base_power >> size);
        expansion.insert(lambda, beta_monomial(size, coefficient));
    }

    expansion
}

fn format_beta_polynomial(polynomial: &BetaPolynomial) -> String {
    polynomial
        .iter()
        .map(|(&degree, &coefficient)| match (coefficient, degree) {
            (1, 0) => "1".to_string(),
            (-1, 0) => "-1".to_string(),
            (c, 0) => c.to_string(),
            (1, 1) => "beta".to_string(),
            (-1, 1) => "-beta".to_string(),
            (c, 1) => format!("{c} beta"),
            (1, d) => format!("beta^{d}"),
            (-1, d) => format!("-beta^{d}"),
            (c, d) => format!("{c} beta^{d}"),
        })
        .collect::<Vec<_>>()
        .join(" + ")
}

fn format_expansion(expansion: &Expansion) -> String {
    let mut terms = Vec::new();
    for (partition, coefficient) in expansion {
        let coefficient = format_beta_polynomial(coefficient);
        let term = format!("{} GP_{partition:?}", coefficient.trim_start_matches('-'));
        if coefficient.starts_with('-') {
            terms.push(format!("- {term}"));
        } else if terms.is_empty() {
            terms.push(term);
        } else {
            terms.push(format!("+ {term}"));
        }
    }
    terms.join(" ")
}

fn main() {
    let mu_32 = q_to_p_coefficients(&[3, 2]);
    assert_eq!(
        mu_32,
        BTreeMap::from([
            (vec![3, 2], beta_monomial(0, 4)),
            (vec![4, 2], beta_monomial(1, 2)),
            (vec![4, 3], beta_monomial(2, -1)),
        ])
    );

    let mu_4 = q_to_p_coefficients(&[4]);
    assert_eq!(
        mu_4,
        BTreeMap::from([
            (vec![4], beta_monomial(0, 2)),
            (vec![5], beta_monomial(1, 1)),
        ])
    );

    let mu_42 = q_to_p_coefficients(&[4, 2]);
    assert!(mu_42
        .values()
        .flat_map(|polynomial| polynomial.values())
        .all(|&coefficient| coefficient >= 0));

    println!("GQ_(3,2) = {}", format_expansion(&mu_32));
    println!("GQ_(4) = {}", format_expansion(&mu_4));
    println!("GQ_(4,2) has only nonnegative beta coefficients:");
    println!("{}", format_expansion(&mu_42));
}
