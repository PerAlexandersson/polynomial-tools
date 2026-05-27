use std::collections::BTreeMap;

use sym_poly_multipoly::{partial_i, MultiPoly};

fn main() {
    let fixed_lambda = vec![4, 3, 1];
    let fixed_flags = vec![1, 3, 5];
    check_wachs_case("fixed-row-2", &fixed_lambda, &fixed_flags, 3, 5);
    check_zero_case("fixed-zero", &fixed_lambda, &fixed_flags, 2, 5);
    check_plateau_case("nonlowerable-zero", &[2, 2], &[1, 2], 1, 1, 2);
    check_plateau_case("nonlowerable-two-terms", &[2, 2], &[1, 4], 1, 1, 4);
    check_plateau_case("nonlowerable-chain", &[3, 3, 3], &[1, 3, 5], 1, 2, 5);

    for n in 0..=3 {
        let lambda = vec![n + 4, n + 2, 1];
        let flags = vec![1, 3, 5];
        check_wachs_case(&format!("affine-N-{n}"), &lambda, &flags, 3, 5);
    }
}

fn check_wachs_case(name: &str, lambda: &[usize], flags: &[usize], i: usize, num_vars: usize) {
    let row = flags
        .iter()
        .position(|&b| b == i)
        .expect("test expects i to be a flag");
    assert!(
        row + 1 == lambda.len() || lambda[row] > lambda[row + 1],
        "test expects the selected row to be lowerable"
    );

    let f = flagged_schur(lambda, flags, num_vars);
    let lhs = partial_i(&f, i - 1);

    let mut rhs_lambda = lambda.to_vec();
    rhs_lambda[row] -= 1;
    let mut rhs_flags = flags.to_vec();
    rhs_flags[row] += 1;
    let rhs = flagged_schur(&rhs_lambda, &rhs_flags, num_vars);

    println!(
        "{name}: lambda={lambda:?}, flags={flags:?}, partial_{i} terms={}, rhs lambda={rhs_lambda:?}, rhs flags={rhs_flags:?}, rhs terms={}, match={}",
        lhs.terms().len(),
        rhs.terms().len(),
        lhs == rhs
    );

    assert_eq!(lhs, rhs, "{name} failed Wachs divided-difference formula");
}

fn check_zero_case(name: &str, lambda: &[usize], flags: &[usize], i: usize, num_vars: usize) {
    assert!(
        !flags.contains(&i),
        "zero test expects i not to occur as a flag"
    );
    let f = flagged_schur(lambda, flags, num_vars);
    let lhs = partial_i(&f, i - 1);
    println!(
        "{name}: lambda={lambda:?}, flags={flags:?}, partial_{i} zero={}",
        lhs.is_zero()
    );
    assert!(lhs.is_zero(), "{name} expected zero divided difference");
}

fn check_plateau_case(
    name: &str,
    lambda: &[usize],
    flags: &[usize],
    i: usize,
    plateau_end: usize,
    num_vars: usize,
) {
    let row = flags
        .iter()
        .position(|&b| b == i)
        .expect("test expects i to be a flag");
    assert!(row < plateau_end, "test expects a non-lowerable row");
    assert!(
        (row..=plateau_end).all(|j| lambda[j] == lambda[row]),
        "test expects a plateau of equal row lengths"
    );
    assert!(
        plateau_end + 1 == lambda.len() || lambda[plateau_end] > lambda[plateau_end + 1],
        "test expects plateau_end to be the last row in the plateau"
    );

    let f = flagged_schur(lambda, flags, num_vars);
    let lhs = partial_i(&f, i - 1);
    let rhs = plateau_formula(lambda, flags, i, row, plateau_end, num_vars);

    println!(
        "{name}: lambda={lambda:?}, flags={flags:?}, partial_{i} terms={}, plateau rhs terms={}, match={}",
        lhs.terms().len(),
        rhs.terms().len(),
        lhs == rhs
    );

    assert_eq!(lhs, rhs, "{name} failed plateau divided-difference formula");
}

fn plateau_formula(
    lambda: &[usize],
    flags: &[usize],
    i: usize,
    row: usize,
    plateau_end: usize,
    num_vars: usize,
) -> MultiPoly<i64> {
    let mut shape = lambda.to_vec();
    for part in &mut shape[row..=plateau_end] {
        *part -= 1;
    }

    let mut rhs = MultiPoly::zero(num_vars);
    let mut chosen = Vec::new();
    accumulate_plateau_terms(
        &shape,
        flags,
        i + 1,
        row + 1,
        plateau_end,
        num_vars,
        &mut chosen,
        &mut rhs,
    );
    rhs
}

fn accumulate_plateau_terms(
    shape: &[usize],
    flags: &[usize],
    previous_flag: usize,
    current_row: usize,
    plateau_end: usize,
    num_vars: usize,
    chosen: &mut Vec<usize>,
    rhs: &mut MultiPoly<i64>,
) {
    if current_row > plateau_end {
        let first_row = plateau_end - chosen.len();
        let mut term_flags = flags.to_vec();
        term_flags[first_row] += 1;
        for (offset, &flag) in chosen.iter().enumerate() {
            term_flags[first_row + 1 + offset] = flag;
        }

        let mut exp = vec![0u32; num_vars];
        for &flag in chosen.iter() {
            exp[flag - 1] += 1;
        }
        let monomial = MultiPoly::x_power(num_vars, exp);
        *rhs = rhs.clone() + monomial * flagged_schur(shape, &term_flags, num_vars);
        return;
    }

    for next_flag in previous_flag + 1..=flags[current_row] {
        chosen.push(next_flag);
        accumulate_plateau_terms(
            shape,
            flags,
            next_flag,
            current_row + 1,
            plateau_end,
            num_vars,
            chosen,
            rhs,
        );
        chosen.pop();
    }
}

fn flagged_schur(lambda: &[usize], flags: &[usize], num_vars: usize) -> MultiPoly<i64> {
    assert_eq!(lambda.len(), flags.len());
    assert!(flags.iter().all(|&b| b <= num_vars));

    let mut terms = BTreeMap::new();
    let mut filling = lambda
        .iter()
        .map(|&row_len| vec![0usize; row_len])
        .collect::<Vec<_>>();
    enumerate_cell(lambda, flags, num_vars, &mut filling, 0, 0, &mut terms);
    MultiPoly::from_terms(num_vars, terms)
}

fn enumerate_cell(
    lambda: &[usize],
    flags: &[usize],
    num_vars: usize,
    filling: &mut [Vec<usize>],
    row: usize,
    col: usize,
    terms: &mut BTreeMap<Vec<u32>, i64>,
) {
    if row == lambda.len() {
        let mut exp = vec![0u32; num_vars];
        for row_values in filling.iter() {
            for &entry in row_values {
                exp[entry - 1] += 1;
            }
        }
        *terms.entry(exp).or_insert(0) += 1;
        return;
    }

    if col == lambda[row] {
        enumerate_cell(lambda, flags, num_vars, filling, row + 1, 0, terms);
        return;
    }

    let row_min = if col == 0 { 1 } else { filling[row][col - 1] };
    let col_min = if row == 0 || col >= lambda[row - 1] {
        1
    } else {
        filling[row - 1][col] + 1
    };
    let min_entry = row_min.max(col_min);

    for entry in min_entry..=flags[row] {
        filling[row][col] = entry;
        enumerate_cell(lambda, flags, num_vars, filling, row, col + 1, terms);
    }
    filling[row][col] = 0;
}
