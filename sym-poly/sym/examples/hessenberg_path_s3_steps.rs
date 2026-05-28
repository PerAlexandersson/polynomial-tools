use std::collections::BTreeMap;

use num_rational::Ratio;
use sym_poly_core::linear_algebra::{
    kernel_basis, matrix_trace, matrix_vector_multiply, quotient_action_matrix, rank,
    solve_linear_system, Matrix, QuotientSpace, Vector,
};
use sym_poly_core::sn_action::{cycle_type, SnIndex};
use sym_poly_core::{symmetric_group_permutation_basis, Ring};

type Q = Ratio<i64>;

fn main() {
    let n = 3;
    let edges = vec![(0, 1), (1, 2)];
    let fixed_points = symmetric_group_permutation_basis(n);
    let sn_index = SnIndex::new(n);

    println!("Area sequence [0,1,1]");
    println!("Graph edges: {{1,2}}, {{2,3}}");
    println!();

    println!("Fixed point order:");
    for (index, point) in fixed_points.iter().enumerate() {
        println!("  w{index} = {}", format_permutation(point));
    }
    println!();

    let monomials_degree_zero = homogeneous_monomials(n, 0);
    let monomials_degree_one = homogeneous_monomials(n, 1);
    let degree_zero_basis = gkm_module_basis(n, &edges, &fixed_points, &sn_index, 0);
    let degree_one_basis = gkm_module_basis(n, &edges, &fixed_points, &sn_index, 1);

    println!("Step 1. Degree 0 equivariant GKM module M_0");
    for (index, basis_vector) in degree_zero_basis.iter().enumerate() {
        println!(
            "  b{index}: {}",
            format_gkm_assignment(basis_vector, &fixed_points, &monomials_degree_zero)
        );
    }
    println!("  dim M_0 = {}", degree_zero_basis.len());
    println!();

    println!("Step 2. Degree 1 equivariant GKM module M_1");
    for (index, basis_vector) in degree_one_basis.iter().enumerate() {
        println!(
            "  m{index}: {}",
            format_gkm_assignment(basis_vector, &fixed_points, &monomials_degree_one)
        );
    }
    println!("  dim M_1 = {}", degree_one_basis.len());
    println!();

    let quotient_relations = (0..n)
        .map(|variable| {
            let ambient =
                global_variable_assignment(variable, &fixed_points, &monomials_degree_one, n);
            coordinates_in_basis(&degree_one_basis, &ambient)
        })
        .collect::<Vec<_>>();

    println!("Step 3. Subspace R_1 = (t_1,t_2,t_3) M_0 to quotient out");
    for (variable, relation) in quotient_relations.iter().enumerate() {
        println!(
            "  r{} = t{} b0 = {}",
            variable + 1,
            variable + 1,
            format_basis_combination(relation, "m")
        );
    }
    println!("  rank R_1 = {}", rank(&quotient_relations));
    println!();

    let ordinary_quotient =
        QuotientSpace::from_relations(degree_one_basis.len(), &quotient_relations);
    println!("Step 4. Ordinary degree-one quotient H^2 = M_1/R_1");
    println!(
        "  quotient free module coordinates: {:?}",
        ordinary_quotient.free_columns
    );
    println!("  so [x1,x2,x3,x4] records normal-form coefficients of m1,m2,m5,m6");
    println!("  dim H^2 = {}", ordinary_quotient.dimension());
    println!();

    let module_actions = class_representatives()
        .into_iter()
        .map(|(name, permutation)| {
            let action = module_dot_action_matrix(
                &degree_one_basis,
                &fixed_points,
                &sn_index,
                &monomials_degree_one,
                &permutation,
            );
            let quotient_action = quotient_action_matrix(&ordinary_quotient, &action);
            (name, permutation, quotient_action)
        })
        .collect::<Vec<_>>();

    let s1_action = module_actions
        .iter()
        .find(|(name, _, _)| *name == "(12)")
        .unwrap()
        .2
        .clone();
    let s2_action = quotient_action_for_permutation(
        &degree_one_basis,
        &ordinary_quotient,
        &fixed_points,
        &sn_index,
        &monomials_degree_one,
        &[0, 2, 1],
    );
    let young_basis = find_young_basis(&s1_action, &s2_action);

    println!("Step 5. A Young-permutation quotient basis");
    println!("  Coordinates are in the quotient basis from Step 4.");
    for (name, vector) in &young_basis {
        println!("  {name} = {}", format_coordinate_vector(vector));
    }
    println!(
        "  span rank = {}",
        rank(
            &young_basis
                .iter()
                .map(|(_, v)| v.clone())
                .collect::<Vec<_>>()
        )
    );
    println!();

    println!("Step 6. Action on the Young basis");
    for (name, permutation, quotient_action) in module_actions {
        let action_in_young_basis = matrix_in_young_basis(&quotient_action, &young_basis);
        println!(
            "  representative {name}, cycle type {}:",
            cycle_type(&permutation)
        );
        println!("  action: {}", format_basis_action(&action_in_young_basis));
        println!("{}", format_matrix(&action_in_young_basis, 4));
        println!(
            "  trace = {}",
            format_q(&matrix_trace(&action_in_young_basis))
        );
    }
    println!();

    println!("Conclusion:");
    println!("  H^2 is one fixed vector plus a three-point permutation orbit.");
    println!("  Frobenius characteristic: h[3] + h[2,1] = 2*s[3] + s[2,1].");
}

fn class_representatives() -> Vec<(&'static str, Vec<usize>)> {
    vec![
        ("id", vec![0, 1, 2]),
        ("(12)", vec![1, 0, 2]),
        ("(123)", vec![1, 2, 0]),
    ]
}

fn gkm_module_basis(
    n: usize,
    hessenberg_edges: &[(usize, usize)],
    fixed_points: &[Vec<usize>],
    sn_index: &SnIndex,
    degree: u32,
) -> Vec<Vector<Q>> {
    let monomials = homogeneous_monomials(n, degree);
    let monomial_count = monomials.len();
    let ambient_dimension = fixed_points.len() * monomial_count;
    let mut constraints = Vec::new();

    for (point_index, point) in fixed_points.iter().enumerate() {
        for &(i, j) in hessenberg_edges {
            let adjacent_point_index = right_transposition_fixed_point_index(sn_index, point, i, j);
            if point_index > adjacent_point_index {
                continue;
            }

            let label_a = point[i];
            let label_b = point[j];
            for group in substitution_groups_for_labels(&monomials, label_a, label_b) {
                let mut row = vec![Q::zero(); ambient_dimension];
                for monomial_idx in group {
                    row[point_index * monomial_count + monomial_idx] =
                        row[point_index * monomial_count + monomial_idx].clone() + Q::one();
                    row[adjacent_point_index * monomial_count + monomial_idx] =
                        row[adjacent_point_index * monomial_count + monomial_idx].clone()
                            - Q::one();
                }
                constraints.push(row);
            }
        }
    }

    kernel_basis(&constraints)
}

fn module_dot_action_matrix(
    basis: &[Vector<Q>],
    fixed_points: &[Vec<usize>],
    sn_index: &SnIndex,
    monomials: &[Vec<u32>],
    permutation: &[usize],
) -> Matrix<Q> {
    let dim = basis.len();
    let mut matrix = vec![vec![Q::zero(); dim]; dim];
    for (col, basis_vector) in basis.iter().enumerate() {
        let image = apply_dot_action_to_ambient_vector(
            fixed_points,
            sn_index,
            monomials,
            permutation,
            basis_vector,
        );
        let coords = coordinates_in_basis(basis, &image);
        for row in 0..dim {
            matrix[row][col] = coords[row].clone();
        }
    }
    matrix
}

fn quotient_action_for_permutation(
    basis: &[Vector<Q>],
    ordinary_quotient: &QuotientSpace<Q>,
    fixed_points: &[Vec<usize>],
    sn_index: &SnIndex,
    monomials: &[Vec<u32>],
    permutation: &[usize],
) -> Matrix<Q> {
    let module_action =
        module_dot_action_matrix(basis, fixed_points, sn_index, monomials, permutation);
    quotient_action_matrix(ordinary_quotient, &module_action)
}

fn apply_dot_action_to_ambient_vector(
    fixed_points: &[Vec<usize>],
    sn_index: &SnIndex,
    monomials: &[Vec<u32>],
    permutation: &[usize],
    vector: &[Q],
) -> Vector<Q> {
    let monomial_index = monomials
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, monomial)| (monomial, index))
        .collect::<BTreeMap<_, _>>();
    let monomial_count = monomials.len();
    let mut result = vec![Q::zero(); vector.len()];

    for (source_point_index, point) in fixed_points.iter().enumerate() {
        let target_point_index = sn_index.rank_left_composition(permutation, point);
        for (source_monomial_index, monomial) in monomials.iter().enumerate() {
            let source_col = source_point_index * monomial_count + source_monomial_index;
            let coeff = &vector[source_col];
            if coeff.is_zero() {
                continue;
            }

            let image_monomial = permute_monomial_variables(monomial, permutation);
            let target_monomial_index = monomial_index[&image_monomial];
            let target_col = target_point_index * monomial_count + target_monomial_index;
            result[target_col] = result[target_col].clone() + coeff.clone();
        }
    }

    result
}

fn global_variable_assignment(
    variable: usize,
    fixed_points: &[Vec<usize>],
    monomials: &[Vec<u32>],
    n: usize,
) -> Vector<Q> {
    let monomial_count = monomials.len();
    let mut result = vec![Q::zero(); fixed_points.len() * monomial_count];
    let mut exponent = vec![0u32; n];
    exponent[variable] = 1;
    let monomial_index = monomials
        .iter()
        .position(|monomial| *monomial == exponent)
        .expect("degree-one monomial missing");
    for point_index in 0..fixed_points.len() {
        result[point_index * monomial_count + monomial_index] = Q::one();
    }
    result
}

fn find_young_basis(
    s1_action: &Matrix<Q>,
    s2_action: &Matrix<Q>,
) -> Vec<(&'static str, Vector<Q>)> {
    let dimension = s1_action.len();
    let candidates = small_integer_vectors(dimension, 1);
    let a1 = candidates
        .iter()
        .find_map(|candidate| {
            if matrix_vector_multiply(s2_action, candidate) != *candidate {
                return None;
            }
            let a2 = matrix_vector_multiply(s1_action, candidate);
            let a3 = matrix_vector_multiply(s2_action, &a2);
            if rank(&[candidate.clone(), a2, a3]) == 3 {
                Some(candidate.clone())
            } else {
                None
            }
        })
        .expect("expected a vector fixed by s2 with a three-point orbit");

    let a2 = matrix_vector_multiply(s1_action, &a1);
    let a3 = matrix_vector_multiply(s2_action, &a2);
    let c = small_integer_vectors(dimension, 2)
        .into_iter()
        .find(|candidate| {
            matrix_vector_multiply(s1_action, candidate) == *candidate
                && matrix_vector_multiply(s2_action, candidate) == *candidate
                && rank(&[candidate.clone(), a1.clone(), a2.clone(), a3.clone()]) == 4
        })
        .expect("expected an invariant vector outside the three-point orbit span");

    vec![("c", c), ("a1", a1), ("a2", a2), ("a3", a3)]
}

fn matrix_in_young_basis(
    action: &Matrix<Q>,
    young_basis: &[(&'static str, Vector<Q>)],
) -> Matrix<Q> {
    let basis = young_basis
        .iter()
        .map(|(_, vector)| vector.clone())
        .collect::<Vec<_>>();
    let mut result = vec![vec![Q::zero(); basis.len()]; basis.len()];
    for (col, vector) in basis.iter().enumerate() {
        let image = matrix_vector_multiply(action, vector);
        let coords = coordinates_in_basis(&basis, &image);
        for row in 0..basis.len() {
            result[row][col] = coords[row].clone();
        }
    }
    result
}

fn coordinates_in_basis(basis: &[Vector<Q>], vector: &[Q]) -> Vector<Q> {
    if basis.is_empty() {
        assert!(vector.iter().all(Ring::is_zero));
        return Vec::new();
    }

    let ambient_dimension = vector.len();
    let matrix = (0..ambient_dimension)
        .map(|row| {
            basis
                .iter()
                .map(|basis_vector| basis_vector[row].clone())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    solve_linear_system(&matrix, vector).expect("vector should lie in the basis span")
}

fn small_integer_vectors(dimension: usize, max_abs_value: i64) -> Vec<Vector<Q>> {
    fn rec(index: usize, current: &mut [Q], max_abs_value: i64, result: &mut Vec<Vector<Q>>) {
        if index == current.len() {
            if current.iter().any(|entry| !entry.is_zero()) {
                result.push(current.to_vec());
            }
            return;
        }

        for value in -max_abs_value..=max_abs_value {
            current[index] = Q::from_integer(value);
            rec(index + 1, current, max_abs_value, result);
        }
        current[index] = Q::zero();
    }

    let mut result = Vec::new();
    let mut current = vec![Q::zero(); dimension];
    rec(0, &mut current, max_abs_value, &mut result);
    result
}

fn right_transposition_fixed_point_index(
    sn_index: &SnIndex,
    point: &[usize],
    i: usize,
    j: usize,
) -> usize {
    let mut adjacent_point = point.to_vec();
    adjacent_point.swap(i, j);
    sn_index.rank(&adjacent_point)
}

fn homogeneous_monomials(num_vars: usize, degree: u32) -> Vec<Vec<u32>> {
    fn rec(remaining: u32, index: usize, current: &mut [u32], result: &mut Vec<Vec<u32>>) {
        if index + 1 == current.len() {
            current[index] = remaining;
            result.push(current.to_vec());
            current[index] = 0;
            return;
        }

        for value in (0..=remaining).rev() {
            current[index] = value;
            rec(remaining - value, index + 1, current, result);
        }
        current[index] = 0;
    }

    let mut result = Vec::new();
    let mut current = vec![0u32; num_vars];
    rec(degree, 0, &mut current, &mut result);
    result
}

fn substitution_groups_for_labels(
    monomials: &[Vec<u32>],
    label_a: usize,
    label_b: usize,
) -> Vec<Vec<usize>> {
    let mut groups_by_key: BTreeMap<Vec<u32>, Vec<usize>> = BTreeMap::new();
    for (monomial_idx, monomial) in monomials.iter().enumerate() {
        groups_by_key
            .entry(substitute_equal_variables(monomial, label_a, label_b))
            .or_default()
            .push(monomial_idx);
    }
    groups_by_key.into_values().collect()
}

fn substitute_equal_variables(monomial: &[u32], label_a: usize, label_b: usize) -> Vec<u32> {
    let mut key = monomial.to_vec();
    key[label_a] += key[label_b];
    key[label_b] = 0;
    key
}

fn permute_monomial_variables(monomial: &[u32], permutation: &[usize]) -> Vec<u32> {
    let mut result = vec![0u32; monomial.len()];
    for (source, &target) in permutation.iter().enumerate() {
        result[target] = monomial[source];
    }
    result
}

fn format_gkm_assignment(
    vector: &[Q],
    fixed_points: &[Vec<usize>],
    monomials: &[Vec<u32>],
) -> String {
    let monomial_count = monomials.len();
    fixed_points
        .iter()
        .enumerate()
        .map(|(point_index, point)| {
            let coeffs = &vector[point_index * monomial_count..(point_index + 1) * monomial_count];
            format!(
                "{} -> {}",
                format_permutation(point),
                format_polynomial(coeffs, monomials)
            )
        })
        .collect::<Vec<_>>()
        .join("; ")
}

fn format_polynomial(coeffs: &[Q], monomials: &[Vec<u32>]) -> String {
    let mut terms = Vec::new();
    for (coeff, monomial) in coeffs.iter().zip(monomials.iter()) {
        if coeff.is_zero() {
            continue;
        }
        let monomial_text = format_monomial(monomial);
        terms.push(format_signed_term(coeff, &monomial_text));
    }

    if terms.is_empty() {
        "0".to_string()
    } else {
        join_signed_terms(&terms)
    }
}

fn format_monomial(monomial: &[u32]) -> String {
    let mut factors = Vec::new();
    for (index, &power) in monomial.iter().enumerate() {
        if power == 0 {
            continue;
        }
        if power == 1 {
            factors.push(format!("t{}", index + 1));
        } else {
            factors.push(format!("t{}^{}", index + 1, power));
        }
    }
    if factors.is_empty() {
        "1".to_string()
    } else {
        factors.join("*")
    }
}

fn format_basis_combination(coords: &[Q], prefix: &str) -> String {
    let mut terms = Vec::new();
    for (index, coeff) in coords.iter().enumerate() {
        if coeff.is_zero() {
            continue;
        }
        terms.push(format_signed_term(coeff, &format!("{prefix}{index}")));
    }
    if terms.is_empty() {
        "0".to_string()
    } else {
        join_signed_terms(&terms)
    }
}

fn format_coordinate_vector(vector: &[Q]) -> String {
    format!(
        "[{}]",
        vector.iter().map(format_q).collect::<Vec<_>>().join(", ")
    )
}

fn format_matrix(matrix: &Matrix<Q>, indent: usize) -> String {
    let spaces = " ".repeat(indent);
    matrix
        .iter()
        .map(|row| format!("{spaces}{}", format_coordinate_vector(row)))
        .collect::<Vec<_>>()
        .join("\n")
}

fn format_basis_action(matrix: &Matrix<Q>) -> String {
    let names = ["c", "a1", "a2", "a3"];
    names
        .iter()
        .enumerate()
        .map(|(col, source)| {
            let image = matrix
                .iter()
                .map(|row| row[col].clone())
                .collect::<Vec<_>>();
            format!("{source} -> {}", format_named_combination(&image, &names))
        })
        .collect::<Vec<_>>()
        .join(", ")
}

fn format_named_combination(coords: &[Q], names: &[&str]) -> String {
    let mut terms = Vec::new();
    for (coeff, name) in coords.iter().zip(names.iter()) {
        if coeff.is_zero() {
            continue;
        }
        terms.push(format_signed_term(coeff, name));
    }
    if terms.is_empty() {
        "0".to_string()
    } else {
        join_signed_terms(&terms)
    }
}

fn format_signed_term(coeff: &Q, term: &str) -> String {
    if term == "1" {
        return format_q(coeff);
    }
    if *coeff == Q::one() {
        term.to_string()
    } else if *coeff == -Q::one() {
        format!("-{term}")
    } else {
        format!("{}*{term}", format_q(coeff))
    }
}

fn join_signed_terms(terms: &[String]) -> String {
    let mut result = String::new();
    for term in terms {
        if result.is_empty() {
            result.push_str(term);
        } else if let Some(stripped) = term.strip_prefix('-') {
            result.push_str(" - ");
            result.push_str(stripped);
        } else {
            result.push_str(" + ");
            result.push_str(term);
        }
    }
    result
}

fn format_permutation(permutation: &[usize]) -> String {
    permutation
        .iter()
        .map(|value| (value + 1).to_string())
        .collect::<Vec<_>>()
        .join("")
}

fn format_q(value: &Q) -> String {
    if *value.denom() == 1 {
        value.numer().to_string()
    } else {
        format!("{}/{}", value.numer(), value.denom())
    }
}
