//! Small GKM-model computations for Hessenberg dot actions.
//!
//! This module is intended for exact small-rank experiments. It computes the
//! homogeneous GKM functions on the regular semisimple Hessenberg moment graph,
//! quotients degree-by-degree by the positive polynomial variables, and
//! descends Tymoczko's dot action to ordinary cohomology.

use std::collections::BTreeMap;

use combinatoric_core::Graph;
use num_rational::Ratio;
use sym_poly_core::linear_algebra::{
    quotient_action_matrix, rref, zero_matrix, Matrix, QuotientSpace, Vector,
};
use sym_poly_core::sn_action::{
    assert_permutation, compose_permutations, conjugacy_class_representatives, inverse_permutation,
};
use sym_poly_core::{Partition, Ring};

use crate::frobenius::graded_frobenius_from_trace_matrices;
use crate::SymmetricFunction;

type Q = Ratio<i64>;

#[derive(Debug, Clone)]
struct HomogeneousGkmComponent {
    fixed_points: Vec<Vec<usize>>,
    fixed_point_index: BTreeMap<Vec<usize>, usize>,
    monomials: Vec<Vec<u32>>,
    monomial_index: BTreeMap<Vec<u32>, usize>,
    module_basis: Vec<Vector<Q>>,
    module_coordinate_columns: Vec<usize>,
    ordinary_quotient: QuotientSpace<Q>,
}

#[derive(Debug, Clone)]
struct GkmModuleBasis {
    vectors: Vec<Vector<Q>>,
    coordinate_columns: Vec<usize>,
}

/// Compute dot-action matrices on ordinary Hessenberg cohomology.
///
/// The input is a Dyck area sequence. The corresponding unit-interval graph is
/// used as the Hessenberg incomparability graph. The output is keyed first by
/// polynomial degree and then by conjugacy class. Each matrix is the ordinary
/// cohomology dot action on that degree.
pub fn hessenberg_gkm_dot_action_matrices(
    area: &[u8],
) -> Option<BTreeMap<u32, BTreeMap<Partition, Matrix<Ratio<i64>>>>> {
    if !is_area_sequence(area) {
        return None;
    }

    let graph = Graph::unit_interval(area);
    let n = area.len();
    let max_degree = graph.num_edges() as u32;
    let fixed_points = sym_poly_core::symmetric_group_permutation_basis(n);
    let fixed_point_index = fixed_points
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, point)| (point, index))
        .collect::<BTreeMap<_, _>>();

    let mut components = Vec::new();
    for degree in 0..=max_degree {
        let module_basis = homogeneous_gkm_module_basis(n, graph.edges(), degree);
        let mut component = HomogeneousGkmComponent::new(
            fixed_points.clone(),
            fixed_point_index.clone(),
            degree,
            module_basis.vectors,
            module_basis.coordinate_columns,
            Vec::new(),
        );
        let relations = if degree == 0 {
            Vec::new()
        } else {
            variable_multiple_relations(&components[(degree - 1) as usize], &component, n)
        };
        component.ordinary_quotient =
            QuotientSpace::from_relations(component.module_basis.len(), &relations);
        components.push(component);
    }

    let mut by_degree = BTreeMap::new();
    for (degree, component) in components.iter().enumerate() {
        if component.ordinary_quotient.dimension() == 0 {
            continue;
        }

        let mut class_matrices = BTreeMap::new();
        for (cycle_type, representative) in conjugacy_class_representatives(n) {
            let module_action = component.module_dot_action_matrix(&representative);
            let ordinary_action =
                quotient_action_matrix(&component.ordinary_quotient, &module_action);
            class_matrices.insert(cycle_type, ordinary_action);
        }
        by_degree.insert(degree as u32, class_matrices);
    }

    Some(by_degree)
}

/// Compute the graded Frobenius characteristic of the GKM dot action.
pub fn hessenberg_gkm_dot_frobenius(
    area: &[u8],
) -> Option<BTreeMap<u32, SymmetricFunction<Ratio<i64>>>> {
    let matrices = hessenberg_gkm_dot_action_matrices(area)?;
    Some(
        graded_frobenius_from_trace_matrices(&matrices)
            .into_iter()
            .filter(|(_, f)| !f.is_zero())
            .collect(),
    )
}

impl HomogeneousGkmComponent {
    fn new(
        fixed_points: Vec<Vec<usize>>,
        fixed_point_index: BTreeMap<Vec<usize>, usize>,
        degree: u32,
        module_basis: Vec<Vector<Q>>,
        module_coordinate_columns: Vec<usize>,
        relations: Vec<Vector<Q>>,
    ) -> Self {
        let monomials = homogeneous_monomials(fixed_points.first().map_or(0, Vec::len), degree);
        let monomial_index = monomials
            .iter()
            .cloned()
            .enumerate()
            .map(|(index, monomial)| (monomial, index))
            .collect::<BTreeMap<_, _>>();
        let ordinary_quotient = QuotientSpace::from_relations(module_basis.len(), &relations);

        Self {
            fixed_points,
            fixed_point_index,
            monomials,
            monomial_index,
            module_basis,
            module_coordinate_columns,
            ordinary_quotient,
        }
    }

    fn ambient_dimension(&self) -> usize {
        self.fixed_points.len() * self.monomials.len()
    }

    fn module_coordinates(&self, vector: &[Q]) -> Vector<Q> {
        assert_eq!(
            vector.len(),
            self.ambient_dimension(),
            "ambient vector has wrong dimension"
        );
        self.module_coordinate_columns
            .iter()
            .map(|&col| vector[col].clone())
            .collect()
    }

    fn module_dot_action_matrix(&self, permutation: &[usize]) -> Matrix<Q> {
        assert_permutation(permutation);
        let dim = self.module_basis.len();
        let mut matrix = zero_matrix::<Q>(dim, dim);

        for (col, basis_vector) in self.module_basis.iter().enumerate() {
            let image = self.apply_dot_action_to_ambient_vector(permutation, basis_vector);
            let coords = self.module_coordinates(&image);
            for row in 0..dim {
                matrix[row][col] = coords[row];
            }
        }

        matrix
    }

    fn apply_dot_action_to_ambient_vector(&self, permutation: &[usize], vector: &[Q]) -> Vector<Q> {
        let mut result = vec![Q::zero(); self.ambient_dimension()];
        let monomial_count = self.monomials.len();
        let inverse = inverse_permutation(permutation);

        for (target_point_index, target_point) in self.fixed_points.iter().enumerate() {
            let source_point = compose_permutations(&inverse, target_point);
            let source_point_index = self.fixed_point_index[&source_point];

            for (source_monomial_index, monomial) in self.monomials.iter().enumerate() {
                let source_col = source_point_index * monomial_count + source_monomial_index;
                let coeff = &vector[source_col];
                if coeff.is_zero() {
                    continue;
                }

                let image_monomial = permute_monomial_variables(monomial, permutation);
                let target_monomial_index = self.monomial_index[&image_monomial];
                let target_row = target_point_index * monomial_count + target_monomial_index;
                result[target_row] = result[target_row].clone() + coeff.clone();
            }
        }

        result
    }
}

fn homogeneous_gkm_module_basis(
    n: usize,
    hessenberg_edges: &[(usize, usize)],
    degree: u32,
) -> GkmModuleBasis {
    let fixed_points = sym_poly_core::symmetric_group_permutation_basis(n);
    let fixed_point_index = fixed_points
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, point)| (point, index))
        .collect::<BTreeMap<_, _>>();
    let monomials = homogeneous_monomials(n, degree);
    let ambient_dimension = fixed_points.len() * monomials.len();
    let mut constraints = Vec::new();

    for (point_index, point) in fixed_points.iter().enumerate() {
        for &(i, j) in hessenberg_edges {
            let transposition = transposition(n, i, j);
            let adjacent_point = compose_permutations(point, &transposition);
            let adjacent_point_index = fixed_point_index[&adjacent_point];
            if point_index > adjacent_point_index {
                continue;
            }

            let label_a = point[i];
            let label_b = point[j];
            for row in divisibility_constraints_for_edge(
                point_index,
                adjacent_point_index,
                label_a,
                label_b,
                &monomials,
                ambient_dimension,
            ) {
                constraints.push(row);
            }
        }
    }

    if constraints.is_empty() {
        GkmModuleBasis {
            vectors: standard_basis(ambient_dimension),
            coordinate_columns: (0..ambient_dimension).collect(),
        }
    } else {
        kernel_basis_with_coordinate_columns(ambient_dimension, &constraints)
    }
}

fn divisibility_constraints_for_edge(
    point_index: usize,
    adjacent_point_index: usize,
    label_a: usize,
    label_b: usize,
    monomials: &[Vec<u32>],
    ambient_dimension: usize,
) -> Vec<Vector<Q>> {
    let monomial_count = monomials.len();
    let mut rows_by_key: BTreeMap<Vec<u32>, Vector<Q>> = BTreeMap::new();

    for (monomial_idx, monomial) in monomials.iter().enumerate() {
        let key = substitute_equal_variables(monomial, label_a, label_b);
        let row = rows_by_key
            .entry(key)
            .or_insert_with(|| vec![Q::zero(); ambient_dimension]);
        row[point_index * monomial_count + monomial_idx] =
            row[point_index * monomial_count + monomial_idx].clone() + Q::one();
        row[adjacent_point_index * monomial_count + monomial_idx] =
            row[adjacent_point_index * monomial_count + monomial_idx].clone() - Q::one();
    }

    rows_by_key
        .into_values()
        .filter(|row| row.iter().any(|entry| !entry.is_zero()))
        .collect()
}

fn variable_multiple_relations(
    previous: &HomogeneousGkmComponent,
    current: &HomogeneousGkmComponent,
    n: usize,
) -> Vec<Vector<Q>> {
    let mut relations = Vec::new();
    let previous_monomial_count = previous.monomials.len();
    let current_monomial_count = current.monomials.len();

    for basis_vector in &previous.module_basis {
        for variable in 0..n {
            let mut ambient = vec![Q::zero(); current.ambient_dimension()];
            for point_index in 0..previous.fixed_points.len() {
                for (monomial_index, monomial) in previous.monomials.iter().enumerate() {
                    let source = point_index * previous_monomial_count + monomial_index;
                    let coeff = &basis_vector[source];
                    if coeff.is_zero() {
                        continue;
                    }

                    let mut product_monomial = monomial.clone();
                    product_monomial[variable] += 1;
                    let target_monomial_index = current.monomial_index[&product_monomial];
                    let target = point_index * current_monomial_count + target_monomial_index;
                    ambient[target] = ambient[target].clone() + coeff.clone();
                }
            }
            relations.push(current.module_coordinates(&ambient));
        }
    }

    relations
}

fn homogeneous_monomials(num_vars: usize, degree: u32) -> Vec<Vec<u32>> {
    if num_vars == 0 {
        return if degree == 0 {
            vec![Vec::new()]
        } else {
            Vec::new()
        };
    }
    let mut result = Vec::new();
    let mut current = vec![0u32; num_vars];
    weak_compositions_rec(degree, 0, &mut current, &mut result);
    result
}

fn weak_compositions_rec(
    remaining: u32,
    index: usize,
    current: &mut [u32],
    result: &mut Vec<Vec<u32>>,
) {
    if index + 1 == current.len() {
        current[index] = remaining;
        result.push(current.to_vec());
        current[index] = 0;
        return;
    }

    for value in 0..=remaining {
        current[index] = value;
        weak_compositions_rec(remaining - value, index + 1, current, result);
    }
    current[index] = 0;
}

fn standard_basis(dimension: usize) -> Vec<Vector<Q>> {
    (0..dimension)
        .map(|index| {
            let mut vector = vec![Q::zero(); dimension];
            vector[index] = Q::one();
            vector
        })
        .collect()
}

fn kernel_basis_with_coordinate_columns(
    ambient_dimension: usize,
    constraints: &[Vector<Q>],
) -> GkmModuleBasis {
    let reduced = rref(constraints);
    let coordinate_columns = complement_columns(ambient_dimension, &reduced.pivot_columns);
    let mut vectors = Vec::with_capacity(coordinate_columns.len());

    for &free_col in &coordinate_columns {
        let mut vector = vec![Q::zero(); ambient_dimension];
        vector[free_col] = Q::one();
        for (pivot_row, &pivot_col) in reduced.pivot_columns.iter().enumerate() {
            vector[pivot_col] = -reduced.matrix[pivot_row][free_col].clone();
        }
        vectors.push(vector);
    }

    GkmModuleBasis {
        vectors,
        coordinate_columns,
    }
}

fn complement_columns(num_cols: usize, pivot_columns: &[usize]) -> Vec<usize> {
    let mut is_pivot = vec![false; num_cols];
    for &col in pivot_columns {
        if col < num_cols {
            is_pivot[col] = true;
        }
    }
    (0..num_cols).filter(|&col| !is_pivot[col]).collect()
}

fn transposition(n: usize, i: usize, j: usize) -> Vec<usize> {
    assert!(i < n && j < n, "transposition index out of range");
    let mut permutation: Vec<_> = (0..n).collect();
    permutation.swap(i, j);
    permutation
}

fn substitute_equal_variables(monomial: &[u32], label_a: usize, label_b: usize) -> Vec<u32> {
    let mut key = monomial.to_vec();
    key[label_a] += key[label_b];
    key[label_b] = 0;
    key
}

fn permute_monomial_variables(monomial: &[u32], permutation: &[usize]) -> Vec<u32> {
    assert_eq!(
        monomial.len(),
        permutation.len(),
        "monomial and permutation have different sizes"
    );
    let mut result = vec![0u32; monomial.len()];
    for (source, &target) in permutation.iter().enumerate() {
        result[target] = monomial[source];
    }
    result
}

fn is_area_sequence(area: &[u8]) -> bool {
    area.iter().enumerate().all(|(i, &v)| v as usize <= i)
        && area
            .windows(2)
            .all(|w| usize::from(w[1]) <= usize::from(w[0]) + 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chromatic::hessenberg_area_dot_frobenius_target;
    use std::collections::BTreeSet;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn assert_matches_shareshian_wachs_target(area: &[u8]) {
        let computed = hessenberg_gkm_dot_frobenius(area).unwrap();
        let target = hessenberg_area_dot_frobenius_target(area).unwrap();

        assert_eq!(
            computed.keys().copied().collect::<Vec<_>>(),
            target.keys().copied().collect::<Vec<_>>()
        );

        for (&degree, computed_function) in &computed {
            let computed_schur = computed_function.to_schur_basis();
            let target_schur = target[&degree].to_schur_basis();
            let partitions = computed_schur
                .terms()
                .keys()
                .chain(target_schur.terms().keys())
                .cloned()
                .collect::<BTreeSet<_>>();

            for partition in partitions {
                assert_eq!(
                    computed_schur.coefficient(&partition),
                    q(target_schur.coefficient(&partition)),
                    "mismatch in degree {degree}, partition {partition:?}"
                );
            }
        }
    }

    #[test]
    fn test_hessenberg_gkm_dot_frobenius_rejects_invalid_area() {
        assert!(hessenberg_gkm_dot_frobenius(&[0, 2]).is_none());
    }

    #[test]
    fn test_hessenberg_gkm_dot_frobenius_edgeless_s3() {
        assert_matches_shareshian_wachs_target(&[0, 0, 0]);
    }

    #[test]
    fn test_hessenberg_gkm_dot_frobenius_complete_graph_s3() {
        assert_matches_shareshian_wachs_target(&[0, 1, 2]);
    }

    #[test]
    fn test_hessenberg_gkm_dot_frobenius_path_s3() {
        assert_matches_shareshian_wachs_target(&[0, 1, 1]);
    }
}
