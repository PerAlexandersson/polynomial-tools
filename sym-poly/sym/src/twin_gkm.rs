//! Small GKM-model computations for twin-manifold dagger actions.
//!
//! This module is intended for exact small-rank experiments. It computes the
//! homogeneous GKM functions for the twin-manifold graph attached to a
//! unit-interval area sequence, quotients degree-by-degree by the positive
//! polynomial variables, and descends the dagger action on fixed points.
//!
//! The convention differs from the regular semisimple Hessenberg GKM model in
//! [`crate::hessenberg_gkm`]: for an edge `(i, j)` the congruence label is the
//! fixed root `t_i - t_j`, while the action sends `(p_v)_v` to
//! `(p_{\sigma^{-1} v})_v` without permuting polynomial variables.

use std::collections::BTreeMap;

use combinatoric_core::Graph;
use num_rational::Ratio;
use sym_poly_core::linear_algebra::{
    matrix_trace, quotient_action_matrix, rref, zero_matrix, Matrix, QuotientSpace, Vector,
};
use sym_poly_core::sn_action::{
    assert_permutation, compose_permutations, conjugacy_class_representatives, inverse_permutation,
};
use sym_poly_core::{Partition, Ring};
use sym_poly_multipoly::{
    elementary_symmetric_generators, quotient_action_matrices_by_multidegree_and_cycle_type,
    quotient_basis, GroebnerBasis, IndexedVariables, MonomialOrder,
};

use crate::frobenius::graded_frobenius_from_character_values;
use crate::SymmetricFunction;

type Q = Ratio<i64>;

#[derive(Debug, Clone)]
struct HomogeneousTwinGkmComponent {
    fixed_points: Vec<Vec<usize>>,
    fixed_point_index: BTreeMap<Vec<usize>, usize>,
    monomials: Vec<Vec<u32>>,
    module_basis: Vec<Vector<Q>>,
    module_coordinate_columns: Vec<usize>,
    ordinary_quotient: QuotientSpace<Q>,
}

#[derive(Debug, Clone)]
struct TwinGkmCombinatorics {
    fixed_points: Vec<Vec<usize>>,
    edge_adjacencies: Vec<TwinEdgeAdjacency>,
    monomials_by_degree: Vec<Vec<Vec<u32>>>,
    substitution_groups_by_degree: Vec<Vec<Vec<Vec<usize>>>>,
    multiplication_maps_by_degree: Vec<Vec<Vec<usize>>>,
}

#[derive(Debug, Clone)]
struct TwinEdgeAdjacency {
    edge_index: usize,
    point_index: usize,
    adjacent_point_index: usize,
}

#[derive(Debug, Clone)]
struct TwinGkmModuleBasis {
    vectors: Vec<Vector<Q>>,
    coordinate_columns: Vec<usize>,
}

/// Compute dagger-action matrices on ordinary twin-manifold cohomology.
///
/// The input is a Dyck area sequence. The corresponding unit-interval graph
/// supplies the twin-manifold GKM edges. The output is keyed first by
/// polynomial degree and then by conjugacy class. Each matrix is the ordinary
/// cohomology dagger action on that degree. For the complete graph, the same
/// action is computed from the Artin coinvariant-ring presentation
/// `Q[x_1,...,x_n]/<e_1,...,e_n>` rather than by dense GKM row reduction.
pub fn twin_gkm_dagger_action_matrices(
    area: &[u8],
) -> Option<BTreeMap<u32, BTreeMap<Partition, Matrix<Ratio<i64>>>>> {
    if !is_area_sequence(area) {
        return None;
    }
    if is_complete_area_sequence(area) {
        return Some(complete_graph_artin_action_matrices(area.len()));
    }

    let n = area.len();
    let components = noncomplete_twin_gkm_components(area);

    let mut by_degree = BTreeMap::new();
    for (degree, component) in components.iter().enumerate() {
        if component.ordinary_quotient.dimension() == 0 {
            continue;
        }

        let mut class_matrices = BTreeMap::new();
        for (cycle_type, representative) in conjugacy_class_representatives(n) {
            let module_action = component.module_dagger_action_matrix(&representative);
            let ordinary_action =
                quotient_action_matrix(&component.ordinary_quotient, &module_action);
            class_matrices.insert(cycle_type, ordinary_action);
        }
        by_degree.insert(degree as u32, class_matrices);
    }

    Some(by_degree)
}

/// Compute graded character values of the twin-manifold dagger action.
pub fn twin_gkm_dagger_character_values_by_degree(
    area: &[u8],
) -> Option<BTreeMap<u32, BTreeMap<Partition, Ratio<i64>>>> {
    if !is_area_sequence(area) {
        return None;
    }
    if is_complete_area_sequence(area) {
        let matrices = complete_graph_artin_action_matrices(area.len());
        return Some(trace_values_from_action_matrices(matrices));
    }

    let n = area.len();
    let components = noncomplete_twin_gkm_components(area);
    let mut by_degree = BTreeMap::new();

    for (degree, component) in components.iter().enumerate() {
        if component.ordinary_quotient.dimension() == 0 {
            continue;
        }

        let mut values = BTreeMap::new();
        for (cycle_type, representative) in conjugacy_class_representatives(n) {
            let trace = component.ordinary_dagger_trace(&representative);
            if !trace.is_zero() {
                values.insert(cycle_type, trace);
            }
        }
        by_degree.insert(degree as u32, values);
    }

    Some(by_degree)
}

/// Compute the graded Frobenius characteristic of the twin GKM dagger action.
pub fn twin_gkm_dagger_frobenius(
    area: &[u8],
) -> Option<BTreeMap<u32, SymmetricFunction<Ratio<i64>>>> {
    let character_values = twin_gkm_dagger_character_values_by_degree(area)?;
    Some(
        graded_frobenius_from_character_values(&character_values)
            .into_iter()
            .filter(|(_, f)| !f.is_zero())
            .collect(),
    )
}

impl HomogeneousTwinGkmComponent {
    fn new(
        fixed_points: Vec<Vec<usize>>,
        fixed_point_index: BTreeMap<Vec<usize>, usize>,
        monomials: Vec<Vec<u32>>,
        module_basis: Vec<Vector<Q>>,
        module_coordinate_columns: Vec<usize>,
        relations: Vec<Vector<Q>>,
    ) -> Self {
        let ordinary_quotient = QuotientSpace::from_relations(module_basis.len(), &relations);

        Self {
            fixed_points,
            fixed_point_index,
            monomials,
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

    fn module_dagger_action_matrix(&self, permutation: &[usize]) -> Matrix<Q> {
        assert_permutation(permutation);
        let dim = self.module_basis.len();
        let mut matrix = zero_matrix::<Q>(dim, dim);
        let source_point_indices = self.dagger_source_point_indices(permutation);

        for (col, basis_vector) in self.module_basis.iter().enumerate() {
            let image =
                self.apply_dagger_action_to_ambient_vector(&source_point_indices, basis_vector);
            let coords = self.module_coordinates(&image);
            for row in 0..dim {
                matrix[row][col] = coords[row].clone();
            }
        }

        matrix
    }

    fn ordinary_dagger_trace(&self, permutation: &[usize]) -> Q {
        assert_permutation(permutation);
        let source_point_indices = self.dagger_source_point_indices(permutation);
        let mut trace = Q::zero();

        for (quotient_col, &module_basis_index) in
            self.ordinary_quotient.free_columns.iter().enumerate()
        {
            let image = self.apply_dagger_action_to_ambient_vector(
                &source_point_indices,
                &self.module_basis[module_basis_index],
            );
            let module_coords = self.module_coordinates(&image);
            let quotient_coords = self.ordinary_quotient.quotient_coordinates(&module_coords);
            trace = trace + quotient_coords[quotient_col].clone();
        }

        trace
    }

    fn dagger_source_point_indices(&self, permutation: &[usize]) -> Vec<usize> {
        let inverse = inverse_permutation(permutation);
        self.fixed_points
            .iter()
            .map(|target_point| {
                let source_point = compose_permutations(&inverse, target_point);
                self.fixed_point_index[&source_point]
            })
            .collect()
    }

    fn apply_dagger_action_to_ambient_vector(
        &self,
        source_point_indices: &[usize],
        vector: &[Q],
    ) -> Vector<Q> {
        let mut result = vec![Q::zero(); self.ambient_dimension()];
        let monomial_count = self.monomials.len();

        for (target_point_index, &source_point_index) in source_point_indices.iter().enumerate() {
            for monomial_index in 0..monomial_count {
                let source_col = source_point_index * monomial_count + monomial_index;
                let coeff = &vector[source_col];
                if coeff.is_zero() {
                    continue;
                }

                let target_row = target_point_index * monomial_count + monomial_index;
                result[target_row] = result[target_row].clone() + coeff.clone();
            }
        }

        result
    }
}

fn noncomplete_twin_gkm_components(area: &[u8]) -> Vec<HomogeneousTwinGkmComponent> {
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
    let combinatorics = TwinGkmCombinatorics::new(
        n,
        graph.edges(),
        fixed_points.clone(),
        fixed_point_index.clone(),
        max_degree,
    );

    let mut components = Vec::new();
    for degree in 0..=max_degree {
        let module_basis = homogeneous_twin_gkm_module_basis(&combinatorics, degree);
        let mut component = HomogeneousTwinGkmComponent::new(
            fixed_points.clone(),
            fixed_point_index.clone(),
            combinatorics.monomials(degree).to_vec(),
            module_basis.vectors,
            module_basis.coordinate_columns,
            Vec::new(),
        );
        let relations = if degree == 0 {
            Vec::new()
        } else {
            variable_multiple_relations(
                &components[(degree - 1) as usize],
                &component,
                combinatorics.multiplication_maps(degree),
            )
        };
        component.ordinary_quotient =
            QuotientSpace::from_relations(component.module_basis.len(), &relations);
        components.push(component);
    }

    components
}

fn trace_values_from_action_matrices(
    matrices: BTreeMap<u32, BTreeMap<Partition, Matrix<Q>>>,
) -> BTreeMap<u32, BTreeMap<Partition, Q>> {
    matrices
        .into_iter()
        .map(|(degree, class_matrices)| {
            let values = class_matrices
                .into_iter()
                .map(|(cycle_type, matrix)| (cycle_type, matrix_trace(&matrix)))
                .filter(|(_, trace)| !trace.is_zero())
                .collect();
            (degree, values)
        })
        .collect()
}

impl TwinGkmCombinatorics {
    fn new(
        n: usize,
        twin_edges: &[(usize, usize)],
        fixed_points: Vec<Vec<usize>>,
        fixed_point_index: BTreeMap<Vec<usize>, usize>,
        max_degree: u32,
    ) -> Self {
        let edge_adjacencies =
            twin_edge_adjacencies(n, twin_edges, &fixed_points, &fixed_point_index);
        let monomials_by_degree = (0..=max_degree)
            .map(|degree| homogeneous_monomials(n, degree))
            .collect::<Vec<_>>();
        let monomial_index_by_degree = monomials_by_degree
            .iter()
            .map(|monomials| {
                monomials
                    .iter()
                    .cloned()
                    .enumerate()
                    .map(|(index, monomial)| (monomial, index))
                    .collect::<BTreeMap<_, _>>()
            })
            .collect::<Vec<_>>();
        let substitution_groups_by_degree = monomials_by_degree
            .iter()
            .map(|monomials| {
                twin_edges
                    .iter()
                    .map(|&(i, j)| substitution_groups_for_edge(monomials, i, j))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut multiplication_maps_by_degree = Vec::with_capacity(max_degree as usize + 1);
        multiplication_maps_by_degree.push(Vec::new());
        for degree in 1..=max_degree {
            multiplication_maps_by_degree.push(multiplication_maps_between_degrees(
                &monomials_by_degree[(degree - 1) as usize],
                &monomial_index_by_degree[degree as usize],
                n,
            ));
        }

        Self {
            fixed_points,
            edge_adjacencies,
            monomials_by_degree,
            substitution_groups_by_degree,
            multiplication_maps_by_degree,
        }
    }

    fn monomials(&self, degree: u32) -> &[Vec<u32>] {
        &self.monomials_by_degree[degree as usize]
    }

    fn substitution_groups(&self, degree: u32, edge_index: usize) -> &[Vec<usize>] {
        &self.substitution_groups_by_degree[degree as usize][edge_index]
    }

    fn multiplication_maps(&self, degree: u32) -> &[Vec<usize>] {
        &self.multiplication_maps_by_degree[degree as usize]
    }
}

fn homogeneous_twin_gkm_module_basis(
    combinatorics: &TwinGkmCombinatorics,
    degree: u32,
) -> TwinGkmModuleBasis {
    let monomials = combinatorics.monomials(degree);
    let ambient_dimension = combinatorics.fixed_points.len() * monomials.len();
    let mut constraints = Vec::new();

    for adjacency in &combinatorics.edge_adjacencies {
        for row in divisibility_constraints_for_edge(
            adjacency.point_index,
            adjacency.adjacent_point_index,
            combinatorics.substitution_groups(degree, adjacency.edge_index),
            monomials.len(),
            ambient_dimension,
        ) {
            constraints.push(row);
        }
    }

    if constraints.is_empty() {
        TwinGkmModuleBasis {
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
    substitution_groups: &[Vec<usize>],
    monomial_count: usize,
    ambient_dimension: usize,
) -> Vec<Vector<Q>> {
    substitution_groups
        .iter()
        .map(|group| {
            let mut row = vec![Q::zero(); ambient_dimension];
            for &monomial_idx in group {
                row[point_index * monomial_count + monomial_idx] =
                    row[point_index * monomial_count + monomial_idx].clone() + Q::one();
                row[adjacent_point_index * monomial_count + monomial_idx] =
                    row[adjacent_point_index * monomial_count + monomial_idx].clone() - Q::one();
            }
            row
        })
        .collect()
}

fn variable_multiple_relations(
    previous: &HomogeneousTwinGkmComponent,
    current: &HomogeneousTwinGkmComponent,
    multiplication_maps: &[Vec<usize>],
) -> Vec<Vector<Q>> {
    let mut relations = Vec::new();
    let previous_monomial_count = previous.monomials.len();
    let current_monomial_count = current.monomials.len();

    for basis_vector in &previous.module_basis {
        for variable_map in multiplication_maps {
            let mut ambient = vec![Q::zero(); current.ambient_dimension()];
            for point_index in 0..previous.fixed_points.len() {
                for (monomial_index, &target_monomial_index) in variable_map.iter().enumerate() {
                    let source = point_index * previous_monomial_count + monomial_index;
                    let coeff = &basis_vector[source];
                    if coeff.is_zero() {
                        continue;
                    }

                    let target = point_index * current_monomial_count + target_monomial_index;
                    ambient[target] = ambient[target].clone() + coeff.clone();
                }
            }
            relations.push(current.module_coordinates(&ambient));
        }
    }

    relations
}

fn twin_edge_adjacencies(
    n: usize,
    twin_edges: &[(usize, usize)],
    fixed_points: &[Vec<usize>],
    fixed_point_index: &BTreeMap<Vec<usize>, usize>,
) -> Vec<TwinEdgeAdjacency> {
    let mut result = Vec::new();
    for (point_index, point) in fixed_points.iter().enumerate() {
        for (edge_index, &(i, j)) in twin_edges.iter().enumerate() {
            let transposition = transposition(n, i, j);
            let adjacent_point = compose_permutations(point, &transposition);
            let adjacent_point_index = fixed_point_index[&adjacent_point];
            if point_index > adjacent_point_index {
                continue;
            }
            result.push(TwinEdgeAdjacency {
                edge_index,
                point_index,
                adjacent_point_index,
            });
        }
    }
    result
}

fn substitution_groups_for_edge(
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

fn multiplication_maps_between_degrees(
    previous_monomials: &[Vec<u32>],
    current_monomial_index: &BTreeMap<Vec<u32>, usize>,
    n: usize,
) -> Vec<Vec<usize>> {
    let mut maps = Vec::with_capacity(n);
    for variable in 0..n {
        let mut variable_map = Vec::with_capacity(previous_monomials.len());
        for monomial in previous_monomials {
            let mut product_monomial = monomial.clone();
            product_monomial[variable] += 1;
            variable_map.push(current_monomial_index[&product_monomial]);
        }
        maps.push(variable_map);
    }
    maps
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
) -> TwinGkmModuleBasis {
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

    TwinGkmModuleBasis {
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

fn is_area_sequence(area: &[u8]) -> bool {
    area.iter().enumerate().all(|(i, &v)| v as usize <= i)
        && area
            .windows(2)
            .all(|w| usize::from(w[1]) <= usize::from(w[0]) + 1)
}

fn is_complete_area_sequence(area: &[u8]) -> bool {
    area.iter()
        .enumerate()
        .all(|(index, &value)| value as usize == index)
}

fn complete_graph_artin_action_matrices(n: usize) -> BTreeMap<u32, BTreeMap<Partition, Matrix<Q>>> {
    let groebner_basis =
        GroebnerBasis::new(elementary_symmetric_generators::<Q>(n), MonomialOrder::Lex);
    let basis = quotient_basis(&groebner_basis).expect("Artin quotient has a finite basis");
    let variables = IndexedVariables::new(1, n);
    quotient_action_matrices_by_multidegree_and_cycle_type(&variables, &groebner_basis, &basis)
        .expect("Artin quotient action preserves degree")
        .into_iter()
        .map(|(degree, matrices)| {
            assert_eq!(degree.len(), 1, "Artin quotient has one grading");
            (degree[0], matrices)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::frobenius::graded_frobenius_from_trace_matrices;
    use crate::llt::{
        unicellular_llt_character_values_by_degree, unicellular_llt_frobenius_target,
    };
    use std::collections::BTreeSet;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn assert_matches_llt_character_target(area: &[u8]) {
        let computed = twin_gkm_dagger_character_values_by_degree(area).unwrap();
        let target = unicellular_llt_character_values_by_degree(area).unwrap();

        assert_eq!(
            computed.keys().copied().collect::<Vec<_>>(),
            target.keys().copied().collect::<Vec<_>>()
        );

        for degree in computed.keys() {
            for cycle_type in Partition::all_of_size(area.len() as u32) {
                let computed_value = computed
                    .get(degree)
                    .and_then(|values| values.get(&cycle_type))
                    .cloned()
                    .unwrap_or_else(Q::zero);
                let target_value = target
                    .get(degree)
                    .and_then(|values| values.get(&cycle_type))
                    .copied()
                    .unwrap_or(0);
                assert_eq!(
                    computed_value,
                    q(target_value),
                    "mismatch in degree {degree}, cycle type {cycle_type:?}"
                );
            }
        }
    }

    fn assert_matches_llt_frobenius_target(area: &[u8]) {
        let computed = twin_gkm_dagger_frobenius(area).unwrap();
        let target = unicellular_llt_frobenius_target(area).unwrap();

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
    fn test_twin_gkm_dagger_rejects_invalid_area() {
        assert!(twin_gkm_dagger_action_matrices(&[0, 2]).is_none());
        assert!(twin_gkm_dagger_character_values_by_degree(&[0, 2]).is_none());
        assert!(twin_gkm_dagger_frobenius(&[0, 2]).is_none());
    }

    #[test]
    fn test_twin_gkm_dagger_edgeless_s3_matches_llt() {
        assert_matches_llt_character_target(&[0, 0, 0]);
        assert_matches_llt_frobenius_target(&[0, 0, 0]);
    }

    #[test]
    fn test_twin_gkm_dagger_path_s3_matches_llt() {
        assert_matches_llt_character_target(&[0, 1, 1]);
        assert_matches_llt_frobenius_target(&[0, 1, 1]);
    }

    #[test]
    fn test_twin_gkm_dagger_path_s3_direct_traces_match_action_matrices() {
        let area = [0, 1, 1];
        let direct = twin_gkm_dagger_character_values_by_degree(&area).unwrap();
        let from_matrices =
            trace_values_from_action_matrices(twin_gkm_dagger_action_matrices(&area).unwrap());

        assert_eq!(direct, from_matrices);
    }

    #[test]
    fn test_twin_gkm_dagger_complete_graph_s3_matches_llt() {
        assert_matches_llt_character_target(&[0, 1, 2]);
        assert_matches_llt_frobenius_target(&[0, 1, 2]);
    }

    #[test]
    fn test_twin_gkm_dagger_complete_graph_s3_artin_matrices_match_llt() {
        let area = [0, 1, 2];
        let computed =
            graded_frobenius_from_trace_matrices(&twin_gkm_dagger_action_matrices(&area).unwrap());
        let target = unicellular_llt_frobenius_target(&area).unwrap();

        for degree in 0..=3 {
            let computed_schur = computed[&degree].to_schur_basis();
            let target_schur = target[&degree].to_schur_basis();
            for partition in Partition::all_of_size(area.len() as u32) {
                assert_eq!(
                    computed_schur.coefficient(&partition),
                    q(target_schur.coefficient(&partition)),
                    "mismatch in degree {degree}, partition {partition:?}"
                );
            }
        }
    }

    #[test]
    fn test_twin_gkm_dagger_complete_graph_s4_artin_presentation_matches_llt() {
        assert_matches_llt_character_target(&[0, 1, 2, 3]);
        assert_matches_llt_frobenius_target(&[0, 1, 2, 3]);
    }
}
