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
    matrix_trace, quotient_action_matrix, rref, zero_matrix, Matrix, QuotientSpace, Vector,
};
use sym_poly_core::packed_sparse_linear_algebra::{
    packed_sparse_kernel_basis_with_free_columns_from_rows, PackedSparseQuotientSpace,
    PackedSparseRow,
};
use sym_poly_core::sn_action::{assert_permutation, conjugacy_class_representatives, SnIndex};
use sym_poly_core::sparse_linear_algebra::{
    sparse_coefficient, sparse_kernel_basis_with_free_columns_from_rows, sparse_vector,
    SparseQuotientSpace, SparseVector,
};
use sym_poly_core::{chinese_remainder, symmetric_residue, Field, Partition, PrimeField, Ring};

use crate::frobenius::graded_frobenius_from_character_values;
use crate::SymmetricFunction;

type Q = Ratio<i64>;
type ResidueCharacterValues = BTreeMap<u32, BTreeMap<Partition, i128>>;

#[derive(Debug, Clone)]
struct HomogeneousGkmComponent {
    fixed_points: Vec<Vec<usize>>,
    sn_index: SnIndex,
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

#[derive(Debug, Clone)]
struct SparseHomogeneousGkmComponent<C: Field> {
    fixed_points: Vec<Vec<usize>>,
    sn_index: SnIndex,
    monomials: Vec<Vec<u32>>,
    monomial_index: BTreeMap<Vec<u32>, usize>,
    module_basis: Vec<SparseVector<C>>,
    module_coordinate_columns: Vec<usize>,
    module_coordinate_index_by_ambient_column: Vec<Option<usize>>,
    ordinary_quotient: SparseQuotientSpace<C>,
}

#[derive(Debug, Clone)]
struct SparseGkmModuleBasis<C: Field> {
    vectors: Vec<SparseVector<C>>,
    coordinate_columns: Vec<usize>,
}

#[derive(Debug, Clone)]
struct PackedHomogeneousGkmComponent {
    fixed_points: Vec<Vec<usize>>,
    sn_index: SnIndex,
    monomials: Vec<Vec<u32>>,
    monomial_index: BTreeMap<Vec<u32>, usize>,
    module_basis: Vec<PackedSparseRow>,
    module_coordinate_columns: Vec<usize>,
    module_coordinate_index_by_ambient_column: Vec<Option<usize>>,
    ordinary_quotient: PackedSparseQuotientSpace,
}

#[derive(Debug, Clone)]
struct PackedGkmModuleBasis {
    vectors: Vec<PackedSparseRow>,
    coordinate_columns: Vec<usize>,
}

#[derive(Debug, Clone)]
struct HessenbergGkmCombinatorics {
    fixed_points: Vec<Vec<usize>>,
    sn_index: SnIndex,
    edge_adjacencies: Vec<GkmEdgeAdjacency>,
    monomials_by_degree: Vec<Vec<Vec<u32>>>,
    monomial_index_by_degree: Vec<BTreeMap<Vec<u32>, usize>>,
    substitution_groups_by_degree: Vec<Vec<Vec<Vec<usize>>>>,
    multiplication_maps_by_degree: Vec<Vec<Vec<usize>>>,
    n: usize,
}

#[derive(Debug, Clone)]
struct GkmEdgeAdjacency {
    point_index: usize,
    adjacent_point_index: usize,
    label_a: usize,
    label_b: usize,
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
    let sn_index = SnIndex::new(n);

    let mut components = Vec::new();
    for degree in 0..=max_degree {
        let module_basis = homogeneous_gkm_module_basis(n, graph.edges(), degree);
        let mut component = HomogeneousGkmComponent::new(
            fixed_points.clone(),
            sn_index.clone(),
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

/// Compute graded character values of the Hessenberg GKM dot action.
///
/// This default high-level path runs sparse linear algebra over prime fields
/// and CRT-lifts the final integral character values.
pub fn hessenberg_gkm_dot_character_values_by_degree(
    area: &[u8],
) -> Option<BTreeMap<u32, BTreeMap<Partition, Ratio<i64>>>> {
    let values = hessenberg_gkm_dot_character_values_crt(area)?;
    Some(
        values
            .into_iter()
            .map(|(degree, degree_values)| {
                (
                    degree,
                    degree_values
                        .into_iter()
                        .map(|(cycle_type, value)| (cycle_type, Q::from_integer(value)))
                        .collect(),
                )
            })
            .collect(),
    )
}

/// Compute graded character values over a single prime field.
pub fn hessenberg_gkm_dot_character_values_mod_prime<const P: u64>(
    area: &[u8],
) -> Option<BTreeMap<u32, BTreeMap<Partition, PrimeField<P>>>> {
    if !is_area_sequence(area) {
        return None;
    }

    let graph = Graph::unit_interval(area);
    let n = graph.num_vertices();
    let components = sparse_gkm_components_from_edges::<PrimeField<P>>(
        n,
        graph.edges(),
        graph.num_edges() as u32,
    );
    let mut by_degree = BTreeMap::new();

    for (degree, component) in components.iter().enumerate() {
        if component.ordinary_quotient.dimension() == 0 {
            continue;
        }

        let mut values = BTreeMap::new();
        for (cycle_type, representative) in conjugacy_class_representatives(n) {
            let trace = component.ordinary_dot_trace(&representative);
            if !trace.is_zero() {
                values.insert(cycle_type, trace);
            }
        }
        by_degree.insert(degree as u32, values);
    }

    Some(by_degree)
}

/// Compute integer graded character values by CRT-lifting prime-field traces.
pub fn hessenberg_gkm_dot_character_values_crt(
    area: &[u8],
) -> Option<BTreeMap<u32, BTreeMap<Partition, i64>>> {
    if !is_area_sequence(area) {
        return None;
    }

    let graph = Graph::unit_interval(area);
    gkm_dot_character_values_crt_for_edges(
        graph.num_vertices(),
        graph.edges(),
        graph.num_edges() as u32,
    )
}

fn gkm_dot_character_values_crt_for_edges(
    n: usize,
    edges: &[(usize, usize)],
    max_degree: u32,
) -> Option<BTreeMap<u32, BTreeMap<Partition, i64>>> {
    let bound = character_trace_bound(n);
    let mut residues = Vec::new();
    let mut modulus = 1i128;

    push_prime_residues_for_edges::<1_000_000_007>(
        n,
        edges,
        max_degree,
        &mut residues,
        &mut modulus,
    )?;
    let required_modulus = bound.checked_mul(2).unwrap_or(i128::MAX);

    if modulus <= required_modulus {
        push_prime_residues_for_edges::<1_000_000_009>(
            n,
            edges,
            max_degree,
            &mut residues,
            &mut modulus,
        )?;
    }
    if modulus <= required_modulus {
        push_prime_residues_for_edges::<998_244_353>(
            n,
            edges,
            max_degree,
            &mut residues,
            &mut modulus,
        )?;
    }

    if modulus <= required_modulus {
        return None;
    }

    Some(lift_character_residues(n, &residues))
}

/// Compute graded character values over a byte-sized prime field.
///
/// This is the same GKM dot-action computation as
/// [`hessenberg_gkm_dot_character_values_mod_prime`], but all sparse rows store
/// coefficients as `u8`. It is an experimental backend for timing small-prime
/// linear algebra against the generic `PrimeField<P>` path.
pub fn hessenberg_gkm_dot_character_values_packed_mod_prime<const P: u8>(
    area: &[u8],
) -> Option<BTreeMap<u32, BTreeMap<Partition, u8>>> {
    if !is_area_sequence(area) {
        return None;
    }

    let n = area.len();
    let components = packed_hessenberg_gkm_components::<P>(area);
    let mut by_degree = BTreeMap::new();

    for (degree, component) in components.iter().enumerate() {
        if component.ordinary_quotient.dimension() == 0 {
            continue;
        }

        let mut values = BTreeMap::new();
        for (cycle_type, representative) in conjugacy_class_representatives(n) {
            let trace = component.ordinary_dot_trace_mod_prime::<P>(&representative);
            if trace != 0 {
                values.insert(cycle_type, trace);
            }
        }
        by_degree.insert(degree as u32, values);
    }

    Some(by_degree)
}

/// Compute integer graded character values by CRT-lifting packed small-prime traces.
pub fn hessenberg_gkm_dot_character_values_packed_crt(
    area: &[u8],
) -> Option<BTreeMap<u32, BTreeMap<Partition, i64>>> {
    if !is_area_sequence(area) {
        return None;
    }

    let bound = character_trace_bound(area.len());
    let mut residues = Vec::new();
    let mut modulus = 1i128;

    push_packed_prime_residues::<251>(area, &mut residues, &mut modulus)?;
    let required_modulus = bound.checked_mul(2).unwrap_or(i128::MAX);

    if modulus <= required_modulus {
        push_packed_prime_residues::<241>(area, &mut residues, &mut modulus)?;
    }
    if modulus <= required_modulus {
        push_packed_prime_residues::<239>(area, &mut residues, &mut modulus)?;
    }
    if modulus <= required_modulus {
        push_packed_prime_residues::<233>(area, &mut residues, &mut modulus)?;
    }
    if modulus <= required_modulus {
        push_packed_prime_residues::<229>(area, &mut residues, &mut modulus)?;
    }

    if modulus <= required_modulus {
        return None;
    }

    Some(lift_character_residues(area.len(), &residues))
}

/// Compute graded character values through the packed small-prime backend.
pub fn hessenberg_gkm_dot_character_values_by_degree_packed(
    area: &[u8],
) -> Option<BTreeMap<u32, BTreeMap<Partition, Ratio<i64>>>> {
    let values = hessenberg_gkm_dot_character_values_packed_crt(area)?;
    Some(
        values
            .into_iter()
            .map(|(degree, degree_values)| {
                (
                    degree,
                    degree_values
                        .into_iter()
                        .map(|(cycle_type, value)| (cycle_type, Q::from_integer(value)))
                        .collect(),
                )
            })
            .collect(),
    )
}

/// Rational reference implementation for small cases and regression checks.
pub fn hessenberg_gkm_dot_character_values_by_degree_rational(
    area: &[u8],
) -> Option<BTreeMap<u32, BTreeMap<Partition, Ratio<i64>>>> {
    let matrices = hessenberg_gkm_dot_action_matrices(area)?;
    Some(trace_values_from_action_matrices(matrices))
}

/// Compute the graded Frobenius characteristic of the GKM dot action.
pub fn hessenberg_gkm_dot_frobenius(
    area: &[u8],
) -> Option<BTreeMap<u32, SymmetricFunction<Ratio<i64>>>> {
    let character_values = hessenberg_gkm_dot_character_values_by_degree(area)?;
    Some(
        graded_frobenius_from_character_values(&character_values)
            .into_iter()
            .filter(|(_, f)| !f.is_zero())
            .collect(),
    )
}

/// Compute the graded Frobenius characteristic using packed small-prime traces.
pub fn hessenberg_gkm_dot_frobenius_packed(
    area: &[u8],
) -> Option<BTreeMap<u32, SymmetricFunction<Ratio<i64>>>> {
    let character_values = hessenberg_gkm_dot_character_values_by_degree_packed(area)?;
    Some(
        graded_frobenius_from_character_values(&character_values)
            .into_iter()
            .filter(|(_, f)| !f.is_zero())
            .collect(),
    )
}

/// Compute graded character values for the naive circular GKM dot-action model.
///
/// The circular area sequence supplies a directed circular unit arc digraph.
/// The naive GKM model forgets the edge orientation and imposes the ordinary
/// Hessenberg-style divisibility relation for each underlying transposition.
pub fn naive_circular_gkm_dot_character_values_by_degree(
    area: &[u8],
) -> Option<BTreeMap<u32, BTreeMap<Partition, Ratio<i64>>>> {
    let graph = Graph::circular_unit_interval(area)?;
    let values = gkm_dot_character_values_crt_for_edges(
        graph.num_vertices(),
        graph.edges(),
        graph.num_edges() as u32,
    )?;
    Some(
        values
            .into_iter()
            .map(|(degree, degree_values)| {
                (
                    degree,
                    degree_values
                        .into_iter()
                        .map(|(cycle_type, value)| (cycle_type, Q::from_integer(value)))
                        .collect(),
                )
            })
            .collect(),
    )
}

/// Compute the graded Frobenius characteristic of the naive circular GKM model.
pub fn naive_circular_gkm_dot_frobenius(
    area: &[u8],
) -> Option<BTreeMap<u32, SymmetricFunction<Ratio<i64>>>> {
    let character_values = naive_circular_gkm_dot_character_values_by_degree(area)?;
    Some(
        graded_frobenius_from_character_values(&character_values)
            .into_iter()
            .filter(|(_, f)| !f.is_zero())
            .collect(),
    )
}

impl HomogeneousGkmComponent {
    fn new(
        fixed_points: Vec<Vec<usize>>,
        sn_index: SnIndex,
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
            sn_index,
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
        let fixed_point_action =
            fixed_point_left_action_table(&self.sn_index, &self.fixed_points, permutation);
        let monomial_action =
            monomial_permutation_action_table(&self.monomials, &self.monomial_index, permutation);

        for (col, basis_vector) in self.module_basis.iter().enumerate() {
            let image = self.apply_dot_action_to_ambient_vector(
                &fixed_point_action,
                &monomial_action,
                basis_vector,
            );
            let coords = self.module_coordinates(&image);
            for row in 0..dim {
                matrix[row][col] = coords[row];
            }
        }

        matrix
    }

    fn apply_dot_action_to_ambient_vector(
        &self,
        fixed_point_action: &[usize],
        monomial_action: &[usize],
        vector: &[Q],
    ) -> Vector<Q> {
        let mut result = vec![Q::zero(); self.ambient_dimension()];
        let monomial_count = self.monomials.len();

        for (source_point_index, &target_point_index) in fixed_point_action.iter().enumerate() {
            for (source_monomial_index, &target_monomial_index) in
                monomial_action.iter().enumerate()
            {
                let source_col = source_point_index * monomial_count + source_monomial_index;
                let coeff = &vector[source_col];
                if coeff.is_zero() {
                    continue;
                }

                let target_row = target_point_index * monomial_count + target_monomial_index;
                result[target_row] = result[target_row].clone() + coeff.clone();
            }
        }

        result
    }
}

impl<C: Field> SparseHomogeneousGkmComponent<C> {
    fn new(
        fixed_points: Vec<Vec<usize>>,
        sn_index: SnIndex,
        monomials: Vec<Vec<u32>>,
        monomial_index: BTreeMap<Vec<u32>, usize>,
        module_basis: Vec<SparseVector<C>>,
        module_coordinate_columns: Vec<usize>,
        relations: Vec<SparseVector<C>>,
    ) -> Self {
        let ambient_dimension = fixed_points.len() * monomials.len();
        let module_coordinate_index_by_ambient_column =
            coordinate_index_by_ambient_column(ambient_dimension, &module_coordinate_columns);
        let ordinary_quotient = SparseQuotientSpace::from_relations(module_basis.len(), &relations);

        Self {
            fixed_points,
            sn_index,
            monomials,
            monomial_index,
            module_basis,
            module_coordinate_columns,
            module_coordinate_index_by_ambient_column,
            ordinary_quotient,
        }
    }

    fn ambient_dimension(&self) -> usize {
        self.fixed_points.len() * self.monomials.len()
    }

    fn module_coordinates_sparse(&self, vector: &SparseVector<C>) -> SparseVector<C> {
        let mut result = Vec::with_capacity(vector.len().min(self.module_coordinate_columns.len()));

        for &(ambient_col, ref coeff) in vector {
            if let Some(coordinate_index) =
                self.module_coordinate_index_by_ambient_column[ambient_col]
            {
                if !coeff.is_zero() {
                    result.push((coordinate_index, coeff.clone()));
                }
            }
        }

        result.sort_by_key(|(coordinate_index, _)| *coordinate_index);
        result
    }

    fn ordinary_dot_trace(&self, permutation: &[usize]) -> C {
        assert_permutation(permutation);
        let mut trace = C::zero();
        let fixed_point_action =
            fixed_point_left_action_table(&self.sn_index, &self.fixed_points, permutation);
        let monomial_action =
            monomial_permutation_action_table(&self.monomials, &self.monomial_index, permutation);

        for (quotient_col, &module_basis_index) in
            self.ordinary_quotient.free_columns.iter().enumerate()
        {
            let image = self.apply_dot_action_to_sparse_ambient_vector(
                &fixed_point_action,
                &monomial_action,
                &self.module_basis[module_basis_index],
            );
            let module_coords = self.module_coordinates_sparse(&image);
            let quotient_coords = self
                .ordinary_quotient
                .quotient_coordinates_sparse(&module_coords);
            trace = trace + sparse_coefficient(&quotient_coords, quotient_col);
        }

        trace
    }

    fn apply_dot_action_to_sparse_ambient_vector(
        &self,
        fixed_point_action: &[usize],
        monomial_action: &[usize],
        vector: &SparseVector<C>,
    ) -> SparseVector<C> {
        let monomial_count = self.monomials.len();
        let mut entries = Vec::with_capacity(vector.len());

        for &(source_col, ref coeff) in vector {
            if coeff.is_zero() {
                continue;
            }

            let source_point_index = source_col / monomial_count;
            let source_monomial_index = source_col % monomial_count;
            let target_point_index = fixed_point_action[source_point_index];
            let target_monomial_index = monomial_action[source_monomial_index];
            let target_col = target_point_index * monomial_count + target_monomial_index;
            entries.push((target_col, coeff.clone()));
        }

        sparse_vector(self.ambient_dimension(), entries)
    }
}

impl PackedHomogeneousGkmComponent {
    fn new<const P: u8>(
        fixed_points: Vec<Vec<usize>>,
        sn_index: SnIndex,
        monomials: Vec<Vec<u32>>,
        monomial_index: BTreeMap<Vec<u32>, usize>,
        module_basis: Vec<PackedSparseRow>,
        module_coordinate_columns: Vec<usize>,
        relations: Vec<PackedSparseRow>,
    ) -> Self {
        let ambient_dimension = fixed_points.len() * monomials.len();
        let module_coordinate_index_by_ambient_column =
            coordinate_index_by_ambient_column(ambient_dimension, &module_coordinate_columns);
        let ordinary_quotient =
            PackedSparseQuotientSpace::from_relations::<P>(module_basis.len(), &relations);

        Self {
            fixed_points,
            sn_index,
            monomials,
            monomial_index,
            module_basis,
            module_coordinate_columns,
            module_coordinate_index_by_ambient_column,
            ordinary_quotient,
        }
    }

    fn ambient_dimension(&self) -> usize {
        self.fixed_points.len() * self.monomials.len()
    }

    fn module_coordinates_sparse<const P: u8>(&self, vector: &PackedSparseRow) -> PackedSparseRow {
        let mut entries =
            Vec::with_capacity(vector.len().min(self.module_coordinate_columns.len()));

        for (ambient_col, coeff) in vector.cols.iter().copied().zip(vector.vals.iter().copied()) {
            let ambient_col = ambient_col as usize;
            if let Some(coordinate_index) =
                self.module_coordinate_index_by_ambient_column[ambient_col]
            {
                if coeff != 0 {
                    entries.push((coordinate_index, coeff));
                }
            }
        }

        PackedSparseRow::new::<P, _>(self.module_coordinate_columns.len(), entries)
    }

    fn ordinary_dot_trace_mod_prime<const P: u8>(&self, permutation: &[usize]) -> u8 {
        assert_permutation(permutation);
        let mut trace = 0u8;
        let fixed_point_action =
            fixed_point_left_action_table(&self.sn_index, &self.fixed_points, permutation);
        let monomial_action =
            monomial_permutation_action_table(&self.monomials, &self.monomial_index, permutation);

        for (quotient_col, &module_basis_index) in
            self.ordinary_quotient.free_columns.iter().enumerate()
        {
            let image = self.apply_dot_action_to_packed_ambient_vector::<P>(
                &fixed_point_action,
                &monomial_action,
                &self.module_basis[module_basis_index],
            );
            let module_coords = self.module_coordinates_sparse::<P>(&image);
            let quotient_coords = self
                .ordinary_quotient
                .quotient_coordinates_sparse::<P>(&module_coords);
            trace = add_mod_u8::<P>(trace, quotient_coords.coefficient(quotient_col));
        }

        trace
    }

    fn apply_dot_action_to_packed_ambient_vector<const P: u8>(
        &self,
        fixed_point_action: &[usize],
        monomial_action: &[usize],
        vector: &PackedSparseRow,
    ) -> PackedSparseRow {
        let monomial_count = self.monomials.len();
        let mut entries = Vec::with_capacity(vector.len());

        for (source_col, coeff) in vector.cols.iter().copied().zip(vector.vals.iter().copied()) {
            if coeff == 0 {
                continue;
            }

            let source_col = source_col as usize;
            let source_point_index = source_col / monomial_count;
            let source_monomial_index = source_col % monomial_count;
            let target_point_index = fixed_point_action[source_point_index];
            let target_monomial_index = monomial_action[source_monomial_index];
            let target_col = target_point_index * monomial_count + target_monomial_index;
            entries.push((target_col, coeff));
        }

        PackedSparseRow::new::<P, _>(self.ambient_dimension(), entries)
    }
}

impl HessenbergGkmCombinatorics {
    fn new(
        n: usize,
        hessenberg_edges: &[(usize, usize)],
        fixed_points: Vec<Vec<usize>>,
        max_degree: u32,
    ) -> Self {
        let sn_index = SnIndex::new(n);
        let edge_adjacencies =
            hessenberg_edge_adjacencies(n, hessenberg_edges, &fixed_points, &sn_index);
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
                (0..n * n)
                    .map(|label_key| {
                        let label_a = label_key / n;
                        let label_b = label_key % n;
                        substitution_groups_for_labels(monomials, label_a, label_b)
                    })
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
            sn_index,
            edge_adjacencies,
            monomials_by_degree,
            monomial_index_by_degree,
            substitution_groups_by_degree,
            multiplication_maps_by_degree,
            n,
        }
    }

    fn monomials(&self, degree: u32) -> &[Vec<u32>] {
        &self.monomials_by_degree[degree as usize]
    }

    fn monomial_index(&self, degree: u32) -> &BTreeMap<Vec<u32>, usize> {
        &self.monomial_index_by_degree[degree as usize]
    }

    fn substitution_groups(&self, degree: u32, label_a: usize, label_b: usize) -> &[Vec<usize>] {
        &self.substitution_groups_by_degree[degree as usize][label_a * self.n + label_b]
    }

    fn multiplication_maps(&self, degree: u32) -> &[Vec<usize>] {
        &self.multiplication_maps_by_degree[degree as usize]
    }
}

fn sparse_gkm_components_from_edges<C: Field>(
    n: usize,
    edges: &[(usize, usize)],
    max_degree: u32,
) -> Vec<SparseHomogeneousGkmComponent<C>> {
    let fixed_points = sym_poly_core::symmetric_group_permutation_basis(n);
    let combinatorics = HessenbergGkmCombinatorics::new(n, edges, fixed_points, max_degree);

    let mut components = Vec::new();
    for degree in 0..=max_degree {
        let module_basis = sparse_homogeneous_gkm_module_basis::<C>(&combinatorics, degree);
        let mut component = SparseHomogeneousGkmComponent::new(
            combinatorics.fixed_points.clone(),
            combinatorics.sn_index.clone(),
            combinatorics.monomials(degree).to_vec(),
            combinatorics.monomial_index(degree).clone(),
            module_basis.vectors,
            module_basis.coordinate_columns,
            Vec::new(),
        );
        let relations = if degree == 0 {
            Vec::new()
        } else {
            sparse_variable_multiple_relations(
                &components[(degree - 1) as usize],
                &component,
                combinatorics.multiplication_maps(degree),
            )
        };
        component.ordinary_quotient =
            SparseQuotientSpace::from_relations(component.module_basis.len(), &relations);
        components.push(component);
    }

    components
}

fn packed_hessenberg_gkm_components<const P: u8>(
    area: &[u8],
) -> Vec<PackedHomogeneousGkmComponent> {
    let graph = Graph::unit_interval(area);
    let n = area.len();
    let max_degree = graph.num_edges() as u32;
    let fixed_points = sym_poly_core::symmetric_group_permutation_basis(n);
    let combinatorics =
        HessenbergGkmCombinatorics::new(n, graph.edges(), fixed_points.clone(), max_degree);

    let mut components = Vec::new();
    for degree in 0..=max_degree {
        let module_basis = packed_homogeneous_gkm_module_basis::<P>(&combinatorics, degree);
        let mut component = PackedHomogeneousGkmComponent::new::<P>(
            combinatorics.fixed_points.clone(),
            combinatorics.sn_index.clone(),
            combinatorics.monomials(degree).to_vec(),
            combinatorics.monomial_index(degree).clone(),
            module_basis.vectors,
            module_basis.coordinate_columns,
            Vec::new(),
        );
        let relations = if degree == 0 {
            Vec::new()
        } else {
            packed_variable_multiple_relations::<P>(
                &components[(degree - 1) as usize],
                &component,
                combinatorics.multiplication_maps(degree),
            )
        };
        component.ordinary_quotient = PackedSparseQuotientSpace::from_relations::<P>(
            component.module_basis.len(),
            &relations,
        );
        components.push(component);
    }

    components
}

fn homogeneous_gkm_module_basis(
    n: usize,
    hessenberg_edges: &[(usize, usize)],
    degree: u32,
) -> GkmModuleBasis {
    let fixed_points = sym_poly_core::symmetric_group_permutation_basis(n);
    let sn_index = SnIndex::new(n);
    let monomials = homogeneous_monomials(n, degree);
    let ambient_dimension = fixed_points.len() * monomials.len();
    let mut constraints = Vec::new();

    for (point_index, point) in fixed_points.iter().enumerate() {
        for &(i, j) in hessenberg_edges {
            let adjacent_point_index =
                right_transposition_fixed_point_index(&sn_index, point, i, j);
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

fn sparse_homogeneous_gkm_module_basis<C: Field>(
    combinatorics: &HessenbergGkmCombinatorics,
    degree: u32,
) -> SparseGkmModuleBasis<C> {
    let monomial_count = combinatorics.monomials(degree).len();
    let ambient_dimension = combinatorics.fixed_points.len() * monomial_count;

    if combinatorics.edge_adjacencies.is_empty() {
        SparseGkmModuleBasis {
            vectors: sparse_standard_basis(ambient_dimension),
            coordinate_columns: (0..ambient_dimension).collect(),
        }
    } else {
        let rows = combinatorics.edge_adjacencies.iter().flat_map(|adjacency| {
            combinatorics
                .substitution_groups(degree, adjacency.label_a, adjacency.label_b)
                .iter()
                .map(|group| {
                    sparse_divisibility_constraint_for_group(
                        adjacency.point_index,
                        adjacency.adjacent_point_index,
                        group,
                        monomial_count,
                        ambient_dimension,
                    )
                })
        });
        let (vectors, coordinate_columns) =
            sparse_kernel_basis_with_free_columns_from_rows(ambient_dimension, rows);
        SparseGkmModuleBasis {
            vectors,
            coordinate_columns,
        }
    }
}

fn sparse_divisibility_constraint_for_group<C: Ring>(
    point_index: usize,
    adjacent_point_index: usize,
    group: &[usize],
    monomial_count: usize,
    ambient_dimension: usize,
) -> SparseVector<C> {
    let mut entries = Vec::with_capacity(2 * group.len());
    for &monomial_idx in group {
        entries.push((point_index * monomial_count + monomial_idx, C::one()));
        entries.push((
            adjacent_point_index * monomial_count + monomial_idx,
            -C::one(),
        ));
    }
    sparse_vector(ambient_dimension, entries)
}

fn sparse_variable_multiple_relations<C: Field>(
    previous: &SparseHomogeneousGkmComponent<C>,
    current: &SparseHomogeneousGkmComponent<C>,
    multiplication_maps: &[Vec<usize>],
) -> Vec<SparseVector<C>> {
    let mut relations = Vec::new();
    let previous_monomial_count = previous.monomials.len();
    let current_monomial_count = current.monomials.len();

    for basis_vector in &previous.module_basis {
        for variable_map in multiplication_maps {
            let mut entries = Vec::new();
            for &(source, ref coeff) in basis_vector {
                if coeff.is_zero() {
                    continue;
                }

                let point_index = source / previous_monomial_count;
                let monomial_index = source % previous_monomial_count;
                let target_monomial_index = variable_map[monomial_index];
                let target = point_index * current_monomial_count + target_monomial_index;
                entries.push((target, coeff.clone()));
            }
            let ambient = sparse_vector(current.ambient_dimension(), entries);
            relations.push(current.module_coordinates_sparse(&ambient));
        }
    }

    relations
}

fn packed_homogeneous_gkm_module_basis<const P: u8>(
    combinatorics: &HessenbergGkmCombinatorics,
    degree: u32,
) -> PackedGkmModuleBasis {
    let monomial_count = combinatorics.monomials(degree).len();
    let ambient_dimension = combinatorics.fixed_points.len() * monomial_count;

    if combinatorics.edge_adjacencies.is_empty() {
        PackedGkmModuleBasis {
            vectors: packed_sparse_standard_basis::<P>(ambient_dimension),
            coordinate_columns: (0..ambient_dimension).collect(),
        }
    } else {
        let rows = combinatorics.edge_adjacencies.iter().flat_map(|adjacency| {
            combinatorics
                .substitution_groups(degree, adjacency.label_a, adjacency.label_b)
                .iter()
                .map(|group| {
                    packed_divisibility_constraint_for_group::<P>(
                        adjacency.point_index,
                        adjacency.adjacent_point_index,
                        group,
                        monomial_count,
                        ambient_dimension,
                    )
                })
        });
        let (vectors, coordinate_columns) =
            packed_sparse_kernel_basis_with_free_columns_from_rows::<P, _>(ambient_dimension, rows);
        PackedGkmModuleBasis {
            vectors,
            coordinate_columns,
        }
    }
}

fn packed_divisibility_constraint_for_group<const P: u8>(
    point_index: usize,
    adjacent_point_index: usize,
    group: &[usize],
    monomial_count: usize,
    ambient_dimension: usize,
) -> PackedSparseRow {
    let mut entries = Vec::with_capacity(2 * group.len());
    for &monomial_idx in group {
        entries.push((point_index * monomial_count + monomial_idx, 1));
        entries.push((adjacent_point_index * monomial_count + monomial_idx, P - 1));
    }
    PackedSparseRow::new::<P, _>(ambient_dimension, entries)
}

fn packed_variable_multiple_relations<const P: u8>(
    previous: &PackedHomogeneousGkmComponent,
    current: &PackedHomogeneousGkmComponent,
    multiplication_maps: &[Vec<usize>],
) -> Vec<PackedSparseRow> {
    let mut relations = Vec::new();
    let previous_monomial_count = previous.monomials.len();
    let current_monomial_count = current.monomials.len();

    for basis_vector in &previous.module_basis {
        for variable_map in multiplication_maps {
            let mut entries = Vec::new();
            for (source, coeff) in basis_vector
                .cols
                .iter()
                .copied()
                .zip(basis_vector.vals.iter().copied())
            {
                if coeff == 0 {
                    continue;
                }

                let source = source as usize;
                let point_index = source / previous_monomial_count;
                let monomial_index = source % previous_monomial_count;
                let target_monomial_index = variable_map[monomial_index];
                let target = point_index * current_monomial_count + target_monomial_index;
                entries.push((target, coeff));
            }
            let ambient = PackedSparseRow::new::<P, _>(current.ambient_dimension(), entries);
            relations.push(current.module_coordinates_sparse::<P>(&ambient));
        }
    }

    relations
}

fn hessenberg_edge_adjacencies(
    _n: usize,
    hessenberg_edges: &[(usize, usize)],
    fixed_points: &[Vec<usize>],
    sn_index: &SnIndex,
) -> Vec<GkmEdgeAdjacency> {
    let mut result = Vec::new();
    for (point_index, point) in fixed_points.iter().enumerate() {
        for &(i, j) in hessenberg_edges {
            let adjacent_point_index = right_transposition_fixed_point_index(sn_index, point, i, j);
            if point_index > adjacent_point_index {
                continue;
            }
            result.push(GkmEdgeAdjacency {
                point_index,
                adjacent_point_index,
                label_a: point[i],
                label_b: point[j],
            });
        }
    }
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

fn sparse_standard_basis<C: Ring>(dimension: usize) -> Vec<SparseVector<C>> {
    (0..dimension)
        .map(|index| vec![(index, C::one())])
        .collect()
}

fn packed_sparse_standard_basis<const P: u8>(dimension: usize) -> Vec<PackedSparseRow> {
    (0..dimension)
        .map(|index| PackedSparseRow::new::<P, _>(dimension, vec![(index, 1)]))
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

fn coordinate_index_by_ambient_column(
    ambient_dimension: usize,
    coordinate_columns: &[usize],
) -> Vec<Option<usize>> {
    let mut lookup = vec![None; ambient_dimension];
    for (coordinate_index, &ambient_col) in coordinate_columns.iter().enumerate() {
        lookup[ambient_col] = Some(coordinate_index);
    }
    lookup
}

fn fixed_point_left_action_table(
    sn_index: &SnIndex,
    fixed_points: &[Vec<usize>],
    permutation: &[usize],
) -> Vec<usize> {
    fixed_points
        .iter()
        .map(|point| sn_index.rank_left_composition(permutation, point))
        .collect()
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

fn monomial_permutation_action_table(
    monomials: &[Vec<u32>],
    monomial_index: &BTreeMap<Vec<u32>, usize>,
    permutation: &[usize],
) -> Vec<usize> {
    monomials
        .iter()
        .map(|monomial| {
            let image = permute_monomial_variables(monomial, permutation);
            monomial_index[&image]
        })
        .collect()
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

fn trace_values_from_action_matrices<C: Ring>(
    matrices: BTreeMap<u32, BTreeMap<Partition, Matrix<C>>>,
) -> BTreeMap<u32, BTreeMap<Partition, C>> {
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

fn push_prime_residues_for_edges<const P: u64>(
    n: usize,
    edges: &[(usize, usize)],
    max_degree: u32,
    residues: &mut Vec<(i128, ResidueCharacterValues)>,
    modulus_product: &mut i128,
) -> Option<()> {
    let components = sparse_gkm_components_from_edges::<PrimeField<P>>(n, edges, max_degree);
    let mut values = BTreeMap::new();

    for (degree, component) in components.iter().enumerate() {
        if component.ordinary_quotient.dimension() == 0 {
            continue;
        }

        let mut degree_values = BTreeMap::new();
        for (cycle_type, representative) in conjugacy_class_representatives(n) {
            let trace = component.ordinary_dot_trace(&representative);
            if !trace.is_zero() {
                degree_values.insert(cycle_type, trace);
            }
        }
        values.insert(degree as u32, degree_values);
    }

    residues.push((P as i128, prime_values_to_residues(values)));
    *modulus_product = modulus_product.checked_mul(P as i128)?;
    Some(())
}

fn prime_values_to_residues<const P: u64>(
    values: BTreeMap<u32, BTreeMap<Partition, PrimeField<P>>>,
) -> ResidueCharacterValues {
    values
        .into_iter()
        .map(|(degree, degree_values)| {
            (
                degree,
                degree_values
                    .into_iter()
                    .map(|(cycle_type, value)| (cycle_type, value.value() as i128))
                    .collect(),
            )
        })
        .collect()
}

fn push_packed_prime_residues<const P: u8>(
    area: &[u8],
    residues: &mut Vec<(i128, ResidueCharacterValues)>,
    modulus_product: &mut i128,
) -> Option<()> {
    let values = hessenberg_gkm_dot_character_values_packed_mod_prime::<P>(area)?;
    residues.push((P as i128, packed_prime_values_to_residues(values)));
    *modulus_product = modulus_product.checked_mul(P as i128)?;
    Some(())
}

fn packed_prime_values_to_residues(
    values: BTreeMap<u32, BTreeMap<Partition, u8>>,
) -> ResidueCharacterValues {
    values
        .into_iter()
        .map(|(degree, degree_values)| {
            (
                degree,
                degree_values
                    .into_iter()
                    .map(|(cycle_type, value)| (cycle_type, value as i128))
                    .collect(),
            )
        })
        .collect()
}

fn lift_character_residues(
    n: usize,
    residues: &[(i128, ResidueCharacterValues)],
) -> BTreeMap<u32, BTreeMap<Partition, i64>> {
    let mut degrees = std::collections::BTreeSet::new();
    for (_, values) in residues {
        degrees.extend(values.keys().copied());
    }

    let cycle_types = Partition::all_of_size(n as u32);
    degrees
        .into_iter()
        .map(|degree| {
            let mut degree_values = BTreeMap::new();
            for cycle_type in &cycle_types {
                let congruences = residues
                    .iter()
                    .map(|(prime, values)| {
                        let residue = values
                            .get(&degree)
                            .and_then(|degree_values| degree_values.get(cycle_type))
                            .copied()
                            .unwrap_or(0);
                        (residue, *prime)
                    })
                    .collect::<Vec<_>>();
                let (residue, modulus) =
                    chinese_remainder(&congruences).expect("prime moduli are coprime");
                let value = symmetric_residue(residue, modulus);
                if value != 0 {
                    degree_values.insert(
                        cycle_type.clone(),
                        i64::try_from(value).expect("lifted character value fits in i64"),
                    );
                }
            }
            (degree, degree_values)
        })
        .collect()
}

fn character_trace_bound(n: usize) -> i128 {
    (1..=n).fold(1i128, |acc, value| acc.saturating_mul(value as i128))
}

fn add_mod_u8<const P: u8>(a: u8, b: u8) -> u8 {
    let sum = a as u16 + b as u16;
    if sum >= P as u16 {
        (sum - P as u16) as u8
    } else {
        sum as u8
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chromatic::{
        circular_area_dot_frobenius_target, hessenberg_area_dot_frobenius_target,
    };
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

    fn assert_crt_matches_rational_reference(area: &[u8]) {
        let crt = hessenberg_gkm_dot_character_values_by_degree(area).unwrap();
        let rational = hessenberg_gkm_dot_character_values_by_degree_rational(area).unwrap();

        assert_eq!(crt, rational);
    }

    fn assert_packed_matches_generic(area: &[u8]) {
        let packed = hessenberg_gkm_dot_character_values_by_degree_packed(area).unwrap();
        let generic = hessenberg_gkm_dot_character_values_by_degree(area).unwrap();

        assert_eq!(packed, generic);
    }

    fn assert_packed_mod_prime_matches_generic_mod_prime(area: &[u8]) {
        let packed = hessenberg_gkm_dot_character_values_packed_mod_prime::<251>(area).unwrap();
        let generic = hessenberg_gkm_dot_character_values_mod_prime::<251>(area).unwrap();
        let generic_as_u8 = generic
            .into_iter()
            .map(|(degree, degree_values)| {
                (
                    degree,
                    degree_values
                        .into_iter()
                        .map(|(cycle_type, value)| (cycle_type, value.value() as u8))
                        .collect::<BTreeMap<_, _>>(),
                )
            })
            .collect::<BTreeMap<_, _>>();

        assert_eq!(packed, generic_as_u8);
    }

    fn assert_naive_circular_matches_target(area: &[u8]) {
        let computed = naive_circular_gkm_dot_frobenius(area).unwrap();
        let target = circular_area_dot_frobenius_target(area).unwrap();

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
        assert!(hessenberg_gkm_dot_character_values_by_degree(&[0, 2]).is_none());
        assert!(hessenberg_gkm_dot_frobenius_packed(&[0, 2]).is_none());
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

    #[test]
    fn test_hessenberg_gkm_dot_path_s3_crt_matches_rational_reference() {
        assert_crt_matches_rational_reference(&[0, 1, 1]);
    }

    #[test]
    fn test_hessenberg_gkm_dot_complete_graph_s3_crt_matches_rational_reference() {
        assert_crt_matches_rational_reference(&[0, 1, 2]);
    }

    #[test]
    fn test_hessenberg_gkm_dot_path_s3_packed_matches_generic_mod_prime() {
        assert_packed_mod_prime_matches_generic_mod_prime(&[0, 1, 1]);
    }

    #[test]
    fn test_hessenberg_gkm_dot_path_s3_packed_matches_generic_crt() {
        assert_packed_matches_generic(&[0, 1, 1]);
    }

    #[test]
    fn test_hessenberg_gkm_dot_complete_graph_s3_packed_matches_generic_crt() {
        assert_packed_matches_generic(&[0, 1, 2]);
    }

    #[test]
    fn test_naive_circular_gkm_extends_unit_interval_case() {
        assert_naive_circular_matches_target(&[0, 1, 1]);
    }

    #[test]
    fn test_naive_circular_gkm_detects_directed_cycle_mismatch() {
        let computed = naive_circular_gkm_dot_frobenius(&[1, 1, 1]).unwrap();
        let target = circular_area_dot_frobenius_target(&[1, 1, 1]).unwrap();

        assert_ne!(computed.keys().copied().collect::<Vec<_>>(), vec![1, 2]);
        assert_eq!(target.keys().copied().collect::<Vec<_>>(), vec![1, 2]);

        let computed_degree_zero = computed[&0].to_schur_basis();
        assert_eq!(
            computed_degree_zero.coefficient(&Partition::new(vec![3])),
            q(1)
        );

        let target_degree_one = target[&1].to_schur_basis();
        assert_eq!(target_degree_one.coefficient(&Partition::new(vec![3])), 3);
    }
}
