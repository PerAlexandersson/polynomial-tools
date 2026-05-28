use std::collections::BTreeMap;
use std::env;

use combinatoric_core::Graph;
use num_traits::ToPrimitive;
use sym_poly_core::sparse_linear_algebra::{
    sparse_kernel_basis_with_free_columns_from_rows, sparse_vector, SparseQuotientSpace,
    SparseVector,
};
use sym_poly_core::{PrimeField, Ring};
use sym_poly_sym::circular_area_dot_frobenius_target;

type F = PrimeField<1_000_000_007>;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct AffinePoint {
    residues: Vec<usize>,
    shifts: Vec<i32>,
}

#[derive(Debug, Clone, Copy)]
struct CircularEdge {
    source: usize,
    target: usize,
    wrap: i32,
}

#[derive(Debug, Clone)]
struct AffineAdjacency {
    point_index: usize,
    adjacent_point_index: usize,
    label_a: usize,
    label_b: usize,
    h_coeff: i32,
    z_shift: i32,
}

#[derive(Debug, Clone)]
struct AffineWindowCombinatorics {
    points: Vec<AffinePoint>,
    adjacencies: Vec<AffineAdjacency>,
    omitted_adjacencies: usize,
    variable_count: usize,
}

#[derive(Debug, Clone)]
struct Component {
    monomials: Vec<Vec<u32>>,
    monomial_index: BTreeMap<Vec<u32>, usize>,
    module_basis: Vec<SparseVector<F>>,
    coordinate_index_by_ambient_column: Vec<Option<usize>>,
    quotient: SparseQuotientSpace<F>,
}

#[derive(Debug, Clone)]
struct WindowResult {
    fixed_points: usize,
    omitted_adjacencies: usize,
    hilbert: BTreeMap<u32, usize>,
}

#[derive(Debug, Clone, Copy)]
enum AffineModel {
    Window { radius: i32 },
    Cyclic { period: i32 },
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let (model, monodromy, input_args) = parse_options(&args).unwrap_or_else(|message| {
        eprintln!("{message}");
        std::process::exit(2);
    });
    let areas = parse_input(input_args).unwrap_or_else(|message| {
        eprintln!("{message}");
        eprintln!("usage:");
        eprintln!("  cargo run -p sym-poly-sym --example circular_dot_affine_window -- 1,1,1");
        eprintln!(
            "  cargo run -p sym-poly-sym --example circular_dot_affine_window -- --radius 2 1,1,1"
        );
        eprintln!(
            "  cargo run -p sym-poly-sym --example circular_dot_affine_window -- --period 2 1,1,1"
        );
        eprintln!(
            "  cargo run -p sym-poly-sym --example circular_dot_affine_window -- --monodromy -1 1,1,1"
        );
        eprintln!("  cargo run -p sym-poly-sym --example circular_dot_affine_window -- 3");
        std::process::exit(2);
    });

    let mut matches = 0usize;
    for area in &areas {
        if print_area_sequence(area, model, monodromy) {
            matches += 1;
        }
    }
    if areas.len() > 1 {
        println!("summary: {matches}/{} Hilbert series matched", areas.len());
    }
}

fn parse_options(args: &[String]) -> Result<(AffineModel, i64, &[String]), String> {
    let mut model = None;
    let mut monodromy = 2;
    let mut index = 0;

    while index < args.len() {
        match args[index].as_str() {
            "--radius" => {
                if index + 1 >= args.len() {
                    return Err("--radius needs a nonnegative integer".to_string());
                }
                let radius = args[index + 1].parse::<i32>().map_err(|_| {
                    format!("expected nonnegative radius, got {:?}", args[index + 1])
                })?;
                if radius < 0 {
                    return Err("--radius must be nonnegative".to_string());
                }
                if model.is_some() {
                    return Err("use at most one of --radius and --period".to_string());
                }
                model = Some(AffineModel::Window { radius });
                index += 2;
            }
            "--period" => {
                if index + 1 >= args.len() {
                    return Err("--period needs a positive integer".to_string());
                }
                let period = args[index + 1]
                    .parse::<i32>()
                    .map_err(|_| format!("expected positive period, got {:?}", args[index + 1]))?;
                if period <= 0 {
                    return Err("--period must be positive".to_string());
                }
                if model.is_some() {
                    return Err("use at most one of --radius and --period".to_string());
                }
                model = Some(AffineModel::Cyclic { period });
                index += 2;
            }
            "--monodromy" => {
                if index + 1 >= args.len() {
                    return Err("--monodromy needs a nonzero integer".to_string());
                }
                monodromy = args[index + 1].parse::<i64>().map_err(|_| {
                    format!("expected nonzero monodromy, got {:?}", args[index + 1])
                })?;
                if monodromy == 0 {
                    return Err("--monodromy must be nonzero".to_string());
                }
                index += 2;
            }
            _ => break,
        }
    }

    Ok((
        model.unwrap_or(AffineModel::Window { radius: 1 }),
        monodromy,
        &args[index..],
    ))
}

fn parse_input(args: &[String]) -> Result<Vec<Vec<u8>>, String> {
    if args.is_empty() {
        return Ok(vec![vec![1, 1, 1]]);
    }

    if args.len() == 1 && !args[0].contains(',') {
        let n = args[0]
            .parse::<usize>()
            .map_err(|_| format!("expected a rank or an area sequence, got {:?}", args[0]))?;
        if n > u8::MAX as usize {
            return Err("rank must be at most 255 for u8 area sequences".to_string());
        }
        return Ok(all_circular_area_sequences(n));
    }

    let area = if args.len() == 1 {
        args[0]
            .split(',')
            .map(parse_area_entry)
            .collect::<Result<Vec<_>, _>>()?
    } else {
        args.iter()
            .map(|arg| parse_area_entry(arg))
            .collect::<Result<Vec<_>, _>>()?
    };

    Ok(vec![area])
}

fn parse_area_entry(text: &str) -> Result<u8, String> {
    text.parse::<u8>()
        .map_err(|_| format!("expected a nonnegative area entry, got {text:?}"))
}

fn print_area_sequence(area: &[u8], model: AffineModel, monodromy: i64) -> bool {
    let Some(result) = affine_window_hilbert(area, model, monodromy) else {
        println!("area {:?}: invalid circular area sequence", area);
        return false;
    };
    let Some(target) = circular_area_dot_frobenius_target(area) else {
        println!("area {:?}: could not compute target", area);
        return false;
    };
    let target_hilbert = frobenius_hilbert(&target);
    let matches_target = result.hilbert == target_hilbert;

    println!(
        "area {:?}, {}, monodromy {monodromy}",
        area,
        format_model(model)
    );
    println!(
        "fixed points: {}, omitted boundary adjacencies: {}",
        result.fixed_points, result.omitted_adjacencies
    );
    println!("affine-window Hilbert matches target dimensions: {matches_target}");
    println!("affine window: {}", format_hilbert(&result.hilbert));
    if !matches_target {
        println!("target:        {}", format_hilbert(&target_hilbert));
    }
    println!();

    matches_target
}

fn affine_window_hilbert(area: &[u8], model: AffineModel, monodromy: i64) -> Option<WindowResult> {
    if monodromy == 0 {
        return None;
    }

    let edges = circular_edges(area)?;
    let combinatorics = match model {
        AffineModel::Window { radius } => affine_window_combinatorics(area.len(), &edges, radius),
        AffineModel::Cyclic { period } => affine_cyclic_combinatorics(area.len(), &edges, period),
    };
    let max_degree = edges.len() as u32;
    let mut components = Vec::new();
    let mut hilbert = BTreeMap::new();

    for degree in 0..=max_degree {
        let mut component = homogeneous_component(&combinatorics, degree, monodromy);
        let relations = if degree == 0 {
            Vec::new()
        } else {
            variable_multiple_relations(
                &components[(degree - 1) as usize],
                &component,
                combinatorics.points.len(),
                combinatorics.variable_count,
            )
        };
        component.quotient =
            SparseQuotientSpace::from_relations(component.module_basis.len(), &relations);
        let dimension = component.quotient.dimension();
        if dimension != 0 {
            hilbert.insert(degree, dimension);
        }
        components.push(component);
    }

    Some(WindowResult {
        fixed_points: combinatorics.points.len(),
        omitted_adjacencies: combinatorics.omitted_adjacencies,
        hilbert,
    })
}

fn affine_window_combinatorics(
    n: usize,
    edges: &[CircularEdge],
    radius: i32,
) -> AffineWindowCombinatorics {
    let points = affine_points(n, radius);
    let point_index = points
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, point)| (point, index))
        .collect::<BTreeMap<_, _>>();
    let mut adjacencies = Vec::new();
    let mut omitted_adjacencies = 0usize;

    for (point_index_value, point) in points.iter().enumerate() {
        for edge in edges {
            let adjacent = reflect_affine_point(point, *edge);
            let Some(&adjacent_point_index) = point_index.get(&adjacent) else {
                omitted_adjacencies += 1;
                continue;
            };
            if point_index_value > adjacent_point_index {
                continue;
            }

            adjacencies.push(AffineAdjacency {
                point_index: point_index_value,
                adjacent_point_index,
                label_a: point.residues[edge.source],
                label_b: point.residues[edge.target],
                h_coeff: point.shifts[edge.source] - point.shifts[edge.target] - edge.wrap,
                z_shift: edge.wrap,
            });
        }
    }

    AffineWindowCombinatorics {
        points,
        adjacencies,
        omitted_adjacencies,
        variable_count: n + 1,
    }
}

fn affine_cyclic_combinatorics(
    n: usize,
    edges: &[CircularEdge],
    period: i32,
) -> AffineWindowCombinatorics {
    debug_assert!(period > 0);
    let points = affine_cyclic_points(n, period);
    let point_index = points
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, point)| (point, index))
        .collect::<BTreeMap<_, _>>();
    let mut adjacencies = Vec::new();

    for (point_index_value, point) in points.iter().enumerate() {
        for edge in edges {
            let adjacent = reflect_affine_point_mod_period(point, *edge, period);
            let adjacent_point_index = point_index[&adjacent];
            if point_index_value > adjacent_point_index {
                continue;
            }

            adjacencies.push(AffineAdjacency {
                point_index: point_index_value,
                adjacent_point_index,
                label_a: point.residues[edge.source],
                label_b: point.residues[edge.target],
                h_coeff: point.shifts[edge.source] - point.shifts[edge.target] - edge.wrap,
                z_shift: edge.wrap,
            });
        }
    }

    AffineWindowCombinatorics {
        points,
        adjacencies,
        omitted_adjacencies: 0,
        variable_count: n + 1,
    }
}

fn homogeneous_component(
    combinatorics: &AffineWindowCombinatorics,
    degree: u32,
    monodromy: i64,
) -> Component {
    let monomials = homogeneous_monomials(combinatorics.variable_count, degree);
    let monomial_index = monomials
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, monomial)| (monomial, index))
        .collect::<BTreeMap<_, _>>();
    let monomial_count = monomials.len();
    let ambient_dimension = combinatorics.points.len() * monomial_count;

    let rows = combinatorics.adjacencies.iter().flat_map(|adjacency| {
        affine_window_constraints_for_edge(
            adjacency,
            combinatorics.variable_count - 1,
            &monomials,
            ambient_dimension,
            monodromy,
        )
    });
    let (module_basis, coordinate_columns) =
        sparse_kernel_basis_with_free_columns_from_rows(ambient_dimension, rows);
    let coordinate_index_by_ambient_column =
        coordinate_index_by_ambient_column(ambient_dimension, &coordinate_columns);
    let quotient = SparseQuotientSpace::from_relations(module_basis.len(), &[]);

    Component {
        monomials,
        monomial_index,
        module_basis,
        coordinate_index_by_ambient_column,
        quotient,
    }
}

fn affine_window_constraints_for_edge(
    adjacency: &AffineAdjacency,
    hbar_index: usize,
    monomials: &[Vec<u32>],
    ambient_dimension: usize,
    monodromy: i64,
) -> Vec<SparseVector<F>> {
    let monomial_count = monomials.len();
    let edge_monodromy = field_power_i64(monodromy, adjacency.z_shift);
    let mut rows_by_key: BTreeMap<Vec<u32>, Vec<(usize, F)>> = BTreeMap::new();

    for (monomial_index, monomial) in monomials.iter().enumerate() {
        for (key, coeff) in affine_substitute_label(
            monomial,
            adjacency.label_a,
            adjacency.label_b,
            adjacency.h_coeff,
            hbar_index,
        ) {
            let row = rows_by_key.entry(key).or_default();
            row.push((
                adjacency.point_index * monomial_count + monomial_index,
                F::from_i64(coeff),
            ));
            row.push((
                adjacency.adjacent_point_index * monomial_count + monomial_index,
                -(F::from_i64(coeff) * edge_monodromy),
            ));
        }
    }

    rows_by_key
        .into_values()
        .map(|entries| sparse_vector(ambient_dimension, entries))
        .filter(|row| !row.is_empty())
        .collect()
}

fn field_power_i64(base: i64, exponent: i32) -> F {
    assert!(exponent >= 0, "negative z-shifts are not supported");
    let mut result = F::from_i64(1);
    let factor = F::from_i64(base);
    for _ in 0..exponent {
        result = result * factor;
    }
    result
}

fn affine_substitute_label(
    monomial: &[u32],
    label_a: usize,
    label_b: usize,
    h_coeff: i32,
    hbar_index: usize,
) -> BTreeMap<Vec<u32>, i64> {
    let mut base = monomial.to_vec();
    let exponent = base[label_a];
    base[label_a] = 0;

    let mut result = BTreeMap::new();
    for hbar_power in 0..=exponent {
        let mut expanded = base.clone();
        expanded[label_b] += exponent - hbar_power;
        expanded[hbar_index] += hbar_power;
        let coeff = binomial_u32(exponent, hbar_power) as i64 * i64::from(-h_coeff).pow(hbar_power);
        *result.entry(expanded).or_insert(0) += coeff;
    }
    result.retain(|_, coeff| *coeff != 0);
    result
}

fn variable_multiple_relations(
    previous: &Component,
    current: &Component,
    point_count: usize,
    variable_count: usize,
) -> Vec<SparseVector<F>> {
    let mut relations = Vec::new();
    let previous_monomial_count = previous.monomials.len();
    let current_monomial_count = current.monomials.len();

    for basis_vector in &previous.module_basis {
        for variable in 0..variable_count {
            let mut entries = Vec::new();
            for &(source, coeff) in basis_vector {
                let point_index = source / previous_monomial_count;
                let monomial_index = source % previous_monomial_count;
                debug_assert!(point_index < point_count);
                let mut product_monomial = previous.monomials[monomial_index].clone();
                product_monomial[variable] += 1;
                let target_monomial_index = current.monomial_index[&product_monomial];
                let target = point_index * current_monomial_count + target_monomial_index;
                if let Some(coordinate_index) = current.coordinate_index_by_ambient_column[target] {
                    entries.push((coordinate_index, coeff));
                }
            }
            relations.push(sparse_vector(current.module_basis.len(), entries));
        }
    }

    relations
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

fn circular_edges(area: &[u8]) -> Option<Vec<CircularEdge>> {
    if !Graph::is_circular_unit_interval_area_sequence(area) {
        return None;
    }

    let n = area.len();
    let mut edges = Vec::new();
    for target in 0..n {
        for gap in 1..=usize::from(area[target]) {
            let source = (target + n - gap) % n;
            edges.push(CircularEdge {
                source,
                target,
                wrap: if source > target { 1 } else { 0 },
            });
        }
    }
    Some(edges)
}

fn reflect_affine_point(point: &AffinePoint, edge: CircularEdge) -> AffinePoint {
    let mut residues = point.residues.clone();
    let mut shifts = point.shifts.clone();

    let source_residue = residues[edge.source];
    let source_shift = shifts[edge.source];
    let target_residue = residues[edge.target];
    let target_shift = shifts[edge.target];

    residues[edge.source] = target_residue;
    shifts[edge.source] = target_shift + edge.wrap;
    residues[edge.target] = source_residue;
    shifts[edge.target] = source_shift - edge.wrap;

    AffinePoint { residues, shifts }
}

fn reflect_affine_point_mod_period(
    point: &AffinePoint,
    edge: CircularEdge,
    period: i32,
) -> AffinePoint {
    let mut reflected = reflect_affine_point(point, edge);
    for shift in &mut reflected.shifts {
        *shift = shift.rem_euclid(period);
    }
    reflected
}

fn affine_points(n: usize, radius: i32) -> Vec<AffinePoint> {
    let permutations = sym_poly_core::symmetric_group_permutation_basis(n);
    let shifts = shift_vectors(n, radius);
    let mut points = Vec::new();
    for residues in permutations {
        for shift in &shifts {
            points.push(AffinePoint {
                residues: residues.clone(),
                shifts: shift.clone(),
            });
        }
    }
    points.sort();
    points
}

fn affine_cyclic_points(n: usize, period: i32) -> Vec<AffinePoint> {
    let permutations = sym_poly_core::symmetric_group_permutation_basis(n);
    let shifts = cyclic_shift_vectors(n, period);
    let mut points = Vec::new();
    for residues in permutations {
        for shift in &shifts {
            points.push(AffinePoint {
                residues: residues.clone(),
                shifts: shift.clone(),
            });
        }
    }
    points.sort();
    points
}

fn shift_vectors(n: usize, radius: i32) -> Vec<Vec<i32>> {
    let mut result = Vec::new();
    let mut current = vec![0; n];
    shift_vectors_rec(0, radius, 0, &mut current, &mut result);
    result
}

fn cyclic_shift_vectors(n: usize, period: i32) -> Vec<Vec<i32>> {
    let mut result = Vec::new();
    let mut current = vec![0; n];
    cyclic_shift_vectors_rec(0, period, 0, &mut current, &mut result);
    result
}

fn cyclic_shift_vectors_rec(
    index: usize,
    period: i32,
    sum: i32,
    current: &mut [i32],
    result: &mut Vec<Vec<i32>>,
) {
    if index == current.len() {
        if sum.rem_euclid(period) == 0 {
            result.push(current.to_vec());
        }
        return;
    }

    for value in 0..period {
        current[index] = value;
        cyclic_shift_vectors_rec(index + 1, period, sum + value, current, result);
    }
}

fn shift_vectors_rec(
    index: usize,
    radius: i32,
    sum: i32,
    current: &mut [i32],
    result: &mut Vec<Vec<i32>>,
) {
    if index == current.len() {
        if sum == 0 {
            result.push(current.to_vec());
        }
        return;
    }

    for value in -radius..=radius {
        current[index] = value;
        shift_vectors_rec(index + 1, radius, sum + value, current, result);
    }
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

fn binomial_u32(n: u32, k: u32) -> u64 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result = 1u64;
    for i in 0..k {
        result = result * u64::from(n - i) / u64::from(i + 1);
    }
    result
}

fn frobenius_hilbert(
    frobenius: &BTreeMap<u32, sym_poly_sym::SymmetricFunction<i64>>,
) -> BTreeMap<u32, usize> {
    frobenius
        .iter()
        .filter_map(|(&degree, function)| {
            let schur = function.to_schur_basis();
            let dimension = schur
                .terms()
                .iter()
                .map(|(partition, &multiplicity)| {
                    let specht_dimension = partition
                        .count_syt()
                        .to_usize()
                        .expect("small-rank Specht dimension fits in usize");
                    usize::try_from(multiplicity).expect("multiplicity is nonnegative")
                        * specht_dimension
                })
                .sum::<usize>();
            (dimension != 0).then_some((degree, dimension))
        })
        .collect()
}

fn format_hilbert(hilbert: &BTreeMap<u32, usize>) -> String {
    if hilbert.is_empty() {
        return "0".to_string();
    }

    hilbert
        .iter()
        .map(|(&degree, &coefficient)| match (coefficient, degree) {
            (1, 0) => "1".to_string(),
            (c, 0) => c.to_string(),
            (1, 1) => "q".to_string(),
            (c, 1) => format!("{c}q"),
            (1, d) => format!("q^{d}"),
            (c, d) => format!("{c}q^{d}"),
        })
        .collect::<Vec<_>>()
        .join(" + ")
}

fn format_model(model: AffineModel) -> String {
    match model {
        AffineModel::Window { radius } => format!("window radius {radius}"),
        AffineModel::Cyclic { period } => format!("cyclic period {period}"),
    }
}

fn all_circular_area_sequences(n: usize) -> Vec<Vec<u8>> {
    if n == 0 {
        return vec![Vec::new()];
    }

    let mut result = Vec::new();
    let mut current = vec![0u8; n];
    circular_area_sequences_rec(0, &mut current, &mut result);
    result
}

fn circular_area_sequences_rec(index: usize, current: &mut [u8], result: &mut Vec<Vec<u8>>) {
    if index == current.len() {
        if Graph::is_circular_unit_interval_area_sequence(current) {
            result.push(current.to_vec());
        }
        return;
    }

    for value in 0..current.len() {
        current[index] = value as u8;
        circular_area_sequences_rec(index + 1, current, result);
    }
}
