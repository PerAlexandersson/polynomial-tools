//! Stembridge peak quasisymmetric functions.
//!
//! The peak function indexed by a peak set `Lambda` in degree `n` is returned
//! in the fundamental basis:
//!
//! `K_Lambda = 2^(|Lambda| + 1) sum F_alpha`,
//!
//! where the sum is over compositions `alpha` of `n` such that
//! `Lambda ⊆ Des(alpha) △ (Des(alpha) + 1)`.

use std::collections::{BTreeMap, BTreeSet};

use sym_poly_core::{Composition, Ring};

use crate::basis::QSymBasis;
use crate::qsym_function::QSymFunction;

/// A standard peak composition tableau.
///
/// Rows are stored bottom-to-top, matching the convention in the peak
/// composition tableau literature.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StandardPeakCompositionTableau {
    pub rows: Vec<Vec<u32>>,
}

/// A standard peak Young composition tableau.
///
/// Rows are stored bottom-to-top, matching the convention in the peak
/// composition tableau literature.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StandardPeakYoungCompositionTableau {
    pub rows: Vec<Vec<u32>>,
}

impl StandardPeakCompositionTableau {
    /// Enumerate standard peak composition tableaux of peak shape `alpha`.
    pub fn enumerate(alpha: &[u32]) -> Vec<StandardPeakCompositionTableau> {
        assert!(
            is_peak_composition(alpha),
            "invalid peak composition {:?}",
            alpha
        );

        let row_count = alpha.len();
        let n: u32 = alpha.iter().sum();
        if n == 0 {
            return vec![StandardPeakCompositionTableau { rows: Vec::new() }];
        }

        let mut grid: Vec<Vec<u32>> = alpha.iter().map(|&a| vec![0; a as usize]).collect();
        let cells: Vec<(usize, usize)> = (0..row_count)
            .flat_map(|r| (0..alpha[r] as usize).map(move |c| (r, c)))
            .collect();
        let mut used = vec![false; n as usize + 1];
        let mut results = Vec::new();

        enumerate_spct(&mut grid, &cells, 0, n, &mut used, &mut results);
        results
    }

    /// Descent set `Des_up(T) = {i : i + 1 is strictly above i}`.
    pub fn upward_descent_set(&self) -> BTreeSet<u32> {
        let n: u32 = self.rows.iter().map(|r| r.len() as u32).sum();
        let mut entry_row = vec![0usize; n as usize + 1];
        for (r, row) in self.rows.iter().enumerate() {
            for &value in row {
                entry_row[value as usize] = r;
            }
        }

        let mut descents = BTreeSet::new();
        for i in 1..n {
            if entry_row[(i + 1) as usize] > entry_row[i as usize] {
                descents.insert(i);
            }
        }
        descents
    }

    /// The peak composition `comp_n(Peak(Des_up(T)))`.
    pub fn upward_peak_composition(&self) -> Composition {
        let n: u32 = self.rows.iter().map(|r| r.len() as u32).sum();
        peak_composition_from_descent_set(&self.upward_descent_set(), n)
    }
}

impl StandardPeakYoungCompositionTableau {
    /// Enumerate standard peak Young composition tableaux of peak shape `alpha`.
    pub fn enumerate(alpha: &[u32]) -> Vec<StandardPeakYoungCompositionTableau> {
        assert!(
            is_peak_composition(alpha),
            "invalid peak composition {:?}",
            alpha
        );

        let row_count = alpha.len();
        let n: u32 = alpha.iter().sum();
        if n == 0 {
            return vec![StandardPeakYoungCompositionTableau { rows: Vec::new() }];
        }

        let mut grid: Vec<Vec<u32>> = alpha.iter().map(|&a| vec![0; a as usize]).collect();
        let cells: Vec<(usize, usize)> = (0..row_count)
            .flat_map(|r| (0..alpha[r] as usize).map(move |c| (r, c)))
            .collect();
        let mut used = vec![false; n as usize + 1];
        let mut results = Vec::new();

        enumerate_spyct(&mut grid, &cells, 0, n, &mut used, &mut results);
        results
    }

    /// Descent set `Des_left(T) = {i : i + 1 is weakly left of i}`.
    pub fn left_descent_set(&self) -> BTreeSet<u32> {
        let n: u32 = self.rows.iter().map(|r| r.len() as u32).sum();
        let mut entry_col = vec![0usize; n as usize + 1];
        for row in &self.rows {
            for (c, &value) in row.iter().enumerate() {
                entry_col[value as usize] = c;
            }
        }

        let mut descents = BTreeSet::new();
        for i in 1..n {
            if entry_col[(i + 1) as usize] <= entry_col[i as usize] {
                descents.insert(i);
            }
        }
        descents
    }

    /// The peak composition `comp_n(Peak(Des_left(T)))`.
    pub fn left_peak_composition(&self) -> Composition {
        let n: u32 = self.rows.iter().map(|r| r.len() as u32).sum();
        peak_composition_from_descent_set(&self.left_descent_set(), n)
    }
}

/// Return whether `peak_set` is a valid peak set in degree `degree`.
///
/// A peak set is a subset of `{2, ..., degree - 1}` with no adjacent entries.
pub fn is_peak_set(peak_set: &[u32], degree: u32) -> bool {
    let set: BTreeSet<u32> = peak_set.iter().copied().collect();
    if set.len() != peak_set.len() {
        return false;
    }
    set.iter().all(|&p| 2 <= p && p < degree) && set.iter().all(|&p| !set.contains(&(p + 1)))
}

/// Return whether `alpha` is a peak composition.
///
/// A peak composition has all parts except possibly the last greater than 1.
pub fn is_peak_composition(alpha: &[u32]) -> bool {
    !alpha.is_empty()
        && alpha.iter().all(|&part| part > 0)
        && alpha
            .iter()
            .take(alpha.len().saturating_sub(1))
            .all(|&part| part > 1)
}

/// Stembridge peak quasisymmetric function in the fundamental basis.
pub fn peak_quasisymmetric<C: Ring>(peak_set: &[u32], degree: u32) -> QSymFunction<C> {
    assert!(
        is_peak_set(peak_set, degree),
        "invalid peak set {:?} in degree {}",
        peak_set,
        degree
    );

    let peak_set: BTreeSet<u32> = peak_set.iter().copied().collect();
    let coeff = peak_coefficient::<C>(peak_set.len());
    let mut terms = BTreeMap::new();

    for alpha in Composition::integer_compositions(degree) {
        if peak_set.is_subset(&descent_symmetric_difference(&alpha)) {
            terms.insert(alpha, coeff.clone());
        }
    }

    QSymFunction::from_terms(QSymBasis::Fundamental, terms)
}

/// Quasisymmetric Schur Q-function in the peak-function basis.
///
/// The keys are peak compositions indexing Stembridge peak functions.
pub fn qsym_schur_q_peak_expansion<C: Ring>(alpha: &[u32]) -> BTreeMap<Composition, C> {
    let tableaux = StandardPeakCompositionTableau::enumerate(alpha);
    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();

    for tableau in tableaux {
        let peak_comp = tableau.upward_peak_composition();
        let entry = terms.entry(peak_comp).or_insert_with(C::zero);
        *entry = entry.clone() + C::one();
    }

    terms
}

/// Quasisymmetric Schur Q-function in the fundamental basis.
pub fn qsym_schur_q<C: Ring>(alpha: &[u32]) -> QSymFunction<C> {
    let degree: u32 = alpha.iter().sum();
    let mut result = QSymFunction::zero(QSymBasis::Fundamental);

    for (peak_comp, coeff) in qsym_schur_q_peak_expansion::<C>(alpha) {
        let peak_set: Vec<u32> = peak_comp.composition_to_descent_set().into_iter().collect();
        result = result + peak_quasisymmetric::<C>(&peak_set, degree).scale(&coeff);
    }

    result
}

/// Peak Young quasisymmetric Schur function in the peak-function basis.
///
/// The keys are peak compositions indexing Stembridge peak functions.
pub fn peak_young_qsym_schur_peak_expansion<C: Ring>(alpha: &[u32]) -> BTreeMap<Composition, C> {
    let tableaux = StandardPeakYoungCompositionTableau::enumerate(alpha);
    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();

    for tableau in tableaux {
        let peak_comp = tableau.left_peak_composition();
        let entry = terms.entry(peak_comp).or_insert_with(C::zero);
        *entry = entry.clone() + C::one();
    }

    terms
}

/// Peak Young quasisymmetric Schur function in the fundamental basis.
pub fn peak_young_qsym_schur<C: Ring>(alpha: &[u32]) -> QSymFunction<C> {
    let degree: u32 = alpha.iter().sum();
    let mut result = QSymFunction::zero(QSymBasis::Fundamental);

    for (peak_comp, coeff) in peak_young_qsym_schur_peak_expansion::<C>(alpha) {
        let peak_set: Vec<u32> = peak_comp.composition_to_descent_set().into_iter().collect();
        result = result + peak_quasisymmetric::<C>(&peak_set, degree).scale(&coeff);
    }

    result
}

/// Convert a descent set to the corresponding peak set.
pub fn peak_set_from_descent_set(descent_set: &BTreeSet<u32>) -> BTreeSet<u32> {
    descent_set
        .iter()
        .filter(|&&d| d != 1 && !descent_set.contains(&(d - 1)))
        .copied()
        .collect()
}

/// Convert a descent set to the corresponding peak composition of degree `n`.
pub fn peak_composition_from_descent_set(descent_set: &BTreeSet<u32>, n: u32) -> Composition {
    let peak_set = peak_set_from_descent_set(descent_set);
    Composition::from_descent_set(&peak_set, n)
}

fn peak_coefficient<C: Ring>(peak_set_size: usize) -> C {
    let exponent = u32::try_from(peak_set_size + 1).expect("peak set size overflow");
    let value = 2_i64
        .checked_pow(exponent)
        .expect("peak coefficient overflow");
    C::from_i64(value)
}

fn descent_symmetric_difference(alpha: &Composition) -> BTreeSet<u32> {
    let descents = alpha.composition_to_descent_set();
    let shifted: BTreeSet<u32> = descents.iter().map(|&d| d + 1).collect();
    descents.symmetric_difference(&shifted).copied().collect()
}

fn enumerate_spct(
    grid: &mut [Vec<u32>],
    cells: &[(usize, usize)],
    idx: usize,
    n: u32,
    used: &mut [bool],
    results: &mut Vec<StandardPeakCompositionTableau>,
) {
    if idx == cells.len() {
        if satisfies_peak_subdiagram_condition(grid, n) {
            results.push(StandardPeakCompositionTableau {
                rows: grid.to_vec(),
            });
        }
        return;
    }

    let (r, c) = cells[idx];

    for value in 1..=n {
        if used[value as usize] {
            continue;
        }
        if c > 0 && value <= grid[r][c - 1] {
            continue;
        }
        if c == 0 && r > 0 && value <= grid[r - 1][0] {
            continue;
        }

        used[value as usize] = true;
        grid[r][c] = value;
        enumerate_spct(grid, cells, idx + 1, n, used, results);
        grid[r][c] = 0;
        used[value as usize] = false;
    }
}

fn enumerate_spyct(
    grid: &mut [Vec<u32>],
    cells: &[(usize, usize)],
    idx: usize,
    n: u32,
    used: &mut [bool],
    results: &mut Vec<StandardPeakYoungCompositionTableau>,
) {
    if idx == cells.len() {
        if satisfies_peak_subdiagram_condition(grid, n) && satisfies_spyct_triple_rule(grid) {
            results.push(StandardPeakYoungCompositionTableau {
                rows: grid.to_vec(),
            });
        }
        return;
    }

    let (r, c) = cells[idx];

    for value in 1..=n {
        if used[value as usize] {
            continue;
        }
        if c > 0 && value <= grid[r][c - 1] {
            continue;
        }
        if c == 0 && r > 0 && value <= grid[r - 1][0] {
            continue;
        }

        used[value as usize] = true;
        grid[r][c] = value;
        enumerate_spyct(grid, cells, idx + 1, n, used, results);
        grid[r][c] = 0;
        used[value as usize] = false;
    }
}

fn satisfies_peak_subdiagram_condition(grid: &[Vec<u32>], n: u32) -> bool {
    (1..=n).all(|k| subdiagram_shape(grid, k).is_some_and(|shape| is_peak_composition(&shape)))
}

fn subdiagram_shape(grid: &[Vec<u32>], max_value: u32) -> Option<Vec<u32>> {
    let mut shape = Vec::with_capacity(grid.len());
    for row in grid {
        let count = row.iter().take_while(|&&value| value <= max_value).count();
        if row.iter().skip(count).any(|&value| value <= max_value) {
            return None;
        }
        shape.push(count as u32);
    }

    while shape.last() == Some(&0) {
        shape.pop();
    }

    if shape.iter().any(|&part| part == 0) {
        return None;
    }

    Some(shape)
}

fn satisfies_spyct_triple_rule(grid: &[Vec<u32>]) -> bool {
    for r in 0..grid.len() {
        for c in 0..grid[r].len() {
            for lower in 0..r {
                if grid[lower].len() <= c + 1 {
                    continue;
                }
                if grid[r][c] < grid[lower][c + 1] {
                    if grid[r].len() <= c + 1 {
                        return false;
                    }
                    if grid[r][c + 1] >= grid[lower][c + 1] {
                        return false;
                    }
                }
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set(values: &[u32]) -> BTreeSet<u32> {
        values.iter().copied().collect()
    }

    #[test]
    fn test_peak_set_validation() {
        assert!(is_peak_set(&[], 4));
        assert!(is_peak_set(&[2], 4));
        assert!(is_peak_set(&[3], 4));
        assert!(!is_peak_set(&[1], 4));
        assert!(!is_peak_set(&[4], 4));
        assert!(!is_peak_set(&[2, 3], 5));
        assert!(!is_peak_set(&[2, 2], 5));
    }

    #[test]
    fn test_peak_composition_validation() {
        assert!(is_peak_composition(&[4]));
        assert!(is_peak_composition(&[2, 2]));
        assert!(is_peak_composition(&[3, 1]));
        assert!(!is_peak_composition(&[]));
        assert!(!is_peak_composition(&[1, 3]));
        assert!(!is_peak_composition(&[2, 1, 1]));
    }

    #[test]
    fn test_peak_quasisymmetric_degree_three_empty() {
        let k_empty: QSymFunction<i64> = peak_quasisymmetric(&[], 3);

        assert_eq!(k_empty.coefficient(&Composition::new(vec![1, 1, 1])), 2);
        assert_eq!(k_empty.coefficient(&Composition::new(vec![1, 2])), 2);
        assert_eq!(k_empty.coefficient(&Composition::new(vec![2, 1])), 2);
        assert_eq!(k_empty.coefficient(&Composition::new(vec![3])), 2);
        assert_eq!(k_empty.terms().len(), 4);
    }

    #[test]
    fn test_peak_quasisymmetric_degree_three_peak_two() {
        let k_two: QSymFunction<i64> = peak_quasisymmetric(&[2], 3);

        assert_eq!(k_two.coefficient(&Composition::new(vec![1, 2])), 4);
        assert_eq!(k_two.coefficient(&Composition::new(vec![2, 1])), 4);
        assert_eq!(k_two.terms().len(), 2);
    }

    #[test]
    fn test_standard_peak_composition_tableaux_source_example() {
        // Searles--Slattery-Holmes, Example 3.17, lists four SPCTs of
        // shape (3,3) with peak compositions (3,3), (2,2,2), (2,3,1),
        // and (2,4).
        let tableaux = StandardPeakCompositionTableau::enumerate(&[3, 3]);
        let peak_comps: Vec<Composition> = tableaux
            .iter()
            .map(StandardPeakCompositionTableau::upward_peak_composition)
            .collect();

        assert_eq!(tableaux.len(), 4);
        assert!(peak_comps.contains(&Composition::new(vec![3, 3])));
        assert!(peak_comps.contains(&Composition::new(vec![2, 2, 2])));
        assert!(peak_comps.contains(&Composition::new(vec![2, 3, 1])));
        assert!(peak_comps.contains(&Composition::new(vec![2, 4])));
    }

    #[test]
    fn test_qsym_schur_q_peak_expansion_source_example() {
        // Q~_(3,3) = K_(3,3) + K_(2,2,2) + K_(2,3,1) + K_(2,4).
        let expansion = qsym_schur_q_peak_expansion::<i64>(&[3, 3]);

        assert_eq!(expansion.get(&Composition::new(vec![3, 3])), Some(&1));
        assert_eq!(expansion.get(&Composition::new(vec![2, 2, 2])), Some(&1));
        assert_eq!(expansion.get(&Composition::new(vec![2, 3, 1])), Some(&1));
        assert_eq!(expansion.get(&Composition::new(vec![2, 4])), Some(&1));
        assert_eq!(expansion.len(), 4);
    }

    #[test]
    fn test_standard_peak_young_composition_tableaux_source_example() {
        // Searles--Slattery-Holmes, Example 4.4, lists three SPYCTs of
        // shape (3,3) with peak compositions (3,3), (2,2,2), and (2,3,1).
        let tableaux = StandardPeakYoungCompositionTableau::enumerate(&[3, 3]);
        let peak_comps: Vec<Composition> = tableaux
            .iter()
            .map(StandardPeakYoungCompositionTableau::left_peak_composition)
            .collect();

        assert_eq!(tableaux.len(), 3);
        assert!(peak_comps.contains(&Composition::new(vec![3, 3])));
        assert!(peak_comps.contains(&Composition::new(vec![2, 2, 2])));
        assert!(peak_comps.contains(&Composition::new(vec![2, 3, 1])));
    }

    #[test]
    fn test_peak_young_qsym_schur_peak_expansion_source_example() {
        // S~_(3,3) = K_(3,3) + K_(2,2,2) + K_(2,3,1).
        let expansion = peak_young_qsym_schur_peak_expansion::<i64>(&[3, 3]);

        assert_eq!(expansion.get(&Composition::new(vec![3, 3])), Some(&1));
        assert_eq!(expansion.get(&Composition::new(vec![2, 2, 2])), Some(&1));
        assert_eq!(expansion.get(&Composition::new(vec![2, 3, 1])), Some(&1));
        assert_eq!(expansion.len(), 3);
    }

    #[test]
    fn test_site_example_peak_shape_431_spct() {
        // Site example: rows are bottom-to-top.
        let tableau = StandardPeakCompositionTableau {
            rows: vec![vec![1, 2, 3, 4], vec![5, 6, 7], vec![8]],
        };

        assert!(StandardPeakCompositionTableau::enumerate(&[4, 3, 1]).contains(&tableau));
        assert_eq!(tableau.upward_descent_set(), set(&[4, 7]));
        assert_eq!(
            tableau.upward_peak_composition(),
            Composition::new(vec![4, 3, 1])
        );
    }

    #[test]
    fn test_site_example_peak_shape_431_spyct() {
        // The same filling also satisfies the peak Young triple rule.
        let tableau = StandardPeakYoungCompositionTableau {
            rows: vec![vec![1, 2, 3, 4], vec![5, 6, 7], vec![8]],
        };

        assert!(StandardPeakYoungCompositionTableau::enumerate(&[4, 3, 1]).contains(&tableau));
        assert_eq!(tableau.left_descent_set(), set(&[4, 7]));
        assert_eq!(
            tableau.left_peak_composition(),
            Composition::new(vec![4, 3, 1])
        );
    }
}
