//! Quasi-symmetric Schur functions and related bases.
//!
//! Implements several composition-indexed bases of QSym:
//!
//! - **Composition tableaux** and their enumeration.
//! - **QSym Schur functions** S_α = Σ F_{Des(T)} over composition tableaux T
//!   (Haglund–Luoto–Mason–van Willigenburg).
//! - **Standard reverse composition tableaux** for the fundamental expansion
//!   of quasisymmetric Schur functions.
//! - **Young quasisymmetric Schur functions** via the rho involution and
//!   reverse-index relation of Luoto--Mykytiuk--van Willigenburg.
//! - **Row-strict quasisymmetric Schur functions** via the omega involution.
//! - **Dual immaculate functions** S\*_α = Σ F_{Des(T)} over standard
//!   immaculate tableaux T (Berg–Bergeron–Saliola–Serrano–Zabrocki).
//! - **Row-strict dual immaculate functions** RS\*_α, with strictly
//!   increasing rows and weakly increasing first column.
//! - **Fundamental slide polynomials** 𝔉_α (Assaf–Searles).

use std::collections::BTreeMap;

use sym_poly_core::{Composition, Ring};

use crate::basis::QSymBasis;
use crate::qsym_function::QSymFunction;

// ── Composition tableaux ────────────────────────────────────────────

/// A composition tableau of shape α is a filling of the composition
/// diagram of α such that:
/// - Entries in row i belong to {1, …, n} for some n.
/// - Each row is weakly increasing left-to-right.
/// - The leftmost column is strictly increasing top-to-bottom.
/// - The "triple rule" holds: for entries in the same column, if the
///   entry below is ≤ the entry to the upper-left, that's forbidden.
///
/// We represent as `rows[i]` = entries in row i (0-indexed, top row = 0).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompositionTableau {
    pub rows: Vec<Vec<u32>>,
}

impl CompositionTableau {
    /// Enumerate all composition tableaux of shape α with entries in `{1, …, max_entry}`.
    pub fn enumerate(alpha: &[u32], max_entry: u32) -> Vec<CompositionTableau> {
        let n = alpha.len();
        if n == 0 {
            return vec![CompositionTableau { rows: Vec::new() }];
        }
        let mut grid: Vec<Vec<u32>> = alpha.iter().map(|&a| vec![0; a as usize]).collect();
        let mut results = Vec::new();

        // Fill cells: row-by-row, left-to-right
        let cells: Vec<(usize, usize)> = (0..n)
            .flat_map(|r| (0..alpha[r] as usize).map(move |c| (r, c)))
            .collect();

        enumerate_ct(&mut grid, &cells, 0, max_entry, alpha, &mut results);
        results
    }

    /// Descent set: position *i* (1-indexed) is a descent when *i+1*
    /// appears in a strictly lower row than *i*.
    ///
    /// For semistandard fillings, positions are entry values (not cell indices).
    pub fn descent_set(&self) -> Vec<u32> {
        let max_entry = self
            .rows
            .iter()
            .flat_map(|r| r.iter())
            .copied()
            .max()
            .unwrap_or(0);
        if max_entry <= 1 {
            return Vec::new();
        }
        // entry → lowest row containing it (last occurrence, since rows are increasing)
        let mut entry_row = vec![0usize; max_entry as usize + 1];
        for (r, row) in self.rows.iter().enumerate() {
            for &v in row {
                entry_row[v as usize] = r;
            }
        }
        let mut descents = Vec::new();
        for i in 1..max_entry {
            if entry_row[i as usize] < entry_row[(i + 1) as usize] {
                descents.push(i);
            }
        }
        descents
    }

    /// Descent composition: converts the descent set into a composition
    /// of max_entry (the largest value in the tableau).
    pub fn descent_composition(&self) -> Composition {
        let max_entry = self
            .rows
            .iter()
            .flat_map(|r| r.iter())
            .copied()
            .max()
            .unwrap_or(0);
        if max_entry == 0 {
            return Composition::new(Vec::new());
        }
        let des = self.descent_set();
        descent_set_to_composition(&des, max_entry)
    }
}

fn enumerate_ct(
    grid: &mut [Vec<u32>],
    cells: &[(usize, usize)],
    idx: usize,
    max_entry: u32,
    alpha: &[u32],
    results: &mut Vec<CompositionTableau>,
) {
    if idx == cells.len() {
        results.push(CompositionTableau {
            rows: grid.to_vec(),
        });
        return;
    }

    let (r, c) = cells[idx];

    // Determine bounds
    let min_val = if c > 0 {
        grid[r][c - 1] // weakly increasing in rows
    } else if r > 0 {
        // First column: strictly increasing top-to-bottom
        grid[r - 1][0] + 1
    } else {
        1
    };

    for v in min_val..=max_entry {
        // Composition tableau triple rule (HLMVW): for each row r' above r
        // that has column c, the filling is invalid if T(r,c) ≤ T(r',c-1).
        // Here T(r,c) = v (the candidate value).
        let mut valid = true;
        if c > 0 {
            for rp in 0..r {
                if (alpha[rp] as usize) > c && grid[rp][c - 1] > 0 && v <= grid[rp][c - 1] {
                    valid = false;
                    break;
                }
            }
        }

        if valid {
            grid[r][c] = v;
            enumerate_ct(grid, cells, idx + 1, max_entry, alpha, results);
        }
    }
    grid[r][c] = 0;
}

/// A standard reverse composition tableau of shape α.
///
/// Rows are weakly decreasing left-to-right, the first column is strictly
/// increasing top-to-bottom, and the reverse-composition triple rule holds.
/// The entries are `1, …, |alpha|`, each exactly once.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReverseCompositionTableau {
    pub rows: Vec<Vec<u32>>,
}

impl ReverseCompositionTableau {
    /// Enumerate all standard reverse composition tableaux of shape α.
    pub fn enumerate_standard(alpha: &[u32]) -> Vec<ReverseCompositionTableau> {
        let k = alpha.len();
        let n: u32 = alpha.iter().sum();
        if n == 0 {
            return vec![ReverseCompositionTableau {
                rows: vec![Vec::new(); k],
            }];
        }

        let mut grid: Vec<Vec<u32>> = alpha.iter().map(|&a| vec![0; a as usize]).collect();
        let cells: Vec<(usize, usize)> = (0..k)
            .flat_map(|r| (0..alpha[r] as usize).map(move |c| (r, c)))
            .collect();
        let mut used = vec![false; n as usize + 1];
        let mut results = Vec::new();

        enumerate_srct(&mut grid, &cells, 0, n, &mut used, alpha, &mut results);
        results
    }

    /// Descent set: `i` is a descent when `i + 1` appears weakly to the
    /// right of `i`.
    pub fn descent_set(&self) -> Vec<u32> {
        let n: u32 = self.rows.iter().map(|r| r.len() as u32).sum();
        if n <= 1 {
            return Vec::new();
        }

        let mut entry_col = vec![0usize; n as usize + 1];
        for row in &self.rows {
            for (c, &v) in row.iter().enumerate() {
                entry_col[v as usize] = c;
            }
        }

        let mut descents = Vec::new();
        for i in 1..n {
            if entry_col[(i + 1) as usize] >= entry_col[i as usize] {
                descents.push(i);
            }
        }
        descents
    }

    /// Convert descent set to composition.
    pub fn descent_composition(&self) -> Composition {
        let n: u32 = self.rows.iter().map(|r| r.len() as u32).sum();
        descent_set_to_composition(&self.descent_set(), n)
    }
}

fn enumerate_srct(
    grid: &mut [Vec<u32>],
    cells: &[(usize, usize)],
    idx: usize,
    n: u32,
    used: &mut [bool],
    alpha: &[u32],
    results: &mut Vec<ReverseCompositionTableau>,
) {
    if idx == cells.len() {
        results.push(ReverseCompositionTableau {
            rows: grid.to_vec(),
        });
        return;
    }

    let (r, c) = cells[idx];

    for v in 1..=n {
        if used[v as usize] {
            continue;
        }
        // Rows are weakly decreasing.
        if c > 0 && v > grid[r][c - 1] {
            continue;
        }
        // First column is strictly increasing top-to-bottom.
        if c == 0 && r > 0 && v <= grid[r - 1][0] {
            continue;
        }
        if !valid_srct_triples_for_cell(grid, alpha, r, c, v) {
            continue;
        }

        used[v as usize] = true;
        grid[r][c] = v;
        enumerate_srct(grid, cells, idx + 1, n, used, alpha, results);
        grid[r][c] = 0;
        used[v as usize] = false;
    }
}

fn valid_srct_triples_for_cell(
    grid: &[Vec<u32>],
    alpha: &[u32],
    r: usize,
    c: usize,
    v: u32,
) -> bool {
    if c == 0 {
        return true;
    }

    // For each upper row with a box b immediately left of this column, the
    // reverse-composition triple rule says: if a <= b, then c must exist and
    // a < c.
    for upper in 0..r {
        if alpha[upper] as usize <= c - 1 {
            continue;
        }
        let b = grid[upper][c - 1];
        if b == 0 || v > b {
            continue;
        }
        if alpha[upper] as usize <= c {
            return false;
        }
        let upper_right = grid[upper][c];
        if upper_right == 0 || v >= upper_right {
            return false;
        }
    }

    true
}

// ── QSym Schur function ─────────────────────────────────────────────

/// Quasi-symmetric Schur function S_α in the fundamental basis.
///
/// S_α = Σ_{T ∈ SRCT(α)} F_{Des(T)}
///
/// where SRCT(α) is the set of standard reverse composition tableaux of
/// shape α. The `_max_entry` argument is retained for backward compatibility.
pub fn qsym_schur<C: Ring>(alpha: &[u32], _max_entry: u32) -> QSymFunction<C> {
    let tableaux = ReverseCompositionTableau::enumerate_standard(alpha);
    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();
    for t in &tableaux {
        let comp = t.descent_composition();
        let entry = terms.entry(comp).or_insert_with(C::zero);
        *entry = entry.clone() + C::one();
    }
    QSymFunction::from_terms(QSymBasis::Fundamental, terms)
}

/// Young quasisymmetric Schur function in the fundamental basis.
///
/// Uses the identity
///
/// `YQS_alpha = rho(S_rev(alpha))`,
///
/// where `S` is the quasisymmetric Schur function and `rho(F_alpha) =
/// F_rev(alpha)`.
pub fn young_qsym_schur<C: Ring>(alpha: &[u32]) -> QSymFunction<C> {
    let rev_alpha: Vec<u32> = alpha.iter().rev().copied().collect();
    qsym_schur::<C>(&rev_alpha, alpha.iter().sum()).rho_involution()
}

/// Row-strict quasisymmetric Schur function in the fundamental basis.
///
/// Uses the identity `RS_alpha = omega(S_alpha)`, where `S_alpha` is the
/// quasisymmetric Schur function.
pub fn row_strict_qsym_schur<C: Ring>(alpha: &[u32]) -> QSymFunction<C> {
    qsym_schur::<C>(alpha, alpha.iter().sum()).omega_involution()
}

// ── Dual immaculate quasisymmetric functions ────────────────────────

/// A standard immaculate tableau of shape α is a filling of the
/// composition diagram with entries 1, …, |α| (each exactly once) such that:
/// - Each row is strictly increasing (left to right).
/// - The leftmost column is strictly increasing (top to bottom).
/// - **No other column constraint** (unlike SYT or composition tableaux).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ImmaculateTableau {
    pub rows: Vec<Vec<u32>>,
}

impl ImmaculateTableau {
    /// Enumerate all standard immaculate tableaux of shape α.
    pub fn enumerate(alpha: &[u32]) -> Vec<ImmaculateTableau> {
        let k = alpha.len();
        let n: u32 = alpha.iter().sum();
        if n == 0 {
            return vec![ImmaculateTableau {
                rows: vec![Vec::new(); k],
            }];
        }

        let mut grid: Vec<Vec<u32>> = alpha.iter().map(|&a| vec![0; a as usize]).collect();
        let cells: Vec<(usize, usize)> = (0..k)
            .flat_map(|r| (0..alpha[r] as usize).map(move |c| (r, c)))
            .collect();
        let mut used = vec![false; n as usize + 1];
        let mut results = Vec::new();

        enumerate_immaculate(&mut grid, &cells, 0, n, &mut used, alpha, &mut results);
        results
    }

    /// Enumerate semistandard immaculate tableaux of shape α with entries in {1,…,max_entry}.
    ///
    /// Rows are weakly increasing; first column is strictly increasing.
    /// No other column constraint.
    pub fn enumerate_semistandard(alpha: &[u32], max_entry: u32) -> Vec<ImmaculateTableau> {
        let k = alpha.len();
        if max_entry == 0 {
            return if alpha.iter().all(|&a| a == 0) {
                vec![ImmaculateTableau {
                    rows: vec![Vec::new(); k],
                }]
            } else {
                Vec::new()
            };
        }

        let mut grid: Vec<Vec<u32>> = alpha.iter().map(|&a| vec![0; a as usize]).collect();
        let cells: Vec<(usize, usize)> = (0..k)
            .flat_map(|r| (0..alpha[r] as usize).map(move |c| (r, c)))
            .collect();
        let mut results = Vec::new();

        enumerate_ss_immaculate(&mut grid, &cells, 0, max_entry, alpha, &mut results);
        results
    }

    /// Enumerate row-strict semistandard immaculate tableaux of shape α.
    ///
    /// Rows are strictly increasing; first column is weakly increasing.
    /// No other column constraint.
    pub fn enumerate_row_strict_semistandard(
        alpha: &[u32],
        max_entry: u32,
    ) -> Vec<ImmaculateTableau> {
        let k = alpha.len();
        if max_entry == 0 {
            return if alpha.iter().all(|&a| a == 0) {
                vec![ImmaculateTableau {
                    rows: vec![Vec::new(); k],
                }]
            } else {
                Vec::new()
            };
        }

        let mut grid: Vec<Vec<u32>> = alpha.iter().map(|&a| vec![0; a as usize]).collect();
        let cells: Vec<(usize, usize)> = (0..k)
            .flat_map(|r| (0..alpha[r] as usize).map(move |c| (r, c)))
            .collect();
        let mut results = Vec::new();

        enumerate_row_strict_ss_immaculate(&mut grid, &cells, 0, max_entry, &mut results);
        results
    }

    /// Descent set of a standard immaculate tableau.
    /// Position i is a descent if i+1 appears in a strictly lower row.
    pub fn descent_set(&self) -> Vec<u32> {
        let n: u32 = self.rows.iter().map(|r| r.len() as u32).sum();
        if n <= 1 {
            return Vec::new();
        }
        let mut entry_row = vec![0usize; n as usize + 1];
        for (r, row) in self.rows.iter().enumerate() {
            for &v in row {
                entry_row[v as usize] = r;
            }
        }
        let mut descents = Vec::new();
        for i in 1..n {
            if entry_row[i as usize] < entry_row[(i + 1) as usize] {
                descents.push(i);
            }
        }
        descents
    }

    /// Convert descent set to composition.
    pub fn descent_composition(&self) -> Composition {
        let n: u32 = self.rows.iter().map(|r| r.len() as u32).sum();
        descent_set_to_composition(&self.descent_set(), n)
    }
}

fn enumerate_immaculate(
    grid: &mut [Vec<u32>],
    cells: &[(usize, usize)],
    idx: usize,
    n: u32,
    used: &mut [bool],
    _alpha: &[u32],
    results: &mut Vec<ImmaculateTableau>,
) {
    if idx == cells.len() {
        results.push(ImmaculateTableau {
            rows: grid.to_vec(),
        });
        return;
    }

    let (r, c) = cells[idx];

    for v in 1..=n {
        if used[v as usize] {
            continue;
        }
        // Row constraint: strictly increasing (v > left neighbor)
        if c > 0 && v <= grid[r][c - 1] {
            continue;
        }
        // First column: strictly increasing top-to-bottom
        if c == 0 && r > 0 && v <= grid[r - 1][0] {
            continue;
        }
        used[v as usize] = true;
        grid[r][c] = v;
        enumerate_immaculate(grid, cells, idx + 1, n, used, _alpha, results);
        grid[r][c] = 0;
        used[v as usize] = false;
    }
}

fn enumerate_ss_immaculate(
    grid: &mut [Vec<u32>],
    cells: &[(usize, usize)],
    idx: usize,
    max_entry: u32,
    _alpha: &[u32],
    results: &mut Vec<ImmaculateTableau>,
) {
    if idx == cells.len() {
        results.push(ImmaculateTableau {
            rows: grid.to_vec(),
        });
        return;
    }

    let (r, c) = cells[idx];

    let min_val = if c > 0 {
        grid[r][c - 1] // weakly increasing rows
    } else if r > 0 {
        grid[r - 1][0] + 1 // strictly increasing first column
    } else {
        1
    };

    for v in min_val..=max_entry {
        grid[r][c] = v;
        enumerate_ss_immaculate(grid, cells, idx + 1, max_entry, _alpha, results);
    }
    grid[r][c] = 0;
}

fn enumerate_row_strict_ss_immaculate(
    grid: &mut [Vec<u32>],
    cells: &[(usize, usize)],
    idx: usize,
    max_entry: u32,
    results: &mut Vec<ImmaculateTableau>,
) {
    if idx == cells.len() {
        results.push(ImmaculateTableau {
            rows: grid.to_vec(),
        });
        return;
    }

    let (r, c) = cells[idx];

    let min_val = if c > 0 {
        grid[r][c - 1] + 1 // strictly increasing rows
    } else if r > 0 {
        grid[r - 1][0] // weakly increasing first column
    } else {
        1
    };

    for v in min_val..=max_entry {
        grid[r][c] = v;
        enumerate_row_strict_ss_immaculate(grid, cells, idx + 1, max_entry, results);
    }
    grid[r][c] = 0;
}

/// Dual immaculate quasisymmetric function S\*_α in the fundamental basis.
///
/// S\*_α = Σ_{T ∈ SIT(α)} F_{Des(T)}
///
/// where SIT(α) is the set of standard immaculate tableaux of shape α.
///
/// Reference: Berg–Bergeron–Saliola–Serrano–Zabrocki,
/// *A lift of the Schur and Hall-Littlewood bases to non-commutative
/// symmetric functions*, Canad. J. Math. 66 (2014), 525–565.
pub fn dual_immaculate<C: Ring>(alpha: &[u32]) -> QSymFunction<C> {
    let tableaux = ImmaculateTableau::enumerate(alpha);
    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();
    for t in &tableaux {
        let comp = t.descent_composition();
        let entry = terms.entry(comp).or_insert_with(C::zero);
        *entry = entry.clone() + C::one();
    }
    QSymFunction::from_terms(QSymBasis::Fundamental, terms)
}

/// Dual immaculate quasisymmetric function S\*_α in the monomial basis.
///
/// Uses semistandard immaculate tableaux: rows weakly increasing and first
/// column strictly increasing. The `max_entry` parameter bounds the packed
/// alphabet length; use `|alpha|` for the full basis element.
pub fn dual_immaculate_monomial<C: Ring>(alpha: &[u32], max_entry: u32) -> QSymFunction<C> {
    let tableaux = ImmaculateTableau::enumerate_semistandard(alpha, max_entry);
    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();
    for t in &tableaux {
        if let Some(comp) = packed_content_composition(&t.rows, max_entry) {
            let entry = terms.entry(comp).or_insert_with(C::zero);
            *entry = entry.clone() + C::one();
        }
    }
    QSymFunction::from_terms(QSymBasis::Monomial, terms)
}

/// Row-strict dual immaculate quasisymmetric function in the monomial basis.
///
/// Uses row-strict semistandard immaculate tableaux: rows strictly increasing,
/// first column weakly increasing. The function is:
///
/// ℜS\*_α = Σ_{T ∈ SSIT(α,n)} x^{weight(T)}
///
/// expressed in the monomial QSym basis by packed content. The `max_entry`
/// parameter bounds the packed alphabet length; use `|alpha|` for the full
/// basis element.
pub fn row_strict_dual_immaculate<C: Ring>(alpha: &[u32], max_entry: u32) -> QSymFunction<C> {
    let tableaux = ImmaculateTableau::enumerate_row_strict_semistandard(alpha, max_entry);
    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();
    for t in &tableaux {
        if let Some(comp) = packed_content_composition(&t.rows, max_entry) {
            let entry = terms.entry(comp).or_insert_with(C::zero);
            *entry = entry.clone() + C::one();
        }
    }
    QSymFunction::from_terms(QSymBasis::Monomial, terms)
}

/// Row-strict dual immaculate quasisymmetric function in the fundamental basis.
///
/// This uses the identity `RS*_alpha = psi(S*_alpha)`, avoiding any finite
/// alphabet bound.
pub fn row_strict_dual_immaculate_fundamental<C: Ring>(alpha: &[u32]) -> QSymFunction<C> {
    dual_immaculate::<C>(alpha).psi_involution()
}

// ── Fundamental slide polynomials ───────────────────────────────────

/// Fundamental slide polynomial 𝔉_α as a QSym function in the fundamental basis.
///
/// 𝔉_α = Σ F_β where the sum is over compositions β that are
/// "slide-compatible" with α: β refines some weak composition γ
/// satisfying the partial sum condition ∑_{i≤k} α_i ≤ ∑_{i≤k} γ_i for all k.
pub fn fundamental_slide<C: Ring>(alpha: &[u32]) -> QSymFunction<C> {
    let n = alpha.len();
    let total: u32 = alpha.iter().sum();
    if total == 0 || n == 0 {
        return QSymFunction::zero(QSymBasis::Fundamental);
    }
    let alpha_cumsum: Vec<u32> = alpha
        .iter()
        .scan(0u32, |acc, &x| {
            *acc += x;
            Some(*acc)
        })
        .collect();

    let flat_alpha: Vec<u32> = alpha.iter().copied().filter(|&x| x > 0).collect();

    // Generate all weak compositions of `total` into n parts
    // that dominate alpha (partial sums ≥ alpha's partial sums)
    // and whose non-zero parts refine flat_alpha.
    let mut terms: BTreeMap<Composition, C> = BTreeMap::new();

    // Enumerate compatible weak compositions
    let mut gamma = vec![0u32; n];
    enumerate_slide_compatible(&mut gamma, 0, total, &alpha_cumsum, &flat_alpha, &mut terms);

    QSymFunction::from_terms(QSymBasis::Fundamental, terms)
}

fn enumerate_slide_compatible<C: Ring>(
    gamma: &mut Vec<u32>,
    pos: usize,
    remaining: u32,
    alpha_cumsum: &[u32],
    flat_alpha: &[u32],
    terms: &mut BTreeMap<Composition, C>,
) {
    let n = gamma.len();
    if pos == n {
        if remaining == 0 {
            let flat_gamma: Vec<u32> = gamma.iter().copied().filter(|&x| x > 0).collect();
            if is_composition_refinement(&flat_gamma, flat_alpha) {
                let comp = Composition::new(flat_gamma);
                let entry = terms.entry(comp).or_insert_with(C::zero);
                *entry = entry.clone() + C::one();
            }
        }
        return;
    }

    let cumsum_so_far: u32 = gamma[..pos].iter().sum();
    let min_needed = alpha_cumsum[pos].saturating_sub(cumsum_so_far);
    let positions_left = n - pos - 1;

    for v in min_needed..=remaining {
        let new_remaining = remaining - v;
        if new_remaining > 0 && positions_left == 0 {
            continue;
        }
        gamma[pos] = v;
        // Check partial sum condition for this position
        let new_cumsum = cumsum_so_far + v;
        if new_cumsum >= alpha_cumsum[pos] {
            enumerate_slide_compatible(
                gamma,
                pos + 1,
                new_remaining,
                alpha_cumsum,
                flat_alpha,
                terms,
            );
        }
        gamma[pos] = 0;
    }
}

/// Check if `fine` is a refinement of `coarse` (as compositions).
/// I.e., one can merge consecutive parts of `fine` to recover `coarse`.
fn is_composition_refinement(fine: &[u32], coarse: &[u32]) -> bool {
    let mut ci = 0;
    let mut partial = 0u32;
    for &part in fine {
        partial += part;
        if ci < coarse.len() && partial == coarse[ci] {
            ci += 1;
            partial = 0;
        } else if ci < coarse.len() && partial > coarse[ci] {
            return false;
        }
    }
    ci == coarse.len() && partial == 0
}

// ── Utility ─────────────────────────────────────────────────────────

/// Convert a descent set (1-indexed positions) to a composition of n.
pub fn descent_set_to_composition(des: &[u32], n: u32) -> Composition {
    if n == 0 {
        return Composition::new(Vec::new());
    }
    let mut parts = Vec::new();
    let mut prev = 0u32;
    for &d in des {
        parts.push(d - prev);
        prev = d;
    }
    parts.push(n - prev);
    Composition::new(parts)
}

/// Convert a composition to its descent set (1-indexed).
pub fn composition_to_descent_set(comp: &Composition) -> Vec<u32> {
    let mut des = Vec::new();
    let mut pos = 0u32;
    let parts = comp.parts();
    for &p in parts.iter().take(parts.len().saturating_sub(1)) {
        pos += p;
        des.push(pos);
    }
    des
}

fn packed_content_composition(rows: &[Vec<u32>], max_entry: u32) -> Option<Composition> {
    let mut content = vec![0u32; max_entry as usize];
    let mut max_used = 0u32;

    for row in rows {
        for &value in row {
            if value == 0 || value > max_entry {
                return None;
            }
            content[value as usize - 1] += 1;
            max_used = max_used.max(value);
        }
    }

    if max_used == 0 {
        return Some(Composition::empty());
    }

    let packed = &content[..max_used as usize];
    if packed.iter().any(|&count| count == 0) {
        return None;
    }

    Some(Composition::new(packed.to_vec()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_composition_tableaux_enumeration() {
        // Shape (2, 1) with max entry 3
        let tabs = CompositionTableau::enumerate(&[2, 1], 3);
        assert!(!tabs.is_empty());
        // Each should have correct shape
        for t in &tabs {
            assert_eq!(t.rows.len(), 2);
            assert_eq!(t.rows[0].len(), 2);
            assert_eq!(t.rows[1].len(), 1);
            // Rows weakly increasing
            assert!(t.rows[0][0] <= t.rows[0][1]);
            // First column strictly increasing
            assert!(t.rows[0][0] < t.rows[1][0]);
        }
    }

    #[test]
    fn test_qsym_schur_in_fundamental() {
        // S_(1) in 2 variables = F_(1) = x1 + x2
        let s: QSymFunction<i64> = qsym_schur(&[1], 2);
        assert_eq!(s.basis(), QSymBasis::Fundamental);
        assert!(!s.is_zero());
    }

    #[test]
    fn test_reverse_composition_tableaux_source_example() {
        // Tewari--van Willigenburg, Example 2.7: SRCT(2,1,3) has
        // descent compositions (2,1,3), (2,2,2), and (1,2,1,2).
        let tabs = ReverseCompositionTableau::enumerate_standard(&[2, 1, 3]);
        let descents: Vec<Composition> = tabs.iter().map(|t| t.descent_composition()).collect();

        assert_eq!(tabs.len(), 3);
        assert!(descents.contains(&Composition::new(vec![2, 1, 3])));
        assert!(descents.contains(&Composition::new(vec![2, 2, 2])));
        assert!(descents.contains(&Composition::new(vec![1, 2, 1, 2])));
    }

    #[test]
    fn test_qsym_schur_source_example() {
        // qSchur_(2,1,3) = F_(2,1,3) + F_(2,2,2) + F_(1,2,1,2).
        let s: QSymFunction<i64> = qsym_schur(&[2, 1, 3], 6);

        assert_eq!(s.coefficient(&Composition::new(vec![2, 1, 3])), 1);
        assert_eq!(s.coefficient(&Composition::new(vec![2, 2, 2])), 1);
        assert_eq!(s.coefficient(&Composition::new(vec![1, 2, 1, 2])), 1);
        assert_eq!(s.terms().len(), 3);
    }

    #[test]
    fn test_young_qsym_schur_source_relation() {
        // Luoto--Mykytiuk--van Willigenburg define the Young quasisymmetric
        // Schur functions by reversing the index and applying rho.
        let alpha = [3, 1, 2];
        let from_relation = qsym_schur::<i64>(&[2, 1, 3], 6).rho_involution();
        let yqs = young_qsym_schur::<i64>(&alpha);

        assert_eq!(yqs, from_relation);
        assert_eq!(yqs.coefficient(&Composition::new(vec![3, 1, 2])), 1);
        assert_eq!(yqs.coefficient(&Composition::new(vec![2, 2, 2])), 1);
        assert_eq!(yqs.coefficient(&Composition::new(vec![2, 1, 2, 1])), 1);
        assert_eq!(yqs.terms().len(), 3);
    }

    #[test]
    fn test_row_strict_qsym_schur_omega_relation() {
        // Mason--Remmel give qSchur_alpha = omega(rsqSchur_alpha).
        // Since omega is an involution, rsqSchur_alpha = omega(qSchur_alpha).
        let alpha = [2, 1, 3];
        let from_relation = qsym_schur::<i64>(&alpha, 6).omega_involution();
        let rsqs = row_strict_qsym_schur::<i64>(&alpha);

        assert_eq!(rsqs, from_relation);
        assert_eq!(rsqs.coefficient(&Composition::new(vec![1, 1, 3, 1])), 1);
        assert_eq!(rsqs.coefficient(&Composition::new(vec![1, 2, 2, 1])), 1);
        assert_eq!(rsqs.coefficient(&Composition::new(vec![1, 3, 2])), 1);
        assert_eq!(rsqs.terms().len(), 3);
    }

    #[test]
    fn test_qsym_schur_partition_is_schur() {
        // For partition shape, QSym Schur should be the Schur function
        // (i.e., s-positive in fundamental basis).
        let s: QSymFunction<i64> = qsym_schur(&[2, 1], 3);
        // All coefficients should be non-negative
        for (_, &c) in s.terms() {
            assert!(c >= 0, "QSym Schur of partition shape should be F-positive");
        }
    }

    #[test]
    fn test_fundamental_slide() {
        // 𝔉_(1) = F_(1) (single part)
        let fs: QSymFunction<i64> = fundamental_slide(&[1]);
        assert!(!fs.is_zero());
    }

    #[test]
    fn test_descent_set_composition_roundtrip() {
        let comp = Composition::new(vec![2, 1, 3]);
        let des = composition_to_descent_set(&comp);
        let back = descent_set_to_composition(&des, comp.size());
        assert_eq!(back, comp);
    }

    #[test]
    fn test_immaculate_tableaux_enumeration() {
        // Shape (2, 1): entries 1,2,3
        // Rows strictly increasing, first column strictly increasing
        // Row 0: [a, b] with a < b. Row 1: [c] with c > a (first col).
        // Valid: [1,2],[3]; [1,3],[2]; [1,3],[3]? No, 3 used twice.
        // Actually entries are 1..3 each once: [1,2],[3] and [1,3],[2]
        // and [2,3],[x] where x>2, x=3 but 3 used → no more.
        // So: 2 standard immaculate tableaux of shape (2,1)
        let tabs = ImmaculateTableau::enumerate(&[2, 1]);
        assert_eq!(tabs.len(), 2);
        for t in &tabs {
            // Rows strictly increasing
            for row in &t.rows {
                for w in row.windows(2) {
                    assert!(w[0] < w[1]);
                }
            }
            // First column strictly increasing
            if t.rows.len() > 1 && !t.rows[0].is_empty() && !t.rows[1].is_empty() {
                assert!(t.rows[0][0] < t.rows[1][0]);
            }
        }
    }

    #[test]
    fn test_immaculate_vs_syt_count() {
        // For partition shapes, immaculate tableaux ≥ SYT
        // (immaculate has fewer constraints than SYT)
        // Shape (3, 2, 1): 16 SYT, but more immaculate tableaux
        let sit = ImmaculateTableau::enumerate(&[3, 2, 1]);
        let syt_count = 16; // known: f^{321} = 16
        assert!(
            sit.len() >= syt_count,
            "SIT({}) should be >= SYT({}): {} vs {}",
            "3,2,1",
            "3,2,1",
            sit.len(),
            syt_count
        );
    }

    #[test]
    fn test_dual_immaculate_f_positive() {
        // Dual immaculate functions are F-positive by definition
        let di: QSymFunction<i64> = dual_immaculate(&[2, 1]);
        assert!(!di.is_zero());
        for (_, &c) in di.terms() {
            assert!(c > 0, "dual immaculate should be F-positive");
        }
    }

    #[test]
    fn test_dual_immaculate_partition_shape() {
        // For shape (n), dual immaculate S*_(n) = F_(n) (single row)
        let di: QSymFunction<i64> = dual_immaculate(&[3]);
        assert_eq!(di.terms().len(), 1);
        let comp = Composition::new(vec![3]);
        assert_eq!(di.coefficient(&comp), 1);
    }

    #[test]
    fn test_dual_immaculate_single_column() {
        // Shape (1,1,1): n=3, single column = strictly increasing first column
        // Only one SIT: [1],[2],[3]. Des = {1,2}. Comp = (1,1,1).
        let di: QSymFunction<i64> = dual_immaculate(&[1, 1, 1]);
        let comp = Composition::new(vec![1, 1, 1]);
        assert_eq!(di.coefficient(&comp), 1);
        assert_eq!(di.terms().len(), 1);
    }

    #[test]
    fn test_dual_immaculate_monomial_small_example() {
        // Since S*_(2,1) = F_(1,2) + F_(2,1), its monomial expansion is
        // M_(1,2) + M_(2,1) + 2 M_(1,1,1).
        let di: QSymFunction<i64> = dual_immaculate_monomial(&[2, 1], 3);
        assert_eq!(di.basis(), QSymBasis::Monomial);
        assert_eq!(di.coefficient(&Composition::new(vec![1, 2])), 1);
        assert_eq!(di.coefficient(&Composition::new(vec![2, 1])), 1);
        assert_eq!(di.coefficient(&Composition::new(vec![1, 1, 1])), 2);
        assert_eq!(di.terms().len(), 3);

        let larger_alphabet: QSymFunction<i64> = dual_immaculate_monomial(&[2, 1], 4);
        assert_eq!(larger_alphabet, di);
    }

    #[test]
    fn test_semistandard_immaculate_count() {
        // SS immaculate tableaux of shape (1,1) with max_entry 3:
        // Rows weakly increasing (trivial: 1 entry each)
        // First col strictly increasing: entry[0] < entry[1]
        // Entries in {1,2,3}: (1,2), (1,3), (2,3) → 3 tableaux
        let tabs = ImmaculateTableau::enumerate_semistandard(&[1, 1], 3);
        assert_eq!(tabs.len(), 3);
    }

    #[test]
    fn test_row_strict_semistandard_immaculate_count() {
        // Row-strict SS immaculate tableaux of shape (1,1) with max_entry 3:
        // Rows strict are trivial, and first column is weakly increasing.
        let tabs = ImmaculateTableau::enumerate_row_strict_semistandard(&[1, 1], 3);
        assert_eq!(tabs.len(), 6);
    }

    #[test]
    fn test_row_strict_dual_immaculate_paper_example_310() {
        // Niese--Sundaram--van Willigenburg--Vega--Wang, Example 3.10:
        // RS*_(1,2) = M_(2,1) + M_(1,1,1).
        let rsdi: QSymFunction<i64> = row_strict_dual_immaculate(&[1, 2], 3);
        assert_eq!(rsdi.basis(), QSymBasis::Monomial);
        assert_eq!(rsdi.coefficient(&Composition::new(vec![2, 1])), 1);
        assert_eq!(rsdi.coefficient(&Composition::new(vec![1, 1, 1])), 1);
        assert_eq!(rsdi.coefficient(&Composition::new(vec![1, 2])), 0);
        assert_eq!(rsdi.terms().len(), 2);

        let larger_alphabet: QSymFunction<i64> = row_strict_dual_immaculate(&[1, 2], 4);
        assert_eq!(larger_alphabet, rsdi);
    }

    #[test]
    fn test_row_strict_dual_immaculate_matches_psi() {
        let from_tableaux: QSymFunction<i64> =
            row_strict_dual_immaculate(&[1, 2], 3).to_fundamental_basis();
        let from_psi = row_strict_dual_immaculate_fundamental::<i64>(&[1, 2]);

        assert_eq!(from_tableaux, from_psi);
        assert_eq!(from_psi.coefficient(&Composition::new(vec![2, 1])), 1);
        assert_eq!(from_psi.terms().len(), 1);
    }

    #[test]
    fn test_row_strict_dual_immaculate_matches_psi_small_degrees() {
        for degree in 1..=4 {
            for alpha in Composition::integer_compositions(degree) {
                let from_tableaux: QSymFunction<i64> =
                    row_strict_dual_immaculate(alpha.parts(), degree).to_fundamental_basis();
                let from_psi = row_strict_dual_immaculate_fundamental::<i64>(alpha.parts());
                assert_eq!(from_tableaux, from_psi, "failed for alpha = {}", alpha);
            }
        }
    }

    #[test]
    fn test_dual_immaculate_monomial_matches_fundamental_small_degrees() {
        for degree in 1..=4 {
            for alpha in Composition::integer_compositions(degree) {
                let from_tableaux: QSymFunction<i64> =
                    dual_immaculate_monomial(alpha.parts(), degree).to_fundamental_basis();
                let from_standard_tableaux = dual_immaculate::<i64>(alpha.parts());
                assert_eq!(
                    from_tableaux, from_standard_tableaux,
                    "failed for alpha = {}",
                    alpha
                );
            }
        }
    }

    #[test]
    fn test_composition_refinement() {
        assert!(is_composition_refinement(&[1, 1, 2], &[2, 2]));
        assert!(is_composition_refinement(&[1, 1, 1, 1], &[2, 2]));
        assert!(!is_composition_refinement(&[1, 2, 1], &[2, 2]));
        assert!(is_composition_refinement(&[3], &[3]));
    }
}
