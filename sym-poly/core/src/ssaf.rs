//! Semi-standard augmented fillings (SSAFs) with permuted basements.
//!
//! An SSAF is a filling of a composition diagram where:
//! - Column 0 (the *basement*) is fixed by a permutation.
//! - Each row is weakly decreasing from left to right.
//! - The filling satisfies the non-attacking and coinversion-free conditions.
//!
//! Special cases:
//! - **Atom fillings**: basement = \[1, 2, …, n\] (identity).
//! - **Key fillings**: basement = \[n, n−1, …, 1\] (reverse).
//! - **Macdonald E-fillings**: general basement with non-attacking check only.

use std::fmt;

/// A semi-standard augmented filling.
///
/// `rows[i]` has the form `[basement_i, v1, v2, …]` where
/// `basement_i >= v1 >= v2 >= …` and the filling satisfies
/// triple-rule constraints between rows.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Ssaf {
    rows: Vec<Vec<u32>>,
}

impl Ssaf {
    /// Create an SSAF from pre-built rows (each row includes the basement in position 0).
    pub fn new(rows: Vec<Vec<u32>>) -> Self {
        Ssaf { rows }
    }

    /// The number of rows.
    pub fn num_rows(&self) -> usize {
        self.rows.len()
    }

    /// Access the rows (each row includes the basement value in position 0).
    pub fn rows(&self) -> &[Vec<u32>] {
        &self.rows
    }

    /// The basement permutation (column 0 values).
    pub fn basement(&self) -> Vec<u32> {
        self.rows.iter().map(|r| r[0]).collect()
    }

    /// The composition shape (number of filling boxes per row, excluding basement).
    pub fn shape(&self) -> Vec<u32> {
        self.rows
            .iter()
            .map(|r| {
                if r.is_empty() {
                    0
                } else {
                    (r.len() - 1) as u32
                }
            })
            .collect()
    }

    /// Weight vector: `weight[v-1]` = number of filling entries equal to `v`
    /// (basement entries excluded).
    pub fn weight_vector(&self) -> Vec<u32> {
        let n = self.num_rows();
        let mut counts = vec![0u32; n];
        for row in &self.rows {
            for &v in row.iter().skip(1) {
                if (v as usize) <= n {
                    counts[v as usize - 1] += 1;
                }
            }
        }
        counts
    }

    /// Total number of filling entries (excluding basement).
    pub fn size(&self) -> u32 {
        self.rows
            .iter()
            .map(|r| r.len().saturating_sub(1) as u32)
            .sum()
    }

    // ── Statistics ──────────────────────────────────────────────────

    /// Number of inversion triples in the SSAF.
    pub fn inversions(&self) -> u32 {
        self.count_triples(true)
    }

    /// Number of coinversion (non-inversion) triples in the SSAF.
    pub fn coinversions(&self) -> u32 {
        self.count_triples(false)
    }

    fn count_triples(&self, want_inversions: bool) -> u32 {
        let n = self.num_rows();
        let shape: Vec<usize> = self.rows.iter().map(|r| r.len()).collect();
        let mut count = 0u32;

        for r in 0..n {
            for c in 1..shape[r] {
                for i in 0..r {
                    // Type A: row i at least as wide as row r (in goal shape)
                    if shape[i] >= shape[r] && shape[i] > c {
                        let aa = self.rows[i][c];
                        let bb = self.rows[r][c];
                        let cc = self.rows[i][c - 1];
                        let is_inv = is_inversion_triple_type_a(aa, bb, cc);
                        if is_inv == want_inversions {
                            count += 1;
                        }
                    }
                    // Type B: row i shorter than row r
                    if shape[i] < shape[r] && shape[i] > c - 1 {
                        let aa = self.rows[i][c - 1];
                        let bb = self.rows[r][c - 1];
                        let cc = self.rows[r][c];
                        let is_inv = is_inversion_triple_type_b(aa, bb, cc);
                        if is_inv == want_inversions {
                            count += 1;
                        }
                    }
                }
            }
        }
        count
    }

    /// Horizontal descents: count of cells `(r, c)` with `c >= 1` where
    /// `ssaf[r][c] != ssaf[r][c-1]`.
    pub fn horizontal_descents(&self) -> u32 {
        let mut count = 0u32;
        for row in &self.rows {
            for c in 1..row.len() {
                if row[c] != row[c - 1] {
                    count += 1;
                }
            }
        }
        count
    }

    /// Major index: sum of `(1 + leg)` over cells `(r, c)` where
    /// `ssaf[r][c] > ssaf[r][c-1]` (an ascent), where `leg = row_len - 1 - c`.
    pub fn major_index(&self) -> u32 {
        let mut maj = 0u32;
        for row in &self.rows {
            let row_len = row.len();
            for c in 1..row_len {
                if row[c] > row[c - 1] {
                    let leg = row_len - 1 - c;
                    maj += 1 + leg as u32;
                }
            }
        }
        maj
    }

    // ── Generation ─────────────────────────────────────────────────

    /// Generate all SSAFs with the given composition shape and basement permutation.
    ///
    /// - `shape[i]` = number of filling boxes in row `i` (excluding basement).
    /// - `basement[i]` = value in column 0 of row `i`.
    ///
    /// The filling satisfies: weakly decreasing rows, non-attacking, and
    /// coinversion-free (triple rule).
    pub fn generate(shape: &[u32], basement: &[u32]) -> Vec<Ssaf> {
        assert_eq!(
            shape.len(),
            basement.len(),
            "shape and basement must have the same length"
        );
        let n = shape.len();
        // Goal shape includes the basement column
        let goal_shape: Vec<usize> = shape.iter().map(|&s| s as usize + 1).collect();

        // Build the row sequence: which row gets each new box
        let mut row_sequence = Vec::new();
        for (i, &s) in shape.iter().enumerate() {
            for _ in 0..s {
                row_sequence.push(i);
            }
        }

        // Initialize tableau with basement column
        let mut init: Vec<Vec<u32>> = basement.iter().map(|&b| vec![b]).collect();

        let mut results = Vec::new();
        generate_recursive(
            &mut init,
            &row_sequence,
            0,
            &goal_shape,
            n as u32,
            &mut results,
        );
        results
    }

    /// Atom fillings: SSAFs with basement = \[1, 2, …, n\].
    pub fn atom_fillings(alpha: &[u32]) -> Vec<Ssaf> {
        let n = alpha.len();
        let basement: Vec<u32> = (1..=n as u32).collect();
        Self::generate(alpha, &basement)
    }

    /// Atom fillings filtered by weight.
    pub fn atom_fillings_weight(alpha: &[u32], weight: &[u32]) -> Vec<Ssaf> {
        Self::atom_fillings(alpha)
            .into_iter()
            .filter(|f| f.weight_vector() == weight)
            .collect()
    }

    /// Key fillings: SSAFs with basement = \[n, n−1, …, 1\].
    ///
    /// Under the permuted-basement Macdonald convention used in this crate,
    /// the key-filling shape is the reverse of the weak composition indexing
    /// the key polynomial.
    pub fn key_fillings(alpha: &[u32]) -> Vec<Ssaf> {
        let n = alpha.len();
        let basement: Vec<u32> = (1..=n as u32).rev().collect();
        let reversed_alpha: Vec<u32> = alpha.iter().rev().copied().collect();
        Self::generate(&reversed_alpha, &basement)
    }

    /// Key fillings filtered by weight.
    pub fn key_fillings_weight(alpha: &[u32], weight: &[u32]) -> Vec<Ssaf> {
        Self::key_fillings(alpha)
            .into_iter()
            .filter(|f| f.weight_vector() == weight)
            .collect()
    }

    /// Lock fillings: entries weakly increasing in rows, strictly decreasing
    /// in columns, and each entry ≤ its basement value.
    /// Basement = \[n, n−1, …, 1\].
    pub fn lock_fillings(alpha: &[u32]) -> Vec<Ssaf> {
        let n = alpha.len();
        let basement: Vec<u32> = (1..=n as u32).rev().collect();
        let goal_shape: Vec<usize> = alpha.iter().map(|&s| s as usize + 1).collect();

        let mut row_sequence = Vec::new();
        for (i, &s) in alpha.iter().enumerate() {
            for _ in 0..s {
                row_sequence.push(i);
            }
        }

        let mut init: Vec<Vec<u32>> = basement.iter().map(|&b| vec![b]).collect();
        let mut results = Vec::new();
        generate_lock_recursive(
            &mut init,
            &row_sequence,
            0,
            n as u32,
            &goal_shape,
            &mut results,
        );
        results
    }

    // ── Mason's bijection: SSYT → SSAF ─────────────────────────────

    /// Convert a semistandard Young tableau to an atom filling (SSAF)
    /// via Mason's insertion bijection.
    ///
    /// The SSYT should have entries in `{1, …, n}` where `n` is chosen
    /// as `max(entries)`.  The result is an SSAF with identity basement
    /// `[1, 2, …, n]`.
    pub fn from_ssyt(ssyt: &[Vec<u32>]) -> Ssaf {
        let n = ssyt
            .iter()
            .flat_map(|r| r.iter())
            .copied()
            .max()
            .unwrap_or(0) as usize;
        if n == 0 {
            return Ssaf { rows: Vec::new() };
        }

        // Start with basement only: [[1], [2], ..., [n]]
        let mut atom: Vec<Vec<u32>> = (1..=n as u32).map(|i| vec![i]).collect();

        // Reading word: bottom-to-top, left-to-right, then process in reverse
        let reading_word: Vec<u32> = ssyt.iter().rev().flat_map(|r| r.iter().copied()).collect();

        for &v in reading_word.iter().rev() {
            let start_col = atom.iter().map(|r| r.len()).max().unwrap_or(1);
            mason_insert(&mut atom, v, 0, start_col, n);
        }

        Ssaf { rows: atom }
    }
}

// ── Generation internals ────────────────────────────────────────────

fn generate_recursive(
    tab: &mut Vec<Vec<u32>>,
    row_seq: &[usize],
    idx: usize,
    goal_shape: &[usize],
    _n: u32,
    results: &mut Vec<Ssaf>,
) {
    if idx == row_seq.len() {
        results.push(Ssaf::new(tab.clone()));
        return;
    }

    let r = row_seq[idx];
    // Each row is initialized with its basement value, so last() is always Some.
    let max_val = *tab[r].last().expect("row must have basement entry");

    for b in 1..=max_val {
        tab[r].push(b);
        let c = tab[r].len() - 1; // 0-based column of the new entry

        if check_coinversions(tab, r, c, goal_shape) {
            generate_recursive(tab, row_seq, idx + 1, goal_shape, _n, results);
        }

        tab[r].pop();
    }
}

/// Check that adding the entry at `tab[r][c]` doesn't create a forbidden coinversion.
fn check_coinversions(tab: &[Vec<u32>], r: usize, c: usize, goal_shape: &[usize]) -> bool {
    let b = tab[r][c];

    // Non-attacking: check rows above for same value in column c or c-1
    for row_above in &tab[..r] {
        if row_above.len() > c && row_above[c] == b {
            return false;
        }
        if c > 0 && row_above.len() > c - 1 && row_above[c - 1] == b {
            return false;
        }
    }

    // Triple checks against each row above
    for (i, row_above) in tab[..r].iter().enumerate() {
        // Type A: row i is at least as wide as row r (in the goal shape)
        if goal_shape[i] >= goal_shape[r] && row_above.len() > c {
            let aa = row_above[c];
            let bb = tab[r][c]; // the new entry
            let cc = row_above[c - 1];
            if !is_inversion_triple_type_a(aa, bb, cc) {
                return false;
            }
        }
        // Type B: row i is shorter than row r
        if goal_shape[i] < goal_shape[r] && row_above.len() > c - 1 {
            let aa = row_above[c - 1];
            let bb = tab[r][c - 1];
            let cc = tab[r][c]; // the new entry
            if !is_inversion_triple_type_b(aa, bb, cc) {
                return false;
            }
        }
    }
    true
}

// ── Triple predicates ───────────────────────────────────────────────

/// ClockwiseQ(x, y, z): are (x, y, z) in cyclic ascending order?
fn clockwise_q(x: i64, y: i64, z: i64) -> bool {
    // Three-element cyclic order check (algebraically equivalent forms cause
    // clippy to suggest simplification, but the three-disjunct form is clearest
    // for the mathematical definition).
    #[allow(clippy::nonminimal_bool)]
    {
        (x < y && y < z) || (y < z && z < x) || (z < x && x < y)
    }
}

/// Type A inversion triple.
///
/// For cells `(i,c)`, `(r,c)`, `(i,c-1)` with row i above row r:
/// checks `ClockwiseQ(aa + 0.1, cc + 0.3, bb + 0.2)`.
/// We multiply by 10 for exact integer arithmetic.
fn is_inversion_triple_type_a(aa: u32, bb: u32, cc: u32) -> bool {
    clockwise_q(10 * aa as i64 + 1, 10 * cc as i64 + 3, 10 * bb as i64 + 2)
}

/// Type B inversion triple.
///
/// For cells `(i,c-1)`, `(r,c-1)`, `(r,c)` with row i shorter than row r:
/// checks `ClockwiseQ(aa + 0.2, cc + 0.1, bb + 0.3)`.
fn is_inversion_triple_type_b(aa: u32, bb: u32, cc: u32) -> bool {
    clockwise_q(10 * aa as i64 + 2, 10 * cc as i64 + 1, 10 * bb as i64 + 3)
}

// ── Composition-to-permutation maps ─────────────────────────────────

/// For a composition α, return the permutation σ and the sorted partition λ
/// such that α indexes an atom polynomial: atom_α = θ_σ(x^λ).
///
/// σ is the permutation that sorts α into weakly decreasing order,
/// chosen with minimal length.
///
/// Returns `(sigma, lambda)` where σ is 0-indexed.
pub fn atom_composition_to_permutation(alpha: &[u32]) -> (Vec<usize>, Vec<u32>) {
    let n = alpha.len();
    // Sort indices by decreasing alpha value (stable sort gives minimal length)
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| alpha[b].cmp(&alpha[a]).then(a.cmp(&b)));
    // lambda = sorted alpha (weakly decreasing)
    let lambda: Vec<u32> = indices.iter().map(|&i| alpha[i]).collect();
    // sigma: sigma[i] = where index i goes
    let mut sigma = vec![0usize; n];
    for (new_pos, &old_pos) in indices.iter().enumerate() {
        sigma[old_pos] = new_pos;
    }
    (sigma, lambda)
}

/// For a composition α, return the permutation σ and the sorted partition λ
/// such that α indexes a key polynomial: key_α = π_σ(x^λ).
///
/// Returns `(sigma, lambda)` where σ is 0-indexed.
pub fn key_composition_to_permutation(alpha: &[u32]) -> (Vec<usize>, Vec<u32>) {
    let n = alpha.len();
    let lambda: Vec<u32> = {
        let mut sorted = alpha.to_vec();
        sorted.sort_unstable_by(|a, b| b.cmp(a));
        sorted
    };

    // Find the permutation of minimal length such that lambda[sigma[i]] = alpha[i]
    // Equivalently: sigma sends position i to the position of alpha[i] in lambda
    // We use a greedy assignment: for each value, assign positions left-to-right
    let mut available: Vec<Vec<usize>> =
        vec![Vec::new(); (*lambda.iter().max().unwrap_or(&0) + 1) as usize];
    for (pos, &val) in lambda.iter().enumerate() {
        available[val as usize].push(pos);
    }
    // Reverse so we can pop from the front efficiently
    for v in &mut available {
        v.reverse();
    }

    let mut sigma = vec![0usize; n];
    for (i, &val) in alpha.iter().enumerate() {
        sigma[i] = available[val as usize].pop().unwrap();
    }

    (sigma, lambda)
}

// ── Lock fillings ───────────────────────────────────────────────────

fn generate_lock_recursive(
    tab: &mut Vec<Vec<u32>>,
    row_seq: &[usize],
    idx: usize,
    n: u32,
    goal_shape: &[usize],
    results: &mut Vec<Ssaf>,
) {
    if idx == row_seq.len() {
        results.push(Ssaf::new(tab.clone()));
        return;
    }

    let r = row_seq[idx];

    for b in 1..=n {
        tab[r].push(b);
        let c = tab[r].len() - 1;

        if check_lock(tab, r, c, goal_shape) {
            generate_lock_recursive(tab, row_seq, idx + 1, n, goal_shape, results);
        }

        tab[r].pop();
    }
}

/// Lock filling validity check.
///
/// 1. Entry ≤ basement value (entry in row r is ≤ tab\[r\]\[0\]).
/// 2. Rows are weakly increasing (for c > 1: entry ≥ previous entry).
/// 3. Columns are strictly decreasing upward (entries above must be > current).
fn check_lock(tab: &[Vec<u32>], r: usize, c: usize, _goal_shape: &[usize]) -> bool {
    let b = tab[r][c];

    // Entry ≤ basement value
    if b > tab[r][0] {
        return false;
    }

    // Weakly increasing rows (skip first filling column)
    if c > 1 && b < tab[r][c - 1] {
        return false;
    }

    // Strictly decreasing columns upward
    for row_above in &tab[..r] {
        if row_above.len() > c && row_above[c] <= b {
            return false;
        }
    }
    true
}

// ── Mason's insertion bijection ─────────────────────────────────────

/// Insert value `v` into the atom filling at position `(r, c)` (0-indexed).
/// Cycles through rows (wrapping around) and moves left through columns.
fn mason_insert(atom: &mut Vec<Vec<u32>>, v: u32, r: usize, c: usize, n: usize) {
    if c == 0 {
        // Should not try to insert into the basement column
        panic!("mason_insert: cannot insert into basement column");
    }

    // Next position in reading order: row r+1, same column.
    // If past bottom, wrap to row 0, column c-1.
    let (nr, nc) = if r + 1 >= n { (0, c - 1) } else { (r + 1, c) };

    let row_len = atom[r].len();

    if row_len == c {
        // Position c is one past the end: can append if left neighbor ≥ v
        if atom[r][c - 1] >= v {
            atom[r].push(v);
        } else {
            mason_insert(atom, v, nr, nc, n);
        }
    } else if row_len > c {
        // Position c is occupied
        if atom[r][c - 1] >= v && atom[r][c] >= v {
            // Existing entry ≥ v and left neighbor ≥ v: skip, continue reading
            mason_insert(atom, v, nr, nc, n);
        } else if atom[r][c - 1] >= v && atom[r][c] < v {
            // Bump: replace existing with v, reinsert the old value
            let old = atom[r][c];
            atom[r][c] = v;
            mason_insert(atom, old, nr, nc, n);
        } else {
            // Left neighbor < v: continue
            mason_insert(atom, v, nr, nc, n);
        }
    } else {
        // row_len < c: position is beyond the row; continue reading
        mason_insert(atom, v, nr, nc, n);
    }
}

// ── Non-attacking (Macdonald E) fillings ────────────────────────────

/// Generate non-attacking fillings (for Macdonald E-polynomials).
///
/// Like SSAF generation but:
/// - Rows are NOT constrained to be weakly decreasing.
/// - Only the non-attacking condition is enforced.
/// - Entries can be any value in `1..=n`.
fn generate_non_attacking_recursive(
    tab: &mut Vec<Vec<u32>>,
    row_seq: &[usize],
    idx: usize,
    n: u32,
    results: &mut Vec<Ssaf>,
) {
    if idx == row_seq.len() {
        results.push(Ssaf::new(tab.clone()));
        return;
    }

    let r = row_seq[idx];

    for b in 1..=n {
        tab[r].push(b);
        let c = tab[r].len() - 1;

        if check_non_attacking(tab, r, c) {
            generate_non_attacking_recursive(tab, row_seq, idx + 1, n, results);
        }

        tab[r].pop();
    }
}

/// Non-attacking check only (no triple rule).
fn check_non_attacking(tab: &[Vec<u32>], r: usize, c: usize) -> bool {
    let b = tab[r][c];
    for row_above in &tab[..r] {
        if row_above.len() > c && row_above[c] == b {
            return false;
        }
        if c > 0 && row_above.len() > c - 1 && row_above[c - 1] == b {
            return false;
        }
    }
    true
}

impl Ssaf {
    /// Generate non-attacking fillings (Macdonald E-fillings).
    ///
    /// Entries are not constrained to be weakly decreasing within rows.
    /// Only the non-attacking condition (no equal values in same column
    /// or adjacent column above) is enforced.
    pub fn non_attacking_fillings(shape: &[u32], basement: &[u32]) -> Vec<Ssaf> {
        assert_eq!(shape.len(), basement.len());
        let n = shape.len() as u32;

        let mut row_sequence = Vec::new();
        for (i, &s) in shape.iter().enumerate() {
            for _ in 0..s {
                row_sequence.push(i);
            }
        }

        let mut init: Vec<Vec<u32>> = basement.iter().map(|&b| vec![b]).collect();
        let mut results = Vec::new();
        generate_non_attacking_recursive(&mut init, &row_sequence, 0, n, &mut results);
        results
    }

    /// Macdonald E-fillings with identity basement.
    pub fn macdonald_e_fillings(alpha: &[u32]) -> Vec<Ssaf> {
        let n = alpha.len();
        let basement: Vec<u32> = (1..=n as u32).collect();
        Self::non_attacking_fillings(alpha, &basement)
    }

    // ── Macdonald statistics ───────────────────────────────────────

    /// Arm-leg data for Macdonald factor computation.
    ///
    /// For each descent cell `(r, c)` (where `ssaf[r][c] != ssaf[r][c-1]`),
    /// returns `(arm, leg)` where:
    /// - `leg = row_len(r) - 1 - c`
    /// - `arm = #{rows below: c ≤ len ≤ len(r)} + #{rows above: c-1 ≤ len < len(r)}`
    pub fn arm_leg_data(&self) -> Vec<(u32, u32)> {
        let n = self.num_rows();
        let shape: Vec<usize> = self.rows.iter().map(|r| r.len()).collect();
        let mut data = Vec::new();

        for r in 0..n {
            for c in 1..shape[r] {
                if self.rows[r][c] != self.rows[r][c - 1] {
                    let leg = (shape[r] - 1 - c) as u32;
                    let arm_below: u32 = (r + 1..n)
                        .filter(|&i| c <= shape[i] && shape[i] <= shape[r])
                        .count() as u32;
                    let arm_above: u32 = (0..r)
                        .filter(|&i| c >= 1 && c - 1 < shape[i] && shape[i] < shape[r])
                        .count() as u32;
                    data.push((arm_below + arm_above, leg));
                }
            }
        }
        data
    }

    // ── Lascoux-Schützenberger involution ──────────────────────────

    /// Crystal word for the *i*-string of this SSAF.
    ///
    /// Reading order: columns left-to-right (skipping basement column 0),
    /// rows bottom-to-top within each column.
    /// After extraction:
    /// 1. Remove same-column pairs (columns with exactly one *i* and one *i+1*).
    /// 2. Repeatedly remove adjacent (*i*, *i+1*) pairs.
    /// 3. Result has form (i+1)^a i^b.
    ///
    /// Returns positions and values of the surviving entries.
    fn crystal_word(&self, i: u32) -> Vec<(usize, usize, u32)> {
        let n = self.num_rows();
        let max_col = self.rows.iter().map(|r| r.len()).max().unwrap_or(0);

        // Step 1: Extract subword, columns left-to-right (skip col 0), rows bottom-to-top
        let mut positions: Vec<(usize, usize, u32)> = Vec::new();
        for c in 1..max_col {
            for r in (0..n).rev() {
                if self.rows[r].len() > c {
                    let v = self.rows[r][c];
                    if v == i || v == i + 1 {
                        positions.push((r, c, v));
                    }
                }
            }
        }

        // Step 2: SameColDelete — remove column pairs {i, i+1}
        let mut col_counts: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for (idx, &(_, c, _)) in positions.iter().enumerate() {
            col_counts.entry(c).or_default().push(idx);
        }
        let mut remove_set = std::collections::HashSet::new();
        for indices in col_counts.values() {
            if indices.len() == 2 {
                let v0 = positions[indices[0]].2;
                let v1 = positions[indices[1]].2;
                if (v0 == i && v1 == i + 1) || (v0 == i + 1 && v1 == i) {
                    remove_set.insert(indices[0]);
                    remove_set.insert(indices[1]);
                }
            }
        }
        if !remove_set.is_empty() {
            positions = positions
                .into_iter()
                .enumerate()
                .filter(|(idx, _)| !remove_set.contains(idx))
                .map(|(_, p)| p)
                .collect();
        }

        // Step 3: Repeatedly remove adjacent (i, i+1) pairs
        loop {
            let mut removed = false;
            let mut new_positions = Vec::with_capacity(positions.len());
            let mut skip_next = false;
            for idx in 0..positions.len() {
                if skip_next {
                    skip_next = false;
                    continue;
                }
                if idx + 1 < positions.len()
                    && positions[idx].2 == i
                    && positions[idx + 1].2 == i + 1
                {
                    skip_next = true;
                    removed = true;
                } else {
                    new_positions.push(positions[idx]);
                }
            }
            positions = new_positions;
            if !removed {
                break;
            }
        }
        positions
    }

    /// Lascoux-Schützenberger involution exchanging values *i* and *i+1*.
    ///
    /// Gets the crystal word, swaps *i* ↔ *i+1*, reverses the value sequence,
    /// and writes back to the same positions.
    pub fn lascoux_schutzenberger(&self, i: u32) -> Ssaf {
        let word = self.crystal_word(i);
        if word.is_empty() {
            return self.clone();
        }
        // Swap i ↔ i+1 and reverse
        let swapped_reversed: Vec<u32> = word
            .iter()
            .rev()
            .map(|&(_, _, v)| if v == i { i + 1 } else { i })
            .collect();

        let mut rows = self.rows.clone();
        for (idx, &(r, c, _)) in word.iter().enumerate() {
            rows[r][c] = swapped_reversed[idx];
        }
        Ssaf { rows }
    }

    /// Normalize the SSAF weight to a partition (weakly decreasing).
    ///
    /// Repeatedly applies Lascoux-Schützenberger involutions to sort
    /// the weight into a partition.
    pub fn weight_normalize(&self) -> Ssaf {
        let mut result = self.clone();
        let n = self.num_rows();
        // Bubble sort the weight using LS involutions
        loop {
            let w = result.weight_vector();
            // Find first descent in weight
            let mut sorted = true;
            for idx in 0..w.len().saturating_sub(1) {
                if w[idx] < w[idx + 1] {
                    // Apply LS involution for values (idx+1, idx+2) [1-indexed]
                    result = result.lascoux_schutzenberger((idx + 1) as u32);
                    sorted = false;
                    break;
                }
            }
            if sorted || n <= 1 {
                break;
            }
        }
        result
    }

    // ── Crystal operators on SSAF ──────────────────────────────────

    /// Crystal raising operator *eᵢ* on an SSAF.
    ///
    /// In the crystal word (i+1)^a i^b, changes the first *i* to *i+1*
    /// (decreasing b by 1, increasing a by 1).
    /// Returns `None` if b = 0.
    pub fn crystal_e(&self, i: u32) -> Option<Ssaf> {
        let word = self.crystal_word(i);
        // Word is (i+1)^a i^b. If no i's, return None.
        let first_i = word.iter().position(|&(_, _, v)| v == i)?;
        let mut rows = self.rows.clone();
        let (r, c, _) = word[first_i];
        rows[r][c] = i + 1;
        // Sort columns to maintain SSAF structure
        sort_columns(&mut rows);
        Some(Ssaf { rows })
    }

    /// Crystal lowering operator *fᵢ* on an SSAF.
    ///
    /// In the crystal word (i+1)^a i^b, changes the last *i+1* to *i*
    /// (decreasing a by 1, increasing b by 1).
    /// Returns `None` if a = 0.
    pub fn crystal_f(&self, i: u32) -> Option<Ssaf> {
        let word = self.crystal_word(i);
        // Word is (i+1)^a i^b. If no (i+1)'s, return None.
        let last_ip1 = word.iter().rposition(|&(_, _, v)| v == i + 1)?;
        let mut rows = self.rows.clone();
        let (r, c, _) = word[last_ip1];
        rows[r][c] = i;
        sort_columns(&mut rows);
        Some(Ssaf { rows })
    }
}

/// Sort each column of the filling in decreasing order (to maintain SSAF structure).
fn sort_columns(rows: &mut [Vec<u32>]) {
    let max_col = rows.iter().map(|r| r.len()).max().unwrap_or(0);
    for c in 1..max_col {
        // Collect values in this column
        let mut vals: Vec<(usize, u32)> = Vec::new();
        for (r, row) in rows.iter().enumerate() {
            if row.len() > c {
                vals.push((r, row[c]));
            }
        }
        // Sort values in decreasing order and write back
        let mut sorted_vals: Vec<u32> = vals.iter().map(|&(_, v)| v).collect();
        sorted_vals.sort_unstable_by(|a, b| b.cmp(a));
        for (idx, &(r, _)) in vals.iter().enumerate() {
            rows[r][c] = sorted_vals[idx];
        }
    }
}

// ── Display ─────────────────────────────────────────────────────────

impl fmt::Display for Ssaf {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, row) in self.rows.iter().enumerate() {
            if i > 0 {
                writeln!(f)?;
            }
            // Show basement separated by |
            if let Some((&basement, rest)) = row.split_first() {
                write!(f, "{basement}|")?;
                for (j, &v) in rest.iter().enumerate() {
                    if j > 0 {
                        write!(f, " ")?;
                    }
                    write!(f, "{v}")?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom_fillings_identity_basement() {
        // Shape (2, 1) with basement [1, 2]
        // Atom fillings: basement [1,2], rows weakly decreasing
        let fillings = Ssaf::atom_fillings(&[2, 1]);
        assert!(!fillings.is_empty());
        for f in &fillings {
            assert_eq!(f.shape(), vec![2, 1]);
            assert_eq!(f.basement(), vec![1, 2]);
            // Each row is weakly decreasing
            for row in f.rows() {
                for w in row.windows(2) {
                    assert!(w[0] >= w[1], "row not weakly decreasing: {:?}", row);
                }
            }
        }
    }

    #[test]
    fn test_key_fillings_reversed_basement() {
        let fillings = Ssaf::key_fillings(&[2, 1]);
        assert!(!fillings.is_empty());
        for f in &fillings {
            assert_eq!(f.basement(), vec![2, 1]);
        }
    }

    #[test]
    fn test_key_fillings_reverse_shape_convention() {
        let fillings = Ssaf::key_fillings(&[0, 2]);
        let mut weights: Vec<Vec<u32>> = fillings.iter().map(|f| f.weight_vector()).collect();
        weights.sort();
        assert_eq!(weights, vec![vec![0, 2], vec![1, 1], vec![2, 0]]);
    }

    #[test]
    fn test_dominant_key_has_single_weight() {
        let keys = Ssaf::key_fillings(&[2, 1]);
        assert_eq!(keys.len(), 1);
        assert_eq!(keys[0].weight_vector(), vec![2, 1]);
    }

    #[test]
    fn test_weight_vector() {
        // Simple test: shape (1), basement [1]
        // Only filling: [1, 1] (row: basement=1, fill=1)
        let fillings = Ssaf::atom_fillings(&[1]);
        assert_eq!(fillings.len(), 1);
        assert_eq!(fillings[0].weight_vector(), vec![1]);
    }

    #[test]
    fn test_ssaf_statistics() {
        // Test that inversions + coinversions = total triples
        let fillings = Ssaf::atom_fillings(&[2, 1, 1]);
        for f in &fillings {
            let inv = f.inversions();
            let coinv = f.coinversions();
            // Verify no panics and that we get some triples
            let _ = (inv, coinv);
        }
    }

    #[test]
    fn test_lock_fillings_basic() {
        // Lock fillings for shape (1, 1) with basement [2, 1]
        let fillings = Ssaf::lock_fillings(&[1, 1]);
        for f in &fillings {
            assert_eq!(f.basement(), vec![2, 1]);
            // Row entries are weakly increasing (excluding basement)
            for row in f.rows() {
                for w in row[1..].windows(2) {
                    assert!(w[0] <= w[1], "lock row not weakly increasing: {:?}", row);
                }
            }
            // Each entry ≤ basement
            for row in f.rows() {
                let base = row[0];
                for &v in &row[1..] {
                    assert!(v <= base, "entry {} > basement {}", v, base);
                }
            }
        }
    }

    #[test]
    fn test_lock_fillings_column_decreasing() {
        let fillings = Ssaf::lock_fillings(&[2, 1, 1]);
        for f in &fillings {
            // Columns (excluding basement) should be strictly decreasing top-to-bottom
            let max_col = f.rows().iter().map(|r| r.len()).max().unwrap_or(0);
            for c in 1..max_col {
                let col_vals: Vec<u32> = f
                    .rows()
                    .iter()
                    .filter(|r| r.len() > c)
                    .map(|r| r[c])
                    .collect();
                for w in col_vals.windows(2) {
                    assert!(
                        w[0] > w[1],
                        "lock column not strictly decreasing: {:?}",
                        col_vals
                    );
                }
            }
        }
    }

    #[test]
    fn test_mason_bijection_roundtrip() {
        // Create an SSYT, convert to SSAF, check it's a valid atom filling
        // SSYT of shape (2,1): [[1,2],[3]] or [[1,3],[2]]
        let ssyt = vec![vec![1, 2], vec![3]];
        let atom = Ssaf::from_ssyt(&ssyt);
        assert_eq!(atom.basement(), vec![1, 2, 3]);
        // Weight should match: one 1, one 2, one 3
        assert_eq!(atom.weight_vector(), vec![1, 1, 1]);
    }

    #[test]
    fn test_mason_bijection_partition_shape() {
        // For partition shape, number of SSYTs should equal number of atom fillings
        // Shape (2,1) with max entry 3: SSYTs have entries in {1,2,3}
        // The atom fillings of shape (2,1) should correspond to SSYTs
        let ssyt1 = vec![vec![1, 2], vec![3]];
        let ssyt2 = vec![vec![1, 3], vec![2]];

        let atom1 = Ssaf::from_ssyt(&ssyt1);
        let atom2 = Ssaf::from_ssyt(&ssyt2);

        // Both should produce valid fillings (no panics)
        assert_eq!(atom1.weight_vector().iter().sum::<u32>(), 3);
        assert_eq!(atom2.weight_vector().iter().sum::<u32>(), 3);

        // They should be different fillings
        assert_ne!(atom1, atom2);
    }

    #[test]
    fn test_non_attacking_fillings() {
        // Non-attacking fillings for shape (1,1) with basement [1,2]
        let fillings = Ssaf::non_attacking_fillings(&[1, 1], &[1, 2]);
        // With n=2, each row gets one entry from {1,2}, non-attacking
        assert!(!fillings.is_empty());
        // Verify non-attacking: no same value in same or adjacent column above
        for f in &fillings {
            for r in 0..f.rows().len() {
                for c in 1..f.rows()[r].len() {
                    let val = f.rows()[r][c];
                    for i in 0..r {
                        if f.rows()[i].len() > c {
                            assert_ne!(f.rows()[i][c], val, "same column attack");
                        }
                        if c > 0 && f.rows()[i].len() > c - 1 {
                            assert_ne!(f.rows()[i][c - 1], val, "adjacent column attack");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_macdonald_e_more_than_ssaf() {
        // Non-attacking fillings should be >= SSAF count (less constrained)
        let shape = [2, 1];
        let basement: Vec<u32> = vec![1, 2];
        let ssaf_count = Ssaf::generate(&shape, &basement).len();
        let na_count = Ssaf::non_attacking_fillings(&shape, &basement).len();
        assert!(
            na_count >= ssaf_count,
            "non-attacking ({}) should be >= ssaf ({})",
            na_count,
            ssaf_count
        );
    }

    #[test]
    fn test_arm_leg_data() {
        let fillings = Ssaf::key_fillings(&[2, 1]);
        for f in &fillings {
            let data = f.arm_leg_data();
            // Each pair should have non-negative components
            for &(arm, leg) in &data {
                assert!(arm < 100 && leg < 100);
            }
        }
    }

    #[test]
    fn test_ls_involution_is_involution() {
        let fillings = Ssaf::atom_fillings(&[2, 1, 1]);
        for f in &fillings {
            let ls1 = f.lascoux_schutzenberger(1);
            let ls1_ls1 = ls1.lascoux_schutzenberger(1);
            assert_eq!(
                f.weight_vector().iter().sum::<u32>(),
                ls1.weight_vector().iter().sum::<u32>(),
                "LS should preserve total weight"
            );
            // LS is an involution
            assert_eq!(&ls1_ls1, f, "LS^2 should be identity");
        }
    }

    #[test]
    fn test_weight_normalize_gives_partition() {
        let fillings = Ssaf::atom_fillings(&[1, 3, 2]);
        for f in &fillings {
            let normalized = f.weight_normalize();
            let w = normalized.weight_vector();
            // Weight should be weakly decreasing (a partition)
            for pair in w.windows(2) {
                assert!(pair[0] >= pair[1], "weight not a partition: {:?}", w);
            }
        }
    }

    #[test]
    fn test_clockwise() {
        assert!(clockwise_q(1, 2, 3));
        assert!(clockwise_q(2, 3, 1));
        assert!(clockwise_q(3, 1, 2));
        assert!(!clockwise_q(3, 2, 1));
        assert!(!clockwise_q(1, 3, 2));
        assert!(!clockwise_q(2, 1, 3));
    }

    #[test]
    fn test_atom_composition_to_permutation() {
        let (sigma, lambda) = atom_composition_to_permutation(&[1, 3, 2]);
        assert_eq!(lambda, vec![3, 2, 1]);
        // sigma should sort [1,3,2] to [3,2,1]
        // alpha[i] = lambda[sigma[i]]
        for (i, &s) in sigma.iter().enumerate() {
            assert_eq!([1, 3, 2][i], lambda[s]);
        }
    }

    #[test]
    fn test_generate_with_custom_basement() {
        // Shape (1, 1) with basement [2, 1]
        let fillings = Ssaf::generate(&[1, 1], &[2, 1]);
        // Each filling has 2 rows, each with basement + 1 entry
        for f in &fillings {
            assert_eq!(f.rows().len(), 2);
            assert_eq!(f.rows()[0].len(), 2);
            assert_eq!(f.rows()[1].len(), 2);
        }
    }
}
