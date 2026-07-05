//! Supersymmetric, or hook, Schur functions via hook tableaux.
//!
//! The alphabet is
//! `1 < 2 < ... < x_count < 1' < 2' < ... < y_count'`.
//! Rows and columns are weakly increasing, unprimed entries strictly increase
//! down columns, and primed entries strictly increase along rows.

use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum HookEntry {
    Unprimed(u32),
    Primed(u32),
}

impl HookEntry {
    pub fn is_unprimed(self) -> bool {
        matches!(self, HookEntry::Unprimed(_))
    }

    pub fn is_primed(self) -> bool {
        matches!(self, HookEntry::Primed(_))
    }

    pub fn index(self) -> u32 {
        match self {
            HookEntry::Unprimed(index) | HookEntry::Primed(index) => index,
        }
    }

    pub fn display(self) -> String {
        match self {
            HookEntry::Unprimed(index) => index.to_string(),
            HookEntry::Primed(index) => format!("{index}'"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HookTableau {
    rows: Vec<Vec<HookEntry>>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HookWeight {
    pub x: Vec<u32>,
    pub y: Vec<u32>,
}

pub type HookSchurExpansion = BTreeMap<HookWeight, u64>;

impl HookTableau {
    pub fn new(rows: Vec<Vec<HookEntry>>) -> Self {
        Self { rows }
    }

    pub fn rows(&self) -> &[Vec<HookEntry>] {
        &self.rows
    }

    pub fn is_hook_tableau(&self) -> bool {
        is_hook_tableau_rows(&self.rows)
    }

    pub fn weight(&self, x_count: u32, y_count: u32) -> HookWeight {
        let mut x = vec![0; x_count as usize];
        let mut y = vec![0; y_count as usize];
        for row in &self.rows {
            for &entry in row {
                match entry {
                    HookEntry::Unprimed(index) => {
                        assert!(
                            (1..=x_count).contains(&index),
                            "unprimed hook-tableau entry outside the x alphabet"
                        );
                        x[index as usize - 1] += 1;
                    }
                    HookEntry::Primed(index) => {
                        assert!(
                            (1..=y_count).contains(&index),
                            "primed hook-tableau entry outside the y alphabet"
                        );
                        y[index as usize - 1] += 1;
                    }
                }
            }
        }
        HookWeight { x, y }
    }

    pub fn display_rows(&self) -> String {
        self.rows
            .iter()
            .map(|row| {
                row.iter()
                    .map(|entry| entry.display())
                    .collect::<Vec<_>>()
                    .join(" & ")
            })
            .collect::<Vec<_>>()
            .join(" \\\\ ")
    }
}

pub fn hook_alphabet(x_count: u32, y_count: u32) -> Vec<HookEntry> {
    let mut alphabet = Vec::with_capacity((x_count + y_count) as usize);
    alphabet.extend((1..=x_count).map(HookEntry::Unprimed));
    alphabet.extend((1..=y_count).map(HookEntry::Primed));
    alphabet
}

pub fn enumerate_hook_tableaux(shape: &[usize], x_count: u32, y_count: u32) -> Vec<HookTableau> {
    if shape.iter().sum::<usize>() == 0 {
        return vec![HookTableau::new(shape.iter().map(|_| Vec::new()).collect())];
    }

    let alphabet = hook_alphabet(x_count, y_count);
    assert!(
        !alphabet.is_empty(),
        "a nonempty hook-tableau shape needs a nonempty alphabet"
    );
    assert!(
        shape.windows(2).all(|pair| pair[0] >= pair[1]),
        "hook-tableau shape must be a partition"
    );

    let mut current = shape
        .iter()
        .map(|&length| vec![alphabet[0]; length])
        .collect::<Vec<_>>();
    let mut results = Vec::new();
    enumerate_hook_tableaux_rec(
        shape,
        &alphabet,
        0,
        0,
        &mut current,
        &mut results,
        x_count,
        y_count,
    );
    results
}

pub fn hook_schur_expansion(shape: &[usize], x_count: u32, y_count: u32) -> HookSchurExpansion {
    let mut expansion = HookSchurExpansion::new();
    for tableau in enumerate_hook_tableaux(shape, x_count, y_count) {
        *expansion
            .entry(tableau.weight(x_count, y_count))
            .or_insert(0) += 1;
    }
    expansion
}

fn enumerate_hook_tableaux_rec(
    shape: &[usize],
    alphabet: &[HookEntry],
    row: usize,
    col: usize,
    current: &mut Vec<Vec<HookEntry>>,
    results: &mut Vec<HookTableau>,
    x_count: u32,
    y_count: u32,
) {
    if row == shape.len() {
        let tableau = HookTableau::new(current.clone());
        if tableau.is_hook_tableau() {
            debug_assert_eq!(tableau.weight(x_count, y_count).x.len(), x_count as usize);
            results.push(tableau);
        }
        return;
    }

    if shape[row] == 0 {
        enumerate_hook_tableaux_rec(
            shape,
            alphabet,
            row + 1,
            0,
            current,
            results,
            x_count,
            y_count,
        );
        return;
    }

    let (next_row, next_col) = if col + 1 == shape[row] {
        (row + 1, 0)
    } else {
        (row, col + 1)
    };

    for &entry in alphabet {
        current[row][col] = entry;
        enumerate_hook_tableaux_rec(
            shape, alphabet, next_row, next_col, current, results, x_count, y_count,
        );
    }
}

fn is_hook_tableau_rows(rows: &[Vec<HookEntry>]) -> bool {
    for row in rows {
        for pair in row.windows(2) {
            if pair[0] > pair[1] {
                return false;
            }
            if pair[0].is_primed() && pair[1].is_primed() && pair[0].index() >= pair[1].index() {
                return false;
            }
        }
    }

    for row in 0..rows.len().saturating_sub(1) {
        for col in 0..rows[row + 1].len() {
            let top = rows[row][col];
            let bottom = rows[row + 1][col];
            if top > bottom {
                return false;
            }
            if top.is_unprimed() && bottom.is_unprimed() && top.index() >= bottom.index() {
                return false;
            }
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn w(x: &[u32], y: &[u32]) -> HookWeight {
        HookWeight {
            x: x.to_vec(),
            y: y.to_vec(),
        }
    }

    #[test]
    fn hook_entry_order_matches_super_alphabet() {
        assert!(HookEntry::Unprimed(2) < HookEntry::Primed(1));
        assert!(HookEntry::Primed(1) < HookEntry::Primed(2));
    }

    #[test]
    fn hook_tableau_rules_match_site_examples() {
        let first = HookTableau::new(vec![
            vec![
                HookEntry::Unprimed(1),
                HookEntry::Unprimed(1),
                HookEntry::Primed(1),
            ],
            vec![HookEntry::Unprimed(2)],
        ]);
        let second = HookTableau::new(vec![
            vec![
                HookEntry::Unprimed(1),
                HookEntry::Unprimed(1),
                HookEntry::Unprimed(2),
            ],
            vec![HookEntry::Primed(1)],
        ]);

        assert!(first.is_hook_tableau());
        assert!(second.is_hook_tableau());
        assert_eq!(first.weight(2, 1), w(&[2, 1], &[1]));
        assert_eq!(second.weight(2, 1), w(&[2, 1], &[1]));
    }

    #[test]
    fn hook_schur_31_site_expansion() {
        let expansion = hook_schur_expansion(&[3, 1], 2, 1);
        let expected = HookSchurExpansion::from([
            (w(&[0, 2], &[2]), 1),
            (w(&[0, 3], &[1]), 1),
            (w(&[1, 1], &[2]), 1),
            (w(&[1, 2], &[1]), 2),
            (w(&[1, 3], &[0]), 1),
            (w(&[2, 0], &[2]), 1),
            (w(&[2, 1], &[1]), 2),
            (w(&[2, 2], &[0]), 1),
            (w(&[3, 0], &[1]), 1),
            (w(&[3, 1], &[0]), 1),
        ]);

        assert_eq!(expansion, expected);
    }
}
