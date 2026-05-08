//! Symbolic obstruction for raw affine branch pair swaps.
//!
//! The packet-row attempt failed for pointwise branch-key preservation.  This
//! diagnostic asks whether the weaker condition needed by an LGV involution,
//! product preservation for the two swapped paths, can hold at the raw affine
//! branch level.
//!
//! We treat every raw branch weight as a symbolic monomial atom:
//! `B_J`, `T_J`, `t Q_J`, or `s Q_J`.  Since the abstract Lean packet data has
//! no algebraic relation between these atoms, a product-preserving raw swap
//! must preserve the multiset of branch atoms.  The check below shows that
//! strict Toeplitz inversions change the required degree multiset, so this is
//! impossible before introducing a finer path-word model.

use std::collections::BTreeMap;
use std::fmt;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
enum Kind {
    Base,
    Tail,
    TCorr,
    SShift,
}

impl Kind {
    fn all() -> [Self; 4] {
        [Self::Base, Self::Tail, Self::TCorr, Self::SShift]
    }

    fn shift(self) -> usize {
        match self {
            Self::Base | Self::TCorr => 0,
            Self::Tail | Self::SShift => 1,
        }
    }

    fn parameter(self) -> Parameter {
        match self {
            Self::Base | Self::Tail => Parameter::None,
            Self::TCorr => Parameter::T,
            Self::SShift => Parameter::S,
        }
    }

    fn atom_family(self) -> AtomFamily {
        match self {
            Self::Base => AtomFamily::Base,
            Self::Tail => AtomFamily::Tail,
            Self::TCorr | Self::SShift => AtomFamily::Correction,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
enum Parameter {
    None,
    S,
    T,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
enum AtomFamily {
    Base,
    Tail,
    Correction,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct Branch {
    kind: Kind,
    mask: u32,
}

impl Branch {
    fn degree(self) -> usize {
        self.mask.count_ones() as usize + self.kind.shift()
    }

    fn atom(self) -> Atom {
        Atom {
            family: self.kind.atom_family(),
            mask: self.mask,
            parameter: self.kind.parameter(),
        }
    }
}

impl fmt::Display for Branch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}({:b})", self.kind, self.mask)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct Atom {
    family: AtomFamily,
    mask: u32,
    parameter: Parameter,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct ProductSig(BTreeMap<Atom, usize>);

impl ProductSig {
    fn of_pair(a: Branch, b: Branch) -> Self {
        let mut map = BTreeMap::new();
        *map.entry(a.atom()).or_insert(0) += 1;
        *map.entry(b.atom()).or_insert(0) += 1;
        Self(map)
    }
}

fn branches(packet_rows: usize) -> Vec<Branch> {
    let mut out = Vec::new();
    for kind in Kind::all() {
        for mask in 0..(1u32 << packet_rows) {
            out.push(Branch { kind, mask });
        }
    }
    out
}

fn by_degree(branches: &[Branch]) -> BTreeMap<usize, Vec<Branch>> {
    let mut out: BTreeMap<usize, Vec<Branch>> = BTreeMap::new();
    for &branch in branches {
        out.entry(branch.degree()).or_default().push(branch);
    }
    out
}

fn product_match_possible(
    by_deg: &BTreeMap<usize, Vec<Branch>>,
    left: Branch,
    right: Branch,
    crossed_left_degree: usize,
    crossed_right_degree: usize,
) -> bool {
    let sig = ProductSig::of_pair(left, right);
    let Some(left_candidates) = by_deg.get(&crossed_left_degree) else {
        return false;
    };
    let Some(right_candidates) = by_deg.get(&crossed_right_degree) else {
        return false;
    };
    left_candidates.iter().any(|&l| {
        right_candidates
            .iter()
            .any(|&r| ProductSig::of_pair(l, r) == sig)
    })
}

#[derive(Clone, Debug)]
struct Failure {
    rows: Vec<usize>,
    cols: Vec<usize>,
    r_hi: usize,
    r_lo: usize,
    c_lo: usize,
    c_hi: usize,
    d_hi: usize,
    d_lo: usize,
    d_cross_lo: usize,
    d_cross_hi: usize,
    left: Branch,
    right: Branch,
}

fn strictly_increasing_tuples(
    len: usize,
    max_value: usize,
    start: usize,
    current: &mut Vec<usize>,
    out: &mut Vec<Vec<usize>>,
) {
    if current.len() == len {
        out.push(current.clone());
        return;
    }
    for x in start..=max_value {
        current.push(x);
        strictly_increasing_tuples(len, max_value, x + 1, current, out);
        current.pop();
    }
}

fn all_strict(len: usize, max_value: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    strictly_increasing_tuples(len, max_value, 0, &mut Vec::new(), &mut out);
    out
}

fn scan(packet_rows: usize, n: usize, max_value: usize) -> (usize, usize, Option<Failure>) {
    let branches = branches(packet_rows);
    let by_deg = by_degree(&branches);
    let rows_list = all_strict(n, max_value);
    let cols_list = all_strict(n, max_value);

    let mut cases = 0usize;
    let mut failures = 0usize;
    let mut first_failure = None;

    for rows in &rows_list {
        for cols in &cols_list {
            for r_hi in 1..n {
                for r_lo in 0..r_hi {
                    for c_lo in 0..n - 1 {
                        for c_hi in c_lo + 1..n {
                            let Some(d_hi) = rows[r_hi].checked_sub(cols[c_lo]) else {
                                continue;
                            };
                            let Some(d_lo) = rows[r_lo].checked_sub(cols[c_hi]) else {
                                continue;
                            };
                            let d_cross_lo = rows[r_lo] - cols[c_lo];
                            let d_cross_hi = rows[r_hi] - cols[c_hi];
                            let Some(lefts) = by_deg.get(&d_hi) else {
                                continue;
                            };
                            let Some(rights) = by_deg.get(&d_lo) else {
                                continue;
                            };
                            for &left in lefts {
                                for &right in rights {
                                    cases += 1;
                                    let ok = product_match_possible(
                                        &by_deg, left, right, d_cross_lo, d_cross_hi,
                                    );
                                    if !ok {
                                        failures += 1;
                                        if first_failure.is_none() {
                                            first_failure = Some(Failure {
                                                rows: rows.clone(),
                                                cols: cols.clone(),
                                                r_hi,
                                                r_lo,
                                                c_lo,
                                                c_hi,
                                                d_hi,
                                                d_lo,
                                                d_cross_lo,
                                                d_cross_hi,
                                                left,
                                                right,
                                            });
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    (cases, failures, first_failure)
}

fn main() {
    let packet_rows = 4;
    let n = 3;
    let max_value = 6;

    let (cases, failures, first_failure) = scan(packet_rows, n, max_value);
    println!("=== Raw affine branch pair-swap diagnostic ===");
    println!("packet_rows={packet_rows}, minor size n={n}, strict values in 0..={max_value}");
    println!("inverted raw branch pairs tested: {cases}");
    println!("product-preserving raw swap failures: {failures}");

    if let Some(f) = first_failure {
        println!("\nfirst failure:");
        println!("  rows={:?}", f.rows);
        println!("  cols={:?}", f.cols);
        println!(
            "  inversion: source {}<{}, sink {}>{}",
            f.c_lo, f.c_hi, f.r_hi, f.r_lo
        );
        println!("  original degrees: left={} right={}", f.d_hi, f.d_lo);
        println!(
            "  crossed degrees: left={} right={}",
            f.d_cross_lo, f.d_cross_hi
        );
        println!("  original branches: left={} right={}", f.left, f.right);
        println!(
            "  degree order check: {} > max({}, {}) and {} < min({}, {})",
            f.d_hi, f.d_cross_lo, f.d_cross_hi, f.d_lo, f.d_cross_lo, f.d_cross_hi
        );
    }

    println!("\ninterpretation:");
    println!("  A raw branch product is a product of independent symbolic atoms.");
    println!("  Preserving such a product forces the same two raw branch atoms.");
    println!("  Strict Toeplitz endpoint crossing changes the degree multiset.");
    println!("  Therefore raw-branch pair swaps cannot be the LGV involution.");
}
