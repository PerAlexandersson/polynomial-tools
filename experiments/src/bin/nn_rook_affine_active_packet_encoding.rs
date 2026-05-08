//! Diagnostics for the proposed active packet-row encoding.
//!
//! This is not a proof of the Ferrers statement.  It is a small finite model of
//! the Lean `DegreePath` row geometry, used to check whether an endpoint-coded
//! active path set can simultaneously provide:
//!
//! - exact entry enumeration by affine raw branches;
//! - the strict inversion-forces-intersection axiom;
//! - local select/select splice closure in the crossed matrix entries;
//! - opposite-branch key preservation under the splice.
//!
//! Run from `/workspace/rust` with:
//!
//! ```text
//! timeout 60s nice -n 10 cargo run -q -p experiments \
//!   --bin nn_rook_affine_active_packet_encoding
//! ```

use std::collections::HashMap;
use std::fmt;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
enum BranchKind {
    Base,
    Tail,
    TCorr,
    SShift,
}

impl BranchKind {
    fn all() -> [Self; 4] {
        [Self::Base, Self::Tail, Self::TCorr, Self::SShift]
    }

    fn natural_is_select(self) -> bool {
        matches!(self, Self::Tail | Self::SShift)
    }

    fn shift(self) -> usize {
        match self {
            Self::Base | Self::TCorr => 0,
            Self::Tail | Self::SShift => 1,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
struct Branch {
    kind: BranchKind,
    mask: u32,
}

impl Branch {
    fn degree(self) -> usize {
        self.mask.count_ones() as usize + self.kind.shift()
    }
}

impl fmt::Display for Branch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}({:b})", self.kind, self.mask)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
enum RowPath {
    Skip { bound: usize },
    Select { lower: usize, upper: usize },
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
struct DegreePath {
    degree: usize,
    row: RowPath,
}

impl DegreePath {
    fn source_state(self) -> (usize, usize) {
        match self.row {
            RowPath::Skip { bound } => (self.degree, bound),
            RowPath::Select { upper, .. } => (self.degree, upper),
        }
    }

    fn target_degree(self) -> usize {
        self.degree
            + match self.row {
                RowPath::Skip { .. } => 0,
                RowPath::Select { .. } => 1,
            }
    }

    fn target_state(self) -> (usize, usize) {
        match self.row {
            RowPath::Skip { bound } => (self.target_degree(), bound),
            RowPath::Select { lower, .. } => (self.target_degree(), lower),
        }
    }

    fn select_parts(self) -> Option<(usize, usize, usize)> {
        match self.row {
            RowPath::Select { lower, upper } => Some((self.degree, lower, upper)),
            RowPath::Skip { .. } => None,
        }
    }

    fn disjoint(self, other: Self) -> bool {
        if self.source_state() == other.source_state() {
            return false;
        }
        if self.target_state() == other.target_state() {
            return false;
        }
        let Some((d1, l1, u1)) = self.select_parts() else {
            return true;
        };
        let Some((d2, l2, u2)) = other.select_parts() else {
            return true;
        };
        d1 != d2 || u1 <= l2 || u2 <= l1
    }
}

#[derive(Clone, Debug)]
struct EncodedPath {
    source: usize,
    sink: usize,
    branch: Branch,
    path: DegreePath,
}

#[derive(Clone, Copy, Debug)]
struct EntryData<'a> {
    rows: &'a [usize],
    cols: &'a [usize],
}

impl EntryData<'_> {
    fn entry_degree(self, r: usize, c: usize) -> Option<usize> {
        self.rows[r].checked_sub(self.cols[c])
    }
}

fn branches(packet_rows: usize) -> Vec<Branch> {
    let mut out = Vec::new();
    for kind in BranchKind::all() {
        for mask in 0..(1u32 << packet_rows) {
            out.push(Branch { kind, mask });
        }
    }
    out
}

fn entry_branches(data: EntryData<'_>, all: &[Branch], r: usize, c: usize) -> Vec<Branch> {
    let Some(degree) = data.entry_degree(r, c) else {
        return Vec::new();
    };
    all.iter()
        .copied()
        .filter(|b| b.degree() == degree)
        .collect()
}

fn build_natural_shift_encoding(data: EntryData<'_>, all: &[Branch]) -> Vec<EncodedPath> {
    let branch_index: HashMap<Branch, usize> =
        all.iter().enumerate().map(|(i, &b)| (b, i)).collect();
    let n = data.rows.len();
    let lower_base = 1usize;
    let lower_slots = all.len() * n;
    let upper_base = lower_base + lower_slots;
    let skip_base = upper_base + n;

    let mut out = Vec::new();
    for r in 0..n {
        for c in 0..n {
            for branch in entry_branches(data, all, r, c) {
                let bidx = branch_index[&branch];
                let path = if branch.kind.natural_is_select() {
                    let lower = lower_base + bidx * n + r;
                    let upper = upper_base + c;
                    DegreePath {
                        degree: branch.degree() - 1,
                        row: RowPath::Select { lower, upper },
                    }
                } else {
                    let bound = skip_base + ((bidx * n + r) * n + c);
                    DegreePath {
                        degree: branch.degree(),
                        row: RowPath::Skip { bound },
                    }
                };
                out.push(EncodedPath {
                    source: c,
                    sink: r,
                    branch,
                    path,
                });
            }
        }
    }
    out
}

fn build_coordinate_all_select_encoding(data: EntryData<'_>, all: &[Branch]) -> Vec<EncodedPath> {
    let branch_index: HashMap<Branch, usize> =
        all.iter().enumerate().map(|(i, &b)| (b, i)).collect();
    let n = data.rows.len();
    let max_value = data
        .rows
        .iter()
        .chain(data.cols.iter())
        .copied()
        .max()
        .unwrap_or(0);
    let lower_slots = all.len() * n;
    let slot = lower_slots + n + 2;
    let upper_payload_base = lower_slots + 1;

    let coord = |value: usize| max_value - value;

    let mut out = Vec::new();
    for r in 0..n {
        for c in 0..n {
            for branch in entry_branches(data, all, r, c) {
                let bidx = branch_index[&branch];
                let lower = coord(data.rows[r]) * slot + bidx * n + r + 1;
                let upper = coord(data.cols[c]) * slot + upper_payload_base + c;
                assert!(lower < upper);
                out.push(EncodedPath {
                    source: c,
                    sink: r,
                    branch,
                    path: DegreePath {
                        degree: 0,
                        row: RowPath::Select { lower, upper },
                    },
                });
            }
        }
    }
    out
}

fn check_entry_bijection(
    name: &str,
    data: EntryData<'_>,
    all: &[Branch],
    paths: &[EncodedPath],
) -> bool {
    let n = data.rows.len();
    let mut ok = true;
    for r in 0..n {
        for c in 0..n {
            let expected = entry_branches(data, all, r, c);
            let actual: Vec<_> = paths
                .iter()
                .filter(|p| p.source == c && p.sink == r)
                .map(|p| p.branch)
                .collect();
            if expected.len() != actual.len()
                || expected.iter().any(|b| !actual.contains(b))
                || actual.iter().any(|b| !expected.contains(b))
            {
                println!(
                    "FAIL {name} entry bijection at r={r} c={c}: expected {} got {}",
                    expected.len(),
                    actual.len()
                );
                ok = false;
            }
        }
    }
    ok
}

fn first_strict_failure(paths: &[EncodedPath]) -> Option<(EncodedPath, EncodedPath)> {
    for p in paths {
        for q in paths {
            if p.source < q.source && q.sink < p.sink && p.path.disjoint(q.path) {
                return Some((p.clone(), q.clone()));
            }
        }
    }
    None
}

fn crossed_paths(p: DegreePath, q: DegreePath) -> Option<(DegreePath, DegreePath)> {
    let (pd, pl, pu) = p.select_parts()?;
    let (qd, ql, qu) = q.select_parts()?;
    if pd != qd || !(ql < pu && pl < qu) {
        return None;
    }
    Some((
        DegreePath {
            degree: pd,
            row: RowPath::Select {
                lower: ql,
                upper: pu,
            },
        },
        DegreePath {
            degree: qd,
            row: RowPath::Select {
                lower: pl,
                upper: qu,
            },
        },
    ))
}

#[derive(Default)]
struct SpliceStats {
    tested: usize,
    closed: usize,
    key_preserved: usize,
    first_failure: Option<String>,
}

fn check_splices(paths: &[EncodedPath]) -> SpliceStats {
    let by_full: HashMap<(usize, usize, DegreePath), &EncodedPath> = paths
        .iter()
        .map(|p| ((p.source, p.sink, p.path), p))
        .collect();

    let mut stats = SpliceStats::default();
    for p in paths {
        for q in paths {
            if !(p.source < q.source && q.sink < p.sink) {
                continue;
            }
            let Some((left, right)) = crossed_paths(p.path, q.path) else {
                continue;
            };
            stats.tested += 1;
            let left_key = (p.source, q.sink, left);
            let right_key = (q.source, p.sink, right);
            let Some(left_path) = by_full.get(&left_key) else {
                if stats.first_failure.is_none() {
                    stats.first_failure = Some(format!(
                        "crossed left path not active: p=({}->{}, {}, {:?}) q=({}->{}, {}, {:?}) left entry {}->{} path {:?}",
                        p.source, p.sink, p.branch, p.path,
                        q.source, q.sink, q.branch, q.path,
                        p.source, q.sink, left
                    ));
                }
                continue;
            };
            let Some(right_path) = by_full.get(&right_key) else {
                if stats.first_failure.is_none() {
                    stats.first_failure = Some(format!(
                        "crossed right path not active: p=({}->{}, {}, {:?}) q=({}->{}, {}, {:?}) right entry {}->{} path {:?}",
                        p.source, p.sink, p.branch, p.path,
                        q.source, q.sink, q.branch, q.path,
                        q.source, p.sink, right
                    ));
                }
                continue;
            };
            stats.closed += 1;
            if left_path.branch == q.branch && right_path.branch == p.branch {
                stats.key_preserved += 1;
            } else if stats.first_failure.is_none() {
                stats.first_failure = Some(format!(
                    "branch keys changed: p=({}->{}, {}) q=({}->{}, {}) left got {} right got {}",
                    p.source,
                    p.sink,
                    p.branch,
                    q.source,
                    q.sink,
                    q.branch,
                    left_path.branch,
                    right_path.branch
                ));
            }
        }
    }
    stats
}

fn print_model(name: &str, data: EntryData<'_>, all: &[Branch], paths: Vec<EncodedPath>) {
    println!("=== {name} ===");
    println!("active paths: {}", paths.len());
    let entry_ok = check_entry_bijection(name, data, all, &paths);
    println!(
        "entry enumeration: {}",
        if entry_ok { "PASS" } else { "FAIL" }
    );

    match first_strict_failure(&paths) {
        None => println!("strict inversion axiom: PASS"),
        Some((p, q)) => {
            println!("strict inversion axiom: FAIL");
            println!(
                "  disjoint inverted pair: p source={} sink={} branch={} path={:?}",
                p.source, p.sink, p.branch, p.path
            );
            println!(
                "                          q source={} sink={} branch={} path={:?}",
                q.source, q.sink, q.branch, q.path
            );
        }
    }

    let splices = check_splices(&paths);
    println!(
        "select/select crossed splices: tested={}, closed={}, key_preserved={}",
        splices.tested, splices.closed, splices.key_preserved
    );
    if splices.tested == 0 {
        println!("  splice closure/key check: VACUOUS (no select/select inverted cases)");
    } else if let Some(failure) = splices.first_failure {
        println!("  first splice failure: {failure}");
    } else {
        println!("  splice closure/key check: PASS");
    }
    println!();
}

fn main() {
    let rows = vec![1, 2, 3];
    let cols = vec![0, 1, 2];
    let data = EntryData {
        rows: &rows,
        cols: &cols,
    };
    let all = branches(3);

    println!("rows={rows:?} cols={cols:?}");
    println!("raw affine branches: {}", all.len());
    println!();

    print_model(
        "natural shift encoding: skip=base/tCorrection, select=tail/sShifted",
        data,
        &all,
        build_natural_shift_encoding(data, &all),
    );

    print_model(
        "coordinate all-select encoding: all branches as intervals",
        data,
        &all,
        build_coordinate_all_select_encoding(data, &all),
    );
}
