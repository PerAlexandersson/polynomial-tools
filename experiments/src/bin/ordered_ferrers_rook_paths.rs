//! Ordered path-component model for rook placements on Ferrers boards.
//!
//! For a Ferrers shape mu, we use the same directed graph D_mu(rho) as in
//! ferrers_rook_paths.rs:
//!   N(mu) = max_i (i + mu_i),
//!   rook (i,c) becomes edge i -> N(mu)-c+1 after global column reversal.
//!
//! Every standard rook placement rho therefore yields a disjoint union of
//! directed paths on [N(mu)].
//!
//! We refine a placement by:
//!   - nest(rho): number of nesting pairs of rooks,
//!   - kappa_mu(rho): number of path components of D_mu(rho),
//!   - rch_j^mu(rho): number of path components having at least j vertices.
//!
//! The "ordered" version counts a total order on the path components, hence a
//! multiplicative factor kappa_mu(rho)!.
//!
//! On the staircase mu = delta_{n-1}, the q=1 specialization reproduces the
//! ordered set-partition big-block polynomials.

use std::collections::BTreeMap;

use polynomial_tools::{format_poly, is_real_rooted};

type Placement = Vec<(usize, usize)>;

#[derive(Clone, Copy, Debug)]
enum QMode {
    Zero,
    One,
}

impl QMode {
    fn label(self) -> &'static str {
        match self {
            Self::Zero => "q=0",
            Self::One => "q=1",
        }
    }

    fn keep(self, nestings: usize) -> bool {
        match self {
            Self::Zero => nestings == 0,
            Self::One => true,
        }
    }
}

fn partitions(n: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut buf = Vec::new();
    partitions_rec(n, n, &mut buf, &mut result);
    result
}

fn partitions_rec(n: usize, max_part: usize, buf: &mut Vec<usize>, result: &mut Vec<Vec<usize>>) {
    if n == 0 {
        result.push(buf.clone());
        return;
    }
    for k in (1..=n.min(max_part)).rev() {
        buf.push(k);
        partitions_rec(n - k, k, buf, result);
        buf.pop();
    }
}

fn staircase(n: usize) -> Vec<usize> {
    (1..n).rev().collect()
}

fn factorials(max_n: usize) -> Vec<i128> {
    let mut fact = vec![1i128; max_n + 1];
    for n in 1..=max_n {
        fact[n] = fact[n - 1] * n as i128;
    }
    fact
}

fn canonical_n(mu: &[usize]) -> usize {
    mu.iter()
        .enumerate()
        .map(|(i, &part)| (i + 1) + part)
        .max()
        .unwrap_or(0)
}

fn standard_placements(mu: &[usize]) -> Vec<Placement> {
    let mut used_cols = vec![false; mu.first().copied().unwrap_or(0) + 1];
    let mut current = Vec::new();
    let mut out = Vec::new();
    rec_placements(mu, 0, &mut used_cols, &mut current, &mut out);
    out
}

fn rec_placements(
    mu: &[usize],
    row_idx: usize,
    used_cols: &mut [bool],
    current: &mut Placement,
    out: &mut Vec<Placement>,
) {
    if row_idx == mu.len() {
        out.push(current.clone());
        return;
    }

    rec_placements(mu, row_idx + 1, used_cols, current, out);

    let row = row_idx + 1;
    for col in 1..=mu[row_idx] {
        if !used_cols[col] {
            used_cols[col] = true;
            current.push((row, col));
            rec_placements(mu, row_idx + 1, used_cols, current, out);
            current.pop();
            used_cols[col] = false;
        }
    }
}

fn count_nestings(placement: &Placement) -> usize {
    let mut count = 0;
    for i in 0..placement.len() {
        for j in (i + 1)..placement.len() {
            let (r1, c1) = placement[i];
            let (r2, c2) = placement[j];
            if r1 < r2 && c1 < c2 {
                count += 1;
            }
        }
    }
    count
}

fn component_lengths(mu: &[usize], placement: &Placement) -> Vec<usize> {
    let n = canonical_n(mu);
    let mut next = vec![None; n + 1];
    let mut prev = vec![None; n + 1];

    for &(row, col) in placement {
        let target = n - col + 1;
        assert!(target > row);
        next[row] = Some(target);
        prev[target] = Some(row);
    }

    let mut lengths = Vec::new();
    for start in 1..=n {
        if prev[start].is_none() {
            let mut len = 1usize;
            let mut v = start;
            while let Some(w) = next[v] {
                len += 1;
                v = w;
            }
            lengths.push(len);
        }
    }
    lengths
}

fn ordered_path_poly(mu: &[usize], j: usize, mode: QMode) -> Vec<i128> {
    let n = canonical_n(mu);
    let fact = factorials(n);
    let placements = standard_placements(mu);
    let mut coeffs: BTreeMap<usize, i128> = BTreeMap::new();

    for placement in &placements {
        let nestings = count_nestings(placement);
        if !mode.keep(nestings) {
            continue;
        }

        let lengths = component_lengths(mu, placement);
        let big = lengths.iter().filter(|&&len| len >= j).count();
        let components = lengths.len();
        assert_eq!(components, n - placement.len());
        *coeffs.entry(big).or_insert(0) += fact[components];
    }

    let degree = coeffs.keys().copied().max().unwrap_or(0);
    (0..=degree)
        .map(|k| coeffs.get(&k).copied().unwrap_or(0))
        .collect()
}

fn binom(n: usize, k: usize) -> i128 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut num: i128 = 1;
    let mut den: i128 = 1;
    for i in 0..k {
        num *= (n - i) as i128;
        den *= (i + 1) as i128;
    }
    num / den
}

fn ordered_set_partition_poly(n: usize, j: usize) -> Vec<i128> {
    let facts = factorials(n);
    let mut dp = vec![Vec::<Vec<i128>>::new(); n + 1];
    dp[0] = vec![vec![1]];

    for cur_n in 0..n {
        let mut next = vec![vec![0i128; cur_n + 2]; cur_n + 2];
        for k in 0..=cur_n {
            let choose = binom(cur_n, k);
            let prev_n = cur_n - k;
            let big_inc = usize::from(k + 1 >= j);
            for m_prev in 0..dp[prev_n].len() {
                for b_prev in 0..dp[prev_n][m_prev].len() {
                    let val = dp[prev_n][m_prev][b_prev];
                    if val == 0 {
                        continue;
                    }
                    next[m_prev + 1][b_prev + big_inc] += choose * val;
                }
            }
        }
        while next.last().is_some_and(|row| row.iter().all(|&x| x == 0)) {
            next.pop();
        }
        for row in &mut next {
            while row.last().is_some_and(|&x| x == 0) {
                row.pop();
            }
        }
        dp[cur_n + 1] = next;
    }

    let mut poly = vec![0i128; n + 1];
    for m in 0..dp[n].len() {
        for b in 0..dp[n][m].len() {
            poly[b] += facts[m] * dp[n][m][b];
        }
    }
    while poly.last().is_some_and(|&x| x == 0) {
        poly.pop();
    }
    poly
}

fn to_i64_poly(poly: &[i128]) -> Option<Vec<i64>> {
    poly.iter()
        .copied()
        .map(|c| i64::try_from(c).ok())
        .collect()
}

fn main() {
    let max_cells = std::env::args()
        .nth(1)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(14);
    let max_j = std::env::args()
        .nth(2)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(4);

    println!("=== Ordered Ferrers rook-path model ===\n");
    println!("For a placement rho on mu:");
    println!("  D_mu(rho) is the reversed-column path digraph on [N(mu)]");
    println!("  kappa_mu(rho) = # path components = N(mu) - |rho|");
    println!("  rch_j^mu(rho) = # path components with at least j vertices");
    println!("  ordered weight = kappa_mu(rho)!");
    println!("  total weight = q^(nest(rho)) * kappa_mu(rho)! * t^(rch_j^mu(rho))\n");

    println!("=== Staircase sanity check: q=1 should match ordered set partitions ===\n");
    for j in 2..=max_j.min(5) {
        let mut all_ok = true;
        for n in 1..=8 {
            let rook_poly = ordered_path_poly(&staircase(n), j, QMode::One);
            let part_poly = ordered_set_partition_poly(n, j);
            if rook_poly != part_poly {
                all_ok = false;
                println!(
                    "FAIL: j={} n={} rook={} part={}",
                    j,
                    n,
                    format_poly(&to_i64_poly(&rook_poly).unwrap_or_default()),
                    format_poly(&to_i64_poly(&part_poly).unwrap_or_default())
                );
                break;
            }
        }
        println!("j={}: {}", j, if all_ok { "all match ✓" } else { "FAIL" });
    }

    println!("\n=== Sample staircase rows ===\n");
    for j in 2..=max_j.min(4) {
        println!("j={j}:");
        for n in 1..=6 {
            let q1 = ordered_path_poly(&staircase(n), j, QMode::One);
            let q0 = ordered_path_poly(&staircase(n), j, QMode::Zero);
            println!(
                "  n={}: q=1 {} ; q=0 {}",
                n,
                format_poly(&to_i64_poly(&q1).unwrap_or_default()),
                format_poly(&to_i64_poly(&q0).unwrap_or_default())
            );
        }
        println!();
    }

    let mut shapes = Vec::new();
    for cells in 1..=max_cells {
        for mu in partitions(cells) {
            shapes.push(mu);
        }
    }

    println!(
        "=== Real-rootedness scan over Ferrers shapes with |mu| <= {} ===\n",
        max_cells
    );
    for mode in [QMode::Zero, QMode::One] {
        println!("--- {} ---", mode.label());
        for j in 2..=max_j {
            let mut tested = 0usize;
            let mut pass = 0usize;
            let mut first_fail = None;
            let mut overflow = 0usize;

            for mu in &shapes {
                let poly = ordered_path_poly(mu, j, mode);
                let Some(poly64) = to_i64_poly(&poly) else {
                    overflow += 1;
                    continue;
                };
                tested += 1;
                if poly64.len() <= 2 || is_real_rooted(&poly64) {
                    pass += 1;
                } else if first_fail.is_none() {
                    first_fail = Some((mu.clone(), poly64));
                }
            }

            println!(
                "j={}: real-rooted {}/{} ({} skipped for size)",
                j, pass, tested, overflow
            );
            if let Some((mu, poly)) = first_fail {
                println!("  first failure: mu={:?} -> {}", mu, format_poly(&poly));
            }
        }
        println!();
    }
}
