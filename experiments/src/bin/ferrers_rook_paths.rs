//! Explore the "union of directed paths" statistic for standard rook placements
//! on general Ferrers boards.
//!
//! For a Ferrers shape mu = (mu_1 >= ... >= mu_ell), let
//!   n(mu) = max_i (i + mu_i).
//! Reverse columns inside the width-(n-1) rectangle. A rook in row i and
//! original column c becomes a directed edge
//!   i -> (n - c + 1).
//! Since c <= mu_i, we have n - c + 1 >= n - mu_i + 1 > i, so every edge points
//! to the right. Thus every standard rook placement becomes a directed graph on
//! [n] with indegree/outdegree at most 1, hence a disjoint union of directed
//! paths.
//!
//! We test the generating polynomial
//!   P_{mu,j}(t) = sum_{rho} t^{rch_j(rho)},
//! where rch_j(rho) is the number of path components with at least j vertices.

use std::collections::BTreeMap;

use combpoly::rook_placements::rook_polynomial;
use polynomial_tools::{check_weak_interlacing, format_poly, is_real_rooted};

type Placement = Vec<(usize, usize)>;

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
            let mut len = 1;
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

fn long_path_polynomial(mu: &[usize], j: usize) -> Vec<i64> {
    let placements = standard_placements(mu);
    let mut coeffs: BTreeMap<usize, i64> = BTreeMap::new();

    for placement in &placements {
        let count = component_lengths(mu, placement)
            .into_iter()
            .filter(|&len| len >= j)
            .count();
        *coeffs.entry(count).or_insert(0) += 1;
    }

    let degree = coeffs.keys().copied().max().unwrap_or(0);
    (0..=degree)
        .map(|k| coeffs.get(&k).copied().unwrap_or(0))
        .collect()
}

fn degree(poly: &[i64]) -> Option<usize> {
    poly.iter().rposition(|&c| c != 0)
}

fn row_delete_interlaces(smaller: &[i64], larger: &[i64]) -> Option<bool> {
    let ds = degree(smaller)?;
    let dl = degree(larger)?;
    if ds <= dl && dl <= ds + 1 {
        check_weak_interlacing(smaller, larger)
    } else {
        None
    }
}

fn has_long_chain_3(mu: &[usize]) -> bool {
    standard_placements(mu)
        .iter()
        .any(|placement| component_lengths(mu, placement).into_iter().any(|len| len >= 3))
}

fn main() {
    let max_cells = std::env::args()
        .nth(1)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(12);
    let max_j = std::env::args()
        .nth(2)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(4);

    println!("=== Ferrers rook path-component scan ===\n");
    println!("Model:");
    println!("  n(mu) = max_i (i + mu_i)");
    println!("  rook (i,c) becomes edge i -> n(mu)-c+1 after global column reversal");
    println!("  P_(mu,j)(t) = sum_rho t^(# path components with at least j vertices)\n");

    // Sanity checks against the staircase/Touchard examples already computed.
    let stair = vec![3, 2, 1];
    assert_eq!(long_path_polynomial(&stair, 2), vec![1, 11, 3]);
    assert_eq!(long_path_polynomial(&stair, 3), vec![10, 5]);
    println!("Sanity check: delta_3 gives");
    println!("  j=2: {}", format_poly(&long_path_polynomial(&stair, 2)));
    println!("  j=3: {}\n", format_poly(&long_path_polynomial(&stair, 3)));

    let mut shapes = Vec::new();
    for cells in 1..=max_cells {
        for mu in partitions(cells) {
            shapes.push(mu);
        }
    }

    let nontrivial_chain_shapes = shapes.iter().filter(|mu| has_long_chain_3(mu)).count();
    println!(
        "Shapes with |mu| <= {}: {} total, {} admit a path of length >= 2 edges\n",
        max_cells,
        shapes.len(),
        nontrivial_chain_shapes
    );

    for j in 2..=max_j {
        let mut rr_total = 0usize;
        let mut rr_pass = 0usize;
        let mut rr_fail = Vec::new();

        let mut interlace_total = 0usize;
        let mut interlace_pass = 0usize;
        let mut interlace_ineligible = 0usize;
        let mut interlace_fail = Vec::new();

        let mut same_as_standard = 0usize;
        let mut constant_polys = 0usize;

        for mu in &shapes {
            let poly = long_path_polynomial(mu, j);
            rr_total += 1;
            if poly.len() <= 2 || is_real_rooted(&poly) {
                rr_pass += 1;
            } else if rr_fail.len() < 12 {
                rr_fail.push((mu.clone(), poly.clone()));
            }

            if poly.iter().skip(1).all(|&c| c == 0) {
                constant_polys += 1;
            }

            if j == 2 && poly == rook_polynomial(mu) {
                same_as_standard += 1;
            }

            if mu.len() >= 2 {
                let smaller = long_path_polynomial(&mu[..mu.len() - 1], j);
                match row_delete_interlaces(&smaller, &poly) {
                    Some(true) => {
                        interlace_total += 1;
                        interlace_pass += 1;
                    }
                    Some(false) => {
                        interlace_total += 1;
                        if interlace_fail.len() < 12 {
                            interlace_fail.push((mu.clone(), smaller.clone(), poly.clone()));
                        }
                    }
                    None => interlace_ineligible += 1,
                }
            }
        }

        println!("--- j = {} ---", j);
        println!("Real-rooted: {}/{}", rr_pass, rr_total);
        if j == 2 {
            println!(
                "Same as standard rook polynomial: {}/{}",
                same_as_standard, rr_total
            );
        }
        println!("Constant polynomials: {}/{}", constant_polys, rr_total);
        println!(
            "Row-deletion interlacing: {}/{} passes ({} ineligible)",
            interlace_pass, interlace_total, interlace_ineligible
        );

        if !rr_fail.is_empty() {
            println!("First non-real-rooted examples:");
            for (mu, poly) in &rr_fail {
                println!("  mu={:?} -> {}", mu, format_poly(poly));
            }
        }

        if !interlace_fail.is_empty() {
            println!("First row-deletion interlacing failures:");
            for (mu, smaller, larger) in &interlace_fail {
                println!(
                    "  mu={:?}: smaller={} ; larger={}",
                    mu,
                    format_poly(smaller),
                    format_poly(larger)
                );
            }
        }

        println!();
    }

    println!("Sample shapes:");
    for mu in [
        vec![4, 3, 2, 1],
        vec![4, 4, 2, 1],
        vec![4, 4, 4],
        vec![5, 3, 2],
        vec![5, 5, 1],
    ] {
        if mu.iter().sum::<usize>() > max_cells {
            continue;
        }
        println!("  mu={:?}, n(mu)={}", mu, canonical_n(&mu));
        for j in 2..=max_j {
            let poly = long_path_polynomial(&mu, j);
            println!("    j={}: {}", j, format_poly(&poly));
        }
    }
}
