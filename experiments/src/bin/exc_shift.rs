//! Explore shifted excedance polynomials F^(s) on Ferrers boards.
//! M_{i,j} = t if j > i+s, 1 if j <= i+s, 0 if j > lambda_i.
//! s = 0: excedance. s = -1: essentially ascent.
//! Row-deletion recursion: F^(s)(lambda) -> F^(s+1)(lambda') where lambda' = lambda minus first row.
use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly};
use std::collections::BTreeSet;

fn pt(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
    v
}
fn pz(p: &[i64]) -> bool {
    p.iter().all(|&c| c == 0)
}
fn pa(a: &[i64], b: &[i64]) -> Vec<i64> {
    let l = a.len().max(b.len());
    let mut r = vec![0i64; l];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] += v;
    }
    pt(&r)
}
fn pmt(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; p.len() + 1];
    for (i, &v) in p.iter().enumerate() {
        r[i + 1] = v;
    }
    pt(&r)
}
fn pmul(a: &[i64], b: &[i64]) -> Vec<i64> {
    if a.is_empty() || b.is_empty() {
        return vec![0];
    }
    let mut r = vec![0i64; a.len() + b.len() - 1];
    for (i, &av) in a.iter().enumerate() {
        for (j, &bv) in b.iter().enumerate() {
            r[i + j] += av * bv;
        }
    }
    pt(&r)
}

/// Compute permanent of an n x m matrix with polynomial entries.
/// mat[i][j] is a polynomial (Vec<i64>).
fn permanent(mat: &Vec<Vec<Vec<i64>>>) -> Vec<i64> {
    let n = mat.len();
    if n == 0 {
        return vec![1];
    }
    let m = mat[0].len();
    // Expand along first row
    let mut result = vec![0i64];
    for j in 0..m {
        if pz(&mat[0][j]) {
            continue;
        }
        // Build submatrix: remove row 0 and column j
        let mut sub = Vec::new();
        for i in 1..n {
            let mut row = Vec::new();
            for jj in 0..m {
                if jj == j {
                    continue;
                }
                row.push(mat[i][jj].clone());
            }
            sub.push(row);
        }
        let sub_perm = permanent(&sub);
        let term = pmul(&mat[0][j], &sub_perm);
        result = pa(&result, &term);
    }
    pt(&result)
}

/// Build the shifted excedance matrix for board lambda with shift s.
/// Rows are 1-indexed (position i), columns are 1-indexed (value j).
/// M_{i,j} = t if j > i+s, 1 if 1 <= j <= i+s, 0 if j > lambda_i.
fn shifted_exc_matrix(lambda: &[u8], s: i32) -> Vec<Vec<Vec<i64>>> {
    let n = lambda.len();
    let m = lambda[0] as usize; // max column
    let mut mat = Vec::new();
    for i in 0..n {
        let mut row = Vec::new();
        let boundary = (i as i32 + 1 + s) as i64; // i+1 because 1-indexed
        for j in 0..m {
            let j1 = (j + 1) as i64; // 1-indexed column
            if j1 > lambda[i] as i64 {
                row.push(vec![0i64]); // outside board
            } else if j1 <= boundary {
                row.push(vec![1i64]); // non-excedance: 1
            } else {
                row.push(vec![0, 1]); // excedance: t
            }
        }
        mat.push(row);
    }
    mat
}

/// Compute F^(s)_j(lambda) = contribution from first entry j to F^(s)(lambda).
fn shifted_exc_by_first(lambda: &[u8], s: i32) -> Vec<Vec<i64>> {
    let n = lambda.len();
    let m = lambda[0] as usize;
    let mut result = vec![vec![0i64]; m + 1]; // result[j] for j=1..m
    if n == 0 {
        return result;
    }

    // Place first entry = j (column j in row 1)
    let boundary = (1i32 + s) as i64; // boundary for row 1
    for j in 1..=m {
        if (j as u8) > lambda[0] {
            continue;
        }
        let factor: Vec<i64> = if (j as i64) <= boundary {
            vec![1]
        } else {
            vec![0, 1]
        }; // 1 or t
           // Build submatrix: rows 2..n, columns {1..m}\{j}, shift s+1
           // For remaining row i (original position i, 2-indexed), column j':
           // Entry = t if j' > i + s, 1 if j' <= i + s.
           // But we removed column j and are now at shift s (NOT s+1) because
           // we haven't relabeled rows! The remaining rows are still at positions 2..n.
           // Actually: the permanent of the submatrix IS F^(s) restricted to rows 2..n.
           // But the row indices are 2..n with boundaries at i+s.
           // If we relabel: new row index i' = i-1, boundary at i'+1+s = i+s. Same shift!
           // Wait: new position i' = i-1, boundary should be i'+1+s_new. We need i'+1+s_new = i+s.
           // i' = i-1, so i'+1 = i, and i + s_new = i + s => s_new = s. NO SHIFT!

        // Hmm, let me reconsider. The permanent is just:
        // sum over assignments sigma(2),...,sigma(n) of (1,...,m)\{j}
        // with weight t^{sum of [sigma(i) > i+s] for i=2..n}.
        // If we relabel position i -> i' = i-1: weight is t^{sum [sigma(i'+1) > (i'+1)+s]}.
        // This is exc with shift s relative to positions 2..n.
        // Relabeled as positions 1..n-1: boundary at i'+1+s for new position i'.
        // So the shift for the subboard is s (NOT s+1). Wait, that's the same...

        // Let me just compute directly.
        let mut sub = Vec::new();
        for i in 1..n {
            // original rows 2..n (0-indexed: 1..n-1)
            let orig_pos = i + 1; // 1-indexed position
            let bnd = (orig_pos as i64) + (s as i64);
            let mut row = Vec::new();
            for jj in 1..=m {
                if jj == j {
                    continue;
                }
                if (jj as u8) > lambda[i] {
                    row.push(vec![0i64]);
                } else if (jj as i64) <= bnd {
                    row.push(vec![1i64]);
                } else {
                    row.push(vec![0, 1]);
                }
            }
            sub.push(row);
        }
        let sub_perm = permanent(&sub);
        result[j] = pmul(&factor, &sub_perm);
    }
    result
}

fn is_rr(p: &[i64]) -> bool {
    let p = pt(p);
    if p.len() <= 2 {
        return true;
    }
    polynomial_tools::real_rootedness::is_real_rooted(&p)
}

fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(5);

    println!("=== Shifted excedance polynomials F^(s)(lambda) ===\n");

    // Test on full boards [n,n,...,n] = S_n
    for n in 2..=max_n {
        let board: Vec<u8> = vec![n as u8; n];
        println!("--- Board {:?} (S_{}) ---", board, n);

        for s in -1..=(n as i32) {
            let mat = shifted_exc_matrix(&board, s);
            let total = permanent(&mat);
            let rr = is_rr(&total);
            println!("  s={:2}: F = {}  RR={}", s, format_poly(&pt(&total)), rr);

            // Also compute F_j (by first entry)
            let fj = shifted_exc_by_first(&board, s);
            let mut all_rr = true;
            for j in 1..=n {
                if !pz(&fj[j]) && !is_rr(&fj[j]) {
                    all_rr = false;
                }
            }
            if !all_rr {
                print!("    F_j not all RR:");
                for j in 1..=n {
                    if !pz(&fj[j]) {
                        print!(" F_{}={}", j, format_poly(&pt(&fj[j])));
                    }
                }
                println!();
            }
        }
        println!();
    }

    // Check: does deleting row 1 with first entry j give shift s or s+1?
    println!(
        "=== Recursion check: does F^(s)(lambda, first=j) use F^(s) or F^(s+1) on sublambda? ==="
    );
    let board: Vec<u8> = vec![4, 4, 4, 4];
    let sublambda: Vec<u8> = vec![4, 4, 4]; // lambda minus first row (for [4,4,4,4])
    for s in -1..=2 {
        let fj = shifted_exc_by_first(&board, s);
        let mat_sub = shifted_exc_matrix(&sublambda, s);
        let total_sub_s = permanent(&mat_sub);
        let mat_sub1 = shifted_exc_matrix(&sublambda, s + 1);
        let total_sub_s1 = permanent(&mat_sub1);

        // F_j^(s)(lambda) should be: (weight for j) * F^(??)(sublambda with col j removed)
        // For the full board, removing col j from [4,4,4] gives another [3,3,3]-ish board.
        // Let's just check totals.
        let mut total_from_fj = vec![0i64];
        for j in 1..=4 {
            total_from_fj = pa(&total_from_fj, &fj[j]);
        }
        let mat_full = shifted_exc_matrix(&board, s);
        let total_full = permanent(&mat_full);
        println!(
            "  s={}: Σ F_j = {}  F^(s)(lam) = {}  match={}",
            s,
            format_poly(&pt(&total_from_fj)),
            format_poly(&pt(&total_full)),
            pt(&total_from_fj) == pt(&total_full)
        );
    }
}
