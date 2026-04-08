//! Ordered big-block coefficient matrices and their planar-network model.
//!
//! For fixed j >= 2, write
//!   O_{n,j}(t) = sum_b M^{(j)}_{n,b} t^b
//! over ordered set partitions of [n].
//!
//! This binary:
//! 1. computes the exact coefficient matrix M^{(j)} from ordered set partitions,
//! 2. computes the same matrix from the planar path model,
//! 3. checks equality, and
//! 4. tests all minors of the finite window rows <= N, cols <= deg(O_{N,j}).

use std::cmp::min;

use num_bigint::BigInt;

fn binom(n: usize, k: usize) -> i128 {
    if k > n {
        return 0;
    }
    let k = min(k, n - k);
    let mut num: i128 = 1;
    let mut den: i128 = 1;
    for i in 0..k {
        num *= (n - i) as i128;
        den *= (i + 1) as i128;
    }
    num / den
}

fn factorials(max_n: usize) -> Vec<i128> {
    let mut fact = vec![1i128; max_n + 1];
    for n in 1..=max_n {
        fact[n] = fact[n - 1] * n as i128;
    }
    fact
}

fn compute_refined(max_n: usize, j: usize) -> Vec<Vec<Vec<i128>>> {
    // dp[n][m][b] = number of set partitions of [n] with m blocks and b big blocks.
    let mut dp = vec![Vec::<Vec<i128>>::new(); max_n + 1];
    dp[0] = vec![vec![1]];

    for n in 0..max_n {
        let mut next = vec![vec![0i128; n + 2]; n + 2];
        for k in 0..=n {
            let choose = binom(n, k);
            let prev_n = n - k;
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
        dp[n + 1] = next;
    }

    dp
}

fn coefficient_matrix(max_n: usize, j: usize) -> Vec<Vec<i128>> {
    let refined = compute_refined(max_n, j);
    let facts = factorials(max_n);

    let mut mat = vec![vec![1]];
    for n in 1..=max_n {
        let mut row = vec![0i128; n + 1];
        for m in 0..refined[n].len() {
            let weight = facts[m];
            for b in 0..refined[n][m].len() {
                row[b] += weight * refined[n][m][b];
            }
        }
        while row.last().is_some_and(|&x| x == 0) {
            row.pop();
        }
        mat.push(row);
    }
    mat
}

fn path_matrix(max_n: usize, j: usize) -> Vec<Vec<i128>> {
    let mut mat = vec![vec![0i128; max_n + 1]; max_n + 1];
    mat[0][0] = 1;

    for n in 1..=max_n {
        for b in 0..=n {
            let mut val = 0i128;
            for r in 1..=n {
                if r < j {
                    val += binom(n, r) * mat[n - r][b];
                } else if b > 0 {
                    val += binom(n, r) * mat[n - r][b - 1];
                }
            }
            mat[n][b] = val;
        }
    }

    let mut trimmed = Vec::with_capacity(max_n + 1);
    for row in mat {
        let mut row = row;
        while row.last().is_some_and(|&x| x == 0) {
            row.pop();
        }
        if row.is_empty() {
            row.push(0);
        }
        trimmed.push(row);
    }
    trimmed
}

fn bareiss_det(mat: &[Vec<i128>]) -> BigInt {
    let n = mat.len();
    if n == 0 {
        return BigInt::from(1);
    }
    let mut a: Vec<Vec<BigInt>> = mat
        .iter()
        .map(|row| row.iter().map(|&x| BigInt::from(x)).collect())
        .collect();
    let mut denom = BigInt::from(1);
    let mut sign = BigInt::from(1);

    for k in 0..(n - 1) {
        let mut pivot_row = k;
        while pivot_row < n && a[pivot_row][k] == BigInt::from(0) {
            pivot_row += 1;
        }
        if pivot_row == n {
            return BigInt::from(0);
        }
        if pivot_row != k {
            a.swap(k, pivot_row);
            sign = -sign;
        }
        let pivot = a[k][k].clone();
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                a[i][j] =
                    ((&a[i][j] * &pivot) - (&a[i][k] * &a[k][j])) / &denom;
            }
        }
        denom = pivot;
    }

    sign * a[n - 1][n - 1].clone()
}

fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut cur = Vec::new();
    combinations_rec(0, n, k, &mut cur, &mut out);
    out
}

fn combinations_rec(
    start: usize,
    n: usize,
    k: usize,
    cur: &mut Vec<usize>,
    out: &mut Vec<Vec<usize>>,
) {
    if cur.len() == k {
        out.push(cur.clone());
        return;
    }
    let remaining = k - cur.len();
    for x in start..=(n - remaining) {
        cur.push(x);
        combinations_rec(x + 1, n, k, cur, out);
        cur.pop();
    }
}

fn first_negative_minor(
    mat: &[Vec<i128>],
    max_row: usize,
    max_col: usize,
    size: usize,
) -> Option<(Vec<usize>, Vec<usize>, BigInt)> {
    let row_sets = combinations(max_row + 1, size);
    let col_sets = combinations(max_col + 1, size);

    for rows in &row_sets {
        for cols in &col_sets {
            let sub: Vec<Vec<i128>> = rows
                .iter()
                .map(|&r| cols.iter().map(|&c| mat[r][c]).collect())
                .collect();
            let det = bareiss_det(&sub);
            if det < BigInt::from(0) {
                return Some((rows.clone(), cols.clone(), det));
            }
        }
    }
    None
}

fn rectangularize(mat: &[Vec<i128>], width: usize) -> Vec<Vec<i128>> {
    mat.iter()
        .map(|row| {
            let mut out = row.clone();
            out.resize(width, 0);
            out
        })
        .collect()
}

fn format_poly(coeffs: &[i128]) -> String {
    let mut terms = Vec::new();
    for (i, &c) in coeffs.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let term = match (c, i) {
            (_, 0) => format!("{c}"),
            (1, 1) => "t".to_string(),
            (_, 1) => format!("{c}t"),
            (1, e) => format!("t^{e}"),
            (_, e) => format!("{c}t^{e}"),
        };
        terms.push(term);
    }
    if terms.is_empty() {
        "0".to_string()
    } else {
        terms.join(" + ")
    }
}

fn same_matrix(a: &[Vec<i128>], b: &[Vec<i128>]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    a.iter().zip(b.iter()).all(|(ra, rb)| ra == rb)
}

fn main() {
    let max_n = std::env::args()
        .nth(1)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(18);
    let max_j = std::env::args()
        .nth(2)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(6);

    println!("=== Ordered big-block coefficient matrices ===\n");
    println!("For fixed j, write O_(n,j)(t) = sum_b M^(j)_(n,b) t^b.");
    println!("We compare M^(j) with the planar path model and test minors.\n");

    for j in 2..=max_j {
        let coeff = coefficient_matrix(max_n, j);
        let path = path_matrix(max_n, j);
        let ok = same_matrix(&coeff, &path);
        let max_col = coeff.iter().map(|row| row.len()).max().unwrap_or(1) - 1;
        let rect = rectangularize(&coeff, max_col + 1);

        println!("--- j = {j} ---");
        println!("Matrix equals planar path model: {}", if ok { "yes" } else { "NO" });
        println!("Window tested: rows <= {max_n}, cols <= {max_col}");
        println!("First rows:");
        for (n, row) in coeff.iter().enumerate().take(min(max_n + 1, 9)) {
            println!("  n={}: {}", n, format_poly(row));
        }

        let max_size = min(max_n + 1, max_col + 1);
        for size in 2..=max_size {
            match first_negative_minor(&rect, max_n, max_col, size) {
                Some((rows, cols, det)) => {
                    println!("{size}x{size} minors: FAIL");
                    println!("  first negative minor rows={rows:?} cols={cols:?} det={det}");
                    break;
                }
                None => println!("{size}x{size} minors: all nonnegative"),
            }
        }
        println!();
    }
}
