//! Explore total nonnegativity of the coefficient matrix for the ordered
//! big-block family with j = 3.
//!
//! We compute O_{n,3}(t) for n <= N and form the coefficient matrix
//!   M_{n,b} = [t^b] O_{n,3}(t).
//! Since dividing row n by n! only rescales rows by a positive factor, the sign
//! pattern of minors is unchanged. Thus this integer matrix is the natural test
//! object for the "matrix pencil" extension suggested by the Branden--Leite proof.
//!
//! We also test the more naive column-production operator coming from
//! multiplication by
//!   f(z) = (e^z - 1 - z - z^2/2) / (1 - z - z^2/2).
//! Since the b-th column has EGF g(z) f(z)^b, the original coefficient matrix is
//! the orbit of column 0 under this operator. If that operator were TN, it would
//! give a simple explanation for the matrix-TN phenomenon. In fact it is not.

use std::cmp::min;

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

fn ordered_big_block_polys(max_n: usize, j: usize) -> Vec<Vec<i128>> {
    let refined = compute_refined(max_n, j);
    let facts = factorials(max_n);

    let mut polys = vec![vec![1]];
    for n in 1..=max_n {
        let mut poly = vec![0i128; n + 1];
        for m in 0..refined[n].len() {
            let weight = facts[m];
            for b in 0..refined[n][m].len() {
                poly[b] += weight * refined[n][m][b];
            }
        }
        while poly.last().is_some_and(|&x| x == 0) {
            poly.pop();
        }
        polys.push(poly);
    }
    polys
}

fn build_matrix(max_n: usize) -> Vec<Vec<i128>> {
    let polys = ordered_big_block_polys(max_n, 3);
    let max_deg = polys.iter().map(Vec::len).max().unwrap_or(0);
    let mut mat = vec![vec![0i128; max_deg]; max_n + 1];
    for (n, poly) in polys.iter().enumerate() {
        for (b, &coeff) in poly.iter().enumerate() {
            mat[n][b] = coeff;
        }
    }
    mat
}

fn f_coeffs(max_n: usize) -> Vec<i128> {
    let mut coeffs = vec![0i128; max_n + 1];
    for n in 0..=max_n {
        let mut val = if n >= 3 { 1 } else { 0 };
        if n >= 1 {
            val += n as i128 * coeffs[n - 1];
        }
        if n >= 2 {
            val += binom(n, 2) * coeffs[n - 2];
        }
        coeffs[n] = val;
    }
    coeffs
}

fn build_column_operator(max_n: usize) -> Vec<Vec<i128>> {
    let coeffs = f_coeffs(max_n);
    let mut mat = vec![vec![0i128; max_n + 1]; max_n + 1];
    for n in 0..=max_n {
        for k in 0..=n {
            mat[n][k] = binom(n, k) * coeffs[n - k];
        }
    }
    mat
}

fn bareiss_det(mat: &[Vec<i128>]) -> i128 {
    let n = mat.len();
    if n == 0 {
        return 1;
    }
    let mut a = mat.to_vec();
    let mut denom = 1i128;
    let mut sign = 1i128;

    for k in 0..(n - 1) {
        let mut pivot_row = k;
        while pivot_row < n && a[pivot_row][k] == 0 {
            pivot_row += 1;
        }
        if pivot_row == n {
            return 0;
        }
        if pivot_row != k {
            a.swap(k, pivot_row);
            sign = -sign;
        }
        let pivot = a[k][k];
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                a[i][j] = (a[i][j] * pivot - a[i][k] * a[k][j]) / denom;
            }
        }
        denom = pivot;
    }

    sign * a[n - 1][n - 1]
}

fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut cur = Vec::new();
    comb_rec(0, n, k, &mut cur, &mut out);
    out
}

fn comb_rec(start: usize, n: usize, k: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
    if cur.len() == k {
        out.push(cur.clone());
        return;
    }
    let remaining = k - cur.len();
    for x in start..=(n - remaining) {
        cur.push(x);
        comb_rec(x + 1, n, k, cur, out);
        cur.pop();
    }
}

fn first_negative_minor(
    mat: &[Vec<i128>],
    max_row: usize,
    max_col: usize,
    size: usize,
) -> Option<(Vec<usize>, Vec<usize>, i128)> {
    let row_sets = combinations(max_row + 1, size);
    let col_sets = combinations(max_col + 1, size);

    for rows in &row_sets {
        for cols in &col_sets {
            let sub: Vec<Vec<i128>> = rows
                .iter()
                .map(|&r| cols.iter().map(|&c| mat[r][c]).collect())
                .collect();
            let det = bareiss_det(&sub);
            if det < 0 {
                return Some((rows.clone(), cols.clone(), det));
            }
        }
    }
    None
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

fn main() {
    let max_n = std::env::args()
        .nth(1)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(12);
    let max_minor = std::env::args()
        .nth(2)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(4);

    let mat = build_matrix(max_n);
    let max_col = mat.iter().map(Vec::len).max().unwrap_or(1) - 1;
    let prod = build_column_operator(max_n);

    println!("=== Ordered j=3 coefficient matrix ===\n");
    println!("M_(n,b) = [t^b] O_(n,3)(t)");
    println!("Testing minors of the integer coefficient matrix.\n");

    println!("First rows:");
    for (n, row) in mat.iter().enumerate().take(min(max_n + 1, 9)) {
        println!("  n={}: {}", n, format_poly(row));
    }
    println!();

    for size in 2..=max_minor {
        let row_bound = max_n.min(match size {
            2 => 12,
            3 => 10,
            4 => 8,
            _ => 7,
        });
        let col_bound = max_col.min(row_bound);
        match first_negative_minor(&mat, row_bound, col_bound, size) {
            Some((rows, cols, det)) => {
                println!(
                    "{}x{} minors up to rows<= {}, cols<= {}: FAIL",
                    size, size, row_bound, col_bound
                );
                println!("  first negative minor rows={rows:?} cols={cols:?} det={det}");
            }
            None => {
                println!(
                    "{}x{} minors up to rows<= {}, cols<= {}: all nonnegative",
                    size, size, row_bound, col_bound
                );
            }
        }
    }

    println!();
    println!("=== Column-production operator ===\n");
    println!("The original columns satisfy C_b = T_f C_(b-1),");
    println!("where T_f is multiplication by f(z) on EGF coefficients.");
    println!("There is no standard row-production matrix: row_0 = row_1 = [1,0,...],");
    println!("but row_2 = [3,0,...], so a fixed P with row_(n+1)=row_n P cannot exist.\n");

    let coeffs = f_coeffs(max_n);
    print!("f(z) coefficients [z^n/n!]:");
    for (n, coeff) in coeffs.iter().enumerate().take(min(max_n + 1, 10)) {
        print!("  f_{n}={coeff}");
    }
    println!("\n");

    for size in 2..=min(max_minor, 4) {
        let row_bound = max_n.min(match size {
            2 => 12,
            3 => 10,
            _ => 8,
        });
        let col_bound = row_bound.min(6);
        match first_negative_minor(&prod, row_bound, col_bound, size) {
            Some((rows, cols, det)) => {
                println!(
                    "T_f {}x{} minors up to rows<= {}, cols<= {}: FAIL",
                    size, size, row_bound, col_bound
                );
                println!("  first negative minor rows={rows:?} cols={cols:?} det={det}");
            }
            None => {
                println!(
                    "T_f {}x{} minors up to rows<= {}, cols<= {}: all nonnegative",
                    size, size, row_bound, col_bound
                );
            }
        }
    }
}
