//! Brenti matrix conjecture: if B(x,t) ∈ ℕ[x,t] with B(0,0)=0,
//! is the coefficient matrix of F(x,t) = 1/(1-B(x,t)) a Brenti matrix?
//!
//! We verify TNN (via Neville elimination) and real-rootedness of each row
//! for concrete examples.
//!
//! The digraph D_B has step types given by the monomials of B:
//!   each monomial c·x^a·t^b gives an edge (m,k) → (m+a, k+b) with weight c.
//! Path counts in D_B equal [x^n t^j] 1/(1-B).
//!
//! Usage: cargo run --release --bin brenti_geometric

use num_bigint::BigInt;
use polynomial_tools::{check_tnn_neville_bigint, is_real_rooted};
use std::time::Instant;

/// A monomial c * x^a * t^b in B(x,t).
#[derive(Clone, Debug)]
struct Step {
    a: usize, // x-exponent (displacement in first coordinate)
    b: usize, // t-exponent (displacement in second coordinate)
    c: i128,  // weight
}

/// Compute the coefficient matrix [x^n t^j] F = 1/(1-B) up to row max_n.
///
/// Uses the recurrence f[n][j] = sum over steps (a,b,c): c * f[n-a][j-b].
fn compute_matrix(steps: &[Step], max_n: usize) -> Vec<Vec<i128>> {
    // Determine max possible j: each step adds at most max_b per step,
    // and we can take at most max_n steps (if min step length is 1).
    let max_b = steps.iter().map(|s| s.b).max().unwrap_or(0);
    let max_j = max_n * max_b;

    let mut f = vec![vec![0i128; max_j + 1]; max_n + 1];
    f[0][0] = 1;

    for n in 1..=max_n {
        for j in 0..=max_j {
            let mut val = 0i128;
            for s in steps {
                if n >= s.a && j >= s.b {
                    val += s.c * f[n - s.a][j - s.b];
                }
            }
            f[n][j] = val;
        }
    }

    // Trim trailing zero columns
    let max_nonzero_j = f
        .iter()
        .map(|row| row.iter().rposition(|&x| x != 0).unwrap_or(0))
        .max()
        .unwrap_or(0);

    f.iter().map(|row| row[..=max_nonzero_j].to_vec()).collect()
}

fn format_poly_t(coeffs: &[i128]) -> String {
    let terms: Vec<String> = coeffs
        .iter()
        .enumerate()
        .filter(|(_, &c)| c != 0)
        .map(|(j, &c)| match j {
            0 => format!("{c}"),
            1 if c == 1 => "t".into(),
            1 => format!("{c}t"),
            _ if c == 1 => format!("t^{j}"),
            _ => format!("{c}t^{j}"),
        })
        .collect();
    if terms.is_empty() {
        "0".into()
    } else {
        terms.join(" + ")
    }
}

fn format_b(steps: &[Step]) -> String {
    let terms: Vec<String> = steps
        .iter()
        .map(|s| {
            let coeff = if s.c == 1 {
                String::new()
            } else {
                format!("{}", s.c)
            };
            let xpart = match s.a {
                0 => String::new(),
                1 => "x".into(),
                a => format!("x^{a}"),
            };
            let tpart = match s.b {
                0 => String::new(),
                1 => "t".into(),
                b => format!("t^{b}"),
            };
            format!("{coeff}{xpart}{tpart}")
        })
        .collect();
    terms.join(" + ")
}

fn analyze(label: &str, steps: &[Step], max_n: usize) {
    println!("═══════════════════════════════════════════════════════");
    println!("B(x,t) = {}", format_b(steps));
    if !label.is_empty() {
        println!("  [{label}]");
    }
    println!();

    let t0 = Instant::now();
    let mat = compute_matrix(steps, max_n);
    let nrows = mat.len();
    let ncols = mat[0].len();

    // Print matrix
    println!("Coefficient matrix of F = 1/(1-B), size {nrows}×{ncols}:");
    // Determine column widths
    let widths: Vec<usize> = (0..ncols)
        .map(|j| {
            mat.iter()
                .map(|row| format!("{}", row[j]).len())
                .max()
                .unwrap_or(1)
                .max(3)
        })
        .collect();

    print!("{:>4} |", "n\\j");
    for (j, &w) in widths.iter().enumerate() {
        print!(" {:>w$}", j, w = w);
    }
    println!();
    print!("-----|");
    for &w in &widths {
        print!("{}", "-".repeat(w + 1));
    }
    println!();
    for (n, row) in mat.iter().enumerate() {
        print!("{:>4} |", n);
        for (j, &w) in widths.iter().enumerate() {
            print!(" {:>w$}", row[j], w = w);
        }
        println!();
    }
    println!();

    // Check real-rootedness of each row polynomial
    println!("Real-rootedness of row polynomials:");
    let mut all_rr = true;
    for (n, row) in mat.iter().enumerate() {
        let coeffs_i64: Vec<i64> = row.iter().map(|&c| c as i64).collect();
        let deg = coeffs_i64.iter().rposition(|&c| c != 0).unwrap_or(0);
        if deg <= 1 {
            println!(
                "  n={n:>3}: {} — trivially real-rooted (deg ≤ 1)",
                format_poly_t(row)
            );
            continue;
        }
        let rr = is_real_rooted(&coeffs_i64);
        if !rr {
            all_rr = false;
        }
        println!(
            "  n={n:>3}: deg {deg} — {}",
            if rr {
                "REAL-ROOTED ✓"
            } else {
                "NOT real-rooted ✗"
            }
        );
    }
    println!(
        "  => All rows real-rooted: {}",
        if all_rr { "YES ✓" } else { "NO ✗" }
    );
    println!();

    // Check TNN via Neville elimination
    let bigmat: Vec<Vec<BigInt>> = mat
        .iter()
        .map(|row| row.iter().map(|&c| BigInt::from(c)).collect())
        .collect();
    let tnn_result = check_tnn_neville_bigint(&bigmat);
    match &tnn_result {
        Ok(()) => println!("TNN check (Neville): PASS ✓"),
        Err(msg) => println!("TNN check (Neville): FAIL ✗ — {msg}"),
    }

    let elapsed = t0.elapsed();
    println!("  (computed in {elapsed:.2?})");
    println!();
}

fn main() {
    let max_n = 30;

    // Example 1: B = x + 3tx^2 + t^2 x^5
    analyze(
        "concrete example",
        &[
            Step { a: 1, b: 0, c: 1 },
            Step { a: 2, b: 1, c: 3 },
            Step { a: 5, b: 2, c: 1 },
        ],
        max_n,
    );

    // Example 2: B = x + xt (Pascal / binomial)
    analyze(
        "Pascal: rows = binomial coefficients",
        &[Step { a: 1, b: 0, c: 1 }, Step { a: 1, b: 1, c: 1 }],
        max_n,
    );

    // Example 3: B = xt (diagonal)
    analyze(
        "diagonal: identity matrix",
        &[Step { a: 1, b: 1, c: 1 }],
        20,
    );

    // Example 4: B = x + 2xt + xt^2 (trinomial)
    analyze(
        "trinomial steps",
        &[
            Step { a: 1, b: 0, c: 1 },
            Step { a: 1, b: 1, c: 2 },
            Step { a: 1, b: 2, c: 1 },
        ],
        max_n,
    );

    // Example 5: B = x + tx^2 (Fibonacci-like)
    analyze(
        "Fibonacci-like: Catalan numbers appear",
        &[Step { a: 1, b: 0, c: 1 }, Step { a: 2, b: 1, c: 1 }],
        max_n,
    );

    // Example 6: B = x + t + xt (includes vertical step)
    analyze(
        "with vertical step t",
        &[
            Step { a: 1, b: 0, c: 1 },
            Step { a: 0, b: 1, c: 1 },
            Step { a: 1, b: 1, c: 1 },
        ],
        max_n,
    );

    // Example 7: steps with different slopes
    analyze(
        "mixed slopes",
        &[
            Step { a: 1, b: 0, c: 1 },
            Step { a: 1, b: 1, c: 1 },
            Step { a: 2, b: 1, c: 1 },
            Step { a: 1, b: 2, c: 1 },
        ],
        max_n,
    );

    // Example 8: larger weights
    analyze(
        "large weight step",
        &[
            Step { a: 1, b: 0, c: 2 },
            Step { a: 1, b: 1, c: 5 },
            Step { a: 3, b: 2, c: 3 },
        ],
        max_n,
    );
}
