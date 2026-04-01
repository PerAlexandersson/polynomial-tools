//! Test the interlacing-sequence conjecture for the (A_m,...,A_1,W_1,...,W_m) ordering,
//! and print the peak recursion matrix M mapping (A,W) at level λ to (A^+,W^+) at level λ+.
//!
//! Recursion:
//!   A_k^{λ+} = Σ_{j<k} A_j^λ + Σ_{j≥k} W_j^λ
//!   W_k^{λ+} = t·Σ_{j<k} A_j^λ + Σ_{j≥k} W_j^λ
//!
//! The input vector ordering is  (A_m, ..., A_1, W_m, ..., W_1)   [length 2m]
//! The output vector ordering is (A_{m'}, ..., A_1, W_1, ..., W_{m'}) [length 2m']
//! where m' = m+1 for a board λ+ = (λ_1+1, ..., λ_n+1, 1).
//!
//! Tests:
//! - "Reversed AA": A_k ≪ A_j for k > j (i.e. higher-index interlaces lower-index)
//! - Full interlacing sequence: (A_m, A_{m-1}, ..., A_1, W_1, W_2, ..., W_m)
//!   is interlacing, meaning f_i ≪ f_j for all i < j in this list.

use combpoly::order::bruhat_lower_ideal;
use experiments::peak_utils::{
    board_to_perm, gen_boards, is_312_avoiding, pa, peak_count, pmt, pt, pz,
};
use num::BigInt;
use num_rational::Ratio;
use polynomial_tools::polynomial::Polynomial;

type Poly = Polynomial<Ratio<BigInt>>;
type BR = Ratio<BigInt>;
fn br(n: i64) -> BR {
    BR::from_integer(BigInt::from(n))
}
fn to_poly(c: &[i64]) -> Poly {
    Polynomial::new(
        c.iter()
            .map(|&x| BR::from_integer(BigInt::from(x)))
            .collect(),
    )
}

fn poly_long_div(f: &Poly, g: &Poly) -> Poly {
    let (q, _) = f.div_rem(g);
    q
}

fn polynomial_gcd(f: &Poly, g: &Poly) -> Poly {
    f.gcd(g)
}

fn lagrange(pts: &[BR], vals: &[BR]) -> Vec<BR> {
    let p = Polynomial::lagrange_interpolation(pts, vals);
    let d = p.degree().unwrap_or(0);
    (0..=d).map(|i| p.coeff(i)).collect()
}

fn exact_interlaces(fc: &[i64], gc: &[i64]) -> bool {
    if let Some(result) = polynomial_tools::check_weak_interlacing(fc, gc) {
        return result;
    }
    polynomial_tools::check_interlacing_sturm(fc, gc).unwrap_or(false)
}

fn pf(p: &[i64]) -> String {
    let p = pt(p);
    if pz(&p) {
        return "0".into();
    }
    let mut t = vec![];
    for (i, &c) in p.iter().enumerate() {
        if c == 0 {
            continue;
        }
        match (c, i) {
            (c, 0) => t.push(format!("{}", c)),
            (1, 1) => t.push("t".into()),
            (c, 1) => t.push(format!("{}t", c)),
            (1, e) => t.push(format!("t^{}", e)),
            (c, e) => t.push(format!("{}t^{}", c, e)),
        }
    }
    t.join(" + ")
}
/// Compute the A and W polynomials for a given board.
fn compute_aw(board: &[u8]) -> (Vec<Vec<i64>>, Vec<Vec<i64>>) {
    let n = board.len();
    let perm = board_to_perm(board);
    let m = (board[0] as usize).min(n);
    let ideal = bruhat_lower_ideal(&perm);
    let mut d = vec![vec![0i64]; m + 1];
    let mut u = vec![vec![0i64]; m + 1];
    for pi in &ideal {
        if pi.len() < 2 {
            let k = pi[0] as usize;
            if k <= m {
                while u[k].len() < 1 {
                    u[k].push(0);
                }
                u[k][0] += 1;
            }
            continue;
        }
        let k = pi[0] as usize;
        if k > m {
            continue;
        }
        let pk = peak_count(pi);
        let poly = if pi[0] > pi[1] { &mut d[k] } else { &mut u[k] };
        while poly.len() <= pk {
            poly.push(0);
        }
        poly[pk] += 1;
    }
    let mut a = vec![vec![0i64]; m + 1];
    let mut w = vec![vec![0i64]; m + 1];
    for k in 1..=m {
        a[k] = pa(&d[k], &u[k]);
        w[k] = pa(&pmt(&d[k]), &u[k]);
    }
    (a, w)
}

/// Print the recursion matrix for a rectangular board (c^r means all rows = c, r rows).
/// The recursion maps (A_m, ..., A_1, W_m, ..., W_1) -> (A_{m+1}, ..., A_1, W_1, ..., W_{m+1}).
///
/// A_k^+ = Σ_{j<k} A_j + Σ_{j≥k} W_j     (for k = 1, ..., m+1)
/// W_k^+ = t * Σ_{j<k} A_j + Σ_{j≥k} W_j  (for k = 1, ..., m+1)
///
/// The matrix entries are in {0, 1, t} indicating the coefficient applied to each input.
fn print_recursion_matrix(board_name: &str, m: usize) {
    let m_prime = m + 1;
    // Input ordering:  A_m, A_{m-1}, ..., A_1, W_m, W_{m-1}, ..., W_1  (length 2m)
    // Output ordering: A_{m'}, A_{m'-1}, ..., A_1, W_1, W_2, ..., W_{m'} (length 2m')

    println!(
        "\n--- Recursion matrix for {} (m={}, m'={}) ---",
        board_name, m, m_prime
    );
    println!("Input  (cols): A_{}, ..., A_1, W_{}, ..., W_1", m, m);
    println!(
        "Output (rows): A_{}, ..., A_1, W_1, ..., W_{}",
        m_prime, m_prime
    );
    println!();

    // Column headers
    let mut col_labels = vec![];
    for j in (1..=m).rev() {
        col_labels.push(format!("A{}", j));
    }
    for j in (1..=m).rev() {
        col_labels.push(format!("W{}", j));
    }
    let ncols = 2 * m;

    // Row labels
    let mut row_labels = vec![];
    for k in (1..=m_prime).rev() {
        row_labels.push(format!("A{}+", k));
    }
    for k in 1..=m_prime {
        row_labels.push(format!("W{}+", k));
    }
    let nrows = 2 * m_prime;

    // Build matrix entries as strings ("0", "1", "t")
    let mut matrix = vec![vec!["0"; ncols]; nrows];

    for k in 1..=m_prime {
        // Row index for A_k^+ : position (m_prime - k) in the first block
        let row_a = m_prime - k;
        // Row index for W_k^+ : position (k - 1) in the second block, offset by m_prime
        let row_w = m_prime + (k - 1);

        // A_k^+ = Σ_{j<k} A_j + Σ_{j≥k} W_j
        // W_k^+ = t * Σ_{j<k} A_j + Σ_{j≥k} W_j

        // Contribution from A_j (j ranges over 1..=m)
        for j in 1..=m {
            let col_a = m - j; // column for A_j in input
            if j < k {
                // A_j contributes with coefficient 1 to A_k^+
                matrix[row_a][col_a] = "1";
                // A_j contributes with coefficient t to W_k^+
                matrix[row_w][col_a] = "t";
            }
            // j >= k: A_j does NOT contribute to A_k^+ or W_k^+
        }

        // Contribution from W_j (j ranges over 1..=m)
        for j in 1..=m {
            let col_w = m + (m - j); // column for W_j in input
            if j >= k {
                // W_j contributes with coefficient 1 to A_k^+
                matrix[row_a][col_w] = "1";
                // W_j contributes with coefficient 1 to W_k^+
                matrix[row_w][col_w] = "1";
            }
        }
    }

    // Print header
    print!("{:>6}", "");
    for l in &col_labels {
        print!("{:>4}", l);
    }
    println!();

    for (i, row) in matrix.iter().enumerate() {
        print!("{:>6}", row_labels[i]);
        for entry in row {
            print!("{:>4}", entry);
        }
        println!();
    }
}

fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);
    println!(
        "=== Interlacing sequence conjecture + recursion matrix, n <= {} ===",
        max_n
    );
    println!("Convention: f << g means f interlaces g (roots of f separate roots of g).\n");

    // ── Part 1: Print recursion matrices for example boards ──
    println!("========== RECURSION MATRICES ==========");
    // (2,2): m = 2
    print_recursion_matrix("lambda = (2,2)", 2);
    // (3,3,3): m = 3
    print_recursion_matrix("lambda = (3,3,3)", 3);
    // (4,4,4,4): m = 4
    print_recursion_matrix("lambda = (4,4,4,4)", 4);

    // ── Part 2: Test interlacing conjectures ──
    println!("\n\n========== INTERLACING TESTS ==========\n");

    let mut total = 0u64;
    let mut valid = 0u64;

    // Counters: [tested, failed]
    let mut aa_rev = [0u64; 2]; // A_k << A_j for k > j ("reversed AA")
    let mut aw_all = [0u64; 2]; // A_j << W_l for all j,l (existing condition)
    let mut ww_ord = [0u64; 2]; // W_j << W_l for j <= l (existing condition)
    let mut full_seq = [0u64; 2]; // Full interlacing sequence check

    for n in 1..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            total += 1;
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            valid += 1;
            let (a, w) = compute_aw(board);
            let m = a.len() - 1; // a and w are indexed 1..=m

            // Test AA reversed: A_k << A_j for k > j
            for k in 1..=m {
                for j in 1..k {
                    if !pz(&a[k]) && !pz(&a[j]) {
                        aa_rev[0] += 1;
                        if !exact_interlaces(&a[k], &a[j]) {
                            aa_rev[1] += 1;
                            if aa_rev[1] <= 5 {
                                println!(
                                    "FAIL AA_rev({},{}) board={:?}: {} << {}",
                                    k,
                                    j,
                                    board,
                                    pf(&a[k]),
                                    pf(&a[j])
                                );
                            }
                        }
                    }
                }
            }

            // Test AW all: A_j << W_l for all j, l
            for j in 1..=m {
                for l in 1..=m {
                    if !pz(&a[j]) && !pz(&w[l]) {
                        aw_all[0] += 1;
                        if !exact_interlaces(&a[j], &w[l]) {
                            aw_all[1] += 1;
                            if aw_all[1] <= 5 {
                                println!(
                                    "FAIL AW({},{}) board={:?}: {} << {}",
                                    j,
                                    l,
                                    board,
                                    pf(&a[j]),
                                    pf(&w[l])
                                );
                            }
                        }
                    }
                }
            }

            // Test WW ordered: W_j << W_l for j <= l
            for j in 1..=m {
                for l in j..=m {
                    if j < l && !pz(&w[j]) && !pz(&w[l]) {
                        ww_ord[0] += 1;
                        if !exact_interlaces(&w[j], &w[l]) {
                            ww_ord[1] += 1;
                            if ww_ord[1] <= 5 {
                                println!(
                                    "FAIL WW({},{}) board={:?}: {} << {}",
                                    j,
                                    l,
                                    board,
                                    pf(&w[j]),
                                    pf(&w[l])
                                );
                            }
                        }
                    }
                }
            }

            // Test FULL interlacing sequence: (A_m, A_{m-1}, ..., A_1, W_1, W_2, ..., W_m)
            // f_i << f_j for all i < j in this ordering
            {
                let mut seq: Vec<&Vec<i64>> = Vec::new();
                let mut seq_labels: Vec<String> = Vec::new();
                for k in (1..=m).rev() {
                    seq.push(&a[k]);
                    seq_labels.push(format!("A{}", k));
                }
                for k in 1..=m {
                    seq.push(&w[k]);
                    seq_labels.push(format!("W{}", k));
                }

                let mut seq_ok = true;
                for i in 0..seq.len() {
                    for j in i + 1..seq.len() {
                        if !pz(seq[i]) && !pz(seq[j]) {
                            if !exact_interlaces(seq[i], seq[j]) {
                                seq_ok = false;
                            }
                        }
                    }
                }
                full_seq[0] += 1;
                if !seq_ok {
                    full_seq[1] += 1;
                    if full_seq[1] <= 5 {
                        println!("FAIL full_seq board={:?}", board);
                        // Print which pairs fail
                        for i in 0..seq.len() {
                            for j in i + 1..seq.len() {
                                if !pz(seq[i]) && !pz(seq[j]) && !exact_interlaces(seq[i], seq[j]) {
                                    println!(
                                        "  {} << {} fails: {} << {}",
                                        seq_labels[i],
                                        seq_labels[j],
                                        pf(seq[i]),
                                        pf(seq[j])
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
        println!("n={}: {}/{} valid 312-avoiding boards", n, valid, total);
    }

    println!("\n\n========== RESULTS ==========\n");
    let show = |name: &str, c: [u64; 2]| {
        if c[0] == 0 {
            println!("  {}: (no tests)", name);
        } else {
            println!(
                "  {}: {}/{} pass {}",
                name,
                c[0] - c[1],
                c[0],
                if c[1] == 0 { "ALL PASS" } else { "FAILURES" }
            );
        }
    };
    println!("Existing conditions (sanity check):");
    show("(a) A_j << W_l (all j,l)", aw_all);
    show("(b) W_j << W_l (j < l)", ww_ord);
    println!("\nNew tests:");
    show("Reversed AA: A_k << A_j (k > j)", aa_rev);
    show("Full interlacing seq (A_m,...,A_1,W_1,...,W_m)", full_seq);

    // ── Part 3: Print example (A,W) data for small boards ──
    println!("\n\n========== EXAMPLE (A,W) DATA ==========\n");
    let example_boards: Vec<Vec<u8>> = vec![
        vec![2, 2],
        vec![3, 3, 3],
        vec![4, 4, 4, 4],
        vec![3, 3], // non-square Ferrers board
        vec![4, 4, 4],
    ];
    for board in &example_boards {
        let perm = board_to_perm(board);
        if !is_312_avoiding(&perm) {
            println!("Board {:?}: not 312-avoiding, skipping", board);
            continue;
        }
        let (a, w) = compute_aw(board);
        let m = a.len() - 1;
        println!("Board {:?} (m={}), perm={:?}", board, m, perm);
        for k in 1..=m {
            println!("  A_{} = {}    W_{} = {}", k, pf(&a[k]), k, pf(&w[k]));
        }
        // Show the interlacing sequence
        print!("  Sequence: ");
        for k in (1..=m).rev() {
            print!("A{}, ", k);
        }
        for k in 1..=m {
            if k < m {
                print!("W{}, ", k);
            } else {
                print!("W{}", k);
            }
        }
        println!();
        // Check pairwise
        let mut seq: Vec<&Vec<i64>> = Vec::new();
        let mut labels: Vec<String> = Vec::new();
        for k in (1..=m).rev() {
            seq.push(&a[k]);
            labels.push(format!("A{}", k));
        }
        for k in 1..=m {
            seq.push(&w[k]);
            labels.push(format!("W{}", k));
        }
        let mut all_ok = true;
        for i in 0..seq.len() {
            for j in i + 1..seq.len() {
                if !pz(seq[i]) && !pz(seq[j]) {
                    let ok = exact_interlaces(seq[i], seq[j]);
                    if !ok {
                        println!(
                            "    {} << {}: FAIL ({} << {})",
                            labels[i],
                            labels[j],
                            pf(seq[i]),
                            pf(seq[j])
                        );
                        all_ok = false;
                    }
                }
            }
        }
        if all_ok {
            println!("  All pairwise interlacings hold.");
        }
        println!();
    }
}
