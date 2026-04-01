/// Check Brändén's Theorem 7.8.5 condition on the transfer matrix.
///
/// A polynomial matrix G maps F_n^+ -> F_m^+ iff:
/// 1. All entries have non-negative coefficients
/// 2. For every 2x2 submatrix [[a,b],[c,d]] and lambda, mu > 0:
///    (lambda*t + mu)*b + d <<  (lambda*t + mu)*a + c
///
/// We check this for:
/// (a) The pure staircase matrix (entries 1 and t)
/// (b) The actual transfer matrix with corrections
///
use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly, is_real_rooted};

/// Check Brändén's condition for a 2x2 submatrix [[a,b],[c,d]]
/// at several (lambda, mu) values.
/// Returns true if the condition holds for all tested values.
fn check_branden_2x2(a: &[i64], b: &[i64], c: &[i64], d: &[i64]) -> bool {
    let test_pairs: Vec<(f64, f64)> = vec![
        (1.0, 1.0),
        (1.0, 2.0),
        (2.0, 1.0),
        (1.0, 0.5),
        (0.5, 1.0),
        (3.0, 1.0),
        (1.0, 3.0),
        (5.0, 1.0),
        (1.0, 5.0),
        (0.1, 1.0),
        (1.0, 0.1),
        (10.0, 1.0),
        (1.0, 10.0),
        (2.0, 3.0),
        (3.0, 2.0),
        (0.01, 1.0),
        (1.0, 0.01),
        (100.0, 1.0),
        (1.0, 100.0),
        (7.0, 3.0),
        (3.0, 7.0),
        (0.5, 0.5),
    ];

    for &(lambda, mu) in &test_pairs {
        // Compute LHS = (lambda*t + mu)*b + d
        // Compute RHS = (lambda*t + mu)*a + c
        let lhs = lin_mul_plus(lambda, mu, b, d);
        let rhs = lin_mul_plus(lambda, mu, a, c);

        // Check LHS << RHS (interleaving)
        // Convert to integers (scale by large factor)
        let scale = 1000000.0;
        let lhs_int: Vec<i64> = lhs.iter().map(|&x| (x * scale).round() as i64).collect();
        let rhs_int: Vec<i64> = rhs.iter().map(|&x| (x * scale).round() as i64).collect();

        // Trim
        let lhs_t = trim(&lhs_int);
        let rhs_t = trim(&rhs_int);

        if lhs_t.len() <= 1 && rhs_t.len() <= 1 {
            continue; // Trivial
        }

        // Check: every convex combination is real-rooted
        let combos: Vec<(i64, i64)> = vec![
            (1, 0),
            (0, 1),
            (1, 1),
            (1, 2),
            (2, 1),
            (1, 3),
            (3, 1),
            (1, 5),
            (5, 1),
            (2, 3),
            (3, 2),
        ];
        for &(alpha, beta) in &combos {
            let maxlen = lhs_t.len().max(rhs_t.len());
            let mut combo = vec![0i64; maxlen];
            for (i, &c) in lhs_t.iter().enumerate() {
                combo[i] += alpha * c;
            }
            for (i, &c) in rhs_t.iter().enumerate() {
                combo[i] += beta * c;
            }
            let combo = trim(&combo);
            if combo.len() > 2 && !is_real_rooted(&combo) {
                return false;
            }
        }
    }
    true
}

fn lin_mul_plus(lambda: f64, mu: f64, b: &[i64], d: &[i64]) -> Vec<f64> {
    // (lambda*t + mu)*b + d
    let deg_b = if b.is_empty() { 0 } else { b.len() - 1 };
    let deg_d = if d.is_empty() { 0 } else { d.len() - 1 };
    let deg_result = (deg_b + 1).max(deg_d);
    let mut result = vec![0.0f64; deg_result + 1];
    // (lambda*t + mu)*b: coefficient of t^k = lambda*b[k-1] + mu*b[k]
    for (k, &bk) in b.iter().enumerate() {
        result[k] += mu * bk as f64;
        if k + 1 < result.len() {
            result[k + 1] += lambda * bk as f64;
        } else {
            result.push(lambda * bk as f64);
        }
    }
    for (k, &dk) in d.iter().enumerate() {
        if k < result.len() {
            result[k] += dk as f64;
        }
    }
    result
}

fn trim(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
    v
}

fn main() {
    println!("=== Checking Brändén's condition on the pure staircase matrix ===\n");
    println!("Staircase: M_{{p,q}} = t if q <= p-2, else 1\n");

    // The possible 2x2 submatrices
    let submatrices: Vec<(&str, [i64; 2], [i64; 2], [i64; 2], [i64; 2])> = vec![
        ("[[t,t],[t,t]]", [0, 1], [0, 1], [0, 1], [0, 1]),
        ("[[t,1],[t,t]]", [0, 1], [1], [0, 1], [0, 1]),
        ("[[t,1],[t,1]]", [0, 1], [1], [0, 1], [1]),
        ("[[1,1],[t,t]]", [1], [1], [0, 1], [0, 1]),
        ("[[1,1],[t,1]]", [1], [1], [0, 1], [1]),
        ("[[1,1],[1,1]]", [1], [1], [1], [1]),
    ];

    for (name, a, b, c, d) in &submatrices {
        let ok = check_branden_2x2(a, b, c, d);
        println!("  {} : {}", name, if ok { "OK" } else { "FAIL" });
    }

    println!("\n=== Now checking with correction entries ===\n");
    println!("When eps1=1, the entry becomes t*t^eps2 instead of t^eps2.");
    println!("So some entries change from 1 to t, or from t to t^2.\n");

    // The corrected entries: for the eps1=1 subset, multiply by extra t
    // So the entry for (p,q) becomes:
    //   t^{eps2} * (L^(q) - C) + t^{eps2+1} * C = t^{eps2} * [L^(q) + (t-1)*C]
    // If we split the source into eps1=0 and eps1=1 parts:
    //   A_q with multiplier t^{eps2}
    //   C_q with multiplier t^{eps2+1}
    //
    // The extended matrix has columns for (A_q, C_q, D_q) per source q,
    // with entries (t^eps2, t^{eps2+1}, t^eps2).

    // Check: does the extended staircase with entries {1, t, t^2} work?
    // Submatrices involving t^2 entries:
    let extended_submatrices: Vec<(&str, Vec<i64>, Vec<i64>, Vec<i64>, Vec<i64>)> = vec![
        // From the A columns (multiplier 1 or t):
        // already checked above
        // From mixing A (mult t^eps2) and C (mult t^{eps2+1}):
        // When eps2=0: A gets 1, C gets t
        ("[[1,t],[1,t]]", vec![1], vec![0, 1], vec![1], vec![0, 1]),
        (
            "[[1,t],[t,t^2]]",
            vec![1],
            vec![0, 1],
            vec![0, 1],
            vec![0, 0, 1],
        ),
        (
            "[[t,t^2],[1,t]]",
            vec![0, 1],
            vec![0, 0, 1],
            vec![1],
            vec![0, 1],
        ),
        (
            "[[t,t^2],[t,t^2]]",
            vec![0, 1],
            vec![0, 0, 1],
            vec![0, 1],
            vec![0, 0, 1],
        ),
        // When eps2=1 for top, eps2=0 for bottom (or vice versa):
        (
            "[[t,t^2],[1,1]]",
            vec![0, 1],
            vec![0, 0, 1],
            vec![1],
            vec![1],
        ),
        (
            "[[1,1],[t,t^2]]",
            vec![1],
            vec![1],
            vec![0, 1],
            vec![0, 0, 1],
        ),
        (
            "[[t,t],[t,t^2]]",
            vec![0, 1],
            vec![0, 1],
            vec![0, 1],
            vec![0, 0, 1],
        ),
        (
            "[[t^2,t],[t^2,t]]",
            vec![0, 0, 1],
            vec![0, 1],
            vec![0, 0, 1],
            vec![0, 1],
        ),
        // Mixing D columns (same as A columns but from desc source):
        // Same as pure staircase — already checked
    ];

    for (name, a, b, c, d) in &extended_submatrices {
        let ok = check_branden_2x2(a, b, c, d);
        println!("  {} : {}", name, if ok { "OK" } else { "FAIL" });
    }

    println!("\n=== Key insight ===");
    println!("If all 2x2 submatrices pass, the transfer matrix preserves");
    println!("interlacing sequences by Branden's Thm 7.8.5.");
}
