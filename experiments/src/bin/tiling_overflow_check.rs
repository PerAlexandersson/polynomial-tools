//! Check if i64 overflow in derivative computation causes false negatives
//! in real-rootedness checks for tiling polynomials.
//!
//! Usage: cargo run --release --bin tiling_overflow_check

use polynomial_tools::{is_real_rooted, is_real_rooted_sturm, format_poly};

/// Compute P_n(t) = sum_{j=0}^{floor(n/3)} C(n-2j, j) * d^j * t^j
/// via the recursion P_n = P_{n-1} + d*t * P_{n-3}.
fn tiling_poly_rect3(n: usize, d: i64) -> Vec<i64> {
    if n == 0 { return vec![1]; }
    if n == 1 { return vec![1]; }
    if n == 2 { return vec![1]; }

    let mut pp = vec![1i64]; // P_{n-3}
    let mut pq = vec![1i64]; // P_{n-2}
    let mut pr = vec![1i64]; // P_{n-1}

    for _ in 3..=n {
        // P_new = P_{n-1} + d*t * P_{n-3}
        let max_deg = pr.len().max(pp.len() + 1);
        let mut pnew = vec![0i64; max_deg];
        for (i, &c) in pr.iter().enumerate() {
            pnew[i] += c;
        }
        for (i, &c) in pp.iter().enumerate() {
            pnew[i + 1] += d * c;
        }
        // Trim
        while pnew.len() > 1 && *pnew.last().unwrap() == 0 {
            pnew.pop();
        }
        pp = pq;
        pq = pr;
        pr = pnew;
    }
    pr
}

fn main() {
    println!("Checking derivative overflow for claw-free rectangle tilings\n");

    // mu=(k,k,k), d=k: P_n(t) = sum_j C(n-2j,j) * k^j * t^j
    // Fails at n=69 for d=3, n=62 for d=4, n=58 for d=5

    let test_cases = vec![
        (3, 69),
        (3, 68),
        (4, 62),
        (4, 61),
        (5, 58),
        (5, 57),
    ];

    for &(d, n) in &test_cases {
        let p = tiling_poly_rect3(n, d);
        let deg = p.len() - 1;

        // Check for derivative overflow
        let mut overflow = false;
        for k in 0..deg {
            let coeff = p[k + 1];
            let mult = (k as i64) + 1;
            if coeff.checked_mul(mult).is_none() {
                println!("d={}, n={}: OVERFLOW at derivative term k={}: {} * {} overflows i64",
                         d, n, k, mult, coeff);
                overflow = true;
                break;
            }
        }
        if !overflow {
            println!("d={}, n={}: No i64 overflow in derivative (deg={})", d, n, deg);
        }

        // Check real-rootedness with both methods
        let rr_bezout = is_real_rooted(&p);
        let rr_sturm = is_real_rooted_sturm(&p);
        println!("  Bezout: {}, Sturm: {}", rr_bezout, rr_sturm);

        if rr_bezout != rr_sturm {
            println!("  *** MISMATCH! ***");
            if deg <= 15 {
                println!("  P(t) = {}", format_poly(&p));
            }
        }
        println!();
    }

    // Also check a few non-claw-free cases where failure is expected
    println!("\n--- Non-claw-free cases (expected failures) ---\n");
    // mu=(3,1,1), d=3: fails at n=19 per the paper
    // This uses a different recursion, so let's compute via transfer matrix
    // For now just check the rectangle cases with d > k
    let extra_cases = vec![
        (3, 56, 5), // [2,2,2], d=3
        (3, 14, 7), // [2,2,2], d=5
    ];
    for &(d, n, _d_actual) in &extra_cases {
        let p = tiling_poly_rect3(n, d);
        let rr_b = is_real_rooted(&p);
        let rr_s = is_real_rooted_sturm(&p);
        println!("d_eff={}, n={}: Bezout={}, Sturm={}, match={}",
                 d, n, rr_b, rr_s, rr_b == rr_s);
    }
}
