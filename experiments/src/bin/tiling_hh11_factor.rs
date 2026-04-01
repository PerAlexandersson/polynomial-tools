//! Check if M^(2)_{n,j} = sum_r C(hr,j)*C(j, n-2j-2r) is TP.
//! If so, then M = (I+S)*M^(2) proves TP of M since (I+S) is TNN.
//!
//! Also check M^(k) for k = 2,3,...,j+1 (peeling off more (1+x) factors).

use polynomial_tools::check_total_positivity;

fn binom(n: i64, k: i64) -> i64 {
    if k < 0 || k > n || n < 0 {
        return 0;
    }
    let k = k.min(n - k) as usize;
    let n = n as usize;
    let mut r = 1i64;
    for i in 0..k {
        r = r * (n - i) as i64 / (i + 1) as i64;
    }
    r
}

/// M^(p)_{n,j} = sum_r C(hr,j) * C(j+1-p, n-2j-2r)
fn m_peeled(h: i64, n: i64, j: i64, p: i64) -> i64 {
    let upper = j + 1 - p; // upper index of second binomial
    let mut total = 0i64;
    for r in 0..=50 {
        let s = n - 2 * j - 2 * r;
        if s < 0 {
            break;
        }
        total += binom(h * r, j) * binom(upper, s);
    }
    total
}

fn main() {
    for h in 2..=4 {
        println!("=== h = {} ===\n", h);

        let max_n = 25;
        let max_j = 8;

        for p in 0..=3 {
            let label = if p == 0 {
                "M (original)".to_string()
            } else {
                format!("M^({})", p)
            };

            // Build matrix
            let mut mat = Vec::new();
            for n in 0..=max_n {
                let mut row = vec![0i64; (max_j + 1) as usize];
                for j in 0..=max_j {
                    row[j as usize] = m_peeled(h as i64, n as i64, j as i64, p);
                }
                mat.push(row);
            }

            // Check TP
            let result = check_total_positivity(&mat, 4, false);
            match &result {
                Ok(()) => println!("  {}: TP holds (minors up to 4x4)", label),
                Err(msg) => println!("  {}: FAILS — {}", label, msg),
            }
        }

        // Check: does (I+S)*M^(2) = M?
        println!("\n  Verifying M = (I+S)*M^(2)...");
        let mut ok = true;
        for n in 0..=max_n {
            for j in 0..=max_j {
                let m_orig = m_peeled(h as i64, n as i64, j as i64, 0);
                let m2_n = m_peeled(h as i64, n as i64, j as i64, 1);
                let m2_n1 = if n > 0 {
                    m_peeled(h as i64, (n - 1) as i64, j as i64, 1)
                } else {
                    0
                };
                if m_orig != m2_n + m2_n1 {
                    println!(
                        "  MISMATCH at n={}, j={}: {} != {} + {}",
                        n, j, m_orig, m2_n, m2_n1
                    );
                    ok = false;
                }
            }
        }
        if ok {
            println!("  M = (I+S)*M^(2) verified. ✓");
        }

        // Check: does M^(2) = (I+S)*M^(3)?
        println!("  Verifying M^(2) = (I+S)*M^(3)...");
        let mut ok2 = true;
        for n in 0..=max_n {
            for j in 0..=max_j {
                let m2 = m_peeled(h as i64, n as i64, j as i64, 1);
                let m3_n = m_peeled(h as i64, n as i64, j as i64, 2);
                let m3_n1 = if n > 0 {
                    m_peeled(h as i64, (n - 1) as i64, j as i64, 2)
                } else {
                    0
                };
                if m2 != m3_n + m3_n1 {
                    println!("  MISMATCH at n={}, j={}", n, j);
                    ok2 = false;
                }
            }
        }
        if ok2 {
            println!("  M^(2) = (I+S)*M^(3) verified. ✓");
        }

        // Print M^(base) = M^(j+1) for a few entries
        println!("\n  Base matrix M^(j+1)_{{n,j}} = C(h(n-2j)/2, j):");
        println!(
            "  {:>4} | {:>5} {:>5} {:>5} {:>5} {:>5}",
            "n", "j=0", "j=1", "j=2", "j=3", "j=4"
        );
        for n in 0..=16i64 {
            let row: Vec<String> = (0..=4)
                .map(|j: i64| {
                    let m = n - 2 * j;
                    if m >= 0 && m % 2 == 0 {
                        format!("{:>5}", binom(h as i64 * m / 2, j))
                    } else {
                        format!("{:>5}", 0)
                    }
                })
                .collect();
            println!("  {:>4} | {}", n, row.join(" "));
        }
        println!();
    }
}
