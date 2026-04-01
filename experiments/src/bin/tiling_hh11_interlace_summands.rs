//! Check if the summands Q_{n,r}(t) = (1+t)^{hr} ⊙ g_{n-2r}(t)
//! interlace each other (for fixed n, varying r).
//! If they form an interlacing family, their sum P_n is real-rooted.

use polynomial_tools::{check_interlacing, check_weak_interlacing, is_real_rooted};

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

/// Compute Q_{n,r}(t) = sum_j C(hr,j)*C(j+1, n-2j-2r) t^j
fn q_poly(h: i64, n: i64, r: i64) -> Vec<i64> {
    let mut coeffs = Vec::new();
    let mut started = false;
    for j in 0..=n {
        let s = n - 2 * j - 2 * r;
        let c = binom(h * r, j) * binom(j + 1, s);
        if c > 0 {
            started = true;
        }
        if started {
            coeffs.push(c);
        }
    }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
        coeffs.pop();
    }
    coeffs
}

fn poly_is_trivial(p: &[i64]) -> bool {
    p.len() <= 1 || p.iter().all(|&c| c == 0)
}

fn main() {
    let max_n = 50;

    for h in 2..=5i64 {
        println!("=== h = {} ===", h);
        let mut all_interlace = true;
        let mut all_common = true;

        for n in 0..=max_n {
            // Collect all non-trivial Q_{n,r}
            let mut qs: Vec<(i64, Vec<i64>)> = Vec::new();
            for r in 0..=n / 2 {
                let q = q_poly(h, n as i64, r as i64);
                if !poly_is_trivial(&q) {
                    qs.push((r as i64, q));
                }
            }

            if qs.len() < 2 {
                continue;
            }

            // Check consecutive interlacing: Q_{n,r} interlaces Q_{n,r+1}?
            for i in 0..qs.len() - 1 {
                let (r1, ref q1) = qs[i];
                let (r2, ref q2) = qs[i + 1];
                let result = check_weak_interlacing(q1, q2);
                if result != Some(true) {
                    if all_interlace {
                        println!(
                            "  Consecutive interlacing FAILS at n={}, r={}->{}: {:?}",
                            n, r1, r2, result
                        );
                    }
                    all_interlace = false;
                }
            }

            // Check if P_n interlaces each Q_{n,r} (common interlacing)
            // P_n = sum_r Q_{n,r}
            let max_deg = qs.iter().map(|(_, q)| q.len()).max().unwrap_or(0);
            // Need to align: Q_{n,r} starts at different t-powers
            // Actually, q_poly already strips leading zeros, so need offset
            // Let me recompute with full coefficient vectors
            let mut p_n = vec![0i64; (n / 2 + 5) as usize];
            let mut q_full: Vec<(i64, Vec<i64>)> = Vec::new();
            for r in 0..=n / 2 {
                let mut q = vec![0i64; p_n.len()];
                for j in 0..p_n.len() {
                    let s = n as i64 - 2 * j as i64 - 2 * r as i64;
                    let c = binom(h * r as i64, j as i64) * binom(j as i64 + 1, s);
                    q[j] = c;
                    p_n[j] += c;
                }
                while q.len() > 1 && *q.last().unwrap() == 0 {
                    q.pop();
                }
                if !poly_is_trivial(&q) {
                    q_full.push((r as i64, q));
                }
            }
            while p_n.len() > 1 && *p_n.last().unwrap() == 0 {
                p_n.pop();
            }

            // Check: does P_n interlace each Q_{n,r}?
            // (This would mean P_n is a common interlacer)
            if p_n.len() > 2 {
                for (r_val, ref q) in &q_full {
                    if q.len() > 2 {
                        // Check if q interlaces p_n (q has smaller degree)
                        let result = check_weak_interlacing(q, &p_n);
                        if result != Some(true) {
                            if all_common {
                                println!(
                                    "  Common interlacing FAILS: Q_{{{},{}}} vs P_{}: {:?}",
                                    n, r_val, n, result
                                );
                            }
                            all_common = false;
                        }
                    }
                }
            }
        }

        if all_interlace {
            println!("  Consecutive Q_{{n,r}} interlacing: PASS (n <= {})", max_n);
        }
        if all_common {
            println!(
                "  Q_{{n,r}} << P_n (common interlacing): PASS (n <= {})",
                max_n
            );
        }
        println!();
    }
}
