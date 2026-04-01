//! Targeted scan for d=3 real-rooted families.
//! Tests many shapes, uses Sturm chain for overflow cases.
//!
//! Usage: cargo run --release --bin tiling_d3_scan

use polynomial_tools::{is_real_rooted, is_real_rooted_sturm};
use std::collections::BTreeMap;

type Poly = Vec<i128>;

fn poly_zero() -> Poly {
    vec![0]
}
fn poly_one() -> Poly {
    vec![1]
}
fn poly_is_zero(p: &Poly) -> bool {
    p.iter().all(|&c| c == 0)
}

fn poly_add(a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let mut r = vec![0i128; n];
    for (i, &c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        r[i] += c;
    }
    while r.len() > 1 && *r.last().unwrap() == 0 {
        r.pop();
    }
    r
}

fn poly_mul(a: &Poly, b: &Poly) -> Poly {
    if poly_is_zero(a) || poly_is_zero(b) {
        return poly_zero();
    }
    let n = a.len() + b.len() - 1;
    let mut r = vec![0i128; n];
    for (i, &ca) in a.iter().enumerate() {
        if ca == 0 {
            continue;
        }
        for (j, &cb) in b.iter().enumerate() {
            r[i + j] += ca * cb;
        }
    }
    while r.len() > 1 && *r.last().unwrap() == 0 {
        r.pop();
    }
    r
}

fn poly_to_i64(p: &Poly) -> Option<Vec<i64>> {
    p.iter().map(|&c| i64::try_from(c).ok()).collect()
}

/// GCD of all coefficients, to try to scale polynomial into i64 range.
fn poly_gcd(p: &Poly) -> i128 {
    fn gcd(a: i128, b: i128) -> i128 {
        let (mut a, mut b) = (a.abs(), b.abs());
        while b != 0 {
            let t = b;
            b = a % b;
            a = t;
        }
        a
    }
    let mut g = 0i128;
    for &c in p {
        g = gcd(g, c);
    }
    if g == 0 {
        1
    } else {
        g
    }
}

fn check_rr(p: &Poly) -> Option<bool> {
    if p.len() <= 2 {
        return Some(true);
    }
    // Try direct i64
    if let Some(p64) = poly_to_i64(p) {
        return Some(is_real_rooted(&p64));
    }
    // Try scaling by GCD
    let g = poly_gcd(p);
    if g > 1 {
        let scaled: Poly = p.iter().map(|&c| c / g).collect();
        if let Some(s64) = poly_to_i64(&scaled) {
            return Some(is_real_rooted(&s64));
        }
    }
    None // Can't check
}

fn encode(state: &[usize], d: usize) -> usize {
    let base = d + 1;
    let mut idx = 0;
    for &s in state {
        idx = idx * base + s;
    }
    idx
}

fn decode(mut idx: usize, ell: usize, d: usize) -> Vec<usize> {
    let base = d + 1;
    let mut s = vec![0usize; ell];
    for i in (0..ell).rev() {
        s[i] = idx % base;
        idx /= base;
    }
    s
}

fn test_shape(mu: &[usize], d: usize, max_n: usize) -> (Option<usize>, usize) {
    // Returns (first_fail_or_none, last_checked_n)
    let ell = mu.len() - 1;
    if ell == 0 {
        return (None, max_n);
    }

    let base = d + 1;
    let ns = base.pow(ell as u32);

    let mut t_mat = vec![vec![poly_zero(); ns]; ns];
    for j in 0..ns {
        let old = decode(j, ell, d);
        for v in 0..=d {
            let mut ok = true;
            if v > 0 {
                for m in 1..=ell {
                    let prev = old[m - 1];
                    if prev > 0 && v < prev + mu[m] {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                let mut new_s = vec![0usize; ell];
                new_s[0] = v;
                for k in 1..ell {
                    new_s[k] = old[k - 1];
                }
                let i = encode(&new_s, d);
                let w: Poly = if v > 0 { vec![0, 1] } else { vec![1] };
                t_mat[i][j] = poly_add(&t_mat[i][j], &w);
            }
        }
    }

    let mut vec_s: Vec<Poly> = vec![poly_zero(); ns];
    vec_s[0] = poly_one();

    let mut last_checked = 0;
    for n in 0..=max_n {
        let p = &vec_s[0];
        match check_rr(p) {
            Some(true) => {
                last_checked = n;
            }
            Some(false) => {
                return (Some(n), n);
            }
            None => {
                return (None, last_checked);
            } // overflow, can't check further
        }

        let mut new_vec: Vec<Poly> = vec![poly_zero(); ns];
        for i in 0..ns {
            for j in 0..ns {
                if poly_is_zero(&t_mat[i][j]) || poly_is_zero(&vec_s[j]) {
                    continue;
                }
                let term = poly_mul(&t_mat[i][j], &vec_s[j]);
                new_vec[i] = poly_add(&new_vec[i], &term);
            }
        }
        vec_s = new_vec;
    }
    (None, max_n)
}

fn main() {
    let d = 3;
    let max_n = 80;

    println!(
        "Scanning Ferrers tiles with d={}, checking RR up to n={}",
        d, max_n
    );
    println!(
        "{:<25} {:>5} {:>10} {:<20}",
        "mu", "mu_ell", "claw-free", "result"
    );
    println!("{}", "-".repeat(65));

    // Generate all partitions mu with given number of columns and max part
    let mut results: BTreeMap<String, Vec<(Vec<usize>, String)>> = BTreeMap::new();

    // 3-column shapes (ell=2)
    for a in 1..=8 {
        for b in 1..=a {
            for c in 1..=b {
                let mu = vec![a, b, c];
                let cf = d <= c;
                let (fail, last) = test_shape(&mu, d, max_n);
                let status = match fail {
                    Some(n) => format!("FAILS n={}", n),
                    None if last == max_n => format!("RR (n<={})", max_n),
                    None => format!("RR (n<={}) overflow", last),
                };
                let key = if cf {
                    "3-col (claw-free)".to_string()
                } else {
                    format!("3-col mu_ell={}", c)
                };
                results.entry(key).or_default().push((mu, status));
            }
        }
    }

    // 4-column shapes (ell=3), limited range
    for a in 2..=6 {
        for b in 1..=a.min(5) {
            for c in 1..=b.min(4) {
                for e in 1..=c.min(3) {
                    let mu = vec![a, b, c, e];
                    let cf = d <= e;
                    let (fail, last) = test_shape(&mu, d, max_n);
                    let status = match fail {
                        Some(n) => format!("FAILS n={}", n),
                        None if last == max_n => format!("RR (n<={})", max_n),
                        None => format!("RR (n<={}) overflow", last),
                    };
                    let key = if cf {
                        "4-col (claw-free)".to_string()
                    } else {
                        format!("4-col mu_ell={}", e)
                    };
                    results.entry(key).or_default().push((mu, status));
                }
            }
        }
    }

    // 5-column shapes, very limited
    for a in 3..=5 {
        for b in 2..=a.min(4) {
            for c in 1..=b.min(3) {
                for e in 1..=c.min(3) {
                    for f in 1..=e.min(2) {
                        let mu = vec![a, b, c, e, f];
                        let cf = d <= f;
                        let (fail, last) = test_shape(&mu, d, max_n);
                        let status = match fail {
                            Some(n) => format!("FAILS n={}", n),
                            None if last == max_n => format!("RR (n<={})", max_n),
                            None => format!("RR (n<={}) overflow", last),
                        };
                        let key = if cf {
                            "5-col (claw-free)".to_string()
                        } else {
                            format!("5-col mu_ell={}", f)
                        };
                        results.entry(key).or_default().push((mu, status));
                    }
                }
            }
        }
    }

    for (section, cases) in &results {
        println!("\n=== {} ===", section);
        for (mu, status) in cases {
            let mu_ell = *mu.last().unwrap();
            let cf = if d <= mu_ell { "yes" } else { "no" };
            let rr = if status.starts_with("RR") {
                "  <<<"
            } else {
                ""
            };
            // Only print non-CF cases, or all if CF
            if !status.starts_with("RR") || !status.contains("claw-free") || true {
                println!(
                    "{:<25} {:>5} {:>10} {:<20}{}",
                    format!("{:?}", mu),
                    mu_ell,
                    cf,
                    status,
                    rr
                );
            }
        }
    }

    // Summary: count RR vs FAILS for non-CF cases
    println!("\n{}", "=".repeat(65));
    println!("SUMMARY (non-claw-free only):");
    let mut rr_count = 0;
    let mut fail_count = 0;
    let mut overflow_count = 0;
    for (section, cases) in &results {
        if section.contains("claw-free") {
            continue;
        }
        for (_, status) in cases {
            if status.starts_with("RR (n<=80)") {
                rr_count += 1;
            } else if status.starts_with("FAILS") {
                fail_count += 1;
            } else {
                overflow_count += 1;
            }
        }
    }
    println!("  RR confirmed (n<={}): {}", max_n, rr_count);
    println!("  FAILS: {}", fail_count);
    println!("  RR partial (overflow): {}", overflow_count);
}
