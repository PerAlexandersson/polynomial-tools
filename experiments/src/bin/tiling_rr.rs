//! Search for Ferrers tiles with 3+ columns giving real-rooted tiling polynomials.
//!
//! Uses the transfer matrix on anchor word states to compute P_{mu,d,n}(t),
//! then checks real-rootedness via Bézout matrix.
//!
//! Usage:
//!   cargo run --release --bin tiling_rr
//!   cargo run --release --bin tiling_rr -- 3,3,3 5    # specific mu and max d
//!   cargo run --release --bin tiling_rr -- scan       # systematic scan

use polynomial_tools::{format_poly, is_real_rooted};
use std::env;

/// Polynomial represented as Vec<i128> to avoid overflow for large n.
/// Converted to Vec<i64> only when calling is_real_rooted.
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

/// Convert i128 poly to i64, returning None if any coefficient overflows.
fn poly_to_i64(p: &Poly) -> Option<Vec<i64>> {
    p.iter().map(|&c| i64::try_from(c).ok()).collect()
}

/// Encode a state (tuple of ell values, each in 0..=d) as a single usize.
/// state[0] = most recent value, state[ell-1] = oldest.
fn encode_state(state: &[usize], d: usize) -> usize {
    let base = d + 1;
    let mut idx = 0;
    for &s in state {
        idx = idx * base + s;
    }
    idx
}

fn decode_state(mut idx: usize, ell: usize, d: usize) -> Vec<usize> {
    let base = d + 1;
    let mut state = vec![0usize; ell];
    for i in (0..ell).rev() {
        state[i] = idx % base;
        idx /= base;
    }
    state
}

/// Build transfer matrix for anchor words.
///
/// State: last `ell` values, each in {0,1,...,d} where 0 = infinity.
/// T[new_state][old_state] = polynomial weight for the transition.
///
/// mu: Ferrers shape (weakly decreasing), mu[0] >= mu[1] >= ... >= mu[ell]
/// d: number of available rows for anchor placement
fn build_transfer_matrix(mu: &[usize], d: usize) -> (Vec<Vec<Poly>>, usize) {
    let ell = mu.len() - 1;
    let base = d + 1;
    let ns = base.pow(ell as u32);

    let mut t_mat = vec![vec![poly_zero(); ns]; ns];

    for j in 0..ns {
        let old = decode_state(j, ell, d);

        for v in 0..=d {
            let mut ok = true;
            if v > 0 {
                // Check: for m=1..ell, if old[m-1] != 0 then v >= old[m-1] + mu[m]
                for m in 1..=ell {
                    let prev = old[m - 1];
                    if prev > 0 && v < prev + mu[m] {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                // New state: (v, old[0], ..., old[ell-2])
                let mut new_state = vec![0usize; ell];
                new_state[0] = v;
                for k in 1..ell {
                    new_state[k] = old[k - 1];
                }
                let i = encode_state(&new_state, d);

                let weight: Poly = if v > 0 { vec![0, 1] } else { vec![1] }; // t or 1
                t_mat[i][j] = poly_add(&t_mat[i][j], &weight);
            }
        }
    }

    (t_mat, ns)
}

/// Compute tiling polynomials P_{mu,d,n}(t) for n=0..max_n.
/// Returns the polynomials and the first n where real-rootedness fails (if any).
fn compute_and_check(mu: &[usize], d: usize, max_n: usize) -> (Vec<Poly>, Option<usize>) {
    let ell = mu.len() - 1;
    if ell == 0 {
        // Single column: P_n(t) = (1 + d*t)^n, always real-rooted
        let mut polys = Vec::new();
        let mut p = poly_one();
        let factor: Poly = vec![1, d as i128]; // 1 + d*t
        for _ in 0..=max_n {
            polys.push(p.clone());
            p = poly_mul(&p, &factor);
        }
        return (polys, None);
    }

    let base = d + 1;
    let ns = base.pow(ell as u32);
    let (t_mat, _) = build_transfer_matrix(mu, d);

    let zero_idx = 0; // state (0,0,...,0) = all infinity

    // State vector
    let mut vec_state: Vec<Poly> = vec![poly_zero(); ns];
    vec_state[zero_idx] = poly_one();

    let mut polys = Vec::new();
    let mut first_fail: Option<usize> = None;

    for n in 0..=max_n {
        let p = vec_state[zero_idx].clone();
        polys.push(p.clone());

        // Check real-rootedness
        if first_fail.is_none() && p.len() > 2 {
            match poly_to_i64(&p) {
                Some(p64) => {
                    if !is_real_rooted(&p64) {
                        first_fail = Some(n);
                    }
                }
                None => {
                    // Coefficients too large for i64 — skip this n
                    // (could use BigInt path, but for now mark as overflow)
                    first_fail = Some(n); // conservative: stop here
                                          // Distinguish overflow from genuine failure
                    return (polys, Some(n + 10000)); // signal overflow with large n
                }
            }
        }

        // Advance: new_vec = T * vec
        let mut new_vec: Vec<Poly> = vec![poly_zero(); ns];
        for i in 0..ns {
            for j in 0..ns {
                if poly_is_zero(&t_mat[i][j]) || poly_is_zero(&vec_state[j]) {
                    continue;
                }
                let term = poly_mul(&t_mat[i][j], &vec_state[j]);
                new_vec[i] = poly_add(&new_vec[i], &term);
            }
        }
        vec_state = new_vec;
    }

    (polys, first_fail)
}

/// Systematic scan over many Ferrers shapes.
fn scan(max_n: usize) {
    println!(
        "{:<22} {:>3} {:>10} {:<30}",
        "mu", "d", "claw-free", "result"
    );
    println!("{}", "-".repeat(72));

    let test_cases: Vec<(&str, Vec<(Vec<usize>, Vec<usize>)>)> = vec![
        (
            "3-col rectangles",
            vec![
                (vec![2, 2, 2], (1..=5).collect()),
                (vec![3, 3, 3], (1..=6).collect()),
                (vec![4, 4, 4], (1..=6).collect()),
                (vec![5, 5, 5], (1..=6).collect()),
            ],
        ),
        (
            "3-col hooks",
            vec![
                (vec![3, 1, 1], (1..=5).collect()),
                (vec![4, 1, 1], (1..=5).collect()),
                (vec![5, 1, 1], (1..=5).collect()),
                (vec![6, 1, 1], (1..=5).collect()),
            ],
        ),
        (
            "3-col staircase/other",
            vec![
                (vec![3, 2, 1], (1..=4).collect()),
                (vec![4, 3, 2], (1..=4).collect()),
                (vec![5, 4, 3], (1..=4).collect()),
                (vec![5, 3, 1], (1..=4).collect()),
                (vec![4, 2, 1], (1..=4).collect()),
            ],
        ),
        (
            "3-col with mu_ell=2",
            vec![
                (vec![3, 2, 2], (1..=4).collect()),
                (vec![3, 3, 2], (1..=4).collect()),
                (vec![4, 2, 2], (1..=4).collect()),
                (vec![4, 3, 2], (1..=4).collect()),
                (vec![4, 4, 2], (1..=4).collect()),
                (vec![4, 3, 3], (1..=4).collect()),
                (vec![4, 4, 3], (1..=5).collect()),
                (vec![5, 3, 2], (1..=4).collect()),
                (vec![5, 4, 2], (1..=4).collect()),
                (vec![5, 5, 2], (1..=4).collect()),
                (vec![5, 4, 3], (1..=5).collect()),
                (vec![5, 5, 3], (1..=5).collect()),
                (vec![5, 5, 4], (1..=5).collect()),
                (vec![6, 4, 2], (1..=4).collect()),
                (vec![6, 5, 3], (1..=5).collect()),
                (vec![6, 6, 4], (1..=5).collect()),
            ],
        ),
        (
            "4-col shapes",
            vec![
                (vec![2, 2, 2, 2], (1..=3).collect()),
                (vec![3, 3, 3, 3], (1..=4).collect()),
                (vec![4, 4, 4, 4], (1..=5).collect()),
                (vec![3, 2, 2, 1], (1..=3).collect()),
                (vec![3, 2, 2, 2], (1..=3).collect()),
                (vec![3, 3, 2, 2], (1..=3).collect()),
                (vec![3, 3, 3, 2], (1..=3).collect()),
                (vec![4, 3, 2, 1], (1..=3).collect()),
                (vec![4, 3, 2, 2], (1..=3).collect()),
                (vec![4, 3, 3, 2], (1..=3).collect()),
                (vec![4, 3, 3, 3], (1..=4).collect()),
                (vec![4, 4, 3, 2], (1..=3).collect()),
                (vec![4, 4, 3, 3], (1..=4).collect()),
                (vec![4, 4, 4, 3], (1..=4).collect()),
                (vec![3, 2, 1, 1], (1..=3).collect()),
                (vec![4, 2, 1, 1], (1..=3).collect()),
                (vec![3, 1, 1, 1], (1..=3).collect()),
                (vec![4, 1, 1, 1], (1..=3).collect()),
                (vec![5, 3, 3, 2], (1..=3).collect()),
                (vec![5, 4, 3, 2], (1..=3).collect()),
                (vec![5, 5, 4, 3], (1..=4).collect()),
            ],
        ),
        (
            "5-col shapes",
            vec![
                (vec![2, 2, 2, 2, 2], (1..=2).collect()),
                (vec![3, 3, 3, 3, 3], (1..=3).collect()),
                (vec![3, 2, 2, 1, 1], (1..=2).collect()),
                (vec![3, 2, 1, 1, 1], (1..=2).collect()),
                (vec![4, 3, 2, 1, 1], (1..=2).collect()),
                (vec![4, 4, 3, 2, 2], (1..=2).collect()),
                (vec![4, 3, 3, 2, 2], (1..=2).collect()),
            ],
        ),
    ];

    for (section, cases) in &test_cases {
        println!("\n--- {} ---", section);
        for (mu, d_values) in cases {
            for &d in d_values {
                let mu_ell = *mu.last().unwrap();
                let cf = if d <= mu_ell { "yes" } else { "no" };
                let mu_str = format!("{:?}", mu);

                let (_, first_fail) = compute_and_check(mu, d, max_n);
                let status = match first_fail {
                    None => format!("REAL-ROOTED (n<={})", max_n),
                    Some(n) if n >= 10000 => format!("OVERFLOW at n={}", n - 10000),
                    Some(n) => format!("FAILS at n={}", n),
                };
                println!("{:<22} {:>3} {:>10} {}", mu_str, d, cf, status);
            }
        }
    }
}

/// Detailed output for a single (mu, d) pair.
fn detailed(mu: &[usize], d: usize, max_n: usize) {
    println!("Ferrers tile mu = {:?}, d = {}", mu, d);
    let mu_ell = *mu.last().unwrap();
    println!(
        "ell = {}, mu_ell = {}, claw-free (d <= mu_ell): {}",
        mu.len() - 1,
        mu_ell,
        if d <= mu_ell { "yes" } else { "no" }
    );
    println!();

    let (polys, first_fail) = compute_and_check(mu, d, max_n);

    for (n, p) in polys.iter().enumerate() {
        let p64 = poly_to_i64(p);
        let rr = if p.len() <= 2 {
            "  (trivial)".to_string()
        } else if let Some(ref p64v) = p64 {
            if is_real_rooted(p64v) {
                "  RR".to_string()
            } else {
                "  NOT RR  <---".to_string()
            }
        } else {
            "  (overflow)".to_string()
        };
        let deg = p.len() - 1;
        if deg <= 8 {
            if let Some(ref p64v) = p64 {
                println!("P_{:>3}(t) = {}{}", n, format_poly(p64v), rr);
            } else {
                println!("P_{:>3}(t): deg={}, (coeffs overflow i64){}", n, deg, rr);
            }
        } else {
            let total: i128 = p.iter().sum();
            println!("P_{:>3}(t): deg={}, P(1)={}{}", n, deg, total, rr);
        }
    }

    println!();
    match first_fail {
        None => println!("All polynomials real-rooted up to n={}.", max_n),
        Some(n) if n >= 10000 => println!("Coefficients overflow i64 at n={}.", n - 10000),
        Some(n) => println!("First failure of real-rootedness at n={}.", n),
    }
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();

    if args.is_empty() || args[0] == "scan" {
        let max_n = if args.len() > 1 {
            args[1].parse().unwrap_or(80)
        } else {
            80
        };
        scan(max_n);
    } else {
        // Parse mu from first arg: "3,3,3" or "3,2,1"
        let mu: Vec<usize> = args[0]
            .split(',')
            .map(|s| s.trim().parse().expect("invalid mu entry"))
            .collect();

        let max_d: usize = if args.len() > 1 {
            args[1].parse().unwrap_or(5)
        } else {
            5
        };

        let max_n: usize = if args.len() > 2 {
            args[2].parse().unwrap_or(80)
        } else {
            80
        };

        for d in 1..=max_d {
            detailed(&mu, d, max_n);
            println!("\n{}\n", "=".repeat(60));
        }
    }
}
