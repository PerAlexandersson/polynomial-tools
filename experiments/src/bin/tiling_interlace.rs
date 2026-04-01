//! Check pairwise interlacing of tiling polynomials P_n(t) and P_{n+1}(t).
//!
//! Usage:
//!   cargo run --release --bin tiling_interlace -- 2,2,2 3
//!   cargo run --release --bin tiling_interlace -- 3,2,1 2

use polynomial_tools::{check_interlacing, check_weak_interlacing, format_poly, is_real_rooted};
use std::env;

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

/// Compute tiling polynomials and check interlacing between consecutive ones.
fn analyze(mu: &[usize], d: usize, max_n: usize) {
    let ell = mu.len() - 1;
    println!("mu = {:?}, d = {}, ell = {}", mu, d, ell);

    if ell == 0 {
        println!("Single column — trivially real-rooted.");
        return;
    }

    let base = d + 1;
    let ns = base.pow(ell as u32);

    // Build transfer matrix
    let mut t_mat = vec![vec![poly_zero(); ns]; ns];
    for j in 0..ns {
        let old = decode_state(j, ell, d);
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
                let mut new_state = vec![0usize; ell];
                new_state[0] = v;
                for k in 1..ell {
                    new_state[k] = old[k - 1];
                }
                let i = encode_state(&new_state, d);
                let w: Poly = if v > 0 { vec![0, 1] } else { vec![1] };
                t_mat[i][j] = poly_add(&t_mat[i][j], &w);
            }
        }
    }

    let zero_idx = 0;
    let mut vec_state: Vec<Poly> = vec![poly_zero(); ns];
    vec_state[zero_idx] = poly_one();

    let mut prev_p64: Option<Vec<i64>> = None;

    for n in 0..=max_n {
        let p = vec_state[zero_idx].clone();
        let p64 = poly_to_i64(&p);
        let deg = p.len() - 1;

        // Check real-rootedness
        let rr = match &p64 {
            Some(v) if deg >= 2 => {
                if is_real_rooted(v) {
                    "RR"
                } else {
                    "NOT-RR"
                }
            }
            Some(_) => "triv",
            None => "overflow",
        };

        // Check interlacing P_{n-1} << P_n
        let interlace = if let (Some(prev), Some(curr)) = (&prev_p64, &p64) {
            if prev.len() >= 2 && curr.len() >= 2 {
                let strict = check_interlacing(curr, prev);
                let weak = check_weak_interlacing(curr, prev);
                match (strict, weak) {
                    (Some(true), _) => "strict-IL",
                    (_, Some(true)) => "weak-IL",
                    (Some(false), Some(false)) => "NOT-IL",
                    _ => "N/A",
                }
            } else {
                "triv"
            }
        } else if prev_p64.is_none() && p64.is_some() {
            "—"
        } else {
            "overflow"
        };

        if deg <= 6 {
            if let Some(ref v) = p64 {
                println!("P_{:>3} = {:40}  {}  {}", n, format_poly(v), rr, interlace);
            } else {
                println!("P_{:>3}: deg={:>3}  {}  {}", n, deg, rr, interlace);
            }
        } else {
            let total: i128 = p.iter().sum();
            println!(
                "P_{:>3}: deg={:>3}, P(1)={:<20}  {}  {}",
                n, deg, total, rr, interlace
            );
        }

        prev_p64 = p64;

        // Advance
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
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.len() < 2 {
        eprintln!("Usage: tiling_interlace <mu> <d> [max_n]");
        eprintln!("  e.g.: tiling_interlace 2,2,2 3 60");
        return;
    }

    let mu: Vec<usize> = args[0]
        .split(',')
        .map(|s| s.trim().parse().expect("bad mu"))
        .collect();
    let d: usize = args[1].parse().expect("bad d");
    let max_n: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(70);

    analyze(&mu, d, max_n);
}
