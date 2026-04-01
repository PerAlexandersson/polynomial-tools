//! Check if tiling polynomials form a generalized Sturm sequence:
//! P_n << P_{n+1} for all n (interlacing or alternating left).
//!
//! Usage: cargo run --release --bin tiling_sturm -- 3,2,2,1 4 80

use num_rational::BigRational as Q;
use polynomial_tools::{is_real_rooted, real_roots};
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

fn deg(p: &Poly) -> usize {
    let mut d = p.len() - 1;
    while d > 0 && p[d] == 0 {
        d -= 1;
    }
    d
}

/// Check f << g (f interleaves g) using exact rational arithmetic.
/// Returns true if f interlaces g (deg f + 1 = deg g) or f alternates left of g (deg f = deg g).
fn check_interleaves(f: &[i64], g: &[i64]) -> Option<bool> {
    let df = {
        let mut d = f.len() - 1;
        while d > 0 && f[d] == 0 {
            d -= 1;
        }
        d
    };
    let dg = {
        let mut d = g.len() - 1;
        while d > 0 && g[d] == 0 {
            d -= 1;
        }
        d
    };

    if df == 0 && dg == 0 {
        return Some(true);
    }
    if df == 0 {
        return Some(true);
    } // 0 or const << anything
    if dg == 0 {
        return Some(true);
    } // anything << 0 or const

    if df + 1 != dg && df != dg {
        return Some(false);
    }

    let rf = real_roots(f)?;
    let rg = real_roots(g)?;

    if rf.len() != df || rg.len() != dg {
        return Some(false);
    }

    // Sort roots (they should already be sorted from real_roots, but be safe)
    let mut rf = rf;
    let mut rg = rg;
    rf.sort();
    rg.sort();

    if df + 1 == dg {
        // f interlaces g: g1 <= f1 <= g2 <= f2 <= ... <= f_{d-1} <= gd
        for i in 0..df {
            if rf[i] < rg[i] {
                return Some(false);
            }
            if rf[i] > rg[i + 1] {
                return Some(false);
            }
        }
        Some(true)
    } else {
        // f alternates left of g: f1 <= g1 <= f2 <= g2 <= ...
        for i in 0..df {
            if rf[i] > rg[i] {
                return Some(false);
            }
            if i + 1 < df && rg[i] > rf[i + 1] {
                return Some(false);
            }
        }
        Some(true)
    }
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

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mu: Vec<usize> = args[0]
        .split(',')
        .map(|s| s.trim().parse().unwrap())
        .collect();
    let d: usize = args[1].parse().unwrap();
    let max_n: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(80);

    let ell = mu.len() - 1;
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

    let mut prev_p64: Option<Vec<i64>> = None;
    let mut all_interleave = true;
    let mut first_fail_il = None;
    let mut first_fail_rr = None;

    println!(
        "mu={:?}, d={}, checking Sturm sequence up to n={}",
        mu, d, max_n
    );
    println!("{:>5} {:>4} {:>8} {:>12}", "n", "deg", "RR", "P_n-1<<P_n");

    for n in 0..=max_n {
        let p = vec_s[0].clone();
        let p64 = poly_to_i64(&p);
        let d_n = deg(&p);

        let rr = match &p64 {
            Some(v) if d_n >= 2 => {
                let r = is_real_rooted(v);
                if !r && first_fail_rr.is_none() {
                    first_fail_rr = Some(n);
                }
                if r {
                    "yes"
                } else {
                    "NO"
                }
            }
            Some(_) => "triv",
            None => "ovfl",
        };

        let il = if let (Some(prev), Some(curr)) = (&prev_p64, &p64) {
            match check_interleaves(prev, curr) {
                Some(true) => "yes",
                Some(false) => {
                    if first_fail_il.is_none() && n > 1 {
                        first_fail_il = Some(n);
                        all_interleave = false;
                    }
                    "NO"
                }
                None => "N/A",
            }
        } else if prev_p64.is_none() && n > 0 {
            "ovfl"
        } else {
            "—"
        };

        if n <= 20 || (rr == "NO" || il == "NO") || n == max_n || n % 10 == 0 {
            println!("{:>5} {:>4} {:>8} {:>12}", n, d_n, rr, il);
        }

        prev_p64 = p64;

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

    println!("\nSummary:");
    if let Some(n) = first_fail_rr {
        println!("  Real-rootedness FAILS at n={}", n);
    } else {
        println!("  All polynomials real-rooted up to n={}", max_n);
    }
    if let Some(n) = first_fail_il {
        println!("  Interlacing FAILS at n={}", n);
    } else if all_interleave {
        println!(
            "  All consecutive pairs interlace (Sturm sequence) up to n={}",
            max_n
        );
    }
}
