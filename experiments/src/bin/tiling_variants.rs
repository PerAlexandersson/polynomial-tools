//! Compute and check real-rootedness of:
//!   1) Dense tiling polynomials Q_n(s) = [x^n] 1/(1 - s*(B(x,1)-x))
//!   2) Row-1-refined polynomials S_n(s): weight s for row-1 tiles, 1 for others
//!   3) Fault-line polynomials R_n(t) = [x^n] 1/(t*(1 - t*B(x,1)))
//!
//! Usage: cargo run --release --bin tiling_variants -- 3,2,1 3

use polynomial_tools::is_real_rooted;
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

fn format_poly(p: &Poly) -> String {
    if poly_is_zero(p) {
        return "0".into();
    }
    let mut terms = vec![];
    for (i, &c) in p.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let s = match i {
            0 => format!("{}", c),
            1 if c == 1 => "s".into(),
            1 if c == -1 => "-s".into(),
            1 => format!("{}s", c),
            _ if c == 1 => format!("s^{}", i),
            _ if c == -1 => format!("-s^{}", i),
            _ => format!("{}s^{}", c, i),
        };
        terms.push(s);
    }
    if terms.is_empty() {
        "0".into()
    } else {
        terms.join(" + ")
    }
}

fn check_rr(p: &Poly) -> Option<bool> {
    if deg(p) <= 1 {
        return Some(true);
    }
    if let Some(p64) = poly_to_i64(p) {
        return Some(is_real_rooted(&p64));
    }
    None
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

/// Compute polynomials using transfer matrix.
/// weight_fn(v) returns the polynomial weight for placing value v.
fn compute_with_weights(
    mu: &[usize],
    d: usize,
    max_n: usize,
    weight_fn: impl Fn(usize) -> Poly,
) -> Vec<Poly> {
    let ell = mu.len() - 1;
    let ns = (d + 1).pow(ell as u32);
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
                t_mat[i][j] = poly_add(&t_mat[i][j], &weight_fn(v));
            }
        }
    }
    let mut vec_s: Vec<Poly> = vec![poly_zero(); ns];
    vec_s[0] = poly_one();
    let mut polys = vec![];
    for _ in 0..=max_n {
        polys.push(vec_s[0].clone());
        let mut nv: Vec<Poly> = vec![poly_zero(); ns];
        for i in 0..ns {
            for j in 0..ns {
                if poly_is_zero(&t_mat[i][j]) || poly_is_zero(&vec_s[j]) {
                    continue;
                }
                nv[i] = poly_add(&nv[i], &poly_mul(&t_mat[i][j], &vec_s[j]));
            }
        }
        vec_s = nv;
    }
    polys
}

fn check_sequence(name: &str, polys: &[Poly], max_print: usize) {
    let mut first_fail = None;
    let mut last_ok = 0;

    for (n, p) in polys.iter().enumerate() {
        if n < max_print && deg(p) <= 8 {
            println!("  {}_{} = {}", name, n, format_poly(p));
        }
        match check_rr(p) {
            Some(true) => last_ok = n,
            Some(false) => {
                if first_fail.is_none() {
                    first_fail = Some(n);
                }
                break;
            }
            None => break,
        }
    }
    match first_fail {
        None => println!(
            "  => RR for all n <= {} (deg {})\n",
            last_ok,
            deg(&polys[last_ok])
        ),
        Some(n) => println!("  => FAILS at n={} (deg {})\n", n, deg(&polys[n])),
    }
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mu: Vec<usize> = args[0]
        .split(',')
        .map(|s| s.trim().parse().unwrap())
        .collect();
    let d: usize = args[1].parse().unwrap();
    let max_n: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(40);

    println!("mu={:?}, d={}\n", mu, d);

    // 1) Standard tiling polynomial P_n(t): weight t for any tile
    println!("--- P_n(t): standard tiling polynomial ---");
    let p_polys = compute_with_weights(&mu, d, max_n, |v| {
        if v > 0 {
            vec![0, 1]
        } else {
            vec![1]
        } // t for tiles, 1 for gaps
    });
    check_sequence("P", &p_polys, 12);

    // 2) Row-1-refined S_n(s): weight s for row-1, weight 1 for other rows
    println!("--- S_n(s): row-1 refined (s marks row-1 tiles) ---");
    let s_polys = compute_with_weights(&mu, d, max_n, |v| {
        match v {
            0 => vec![1],    // gap
            1 => vec![0, 1], // row 1: weight s
            _ => vec![1],    // other rows: weight 1
        }
    });
    check_sequence("S", &s_polys, 12);

    // 3) Bottom-row-refined: weight s for row-d (bottom), weight 1 for others
    if d > 1 {
        println!(
            "--- T_n(s): bottom-row refined (s marks row-{} tiles) ---",
            d
        );
        let t_polys = compute_with_weights(&mu, d, max_n, |v| {
            if v == 0 {
                vec![1]
            } else if v == d {
                vec![0, 1]
            }
            // bottom row: weight s
            else {
                vec![1]
            }
        });
        check_sequence("T", &t_polys, 12);
    }

    // 4) Dense tiling: Q_n(s) from B(x,1).
    // Compute B coefficients at t=1 first
    println!("--- Q_n(s): dense tiling (s marks fault-free blocks) ---");
    // B(x,1) = F^{-1} at t=1, which we get from p_polys evaluated at t=1
    let p_at_1: Vec<i128> = p_polys.iter().map(|p| p.iter().sum::<i128>()).collect();
    // Extract B_j from: P_n(1) = sum B_j(1) P_{n-j}(1)
    let mut b_at_1 = vec![0i128; max_n + 1];
    b_at_1[1] = 1; // empty column
    for n in 2..=max_n {
        let mut rhs = p_at_1[n];
        for j in 1..n {
            rhs -= b_at_1[j] * p_at_1[n - j];
        }
        b_at_1[n] = rhs;
    }
    // b(x) = B(x,1) - x: the fault-free block part
    let mut b_ff = b_at_1.clone();
    b_ff[1] = 0; // remove the empty column

    // Q_n(s) = sum over compositions into parts from b_ff, each part weighted by s
    // Recurrence: Q_n = s * sum_j b_ff[j] * Q_{n-j}
    let mut q_polys: Vec<Poly> = vec![poly_zero(); max_n + 1];
    q_polys[0] = poly_one();
    for n in 1..=max_n {
        let mut qn = poly_zero();
        for j in 1..=n {
            if b_ff[j] == 0 {
                continue;
            }
            // s * b_ff[j] * Q_{n-j}: multiply Q_{n-j} by s, then scale by b_ff[j]
            let mut term = vec![0i128; q_polys[n - j].len() + 1];
            for (i, &c) in q_polys[n - j].iter().enumerate() {
                term[i + 1] += c * b_ff[j]; // shift by 1 = multiply by s
            }
            qn = poly_add(&qn, &term);
        }
        while qn.len() > 1 && *qn.last().unwrap() == 0 {
            qn.pop();
        }
        q_polys[n] = qn;
    }
    check_sequence("Q", &q_polys, 15);

    // 5) Middle-row refined for d>=3
    if d >= 3 {
        let mid = (d + 1) / 2;
        println!(
            "--- M_n(s): mid-row refined (s marks row-{} tiles) ---",
            mid
        );
        let m_polys = compute_with_weights(&mu, d, max_n, |v| {
            if v == 0 {
                vec![1]
            } else if v == mid {
                vec![0, 1]
            } else {
                vec![1]
            }
        });
        check_sequence("M", &m_polys, 10);
    }
}
