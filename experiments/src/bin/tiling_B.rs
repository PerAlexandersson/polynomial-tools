//! Compute the fault-free generating function B(x,t) for a given mu and d
//! by extracting the recursion coefficients from the tiling polynomials.
//! Usage: cargo run --release --bin tiling_B -- 3,2,2,1 4

use std::env;

type Poly = Vec<i128>;

fn poly_zero() -> Poly { vec![0] }
fn poly_one() -> Poly { vec![1] }
fn poly_is_zero(p: &Poly) -> bool { p.iter().all(|&c| c == 0) }
fn poly_trim(p: &mut Poly) { while p.len() > 1 && *p.last().unwrap() == 0 { p.pop(); } }

fn poly_add(a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let mut r = vec![0i128; n];
    for (i, &c) in a.iter().enumerate() { r[i] += c; }
    for (i, &c) in b.iter().enumerate() { r[i] += c; }
    poly_trim(&mut r); r
}

fn poly_sub(a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let mut r = vec![0i128; n];
    for (i, &c) in a.iter().enumerate() { r[i] += c; }
    for (i, &c) in b.iter().enumerate() { r[i] -= c; }
    poly_trim(&mut r); r
}

fn poly_mul(a: &Poly, b: &Poly) -> Poly {
    if poly_is_zero(a) || poly_is_zero(b) { return poly_zero(); }
    let n = a.len() + b.len() - 1;
    let mut r = vec![0i128; n];
    for (i, &ca) in a.iter().enumerate() {
        if ca == 0 { continue; }
        for (j, &cb) in b.iter().enumerate() { r[i + j] += ca * cb; }
    }
    poly_trim(&mut r); r
}

fn poly_scale(p: &Poly, s: i128) -> Poly {
    p.iter().map(|&c| c * s).collect()
}

fn format_poly(p: &Poly) -> String {
    if poly_is_zero(p) { return "0".into(); }
    let mut terms = vec![];
    for (i, &c) in p.iter().enumerate() {
        if c == 0 { continue; }
        let t_part = match i {
            0 => format!("{}", c),
            1 => if c == 1 { "t".into() } else { format!("{}t", c) },
            _ => if c == 1 { format!("t^{}", i) } else { format!("{}t^{}", c, i) },
        };
        terms.push(t_part);
    }
    terms.join(" + ")
}

fn encode(state: &[usize], d: usize) -> usize {
    let base = d + 1;
    let mut idx = 0;
    for &s in state { idx = idx * base + s; }
    idx
}

fn decode(mut idx: usize, ell: usize, d: usize) -> Vec<usize> {
    let base = d + 1;
    let mut s = vec![0usize; ell];
    for i in (0..ell).rev() { s[i] = idx % base; idx /= base; }
    s
}

fn compute_polys(mu: &[usize], d: usize, max_n: usize) -> Vec<Poly> {
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
                    if prev > 0 && v < prev + mu[m] { ok = false; break; }
                }
            }
            if ok {
                let mut new_s = vec![0usize; ell];
                new_s[0] = v;
                for k in 1..ell { new_s[k] = old[k - 1]; }
                let i = encode(&new_s, d);
                let w: Poly = if v > 0 { vec![0, 1] } else { vec![1] };
                t_mat[i][j] = poly_add(&t_mat[i][j], &w);
            }
        }
    }

    let mut vec_s: Vec<Poly> = vec![poly_zero(); ns];
    vec_s[0] = poly_one();
    let mut polys = vec![];

    for _n in 0..=max_n {
        polys.push(vec_s[0].clone());
        let mut new_vec: Vec<Poly> = vec![poly_zero(); ns];
        for i in 0..ns {
            for j in 0..ns {
                if poly_is_zero(&t_mat[i][j]) || poly_is_zero(&vec_s[j]) { continue; }
                let term = poly_mul(&t_mat[i][j], &vec_s[j]);
                new_vec[i] = poly_add(&new_vec[i], &term);
            }
        }
        vec_s = new_vec;
    }
    polys
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mu: Vec<usize> = args[0].split(',').map(|s| s.trim().parse().unwrap()).collect();
    let d: usize = args[1].parse().unwrap();
    let max_n: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(40);

    println!("mu = {:?}, d = {}", mu, d);
    let polys = compute_polys(&mu, d, max_n);

    // Extract B coefficients: B_j such that P_n = sum_j B_j * P_{n-j}
    // B_1 = 1 (empty column), B_2 = ... = B_{ell} = 0
    let ell = mu.len() - 1;
    let mut b_coeffs: Vec<Poly> = vec![poly_zero(); max_n + 1];
    b_coeffs[1] = poly_one();

    println!("\nFault-free generating function B(x,t):");
    println!("  B_1 = 1  (empty column, x)");

    for n in 2..=max_n {
        // P_n = sum_{j=1}^{n} B_j * P_{n-j}
        // B_n = P_n - sum_{j=1}^{n-1} B_j * P_{n-j}
        let mut rhs = polys[n].clone();
        for j in 1..n {
            if poly_is_zero(&b_coeffs[j]) { continue; }
            let term = poly_mul(&b_coeffs[j], &polys[n - j]);
            rhs = poly_sub(&rhs, &term);
        }
        poly_trim(&mut rhs);
        if !poly_is_zero(&rhs) {
            b_coeffs[n] = rhs.clone();
            println!("  B_{} = {}  (x^{})", n, format_poly(&rhs), n);
        }
    }

    // Print the full B
    println!("\nB(x,t) = x");
    for n in 2..=max_n {
        if !poly_is_zero(&b_coeffs[n]) {
            println!("       + ({}) x^{}", format_poly(&b_coeffs[n]), n);
        }
    }
}
