//! Print a specific tiling polynomial P_n(t) for given mu, d, n.
//! Usage: cargo run --release --bin tiling_print_poly -- 2,2,2 3 56

use std::env;

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mu: Vec<usize> = args[0].split(',').map(|s| s.trim().parse().unwrap()).collect();
    let d: usize = args[1].parse().unwrap();
    let target_n: usize = args[2].parse().unwrap();

    let ell = mu.len() - 1;
    let base = d + 1;
    let ns = base.pow(ell as u32);

    // Build transfer matrix (i128 polynomials)
    type Poly = Vec<i128>;
    let poly_zero = || -> Poly { vec![0] };
    let poly_add = |a: &Poly, b: &Poly| -> Poly {
        let n = a.len().max(b.len());
        let mut r = vec![0i128; n];
        for (i, &c) in a.iter().enumerate() { r[i] += c; }
        for (i, &c) in b.iter().enumerate() { r[i] += c; }
        while r.len() > 1 && *r.last().unwrap() == 0 { r.pop(); }
        r
    };
    let poly_mul = |a: &Poly, b: &Poly| -> Poly {
        if a.iter().all(|&c| c == 0) || b.iter().all(|&c| c == 0) { return vec![0]; }
        let n = a.len() + b.len() - 1;
        let mut r = vec![0i128; n];
        for (i, &ca) in a.iter().enumerate() {
            if ca == 0 { continue; }
            for (j, &cb) in b.iter().enumerate() { r[i + j] += ca * cb; }
        }
        while r.len() > 1 && *r.last().unwrap() == 0 { r.pop(); }
        r
    };

    let decode = |mut idx: usize| -> Vec<usize> {
        let mut s = vec![0usize; ell];
        for i in (0..ell).rev() { s[i] = idx % base; idx /= base; }
        s
    };
    let encode = |state: &[usize]| -> usize {
        let mut idx = 0;
        for &s in state { idx = idx * base + s; }
        idx
    };

    let mut t_mat = vec![vec![poly_zero(); ns]; ns];
    for j in 0..ns {
        let old = decode(j);
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
                let i = encode(&new_s);
                let w: Poly = if v > 0 { vec![0, 1] } else { vec![1] };
                t_mat[i][j] = poly_add(&t_mat[i][j], &w);
            }
        }
    }

    let mut vec_s: Vec<Poly> = vec![poly_zero(); ns];
    vec_s[0] = vec![1];

    for n in 0..=target_n {
        if n == target_n {
            let p = &vec_s[0];
            let deg = p.len() - 1;
            println!("mu={:?}, d={}, n={}, deg={}", mu, d, n, deg);
            // Print coefficients
            for (j, &c) in p.iter().enumerate() {
                if c != 0 {
                    println!("  [t^{}] = {}", j, c);
                }
            }
            // Print in polynomial form (reversed, leading term first)
            let mut terms = Vec::new();
            for j in (0..=deg).rev() {
                if p[j] != 0 {
                    if j == 0 { terms.push(format!("{}", p[j])); }
                    else if j == 1 { terms.push(format!("{}t", p[j])); }
                    else { terms.push(format!("{}t^{{{}}}", p[j], j)); }
                }
            }
            println!("\nP_{{{}}}(t) = {}", n, terms.join(" + "));
        }

        let mut new_vec: Vec<Poly> = vec![poly_zero(); ns];
        for i in 0..ns {
            for j in 0..ns {
                if t_mat[i][j].iter().all(|&c| c == 0) || vec_s[j].iter().all(|&c| c == 0) { continue; }
                let term = poly_mul(&t_mat[i][j], &vec_s[j]);
                new_vec[i] = poly_add(&new_vec[i], &term);
            }
        }
        vec_s = new_vec;
    }
}
