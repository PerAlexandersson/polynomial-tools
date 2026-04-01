//! Generate integer sequences from tiling counts for OEIS lookup.
//! Outputs: P_n(1), dense Q_n(1), B(x,1) coefficients, and row-1 S_n(1).
//!
//! Usage: cargo run --release --bin tiling_oeis -- 3,2,1 3

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

fn compute_polys(mu: &[usize], d: usize, max_n: usize) -> Vec<Poly> {
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
                let w: Poly = if v > 0 { vec![0, 1] } else { vec![1] };
                t_mat[i][j] = poly_add(&t_mat[i][j], &w);
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

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mu: Vec<usize> = args[0]
        .split(',')
        .map(|s| s.trim().parse().unwrap())
        .collect();
    let d: usize = args[1].parse().unwrap();
    let max_n = 25;

    let polys = compute_polys(&mu, d, max_n);

    // P_n(1): total tilings
    let p1: Vec<i128> = polys.iter().map(|p| p.iter().sum()).collect();
    println!("mu={:?}, d={}", mu, d);
    println!(
        "P_n(1) [total tilings]: {}",
        p1.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    );

    // B coefficients at t=1
    let mut b1 = vec![0i128; max_n + 1];
    b1[1] = 1;
    for n in 2..=max_n {
        let mut rhs = p1[n];
        for j in 1..n {
            rhs -= b1[j] * p1[n - j];
        }
        b1[n] = rhs;
    }
    let b_ff: Vec<i128> = b1.iter().skip(2).copied().take_while(|&x| x != 0).collect();
    println!(
        "B(x,1)-x coeffs [fault-free counts]: {}",
        b_ff.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    );

    // Dense Q_n(1): total dense tilings
    let mut b_dense = b1.clone();
    b_dense[1] = 0;
    let mut q1 = vec![0i128; max_n + 1];
    q1[0] = 1;
    for n in 1..=max_n {
        for j in 1..=n {
            if b_dense[j] != 0 {
                q1[n] += b_dense[j] * q1[n - j];
            }
        }
    }
    println!(
        "Q_n(1) [dense tilings]: {}",
        q1.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    );

    // Fault-line counts R_n(1) = [x^n] 1/(1 - B(x,1))... wait, R counts fault lines.
    // R_n(1) = number of (fault-free blocks) in all tilings = sum over tilings of (fl+1)
    // Actually R_n(1) = P_n(1) since at t=1 the GF is the same.
    // The fault-line polynomial R_n(t) uses different variable. Let me compute it differently.
    // R_n(t) = [x^n] t^{-1}/(1 - t*B(x,1)). So R_n(1) = [x^n] 1/(1-B(x,1)) = P_n(1).
    // Not interesting at t=1.

    // S_n(1) [row-1 count at s=1]: same as P_n(1) since s=1 is uniform.
    // More interesting: coefficient of s^1 in S_n(s) = number of row-1 tiles across all tilings
    // Let me compute S_n(s) and extract coefficients

    // Row-1 refined
    let s_polys = {
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
                    let w: Poly = match v {
                        0 => vec![1],
                        1 => vec![0, 1], // s for row 1
                        _ => vec![1],    // 1 for other rows
                    };
                    t_mat[i][j] = poly_add(&t_mat[i][j], &w);
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
    };

    // Coefficient of s^1 in S_n(s): counts total row-1 tiles across all tilings
    let s1_coeffs: Vec<i128> = s_polys
        .iter()
        .map(|p| if p.len() > 1 { p[1] } else { 0 })
        .collect();
    println!(
        "S_n'(1) [total row-1 tiles]: {}",
        s1_coeffs
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    );

    // Constant term of S_n(s): tilings with NO row-1 tiles (= tilings using only rows 2..d)
    let s0_coeffs: Vec<i128> = s_polys.iter().map(|p| p[0]).collect();
    println!(
        "S_n(0) [tilings avoiding row 1]: {}",
        s0_coeffs
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    );
}
