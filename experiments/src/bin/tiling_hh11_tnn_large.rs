//! Check TNN via Neville elimination for mu=(h,h,1,1), d=h
//! with large h and large n.
//!
//! Usage: cargo run --release --bin tiling_hh11_tnn_large [max_h] [max_n]

use num_bigint::BigInt;
use polynomial_tools::check_tnn_neville_bigint;
use std::time::Instant;

type Poly = Vec<i128>;

fn poly_zero() -> Poly {
    vec![0]
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

fn compute_pn(mu: &[usize], d: usize, max_n: usize) -> Vec<Poly> {
    let ell = mu.len() - 1;
    let base = d + 1;
    let ns = base.pow(ell as u32);
    let decode = |mut idx: usize| -> Vec<usize> {
        let mut s = vec![0usize; ell];
        for i in (0..ell).rev() {
            s[i] = idx % base;
            idx /= base;
        }
        s
    };
    let encode = |state: &[usize]| -> usize {
        let mut idx = 0;
        for &s in state {
            idx = idx * base + s;
        }
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
                let i = encode(&new_s);
                let w: Poly = if v > 0 { vec![0, 1] } else { vec![1] };
                t_mat[i][j] = poly_add(&t_mat[i][j], &w);
            }
        }
    }
    let mut results = Vec::with_capacity(max_n + 1);
    let mut vec_s: Vec<Poly> = vec![poly_zero(); ns];
    vec_s[0] = vec![1];
    for n in 0..=max_n {
        results.push(vec_s[0].clone());
        if n < max_n {
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
    }
    results
}

fn to_bigint_matrix(pn: &[Poly]) -> Vec<Vec<BigInt>> {
    let max_j = pn.iter().map(|p| p.len()).max().unwrap_or(1);
    pn.iter()
        .map(|p| {
            let mut row = vec![BigInt::from(0); max_j];
            for (j, &c) in p.iter().enumerate() {
                row[j] = BigInt::from(c);
            }
            row
        })
        .collect()
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let max_h: usize = args.get(0).and_then(|s| s.parse().ok()).unwrap_or(10);
    let max_n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(200);

    for h in 2..=max_h {
        let mu = vec![h, h, 1, 1];
        let t0 = Instant::now();
        let pn = compute_pn(&mu, h, max_n);
        let compute_time = t0.elapsed();

        let mat = to_bigint_matrix(&pn);
        let nrows = mat.len();
        let ncols = mat[0].len();

        let t1 = Instant::now();
        let result = check_tnn_neville_bigint(&mat);
        let check_time = t1.elapsed();

        match &result {
            Ok(()) => println!(
                "h={:>2}: TNN PASS  ({}x{} matrix, compute {:.2?}, Neville {:.2?})",
                h, nrows, ncols, compute_time, check_time
            ),
            Err(msg) => println!(
                "h={:>2}: TNN FAIL  ({}x{} matrix, compute {:.2?}, Neville {:.2?}) — {}",
                h, nrows, ncols, compute_time, check_time, msg
            ),
        }
    }
}
