//! Check if the fault-free coefficient matrix [x^n t^j] B_{mu,d}(x,t)
//! is TNN (via Neville elimination) for various mu with mu_1 >= d.
//!
//! Usage: cargo run --release --bin tiling_B_tnn [max_n]

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
fn poly_sub(a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let mut r = vec![0i128; n];
    for (i, &c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        r[i] -= c;
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

fn extract_b(pn: &[Poly]) -> Vec<Poly> {
    let max_n = pn.len() - 1;
    let mut bn: Vec<Poly> = vec![poly_zero(); max_n + 1];
    for n in 1..=max_n {
        let mut sum = poly_zero();
        for k in 1..n {
            if poly_is_zero(&bn[n - k]) {
                continue;
            }
            let term = poly_mul(&pn[k], &bn[n - k]);
            sum = poly_add(&sum, &term);
        }
        bn[n] = poly_sub(&pn[n], &sum);
    }
    bn
}

fn to_bigint_matrix(bn: &[Poly]) -> Vec<Vec<BigInt>> {
    let max_j = bn.iter().map(|p| p.len()).max().unwrap_or(1);
    bn.iter()
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
    let max_n: usize = args.get(0).and_then(|s| s.parse().ok()).unwrap_or(50);

    // Generate test cases: partitions mu with mu_1 >= d, various ell
    let mut cases: Vec<(Vec<usize>, usize, String)> = Vec::new();

    // (h,h,1,1) with d=h — the proved case
    for h in 2..=6 {
        cases.push((vec![h, h, 1, 1], h, format!("({h},{h},1,1) d={h}")));
    }
    // (h,h,1) with d=h
    for h in 2..=7 {
        cases.push((vec![h, h, 1], h, format!("({h},{h},1) d={h}")));
    }
    // (h,h,1,1,1) with d=h
    for h in 2..=4 {
        cases.push((vec![h, h, 1, 1, 1], h, format!("({h},{h},1,1,1) d={h}")));
    }
    // (h,h,2,1) with d=h
    for h in 2..=5 {
        cases.push((vec![h, h, 2, 1], h, format!("({h},{h},2,1) d={h}")));
    }
    // Rectangles (h,h,h) with d=h
    for h in 2..=4 {
        cases.push((vec![h, h, h], h, format!("({h},{h},{h}) d={h}")));
    }
    // Mixed shapes with mu_1 >= d
    cases.push((vec![3, 3, 2], 3, "(3,3,2) d=3".into()));
    cases.push((vec![4, 3, 2, 1], 3, "(4,3,2,1) d=3".into()));
    cases.push((vec![5, 4, 3, 1], 4, "(5,4,3,1) d=4".into()));
    cases.push((vec![4, 4, 2, 1], 4, "(4,4,2,1) d=4".into()));
    cases.push((vec![3, 3, 1, 1, 1], 3, "(3,3,1,1,1) d=3".into()));
    cases.push((vec![4, 3, 1, 1], 3, "(4,3,1,1) d=3".into()));
    cases.push((vec![3, 2, 1], 2, "(3,2,1) d=2".into()));
    cases.push((vec![4, 3, 2], 3, "(4,3,2) d=3".into()));
    cases.push((vec![4, 4, 3, 2, 1], 4, "(4,4,3,2,1) d=4".into()));
    cases.push((vec![5, 5, 2, 1, 1], 5, "(5,5,2,1,1) d=5".into()));
    cases.push((vec![3, 3, 3, 1], 3, "(3,3,3,1) d=3".into()));
    cases.push((vec![4, 4, 4, 1], 4, "(4,4,4,1) d=4".into()));

    // Cases with mu_1 < d (to see if property still holds)
    cases.push((vec![4, 1, 1], 3, "(4,1,1) d=3 [mu1<d]".into()));
    cases.push((vec![5, 1, 1], 4, "(5,1,1) d=4 [mu1<d]".into()));
    cases.push((vec![3, 2, 1, 1], 3, "(3,2,1,1) d=3 [mu1<d]".into()));
    cases.push((vec![4, 2, 2, 1], 3, "(4,2,2,1) d=3 [mu1<d]".into()));
    cases.push((vec![4, 1, 1, 1], 3, "(4,1,1,1) d=3 [mu1<d]".into()));
    cases.push((vec![5, 2, 1, 1], 4, "(5,2,1,1) d=4 [mu1<d]".into()));

    println!(
        "{:<30} {:>6} {:>8}  {}",
        "Case", "mu1>=d", "Size", "TNN result"
    );
    println!("{}", "-".repeat(75));

    for (mu, d, desc) in &cases {
        let ell = mu.len() - 1;
        let state_size = (d + 1).pow(ell as u32);
        if state_size > 100_000 {
            println!(
                "{:<30} {:>6} {:>8}  skipped (state space {state_size})",
                desc,
                if mu[1] >= *d { "yes" } else { "no" },
                ""
            );
            continue;
        }
        let t0 = Instant::now();
        let pn = compute_pn(mu, *d, max_n);
        let bn = extract_b(&pn);
        let mat = to_bigint_matrix(&bn);
        let nrows = mat.len();
        let ncols = mat[0].len();
        let result = check_tnn_neville_bigint(&mat);
        let elapsed = t0.elapsed();
        let cond = if mu[1] >= *d { "yes" } else { "no" };
        let size_str = format!("{}x{}", nrows, ncols);
        match &result {
            Ok(()) => println!(
                "{:<30} {:>6} {:>8}  TNN PASS  ({:.2?})",
                desc, cond, size_str, elapsed
            ),
            Err(msg) => println!(
                "{:<30} {:>6} {:>8}  TNN FAIL  ({:.2?}) — {}",
                desc, cond, size_str, elapsed, msg
            ),
        }
    }
}
