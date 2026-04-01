//! Compute B_{mu,d}(x,t) - x for mu = (h,h,1,1) with d = h.
//! Usage: cargo run --release --bin tiling_B_hh11

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
            let term = poly_mul(&pn[k], &bn[n - k]);
            sum = poly_add(&sum, &term);
        }
        bn[n] = poly_sub(&pn[n], &sum);
    }
    bn
}

fn format_b_latex(bn: &[Poly]) -> String {
    let mut terms: Vec<String> = Vec::new();
    for n in 2..bn.len() {
        let p = &bn[n];
        if poly_is_zero(p) {
            continue;
        }
        for (m, &c) in p.iter().enumerate() {
            if c == 0 {
                continue;
            }
            let coeff_str = if c == 1 {
                String::new()
            } else {
                format!("{}", c)
            };
            let t_str = match m {
                0 => String::new(),
                1 => "t".to_string(),
                _ => format!("t^{}", m),
            };
            let x_str = format!("x^{{{}}}", n);
            let term = if coeff_str.is_empty() && t_str.is_empty() {
                x_str
            } else if coeff_str.is_empty() {
                format!("{}{}", t_str, x_str)
            } else if t_str.is_empty() {
                format!("{}{}", coeff_str, x_str)
            } else {
                format!("{}{}{}", coeff_str, t_str, x_str)
            };
            terms.push(term);
        }
    }
    if terms.is_empty() {
        "0".to_string()
    } else {
        terms.join(" + ")
    }
}

fn main() {
    println!("B_{{mu,d}}(x,t) - x for mu = (h,h,1,1), d = h\n");

    for h in 2..=8 {
        let mu = vec![h, h, 1, 1];
        let d = h;
        let ell = mu.len() - 1; // = 3
        let max_width = d * ell + 1;
        let pn = compute_pn(&mu, d, max_width);
        let bn = extract_b(&pn);
        let latex = format_b_latex(&bn);

        println!("h = {}: B - x = {}", h, latex);

        // Also print grouped by t-power for clarity
        let mut by_t: std::collections::BTreeMap<usize, Vec<(usize, i128)>> =
            std::collections::BTreeMap::new();
        for n in 2..bn.len() {
            let p = &bn[n];
            for (m, &c) in p.iter().enumerate() {
                if c != 0 {
                    by_t.entry(m).or_default().push((n, c));
                }
            }
        }
        for (m, terms) in &by_t {
            let parts: Vec<String> = terms
                .iter()
                .map(|(n, c)| {
                    if *c == 1 {
                        format!("x^{}", n)
                    } else {
                        format!("{}x^{}", c, n)
                    }
                })
                .collect();
            println!("  [t^{}]: {}", m, parts.join(" + "));
        }
        println!();
    }
}
