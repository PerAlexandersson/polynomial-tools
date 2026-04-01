//! Compute B_{mu,d}(x,t) - x for non-rigid strictly decreasing tiles.
//! Extracts B from the transfer matrix via B_n = P_n - sum_{k=1}^{n-1} P_k B_{n-k}.
//!
//! Usage: cargo run --release --bin tiling_B_table

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

/// Compute P_n(t) for n = 0..=max_n using the transfer matrix approach.
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

    // Build transfer matrix
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
    vec_s[0] = vec![1]; // initial state

    for n in 0..=max_n {
        results.push(vec_s[0].clone()); // P_n(t) = vec_s[0]

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

/// Extract B_n(t) from P_n(t) using B_n = P_n - sum_{k=1}^{n-1} P_k * B_{n-k}
fn extract_b(pn: &[Poly]) -> Vec<Poly> {
    let max_n = pn.len() - 1;
    let mut bn: Vec<Poly> = vec![poly_zero(); max_n + 1];
    // B_0 = 0
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

/// Find smallest feasible offset j0
fn find_j0(mu: &[usize], d: usize) -> usize {
    let ell = mu.len() - 1;
    for j in 1..=ell {
        if mu[j] <= d - 1 {
            return j;
        }
    }
    ell + 1 // shouldn't happen if tile is valid
}

/// Check if tile is rigid: mu_ell = 1 and mu_{ell-1} >= d
fn is_rigid(mu: &[usize], d: usize) -> bool {
    let ell = mu.len() - 1;
    if ell == 0 {
        return false;
    }
    mu[ell] == 1 && (ell == 1 || mu[ell - 1] >= d)
}

/// Check claw-free: d <= mu_{ell-1}
fn is_claw_free(mu: &[usize], d: usize) -> bool {
    let ell = mu.len() - 1;
    if ell <= 1 {
        return true;
    }
    d <= mu[ell - 1]
}

/// Format B(x,t) - x as LaTeX
fn format_b_latex(bn: &[Poly]) -> String {
    // bn[n] is the coefficient of x^n in B, which is a polynomial in t
    // We want to print sum of c * t^m * x^n terms
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
            let x_str = match n {
                1 => "x".to_string(),
                _ => format!("x^{{{}}}", n),
            };
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
        terms.join("+")
    }
}

fn main() {
    // Generate non-rigid strictly decreasing tiles with mu_1 >= d
    // Strictly decreasing: mu_1 > mu_2 > ... > mu_ell >= 1
    // Non-rigid: NOT (mu_ell = 1 AND mu_{ell-1} >= d)

    println!("Non-rigid strictly decreasing tiles with mu_1 >= d:\n");
    println!(
        "{:<20} {:<4} {:<4} {:<6} {}",
        "mu", "d", "j0", "claw?", "B(x,t) - x"
    );
    println!("{}", "-".repeat(90));

    // We'll enumerate tiles with ell = 2,3,4 and various d
    let mut entries: Vec<(Vec<usize>, usize)> = Vec::new();

    // ell=2: mu = (h, a, b) with a > b >= 1, a >= d
    for a in 2..=7 {
        for b in 1..a {
            for d in 2..=a.min(6) {
                let mu = vec![a + d, a, b]; // mu_0 >= d, mu_1 = a >= d
                if !is_rigid(&mu, d) {
                    entries.push((mu, d));
                }
            }
        }
    }

    // ell=3: mu = (h, a, b, c) with a > b > c >= 1, a >= d
    for a in 3..=6 {
        for b in 2..a {
            for c in 1..b {
                for d in 2..=a.min(5) {
                    let mu = vec![a + d, a, b, c];
                    if !is_rigid(&mu, d) {
                        entries.push((mu, d));
                    }
                }
            }
        }
    }

    // ell=4: mu = (h, a, b, c, e) with a > b > c > e >= 1, a >= d
    for a in 4..=6 {
        for b in 3..a {
            for c in 2..b {
                for e in 1..c {
                    for d in 2..=a.min(5) {
                        let mu = vec![a + d, a, b, c, e];
                        if !is_rigid(&mu, d) {
                            entries.push((mu, d));
                        }
                    }
                }
            }
        }
    }

    // ell=5: mu = (h, a, b, c, e, f) with a > b > c > e > f >= 1, a >= d
    for a in 5..=6 {
        for b in 4..a {
            for c in 3..b {
                for e in 2..c {
                    for f in 1..e {
                        for d in 2..=a.min(4) {
                            let mu = vec![a + d, a, b, c, e, f];
                            if !is_rigid(&mu, d) {
                                entries.push((mu, d));
                            }
                        }
                    }
                }
            }
        }
    }

    for (mu, d) in &entries {
        let ell = mu.len() - 1;
        let max_width = d * ell + 1;
        let pn = compute_pn(mu, *d, max_width);
        let bn = extract_b(&pn);
        let j0 = find_j0(mu, *d);
        let claw_free = is_claw_free(mu, *d);
        let latex = format_b_latex(&bn);

        // Format mu for display (using h for mu_0)
        let mu_display: Vec<String> = mu
            .iter()
            .enumerate()
            .map(|(i, &v)| {
                if i == 0 {
                    "h".to_string()
                } else {
                    v.to_string()
                }
            })
            .collect();
        let mu_str = format!("({})", mu_display.join(","));

        println!(
            "{:<20} {:<4} {:<4} {:<6} {}",
            mu_str,
            d,
            j0,
            if claw_free { "CF" } else { "claw" },
            latex
        );
    }
}
