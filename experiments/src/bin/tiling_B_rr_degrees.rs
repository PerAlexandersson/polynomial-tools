//! Check effective degrees of fault-free row polynomials B_n(t)
//! after factoring out t^k. How many have effective degree >= 3?
//!
//! Usage: cargo run --release --bin tiling_B_rr_degrees [max_n] [max_parts] [max_val] [max_d]

use polynomial_tools::is_real_rooted;

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

/// Effective degree: degree minus valuation (lowest nonzero power).
fn effective_degree(p: &Poly) -> usize {
    if poly_is_zero(p) {
        return 0;
    }
    let lo = p.iter().position(|&c| c != 0).unwrap_or(0);
    let hi = p.len() - 1 - p.iter().rev().position(|&c| c != 0).unwrap_or(0);
    hi - lo
}

fn partitions(k: usize, max_val: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    fn helper(k: usize, max_val: usize, current: &mut Vec<usize>, result: &mut Vec<Vec<usize>>) {
        if k == 0 {
            result.push(current.clone());
            return;
        }
        let upper = if current.is_empty() {
            max_val
        } else {
            *current.last().unwrap()
        };
        for v in (1..=upper).rev() {
            current.push(v);
            helper(k - 1, max_val, current, result);
            current.pop();
        }
    }
    helper(k, max_val, &mut Vec::new(), &mut result);
    result
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let max_n: usize = args.get(0).and_then(|s| s.parse().ok()).unwrap_or(50);
    let max_parts: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(5);
    let max_val: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(6);
    let max_d: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(7);

    let mut total_polys = 0usize;
    let mut by_eff_deg: Vec<usize> = vec![0; 30];
    let mut max_eff_deg_seen = 0usize;
    let mut rr_fail_count = 0usize;
    let mut nontrivial_rr_count = 0usize; // eff_deg >= 3 and RR
    let mut nontrivial_total = 0usize;
    let mut total_cases = 0usize;

    for num_parts in 2..=max_parts {
        for mu in partitions(num_parts, max_val) {
            let ell = mu.len() - 1;
            for d in 1..=max_d.min(mu[0]) {
                let ns = (d + 1).pow(ell as u32);
                if ns > 50_000 {
                    continue;
                }

                total_cases += 1;
                let pn = compute_pn(&mu, d, max_n);
                let bn = extract_b(&pn);

                for p in &bn {
                    if poly_is_zero(p) {
                        continue;
                    }
                    total_polys += 1;
                    let ed = effective_degree(p);
                    if ed < by_eff_deg.len() {
                        by_eff_deg[ed] += 1;
                    }
                    if ed > max_eff_deg_seen {
                        max_eff_deg_seen = ed;
                    }

                    if ed >= 3 {
                        nontrivial_total += 1;
                        let coeffs64: Vec<i64> = p.iter().map(|&c| c as i64).collect();
                        if is_real_rooted(&coeffs64) {
                            nontrivial_rr_count += 1;
                        } else {
                            rr_fail_count += 1;
                            // Print the failure
                            let lo = p.iter().position(|&c| c != 0).unwrap_or(0);
                            let hi = p.len() - 1;
                            let trimmed: Vec<i128> = p[lo..=hi].to_vec();
                            eprintln!(
                                "RR FAIL: mu={:?} d={} eff_deg={} poly={:?}",
                                mu, d, ed, trimmed
                            );
                        }
                    }
                }
            }
        }
    }

    println!("Scanned {} (mu,d) cases, max_n={}", total_cases, max_n);
    println!("Total nonzero B_n polynomials: {}", total_polys);
    println!("\nEffective degree distribution:");
    for deg in 0..=max_eff_deg_seen {
        if by_eff_deg[deg] > 0 {
            println!("  eff_deg={}: {} polys", deg, by_eff_deg[deg]);
        }
    }
    println!("\nNontrivial (eff_deg >= 3): {} total", nontrivial_total);
    println!("  Real-rooted: {}", nontrivial_rr_count);
    println!("  NOT real-rooted: {}", rr_fail_count);
    println!("\nMax effective degree seen: {}", max_eff_deg_seen);
}
