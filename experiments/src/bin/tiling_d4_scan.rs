//! Targeted scan for d=4 real-rooted families.
//! Usage: cargo run --release --bin tiling_d4_scan

use polynomial_tools::is_real_rooted;

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

fn poly_gcd(p: &Poly) -> i128 {
    fn gcd(a: i128, b: i128) -> i128 {
        let (mut a, mut b) = (a.abs(), b.abs());
        while b != 0 {
            let t = b;
            b = a % b;
            a = t;
        }
        a
    }
    let mut g = 0i128;
    for &c in p {
        g = gcd(g, c);
    }
    if g == 0 {
        1
    } else {
        g
    }
}

fn check_rr(p: &Poly) -> Option<bool> {
    if p.len() <= 2 {
        return Some(true);
    }
    if let Some(p64) = poly_to_i64(p) {
        return Some(is_real_rooted(&p64));
    }
    let g = poly_gcd(p);
    if g > 1 {
        let scaled: Poly = p.iter().map(|&c| c / g).collect();
        if let Some(s64) = poly_to_i64(&scaled) {
            return Some(is_real_rooted(&s64));
        }
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

fn test_shape(mu: &[usize], d: usize, max_n: usize) -> (Option<usize>, usize) {
    let ell = mu.len() - 1;
    if ell == 0 {
        return (None, max_n);
    }

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

    let mut last_checked = 0;
    for n in 0..=max_n {
        let p = &vec_s[0];
        match check_rr(p) {
            Some(true) => {
                last_checked = n;
            }
            Some(false) => {
                return (Some(n), n);
            }
            None => {
                return (None, last_checked);
            }
        }

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
    (None, max_n)
}

fn main() {
    let d = 4;
    let max_n = 60;

    println!(
        "Scanning Ferrers tiles with d={}, checking RR up to n={}",
        d, max_n
    );
    println!(
        "{:<25} {:>3} {:>3} {:>5} {:<25}",
        "mu", "m1", "m2", "CF", "result"
    );
    println!("{}", "-".repeat(70));

    // 3-column shapes (ell=2)
    println!("\n=== 3-column tiles ===");
    for a in 1..=8 {
        for b in 1..=a {
            for c in 1..=b {
                let mu = vec![a, b, c];
                let cf = d <= c;
                let (fail, last) = test_shape(&mu, d, max_n);
                let status = match fail {
                    Some(n) => format!("FAILS n={}", n),
                    None if last >= max_n => format!("RR (n<={})", max_n),
                    None => format!("RR (n<={}) ovfl", last),
                };
                // Only print non-CF cases (CF cases are boring)
                if !cf {
                    let marker = if status.starts_with("RR") && !status.contains("FAILS") {
                        " <<<"
                    } else {
                        ""
                    };
                    println!(
                        "{:<25} {:>3} {:>3} {:>5} {:<25}{}",
                        format!("{:?}", mu),
                        b,
                        c,
                        if cf { "yes" } else { "no" },
                        status,
                        marker
                    );
                }
            }
        }
    }

    // 4-column shapes (ell=3), limited range
    println!("\n=== 4-column tiles ===");
    for a in 2..=6 {
        for b in 1..=a.min(5) {
            for c in 1..=b.min(4) {
                for e in 1..=c.min(3) {
                    let mu = vec![a, b, c, e];
                    let cf = d <= e;
                    let (fail, last) = test_shape(&mu, d, max_n);
                    let status = match fail {
                        Some(n) => format!("FAILS n={}", n),
                        None if last >= max_n => format!("RR (n<={})", max_n),
                        None => format!("RR (n<={}) ovfl", last),
                    };
                    if !cf {
                        let marker = if status.starts_with("RR") && !status.contains("FAILS") {
                            " <<<"
                        } else {
                            ""
                        };
                        println!(
                            "{:<25} {:>3} {:>3} {:>5} {:<25}{}",
                            format!("{:?}", mu),
                            b,
                            c,
                            if cf { "yes" } else { "no" },
                            status,
                            marker
                        );
                    }
                }
            }
        }
    }

    // Summary
    println!("\n{}", "=".repeat(70));
    println!("Pattern check: does mu_1 > mu_2 predict RR?");
    println!(
        "{:<12} {:>8} {:>8} {:>8}",
        "condition", "RR", "FAILS", "overflow"
    );

    let mut counts = std::collections::BTreeMap::new();
    // Re-run to collect stats
    for ncols in 3..=4 {
        let ranges: Vec<Vec<usize>> = if ncols == 3 {
            let mut v = vec![];
            for a in 1..=8 {
                for b in 1..=a {
                    for c in 1..=b {
                        v.push(vec![a, b, c]);
                    }
                }
            }
            v
        } else {
            let mut v = vec![];
            for a in 2..=6 {
                for b in 1..=a.min(5) {
                    for c in 1..=b.min(4) {
                        for e in 1..=c.min(3) {
                            v.push(vec![a, b, c, e]);
                        }
                    }
                }
            }
            v
        };
        for mu in &ranges {
            let cf = d <= *mu.last().unwrap();
            if cf {
                continue;
            }
            let m1 = mu[1];
            let m2 = if mu.len() > 2 { mu[2] } else { 0 };
            let cond = if m1 > m2 { "m1>m2" } else { "m1=m2" };
            let (fail, last) = test_shape(mu, d, max_n);
            let cat = match fail {
                Some(_) => "FAILS",
                None if last >= max_n => "RR",
                None => "overflow",
            };
            let entry = counts
                .entry((cond.to_string(), cat.to_string()))
                .or_insert(0u32);
            *entry += 1;
        }
    }
    for cond in &["m1>m2", "m1=m2"] {
        let rr = counts
            .get(&(cond.to_string(), "RR".to_string()))
            .unwrap_or(&0);
        let fail = counts
            .get(&(cond.to_string(), "FAILS".to_string()))
            .unwrap_or(&0);
        let ovfl = counts
            .get(&(cond.to_string(), "overflow".to_string()))
            .unwrap_or(&0);
        println!("{:<12} {:>8} {:>8} {:>8}", cond, rr, fail, ovfl);
    }
}
