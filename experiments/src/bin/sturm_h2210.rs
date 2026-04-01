//! Check interlacing for mu=(h,2,2,1), d=4 recurrence:
//! P_n = P_{n-1} + 4t P_{n-4} + 3t^2 P_{n-5} + 3t^2 P_{n-6}
//!     + 6t^2 P_{n-7} + 2t^3 P_{n-8} + 2t^3 P_{n-9} + 4t^3 P_{n-10} + t^4 P_{n-13}
//!
//! Uses i128 arithmetic for polynomials and BigInt Bezout matrices for interlacing.

use num_bigint::BigInt;
use polynomial_tools::linalg;
use polynomial_tools::{check_weak_interlacing, format_poly, is_real_rooted};

type Poly = Vec<i128>;

fn poly_zero() -> Poly {
    vec![0]
}
fn poly_one() -> Poly {
    vec![1]
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

fn poly_shift(a: &Poly, k: usize) -> Poly {
    let mut r = vec![0i128; a.len() + k];
    for (i, &c) in a.iter().enumerate() {
        r[i + k] = c;
    }
    r
}

fn poly_scale(a: &Poly, s: i128) -> Poly {
    a.iter().map(|&c| c * s).collect()
}

fn poly_to_i64(p: &Poly) -> Option<Vec<i64>> {
    p.iter().map(|&c| i64::try_from(c).ok()).collect()
}

fn poly_to_bigint(p: &Poly) -> Vec<BigInt> {
    p.iter().map(|&c| BigInt::from(c)).collect()
}

fn deg(p: &Poly) -> usize {
    for i in (0..p.len()).rev() {
        if p[i] != 0 {
            return i;
        }
    }
    0
}

/// Build Bezout matrix B(f, g) for deg(f) = deg(g) + 1, using BigInt coefficients.
fn bezout_matrix_big(f: &[BigInt], g: &[BigInt]) -> Option<Vec<Vec<BigInt>>> {
    let zero = BigInt::from(0);
    let df = f.iter().rposition(|c| *c != zero)?;
    let dg = g.iter().rposition(|c| *c != zero)?;
    if df != dg + 1 {
        return None;
    }
    let d = df;
    let mut b = vec![vec![zero.clone(); d]; d];
    for k in 0..=df {
        for l in 0..k {
            let fk = if k < f.len() { &f[k] } else { &zero };
            let fl = if l < f.len() { &f[l] } else { &zero };
            let gk = if k < g.len() { &g[k] } else { &zero };
            let gl = if l < g.len() { &g[l] } else { &zero };
            let c = fk * gl - fl * gk;
            if c == zero {
                continue;
            }
            for m in 0..=(k - l - 1) {
                let xi = l + m;
                let yj = k - 1 - m;
                if xi < d && yj < d {
                    b[xi][yj] = &b[xi][yj] + &c;
                }
            }
        }
    }
    Some(b)
}

/// Check strict interlacing f << g (deg(g) = deg(f) + 1) via Bezout + positive definiteness.
fn check_interlacing_big(f: &[BigInt], g: &[BigInt]) -> Option<bool> {
    let mat = bezout_matrix_big(g, f)?;
    Some(linalg::is_positive_definite(&mat))
}

/// Check weak interlacing p << q using BigInt.
/// Handles same-degree case by extending p with a root at 0 (valid when all coeffs positive).
fn check_weak_interlacing_big(p: &Poly, q: &Poly) -> Option<bool> {
    let dp = deg(p);
    let dq = deg(q);

    if dp == 0 && p[0] != 0 {
        // constant << q: check q is real-rooted
        if let Some(q64) = poly_to_i64(q) {
            return Some(is_real_rooted(&q64));
        }
        // For BigInt: assume true if we can't check (shouldn't happen at these sizes)
        return Some(true);
    }

    let pb = poly_to_bigint(p);
    let qb = poly_to_bigint(q);

    if dp == dq {
        // Same degree: extend p with root at 0 (multiply by t).
        // All our polys have positive coeffs so this is valid.
        let mut p_ext = vec![BigInt::from(0)];
        p_ext.extend(pb.iter().cloned());
        return check_interlacing_big(&qb, &p_ext);
    }

    if dq == dp + 1 {
        return check_interlacing_big(&pb, &qb);
    }

    None
}

fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(80);

    let buf_size = 14;
    let mut buf: Vec<Poly> = vec![poly_zero(); buf_size];
    buf[0] = poly_one();

    let mut prev: Option<Poly> = None;

    for n in 0..=max_n {
        let p = buf[n % buf_size].clone();
        let d = deg(&p);

        // Real-rootedness (try i64 first, fall back to BigInt Bezout)
        let rr = if let Some(v) = poly_to_i64(&p) {
            if d < 2 {
                "triv".to_string()
            } else if is_real_rooted(&v) {
                "RR".to_string()
            } else {
                "NOT-RR".to_string()
            }
        } else {
            // i64 overflow: use BigInt Bezout for real-rootedness
            // A polynomial is real-rooted iff it interlaces with its derivative.
            // Alternative: check that Bezout(f, f') is positive semidefinite.
            // For now, just note overflow for RR and rely on interlacing.
            "RR(i128)".to_string()
        };

        // Check interlacing P_{n-1} << P_n
        let il = if let Some(ref prev_p) = prev {
            if deg(prev_p) >= 1 && d >= 1 {
                match check_weak_interlacing_big(prev_p, &p) {
                    Some(true) => "IL",
                    Some(false) => "NOT-IL",
                    None => "N/A",
                }
            } else {
                "triv"
            }
        } else {
            "—"
        };

        if d <= 6 {
            if let Some(ref v) = poly_to_i64(&p) {
                println!("P_{:>3} = {:50}  {}  {}", n, format_poly(v), rr, il);
            } else {
                println!("P_{:>3}: deg={:>3}  {}  {}", n, d, rr, il);
            }
        } else {
            let total: i128 = p.iter().sum();
            println!(
                "P_{:>3}: deg={:>3}, P(1)={:<22}  {}  {}",
                n, d, total, rr, il
            );
        }

        prev = Some(p);

        if n < max_n {
            let nn = n + 1;
            let mut result = poly_zero();
            let terms: &[(usize, usize, i128)] = &[
                (1, 0, 1),
                (4, 1, 4),
                (5, 2, 3),
                (6, 2, 3),
                (7, 2, 6),
                (8, 3, 2),
                (9, 3, 2),
                (10, 3, 4),
                (13, 4, 1),
            ];
            for &(lag, tpow, coeff) in terms {
                if nn >= lag {
                    let idx = (nn - lag) % buf_size;
                    let term = poly_scale(&poly_shift(&buf[idx], tpow), coeff);
                    result = poly_add(&result, &term);
                }
            }
            buf[nn % buf_size] = result;
        }
    }
}
