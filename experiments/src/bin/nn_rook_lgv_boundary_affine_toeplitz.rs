//! Symbolic Toeplitz-minor diagnostic for affine same-degree boundary packets.
//!
//! The hOrient-free Lean endpoint asks for Toeplitz total nonnegativity of
//! the coefficient sequence of
//!
//!   A(x) + (s*x + t) Q(x)
//!
//! for all positive `s,t`, with `Q` equal to the upper or lower boundary
//! packet polynomial.  This binary checks selected Toeplitz minors of that
//! sequence coefficientwise as polynomials in the two variables `s,t`.
//! Coefficientwise nonnegativity is stronger than positivity for positive
//! parameter values, and failure would identify a bad candidate packet-row
//! certificate target.

use std::collections::{BTreeMap, HashMap};
use std::env;

use experiments::nn_rook_utils::{
    add, degree, fixed_row_count_lgv, monomial, partitions, subset_size, trim,
};
use polynomial_tools::format_poly;

type BiPoly = BTreeMap<(usize, usize), i128>;

fn fixed_count_cached(
    cache: &mut HashMap<(usize, u64), i64>,
    eta: &[usize],
    strip: usize,
    mask: u64,
) -> i64 {
    if let Some(&value) = cache.get(&(strip, mask)) {
        return value;
    }
    let value = fixed_row_count_lgv(eta, strip, mask);
    cache.insert((strip, mask), value);
    value
}

fn rho_count(eta: &[usize], mask: u64, c: usize) -> i64 {
    let highest = (0..eta.len()).filter(|&r| (mask & (1u64 << r)) != 0).max();
    let first_forced = highest.map_or(0, |r| r + 1);
    eta.iter()
        .enumerate()
        .skip(first_forced)
        .filter(|&(_, &width)| width > c)
        .count() as i64
}

fn prefix_packet(
    cache: &mut HashMap<(usize, u64), i64>,
    eta: &[usize],
    h: usize,
    c: usize,
    p: usize,
) -> (Vec<i64>, Vec<i64>, Vec<i64>) {
    let mut a_poly = vec![0];
    let mut u_poly = vec![0];
    let mut l_poly = vec![0];

    for mask in 0..(1u64 << (h + 1)) {
        let k = subset_size(mask);
        let d_c1 = fixed_count_cached(cache, eta, c + 1, mask);
        let d_p1 = fixed_count_cached(cache, eta, p + 1, mask);
        let tail_sum: i64 = ((c + 2)..=(p + 1))
            .map(|strip| fixed_count_cached(cache, eta, strip, mask))
            .sum();
        let rho = rho_count(eta, mask, c);
        let window = d_c1 - d_p1;

        a_poly = add(&a_poly, &monomial(k, d_c1));
        a_poly = add(&a_poly, &monomial(k + 1, tail_sum));
        u_poly = add(&u_poly, &monomial(k, rho * d_c1 + d_c1 + d_p1));
        l_poly = add(&l_poly, &monomial(k, rho * d_c1 + window));
    }

    (trim(a_poly), trim(u_poly), trim(l_poly))
}

fn prefix_strip_capacity(
    cache: &mut HashMap<(usize, u64), i64>,
    eta: &[usize],
    h: usize,
    strip: usize,
) -> usize {
    let mut cap = 0usize;
    for mask in 0..(1u64 << (h + 1)) {
        if fixed_count_cached(cache, eta, strip, mask) != 0 {
            cap = cap.max(subset_size(mask));
        }
    }
    cap
}

fn bi_zero() -> BiPoly {
    BTreeMap::new()
}

fn bi_const(c: i128) -> BiPoly {
    if c == 0 {
        return bi_zero();
    }
    BTreeMap::from([((0, 0), c)])
}

fn bi_monom(s_exp: usize, t_exp: usize, c: i128) -> BiPoly {
    if c == 0 {
        return bi_zero();
    }
    BTreeMap::from([((s_exp, t_exp), c)])
}

fn bi_add(a: &BiPoly, b: &BiPoly) -> BiPoly {
    let mut out = a.clone();
    for (&key, &value) in b {
        let entry = out.entry(key).or_insert(0);
        *entry += value;
        if *entry == 0 {
            out.remove(&key);
        }
    }
    out
}

fn bi_sub(a: &BiPoly, b: &BiPoly) -> BiPoly {
    let mut out = a.clone();
    for (&key, &value) in b {
        let entry = out.entry(key).or_insert(0);
        *entry -= value;
        if *entry == 0 {
            out.remove(&key);
        }
    }
    out
}

fn bi_mul(a: &BiPoly, b: &BiPoly) -> BiPoly {
    let mut out = bi_zero();
    for (&(sa, ta), &ca) in a {
        for (&(sb, tb), &cb) in b {
            let key = (sa + sb, ta + tb);
            *out.entry(key).or_insert(0) += ca * cb;
        }
    }
    out.retain(|_, c| *c != 0);
    out
}

fn det_bipoly(mat: &[Vec<BiPoly>]) -> BiPoly {
    let n = mat.len();
    if n == 0 {
        return bi_const(1);
    }
    if n == 1 {
        return mat[0][0].clone();
    }

    let mut total = bi_zero();
    for col in 0..n {
        if mat[0][col].is_empty() {
            continue;
        }
        let mut submat = Vec::with_capacity(n - 1);
        for row in mat.iter().skip(1) {
            let mut subrow = Vec::with_capacity(n - 1);
            for (j, entry) in row.iter().enumerate() {
                if j != col {
                    subrow.push(entry.clone());
                }
            }
            submat.push(subrow);
        }
        let term = bi_mul(&mat[0][col], &det_bipoly(&submat));
        total = if col % 2 == 0 {
            bi_add(&total, &term)
        } else {
            bi_sub(&total, &term)
        };
    }
    total
}

fn coeff_at(poly: &[i64], k: usize) -> i128 {
    poly.get(k).copied().unwrap_or(0) as i128
}

fn affine_sequence(q: &[i64], a: &[i64]) -> Vec<BiPoly> {
    let len = (q.len() + 1).max(a.len());
    (0..len)
        .map(|k| {
            let mut out = bi_const(coeff_at(a, k));
            if k > 0 {
                out = bi_add(&out, &bi_monom(1, 0, coeff_at(q, k - 1)));
            }
            bi_add(&out, &bi_monom(0, 1, coeff_at(q, k)))
        })
        .collect()
}

fn toeplitz_entry(seq: &[BiPoly], row: usize, col: usize) -> BiPoly {
    if col < row {
        return bi_zero();
    }
    seq.get(col - row).cloned().unwrap_or_default()
}

fn solid_minor(seq: &[BiPoly], row0: usize, col0: usize, size: usize) -> BiPoly {
    let mut mat = vec![vec![bi_zero(); size]; size];
    for i in 0..size {
        for j in 0..size {
            mat[i][j] = toeplitz_entry(seq, row0 + i, col0 + j);
        }
    }
    det_bipoly(&mat)
}

fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    fn rec(n: usize, k: usize, start: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if cur.len() == k {
            out.push(cur.clone());
            return;
        }
        let need = k - cur.len();
        for value in start..=n - need {
            cur.push(value);
            rec(n, k, value + 1, cur, out);
            cur.pop();
        }
    }

    let mut out = Vec::new();
    rec(n, k, 0, &mut Vec::new(), &mut out);
    out
}

fn arbitrary_minor(seq: &[BiPoly], rows: &[usize], cols: &[usize]) -> BiPoly {
    let size = rows.len();
    let mut mat = vec![vec![bi_zero(); size]; size];
    for (i, &row) in rows.iter().enumerate() {
        for (j, &col) in cols.iter().enumerate() {
            mat[i][j] = toeplitz_entry(seq, row, col);
        }
    }
    det_bipoly(&mat)
}

fn first_negative_coeff(p: &BiPoly) -> Option<((usize, usize), i128)> {
    p.iter()
        .find_map(|(&key, &value)| (value < 0).then_some((key, value)))
}

fn format_bipoly(p: &BiPoly) -> String {
    if p.is_empty() {
        return "0".to_string();
    }
    p.iter()
        .map(|(&(s, t), &c)| {
            let mut out = c.to_string();
            if s > 0 {
                out.push_str(&format!("*s^{s}"));
            }
            if t > 0 {
                out.push_str(&format!("*t^{t}"));
            }
            out
        })
        .collect::<Vec<_>>()
        .join(" + ")
}

struct FailCtx<'a> {
    label: &'a str,
    eta: &'a [usize],
    h: usize,
    c: usize,
    p: usize,
    size: usize,
    rows: Vec<usize>,
    cols: Vec<usize>,
    a: &'a [i64],
    q: &'a [i64],
}

fn fail(ctx: FailCtx<'_>, minor: &BiPoly, bad: ((usize, usize), i128)) -> ! {
    println!("FAIL affine boundary Toeplitz minor:");
    println!("  label={} eta={:?}", ctx.label, ctx.eta);
    println!("  h={} c={} p={}", ctx.h, ctx.c, ctx.p);
    println!("  deg A={} deg Q={}", degree(ctx.a), degree(ctx.q));
    println!(
        "  size={} rows={:?} cols={:?}",
        ctx.size, ctx.rows, ctx.cols
    );
    println!("  negative coeff s^{} t^{} = {}", bad.0 .0, bad.0 .1, bad.1);
    println!("  A={}", format_poly(ctx.a));
    println!("  Q={}", format_poly(ctx.q));
    println!("  minor={}", format_bipoly(minor));
    std::process::exit(1);
}

fn main() {
    let max_n = env::args()
        .nth(1)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(14);
    let max_minor_size = env::args()
        .nth(2)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(4);
    let matrix_pad = env::args()
        .nth(3)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(1);
    let mode = env::args().nth(4).unwrap_or_else(|| "solid".to_string());

    println!("=== LGV boundary affine symbolic Toeplitz scan ===");
    println!(
        "Checking Ferrers partitions of size <= {max_n}, {mode} minors through size {max_minor_size}.\n"
    );

    let mut boundary_packets = 0usize;
    let mut minors_checked = 0usize;

    for n in 1..=max_n {
        for eta in partitions(n) {
            if eta.len() >= 63 {
                continue;
            }
            let min_width = eta[eta.len() - 1];
            if min_width < 2 {
                continue;
            }
            let mut cache = HashMap::new();

            for h in 0..eta.len() {
                for c in 0..min_width {
                    let cap_c1 = prefix_strip_capacity(&mut cache, &eta, h, c + 1);
                    let cap_c2 = prefix_strip_capacity(&mut cache, &eta, h, c + 2);
                    if cap_c2 >= cap_c1 {
                        continue;
                    }
                    for p in (c + 1)..min_width {
                        let (a, u, l) = prefix_packet(&mut cache, &eta, h, c, p);
                        if degree(&a) != degree(&u) || degree(&a) != degree(&l) {
                            continue;
                        }
                        boundary_packets += 1;

                        for (label, q) in [("U", &u), ("L", &l)] {
                            let seq = affine_sequence(q, &a);
                            let matrix_n = seq.len() + matrix_pad;
                            if mode == "all" {
                                for size in 1..=max_minor_size.min(matrix_n) {
                                    let choices = combinations(matrix_n, size);
                                    for rows in &choices {
                                        for cols in &choices {
                                            minors_checked += 1;
                                            let minor = arbitrary_minor(&seq, rows, cols);
                                            if let Some(bad) = first_negative_coeff(&minor) {
                                                fail(
                                                    FailCtx {
                                                        label,
                                                        eta: &eta,
                                                        h,
                                                        c,
                                                        p,
                                                        size,
                                                        rows: rows.clone(),
                                                        cols: cols.clone(),
                                                        a: &a,
                                                        q,
                                                    },
                                                    &minor,
                                                    bad,
                                                );
                                            }
                                        }
                                    }
                                }
                            } else {
                                for size in 1..=max_minor_size.min(matrix_n) {
                                    for row0 in 0..=(matrix_n - size) {
                                        for col0 in 0..=(matrix_n - size) {
                                            minors_checked += 1;
                                            let minor = solid_minor(&seq, row0, col0, size);
                                            if let Some(bad) = first_negative_coeff(&minor) {
                                                fail(
                                                    FailCtx {
                                                        label,
                                                        eta: &eta,
                                                        h,
                                                        c,
                                                        p,
                                                        size,
                                                        rows: (row0..row0 + size).collect(),
                                                        cols: (col0..col0 + size).collect(),
                                                        a: &a,
                                                        q,
                                                    },
                                                    &minor,
                                                    bad,
                                                );
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    println!("boundary packets checked: {boundary_packets}");
    println!("Toeplitz minors checked: {minors_checked}");
    println!(
        "All checked affine boundary Toeplitz minors were coefficientwise nonnegative in s,t."
    );
}
