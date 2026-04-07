//! Explore the 2x2 interlacing squares
//!
//!   A   = G(c+1,s)
//!   B   = G(c,s)   = A + t Delta
//!   A'  = G(c+1,s+1) = A + t R
//!   B'  = G(c,s+1)   = B + t R = A + t (Delta + R)
//!
//! where R = F_{s+1}^{mu'}.
//!
//! This isolates which relations in the square are always true and which
//! candidate side-lemmas might explain the diagonal/anti-diagonal behavior.

use combpoly::rook_placements::non_nesting_rook_polynomial;
use polynomial_tools::{check_weak_interlacing, format_poly};

fn strip_columns(mu: &[usize], c: usize) -> Vec<usize> {
    mu.iter()
        .map(|&x| x.saturating_sub(c))
        .filter(|&x| x > 0)
        .collect()
}

fn partitions(n: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut buf = Vec::new();
    part_rec(n, n, &mut buf, &mut result);
    result
}

fn part_rec(n: usize, max: usize, buf: &mut Vec<usize>, result: &mut Vec<Vec<usize>>) {
    if n == 0 {
        result.push(buf.clone());
        return;
    }
    for k in (1..=n.min(max)).rev() {
        buf.push(k);
        part_rec(n - k, k, buf, result);
        buf.pop();
    }
}

fn poly_degree(p: &[i64]) -> usize {
    p.iter().rposition(|&c| c != 0).unwrap_or(0)
}

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let n = a.len().max(b.len());
    let mut r = vec![0i64; n];
    for (i, &c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        r[i] += c;
    }
    r
}

fn poly_sub(a: &[i64], b: &[i64]) -> Vec<i64> {
    let n = a.len().max(b.len());
    let mut r = vec![0i64; n];
    for (i, &c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        r[i] -= c;
    }
    r
}

fn poly_mul_t(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; p.len() + 1];
    for (i, &c) in p.iter().enumerate() {
        r[i + 1] = c;
    }
    r
}

fn poly_div_t(p: &[i64]) -> Vec<i64> {
    assert_eq!(p.first().copied().unwrap_or(0), 0, "constant term must be 0");
    if p.len() <= 1 {
        return vec![0];
    }
    p[1..].to_vec()
}

fn interlaces(f: &[i64], g: &[i64]) -> Option<bool> {
    let df = poly_degree(f);
    let dg = poly_degree(g);
    if df == 0 && dg == 0 {
        return Some(true);
    }
    if df > dg + 1 || dg > df + 1 {
        return None;
    }
    if df <= dg {
        check_weak_interlacing(f, g)
    } else {
        check_weak_interlacing(g, f)
    }
}

struct TestResult {
    name: String,
    pass: usize,
    total: usize,
    fails: Vec<String>,
}

impl TestResult {
    fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            pass: 0,
            total: 0,
            fails: Vec::new(),
        }
    }

    fn check(&mut self, ok: bool, fail_msg: impl FnOnce() -> String) {
        self.total += 1;
        if ok {
            self.pass += 1;
        } else if self.fails.len() < 8 {
            self.fails.push(fail_msg());
        }
    }

    fn print(&self) {
        println!("{}: {}/{}", self.name, self.pass, self.total);
        if self.fails.is_empty() {
            println!("  All pass! ✓\n");
        } else {
            for msg in &self.fails {
                println!("  {msg}");
            }
            println!();
        }
    }
}

fn main() {
    let max_n = 14;

    let mut row = TestResult::new("Row: A << B");
    let mut row_next = TestResult::new("Row+: A' << B'");
    let mut col_left = TestResult::new("Column left: A << A'");
    let mut col_right = TestResult::new("Column right: B << B'");
    let mut diagonal = TestResult::new("Diagonal: A << B'");
    let mut anti = TestResult::new("Anti-diagonal: A' << B");

    let mut r_left = TestResult::new("Input: R << A");
    let mut r_right = TestResult::new("Input: R << B");
    let mut delta_left = TestResult::new("Input: Delta << A");
    let mut delta_right = TestResult::new("Input: Delta << B");
    let mut sum_left = TestResult::new("Candidate: Delta + R << A");
    let mut r_delta = TestResult::new("Candidate: R << Delta");
    let mut diff_nonneg = TestResult::new("Candidate: Delta - R has nonnegative coefficients");
    let mut diff_left = TestResult::new("Candidate: Delta - R << A+tR");
    let mut diff_right = TestResult::new("Candidate: Delta - R << B");

    for n in 2..=max_n {
        for mu in partitions(n) {
            if mu.len() < 2 {
                continue;
            }
            let ell = mu.len();
            let m = mu[ell - 1];
            let mp = &mu[..ell - 1];

            let f: Vec<Vec<i64>> = (0..=mu[0])
                .map(|c| non_nesting_rook_polynomial(&strip_columns(mp, c)))
                .collect();

            for c in 0..m {
                for s in (c + 1)..m {
                    let mut a = f[c + 1].clone();
                    for term in f.iter().take(s + 1).skip(c + 2) {
                        a = poly_add(&a, &poly_mul_t(term));
                    }

                    let mut b = f[c].clone();
                    for term in f.iter().take(s + 1).skip(c + 1) {
                        b = poly_add(&b, &poly_mul_t(term));
                    }

                    let r = f[s + 1].clone();
                    let t_r = poly_mul_t(&r);
                    let a_next = poly_add(&a, &t_r);
                    let b_next = poly_add(&b, &t_r);

                    let diff = poly_sub(&b, &a);
                    let delta = poly_div_t(&diff);
                    let delta_plus_r = poly_add(&delta, &r);
                    let delta_minus_r = poly_sub(&delta, &r);

                    row.check(interlaces(&a, &b) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: {} << {}",
                            mu,
                            c,
                            s,
                            format_poly(&a),
                            format_poly(&b)
                        )
                    });
                    row_next.check(interlaces(&a_next, &b_next) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: {} << {}",
                            mu,
                            c,
                            s,
                            format_poly(&a_next),
                            format_poly(&b_next)
                        )
                    });
                    col_left.check(interlaces(&a, &a_next) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: {} << {}",
                            mu,
                            c,
                            s,
                            format_poly(&a),
                            format_poly(&a_next)
                        )
                    });
                    col_right.check(interlaces(&b, &b_next) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: {} << {}",
                            mu,
                            c,
                            s,
                            format_poly(&b),
                            format_poly(&b_next)
                        )
                    });
                    diagonal.check(interlaces(&a, &b_next) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: {} << {}",
                            mu,
                            c,
                            s,
                            format_poly(&a),
                            format_poly(&b_next)
                        )
                    });
                    anti.check(interlaces(&a_next, &b) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: {} << {}",
                            mu,
                            c,
                            s,
                            format_poly(&a_next),
                            format_poly(&b)
                        )
                    });

                    r_left.check(interlaces(&r, &a) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: R={} << A={:?}",
                            mu,
                            c,
                            s,
                            format_poly(&r),
                            interlaces(&r, &a)
                        )
                    });
                    r_right.check(interlaces(&r, &b) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: R={} << B={:?}",
                            mu,
                            c,
                            s,
                            format_poly(&r),
                            interlaces(&r, &b)
                        )
                    });
                    delta_left.check(interlaces(&delta, &a) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: Delta={} << A={}",
                            mu,
                            c,
                            s,
                            format_poly(&delta),
                            format_poly(&a)
                        )
                    });
                    delta_right.check(interlaces(&delta, &b) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: Delta={} << B={}",
                            mu,
                            c,
                            s,
                            format_poly(&delta),
                            format_poly(&b)
                        )
                    });
                    sum_left.check(interlaces(&delta_plus_r, &a) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: Delta+R={} << A={:?}",
                            mu,
                            c,
                            s,
                            format_poly(&delta_plus_r),
                            interlaces(&delta_plus_r, &a)
                        )
                    });
                    r_delta.check(interlaces(&r, &delta) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: R={} << Delta={:?}",
                            mu,
                            c,
                            s,
                            format_poly(&r),
                            interlaces(&r, &delta)
                        )
                    });
                    diff_nonneg.check(delta_minus_r.iter().all(|&v| v >= 0), || {
                        format!(
                            "mu={:?} c={} s={}: Delta-R={}",
                            mu,
                            c,
                            s,
                            format_poly(&delta_minus_r)
                        )
                    });
                    diff_left.check(interlaces(&delta_minus_r, &a_next) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: Delta-R={} << A'={:?}",
                            mu,
                            c,
                            s,
                            format_poly(&delta_minus_r),
                            interlaces(&delta_minus_r, &a_next)
                        )
                    });
                    diff_right.check(interlaces(&delta_minus_r, &b) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: Delta-R={} << B={:?}",
                            mu,
                            c,
                            s,
                            format_poly(&delta_minus_r),
                            interlaces(&delta_minus_r, &b)
                        )
                    });
                }
            }
        }
    }

    row.print();
    row_next.print();
    col_left.print();
    col_right.print();
    diagonal.print();
    anti.print();

    r_left.print();
    r_right.print();
    delta_left.print();
    delta_right.print();
    sum_left.print();
    r_delta.print();
    diff_nonneg.print();
    diff_left.print();
    diff_right.print();
}
