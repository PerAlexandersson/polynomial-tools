//! Focused evidence for the two-state / grid recurrence package.
//!
//! For a board mu with last row width m and mu' = all-but-last:
//!   f_c = R_{mu'_<c>}
//!   G(c,s) = f_c + t * sum_{j=c+1}^s f_j
//!   Delta_c(mu) = (R_{mu_<c>} - R_{mu_<c+1>}) / t
//!               = f_{c+1} + (f_c - f_{c+1}) / t
//!
//! Then:
//!   R_{mu_<c>} = G(c,m)
//!   G(c,s) = G(c+1,s) + t * Delta_c(mu)
//!   G(c,s+1) = G(c,s) + t * f_{s+1}
//!
//! This binary checks the interlacing package suggested by these recurrences.

use std::env;

use combpoly::rook_placements::non_nesting_rook_polynomial;
use polynomial_tools::{check_weak_interlacing, format_poly, is_real_rooted};

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
    assert_eq!(p.first().copied().unwrap_or(0), 0);
    if p.len() <= 1 {
        return vec![0];
    }
    p[1..].to_vec()
}

fn interlaces(f: &[i64], g: &[i64]) -> Option<bool> {
    let df = poly_degree(f);
    let dg = poly_degree(g);
    if df == 0 && f.first().copied().unwrap_or(0) == 0 {
        return Some(is_real_rooted(g));
    }
    if dg == 0 && g.first().copied().unwrap_or(0) == 0 {
        return Some(df == 0 && f.first().copied().unwrap_or(0) == 0);
    }
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
    name: &'static str,
    pass: usize,
    total: usize,
    fails: Vec<String>,
}

impl TestResult {
    fn new(name: &'static str) -> Self {
        Self {
            name,
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
            for fail in &self.fails {
                println!("  {fail}");
            }
            println!();
        }
    }
}

fn main() {
    let max_n = env::args()
        .nth(1)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(16);

    println!("=== Non-nesting rook grid package ===");
    println!("Checking all partitions of size <= {max_n}.\n");

    let mut prev_delta_rr = TestResult::new("Prev delta = (f_c-f_{c+1})/t is real-rooted");
    let mut prev_delta_nonneg = TestResult::new("Prev delta = (f_c-f_{c+1})/t has nonnegative coefficients");
    let mut delta_rr = TestResult::new("Delta_c(mu) is real-rooted");
    let mut delta_nonneg = TestResult::new("Delta_c(mu) has nonnegative coefficients");
    let mut delta_base = TestResult::new("Delta_c(mu) << f_{c+1}");
    let mut delta_left = TestResult::new("Delta_c(mu) << G(c+1,s)");
    let mut delta_right = TestResult::new("Delta_c(mu) << G(c,s)");
    let mut adjacent = TestResult::new("G(c+1,s) << G(c,s)");
    let mut diagonal = TestResult::new("G(c+1,s) << G(c,s+1)");
    let mut anti_diagonal = TestResult::new("G(c+1,s+1) << G(c,s)");

    for n in 2..=max_n {
        for mu in partitions(n) {
            if mu.len() < 2 {
                continue;
            }
            let ell = mu.len();
            let m = mu[ell - 1];
            let mp = &mu[..ell - 1];
            let big_m = mp[0];

            let f: Vec<Vec<i64>> = (0..=big_m)
                .map(|c| non_nesting_rook_polynomial(&strip_columns(mp, c)))
                .collect();

            let prev_delta: Vec<Vec<i64>> = (0..big_m)
                .map(|c| poly_div_t(&poly_sub(&f[c], &f[c + 1])))
                .collect();

            let delta_mu: Vec<Vec<i64>> = (0..m)
                .map(|c| poly_add(&prev_delta[c], &f[c + 1]))
                .collect();

            for c in 0..big_m {
                prev_delta_rr.check(is_real_rooted(&prev_delta[c]), || {
                    format!(
                        "mu={:?} c={}: prev_delta = {}",
                        mu,
                        c,
                        format_poly(&prev_delta[c])
                    )
                });
                prev_delta_nonneg.check(prev_delta[c].iter().all(|&x| x >= 0), || {
                    format!(
                        "mu={:?} c={}: prev_delta = {}",
                        mu,
                        c,
                        format_poly(&prev_delta[c])
                    )
                });
            }

            for c in 0..m {
                delta_rr.check(is_real_rooted(&delta_mu[c]), || {
                    format!(
                        "mu={:?} c={}: Delta = {}",
                        mu,
                        c,
                        format_poly(&delta_mu[c])
                    )
                });
                delta_nonneg.check(delta_mu[c].iter().all(|&x| x >= 0), || {
                    format!(
                        "mu={:?} c={}: Delta = {}",
                        mu,
                        c,
                        format_poly(&delta_mu[c])
                    )
                });
                delta_base.check(interlaces(&delta_mu[c], &f[c + 1]) == Some(true), || {
                    format!(
                        "mu={:?} c={}: Delta << f_{{c+1}} failed: {} << {}",
                        mu,
                        c,
                        format_poly(&delta_mu[c]),
                        format_poly(&f[c + 1])
                    )
                });

                let mut g_c1 = f[c + 1].clone();
                let mut g_c = poly_add(&f[c], &poly_mul_t(&f[c + 1]));

                delta_left.check(interlaces(&delta_mu[c], &g_c1) == Some(true), || {
                    format!(
                        "mu={:?} c={} s={}: Delta << G(c+1,s) failed: {} << {}",
                        mu,
                        c,
                        c + 1,
                        format_poly(&delta_mu[c]),
                        format_poly(&g_c1)
                    )
                });
                delta_right.check(interlaces(&delta_mu[c], &g_c) == Some(true), || {
                    format!(
                        "mu={:?} c={} s={}: Delta << G(c,s) failed: {} << {}",
                        mu,
                        c,
                        c + 1,
                        format_poly(&delta_mu[c]),
                        format_poly(&g_c)
                    )
                });
                adjacent.check(interlaces(&g_c1, &g_c) == Some(true), || {
                    format!(
                        "mu={:?} c={} s={}: G(c+1,s) << G(c,s) failed: {} << {}",
                        mu,
                        c,
                        c + 1,
                        format_poly(&g_c1),
                        format_poly(&g_c)
                    )
                });

                for s in (c + 1)..m {
                    let t_next = poly_mul_t(&f[s + 1]);
                    let g_c1_next = poly_add(&g_c1, &t_next);
                    let g_c_next = poly_add(&g_c, &t_next);

                    diagonal.check(interlaces(&g_c1, &g_c_next) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: diagonal failed: {} << {}",
                            mu,
                            c,
                            s,
                            format_poly(&g_c1),
                            format_poly(&g_c_next)
                        )
                    });
                    anti_diagonal.check(interlaces(&g_c1_next, &g_c) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: anti-diagonal failed: {} << {}",
                            mu,
                            c,
                            s,
                            format_poly(&g_c1_next),
                            format_poly(&g_c)
                        )
                    });

                    g_c1 = g_c1_next;
                    g_c = g_c_next;

                    delta_left.check(interlaces(&delta_mu[c], &g_c1) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: Delta << G(c+1,s) failed: {} << {}",
                            mu,
                            c,
                            s + 1,
                            format_poly(&delta_mu[c]),
                            format_poly(&g_c1)
                        )
                    });
                    delta_right.check(interlaces(&delta_mu[c], &g_c) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: Delta << G(c,s) failed: {} << {}",
                            mu,
                            c,
                            s + 1,
                            format_poly(&delta_mu[c]),
                            format_poly(&g_c)
                        )
                    });
                    adjacent.check(interlaces(&g_c1, &g_c) == Some(true), || {
                        format!(
                            "mu={:?} c={} s={}: G(c+1,s) << G(c,s) failed: {} << {}",
                            mu,
                            c,
                            s + 1,
                            format_poly(&g_c1),
                            format_poly(&g_c)
                        )
                    });
                }
            }
        }
    }

    prev_delta_rr.print();
    prev_delta_nonneg.print();
    delta_rr.print();
    delta_nonneg.print();
    delta_base.print();
    delta_left.print();
    delta_right.print();
    adjacent.print();
    diagonal.print();
    anti_diagonal.print();
}
