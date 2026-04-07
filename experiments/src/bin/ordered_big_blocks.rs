//! Ordered set partitions with the big-block statistic.
//!
//! Let O_{n,j}(t) = sum t^{bb_j(pi)}, where the sum ranges over ordered set
//! partitions of [n] and bb_j(pi) counts blocks of size at least j.
//!
//! We compute O_{n,j}(t) by first refining ordinary set partitions by
//! (number of blocks, number of big blocks), then weighting each partition
//! by m! according to its number of blocks.

use polynomial_tools::{check_weak_interlacing, is_real_rooted};

fn binom(n: usize, k: usize) -> i128 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut num: i128 = 1;
    let mut den: i128 = 1;
    for i in 0..k {
        num *= (n - i) as i128;
        den *= (i + 1) as i128;
    }
    num / den
}

fn factorials(max_n: usize) -> Vec<i128> {
    let mut fact = vec![1i128; max_n + 1];
    for n in 1..=max_n {
        fact[n] = fact[n - 1] * n as i128;
    }
    fact
}

fn compute_refined(max_n: usize, j: usize) -> Vec<Vec<Vec<i128>>> {
    // dp[n][m][b] = number of set partitions of [n] with m blocks and b big blocks.
    let mut dp = vec![Vec::<Vec<i128>>::new(); max_n + 1];
    dp[0] = vec![vec![1]];

    for n in 0..max_n {
        let mut next = vec![vec![0i128; n + 2]; n + 2];

        for k in 0..=n {
            let choose = binom(n, k);
            let prev_n = n - k;
            let big_inc = usize::from(k + 1 >= j);

            for m_prev in 0..dp[prev_n].len() {
                for b_prev in 0..dp[prev_n][m_prev].len() {
                    let val = dp[prev_n][m_prev][b_prev];
                    if val == 0 {
                        continue;
                    }
                    next[m_prev + 1][b_prev + big_inc] += choose * val;
                }
            }
        }

        while next.last().is_some_and(|row| row.iter().all(|&x| x == 0)) {
            next.pop();
        }
        for row in &mut next {
            while row.last().is_some_and(|&x| x == 0) {
                row.pop();
            }
        }
        dp[n + 1] = next;
    }

    dp
}

fn ordered_big_block_polys(max_n: usize, j: usize) -> Vec<Vec<i128>> {
    let refined = compute_refined(max_n, j);
    let facts = factorials(max_n);

    let mut polys = vec![vec![1]];
    for n in 1..=max_n {
        let mut poly = vec![0i128; n + 1];
        for m in 0..refined[n].len() {
            let weight = facts[m];
            for b in 0..refined[n][m].len() {
                poly[b] += weight * refined[n][m][b];
            }
        }
        while poly.last().is_some_and(|&x| x == 0) {
            poly.pop();
        }
        polys.push(poly);
    }
    polys
}

fn degree_i128(poly: &[i128]) -> Option<usize> {
    poly.iter().rposition(|&c| c != 0)
}

fn to_i64_poly(poly: &[i128]) -> Option<Vec<i64>> {
    poly.iter()
        .copied()
        .map(|c| i64::try_from(c).ok())
        .collect()
}

fn format_poly_i128(coeffs: &[i128]) -> String {
    let mut terms = Vec::new();
    for (i, &c) in coeffs.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let term = match (c, i) {
            (_, 0) => format!("{}", c),
            (1, 1) => "t".to_string(),
            (-1, 1) => "-t".to_string(),
            (_, 1) => format!("{}t", c),
            (1, e) => format!("t^{}", e),
            (-1, e) => format!("-t^{}", e),
            (_, e) => format!("{}t^{}", c, e),
        };
        terms.push(term);
    }
    if terms.is_empty() {
        return "0".to_string();
    }
    let mut result = terms[0].clone();
    for term in &terms[1..] {
        if let Some(rest) = term.strip_prefix('-') {
            result.push_str(" - ");
            result.push_str(rest);
        } else {
            result.push_str(" + ");
            result.push_str(term);
        }
    }
    result
}

fn interlaces_consecutive(a: &[i128], b: &[i128]) -> Option<bool> {
    let a64 = to_i64_poly(a)?;
    let b64 = to_i64_poly(b)?;
    let da = degree_i128(a)?;
    let db = degree_i128(b)?;
    if da <= db && db <= da + 1 {
        check_weak_interlacing(&a64, &b64)
    } else {
        None
    }
}

fn main() {
    let max_n = std::env::args()
        .nth(1)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(18);
    let max_j = std::env::args()
        .nth(2)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(6);

    println!("=== Ordered set partitions with big blocks ===\n");
    println!("O_(n,j)(t) = sum t^(# blocks of size >= j) over ordered set partitions.\n");

    for j in 2..=max_j {
        let polys = ordered_big_block_polys(max_n, j);
        let mut rr_pass = 0usize;
        let mut rr_testable = 0usize;
        let mut rr_fail = Vec::new();

        let mut interlace_pass = 0usize;
        let mut interlace_total = 0usize;
        let mut interlace_ineligible = 0usize;
        let mut interlace_fail = Vec::new();

        for n in 1..=max_n {
            if let Some(p64) = to_i64_poly(&polys[n]) {
                rr_testable += 1;
                if p64.len() <= 2 || is_real_rooted(&p64) {
                    rr_pass += 1;
                } else if rr_fail.len() < 10 {
                    rr_fail.push((n, polys[n].clone()));
                }
            }

            if n >= 2 {
                match interlaces_consecutive(&polys[n - 1], &polys[n]) {
                    Some(true) => {
                        interlace_total += 1;
                        interlace_pass += 1;
                    }
                    Some(false) => {
                        interlace_total += 1;
                        if interlace_fail.len() < 10 {
                            interlace_fail.push((n - 1, polys[n - 1].clone(), n, polys[n].clone()));
                        }
                    }
                    None => interlace_ineligible += 1,
                }
            }
        }

        println!("--- j = {} ---", j);
        println!("Real-rooted: {}/{} tested rows", rr_pass, rr_testable);
        println!(
            "Consecutive interlacing: {}/{} passes ({} ineligible)",
            interlace_pass, interlace_total, interlace_ineligible
        );

        println!("First rows:");
        for n in 0..=max_n.min(8) {
            println!("  n={}: {}", n, format_poly_i128(&polys[n]));
        }

        if !rr_fail.is_empty() {
            println!("First non-real-rooted examples:");
            for (n, p) in &rr_fail {
                println!("  n={}: {}", n, format_poly_i128(p));
            }
        }

        if !interlace_fail.is_empty() {
            println!("First interlacing failures:");
            for (n1, p1, n2, p2) in &interlace_fail {
                println!(
                    "  n={} -> {}, n={} -> {}",
                    n1,
                    format_poly_i128(p1),
                    n2,
                    format_poly_i128(p2)
                );
            }
        }
        println!();
    }
}
