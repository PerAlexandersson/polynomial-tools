//! Ordered set partitions with the big-block statistic, refined by the
//! position of the block containing the maximal element.
//!
//! For fixed j >= 2, let H_{n,j,p}(t) count ordered set partitions of [n]
//! where the block containing n is in position p, weighted by t^{bb_j}.
//! The corresponding bivariate polynomial is
//!   H_{n,j}(x,t) = sum_{p >= 1} H_{n,j,p}(t) x^{p-1}.
//!
//! This scanner checks whether the coefficient polynomials H_{n,j,p}(t)
//! appear to form an interlacing sequence in p, and whether H_{n,j}(x,t)
//! remains real-rooted after positive specializations x = c.

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

fn ordered_exact(max_n: usize, j: usize) -> Vec<Vec<Vec<i128>>> {
    let refined = compute_refined(max_n, j);
    let facts = factorials(max_n);
    let mut out = vec![Vec::<Vec<i128>>::new(); max_n + 1];

    for n in 0..=max_n {
        let mut rows = vec![Vec::<i128>::new(); refined[n].len()];
        for m in 0..refined[n].len() {
            rows[m] = vec![0i128; refined[n][m].len()];
            for b in 0..refined[n][m].len() {
                rows[m][b] = facts[m] * refined[n][m][b];
            }
        }
        out[n] = rows;
    }

    out
}

fn ordered_total(exact: &[Vec<Vec<i128>>]) -> Vec<Vec<i128>> {
    let max_n = exact.len() - 1;
    let mut out = vec![vec![1]];
    for rows in exact.iter().take(max_n + 1).skip(1) {
        let mut poly = vec![0i128; rows.len() + 1];
        for row in rows {
            for (b, &val) in row.iter().enumerate() {
                poly[b] += val;
            }
        }
        while poly.last().is_some_and(|&x| x == 0) {
            poly.pop();
        }
        out.push(poly);
    }
    out
}

fn position_refined(max_n: usize, j: usize) -> Vec<Vec<Vec<i128>>> {
    let exact = ordered_exact(max_n, j);
    let total = ordered_total(&exact);
    let mut out = vec![Vec::<Vec<i128>>::new(); max_n + 1];
    out[0] = vec![vec![1]];

    for n in 1..=max_n {
        let mut rows = vec![vec![0i128; n + 1]; n + 1];
        let ordinary = n - 1;

        for k in 0..=ordinary {
            let choose_block = binom(ordinary, k);
            let big_inc = usize::from(k + 1 >= j);
            let rem = ordinary - k;

            for left_size in 0..=rem {
                let choose_left = binom(rem, left_size);
                let right_size = rem - left_size;

                for p in 1..=left_size + 1 {
                    if p - 1 >= exact[left_size].len() {
                        continue;
                    }
                    let left_row = &exact[left_size][p - 1];
                    let right_poly = &total[right_size];

                    for (b_left, &left_val) in left_row.iter().enumerate() {
                        if left_val == 0 {
                            continue;
                        }
                        for (b_right, &right_val) in right_poly.iter().enumerate() {
                            if right_val == 0 {
                                continue;
                            }
                            rows[p][b_left + b_right + big_inc] +=
                                choose_block * choose_left * left_val * right_val;
                        }
                    }
                }
            }
        }

        while rows.last().is_some_and(|poly| poly.iter().all(|&x| x == 0)) {
            rows.pop();
        }
        for poly in &mut rows {
            while poly.last().is_some_and(|&x| x == 0) {
                poly.pop();
            }
        }
        out[n] = rows;
    }

    out
}

fn to_i64_poly(poly: &[i128]) -> Option<Vec<i64>> {
    poly.iter()
        .copied()
        .map(|c| i64::try_from(c).ok())
        .collect()
}

fn degree_i128(poly: &[i128]) -> Option<usize> {
    poly.iter().rposition(|&c| c != 0)
}

fn interlaces(a: &[i128], b: &[i128]) -> Option<bool> {
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

fn eval_x(position_rows: &[Vec<i128>], x: i128) -> Vec<i128> {
    let mut out = vec![0i128; position_rows.len() + 2];
    let mut x_pow = 1i128;
    for poly in position_rows.iter().skip(1) {
        for (b, &coeff) in poly.iter().enumerate() {
            out[b] += x_pow * coeff;
        }
        x_pow *= x;
    }
    while out.last().is_some_and(|&c| c == 0) {
        out.pop();
    }
    out
}

fn eval_block_line(exact_rows: &[Vec<i128>], a: i128, b: i128) -> Vec<i128> {
    let mut out = vec![0i128; exact_rows.len() * 2 + 2];
    let mut a_pow = 1i128;
    for (m, row) in exact_rows.iter().enumerate().skip(1) {
        a_pow *= a;
        let mut b_pow = 1i128;
        for (big, &coeff) in row.iter().enumerate() {
            out[m + big] += coeff * a_pow * b_pow;
            b_pow *= b;
        }
    }
    while out.last().is_some_and(|&c| c == 0) {
        out.pop();
    }
    out
}

fn eval_block_u(exact_rows: &[Vec<i128>], u: i128) -> Vec<i128> {
    let mut out = vec![0i128; exact_rows.len() + 2];
    let mut u_pow = 1i128;
    for row in exact_rows.iter().skip(1) {
        u_pow *= u;
        for (big, &coeff) in row.iter().enumerate() {
            out[big] += coeff * u_pow;
        }
    }
    while out.last().is_some_and(|&c| c == 0) {
        out.pop();
    }
    out
}

fn eval_position_line(position_rows: &[Vec<i128>], a: i128, b: i128) -> Vec<i128> {
    let mut out = vec![0i128; position_rows.len() * 2 + 2];
    let mut a_pow = 1i128;
    for (idx, row) in position_rows.iter().enumerate().skip(1) {
        if idx > 1 {
            a_pow *= a;
        }
        let mut b_pow = 1i128;
        for (big, &coeff) in row.iter().enumerate() {
            out[(idx - 1) + big] += coeff * a_pow * b_pow;
            b_pow *= b;
        }
    }
    while out.last().is_some_and(|&c| c == 0) {
        out.pop();
    }
    out
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

fn main() {
    let max_n = std::env::args()
        .nth(1)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(14);
    let max_j = std::env::args()
        .nth(2)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(6);

    println!("=== Ordered big blocks refined by position of n ===\n");

    for j in 2..=max_j {
        let exact = ordered_exact(max_n, j);
        let pos = position_refined(max_n, j);
        let mut block_interlacing_pass = 0usize;
        let mut block_interlacing_total = 0usize;
        let mut block_ineligible = 0usize;
        let mut block_fail = Vec::new();

        let mut block_reverse_interlacing_pass = 0usize;
        let mut block_reverse_interlacing_total = 0usize;
        let mut block_reverse_ineligible = 0usize;
        let mut block_reverse_fail = Vec::new();

        let mut coeff_interlacing_pass = 0usize;
        let mut coeff_interlacing_total = 0usize;
        let mut coeff_ineligible = 0usize;
        let mut coeff_fail = Vec::new();

        let mut coeff_reverse_interlacing_pass = 0usize;
        let mut coeff_reverse_interlacing_total = 0usize;
        let mut coeff_reverse_ineligible = 0usize;
        let mut coeff_reverse_fail = Vec::new();

        let mut eval_rr_pass = 0usize;
        let mut eval_rr_total = 0usize;
        let mut eval_rr_fail = Vec::new();

        let mut block_eval_rr_pass = 0usize;
        let mut block_eval_rr_total = 0usize;
        let mut block_eval_rr_fail = Vec::new();

        let mut block_line_rr_pass = 0usize;
        let mut block_line_rr_total = 0usize;
        let mut block_line_rr_fail = Vec::new();

        let mut position_line_rr_pass = 0usize;
        let mut position_line_rr_total = 0usize;
        let mut position_line_rr_fail = Vec::new();

        for (n, rows) in pos.iter().enumerate().take(max_n + 1).skip(1) {
            for m in 1..exact[n].len().saturating_sub(1) {
                match interlaces(&exact[n][m], &exact[n][m + 1]) {
                    Some(true) => {
                        block_interlacing_total += 1;
                        block_interlacing_pass += 1;
                    }
                    Some(false) => {
                        block_interlacing_total += 1;
                        if block_fail.len() < 10 {
                            block_fail.push((n, m, exact[n][m].clone(), exact[n][m + 1].clone()));
                        }
                    }
                    None => block_ineligible += 1,
                }
            }
            for m in (2..exact[n].len()).rev() {
                match interlaces(&exact[n][m], &exact[n][m - 1]) {
                    Some(true) => {
                        block_reverse_interlacing_total += 1;
                        block_reverse_interlacing_pass += 1;
                    }
                    Some(false) => {
                        block_reverse_interlacing_total += 1;
                        if block_reverse_fail.len() < 10 {
                            block_reverse_fail
                                .push((n, m, exact[n][m].clone(), exact[n][m - 1].clone()));
                        }
                    }
                    None => block_reverse_ineligible += 1,
                }
            }

            for &u in &[1i128, 2, 3] {
                let poly = eval_block_u(&exact[n], u);
                if let Some(p64) = to_i64_poly(&poly) {
                    block_eval_rr_total += 1;
                    if p64.len() <= 2 || is_real_rooted(&p64) {
                        block_eval_rr_pass += 1;
                    } else if block_eval_rr_fail.len() < 10 {
                        block_eval_rr_fail.push((n, u, poly.clone()));
                    }
                }
            }

            for p in 1..rows.len().saturating_sub(1) {
                match interlaces(&rows[p], &rows[p + 1]) {
                    Some(true) => {
                        coeff_interlacing_total += 1;
                        coeff_interlacing_pass += 1;
                    }
                    Some(false) => {
                        coeff_interlacing_total += 1;
                        if coeff_fail.len() < 10 {
                            coeff_fail.push((n, p, rows[p].clone(), rows[p + 1].clone()));
                        }
                    }
                    None => coeff_ineligible += 1,
                }
            }
            for p in (2..rows.len()).rev() {
                match interlaces(&rows[p], &rows[p - 1]) {
                    Some(true) => {
                        coeff_reverse_interlacing_total += 1;
                        coeff_reverse_interlacing_pass += 1;
                    }
                    Some(false) => {
                        coeff_reverse_interlacing_total += 1;
                        if coeff_reverse_fail.len() < 10 {
                            coeff_reverse_fail.push((n, p, rows[p].clone(), rows[p - 1].clone()));
                        }
                    }
                    None => coeff_reverse_ineligible += 1,
                }
            }

            for &x in &[1i128, 2, 3] {
                let poly = eval_x(rows, x);
                if let Some(p64) = to_i64_poly(&poly) {
                    eval_rr_total += 1;
                    if p64.len() <= 2 || is_real_rooted(&p64) {
                        eval_rr_pass += 1;
                    } else if eval_rr_fail.len() < 10 {
                        eval_rr_fail.push((n, x, poly.clone()));
                    }
                }
            }

            for &(a, b) in &[(1i128, 1i128), (2, 1), (1, 2), (2, 2)] {
                let poly = eval_block_line(&exact[n], a, b);
                if let Some(p64) = to_i64_poly(&poly) {
                    block_line_rr_total += 1;
                    if p64.len() <= 2 || is_real_rooted(&p64) {
                        block_line_rr_pass += 1;
                    } else if block_line_rr_fail.len() < 10 {
                        block_line_rr_fail.push((n, a, b, poly.clone()));
                    }
                }
            }

            for &(a, b) in &[(1i128, 1i128), (2, 1), (1, 2), (2, 2)] {
                let poly = eval_position_line(rows, a, b);
                if let Some(p64) = to_i64_poly(&poly) {
                    position_line_rr_total += 1;
                    if p64.len() <= 2 || is_real_rooted(&p64) {
                        position_line_rr_pass += 1;
                    } else if position_line_rr_fail.len() < 10 {
                        position_line_rr_fail.push((n, a, b, poly.clone()));
                    }
                }
            }
        }

        println!("--- j = {} ---", j);
        println!(
            "Block-coefficient interlacing: {}/{} passes ({} ineligible)",
            block_interlacing_pass, block_interlacing_total, block_ineligible
        );
        println!(
            "Block-coefficient interlacing, reversed: {}/{} passes ({} ineligible)",
            block_reverse_interlacing_pass,
            block_reverse_interlacing_total,
            block_reverse_ineligible
        );
        println!(
            "Positive u specialization real-rootedness: {}/{} tested",
            block_eval_rr_pass, block_eval_rr_total
        );
        println!(
            "Same-phase line tests for G_n(u,t): {}/{} tested",
            block_line_rr_pass, block_line_rr_total
        );
        println!(
            "Position-coefficient interlacing: {}/{} passes ({} ineligible)",
            coeff_interlacing_pass, coeff_interlacing_total, coeff_ineligible
        );
        println!(
            "Position-coefficient interlacing, reversed: {}/{} passes ({} ineligible)",
            coeff_reverse_interlacing_pass,
            coeff_reverse_interlacing_total,
            coeff_reverse_ineligible
        );
        println!(
            "Positive x specialization real-rootedness: {}/{} tested",
            eval_rr_pass, eval_rr_total
        );
        println!(
            "Same-phase line tests for H_n(x,t): {}/{} tested",
            position_line_rr_pass, position_line_rr_total
        );

        let sample_n = max_n.min(j + 3);
        println!("Sample H_{{{}, {}}}(x,t) coefficients in x:", sample_n, j);
        for p in 1..pos[sample_n].len() {
            println!("  p={}: {}", p, format_poly_i128(&pos[sample_n][p]));
        }
        println!("  x=1: {}", format_poly_i128(&eval_x(&pos[sample_n], 1)));
        println!("  x=2: {}", format_poly_i128(&eval_x(&pos[sample_n], 2)));
        println!("Sample G_{{{}, {}}}(u,t) coefficients in u:", sample_n, j);
        for m in 1..exact[sample_n].len() {
            println!("  m={}: {}", m, format_poly_i128(&exact[sample_n][m]));
        }

        if !block_fail.is_empty() {
            println!("First block-coefficient interlacing failures:");
            for (n, m, a, b) in &block_fail {
                println!(
                    "  n={}, m={} -> {}, m+1={} -> {}",
                    n,
                    m,
                    format_poly_i128(a),
                    m + 1,
                    format_poly_i128(b)
                );
            }
        }

        if !block_reverse_fail.is_empty() {
            println!("First reversed block-coefficient interlacing failures:");
            for (n, m, a, b) in &block_reverse_fail {
                println!(
                    "  n={}, m={} -> {}, m-1={} -> {}",
                    n,
                    m,
                    format_poly_i128(a),
                    m - 1,
                    format_poly_i128(b)
                );
            }
        }

        if !block_eval_rr_fail.is_empty() {
            println!("First positive-u real-rootedness failures:");
            for (n, u, poly) in &block_eval_rr_fail {
                println!("  n={}, u={} -> {}", n, u, format_poly_i128(poly));
            }
        }

        if !block_line_rr_fail.is_empty() {
            println!("First same-phase line failures for G_n(u,t):");
            for (n, a, b, poly) in &block_line_rr_fail {
                println!(
                    "  n={}, (u,t)=({}s,{}s) -> {}",
                    n,
                    a,
                    b,
                    format_poly_i128(poly)
                );
            }
        }

        if !coeff_fail.is_empty() {
            println!("First coefficient-interlacing failures:");
            for (n, p, a, b) in &coeff_fail {
                println!(
                    "  n={}, p={} -> {}, p+1={} -> {}",
                    n,
                    p,
                    format_poly_i128(a),
                    p + 1,
                    format_poly_i128(b)
                );
            }
        }

        if !coeff_reverse_fail.is_empty() {
            println!("First reversed coefficient-interlacing failures:");
            for (n, p, a, b) in &coeff_reverse_fail {
                println!(
                    "  n={}, p={} -> {}, p-1={} -> {}",
                    n,
                    p,
                    format_poly_i128(a),
                    p - 1,
                    format_poly_i128(b)
                );
            }
        }

        if !eval_rr_fail.is_empty() {
            println!("First positive-x real-rootedness failures:");
            for (n, x, poly) in &eval_rr_fail {
                println!("  n={}, x={} -> {}", n, x, format_poly_i128(poly));
            }
        }

        if !position_line_rr_fail.is_empty() {
            println!("First same-phase line failures for H_n(x,t):");
            for (n, a, b, poly) in &position_line_rr_fail {
                println!(
                    "  n={}, (x,t)=({}s,{}s) -> {}",
                    n,
                    a,
                    b,
                    format_poly_i128(poly)
                );
            }
        }

        println!();
    }
}
