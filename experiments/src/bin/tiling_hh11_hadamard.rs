//! Check if g_{n,r}(t) = sum_j C(j+1, n-2j-2r) t^j is real-rooted
//! for all valid (n, r) pairs. If so, then Q_{n,r} = f_r ⊙ g_{n,r}
//! is real-rooted by the Hadamard product theorem.
//!
//! Also check if Q_{n,r} polynomials have a common interlacing
//! so that P_n = sum_r Q_{n,r} is real-rooted.

use polynomial_tools::is_real_rooted;

fn binom(n: i64, k: i64) -> i64 {
    if k < 0 || k > n || n < 0 {
        return 0;
    }
    let k = k.min(n - k) as usize;
    let n = n as usize;
    let mut r = 1i64;
    for i in 0..k {
        r = r * (n - i) as i64 / (i + 1) as i64;
    }
    r
}

fn main() {
    println!("=== Checking real-rootedness of g_{{n,r}}(t) = sum_j C(j+1, n-2j-2r) t^j ===\n");

    let max_n = 60;
    let mut g_failures = 0;
    let mut g_total = 0;

    for n in 0..=max_n {
        for r in 0..=n / 2 {
            // Compute g_{n,r}(t)
            let mut coeffs = Vec::new();
            let mut min_j = None;
            for j in 0..=n as i64 {
                let s = n as i64 - 2 * j - 2 * r as i64;
                let c = binom(j + 1, s);
                if c > 0 && min_j.is_none() {
                    min_j = Some(j);
                }
                if c > 0 || min_j.is_some() {
                    coeffs.push(c);
                }
            }
            // Trim trailing zeros
            while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
                coeffs.pop();
            }

            if coeffs.len() <= 1 || coeffs.iter().all(|&c| c == 0) {
                continue;
            }
            g_total += 1;

            if !is_real_rooted(&coeffs) {
                println!("  g_{{{},{}}} NOT real-rooted: {:?}", n, r, coeffs);
                g_failures += 1;
            }
        }
    }

    if g_failures == 0 {
        println!(
            "  All g_{{n,r}} real-rooted for n <= {} ({} polynomials checked)",
            max_n, g_total
        );
    } else {
        println!("  {} failures out of {}", g_failures, g_total);
    }

    // Now check Q_{n,r}(t) = sum_j C(hr,j)*C(j+1,n-2j-2r) t^j for various h
    println!("\n=== Checking Q_{{n,r}} = f_r ⊙ g_{{n,r}} for h=2..6 ===\n");

    for h in 2..=6 {
        let mut q_failures = 0;
        let mut q_total = 0;
        for n in 0..=max_n {
            for r in 0..=n / 2 {
                let mut coeffs = Vec::new();
                let mut min_j = None;
                for j in 0..=n as i64 {
                    let s = n as i64 - 2 * j - 2 * r as i64;
                    let c = binom(h as i64 * r as i64, j) * binom(j + 1, s);
                    if c > 0 && min_j.is_none() {
                        min_j = Some(j);
                    }
                    if c > 0 || min_j.is_some() {
                        coeffs.push(c);
                    }
                }
                while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
                    coeffs.pop();
                }
                if coeffs.len() <= 2 || coeffs.iter().all(|&c| c == 0) {
                    continue;
                }
                q_total += 1;
                if !is_real_rooted(&coeffs) {
                    println!("  h={}: Q_{{{},{}}} NOT real-rooted: {:?}", h, n, r, coeffs);
                    q_failures += 1;
                }
            }
        }
        if q_failures == 0 {
            println!("  h={}: All Q_{{n,r}} real-rooted ({} checked)", h, q_total);
        } else {
            println!("  h={}: {} failures out of {}", h, q_failures, q_total);
        }
    }

    // Check: does P_n = sum_r Q_{n,r} remain real-rooted? (We already know this, but verify)
    println!("\n=== Verifying P_n = sum_r Q_{{n,r}} is real-rooted for h=2..6 ===\n");
    for h in 2..=6 {
        let mut ok = true;
        for n in 0..=max_n {
            let max_j = n / 2 + 5;
            let mut p = vec![0i64; max_j + 1];
            for r in 0..=n / 2 {
                for j in 0..=max_j {
                    let s = n as i64 - 2 * j as i64 - 2 * r as i64;
                    let c = binom(h as i64 * r as i64, j as i64) * binom(j as i64 + 1, s);
                    p[j] += c;
                }
            }
            while p.len() > 1 && *p.last().unwrap() == 0 {
                p.pop();
            }
            if p.len() > 2 && !is_real_rooted(&p) {
                println!("  h={}: P_{} NOT real-rooted", h, n);
                ok = false;
                break;
            }
        }
        if ok {
            println!("  h={}: All P_n real-rooted (n <= {})", h, max_n);
        }
    }
}
