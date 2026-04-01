/// Check: does the interlacing chain of L^{(p)} DIRECTLY imply stability
/// of Phi(x,t) = sum x^p L^{(p)}(t)?
///
/// The Borcea-Brändén characterization: for consecutive exponents,
/// f_0 + f_1 x + ... + f_d x^d is stable iff {f_i} is an interlacing sequence.
///
/// For non-consecutive exponents, this doesn't directly apply.
/// But we can check: if we RENUMBER the positions to be consecutive
/// (i.e., define Psi(y,t) = sum_{i=1}^m y^i L^{(p_i)}(t)),
/// is Psi stable?
///
/// Key test: Psi stable iff the interlacing chain holds.
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{check_weak_interlacing, is_real_rooted};
use std::collections::BTreeMap;

fn build_poly(vals: &[usize]) -> Vec<i64> {
    if vals.is_empty() {
        return vec![0];
    }
    let max_s = *vals.iter().max().unwrap();
    let mut coeffs = vec![0i64; max_s + 1];
    for &s in vals {
        coeffs[s] += 1;
    }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
        coeffs.pop();
    }
    coeffs
}

fn valid_positions(s_mask: u64, n: u8) -> Vec<u8> {
    let mut positions = Vec::new();
    for p in 1..n {
        if (s_mask >> (p - 1)) & 1 == 1 && (p < 2 || (s_mask >> (p - 2)) & 1 == 0) {
            positions.push(p);
        }
    }
    if n >= 2 && (s_mask >> (n - 2)) & 1 == 0 {
        positions.push(n);
    }
    positions
}

fn cmul(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}
fn cadd(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 + b.0, a.1 + b.1)
}
fn cpow(base: (f64, f64), exp: u32) -> (f64, f64) {
    let mut r = (1.0, 0.0);
    for _ in 0..exp {
        r = cmul(r, base);
    }
    r
}
fn cabs(a: (f64, f64)) -> f64 {
    (a.0 * a.0 + a.1 * a.1).sqrt()
}

/// Evaluate Psi(y, t) = sum_{i=0}^{m-1} y^{i+1} L^{(p_i)}(t)
/// with CONSECUTIVE y-exponents (renumbered positions).
fn eval_psi(polys: &[Vec<i64>], y: (f64, f64), t: (f64, f64)) -> (f64, f64) {
    let mut result = (0.0, 0.0);
    for (i, poly) in polys.iter().enumerate() {
        // Evaluate L^{(p_i)}(t)
        let mut lt = (0.0, 0.0);
        for (k, &c) in poly.iter().enumerate() {
            lt = cadd(lt, cmul((c as f64, 0.0), cpow(t, k as u32)));
        }
        // Multiply by y^{i+1}
        let yi = cpow(y, (i + 1) as u32);
        result = cadd(result, cmul(yi, lt));
    }
    result
}

/// Check stability of Psi(y,t) = sum y^{i+1} f_i(t)
fn check_psi_stable(polys: &[Vec<i64>], num_tests: usize) -> bool {
    let mut rng: u64 = 12345;
    let next = |s: &mut u64| -> f64 {
        *s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*s >> 33) as f64) / (1u64 << 31) as f64 * 4.0 - 2.0
    };
    let next_pos = |s: &mut u64| -> f64 {
        *s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*s >> 33) as f64) / (1u64 << 31) as f64 * 2.0 + 0.01
    };

    for _ in 0..num_tests {
        let y = (next(&mut rng), next_pos(&mut rng)); // Im(y) > 0
        let t = (next(&mut rng), next_pos(&mut rng)); // Im(t) > 0
        let val = eval_psi(polys, y, t);
        if cabs(val) < 1e-8 {
            return false;
        }
    }
    true
}

fn check_chain(polys: &[Vec<i64>]) -> bool {
    for i in 0..polys.len() - 1 {
        let f = &polys[i];
        let g = &polys[i + 1];
        if f.len() <= 1 || g.len() <= 1 {
            continue;
        }
        let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
        match check_weak_interlacing(small, large) {
            Some(true) => {}
            _ => {
                return false;
            }
        }
    }
    true
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    println!("=== Does interlacing chain ↔ stability of Psi(y,t)? ===\n");

    let mut chain_ok = 0u32;
    let mut psi_stable = 0u32;
    let mut both = 0u32;
    let mut total = 0u32;

    for n in 4..=max_n {
        let perms = all_permutations(n);
        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        for (&mask, class) in &by_descent {
            if mask & 1 != 0 {
                continue;
            }
            let vp = valid_positions(mask, n);
            if vp.len() < 2 {
                continue;
            }

            total += 1;

            // Build L^{(p)} polynomials ordered by p
            let mut polys: Vec<Vec<i64>> = Vec::new();
            for &p in &vp {
                let vals: Vec<usize> = class
                    .iter()
                    .filter(|pi| pi.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                    .map(|pi| compute(pi, Stat::Swaps))
                    .collect();
                polys.push(build_poly(&vals));
            }

            let is_chain = check_chain(&polys);
            let is_stable = check_psi_stable(&polys, 500);

            if is_chain {
                chain_ok += 1;
            }
            if is_stable {
                psi_stable += 1;
            }
            if is_chain && is_stable {
                both += 1;
            }

            if is_chain != is_stable {
                let s_str = {
                    let mut s = String::from("{");
                    let mut first = true;
                    for i in 1..n {
                        if (mask >> (i - 1)) & 1 == 1 {
                            if !first {
                                s.push(',');
                            }
                            s.push_str(&i.to_string());
                            first = false;
                        }
                    }
                    s.push('}');
                    s
                };
                println!(
                    "  MISMATCH n={} S={}: chain={} stable={}",
                    n, s_str, is_chain, is_stable
                );
            }
        }

        println!(
            "n={}: chain={}/{} stable={}/{} both={}/{}",
            n, chain_ok, total, psi_stable, total, both, total
        );
    }

    println!(
        "\nchain ↔ Psi stable: {}",
        if chain_ok == psi_stable && chain_ok == both {
            "EQUIVALENT"
        } else {
            "NOT equivalent"
        }
    );
    println!(
        "Total: chain={}/{} Psi_stable={}/{}",
        chain_ok, total, psi_stable, total
    );
}
