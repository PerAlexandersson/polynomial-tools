use combpoly::permutation::all_permutations;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};
use polynomial_tools::real_rootedness::{
    check_weak_interlacing, ehrhart_to_hstar, format_poly, gamma_coefficients, is_palindromic,
    is_real_rooted,
};
/// h*-vectors of key polytopes for staircase partitions ρ_n = (n-1,...,1,0).
///
/// Bespoke algorithm: compute κ_{σ(kρ_n)}(1,...,1) using Demazure operators
/// π_i, tracking exponent vectors only at "live" variables.
use std::collections::{BTreeMap, HashMap};
use std::time::Instant;

/// Bubble-sort reduced word: 0-indexed positions where swaps occur.
fn sorting_word(alpha: &[u32]) -> Vec<usize> {
    let mut perm = alpha.to_vec();
    let n = perm.len();
    let mut word = Vec::new();
    loop {
        let mut swapped = false;
        for i in 0..n - 1 {
            if perm[i] < perm[i + 1] {
                perm.swap(i, i + 1);
                word.push(i);
                swapped = true;
            }
        }
        if !swapped {
            break;
        }
    }
    word
}

/// Compute κ_α(1,...,1) via Demazure operators with live-variable compression.
fn key_eval_ones(n: usize, alpha: &[u32]) -> u64 {
    let word = sorting_word(alpha);
    if word.is_empty() {
        return 1;
    }
    let mut lambda = alpha.to_vec();
    lambda.sort_unstable_by(|a, b| b.cmp(a));
    let l = word.len();

    // live_at[t] = variables whose exponents matter at step t
    let mut live_at: Vec<Vec<usize>> = vec![vec![]; l];
    let mut needed = vec![false; n];
    for t in (0..l).rev() {
        needed.fill(false);
        for s in t..l {
            needed[word[s]] = true;
            needed[word[s] + 1] = true;
        }
        live_at[t] = (0..n).filter(|&v| needed[v]).collect();
    }

    // Signed counts: intermediate steps can have negative coefficients
    let mut state: HashMap<Vec<u32>, i64> = HashMap::new();
    let init_key: Vec<u32> = live_at[0].iter().map(|&v| lambda[v]).collect();
    state.insert(init_key, 1);

    for t in 0..l {
        let op = word[t];
        let live_now = &live_at[t];
        let pos_i = live_now.iter().position(|&v| v == op).unwrap();
        let pos_ip1 = live_now.iter().position(|&v| v == op + 1).unwrap();

        let live_next: &[usize] = if t + 1 < l { &live_at[t + 1] } else { &[] };
        let proj: Vec<(usize, usize)> = live_next
            .iter()
            .enumerate()
            .map(|(np, &v)| (np, live_now.iter().position(|&u| u == v).unwrap()))
            .collect();
        let next_len = live_next.len();
        let mut new_state: HashMap<Vec<u32>, i64> = HashMap::new();

        for (exps, &count) in &state {
            if count == 0 {
                continue;
            }
            let ai = exps[pos_i] as i64;
            let aip1 = exps[pos_ip1] as i64;
            let p = ai + 1;
            let q = aip1;

            if p == q {
                continue;
            } // π_i(x^a) = 0

            let (sign, hi, lo, spread) = if p > q {
                (1i64, ai, aip1, (p - q - 1) as u32)
            } else {
                (-1i64, aip1 - 1, ai + 1, (q - p - 1) as u32)
            };

            for j in 0..=(spread as i64) {
                let new_ai = (hi - j) as u32;
                let new_aip1 = (lo + j) as u32;

                let mut new_key = vec![0u32; next_len];
                for &(np, op_idx) in &proj {
                    let v = live_next[np];
                    new_key[np] = if v == op {
                        new_ai
                    } else if v == op + 1 {
                        new_aip1
                    } else {
                        exps[op_idx]
                    };
                }
                *new_state.entry(new_key).or_insert(0) += sign * count;
            }
        }

        new_state.retain(|_, v| *v != 0);
        state = new_state;
    }

    let total: i64 = state.values().sum();
    assert!(
        total >= 0,
        "key_eval negative: {} for alpha {:?}",
        total,
        alpha
    );
    total as u64
}

fn staircase_ehrhart(n: usize, sigma: &[u8]) -> (Vec<BigRational>, Vec<i64>) {
    let ell = inv_count(sigma);
    let staircase: Vec<u32> = (0..n).map(|j| (n - 1 - j) as u32).collect();
    let id: Vec<u8> = (1..=n as u8).collect();
    let w0: Vec<u8> = (1..=n as u8).rev().collect();

    if *sigma == *id {
        return (vec![BigRational::one()], vec![1]);
    }
    if *sigma == *w0 {
        let d = n * (n - 1) / 2;
        let e = w0_ehrhart(d);
        let h = ehrhart_to_hstar(&e);
        return (e, h);
    }

    let mut points: Vec<(i64, BigRational)> = Vec::new();
    for k in 0..=(ell as u32) {
        let alpha: Vec<u32> = (0..n)
            .map(|j| staircase[(sigma[j] as usize) - 1] * k)
            .collect();
        let val = key_eval_ones(n, &alpha);
        points.push((k as i64, BigRational::from(BigInt::from(val))));
    }

    let e = lagrange_interpolation(&points);
    let h = ehrhart_to_hstar(&e);
    (e, h)
}

fn lagrange_interpolation(points: &[(i64, BigRational)]) -> Vec<BigRational> {
    let n = points.len();
    if n == 0 {
        return vec![];
    }
    let mut result = vec![BigRational::zero(); n];
    for i in 0..n {
        let (xi, ref yi) = points[i];
        let mut basis = vec![BigRational::zero(); n];
        basis[0] = BigRational::one();
        let mut denom = BigRational::one();
        let mut deg = 0;
        for j in 0..n {
            if j == i {
                continue;
            }
            let (xj, _) = points[j];
            denom *= BigRational::from(BigInt::from(xi - xj));
            let mut nb = vec![BigRational::zero(); n];
            for d in 0..=deg {
                nb[d + 1] += &basis[d];
                nb[d] -= &basis[d] * BigRational::from(BigInt::from(xj));
            }
            basis = nb;
            deg += 1;
        }
        for d in 0..n {
            result[d] += yi * &basis[d] / &denom;
        }
    }
    while result.len() > 1 && result.last().map_or(false, |c| c.is_zero()) {
        result.pop();
    }
    result
}

fn inv_count(perm: &[u8]) -> usize {
    let n = perm.len();
    let mut c = 0;
    for i in 0..n {
        for j in i + 1..n {
            if perm[i] > perm[j] {
                c += 1;
            }
        }
    }
    c
}
fn descent_set(perm: &[u8]) -> Vec<usize> {
    (0..perm.len() - 1)
        .filter(|&i| perm[i] > perm[i + 1])
        .map(|i| i + 1)
        .collect()
}
fn perm_str(perm: &[u8]) -> String {
    perm.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join("")
}
fn w0_conjugate(perm: &[u8]) -> Vec<u8> {
    let n = perm.len();
    (0..n).map(|i| (n as u8) + 1 - perm[n - 1 - i]).collect()
}
fn inverse_perm(perm: &[u8]) -> Vec<u8> {
    let n = perm.len();
    let mut inv = vec![0u8; n];
    for i in 0..n {
        inv[(perm[i] - 1) as usize] = (i + 1) as u8;
    }
    inv
}
fn w0_ehrhart(d: usize) -> Vec<BigRational> {
    let mut c = vec![BigRational::zero(); d + 1];
    for i in 0..=d {
        let mut b = BigRational::one();
        for j in 0..i {
            b *= BigRational::from(BigInt::from(d - j));
            b /= BigRational::from(BigInt::from(j + 1));
        }
        c[i] = b;
    }
    c
}
fn bruhat_covers(sigma: &[u8]) -> Vec<Vec<u8>> {
    let n = sigma.len();
    let mut covers = Vec::new();
    for a in 0..n {
        for b in a + 1..n {
            if sigma[a] > sigma[b] {
                let lo = sigma[b];
                let hi = sigma[a];
                if !(a + 1..b).any(|c| sigma[c] > lo && sigma[c] < hi) {
                    let mut tau = sigma.to_vec();
                    tau.swap(a, b);
                    covers.push(tau);
                }
            }
        }
    }
    covers
}

fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(5);

    for n in 3..=max_n {
        let d = n * (n - 1) / 2;
        let parts: Vec<u32> = (1..n as u32).rev().collect();
        println!("================================================================");
        println!("STAIRCASE ρ_{} = {:?},  n = {},  d = {}", n, parts, n, d);
        println!("================================================================");

        let perms = all_permutations(n as u8);
        let t0 = Instant::now();
        let mut hstar_map: BTreeMap<Vec<u8>, Vec<i64>> = BTreeMap::new();

        for (idx, perm) in perms.iter().enumerate() {
            let (_, hstar) = staircase_ehrhart(n, perm);
            hstar_map.insert(perm.clone(), hstar);
            if (idx + 1) % 20 == 0 || idx + 1 == perms.len() {
                eprint!(
                    "\r  computed {}/{} ({:.1?})   ",
                    idx + 1,
                    perms.len(),
                    t0.elapsed()
                );
            }
        }
        eprintln!();

        let mut by_length: BTreeMap<usize, Vec<&Vec<u8>>> = BTreeMap::new();
        for perm in &perms {
            by_length.entry(inv_count(perm)).or_default().push(perm);
        }

        println!("\nh*-vectors by Bruhat length:");
        for (len, ps) in &by_length {
            for perm in ps {
                let h = &hstar_map[*perm];
                let conj = w0_conjugate(perm);
                let inv = inverse_perm(perm);
                let rr = if h.len() <= 1 {
                    true
                } else {
                    is_real_rooted(h)
                };
                let pal = is_palindromic(h);
                let gam = if pal { gamma_coefficients(h) } else { None };
                let gp = gam.as_ref().map_or(false, |g| g.iter().all(|&c| c >= 0));
                let mut fl = Vec::new();
                if rr {
                    fl.push("RR");
                }
                if pal {
                    fl.push("pal");
                }
                if gp {
                    fl.push("γ+");
                }
                let cs = hstar_map.get(&conj).map_or(false, |ch| ch == h) && conj != **perm;
                let is = hstar_map.get(&inv).map_or(false, |ih| ih == h) && inv != **perm;
                println!(
                    "  ℓ={:2} σ={} des={:?}  h*={:<50} [{}]{}{}",
                    len,
                    perm_str(perm),
                    descent_set(perm),
                    format_poly(h),
                    fl.join(","),
                    if cs { "  =w₀conj" } else { "" },
                    if is { "  =inv" } else { "" }
                );
            }
        }

        println!("\nIdentical h*-vectors:");
        let mut h2p: BTreeMap<Vec<i64>, Vec<&Vec<u8>>> = BTreeMap::new();
        for perm in &perms {
            h2p.entry(hstar_map[perm].clone()).or_default().push(perm);
        }
        for (h, g) in &h2p {
            if g.len() > 1 {
                let ns: Vec<String> = g.iter().map(|p| perm_str(p)).collect();
                println!("  {} ← {} ({})", format_poly(h), ns.join(", "), g.len());
            }
        }

        println!("\nBruhat cover interlacing:");
        let (mut total, mut pass, mut fail) = (0, 0, 0);
        for sigma in &perms {
            for tau in &bruhat_covers(sigma) {
                let hs = &hstar_map[sigma];
                let ht = &hstar_map[tau];
                total += 1;
                if hs.len() <= 1 && ht.len() <= 1 {
                    pass += 1;
                    continue;
                }
                if check_weak_interlacing(ht, hs) == Some(true) {
                    pass += 1;
                } else {
                    fail += 1;
                    if fail <= 10 {
                        println!(
                            "  FAIL: {} ⋗ {}  h*(σ)={} vs h*(τ)={}",
                            perm_str(sigma),
                            perm_str(tau),
                            format_poly(hs),
                            format_poly(ht)
                        );
                    }
                }
            }
        }
        if fail == 0 {
            println!("  {}/{} — ALL PASS", pass, total);
        } else {
            println!("  {}/{} — {} failures", pass, total, fail);
        }
        println!();
    }
}
