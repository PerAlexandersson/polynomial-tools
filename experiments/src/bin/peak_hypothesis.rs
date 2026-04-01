//! Test inductive hypotheses for peak polynomial real-rootedness.
//! Tests ALL Ferrers boards with n rows up to a given size.

use std::collections::BTreeSet;

fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> {
    let n = perm.len();
    let mut visited: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut queue: BTreeSet<Vec<u8>> = BTreeSet::new();
    queue.insert(perm.to_vec());
    while let Some(current) = queue.pop_last() {
        for i in 0..n {
            for j in i + 1..n {
                if current[i] > current[j] {
                    let mut child = current.clone();
                    child.swap(i, j);
                    if !visited.contains(&child) {
                        queue.insert(child);
                    }
                }
            }
        }
        visited.insert(current);
    }
    visited.into_iter().collect()
}

fn board_to_312_perm(board: &[u8]) -> Vec<u8> {
    let n = board.len();
    let mut perm = vec![0u8; n];
    let mut used = vec![false; n + 1];
    for i in 0..n {
        let max_col = (board[i] as usize).min(n);
        for c in (1..=max_col).rev() {
            if !used[c] {
                perm[i] = c as u8;
                used[c] = true;
                break;
            }
        }
    }
    perm
}

fn peaks(w: &[u8]) -> usize {
    if w.len() < 3 {
        return 0;
    }
    (1..w.len() - 1)
        .filter(|&i| w[i - 1] < w[i] && w[i] > w[i + 1])
        .count()
}

fn poly_trim(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
    v
}
fn poly_is_zero(p: &[i64]) -> bool {
    p.iter().all(|&c| c == 0)
}
fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] += v;
    }
    poly_trim(&r)
}
fn poly_mul_t(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; p.len() + 1];
    for (i, &v) in p.iter().enumerate() {
        r[i + 1] = v;
    }
    poly_trim(&r)
}
fn poly_deg(p: &[i64]) -> usize {
    let t = poly_trim(p);
    if poly_is_zero(&t) {
        return 0;
    }
    t.len() - 1
}
fn poly_fmt(p: &[i64]) -> String {
    let p = poly_trim(p);
    if poly_is_zero(&p) {
        return "0".into();
    }
    let mut terms = vec![];
    for (i, &c) in p.iter().enumerate() {
        if c == 0 {
            continue;
        }
        match (c, i) {
            (c, 0) => terms.push(format!("{}", c)),
            (1, 1) => terms.push("t".into()),
            (c, 1) => terms.push(format!("{}t", c)),
            (1, e) => terms.push(format!("t^{}", e)),
            (c, e) => terms.push(format!("{}t^{}", c, e)),
        }
    }
    terms.join(" + ")
}

fn find_real_roots(coeffs: &[i64]) -> Option<Vec<f64>> {
    let roots = polynomial_tools::real_roots(coeffs)?;
    let mut float_roots: Vec<f64> = roots
        .iter()
        .map(|r| {
            let n = r.numer().to_str_radix(10).parse::<f64>().unwrap();
            let d = r.denom().to_str_radix(10).parse::<f64>().unwrap();
            n / d
        })
        .collect();
    float_roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Some(float_roots)
}

/// Check f ≼ g: roots of f interlace roots of g.
/// Convention: ... ≤ a_2 ≤ b_2 ≤ a_1 ≤ b_1 ≤ 0
/// where a_i = roots of f (sorted ascending), b_i = roots of g.
/// Requires deg(g) = deg(f) or deg(g) = deg(f)+1.
fn interlaces(f: &[i64], g: &[i64]) -> bool {
    let f = poly_trim(f);
    let g = poly_trim(g);
    if poly_is_zero(&f) {
        return true;
    }
    if poly_is_zero(&g) {
        return false;
    }
    let fr = match find_real_roots(&f) {
        Some(r) => r,
        None => return false,
    };
    let gr = match find_real_roots(&g) {
        Some(r) => r,
        None => return false,
    };
    let df = fr.len();
    let dg = gr.len();
    if dg != df && dg != df + 1 {
        return false;
    }
    let eps = 1e-6;
    if df == 0 {
        return true;
    }
    if dg == df + 1 {
        // g has one more root. Check: g[0] ≤ f[0] ≤ g[1] ≤ f[1] ≤ ... ≤ f[df-1] ≤ g[dg-1]
        for i in 0..df {
            if gr[i] > fr[i] + eps || fr[i] > gr[i + 1] + eps {
                return false;
            }
        }
        true
    } else {
        // Same degree. Check: g[0] ≤ f[0] ≤ g[1] ≤ f[1] ≤ ... ≤ g[df-1] ≤ f[df-1]
        for i in 0..df {
            if gr[i] > fr[i] + eps {
                return false;
            }
            if i + 1 < df && fr[i] > gr[i + 1] + eps {
                return false;
            }
        }
        true
    }
}

fn is_real_rooted(coeffs: &[i64]) -> bool {
    polynomial_tools::is_real_rooted(coeffs)
}

struct PeakData {
    m: usize,
    a: Vec<Vec<i64>>, // 1-indexed
    w: Vec<Vec<i64>>, // 1-indexed
}

fn compute_peak_data(board: &[u8]) -> PeakData {
    let n = board.len();
    let m = (board[0] as usize).min(n);
    let perm = board_to_312_perm(board);
    let ideal = bruhat_lower_ideal(&perm);

    let mut d = vec![vec![0i64]; m + 1];
    let mut u = vec![vec![0i64]; m + 1];

    for pi in &ideal {
        if pi.len() < 2 {
            // n=1: single element, no descent/ascent
            let k = pi[0] as usize;
            if k <= m {
                while u[k].len() < 1 {
                    u[k].push(0);
                }
                u[k][0] += 1;
            }
            continue;
        }
        let k = pi[0] as usize;
        if k > m {
            continue;
        }
        let pk = peaks(pi);
        let is_ascent = pi[0] < pi[1];
        let poly = if is_ascent { &mut u[k] } else { &mut d[k] };
        while poly.len() <= pk {
            poly.push(0);
        }
        poly[pk] += 1;
    }

    let mut a = vec![vec![0i64]; m + 1];
    let mut w = vec![vec![0i64]; m + 1];
    for k in 1..=m {
        a[k] = poly_add(&d[k], &u[k]);
        w[k] = poly_add(&poly_mul_t(&d[k]), &u[k]);
    }

    PeakData { m, a, w }
}

fn generate_boards(n: usize) -> Vec<Vec<u8>> {
    let mut results = vec![];
    let mut current = vec![];
    gen_boards_rec(n, n, 0, &mut current, &mut results);
    results
}

fn gen_boards_rec(
    n: usize,
    max_col: usize,
    depth: usize,
    current: &mut Vec<u8>,
    results: &mut Vec<Vec<u8>>,
) {
    if depth == n {
        results.push(current.clone());
        return;
    }
    let min_val = (depth + 1).max(if depth > 0 {
        current[depth - 1] as usize
    } else {
        1
    });
    for v in min_val..=max_col {
        current.push(v as u8);
        gen_boards_rec(n, max_col, depth + 1, current, results);
        current.pop();
    }
}

fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);

    println!("Testing on all Ferrers boards with n <= {}", max_n);
    println!("=========================================================\n");

    let mut total = 0;
    let mut f_w_rev = 0; // (a) W_m ≼ ... ≼ W_1
    let mut f_a_fwd = 0; // (b) A_1 ≼ ... ≼ A_m
    let mut f_cross = 0; // (c) W_k ≼ A_j for k >= j (equivalent: A_j ≼ ?, rethink)
    let mut f_total_rr = 0; // total peak poly real-rooted

    for n in 1..=max_n {
        let boards = generate_boards(n);
        for board in &boards {
            total += 1;
            let pd = compute_peak_data(board);
            let m = pd.m;
            if m <= 1 {
                continue;
            }

            // Check total peak poly real-rooted
            let mut total_peak = vec![0i64];
            for k in 1..=m {
                total_peak = poly_add(&total_peak, &pd.a[k]);
            }
            if !is_real_rooted(&total_peak) {
                f_total_rr += 1;
                println!(
                    "FAIL: total peak poly not RR for board {:?}: {}",
                    board,
                    poly_fmt(&total_peak)
                );
            }

            // (a) W reversed interlacing: W_m ≼ W_{m-1} ≼ ... ≼ W_1
            for i in (2..=m).rev() {
                let k1 = i;
                let k2 = i - 1;
                if !poly_is_zero(&pd.w[k1]) && !poly_is_zero(&pd.w[k2]) {
                    if !interlaces(&pd.w[k1], &pd.w[k2]) {
                        f_w_rev += 1;
                        if f_w_rev <= 3 {
                            println!(
                                "FAIL W_{} ≼ W_{} for {:?}: {} vs {}",
                                k1,
                                k2,
                                board,
                                poly_fmt(&pd.w[k1]),
                                poly_fmt(&pd.w[k2])
                            );
                        }
                    }
                }
            }

            // (b) A forward interlacing: A_1 ≼ A_2 ≼ ... ≼ A_m
            for i in 1..m {
                if !poly_is_zero(&pd.a[i]) && !poly_is_zero(&pd.a[i + 1]) {
                    if !interlaces(&pd.a[i], &pd.a[i + 1]) {
                        f_a_fwd += 1;
                        if f_a_fwd <= 3 {
                            println!(
                                "FAIL A_{} ≼ A_{} for {:?}: {} vs {}",
                                i,
                                i + 1,
                                board,
                                poly_fmt(&pd.a[i]),
                                poly_fmt(&pd.a[i + 1])
                            );
                        }
                    }
                }
            }

            // (c) Cross: W_k ≼ A_j for k ≥ j (W has smaller roots)
            for k in 1..=m {
                for j in 1..=k {
                    if !poly_is_zero(&pd.w[k]) && !poly_is_zero(&pd.a[j]) {
                        if !interlaces(&pd.w[k], &pd.a[j]) {
                            f_cross += 1;
                            if f_cross <= 3 {
                                println!(
                                    "FAIL W_{} ≼ A_{} for {:?}: {} vs {}",
                                    k,
                                    j,
                                    board,
                                    poly_fmt(&pd.w[k]),
                                    poly_fmt(&pd.a[j])
                                );
                            }
                        }
                    }
                }
            }
        }
        println!(
            "n={}: {} boards, W_rev={}, A_fwd={}, cross={}, total_rr={}",
            n, total, f_w_rev, f_a_fwd, f_cross, f_total_rr
        );
    }

    println!("\n=========================================================");
    println!(
        "Total: {}, W_rev fails: {}, A_fwd fails: {}, cross fails: {}, total_rr fails: {}",
        total, f_w_rev, f_a_fwd, f_cross, f_total_rr
    );
}
