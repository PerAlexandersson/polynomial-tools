//! Explore peak-generating polynomials D_k^lambda(t) and U_k^lambda(t)
//! on Ferrers boards.
//!
//! Usage: cargo run --bin peak_explore -- 2,3,3
//!        cargo run --bin peak_explore -- 3,4,4
//!        cargo run --bin peak_explore -- 2,3,4,4
//!        cargo run --bin peak_explore -- 4,4,4,4

use std::collections::BTreeSet;
use std::env;

// ── Permutation utilities ──────────────────────────────────────────

#[allow(dead_code)]
fn all_permutations(n: u8) -> Vec<Vec<u8>> {
    if n == 0 {
        return vec![vec![]];
    }
    let mut perms = Vec::new();
    let mut current: Vec<u8> = (1..=n).collect();
    loop {
        perms.push(current.clone());
        if !next_permutation(&mut current) {
            break;
        }
    }
    perms
}

#[allow(dead_code)]
fn next_permutation(perm: &mut [u8]) -> bool {
    let n = perm.len();
    if n <= 1 {
        return false;
    }
    let mut i = n - 2;
    loop {
        if perm[i] < perm[i + 1] {
            break;
        }
        if i == 0 {
            return false;
        }
        i -= 1;
    }
    let mut j = n - 1;
    while perm[j] <= perm[i] {
        j -= 1;
    }
    perm.swap(i, j);
    perm[i + 1..].reverse();
    true
}

// ── Bruhat order ───────────────────────────────────────────────────

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

// ── Ferrers board -> 312-avoiding permutation ──────────────────────

/// Convert a Ferrers board lambda = (l1, ..., ln) to its maximal
/// 312-avoiding permutation.  Greedy: from row 1 to n, place rook
/// in the rightmost available column <= min(lambda_i, n).
fn ferrers_to_perm(lambda: &[usize]) -> Vec<u8> {
    let n = lambda.len();
    let mut used = vec![false; n + 1]; // columns 1..n
    let mut perm = vec![0u8; n];

    for i in 0..n {
        let cap = lambda[i].min(n);
        // Find rightmost unused column in 1..=cap
        let mut best = 0;
        for c in (1..=cap).rev() {
            if !used[c] {
                best = c;
                break;
            }
        }
        assert!(best > 0, "No available column for row {} (lambda={:?})", i + 1, lambda);
        perm[i] = best as u8;
        used[best] = true;
    }
    perm
}

// ── Peak counting ──────────────────────────────────────────────────

fn count_peaks(w: &[u8]) -> usize {
    if w.len() < 3 {
        return 0;
    }
    (1..w.len() - 1)
        .filter(|&i| w[i - 1] < w[i] && w[i] > w[i + 1])
        .count()
}

// ── Polynomial utilities ───────────────────────────────────────────

/// Add two polynomials (coefficient vectors).
fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut result = vec![0i64; len];
    for (i, &v) in a.iter().enumerate() {
        result[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        result[i] += v;
    }
    result
}

/// Multiply polynomial by t (shift coefficients up by 1).
fn poly_mul_t(a: &[i64]) -> Vec<i64> {
    if a.is_empty() || (a.len() == 1 && a[0] == 0) {
        return vec![0];
    }
    let mut result = vec![0i64; a.len() + 1];
    for (i, &v) in a.iter().enumerate() {
        result[i + 1] = v;
    }
    result
}

fn poly_is_zero(a: &[i64]) -> bool {
    a.iter().all(|&c| c == 0)
}

fn poly_trim(a: &[i64]) -> Vec<i64> {
    let end = a.iter().rposition(|&c| c != 0);
    match end {
        Some(e) => a[..=e].to_vec(),
        None => vec![0],
    }
}

fn format_poly(coeffs: &[i64]) -> String {
    if coeffs.is_empty() {
        return "0".to_string();
    }
    let mut terms = Vec::new();
    for (i, &c) in coeffs.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let term = match (c, i) {
            (c, 0) => format!("{}", c),
            (1, 1) => "t".to_string(),
            (-1, 1) => "-t".to_string(),
            (c, 1) => format!("{}t", c),
            (1, e) => format!("t^{}", e),
            (-1, e) => format!("-t^{}", e),
            (c, e) => format!("{}t^{}", c, e),
        };
        terms.push(term);
    }
    if terms.is_empty() {
        return "0".to_string();
    }
    let mut result = terms[0].clone();
    for term in &terms[1..] {
        if term.starts_with('-') {
            result.push_str(" - ");
            result.push_str(&term[1..]);
        } else {
            result.push_str(" + ");
            result.push_str(term);
        }
    }
    result
}

// ── Real root finding ──────────────────────────────────────────────

/// Find all real roots of a polynomial using Sturm chains.
/// Returns sorted roots (as f64 approximations).
fn find_real_roots(coeffs: &[i64]) -> Option<Vec<f64>> {
    let roots = polynomial_tools::real_roots(coeffs)?;
    let mut float_roots: Vec<f64> = roots.iter().map(|r| {
        let n = r.numer().to_str_radix(10).parse::<f64>().unwrap();
        let d = r.denom().to_str_radix(10).parse::<f64>().unwrap();
        n / d
    }).collect();
    float_roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Some(float_roots)
}


fn is_real_rooted(coeffs: &[i64]) -> bool {
    polynomial_tools::is_real_rooted(coeffs)
}

/// Check if f interlaces g (written f << g).
///
/// Convention: f << g means
///   - both f and g are real-rooted with all roots <= 0
///   - deg(g) = deg(f) or deg(g) = deg(f) + 1
///   - the roots interleave: if deg(g) = deg(f) + 1, then
///       s_1 <= r_1 <= s_2 <= r_2 <= ... <= r_k <= s_{k+1}
///     if deg(g) = deg(f), then
///       s_1 <= r_1 <= s_2 <= r_2 <= ... <= s_k <= r_k
///     (where r_i are roots of f and s_j are roots of g)
///
/// Special cases:
///   - 0 << g for any real-rooted g (the zero polynomial interlaces everything)
///   - c << g for any nonzero constant c and real-rooted g with deg(g) <= 1
///     (a constant has no roots to interleave, so we only need deg(g) - deg(f) <= 1)
fn check_interlacing(f: &[i64], g: &[i64]) -> InterlaceResult {
    let f = poly_trim(f);
    let g = poly_trim(g);

    if poly_is_zero(&f) {
        return InterlaceResult::Yes; // zero interlaces everything
    }
    if poly_is_zero(&g) {
        return InterlaceResult::No("g is zero but f is not".into());
    }

    let f_roots = match find_real_roots(&f) {
        Some(r) => r,
        None => return InterlaceResult::No("f not real-rooted".into()),
    };
    let g_roots = match find_real_roots(&g) {
        Some(r) => r,
        None => return InterlaceResult::No("g not real-rooted".into()),
    };

    let df = f_roots.len();
    let dg = g_roots.len();

    // Check roots are non-positive
    let eps = 1e-6;
    for &r in &f_roots {
        if r > eps {
            return InterlaceResult::No(format!("f has positive root {:.6}", r));
        }
    }
    for &r in &g_roots {
        if r > eps {
            return InterlaceResult::No(format!("g has positive root {:.6}", r));
        }
    }

    // For f << g, we need dg = df or dg = df + 1.
    if dg != df && dg != df + 1 {
        return InterlaceResult::No(format!(
            "degree mismatch: deg(f)={}, deg(g)={}", df, dg
        ));
    }

    // If f is constant (df == 0), then dg must be 0 or 1, already checked above.
    // No roots to interleave, so it's automatic.
    if df == 0 {
        return InterlaceResult::Yes;
    }

    if dg == df + 1 {
        // s_1 <= r_1 <= s_2 <= r_2 <= ... <= r_{df} <= s_{df+1}
        for i in 0..df {
            if g_roots[i] > f_roots[i] + eps || f_roots[i] > g_roots[i + 1] + eps {
                return InterlaceResult::No(format!(
                    "interleaving fails at i={}: g[{}]={:.6} f[{}]={:.6} g[{}]={:.6}",
                    i, i, g_roots[i], i, f_roots[i], i + 1, g_roots[i + 1]
                ));
            }
        }
        InterlaceResult::Yes
    } else {
        // dg == df
        // s_1 <= r_1 <= s_2 <= r_2 <= ... <= s_k <= r_k
        let mut ok = true;
        for i in 0..df {
            if g_roots[i] > f_roots[i] + eps {
                ok = false;
                break;
            }
            if i + 1 < df && f_roots[i] > g_roots[i + 1] + eps {
                ok = false;
                break;
            }
        }
        if ok {
            return InterlaceResult::Yes;
        }

        // Also try the other interleaving: r_1 <= s_1 <= r_2 <= s_2 <= ...
        let mut ok2 = true;
        for i in 0..df {
            if f_roots[i] > g_roots[i] + eps {
                ok2 = false;
                break;
            }
            if i + 1 < df && g_roots[i] > f_roots[i + 1] + eps {
                ok2 = false;
                break;
            }
        }
        if ok2 {
            return InterlaceResult::Yes;
        }

        InterlaceResult::No("roots don't interleave in equal-degree case".into())
    }
}

enum InterlaceResult {
    Yes,
    No(String),
}

impl std::fmt::Display for InterlaceResult {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            InterlaceResult::Yes => write!(f, "YES"),
            InterlaceResult::No(reason) => write!(f, "NO ({})", reason),
        }
    }
}

/// Check if a sequence of polynomials is interlacing:
/// p_1, p_2, ..., p_m is interlacing if p_i interlaces p_j for all i < j.
fn check_sequence_interlacing(polys: &[(usize, Vec<i64>)]) -> bool {
    let mut all_ok = true;
    for i in 0..polys.len() {
        for j in i + 1..polys.len() {
            if poly_is_zero(&polys[i].1) || poly_is_zero(&polys[j].1) {
                continue;
            }
            match check_interlacing(&polys[i].1, &polys[j].1) {
                InterlaceResult::Yes => {}
                InterlaceResult::No(_) => {
                    all_ok = false;
                }
            }
        }
    }
    all_ok
}

// ── Compute D/U maps from a board ──────────────────────────────────

type PolyMap = std::collections::BTreeMap<u8, Vec<i64>>;

fn compute_du(lambda: &[usize]) -> (PolyMap, PolyMap) {
    let pi_lambda = ferrers_to_perm(lambda);
    let interval = bruhat_lower_ideal(&pi_lambda);

    let mut d_polys: PolyMap = std::collections::BTreeMap::new();
    let mut u_polys: PolyMap = std::collections::BTreeMap::new();

    for pi in &interval {
        let k = pi[0];
        let peaks = count_peaks(pi);

        if pi.len() < 2 {
            let entry = u_polys.entry(k).or_insert_with(Vec::new);
            while entry.len() <= peaks {
                entry.push(0);
            }
            entry[peaks] += 1;
            continue;
        }

        let is_descent = k > pi[1];
        let map = if is_descent { &mut d_polys } else { &mut u_polys };
        let entry = map.entry(k).or_insert_with(Vec::new);
        while entry.len() <= peaks {
            entry.push(0);
        }
        entry[peaks] += 1;
    }

    (d_polys, u_polys)
}

/// Verify the recursion:
///   D_k^{lambda+}(t) = sum_{j=1}^{k-1} [D_j^lambda(t) + U_j^lambda(t)]
///   U_k^{lambda+}(t) = sum_{j=k}^{lambda_1} [t*D_j^lambda(t) + U_j^lambda(t)]
/// where lambda+ = (m, lambda_1+1, ..., lambda_n+1).
fn verify_recursion(lambda: &[usize], lambda_plus: &[usize]) {
    let n = lambda.len();
    let m = lambda_plus[0];

    // Check that lambda+ has the right shape
    assert_eq!(lambda_plus.len(), n + 1, "lambda+ must have one more row");
    for i in 0..n {
        assert_eq!(lambda_plus[i + 1], lambda[i] + 1,
            "lambda+[{}] = {} should be lambda[{}]+1 = {}",
            i + 1, lambda_plus[i + 1], i, lambda[i] + 1);
    }

    println!("\nRecursion: lambda={:?} -> lambda+={:?}", lambda, lambda_plus);

    let (d_lam, u_lam) = compute_du(lambda);
    let (d_plus, u_plus) = compute_du(lambda_plus);

    let lambda_1 = lambda[0]; // max column in lambda
    let n_plus = lambda_plus.len(); // = n + 1
    let _ = n_plus;

    // For each k in 1..=m, compute the recursion prediction and compare
    let mut all_ok = true;

    for k in 1..=m as u8 {
        // D_k^{lambda+} = sum_{j=1}^{k-1} [D_j^lambda + U_j^lambda]
        let mut d_predicted = vec![0i64];
        for j in 1..k {
            let dj = d_lam.get(&j).map_or(vec![0], |v| v.clone());
            let uj = u_lam.get(&j).map_or(vec![0], |v| v.clone());
            d_predicted = poly_add(&d_predicted, &poly_add(&dj, &uj));
        }
        let d_predicted = poly_trim(&d_predicted);

        let d_actual = d_plus.get(&k).map_or(vec![0], |v| poly_trim(v));

        let d_match = d_predicted == d_actual;
        if !d_match {
            println!("  D_{}^+ MISMATCH: predicted={} actual={}",
                k, format_poly(&d_predicted), format_poly(&d_actual));
            all_ok = false;
        }

        // U_k^{lambda+} = sum_{j=k}^{lambda_1} [t*D_j^lambda + U_j^lambda]
        let mut u_predicted = vec![0i64];
        for j in k..=lambda_1 as u8 {
            let dj = d_lam.get(&j).map_or(vec![0], |v| v.clone());
            let uj = u_lam.get(&j).map_or(vec![0], |v| v.clone());
            u_predicted = poly_add(&u_predicted, &poly_add(&poly_mul_t(&dj), &uj));
        }
        let u_predicted = poly_trim(&u_predicted);

        let u_actual = u_plus.get(&k).map_or(vec![0], |v| poly_trim(v));

        let u_match = u_predicted == u_actual;
        if !u_match {
            println!("  U_{}^+ MISMATCH: predicted={} actual={}",
                k, format_poly(&u_predicted), format_poly(&u_actual));
            all_ok = false;
        }
    }

    if all_ok {
        println!("  All D_k and U_k match the recursion: OK");
    }
}

// ── Main logic ─────────────────────────────────────────────────────

fn main() {
    let args: Vec<String> = env::args().collect();

    let boards: Vec<Vec<usize>> = if args.len() > 1 {
        args[1..]
            .iter()
            .map(|s| {
                s.split(',')
                    .map(|x| x.trim().parse::<usize>().expect("invalid number"))
                    .collect()
            })
            .collect()
    } else {
        // Default test boards
        vec![
            vec![2, 3, 3],
            vec![3, 4, 4],
            vec![2, 3, 4, 4],
            vec![4, 4, 4, 4],
        ]
    };

    for lambda in &boards {
        process_board(lambda);
    }

    // Verify recursion: for pairs (lambda, lambda+)
    println!("\n{}", "=".repeat(72));
    println!("RECURSION VERIFICATION");
    println!("{}", "=".repeat(72));
    let recursion_pairs: Vec<(Vec<usize>, Vec<usize>)> = vec![
        // lambda -> lambda+ = (m, l1+1, ..., ln+1)
        // (2,3) -> (m, 3, 4) with m = 3 or 4 etc.
        (vec![2, 3], vec![2, 3, 4]),
        (vec![2, 3], vec![3, 3, 4]),
        (vec![2, 3, 3], vec![2, 3, 4, 4]),
        (vec![2, 3, 3], vec![3, 3, 4, 4]),
        (vec![3, 3, 3], vec![3, 4, 4, 4]),
        (vec![3, 3, 3], vec![4, 4, 4, 4]),
    ];
    for (lam, lam_plus) in &recursion_pairs {
        verify_recursion(lam, lam_plus);
    }
}

fn process_board(lambda: &[usize]) {
    let _n = lambda.len();
    println!("{}",  "=".repeat(72));
    println!("Ferrers board lambda = {:?}", lambda);

    // Convert to 312-avoiding permutation
    let pi_lambda = ferrers_to_perm(lambda);
    println!(
        "312-avoiding perm pi_lambda = [{}]",
        pi_lambda
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(", ")
    );

    // Compute Bruhat interval [id, pi_lambda]
    let interval = bruhat_lower_ideal(&pi_lambda);
    println!("Bruhat interval size: {}", interval.len());

    // Classify each permutation by (first entry, ascent/descent at pos 1, peaks)
    // D_k: pi_1 = k and k > pi_2  (position 1 is a descent)
    // U_k: pi_1 = k and k < pi_2  (position 1 is an ascent)
    // For length-1 permutations, there is no position 2, so we skip D/U classification.

    let mut d_polys: std::collections::BTreeMap<u8, Vec<i64>> = std::collections::BTreeMap::new();
    let mut u_polys: std::collections::BTreeMap<u8, Vec<i64>> = std::collections::BTreeMap::new();

    for pi in &interval {
        let k = pi[0];
        let peaks = count_peaks(pi);

        if pi.len() < 2 {
            // Single element: no ascent/descent. Count it as U_k (trivially an ascent).
            let entry = u_polys.entry(k).or_insert_with(Vec::new);
            while entry.len() <= peaks {
                entry.push(0);
            }
            entry[peaks] += 1;
            continue;
        }

        let is_descent = k > pi[1]; // position 1 is a descent

        if is_descent {
            let entry = d_polys.entry(k).or_insert_with(Vec::new);
            while entry.len() <= peaks {
                entry.push(0);
            }
            entry[peaks] += 1;
        } else {
            let entry = u_polys.entry(k).or_insert_with(Vec::new);
            while entry.len() <= peaks {
                entry.push(0);
            }
            entry[peaks] += 1;
        }
    }

    // Collect all keys
    let all_keys: BTreeSet<u8> = d_polys.keys().chain(u_polys.keys()).copied().collect();

    println!("\n--- D_k polynomials (pi_1=k, descent at pos 1) ---");
    for &k in &all_keys {
        let d = d_polys.get(&k).map_or(vec![0], |v| v.clone());
        println!("  D_{} = {}", k, format_poly(&d));
    }

    println!("\n--- U_k polynomials (pi_1=k, ascent at pos 1) ---");
    for &k in &all_keys {
        let u = u_polys.get(&k).map_or(vec![0], |v| v.clone());
        println!("  U_{} = {}", k, format_poly(&u));
    }

    // Build A_k = D_k + U_k
    println!("\n--- A_k = D_k + U_k ---");
    let mut a_polys: std::collections::BTreeMap<u8, Vec<i64>> = std::collections::BTreeMap::new();
    for &k in &all_keys {
        let d = d_polys.get(&k).map_or(vec![0], |v| v.clone());
        let u = u_polys.get(&k).map_or(vec![0], |v| v.clone());
        let a = poly_add(&d, &u);
        println!("  A_{} = {}", k, format_poly(&a));
        a_polys.insert(k, a);
    }

    // Build W_k = t*D_k + U_k
    println!("\n--- W_k = t*D_k + U_k ---");
    let mut w_polys: std::collections::BTreeMap<u8, Vec<i64>> = std::collections::BTreeMap::new();
    for &k in &all_keys {
        let d = d_polys.get(&k).map_or(vec![0], |v| v.clone());
        let u = u_polys.get(&k).map_or(vec![0], |v| v.clone());
        let w = poly_add(&poly_mul_t(&d), &u);
        println!("  W_{} = {}", k, format_poly(&w));
        w_polys.insert(k, w);
    }

    // Check real-rootedness
    println!("\n--- Real-rootedness ---");
    for &k in &all_keys {
        let d = d_polys.get(&k).map_or(vec![0], |v| v.clone());
        let u = u_polys.get(&k).map_or(vec![0], |v| v.clone());
        let a = a_polys.get(&k).unwrap();
        let w = w_polys.get(&k).unwrap();

        let d_rr = is_real_rooted(&d);
        let u_rr = is_real_rooted(&u);
        let a_rr = is_real_rooted(a);
        let w_rr = is_real_rooted(w);

        println!(
            "  k={}: D_k={} U_k={} A_k={} W_k={}",
            k,
            if d_rr { "RR" } else { "NOT-RR" },
            if u_rr { "RR" } else { "NOT-RR" },
            if a_rr { "RR" } else { "NOT-RR" },
            if w_rr { "RR" } else { "NOT-RR" },
        );
    }

    // Check interlacing: D_i interlaces U_j for i <= j
    println!("\n--- Interlacing: D_i interlaces U_j (i <= j) ---");
    let keys: Vec<u8> = all_keys.iter().copied().collect();
    for &i in &keys {
        for &j in &keys {
            if i > j {
                continue;
            }
            let d = d_polys.get(&i).map_or(vec![0], |v| v.clone());
            let u = u_polys.get(&j).map_or(vec![0], |v| v.clone());
            if poly_is_zero(&d) && poly_is_zero(&u) {
                continue;
            }
            let result = check_interlacing(&d, &u);
            println!("  D_{} interlaces U_{}: {}", i, j, result);
        }
    }

    // Check {D_k} is an interlacing sequence (D_1, D_2, ..., D_n: pairwise interlacing)
    println!("\n--- {{D_k}} interlacing sequence ---");
    let d_seq: Vec<(usize, Vec<i64>)> = keys
        .iter()
        .map(|&k| (k as usize, d_polys.get(&k).map_or(vec![0], |v| v.clone())))
        .filter(|(_, p)| !poly_is_zero(p))
        .collect();
    for i in 0..d_seq.len() {
        for j in i + 1..d_seq.len() {
            let result = check_interlacing(&d_seq[i].1, &d_seq[j].1);
            println!(
                "  D_{} interlaces D_{}: {}",
                d_seq[i].0, d_seq[j].0, result
            );
        }
    }
    let d_interlacing = check_sequence_interlacing(&d_seq);
    println!("  => {{{{D_k}}}} is interlacing: {}", if d_interlacing { "YES" } else { "NO" });

    // Check {U_k} is an interlacing sequence
    println!("\n--- {{U_k}} interlacing sequence ---");
    let u_seq: Vec<(usize, Vec<i64>)> = keys
        .iter()
        .map(|&k| (k as usize, u_polys.get(&k).map_or(vec![0], |v| v.clone())))
        .filter(|(_, p)| !poly_is_zero(p))
        .collect();
    for i in 0..u_seq.len() {
        for j in i + 1..u_seq.len() {
            let result = check_interlacing(&u_seq[i].1, &u_seq[j].1);
            println!(
                "  U_{} interlaces U_{}: {}",
                u_seq[i].0, u_seq[j].0, result
            );
        }
    }
    let u_interlacing = check_sequence_interlacing(&u_seq);
    println!("  => {{{{U_k}}}} is interlacing: {}", if u_interlacing { "YES" } else { "NO" });

    // Check {A_k} interlacing sequence
    println!("\n--- {{A_k}} interlacing sequence ---");
    let a_seq: Vec<(usize, Vec<i64>)> = keys
        .iter()
        .map(|&k| (k as usize, a_polys.get(&k).unwrap().clone()))
        .filter(|(_, p)| !poly_is_zero(p))
        .collect();
    for i in 0..a_seq.len() {
        for j in i + 1..a_seq.len() {
            let result = check_interlacing(&a_seq[i].1, &a_seq[j].1);
            println!(
                "  A_{} interlaces A_{}: {}",
                a_seq[i].0, a_seq[j].0, result
            );
        }
    }
    let a_interlacing = check_sequence_interlacing(&a_seq);
    println!("  => {{{{A_k}}}} is interlacing: {}", if a_interlacing { "YES" } else { "NO" });

    // Check {W_k} interlacing sequence
    println!("\n--- {{W_k}} interlacing sequence ---");
    let w_seq: Vec<(usize, Vec<i64>)> = keys
        .iter()
        .map(|&k| (k as usize, w_polys.get(&k).unwrap().clone()))
        .filter(|(_, p)| !poly_is_zero(p))
        .collect();
    for i in 0..w_seq.len() {
        for j in i + 1..w_seq.len() {
            let result = check_interlacing(&w_seq[i].1, &w_seq[j].1);
            println!(
                "  W_{} interlaces W_{}: {}",
                w_seq[i].0, w_seq[j].0, result
            );
        }
    }
    let w_interlacing = check_sequence_interlacing(&w_seq);
    println!("  => {{{{W_k}}}} is interlacing: {}", if w_interlacing { "YES" } else { "NO" });

    // Check ladder pair: A_k << W_k and W_k << A_{k+1}
    // AND the reverse: A_{k+1} << W_k and W_k << A_k
    println!("\n--- Ladder pair (direction 1): A_k << W_k << A_{{k+1}} ---");
    for i in 0..keys.len() {
        let k = keys[i];
        let a_k = a_polys.get(&k).unwrap();
        let w_k = w_polys.get(&k).unwrap();
        let r1 = check_interlacing(a_k, w_k);
        println!("  A_{} << W_{}: {}", k, k, r1);

        if i + 1 < keys.len() {
            let k_next = keys[i + 1];
            let a_next = a_polys.get(&k_next).unwrap();
            let r2 = check_interlacing(w_k, a_next);
            println!("  W_{} << A_{}: {}", k, k_next, r2);
        }
    }

    println!("\n--- Ladder pair (direction 2): A_{{k+1}} << W_k << A_k ---");
    for i in 0..keys.len() {
        let k = keys[i];
        let a_k = a_polys.get(&k).unwrap();
        let w_k = w_polys.get(&k).unwrap();
        let r1 = check_interlacing(w_k, a_k);
        println!("  W_{} << A_{}: {}", k, k, r1);

        if i + 1 < keys.len() {
            let k_next = keys[i + 1];
            let a_next = a_polys.get(&k_next).unwrap();
            let r2 = check_interlacing(a_next, w_k);
            println!("  A_{} << W_{}: {}", k_next, k, r2);
        }
    }

    // Also check W_k << W_{k+1} and A_k << W_k for the "common interlacing" property
    println!("\n--- Cross: W_k << A_{{k+1}} and A_k << W_{{k+1}} ---");
    for i in 0..keys.len() {
        if i + 1 < keys.len() {
            let k = keys[i];
            let k_next = keys[i + 1];
            let a_k = a_polys.get(&k).unwrap();
            let w_k = w_polys.get(&k).unwrap();
            let a_next = a_polys.get(&k_next).unwrap();
            let w_next = w_polys.get(&k_next).unwrap();

            let r1 = check_interlacing(a_k, w_next);
            println!("  A_{} << W_{}: {}", k, k_next, r1);
            let r2 = check_interlacing(w_k, a_next);
            println!("  W_{} << A_{}: {}", k, k_next, r2);
        }
    }

    // Also check the recursion:
    // A_k = D_k + U_k  (sum over all pi with pi_1=k)
    // W_k = t*D_k + U_k
    // The recursion relates lambda+ to lambda.
    // For reference, print the total peak polynomial.
    let total = interval.iter().fold(vec![0i64], |acc, pi| {
        let peaks = count_peaks(pi);
        let mut p = vec![0i64; peaks + 1];
        p[peaks] = 1;
        poly_add(&acc, &p)
    });
    println!("\n--- Total peak polynomial P^lambda(t) = sum_k A_k ---");
    println!("  P(t) = {}", format_poly(&total));
    println!("  Real-rooted: {}", if is_real_rooted(&total) { "YES" } else { "NO" });

    // Print roots of each polynomial for inspection
    println!("\n--- Roots (approximate) ---");
    for &k in &all_keys {
        let d = d_polys.get(&k).map_or(vec![0], |v| v.clone());
        let u = u_polys.get(&k).map_or(vec![0], |v| v.clone());
        if !poly_is_zero(&d) {
            if let Some(roots) = find_real_roots(&d) {
                println!("  D_{} roots: {:?}", k, roots.iter().map(|r| format!("{:.4}", r)).collect::<Vec<_>>());
            } else {
                println!("  D_{} roots: NOT ALL REAL", k);
            }
        }
        if !poly_is_zero(&u) {
            if let Some(roots) = find_real_roots(&u) {
                println!("  U_{} roots: {:?}", k, roots.iter().map(|r| format!("{:.4}", r)).collect::<Vec<_>>());
            } else {
                println!("  U_{} roots: NOT ALL REAL", k);
            }
        }
    }

    println!();
}
