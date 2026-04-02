//! Priority-queue search for GT polytopes with negative Ehrhart coefficients.
//!
//! Explores (λ/μ, w) triples in dimension range [MIN_DIM, MAX_DIM] using
//! best-first search guided by the normalized linear coefficient d! · c₁.
//! Smaller (more negative) values are explored first.
//!
//! Records are saved to `records.jsonl` and loaded on restart.

use hashbrown::HashSet;
use kostka::ehrhart::{compute_ehrhart, EhrhartPoly};
use kostka::gt_dim::gt_polytope_dim;
use kostka::partition::Partition;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, ToPrimitive, Zero};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::io::Write;
use std::time::Instant;

const MIN_DIM: usize = 18;
const MAX_DIM: usize = 23;
const MAX_TOTAL_SIZE: u32 = 30;
const MAX_PART: u32 = 12;
const MAX_STATES: usize = 5_000_000;

fn records_path() -> std::path::PathBuf {
    // Place records.jsonl next to our Cargo.toml.
    let mut p = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("records.jsonl");
    p
}

// ── Search entry ────────────────────────────────────────────────────────────

#[derive(Clone)]
struct SearchEntry {
    lambda: Vec<u32>,
    mu: Vec<u32>,
    weight: Vec<u32>,
    /// Parent's raw linear coefficient c₁.  Used as priority (smaller first).
    parent_score: BigRational,
}

impl Eq for SearchEntry {}
impl PartialEq for SearchEntry {
    fn eq(&self, other: &Self) -> bool {
        self.parent_score == other.parent_score
    }
}
impl Ord for SearchEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        other.parent_score.cmp(&self.parent_score) // min-heap
    }
}
impl PartialOrd for SearchEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// ── Canonical key ───────────────────────────────────────────────────────────

/// Canonical form: λ sorted descending (already), μ sorted descending (already),
/// w sorted descending (K is symmetric in w for non-flagged).
fn canonical_key(lambda: &[u32], mu: &[u32], weight: &[u32]) -> (Vec<u32>, Vec<u32>, Vec<u32>) {
    let mut w = weight.to_vec();
    w.sort_unstable_by(|a, b| b.cmp(a));
    while w.last() == Some(&0) {
        w.pop();
    }
    let mut m = mu.to_vec();
    while m.last() == Some(&0) {
        m.pop();
    }
    (lambda.to_vec(), m, w)
}

// ── Score extraction ────────────────────────────────────────────────────────

/// Raw rational linear coefficient c₁ of the Ehrhart polynomial.
fn linear_coeff(poly: &EhrhartPoly) -> Option<BigRational> {
    if poly.degree == 0 {
        return None;
    }
    Some(poly.coeffs[1].clone())
}

/// d! * c_k (integer) for display of the factored polynomial.
fn normalized_coeff(poly: &EhrhartPoly, k: usize) -> BigInt {
    if k > poly.degree {
        return BigInt::zero();
    }
    let d_fact_r = BigRational::from(factorial(poly.degree));
    (&poly.coeffs[k] * &d_fact_r).to_integer()
}

fn factorial(n: usize) -> BigInt {
    (1..=n as u64).fold(BigInt::one(), |acc, i| acc * BigInt::from(i))
}

/// Format a BigRational as a short decimal approximation for display.
fn fmt_rational(r: &BigRational) -> String {
    let num_f = r.numer().to_f64().unwrap_or(f64::INFINITY);
    let den_f = r.denom().to_f64().unwrap_or(1.0);
    let f = num_f / den_f;
    if f.is_finite() {
        format!("{:.6}", f)
    } else {
        format!("{}/{}", r.numer(), r.denom())
    }
}

// ── Formatting ──────────────────────────────────────────────────────────────

fn fmt_p(p: &[u32]) -> String {
    if p.is_empty() {
        "∅".to_string()
    } else {
        p.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
    }
}

fn fmt_shape(lam: &[u32], mu: &[u32]) -> String {
    if mu.is_empty() {
        format!("({})", fmt_p(lam))
    } else {
        format!("({})/({})", fmt_p(lam), fmt_p(mu))
    }
}

// ── Containment check ───────────────────────────────────────────────────────

fn mu_contained_in_lambda(mu: &[u32], lambda: &[u32]) -> bool {
    for (i, &m) in mu.iter().enumerate() {
        let l = lambda.get(i).copied().unwrap_or(0);
        if m > l {
            return false;
        }
    }
    true
}

// ── Mutations ───────────────────────────────────────────────────────────────

fn generate_mutations(
    lambda: &[u32],
    mu: &[u32],
    weight: &[u32],
) -> Vec<(Vec<u32>, Vec<u32>, Vec<u32>)> {
    let mut results = Vec::new();
    let n = lambda.len();
    let p = mu.len();
    let m = weight.len();

    // ── A. Reshape λ (keep |λ|, μ, w fixed) ────────────────────────────────
    for i in 0..n {
        for j in 0..n {
            if i == j || lambda[j] == 0 {
                continue;
            }
            let mut l2 = lambda.to_vec();
            l2[i] += 1;
            l2[j] -= 1;
            l2.sort_unstable_by(|a, b| b.cmp(a));
            trim_zeros(&mut l2);
            if l2 != lambda {
                try_add(&mut results, l2, mu.to_vec(), weight.to_vec());
            }
        }
    }

    // ── B. Reshape μ (keep |μ|, λ, w fixed) ────────────────────────────────
    if p >= 2 {
        for i in 0..p {
            for j in 0..p {
                if i == j || mu[j] == 0 {
                    continue;
                }
                let mut m2 = mu.to_vec();
                m2[i] += 1;
                m2[j] -= 1;
                m2.sort_unstable_by(|a, b| b.cmp(a));
                trim_zeros(&mut m2);
                if m2 != mu {
                    try_add(&mut results, lambda.to_vec(), m2, weight.to_vec());
                }
            }
        }
    }

    // ── C. Reshape w (keep |w|, λ, μ fixed) ────────────────────────────────
    for i in 0..m {
        for j in 0..m {
            if i == j || weight[j] == 0 {
                continue;
            }
            let mut w2 = weight.to_vec();
            w2[i] += 1;
            w2[j] -= 1;
            w2.sort_unstable_by(|a, b| b.cmp(a));
            trim_zeros(&mut w2);
            if w2 != weight {
                try_add(&mut results, lambda.to_vec(), mu.to_vec(), w2);
            }
        }
    }

    // ── D. Grow skew size by 1: λ_i += 1, w_j += 1 ────────────────────────
    for i in 0..n {
        let mut l2 = lambda.to_vec();
        l2[i] += 1;
        l2.sort_unstable_by(|a, b| b.cmp(a));

        for j in 0..m {
            let mut w2 = weight.to_vec();
            w2[j] += 1;
            w2.sort_unstable_by(|a, b| b.cmp(a));
            try_add(&mut results, l2.clone(), mu.to_vec(), w2);
        }
        // Add new w part = 1.
        let mut w2 = weight.to_vec();
        w2.push(1);
        w2.sort_unstable_by(|a, b| b.cmp(a));
        try_add(&mut results, l2, mu.to_vec(), w2);
    }
    // Add new λ part = 1.
    {
        let mut l2 = lambda.to_vec();
        l2.push(1);
        for j in 0..m {
            let mut w2 = weight.to_vec();
            w2[j] += 1;
            w2.sort_unstable_by(|a, b| b.cmp(a));
            try_add(&mut results, l2.clone(), mu.to_vec(), w2);
        }
        let mut w2 = weight.to_vec();
        w2.push(1);
        w2.sort_unstable_by(|a, b| b.cmp(a));
        try_add(&mut results, l2, mu.to_vec(), w2);
    }

    // ── E. Grow skew size by 1: μ_i -= 1, w_j += 1 ────────────────────────
    for i in 0..p {
        if mu[i] == 0 {
            continue;
        }
        let mut m2 = mu.to_vec();
        m2[i] -= 1;
        m2.sort_unstable_by(|a, b| b.cmp(a));
        trim_zeros(&mut m2);

        for j in 0..m {
            let mut w2 = weight.to_vec();
            w2[j] += 1;
            w2.sort_unstable_by(|a, b| b.cmp(a));
            try_add(&mut results, lambda.to_vec(), m2.clone(), w2);
        }
        let mut w2 = weight.to_vec();
        w2.push(1);
        w2.sort_unstable_by(|a, b| b.cmp(a));
        try_add(&mut results, lambda.to_vec(), m2, w2);
    }

    // ── F. Shrink skew size by 1: λ_i -= 1, w_j -= 1 ──────────────────────
    for i in 0..n {
        if lambda[i] == 0 {
            continue;
        }
        let mut l2 = lambda.to_vec();
        l2[i] -= 1;
        l2.sort_unstable_by(|a, b| b.cmp(a));
        trim_zeros(&mut l2);
        if l2.is_empty() {
            continue;
        }

        for j in 0..m {
            if weight[j] == 0 {
                continue;
            }
            let mut w2 = weight.to_vec();
            w2[j] -= 1;
            w2.sort_unstable_by(|a, b| b.cmp(a));
            trim_zeros(&mut w2);
            if w2.is_empty() {
                continue;
            }
            try_add(&mut results, l2.clone(), mu.to_vec(), w2);
        }
    }

    // ── G. Shrink skew size by 1: μ_i += 1, w_j -= 1 ──────────────────────
    for j in 0..m {
        if weight[j] == 0 {
            continue;
        }
        let mut w2 = weight.to_vec();
        w2[j] -= 1;
        w2.sort_unstable_by(|a, b| b.cmp(a));
        trim_zeros(&mut w2);
        if w2.is_empty() {
            continue;
        }

        // Increment existing μ part.
        for i in 0..p {
            let mut m2 = mu.to_vec();
            m2[i] += 1;
            m2.sort_unstable_by(|a, b| b.cmp(a));
            try_add(&mut results, lambda.to_vec(), m2, w2.clone());
        }
        // Add new μ part = 1.
        let mut m2 = mu.to_vec();
        m2.push(1);
        m2.sort_unstable_by(|a, b| b.cmp(a));
        try_add(&mut results, lambda.to_vec(), m2, w2);
    }

    results
}

fn trim_zeros(v: &mut Vec<u32>) {
    while v.last() == Some(&0) {
        v.pop();
    }
}

fn try_add(
    results: &mut Vec<(Vec<u32>, Vec<u32>, Vec<u32>)>,
    lambda: Vec<u32>,
    mu: Vec<u32>,
    weight: Vec<u32>,
) {
    if lambda.is_empty() || weight.is_empty() {
        return;
    }
    let lam_sum: u32 = lambda.iter().sum();
    let mu_sum: u32 = mu.iter().sum();
    let w_sum: u32 = weight.iter().sum();
    if lam_sum < mu_sum || lam_sum - mu_sum != w_sum {
        return;
    }
    if lam_sum > MAX_TOTAL_SIZE {
        return;
    }
    if lambda.iter().any(|&p| p > MAX_PART) || weight.iter().any(|&p| p > MAX_PART) {
        return;
    }
    if !mu_contained_in_lambda(&mu, &lambda) {
        return;
    }
    if let Some(d) = gt_polytope_dim(&lambda, &mu, &weight) {
        if d >= MIN_DIM && d <= MAX_DIM {
            results.push((lambda, mu, weight));
        }
    }
}

// ── Records persistence ─────────────────────────────────────────────────────

fn save_record(
    lam: &[u32],
    mu: &[u32],
    w: &[u32],
    degree: usize,
    score: &BigRational,
    poly_str: &str,
    has_negative: bool,
) {
    let line = format!(
        "{{\"lambda\":[{}],\"mu\":[{}],\"weight\":[{}],\"degree\":{},\"c1_num\":\"{}\",\"c1_den\":\"{}\",\"c1_approx\":{},\"negative\":{},\"polynomial\":\"{}\"}}\n",
        fmt_p(lam), fmt_p(mu), fmt_p(w), degree,
        score.numer(), score.denom(), fmt_rational(score),
        has_negative, poly_str,
    );
    if let Ok(mut f) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(records_path())
    {
        let _ = f.write_all(line.as_bytes());
    }
}

fn load_best_score() -> Option<BigRational> {
    let content = std::fs::read_to_string(records_path()).ok()?;
    let mut best: Option<BigRational> = None;
    for line in content.lines() {
        let num_str = extract_json_str(line, "c1_num")?;
        let den_str = extract_json_str(line, "c1_den")?;
        let num: BigInt = num_str.parse().ok()?;
        let den: BigInt = den_str.parse().ok()?;
        if den.is_zero() {
            continue;
        }
        let val = BigRational::new(num, den);
        let is_better = match &best {
            None => true,
            Some(b) => val < *b,
        };
        if is_better {
            best = Some(val);
        }
    }
    best
}

fn extract_json_str<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let needle = format!("\"{}\":\"", key);
    let start = line.find(&needle)? + needle.len();
    let rest = &line[start..];
    let end = rest.find('"')?;
    Some(&rest[..end])
}

// ── Seeds ───────────────────────────────────────────────────────────────────

fn seed_states() -> Vec<(Vec<u32>, Vec<u32>, Vec<u32>)> {
    let mut seeds: Vec<(Vec<u32>, Vec<u32>, Vec<u32>)> = Vec::new();

    // ── Non-skew seeds ──────────────────────────────────────────────────────
    let straight_lambdas: Vec<Vec<u32>> = vec![
        vec![3, 3, 2, 1, 1],
        vec![4, 3, 2, 1, 1],
        vec![4, 3, 2, 2, 1],
        vec![4, 3, 3, 2, 1],
        vec![4, 4, 3, 2, 1],
        vec![5, 4, 3, 2, 1],
        vec![3, 3, 2, 2, 1, 1],
        vec![3, 2, 2, 1, 1],
        vec![4, 3, 2, 1],
        vec![6, 5, 4, 3, 2, 1],
        vec![5, 4, 3, 2, 1, 1],
        vec![6, 5, 4, 3, 2],
        vec![7, 5, 4, 3, 2, 1],
        vec![4, 4, 4, 3, 3, 2, 1],
        vec![3, 3, 3, 3, 2, 2, 1, 1],
        vec![5, 4, 4, 3, 2, 1, 1],
    ];

    for lam in &straight_lambdas {
        let s: u32 = lam.iter().sum();
        let mu = vec![];

        // Unit weight.
        let w_unit: Vec<u32> = vec![1; s as usize];
        try_seed(&mut seeds, lam.clone(), mu.clone(), w_unit);

        // Weight with some 2s.
        if s >= 4 {
            let mut w2 = vec![2, 2];
            for _ in 0..(s - 4) {
                w2.push(1);
            }
            w2.sort_unstable_by(|a, b| b.cmp(a));
            try_seed(&mut seeds, lam.clone(), mu.clone(), w2);
        }
    }

    // ── Skew seeds ──────────────────────────────────────────────────────────
    // For each straight seed λ, try small inner shapes μ.
    let skew_pairs: Vec<(Vec<u32>, Vec<u32>)> = vec![
        (vec![5, 4, 3, 2, 1], vec![1]),
        (vec![5, 4, 3, 2, 1], vec![2, 1]),
        (vec![5, 4, 3, 2, 1], vec![3, 2, 1]),
        (vec![6, 5, 4, 3, 2, 1], vec![1]),
        (vec![6, 5, 4, 3, 2, 1], vec![2, 1]),
        (vec![6, 5, 4, 3, 2, 1], vec![3, 2, 1]),
        (vec![6, 5, 4, 3, 2, 1], vec![4, 3, 2, 1]),
        (vec![7, 5, 4, 3, 2, 1], vec![2, 1]),
        (vec![7, 5, 4, 3, 2, 1], vec![3, 2, 1]),
        (vec![7, 6, 5, 4, 3, 2, 1], vec![2, 1]),
        (vec![7, 6, 5, 4, 3, 2, 1], vec![3, 2, 1]),
        (vec![7, 6, 5, 4, 3, 2, 1], vec![4, 3, 2, 1]),
        (vec![5, 4, 4, 3, 2, 1, 1], vec![2, 1]),
        (vec![5, 4, 4, 3, 2, 1, 1], vec![3, 2, 1]),
        (vec![4, 4, 3, 2, 1], vec![1]),
        (vec![4, 4, 3, 2, 1], vec![2, 1]),
        (vec![6, 5, 4, 3, 2], vec![2, 1]),
        (vec![6, 5, 4, 3, 2], vec![3, 2]),
        (vec![8, 6, 5, 4, 3, 2, 1], vec![3, 2, 1]),
        (vec![8, 7, 6, 5, 4, 3, 2, 1], vec![4, 3, 2, 1]),
    ];

    for (lam, mu) in &skew_pairs {
        let s: u32 = lam.iter().sum::<u32>() - mu.iter().sum::<u32>();
        let w_unit: Vec<u32> = vec![1; s as usize];
        try_seed(&mut seeds, lam.clone(), mu.clone(), w_unit);

        if s >= 4 {
            let mut w2 = vec![2, 2];
            for _ in 0..(s - 4) {
                w2.push(1);
            }
            w2.sort_unstable_by(|a, b| b.cmp(a));
            try_seed(&mut seeds, lam.clone(), mu.clone(), w2);
        }
    }

    // Deduplicate.
    let mut seen: HashSet<(Vec<u32>, Vec<u32>, Vec<u32>)> = HashSet::new();
    seeds.retain(|(l, m, w)| {
        let key = canonical_key(l, m, w);
        seen.insert(key)
    });

    seeds
}

fn try_seed(
    seeds: &mut Vec<(Vec<u32>, Vec<u32>, Vec<u32>)>,
    lam: Vec<u32>,
    mu: Vec<u32>,
    w: Vec<u32>,
) {
    if !mu_contained_in_lambda(&mu, &lam) {
        return;
    }
    let lam_sum: u32 = lam.iter().sum();
    if lam_sum > MAX_TOTAL_SIZE || lam.iter().any(|&p| p > MAX_PART) {
        return;
    }
    if let Some(d) = gt_polytope_dim(&lam, &mu, &w) {
        if d >= MIN_DIM && d <= MAX_DIM {
            seeds.push((lam, mu, w));
        }
    }
}

// ── Main ────────────────────────────────────────────────────────────────────

fn main() {
    let mut heap: BinaryHeap<SearchEntry> = BinaryHeap::new();
    let mut visited: HashSet<(Vec<u32>, Vec<u32>, Vec<u32>)> = HashSet::new();

    // Load previous best from records file.
    let mut best_score: Option<BigRational> = load_best_score();
    let mut best_shape = String::new();
    let mut best_degree: usize = 0;
    let mut best_poly_str = String::new();
    let mut found_negative = false;
    let mut computed_count: u64 = 0;
    let start = Instant::now();

    if let Some(ref bs) = best_score {
        eprintln!("Loaded previous best c₁ ≈ {} from {}", fmt_rational(bs), records_path().display());
    }

    // Seed.
    let seeds = seed_states();
    eprintln!("Seeded with {} initial (λ/μ, w) triples", seeds.len());
    for (lam, mu, w) in seeds {
        let key = canonical_key(&lam, &mu, &w);
        if visited.insert(key) {
            heap.push(SearchEntry {
                lambda: lam,
                mu,
                weight: w,
                parent_score: BigRational::zero(),
            });
        }
    }

    // Main loop.
    while let Some(entry) = heap.pop() {
        let lam = &entry.lambda;
        let mu = &entry.mu;
        let w = &entry.weight;

        let shape = fmt_shape(lam, mu);
        eprint!(
            "\r\x1b[2K[{:>6} | Q:{:>7} | V:{:>7} | {:.0}s] {}  w=({})",
            computed_count,
            heap.len(),
            visited.len(),
            start.elapsed().as_secs_f64(),
            shape,
            fmt_p(w),
        );
        let _ = std::io::stderr().flush();

        let lam_p = Partition::new(lam.clone());
        let mu_p = if mu.is_empty() {
            Partition::empty()
        } else {
            Partition::new(mu.clone())
        };

        let poly = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            compute_ehrhart(&lam_p, &mu_p, w, None, None, false, Some(MAX_STATES), true)
        }));

        computed_count += 1;

        let poly = match poly {
            Ok(p) => p,
            Err(_) => continue,
        };

        if poly.degree < MIN_DIM {
            continue;
        }

        let score = match linear_coeff(&poly) {
            Some(s) => s,
            None => continue,
        };

        // Check for negative coefficients.
        let has_negative = poly.has_negative_coefficient();
        if has_negative {
            found_negative = true;
            eprintln!();
            eprintln!("════════════════════════════════════════════════════════════════");
            eprintln!("  NEGATIVE COEFFICIENT FOUND!");
            eprintln!("  shape = {}", shape);
            eprintln!("  w = ({})", fmt_p(w));
            eprintln!("  degree = {}", poly.degree);
            eprintln!("  polynomial = {}", poly.display_factored());
            eprintln!("  normalized coefficients (d! · c_k):");
            for k in 0..=poly.degree {
                let ck = normalized_coeff(&poly, k);
                let marker = if ck < BigInt::zero() {
                    " ← NEGATIVE"
                } else {
                    ""
                };
                eprintln!("    c_{} = {}{}", k, ck, marker);
            }
            eprintln!("════════════════════════════════════════════════════════════════");
            save_record(lam, mu, w, poly.degree, &score, &poly.display_factored(), true);
        }

        // Update best record.
        let is_new_best = match &best_score {
            None => true,
            Some(bs) => score < *bs,
        };
        if is_new_best {
            best_score = Some(score.clone());
            best_shape = shape.clone();
            best_degree = poly.degree;
            best_poly_str = poly.display_factored();

            eprintln!();
            eprintln!(
                "  ★ New record: c₁ ≈ {}  (deg {})  {}  w=({})",
                fmt_rational(&score),
                poly.degree,
                best_shape,
                fmt_p(w),
            );
            eprintln!("    P(n) = {}", best_poly_str);

            save_record(lam, mu, w, poly.degree, &score, &best_poly_str, has_negative);
        }

        // Enqueue mutations.
        let mutations = generate_mutations(lam, mu, w);
        for (ml, mm, mw) in mutations {
            let key = canonical_key(&ml, &mm, &mw);
            if visited.insert(key) {
                heap.push(SearchEntry {
                    lambda: ml,
                    mu: mm,
                    weight: mw,
                    parent_score: score.clone(),
                });
            }
        }
    }

    // Summary.
    eprintln!();
    eprintln!("═══════════════════════════ SEARCH COMPLETE ═══════════════════════════");
    eprintln!("  Computed: {}", computed_count);
    eprintln!("  Visited:  {}", visited.len());
    eprintln!("  Elapsed:  {:.1}s", start.elapsed().as_secs_f64());
    if found_negative {
        eprintln!("  Result:   NEGATIVE COEFFICIENTS FOUND (see above)");
    } else {
        eprintln!("  Result:   No negative coefficients found");
    }
    if let Some(bs) = &best_score {
        eprintln!("  Best c₁ ≈ {}  (deg {})  {}", fmt_rational(bs), best_degree, best_shape);
        eprintln!("    P(n) = {}", best_poly_str);
    }
    eprintln!("═══════════════════════════════════════════════════════════════════════");
}
