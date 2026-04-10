//! Priority-queue search for GT polytopes with negative Ehrhart coefficients.
//!
//! Explores (λ/μ, w) triples in dimension range [MIN_DIM, MAX_DIM] using
//! best-first search guided by the normalized linear coefficient d! · c₁.
//! Smaller (more negative) values are explored first.
//!
//! Records are saved to `records.jsonl` and loaded on restart.

use clap::Parser;
use dotenv::from_path;
use hashbrown::HashSet;
use kostka::ehrhart::{compute_ehrhart, compute_hstar, EhrhartPoly};
use kostka::gt_dim::gt_polytope_dim_full;
use kostka::Partition;
use mysql::prelude::*;
use mysql::*;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, ToPrimitive, Zero};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::io::Write;
use std::path::Path;
use std::sync::OnceLock;
use std::time::Instant;

const DEFAULT_MIN_DIM: usize = 14;
const DEFAULT_MAX_DIM: usize = 22;
const DEFAULT_MAX_LAM_PART: u32 = 5;
const DEFAULT_MAX_W_PART: u32 = 5;
const DEFAULT_MAX_TOTAL_SIZE: u32 = 30;
const DEFAULT_MAX_STATES: usize = 200_000;

const CREATE_TABLE: &str = r"
CREATE TABLE IF NOT EXISTS ktt_search_records (
    id              BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    lambda          VARCHAR(255) NOT NULL,
    mu              VARCHAR(255) NOT NULL DEFAULT '',
    weight          VARCHAR(255) NOT NULL,
    degree          INT UNSIGNED NOT NULL,
    c1_num          TEXT NOT NULL,
    c1_den          TEXT NOT NULL,
    c1_approx       DOUBLE,
    negative        BOOLEAN NOT NULL DEFAULT FALSE,
    reason          VARCHAR(32) NOT NULL,
    polynomial      TEXT NOT NULL,
    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE KEY uq_ktt_record (lambda, mu, weight, reason)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
";

const CREATE_STATES_TABLE: &str = r"
CREATE TABLE IF NOT EXISTS ktt_search_states (
    id              BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    lambda          VARCHAR(255) NOT NULL,
    mu              VARCHAR(255) NOT NULL DEFAULT '',
    weight          VARCHAR(255) NOT NULL,
    degree          INT UNSIGNED,
    c1_num          TEXT,
    c1_den          TEXT,
    c1_approx       DOUBLE,
    negative        BOOLEAN NOT NULL DEFAULT FALSE,
    status          VARCHAR(32) NOT NULL,
    polynomial      TEXT,
    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    UNIQUE KEY uq_ktt_state (lambda, mu, weight)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
";

static CONFIG: OnceLock<SearchConfig> = OnceLock::new();

#[derive(Parser, Debug, Clone)]
#[command(
    name = "ktt-search",
    about = "Priority-queue search for GT polytopes with negative Ehrhart coefficients"
)]
struct Args {
    /// Minimum GT-polytope dimension to keep in the search.
    #[arg(long, default_value_t = DEFAULT_MIN_DIM)]
    min_dim: usize,

    /// Maximum GT-polytope dimension to keep in the search.
    #[arg(long, default_value_t = DEFAULT_MAX_DIM)]
    max_dim: usize,

    /// Maximum allowed part size in lambda.
    #[arg(long, default_value_t = DEFAULT_MAX_LAM_PART)]
    max_lam_part: u32,

    /// Maximum allowed part size in the weight.
    #[arg(long, default_value_t = DEFAULT_MAX_W_PART)]
    max_w_part: u32,

    /// Maximum total size |lambda|.
    #[arg(long, default_value_t = DEFAULT_MAX_TOTAL_SIZE)]
    max_total_size: u32,

    /// Abort a DP call if any level exceeds this many states.
    #[arg(long, default_value_t = DEFAULT_MAX_STATES)]
    max_states: usize,

    /// Allow skew shapes mu != empty during the search.
    #[arg(long)]
    skew: bool,

    /// MariaDB connection URL. Also read from KOSTKA_DB_URL if set.
    #[arg(long, env = "KOSTKA_DB_URL")]
    db_url: Option<String>,
}

#[derive(Debug, Clone)]
struct SearchConfig {
    min_dim: usize,
    max_dim: usize,
    max_lam_part: u32,
    max_w_part: u32,
    max_total_size: u32,
    max_states: usize,
    skew: bool,
}

impl From<&Args> for SearchConfig {
    fn from(args: &Args) -> Self {
        Self {
            min_dim: args.min_dim,
            max_dim: args.max_dim,
            max_lam_part: args.max_lam_part,
            max_w_part: args.max_w_part,
            max_total_size: args.max_total_size,
            max_states: args.max_states,
            skew: args.skew,
        }
    }
}

fn config() -> &'static SearchConfig {
    CONFIG.get().expect("search config not initialized")
}

fn records_path() -> std::path::PathBuf {
    // Place records.jsonl next to our Cargo.toml.
    let mut p = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("records.jsonl");
    p
}

fn load_db_url(cli_url: Option<String>) -> Option<String> {
    if cli_url.is_some() {
        return cli_url;
    }

    let env_path = Path::new("/home/paxinum/Dropbox/projects/rust/kostka/.env");
    let _ = from_path(env_path);

    let host = match std::env::var("DB_HOST") {
        Ok(host) if host == "localhost" => "127.0.0.1".to_string(),
        Ok(host) => host,
        Err(_) => "127.0.0.1".to_string(),
    };
    let user = std::env::var("DB_USER").ok()?;
    let password = std::env::var("DB_PASSWORD").ok()?;
    let db = std::env::var("DB_NAME").ok()?;

    Some(format!("mysql://{}:{}@{}/{}", user, password, host, db))
}

fn init_db_pool(db_url: Option<String>) -> Option<Pool> {
    let Some(db_url) = load_db_url(db_url) else {
        eprintln!("Database persistence disabled: no --db-url, KOSTKA_DB_URL, or kostka/.env");
        return None;
    };

    let opts = Opts::from_url(&db_url).unwrap_or_else(|e| {
        eprintln!("Invalid database URL: {}", e);
        std::process::exit(1);
    });
    let pool = Pool::new(opts).unwrap_or_else(|e| {
        eprintln!("Failed to connect to database: {}", e);
        std::process::exit(1);
    });
    let mut conn = pool.get_conn().expect("Failed to get DB connection");
    conn.query_drop(CREATE_TABLE)
        .expect("Failed to create ktt_search_records table");
    conn.query_drop(CREATE_STATES_TABLE)
        .expect("Failed to create ktt_search_states table");
    // Add flag columns if missing (idempotent migration).
    for tbl in ["ktt_search_records", "ktt_search_states"] {
        for col in ["upper_flags", "lower_flags"] {
            let _ = conn.query_drop(format!(
                "ALTER TABLE {} ADD COLUMN {} VARCHAR(255) NOT NULL DEFAULT ''",
                tbl, col
            ));
        }
    }
    // Widen the unique key to include flags (drop + re-add, ignore errors if already done).
    let _ = conn.query_drop(
        "ALTER TABLE ktt_search_states DROP INDEX uq_ktt_state, \
         ADD UNIQUE KEY uq_ktt_state (lambda, mu, weight, upper_flags, lower_flags)",
    );
    let _ = conn.query_drop(
        "ALTER TABLE ktt_search_records DROP INDEX uq_ktt_record, \
         ADD UNIQUE KEY uq_ktt_record (lambda, mu, weight, upper_flags, lower_flags, reason)",
    );
    eprintln!("Connected to database tables ktt_search_records and ktt_search_states.");
    Some(pool)
}

// ── Search entry ────────────────────────────────────────────────────────────

/// 5-tuple key: (λ, μ, w, upper_flags, lower_flags).
type SearchKey = (Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>);

#[derive(Clone)]
struct SearchEntry {
    lambda: Vec<u32>,
    mu: Vec<u32>,
    weight: Vec<u32>,
    upper_flags: Vec<u32>,
    lower_flags: Vec<u32>,
    /// Search priority inherited from the parent state. Smaller is explored first.
    priority: BigRational,
}

impl Eq for SearchEntry {}
impl PartialEq for SearchEntry {
    fn eq(&self, other: &Self) -> bool {
        self.priority == other.priority
    }
}
impl Ord for SearchEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        other.priority.cmp(&self.priority) // min-heap
    }
}
impl PartialOrd for SearchEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// ── Canonical key ───────────────────────────────────────────────────────────

/// Canonical form: λ sorted desc, μ sorted desc.
/// w sorted desc only when unflagged (K is symmetric in w for non-flagged).
/// Flags are positional (tied to w), so w order is fixed when flagged.
fn canonical_key(lambda: &[u32], mu: &[u32], weight: &[u32], uf: &[u32], lf: &[u32]) -> SearchKey {
    let flagged = !uf.is_empty() || !lf.is_empty();
    let mut w = weight.to_vec();
    if !flagged {
        w.sort_unstable_by(|a, b| b.cmp(a));
    }
    while w.last() == Some(&0) {
        w.pop();
    }
    let mut m = mu.to_vec();
    while m.last() == Some(&0) {
        m.pop();
    }
    (lambda.to_vec(), m, w, uf.to_vec(), lf.to_vec())
}

// ── Score extraction ────────────────────────────────────────────────────────

/// Score: c₁ − H_d (deviation from simplex baseline).
/// The simplex has c₁ = H_d at every degree, so this score is 0 for all simplices.
/// Degree-neutral: low and high dimensions compete equally.
/// Negative score means c₁ is below the simplex baseline (promising).
fn score_c1_deviation(poly: &EhrhartPoly) -> Option<BigRational> {
    if poly.degree == 0 {
        return None;
    }
    let h_d = harmonic(poly.degree);
    Some(&poly.coeffs[1] - &h_d)
}

/// H_d = 1 + 1/2 + ... + 1/d (exact rational).
fn harmonic(d: usize) -> BigRational {
    let mut h = BigRational::zero();
    for k in 1..=d as u64 {
        h += BigRational::new(BigInt::one(), BigInt::from(k));
    }
    h
}

fn run_repeat_score(parts: &[u32]) -> u64 {
    if parts.is_empty() {
        return 0;
    }

    let mut score = 0u64;
    let mut run_len = 1u64;
    for i in 1..parts.len() {
        if parts[i] == parts[i - 1] {
            run_len += 1;
        } else {
            if run_len >= 2 {
                score += run_len * run_len;
            }
            run_len = 1;
        }
    }
    if run_len >= 2 {
        score += run_len * run_len;
    }
    score
}

fn structural_complexity(lam: &[u32], mu: &[u32], w: &[u32]) -> u64 {
    run_repeat_score(lam) + run_repeat_score(mu) + run_repeat_score(w)
}

fn seed_priority(lam: &[u32], mu: &[u32], w: &[u32]) -> BigRational {
    let complexity = structural_complexity(lam, mu, w).max(1);
    BigRational::new(BigInt::one(), BigInt::from(complexity))
}

fn child_priority(
    parent_priority: &BigRational,
    lam: &[u32],
    mu: &[u32],
    w: &[u32],
) -> BigRational {
    let complexity = structural_complexity(lam, mu, w).max(1);
    parent_priority / BigRational::from(BigInt::from(complexity))
}

fn search_priority(
    raw_score: &BigRational,
    poly: &EhrhartPoly,
    _lam: &[u32],
    _mu: &[u32],
    _w: &[u32],
) -> BigRational {
    // Use h*-center-of-mass to guide search.
    // Top-heavy h* → large COM → smaller priority → explored first in min-heap.
    // Simplex has h* = (1,0,...,0) → COM = 0 → priority = raw_score (no boost).
    // For negative c₁, we need large h*_i at high indices (where binomial basis
    // coefficients contribute negatively to the linear monomial term).
    let hstar = compute_hstar(poly);

    let total: BigInt = hstar.iter().sum();
    if total <= BigInt::zero() {
        return raw_score.clone();
    }
    let weighted: BigInt = hstar
        .iter()
        .enumerate()
        .map(|(i, h)| BigInt::from(i as u64) * h)
        .sum();
    // priority = raw_score / (h_com + 1) = raw_score * total / (weighted + total)
    let denom = &weighted + &total;
    if denom <= BigInt::zero() {
        return raw_score.clone();
    }

    raw_score * &BigRational::new(total, denom)
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
        p.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    }
}

fn csv_p(p: &[u32]) -> String {
    p.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn json_p(p: &[u32]) -> String {
    p.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
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
    uf: &[u32],
    lf: &[u32],
) -> Vec<(Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>)> {
    let mut results = Vec::new();
    let flagged = !uf.is_empty() || !lf.is_empty();
    let n = lambda.len();
    let p = mu.len();
    let m = weight.len();

    let no_f = vec![];

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
                try_add(
                    &mut results,
                    l2,
                    mu.to_vec(),
                    weight.to_vec(),
                    uf.to_vec(),
                    lf.to_vec(),
                );
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
                    try_add(
                        &mut results,
                        lambda.to_vec(),
                        m2,
                        weight.to_vec(),
                        uf.to_vec(),
                        lf.to_vec(),
                    );
                }
            }
        }
    }

    // ── C. Reshape w (keep |w|, λ, μ fixed; drops flags since w order changes)
    if !flagged {
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
                    try_add(
                        &mut results,
                        lambda.to_vec(),
                        mu.to_vec(),
                        w2,
                        no_f.clone(),
                        no_f.clone(),
                    );
                }
            }
        }
    }

    // ── D. Grow skew size by 1: λ_i += 1, w_j += 1 (drops flags) ──────────
    for i in 0..n {
        let mut l2 = lambda.to_vec();
        l2[i] += 1;
        l2.sort_unstable_by(|a, b| b.cmp(a));

        for j in 0..m {
            let mut w2 = weight.to_vec();
            w2[j] += 1;
            w2.sort_unstable_by(|a, b| b.cmp(a));
            try_add(
                &mut results,
                l2.clone(),
                mu.to_vec(),
                w2,
                no_f.clone(),
                no_f.clone(),
            );
        }
        let mut w2 = weight.to_vec();
        w2.push(1);
        w2.sort_unstable_by(|a, b| b.cmp(a));
        try_add(
            &mut results,
            l2,
            mu.to_vec(),
            w2,
            no_f.clone(),
            no_f.clone(),
        );
    }
    // Add new λ part = 1.
    {
        let mut l2 = lambda.to_vec();
        l2.push(1);
        for j in 0..m {
            let mut w2 = weight.to_vec();
            w2[j] += 1;
            w2.sort_unstable_by(|a, b| b.cmp(a));
            try_add(
                &mut results,
                l2.clone(),
                mu.to_vec(),
                w2,
                no_f.clone(),
                no_f.clone(),
            );
        }
        let mut w2 = weight.to_vec();
        w2.push(1);
        w2.sort_unstable_by(|a, b| b.cmp(a));
        try_add(
            &mut results,
            l2,
            mu.to_vec(),
            w2,
            no_f.clone(),
            no_f.clone(),
        );
    }

    // ── E. Grow skew size by 1: μ_i -= 1, w_j += 1 (drops flags) ──────────
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
            try_add(
                &mut results,
                lambda.to_vec(),
                m2.clone(),
                w2,
                no_f.clone(),
                no_f.clone(),
            );
        }
        let mut w2 = weight.to_vec();
        w2.push(1);
        w2.sort_unstable_by(|a, b| b.cmp(a));
        try_add(
            &mut results,
            lambda.to_vec(),
            m2,
            w2,
            no_f.clone(),
            no_f.clone(),
        );
    }

    // ── F. Shrink skew size by 1: λ_i -= 1, w_j -= 1 (drops flags) ────────
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
            try_add(
                &mut results,
                l2.clone(),
                mu.to_vec(),
                w2,
                no_f.clone(),
                no_f.clone(),
            );
        }
    }

    // ── G. Shrink skew size by 1: μ_i += 1, w_j -= 1 (drops flags) ────────
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

        for i in 0..p {
            let mut m2 = mu.to_vec();
            m2[i] += 1;
            m2.sort_unstable_by(|a, b| b.cmp(a));
            try_add(
                &mut results,
                lambda.to_vec(),
                m2,
                w2.clone(),
                no_f.clone(),
                no_f.clone(),
            );
        }
        let mut m2 = mu.to_vec();
        m2.push(1);
        m2.sort_unstable_by(|a, b| b.cmp(a));
        try_add(
            &mut results,
            lambda.to_vec(),
            m2,
            w2,
            no_f.clone(),
            no_f.clone(),
        );
    }

    // ── H. Flag mutations: shift upper/lower flags ±1, add/remove flags ────
    // Keep λ, μ, w fixed; only change flags.
    {
        let n_rows = n as u32;
        // Mutate existing upper flags.
        for i in 0..uf.len() {
            if uf[i] > 1 {
                let mut uf2 = uf.to_vec();
                uf2[i] -= 1;
                try_add(
                    &mut results,
                    lambda.to_vec(),
                    mu.to_vec(),
                    weight.to_vec(),
                    uf2,
                    lf.to_vec(),
                );
            }
            if uf[i] < n_rows {
                let mut uf2 = uf.to_vec();
                uf2[i] += 1;
                try_add(
                    &mut results,
                    lambda.to_vec(),
                    mu.to_vec(),
                    weight.to_vec(),
                    uf2,
                    lf.to_vec(),
                );
            }
            // Remove this flag (set to n = unrestricted).
            let mut uf2 = uf.to_vec();
            uf2[i] = n_rows;
            try_add(
                &mut results,
                lambda.to_vec(),
                mu.to_vec(),
                weight.to_vec(),
                uf2,
                lf.to_vec(),
            );
        }
        // Mutate existing lower flags.
        for i in 0..lf.len() {
            if lf[i] > 1 {
                let mut lf2 = lf.to_vec();
                lf2[i] -= 1;
                try_add(
                    &mut results,
                    lambda.to_vec(),
                    mu.to_vec(),
                    weight.to_vec(),
                    uf.to_vec(),
                    lf2,
                );
            }
            if lf[i] < n_rows {
                let mut lf2 = lf.to_vec();
                lf2[i] += 1;
                try_add(
                    &mut results,
                    lambda.to_vec(),
                    mu.to_vec(),
                    weight.to_vec(),
                    uf.to_vec(),
                    lf2,
                );
            }
        }
        // Add a new upper flag: restrict one weight step to rows 1..=r for r in {1, n/2, n-1}.
        if uf.is_empty() && n_rows >= 2 {
            for i in 0..m {
                for &r in &[1, n_rows / 2, n_rows - 1] {
                    if r == 0 || r >= n_rows {
                        continue;
                    }
                    let mut uf2 = vec![n_rows; m]; // all unrestricted
                    uf2[i] = r;
                    try_add(
                        &mut results,
                        lambda.to_vec(),
                        mu.to_vec(),
                        weight.to_vec(),
                        uf2,
                        no_f.clone(),
                    );
                }
            }
        }
    }

    results
}

fn trim_zeros(v: &mut Vec<u32>) {
    while v.last() == Some(&0) {
        v.pop();
    }
}

fn try_add(
    results: &mut Vec<(Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>)>,
    lambda: Vec<u32>,
    mu: Vec<u32>,
    weight: Vec<u32>,
    uf: Vec<u32>,
    lf: Vec<u32>,
) {
    let cfg = config();
    if lambda.is_empty() || weight.is_empty() {
        return;
    }
    if !cfg.skew && !mu.is_empty() {
        return;
    }
    let lam_sum: u32 = lambda.iter().sum();
    let mu_sum: u32 = mu.iter().sum();
    let w_sum: u32 = weight.iter().sum();
    if lam_sum < mu_sum || lam_sum - mu_sum != w_sum {
        return;
    }
    if lam_sum > cfg.max_total_size {
        return;
    }
    if lambda.iter().any(|&p| p > cfg.max_lam_part) || weight.iter().any(|&p| p > cfg.max_w_part) {
        return;
    }
    if !mu.is_empty() && !mu_contained_in_lambda(&mu, &lambda) {
        return;
    }
    let uf_opt = if uf.is_empty() {
        None
    } else {
        Some(uf.as_slice())
    };
    let lf_opt = if lf.is_empty() {
        None
    } else {
        Some(lf.as_slice())
    };
    if let Some(d) = gt_polytope_dim_full(&lambda, &mu, &weight, uf_opt, lf_opt) {
        if d >= cfg.min_dim && d <= cfg.max_dim {
            results.push((lambda, mu, weight, uf, lf));
        }
    }
}

// ── Records persistence ─────────────────────────────────────────────────────

fn save_record(
    lam: &[u32],
    mu: &[u32],
    w: &[u32],
    uf: &[u32],
    lf: &[u32],
    degree: usize,
    score: &BigRational,
    poly_str: &str,
    has_negative: bool,
    reason: &str,
    db_pool: Option<&Pool>,
) {
    let line = format!(
        "{{\"lambda\":[{}],\"mu\":[{}],\"weight\":[{}],\"upper_flags\":[{}],\"lower_flags\":[{}],\"degree\":{},\"c1_num\":\"{}\",\"c1_den\":\"{}\",\"c1_approx\":{},\"negative\":{},\"reason\":\"{}\",\"polynomial\":\"{}\"}}\n",
        json_p(lam), json_p(mu), json_p(w), json_p(uf), json_p(lf), degree,
        score.numer(), score.denom(), fmt_rational(score),
        has_negative, reason, poly_str,
    );
    if let Ok(mut f) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(records_path())
    {
        let _ = f.write_all(line.as_bytes());
    }

    if let Some(pool) = db_pool {
        let mut conn = pool.get_conn().expect("Failed to get DB connection");
        conn.exec_drop(
            r"INSERT IGNORE INTO ktt_search_records
              (lambda, mu, weight, upper_flags, lower_flags, degree, c1_num, c1_den, c1_approx, negative, reason, polynomial)
              VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                csv_p(lam),
                csv_p(mu),
                csv_p(w),
                csv_p(uf),
                csv_p(lf),
                degree as u32,
                score.numer().to_string(),
                score.denom().to_string(),
                score.numer().to_f64().unwrap_or(f64::INFINITY)
                    / score.denom().to_f64().unwrap_or(1.0),
                has_negative,
                reason,
                poly_str,
            ),
        )
        .unwrap_or_else(|e| eprintln!("DB insert error: {}", e));
    }
}

fn save_state(
    lam: &[u32],
    mu: &[u32],
    w: &[u32],
    uf: &[u32],
    lf: &[u32],
    degree: Option<usize>,
    score: Option<&BigRational>,
    poly_str: Option<&str>,
    has_negative: bool,
    status: &str,
    db_pool: Option<&Pool>,
) {
    let Some(pool) = db_pool else {
        return;
    };

    let c1_num = score.map(|s| s.numer().to_string());
    let c1_den = score.map(|s| s.denom().to_string());
    let c1_approx = score
        .map(|s| s.numer().to_f64().unwrap_or(f64::INFINITY) / s.denom().to_f64().unwrap_or(1.0));
    let degree_db = degree.map(|d| d as u32);

    let mut conn = pool.get_conn().expect("Failed to get DB connection");
    conn.exec_drop(
        r"INSERT INTO ktt_search_states
          (lambda, mu, weight, upper_flags, lower_flags, degree, c1_num, c1_den, c1_approx, negative, status, polynomial)
          VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
          ON DUPLICATE KEY UPDATE
            degree = VALUES(degree),
            c1_num = VALUES(c1_num),
            c1_den = VALUES(c1_den),
            c1_approx = VALUES(c1_approx),
            negative = VALUES(negative),
            status = VALUES(status),
            polynomial = VALUES(polynomial)",
        (
            csv_p(lam),
            csv_p(mu),
            csv_p(w),
            csv_p(uf),
            csv_p(lf),
            degree_db,
            c1_num,
            c1_den,
            c1_approx,
            has_negative,
            status,
            poly_str,
        ),
    )
    .unwrap_or_else(|e| eprintln!("DB state insert error: {}", e));
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

fn seed_states() -> Vec<(Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>)> {
    let cfg = config();
    let mut seeds: Vec<(Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>)> = Vec::new();
    let mu = vec![];
    let no_f = vec![];

    // Partitions with many repeated parts, max entry ≤ 5.
    // (v^a) for various v, a:
    for v in 2..=5u32 {
        for a in 3..=10usize {
            let lam = vec![v; a];
            let s = v * a as u32;
            if s > cfg.max_total_size {
                continue;
            }

            // Unit weight.
            let w_unit: Vec<u32> = vec![1; s as usize];
            try_seed(
                &mut seeds,
                lam.clone(),
                mu.clone(),
                w_unit,
                no_f.clone(),
                no_f.clone(),
            );

            // Weight with 2s.
            for n2 in 1..=(s / 2) as usize {
                let ones = s as usize - 2 * n2;
                let mut w: Vec<u32> = vec![2; n2];
                w.extend(vec![1u32; ones]);
                try_seed(
                    &mut seeds,
                    lam.clone(),
                    mu.clone(),
                    w,
                    no_f.clone(),
                    no_f.clone(),
                );
            }
        }
    }

    // (v^a, u^b) for v > u, many repeated:
    for v in 3..=5u32 {
        for u in 1..v {
            for a in 2..=8usize {
                for b in 2..=8usize {
                    let mut lam = vec![v; a];
                    lam.extend(vec![u; b]);
                    let s: u32 = lam.iter().sum();
                    if s > cfg.max_total_size {
                        continue;
                    }

                    let w_unit: Vec<u32> = vec![1; s as usize];
                    try_seed(
                        &mut seeds,
                        lam.clone(),
                        mu.clone(),
                        w_unit,
                        no_f.clone(),
                        no_f.clone(),
                    );

                    // Mixed weight with some 2s.
                    if s >= 6 {
                        let n2 = (s / 4) as usize;
                        let ones = s as usize - 2 * n2;
                        let mut w: Vec<u32> = vec![2; n2];
                        w.extend(vec![1u32; ones]);
                        try_seed(
                            &mut seeds,
                            lam.clone(),
                            mu.clone(),
                            w,
                            no_f.clone(),
                            no_f.clone(),
                        );
                    }

                    // Weight matching λ shape.
                    let w_match = lam.clone();
                    try_seed(
                        &mut seeds,
                        lam.clone(),
                        mu.clone(),
                        w_match,
                        no_f.clone(),
                        no_f.clone(),
                    );
                }
            }
        }
    }

    // Staircase shapes (top, top-1, ..., 1).
    for top in 3..=cfg.max_lam_part {
        let lam: Vec<u32> = (1..=top).rev().collect();
        let s: u32 = lam.iter().sum();
        if s <= cfg.max_total_size {
            let w_unit: Vec<u32> = vec![1; s as usize];
            try_seed(
                &mut seeds,
                lam.clone(),
                mu.clone(),
                w_unit,
                no_f.clone(),
                no_f.clone(),
            );
            for n2 in 1..=(s / 2) as usize {
                let ones = s as usize - 2 * n2;
                let mut w: Vec<u32> = vec![2; n2];
                w.extend(vec![1u32; ones]);
                try_seed(
                    &mut seeds,
                    lam.clone(),
                    mu.clone(),
                    w,
                    no_f.clone(),
                    no_f.clone(),
                );
            }
        }
        // Repeated staircases: (top, top-1, ..., 1)^r
        for r in 2..=4usize {
            let lam_r: Vec<u32> = (0..r).flat_map(|_| (1..=top).rev()).collect();
            let mut lam_r_sorted = lam_r;
            lam_r_sorted.sort_unstable_by(|a, b| b.cmp(a));
            let s: u32 = lam_r_sorted.iter().sum();
            if s <= cfg.max_total_size {
                let w_unit: Vec<u32> = vec![1; s as usize];
                try_seed(
                    &mut seeds,
                    lam_r_sorted.clone(),
                    mu.clone(),
                    w_unit,
                    no_f.clone(),
                    no_f.clone(),
                );
            }
        }
    }

    // Hook shapes (v, 1^a).
    for v in 3..=cfg.max_lam_part {
        for a in 2..=15usize {
            let mut lam = vec![v];
            lam.extend(vec![1u32; a]);
            let s: u32 = lam.iter().sum();
            if s > cfg.max_total_size {
                continue;
            }
            let w_unit: Vec<u32> = vec![1; s as usize];
            try_seed(
                &mut seeds,
                lam.clone(),
                mu.clone(),
                w_unit,
                no_f.clone(),
                no_f.clone(),
            );
            if s >= 6 {
                let n2 = (s / 4) as usize;
                let ones = s as usize - 2 * n2;
                let mut w: Vec<u32> = vec![2; n2];
                w.extend(vec![1u32; ones]);
                try_seed(&mut seeds, lam, mu.clone(), w, no_f.clone(), no_f.clone());
            }
        }
    }

    // Near-staircase (top, top-1, ..., 1, 1, ..., 1).
    for top in 3..=cfg.max_lam_part {
        for extra in 1..=8usize {
            let mut lam: Vec<u32> = (1..=top).rev().collect();
            lam.extend(vec![1u32; extra]);
            let s: u32 = lam.iter().sum();
            if s > cfg.max_total_size {
                continue;
            }
            let w_unit: Vec<u32> = vec![1; s as usize];
            try_seed(
                &mut seeds,
                lam,
                mu.clone(),
                w_unit,
                no_f.clone(),
                no_f.clone(),
            );
        }
    }

    // Skew seeds: λ/μ pairs (only when --skew is enabled).
    if cfg.skew {
        for v in 3..=cfg.max_lam_part {
            for a in 3..=8usize {
                let lam = vec![v; a];
                let lam_sum = v * a as u32;
                // μ = (u^b) for u < v, b < a.
                for u in 1..v {
                    for b in 1..a {
                        let skew_mu = vec![u; b];
                        let mu_sum = u * b as u32;
                        if mu_sum >= lam_sum {
                            continue;
                        }
                        let skew_size = lam_sum - mu_sum;
                        let w_unit: Vec<u32> = vec![1; skew_size as usize];
                        try_seed(
                            &mut seeds,
                            lam.clone(),
                            skew_mu.clone(),
                            w_unit,
                            no_f.clone(),
                            no_f.clone(),
                        );
                        if skew_size >= 6 {
                            let n2 = (skew_size / 4) as usize;
                            let ones = skew_size as usize - 2 * n2;
                            let mut w: Vec<u32> = vec![2; n2];
                            w.extend(vec![1u32; ones]);
                            try_seed(
                                &mut seeds,
                                lam.clone(),
                                skew_mu,
                                w,
                                no_f.clone(),
                                no_f.clone(),
                            );
                        }
                    }
                }
            }
        }
    }

    // Deduplicate.
    let mut seen: HashSet<SearchKey> = HashSet::new();
    seeds.retain(|(l, m, w, uf, lf)| {
        let key = canonical_key(l, m, w, uf, lf);
        seen.insert(key)
    });

    seeds
}

fn try_seed(
    seeds: &mut Vec<(Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>)>,
    lam: Vec<u32>,
    mu: Vec<u32>,
    w: Vec<u32>,
    uf: Vec<u32>,
    lf: Vec<u32>,
) {
    let cfg = config();
    if !mu_contained_in_lambda(&mu, &lam) {
        return;
    }
    let lam_sum: u32 = lam.iter().sum();
    if lam_sum > cfg.max_total_size || lam.iter().any(|&p| p > cfg.max_lam_part) {
        return;
    }
    let uf_opt = if uf.is_empty() {
        None
    } else {
        Some(uf.as_slice())
    };
    let lf_opt = if lf.is_empty() {
        None
    } else {
        Some(lf.as_slice())
    };
    if let Some(d) = gt_polytope_dim_full(&lam, &mu, &w, uf_opt, lf_opt) {
        if d >= cfg.min_dim && d <= cfg.max_dim {
            seeds.push((lam, mu, w, uf, lf));
        }
    }
}

// ── Main ────────────────────────────────────────────────────────────────────

fn main() {
    let args = Args::parse();
    if args.min_dim > args.max_dim {
        eprintln!("Invalid search bounds: --min-dim must be <= --max-dim");
        std::process::exit(2);
    }
    CONFIG
        .set(SearchConfig::from(&args))
        .expect("search config already initialized");
    let db_pool = init_db_pool(args.db_url.clone());

    let mut heap: BinaryHeap<SearchEntry> = BinaryHeap::new();
    let mut visited: HashSet<SearchKey> = HashSet::new();
    let mut seen_polys: HashSet<String> = HashSet::new();

    // Load previous best from records file.
    let mut best_score: Option<BigRational> = load_best_score();
    let mut best_shape = String::new();
    let mut best_degree: usize = 0;
    let mut best_poly_str = String::new();
    let mut found_negative = false;
    let mut computed_count: u64 = 0;
    let mut dedup_count: u64 = 0;
    let start = Instant::now();
    let cfg = config();

    eprintln!(
        "Search config: dim {}..={}, max_lam_part={}, max_w_part={}, max_total_size={}, max_states={}, skew={}",
        cfg.min_dim,
        cfg.max_dim,
        cfg.max_lam_part,
        cfg.max_w_part,
        cfg.max_total_size,
        cfg.max_states,
        cfg.skew,
    );

    if let Some(ref bs) = best_score {
        eprintln!(
            "Loaded previous best c₁−H_d ≈ {} from {}",
            fmt_rational(bs),
            records_path().display()
        );
    }

    // Seed.
    let seeds = seed_states();
    eprintln!("Seeded with {} initial (λ/μ, w) triples", seeds.len());
    for (lam, mu, w, uf, lf) in seeds {
        let priority = seed_priority(&lam, &mu, &w);
        let key = canonical_key(&lam, &mu, &w, &uf, &lf);
        if visited.insert(key) {
            heap.push(SearchEntry {
                lambda: lam,
                mu,
                weight: w,
                upper_flags: uf,
                lower_flags: lf,
                priority,
            });
        }
    }

    // Main loop.
    while let Some(entry) = heap.pop() {
        let lam = &entry.lambda;
        let mu = &entry.mu;
        let w = &entry.weight;
        let uf = &entry.upper_flags;
        let lf = &entry.lower_flags;

        let shape = fmt_shape(lam, mu);
        let flag_str = if uf.is_empty() && lf.is_empty() {
            String::new()
        } else {
            format!(" uf=({}) lf=({})", fmt_p(uf), fmt_p(lf))
        };
        eprint!(
            "\r\x1b[2K[{:>6} | Q:{:>7} | V:{:>7} | D:{:>5} | {:.0}s] {}  w=({}){}",
            computed_count,
            heap.len(),
            visited.len(),
            dedup_count,
            start.elapsed().as_secs_f64(),
            shape,
            fmt_p(w),
            flag_str,
        );
        let _ = std::io::stderr().flush();

        let lam_p = Partition::new(lam.clone());
        let mu_p = if mu.is_empty() {
            Partition::empty()
        } else {
            Partition::new(mu.clone())
        };
        let uf_opt = if uf.is_empty() {
            None
        } else {
            Some(uf.as_slice())
        };
        let lf_opt = if lf.is_empty() {
            None
        } else {
            Some(lf.as_slice())
        };

        let poly = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            compute_ehrhart(
                &lam_p,
                &mu_p,
                w,
                uf_opt,
                lf_opt,
                false,
                Some(cfg.max_states),
                true,
            )
        }));

        computed_count += 1;

        let poly = match poly {
            Ok(p) => p,
            Err(_) => {
                save_state(
                    lam,
                    mu,
                    w,
                    uf,
                    lf,
                    None,
                    None,
                    None,
                    false,
                    "panic",
                    db_pool.as_ref(),
                );
                continue;
            }
        };

        let poly_str = poly.display_factored();
        let is_new_poly = seen_polys.insert(poly_str.clone());
        if !is_new_poly {
            dedup_count += 1;
        }

        let raw_score = match score_c1_deviation(&poly) {
            Some(s) => s,
            None => {
                save_state(
                    lam,
                    mu,
                    w,
                    uf,
                    lf,
                    Some(poly.degree),
                    None,
                    Some(&poly_str),
                    false,
                    "computed",
                    db_pool.as_ref(),
                );
                continue;
            }
        };
        let priority = search_priority(&raw_score, &poly, lam, mu, w);

        // Check for negative coefficients.
        let has_negative = poly.has_negative_coefficient();
        save_state(
            lam,
            mu,
            w,
            uf,
            lf,
            Some(poly.degree),
            Some(&raw_score),
            Some(&poly_str),
            has_negative,
            "computed",
            db_pool.as_ref(),
        );
        if has_negative {
            found_negative = true;
            eprintln!();
            eprintln!("════════════════════════════════════════════════════════════════");
            eprintln!("  NEGATIVE COEFFICIENT FOUND!");
            eprintln!("  shape = {}", shape);
            eprintln!("  w = ({})", fmt_p(w));
            eprintln!("  degree = {}", poly.degree);
            eprintln!("  polynomial = {}", poly_str);
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
            save_record(
                lam,
                mu,
                w,
                uf,
                lf,
                poly.degree,
                &raw_score,
                &poly_str,
                true,
                "negative",
                db_pool.as_ref(),
            );
        }

        // Update best record.
        let is_new_best = match &best_score {
            None => true,
            Some(bs) => raw_score < *bs,
        };
        if is_new_best {
            best_score = Some(raw_score.clone());
            best_shape = shape.clone();
            best_degree = poly.degree;
            best_poly_str = poly_str.clone();

            // Compute h*-COM for display.
            let hstar = compute_hstar(&poly);
            let h_total: BigInt = hstar.iter().sum();
            let h_weighted: BigInt = hstar
                .iter()
                .enumerate()
                .map(|(i, h)| BigInt::from(i as u64) * h)
                .sum();
            let h_com_str = if !h_total.is_zero() {
                let com = BigRational::new(h_weighted, h_total);
                fmt_rational(&com)
            } else {
                "0".to_string()
            };

            eprintln!();
            eprintln!(
                "  ★ New record: c₁−H_d ≈ {}  (h*-COM ≈ {})  (deg {})  {}  w=({})",
                fmt_rational(&raw_score),
                h_com_str,
                poly.degree,
                best_shape,
                fmt_p(w),
            );
            eprintln!("    P(n) = {}", best_poly_str);

            save_record(
                lam,
                mu,
                w,
                uf,
                lf,
                poly.degree,
                &raw_score,
                &best_poly_str,
                has_negative,
                "new_best",
                db_pool.as_ref(),
            );
        }

        // Enqueue mutations (skip if this polynomial was already seen from another shape).
        if is_new_poly {
            let mutations = generate_mutations(lam, mu, w, uf, lf);
            for (ml, mm, mw, muf, mlf) in mutations {
                let next_priority = child_priority(&priority, &ml, &mm, &mw);
                let key = canonical_key(&ml, &mm, &mw, &muf, &mlf);
                if visited.insert(key) {
                    heap.push(SearchEntry {
                        lambda: ml,
                        mu: mm,
                        weight: mw,
                        upper_flags: muf,
                        lower_flags: mlf,
                        priority: next_priority,
                    });
                }
            }
        }
    }

    // Summary.
    eprintln!();
    eprintln!("═══════════════════════════ SEARCH COMPLETE ═══════════════════════════");
    eprintln!("  Computed: {}", computed_count);
    eprintln!(
        "  Deduped:  {} (same polynomial, mutations skipped)",
        dedup_count
    );
    eprintln!("  Distinct: {} polynomials", seen_polys.len());
    eprintln!("  Visited:  {}", visited.len());
    eprintln!("  Elapsed:  {:.1}s", start.elapsed().as_secs_f64());
    if found_negative {
        eprintln!("  Result:   NEGATIVE COEFFICIENTS FOUND (see above)");
    } else {
        eprintln!("  Result:   No negative coefficients found");
    }
    if let Some(bs) = &best_score {
        if best_shape.is_empty() {
            eprintln!(
                "  Best c₁−H_d ≈ {}  (loaded from records.jsonl; details omitted)",
                fmt_rational(bs)
            );
            eprintln!("═══════════════════════════════════════════════════════════════════════");
            return;
        }
        eprintln!(
            "  Best c₁−H_d ≈ {}  (deg {})  {}",
            fmt_rational(bs),
            best_degree,
            best_shape
        );
        eprintln!("    P(n) = {}", best_poly_str);
    }
    eprintln!("═══════════════════════════════════════════════════════════════════════");
}
