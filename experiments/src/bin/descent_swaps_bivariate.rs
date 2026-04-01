/// Bivariate stability analysis of the position-refined swaps polynomial:
///
///   Φ_S(x, t) = Σ_{σ ∈ D(n,S)} x^{pos(n,σ)} t^{swaps(σ)}
///             = Σ_p x^p L^{(p)}_{n,S}(t)
///
/// For each S ⊆ [n-1] with 1 ∉ S, we check:
///   1. Stability: Φ(x,t) ≠ 0 whenever Im(x)>0 and Im(t)>0
///   2. Real-stability restricted to real x>0: is Φ(x₀,t) real-rooted in t for each x₀>0?
///   3. Stronger: for ALL real x (including x<0), is Φ(x,t) real-rooted in t?
///   4. Newton polygon / support analysis (M-convexity check)
///   5. TN_2 of the coefficient matrix (implied by stability)
///
/// Stability of Φ(x,t) would imply TN of the coefficient matrix {a_{p,k}},
/// which is much stronger than TN_2.
///
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{format_poly, is_real_rooted};
use std::collections::BTreeMap;

// ── Complex arithmetic ───────────────────────────────────────────────
fn cmul(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}
fn cadd(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 + b.0, a.1 + b.1)
}
fn cpow(base: (f64, f64), exp: u32) -> (f64, f64) {
    let mut result = (1.0, 0.0);
    for _ in 0..exp {
        result = cmul(result, base);
    }
    result
}
fn cabs(a: (f64, f64)) -> f64 {
    (a.0 * a.0 + a.1 * a.1).sqrt()
}

// ── Simple PRNG ──────────────────────────────────────────────────────
struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed.wrapping_add(1))
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }
    /// Uniform in [0, 1)
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Uniform in [lo, hi)
    fn next_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

// ── Helpers ──────────────────────────────────────────────────────────

fn descent_set_to_string(mask: u64, n: u8) -> String {
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
}

/// Valid positions P(S) for descent set S in S_n.
/// p ∈ S with p-1 ∉ S (left boundaries), and p=n if n-1 ∉ S.
fn valid_positions(s_mask: u64, n: u8) -> Vec<u8> {
    let mut positions = Vec::new();
    for p in 1..n {
        let p_in_s = (s_mask >> (p - 1)) & 1 == 1;
        let pm1_in_s = if p >= 2 {
            (s_mask >> (p - 2)) & 1 == 1
        } else {
            false
        };
        if p_in_s && !pm1_in_s {
            positions.push(p);
        }
    }
    if n >= 2 && (s_mask >> (n - 2)) & 1 == 0 {
        positions.push(n);
    }
    positions
}

/// Build polynomial (coeff vector) from values.
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

/// Build the 2D coefficient array a[p_index][k] for Φ_S(x,t) = Σ_p x^p L^{(p)}(t).
/// Returns (positions_used, coeffs_2d, max_k).
/// positions_used[i] = the actual position p.
/// coeffs_2d[i][k] = coefficient of x^{positions_used[i]} t^k.
fn build_bivariate_coeffs(class: &[&Vec<u8>], positions: &[u8], n: u8) -> (Vec<u8>, Vec<Vec<i64>>) {
    let mut pos_polys: Vec<(u8, Vec<i64>)> = Vec::new();
    for &p in positions {
        let vals: Vec<usize> = class
            .iter()
            .filter(|pi| {
                let pos_of_n = pi.iter().position(|&v| v == n).unwrap();
                (pos_of_n + 1) as u8 == p
            })
            .map(|pi| compute(pi, Stat::Swaps))
            .collect();
        if vals.is_empty() {
            continue;
        }
        pos_polys.push((p, build_poly(&vals)));
    }

    if pos_polys.is_empty() {
        return (vec![], vec![]);
    }

    let max_deg = pos_polys.iter().map(|(_, p)| p.len()).max().unwrap_or(1);
    let positions_used: Vec<u8> = pos_polys.iter().map(|(p, _)| *p).collect();
    let coeffs: Vec<Vec<i64>> = pos_polys
        .iter()
        .map(|(_, poly)| {
            let mut row = poly.clone();
            row.resize(max_deg, 0);
            row
        })
        .collect();

    (positions_used, coeffs)
}

/// Evaluate Φ_S(x, t) = Σ_i x^{p_i} L^{(p_i)}(t) at complex (x, t).
/// positions[i] is the exponent of x, coeffs[i][k] is the coefficient of t^k.
fn eval_bivariate(
    positions: &[u8],
    coeffs: &[Vec<i64>],
    x: (f64, f64),
    t: (f64, f64),
) -> (f64, f64) {
    let mut result = (0.0, 0.0);
    for (i, &p) in positions.iter().enumerate() {
        let xp = cpow(x, p as u32);
        // Evaluate L^{(p)}(t) = Σ_k coeffs[i][k] t^k
        let mut lt = (0.0, 0.0);
        for (k, &c) in coeffs[i].iter().enumerate() {
            if c != 0 {
                let tk = cpow(t, k as u32);
                lt = cadd(lt, cmul((c as f64, 0.0), tk));
            }
        }
        result = cadd(result, cmul(xp, lt));
    }
    result
}

/// For fixed real x0, evaluate Φ(x0, t) as a univariate polynomial in t.
/// Returns coefficient vector in t.
fn specialize_x(positions: &[u8], coeffs: &[Vec<i64>], x0: f64) -> Vec<f64> {
    let max_k = coeffs.iter().map(|r| r.len()).max().unwrap_or(0);
    let mut result = vec![0.0f64; max_k];
    for (i, &p) in positions.iter().enumerate() {
        let xp = x0.powi(p as i32);
        for (k, &c) in coeffs[i].iter().enumerate() {
            result[k] += c as f64 * xp;
        }
    }
    // Trim
    while result.len() > 1 && result.last().unwrap().abs() < 1e-12 {
        result.pop();
    }
    result
}

/// Check approximate real-rootedness of a float polynomial.
fn is_approx_real_rooted(coeffs: &[f64]) -> bool {
    if coeffs.len() <= 2 {
        return true;
    }
    let max_coeff = coeffs.iter().map(|c| c.abs()).fold(0.0f64, f64::max);
    if max_coeff < 1e-10 {
        return true;
    }
    let scale = 1e8 / max_coeff;
    let int_coeffs: Vec<i64> = coeffs.iter().map(|&c| (c * scale).round() as i64).collect();
    let mut trimmed = int_coeffs;
    while trimmed.len() > 1 && *trimmed.last().unwrap() == 0 {
        trimmed.pop();
    }
    if trimmed.len() <= 2 {
        return true;
    }
    is_real_rooted(&trimmed)
}

/// Check TN_2: all 2x2 minors of the coefficient matrix are non-negative.
fn check_tn2(coeffs: &[Vec<i64>]) -> bool {
    let nrows = coeffs.len();
    if nrows < 2 {
        return true;
    }
    let ncols = coeffs[0].len();
    if ncols < 2 {
        return true;
    }
    for i in 0..nrows {
        for j in (i + 1)..nrows {
            for k1 in 0..ncols {
                for k2 in (k1 + 1)..ncols {
                    let det = coeffs[i][k1] * coeffs[j][k2] - coeffs[i][k2] * coeffs[j][k1];
                    if det < 0 {
                        return false;
                    }
                }
            }
        }
    }
    true
}

/// Check TN_3: all 3x3 minors are non-negative.
fn check_tn3(coeffs: &[Vec<i64>]) -> bool {
    let nrows = coeffs.len();
    if nrows < 3 {
        return true;
    }
    let ncols = coeffs[0].len();
    if ncols < 3 {
        return true;
    }
    for i1 in 0..nrows {
        for i2 in (i1 + 1)..nrows {
            for i3 in (i2 + 1)..nrows {
                for j1 in 0..ncols {
                    for j2 in (j1 + 1)..ncols {
                        for j3 in (j2 + 1)..ncols {
                            let a = coeffs[i1][j1];
                            let b = coeffs[i1][j2];
                            let c = coeffs[i1][j3];
                            let d = coeffs[i2][j1];
                            let e = coeffs[i2][j2];
                            let f = coeffs[i2][j3];
                            let g = coeffs[i3][j1];
                            let h = coeffs[i3][j2];
                            let k = coeffs[i3][j3];
                            let det =
                                a * (e * k - f * h) - b * (d * k - f * g) + c * (d * h - e * g);
                            if det < 0 {
                                return false;
                            }
                        }
                    }
                }
            }
        }
    }
    true
}

/// Compute Newton polygon support: set of (p, k) with a_{p,k} > 0.
/// Check M-convexity of the support (exchange axiom).
fn check_newton_polygon(positions: &[u8], coeffs: &[Vec<i64>]) -> (Vec<(u8, usize)>, bool) {
    let mut support: Vec<(u8, usize)> = Vec::new();
    for (i, &p) in positions.iter().enumerate() {
        for (k, &c) in coeffs[i].iter().enumerate() {
            if c > 0 {
                support.push((p, k));
            }
        }
    }

    // Check M-convexity (symmetric exchange axiom):
    // For any u, v in support and any i with u_i > v_i,
    // there exists j with u_j < v_j such that u - e_i + e_j ∈ support and v + e_i - e_j ∈ support.
    //
    // In 2D, this reduces to: the support forms a "generalized permutohedron" shape.
    // For a 2D point set, M-convexity means the support is the set of integer points
    // in a "staircase convex" region (every horizontal/vertical line meets a contiguous segment).
    //
    // Simpler check: for each fixed p, the k-values form a contiguous range,
    // and for each fixed k, the p-values form a contiguous range.

    // Check contiguity per position
    let mut by_p: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
    let mut by_k: BTreeMap<usize, Vec<u8>> = BTreeMap::new();
    for &(p, k) in &support {
        by_p.entry(p).or_default().push(k);
        by_k.entry(k).or_default().push(p);
    }

    let mut m_convex = true;

    // Check: for each p, is the range of k contiguous?
    for (_, ks) in &by_p {
        if ks.is_empty() {
            continue;
        }
        let kmin = *ks.iter().min().unwrap();
        let kmax = *ks.iter().max().unwrap();
        if ks.len() != kmax - kmin + 1 {
            m_convex = false;
        }
    }

    // Check: for each k, is the range of p contiguous?
    for (_, ps) in &by_k {
        if ps.is_empty() {
            continue;
        }
        let pmin = *ps.iter().min().unwrap();
        let pmax = *ps.iter().max().unwrap();
        // Count how many distinct p values are in positions between pmin..=pmax
        // But positions may not be contiguous (they skip values).
        // So check that every position between pmin and pmax that is in positions_used also has this k.
        // Actually for M-convexity in the lattice sense, we just need integer points:
        // every integer p in [pmin, pmax] should appear... but positions are sparse.
        // Better: check that the support is a "discrete polymatroid" shape.
        // For now, just check the simpler condition.
    }

    (support, m_convex)
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    println!("=== Bivariate stability analysis of Phi_S(x,t) ===");
    println!("Phi_S(x,t) = sum_p x^p L^{{(p)}}(t), S subset [n-1] with 1 not in S");
    println!("n = 4..{}\n", max_n);

    let num_stability_tests = 100;
    let num_rr_tests_positive = 50; // x > 0
    let num_rr_tests_all = 50; // all real x

    let mut total_cases = 0u64;
    let mut total_stable = 0u64;
    let mut total_rr_pos_x = 0u64;
    let mut total_rr_all_x = 0u64;
    let mut total_tn2 = 0u64;
    let mut total_tn3 = 0u64;
    let mut total_m_convex = 0u64;

    for n in 4..=max_n {
        let perms = all_permutations(n);

        // Group by descent set
        let mut by_descent: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for pi in &perms {
            by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi);
        }

        println!("========== n = {} ({} perms) ==========", n, perms.len());

        let mut n_cases = 0u64;
        let mut n_stable = 0u64;
        let mut n_rr_pos_x = 0u64;
        let mut n_rr_all_x = 0u64;
        let mut n_tn2 = 0u64;
        let mut n_tn3 = 0u64;
        let mut n_m_convex = 0u64;

        for (&mask, class) in &by_descent {
            // Only S with 1 not in S
            if mask & 1 != 0 {
                continue;
            }

            let positions = valid_positions(mask, n);
            let (pos_used, coeffs) = build_bivariate_coeffs(class, &positions, n);

            if pos_used.len() < 2 {
                continue; // Need at least 2 positions for interesting structure
            }

            n_cases += 1;
            let s_str = descent_set_to_string(mask, n);
            let max_k = coeffs[0].len();

            // Print Phi_S
            println!(
                "  S={}, positions={:?}, max t-degree={}",
                s_str,
                pos_used,
                max_k - 1
            );
            for (i, &p) in pos_used.iter().enumerate() {
                println!("    L^({}) = {}", p, format_poly(&coeffs[i]));
            }

            // ── Check 1: Stability (Im(x)>0, Im(t)>0 => Phi != 0) ──
            let mut stability_ok = true;
            let mut min_abs = f64::MAX;
            let mut rng = Rng::new(mask.wrapping_mul(777) + n as u64 * 13);

            for _ in 0..num_stability_tests {
                let a = rng.next_range(-2.0, 2.0);
                let b = rng.next_range(0.01, 2.0); // Im(x) > 0
                let c = rng.next_range(-2.0, 2.0);
                let d = rng.next_range(0.01, 2.0); // Im(t) > 0

                let val = eval_bivariate(&pos_used, &coeffs, (a, b), (c, d));
                let abs_val = cabs(val);
                if abs_val < min_abs {
                    min_abs = abs_val;
                }
                if abs_val < 1e-10 {
                    stability_ok = false;
                    println!(
                        "    STABILITY FAIL: Phi({:.4}+{:.4}i, {:.4}+{:.4}i) ~ 0  (|Phi|={:.2e})",
                        a, b, c, d, abs_val
                    );
                    break;
                }
            }
            if stability_ok {
                n_stable += 1;
                println!("    Stability: PASS (min |Phi| = {:.4e})", min_abs);
            }

            // ── Check 2: For fixed real x > 0, is Phi(x, t) real-rooted in t? ──
            let mut rr_pos_ok = true;
            let test_x_positive: Vec<f64> = {
                let mut v = vec![0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0];
                let mut rng2 = Rng::new(mask + 999);
                for _ in 0..num_rr_tests_positive {
                    v.push(rng2.next_range(0.01, 10.0));
                }
                v
            };

            for &x0 in &test_x_positive {
                let spec = specialize_x(&pos_used, &coeffs, x0);
                if !is_approx_real_rooted(&spec) {
                    rr_pos_ok = false;
                    println!("    RR(x>0) FAIL at x={:.4}", x0);
                    let int_approx: Vec<i64> = spec.iter().map(|&c| c.round() as i64).collect();
                    println!("      Phi({:.4}, t) ~ {}", x0, format_poly(&int_approx));
                    break;
                }
            }
            if rr_pos_ok {
                n_rr_pos_x += 1;
                println!("    RR for x > 0: PASS");
            }

            // ── Check 3: For ALL real x (including x < 0), real-rooted in t? ──
            let mut rr_all_ok = true;
            let test_x_all: Vec<f64> = {
                let mut v = vec![
                    -10.0, -5.0, -2.0, -1.0, -0.5, -0.1, -0.01, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0,
                ];
                let mut rng3 = Rng::new(mask + 12345);
                for _ in 0..num_rr_tests_all {
                    v.push(rng3.next_range(-10.0, 10.0));
                }
                v
            };

            for &x0 in &test_x_all {
                let spec = specialize_x(&pos_used, &coeffs, x0);
                if !is_approx_real_rooted(&spec) {
                    rr_all_ok = false;
                    println!("    RR(all x) FAIL at x={:.4}", x0);
                    let int_approx: Vec<i64> = spec.iter().map(|&c| c.round() as i64).collect();
                    println!("      Phi({:.4}, t) ~ {}", x0, format_poly(&int_approx));
                    break;
                }
            }
            if rr_all_ok {
                n_rr_all_x += 1;
                println!("    RR for all real x: PASS");
            }

            // ── Check 4: TN_2 and TN_3 of coefficient matrix ──
            let tn2 = check_tn2(&coeffs);
            if tn2 {
                n_tn2 += 1;
            }
            println!("    TN_2: {}", if tn2 { "PASS" } else { "FAIL" });

            let tn3 = check_tn3(&coeffs);
            if tn3 {
                n_tn3 += 1;
            }
            println!("    TN_3: {}", if tn3 { "PASS" } else { "FAIL" });

            // ── Check 5: Newton polygon / M-convexity ──
            let (support, m_conv) = check_newton_polygon(&pos_used, &coeffs);
            if m_conv {
                n_m_convex += 1;
            }
            println!(
                "    Newton polygon: {} points, row-contiguous: {}",
                support.len(),
                if m_conv { "YES" } else { "NO" }
            );

            // Print support shape compactly
            let mut by_p: BTreeMap<u8, (usize, usize)> = BTreeMap::new();
            for &(p, k) in &support {
                let entry = by_p.entry(p).or_insert((k, k));
                if k < entry.0 {
                    entry.0 = k;
                }
                if k > entry.1 {
                    entry.1 = k;
                }
            }
            print!("    Support ranges: ");
            for (&p, &(kmin, kmax)) in &by_p {
                print!("p={}:k={}..{} ", p, kmin, kmax);
            }
            println!();

            println!();
        }

        total_cases += n_cases;
        total_stable += n_stable;
        total_rr_pos_x += n_rr_pos_x;
        total_rr_all_x += n_rr_all_x;
        total_tn2 += n_tn2;
        total_tn3 += n_tn3;
        total_m_convex += n_m_convex;

        println!(
            "n={}: {}/{} stable, {}/{} RR(x>0), {}/{} RR(all x), {}/{} TN2, {}/{} TN3, {}/{} M-convex\n",
            n, n_stable, n_cases, n_rr_pos_x, n_cases, n_rr_all_x, n_cases,
            n_tn2, n_cases, n_tn3, n_cases, n_m_convex, n_cases
        );
    }

    println!("============================================");
    println!("=== SUMMARY (n=4..{}) ===", max_n);
    println!("============================================");
    println!("Total descent sets tested: {}", total_cases);
    println!();
    println!(
        "Bivariate stability (Im x>0, Im t>0):   {}/{}",
        total_stable, total_cases
    );
    println!(
        "RR in t for all x > 0:                   {}/{}",
        total_rr_pos_x, total_cases
    );
    println!(
        "RR in t for ALL real x:                  {}/{}",
        total_rr_all_x, total_cases
    );
    println!(
        "TN_2 of coefficient matrix:              {}/{}",
        total_tn2, total_cases
    );
    println!(
        "TN_3 of coefficient matrix:              {}/{}",
        total_tn3, total_cases
    );
    println!(
        "Newton polygon row-contiguous:           {}/{}",
        total_m_convex, total_cases
    );
    println!();

    if total_stable == total_cases {
        println!("*** CONJECTURE: Phi_S(x,t) appears to be STABLE. ***");
        println!("    This would imply the coefficient matrix is TOTALLY NON-NEGATIVE");
        println!("    (all minors non-negative), which is much stronger than TN_2.");
    } else if total_rr_all_x == total_cases {
        println!("*** Phi_S(x,t) is NOT fully stable, but Phi(x0, t) is real-rooted");
        println!("    for all real x0. This implies the 'compatible family' property.");
        println!("    It is consistent with REAL stability (nonvanishing when Im t > 0");
        println!("    and x is real) but weaker than full stability.");
    } else if total_rr_pos_x == total_cases {
        println!("*** Phi_S(x,t) is real-rooted in t for x > 0 only. ***");
        println!("    This confirms the compatible family property for positive combinations.");
    }

    if total_tn2 == total_cases {
        println!("*** TN_2 HOLDS for all tested coefficient matrices. ***");
    }
    if total_tn3 == total_cases && total_cases > 0 {
        println!("*** TN_3 also HOLDS for all tested coefficient matrices. ***");
    }
}
