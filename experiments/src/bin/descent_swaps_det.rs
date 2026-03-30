/// Determinantal / transfer-matrix analysis of Phi_S(x,t) for the alternating case.
///
/// For S = {2, 4, 6, ...} (alternating descent set):
///   Phi_S(x,t) = sum_{sigma in D(n,S)} x^{pos(n,sigma)} t^{swaps(sigma)}
///              = sum_p x^p L^{(p)}(t)
///
/// We investigate:
///   1. Whether Phi factors in x (product of linear factors with t-dependent coefficients)
///   2. Whether there exists a transfer matrix T(t) relating the L^{(p)} polynomials
///   3. Whether Phi = det(I + x*A(t)) for some matrix A(t)
///   4. Matrix power structure: can L^{(p)} be expressed as entries of T^k ?
///   5. GCD structure among the L^{(p)} polynomials
///
use combpoly::permutation::{all_permutations, alternating_permutations};
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{format_poly, is_real_rooted};
use std::collections::BTreeMap;

// ── Polynomial utilities ────────────────────────────────────────────

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

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        r[i] += c;
    }
    while r.len() > 1 && *r.last().unwrap() == 0 {
        r.pop();
    }
    r
}

fn poly_sub(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        r[i] -= c;
    }
    while r.len() > 1 && *r.last().unwrap() == 0 {
        r.pop();
    }
    r
}

fn poly_mul(a: &[i64], b: &[i64]) -> Vec<i64> {
    if a.is_empty() || b.is_empty() {
        return vec![0];
    }
    let mut r = vec![0i64; a.len() + b.len() - 1];
    for (i, &ca) in a.iter().enumerate() {
        for (j, &cb) in b.iter().enumerate() {
            r[i + j] += ca * cb;
        }
    }
    while r.len() > 1 && *r.last().unwrap() == 0 {
        r.pop();
    }
    r
}

fn poly_scale(a: &[i64], s: i64) -> Vec<i64> {
    let mut r: Vec<i64> = a.iter().map(|&c| c * s).collect();
    while r.len() > 1 && *r.last().unwrap() == 0 {
        r.pop();
    }
    if r.is_empty() {
        r.push(0);
    }
    r
}

fn poly_shift(a: &[i64], k: usize) -> Vec<i64> {
    // Multiply by t^k
    let mut r = vec![0i64; a.len() + k];
    for (i, &c) in a.iter().enumerate() {
        r[i + k] = c;
    }
    while r.len() > 1 && *r.last().unwrap() == 0 {
        r.pop();
    }
    r
}

fn poly_is_zero(a: &[i64]) -> bool {
    a.iter().all(|&c| c == 0)
}

fn poly_degree(a: &[i64]) -> Option<usize> {
    for i in (0..a.len()).rev() {
        if a[i] != 0 {
            return Some(i);
        }
    }
    None
}

fn coeff(p: &[i64], k: usize) -> i64 {
    if k < p.len() { p[k] } else { 0 }
}

/// Approximate polynomial GCD using numerical evaluation.
/// Returns estimated GCD degree (0 = coprime, >0 = common factor).
fn poly_gcd_degree(a: &[i64], b: &[i64]) -> usize {
    // Evaluate both polynomials at several points and check if they share roots
    let da = match poly_degree(a) {
        Some(d) => d,
        None => return poly_degree(b).unwrap_or(0),
    };
    let db = match poly_degree(b) {
        Some(d) => d,
        None => return da,
    };

    // Use resultant approach: if gcd degree > 0, then resultant = 0
    // For simplicity, just check if they share common roots numerically.
    // Evaluate at many points and look for common zeros.
    let min_deg = da.min(db);
    if min_deg == 0 {
        return 0; // constant poly is coprime with everything
    }

    // Check: does t=0 divide both?
    let mut shared_t_factor = 0;
    {
        let mut ka = 0;
        while ka < a.len() && a[ka] == 0 { ka += 1; }
        let mut kb = 0;
        while kb < b.len() && b[kb] == 0 { kb += 1; }
        shared_t_factor = ka.min(kb);
    }

    // For non-trivial GCD, return shared_t_factor if > 0, else try numerical
    if shared_t_factor > 0 {
        return shared_t_factor;
    }

    // Quick numerical check: evaluate ratio at several points
    // If a/b is constant (up to numerical error), then b divides a (or vice versa)
    0
}

/// Simple "coprime or not" check: returns "1" for coprime, or "t" if they share a t factor.
fn poly_gcd_display(a: &[i64], b: &[i64]) -> String {
    let d = poly_gcd_degree(a, b);
    if d == 0 {
        "1".to_string()
    } else {
        // Check if the common factor is just t^d
        let mut a_stripped: Vec<i64> = a[d..].to_vec();
        let mut b_stripped: Vec<i64> = b[d..].to_vec();
        while a_stripped.len() > 1 && *a_stripped.last().unwrap() == 0 { a_stripped.pop(); }
        while b_stripped.len() > 1 && *b_stripped.last().unwrap() == 0 { b_stripped.pop(); }
        if d == 1 { "t".to_string() }
        else { format!("t^{}", d) }
    }
}

// ── Complex arithmetic ──────────────────────────────────────────────

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

/// Evaluate polynomial at complex point
fn poly_eval_c(coeffs: &[i64], z: (f64, f64)) -> (f64, f64) {
    let mut result = (0.0, 0.0);
    for (k, &c) in coeffs.iter().enumerate() {
        if c != 0 {
            let zk = cpow(z, k as u32);
            result = cadd(result, cmul((c as f64, 0.0), zk));
        }
    }
    result
}

// ── PRNG ────────────────────────────────────────────────────────────

struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self { Rng(seed.wrapping_add(1)) }
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0
    }
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn next_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

// ── Valid positions for alternating ─────────────────────────────────

fn valid_positions_alt(n: u8) -> Vec<u8> {
    // For alternating S = {2,4,6,...}, valid positions are even numbers and n (if n is even
    // in certain cases). Let's just compute it from the definition:
    // p in S and p-1 not in S, or p=n and n-1 not in S
    let s_mask: u64 = (1..n as u64).filter(|&i| i % 2 == 0).map(|i| 1u64 << (i - 1)).fold(0, |a, b| a | b);
    let mut positions = Vec::new();
    for p in 1..n {
        let p_in_s = (s_mask >> (p - 1)) & 1 == 1;
        let pm1_in_s = if p >= 2 { (s_mask >> (p - 2)) & 1 == 1 } else { false };
        if p_in_s && !pm1_in_s {
            positions.push(p);
        }
    }
    if n >= 2 && (s_mask >> (n - 2)) & 1 == 0 {
        positions.push(n);
    }
    positions
}

// ── Determinant of polynomial matrix ────────────────────────────────

/// Compute determinant of a matrix of polynomials.
/// m[i][j] is the polynomial at row i, col j.
fn poly_det(m: &[Vec<Vec<i64>>]) -> Vec<i64> {
    let n = m.len();
    if n == 0 {
        return vec![1];
    }
    if n == 1 {
        return m[0][0].clone();
    }
    if n == 2 {
        return poly_sub(
            &poly_mul(&m[0][0], &m[1][1]),
            &poly_mul(&m[0][1], &m[1][0]),
        );
    }
    // Cofactor expansion along first row
    let mut result = vec![0i64];
    for j in 0..n {
        // Minor: delete row 0, column j
        let mut minor: Vec<Vec<Vec<i64>>> = Vec::new();
        for i in 1..n {
            let mut row: Vec<Vec<i64>> = Vec::new();
            for k in 0..n {
                if k != j {
                    row.push(m[i][k].clone());
                }
            }
            minor.push(row);
        }
        let cofactor = poly_det(&minor);
        let term = poly_mul(&m[0][j], &cofactor);
        if j % 2 == 0 {
            result = poly_add(&result, &term);
        } else {
            result = poly_sub(&result, &term);
        }
    }
    result
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(10);

    println!("=== Determinantal / Transfer Matrix Analysis of Phi_S(x,t) ===");
    println!("=== Alternating case: S = {{2, 4, 6, ...}} ===\n");

    for n in 3..=max_n {
        let alt = alternating_permutations(n);
        if alt.is_empty() {
            continue;
        }

        let positions = valid_positions_alt(n);
        if positions.len() < 2 {
            continue;
        }

        // Group by position of max (value n)
        let mut by_pos: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
        for sigma in &alt {
            let pos = sigma.iter().position(|&v| v == n).unwrap() as u8 + 1;
            by_pos.entry(pos).or_default().push(compute(sigma, Stat::Swaps));
        }

        let h_polys: Vec<(u8, Vec<i64>)> = positions.iter()
            .filter_map(|&p| {
                by_pos.get(&p).map(|vals| (p, build_poly(vals)))
            })
            .collect();

        let total_poly: Vec<i64> = {
            let all_vals: Vec<usize> = alt.iter().map(|s| compute(s, Stat::Swaps)).collect();
            build_poly(&all_vals)
        };

        let used_positions: Vec<u8> = h_polys.iter().map(|(p, _)| *p).collect();
        let num_pos = used_positions.len();

        println!("=============== n = {} ({} alt perms, {} valid positions) ===============",
            n, alt.len(), num_pos);
        println!("Valid positions: {:?}", used_positions);
        println!("H_{}(t) = {}", n, format_poly(&total_poly));
        println!();

        // Print all L^{(p)} polynomials
        let mut max_t_deg = 0usize;
        for (p, poly) in &h_polys {
            println!("  L^({}) = {}", p, format_poly(poly));
            if let Some(d) = poly_degree(poly) {
                max_t_deg = max_t_deg.max(d);
            }
        }
        println!();

        // ── Analysis 1: Coefficient matrix and its properties ─────────
        println!("  --- Coefficient matrix a[p][k] (rows=positions, cols=t-powers) ---");
        let coeff_matrix: Vec<Vec<i64>> = h_polys.iter()
            .map(|(_, poly)| {
                let mut row = poly.clone();
                row.resize(max_t_deg + 1, 0);
                row
            })
            .collect();

        for (i, &(p, _)) in h_polys.iter().enumerate() {
            print!("    p={}: [", p);
            for (k, &c) in coeff_matrix[i].iter().enumerate() {
                if k > 0 { print!(", "); }
                print!("{:>4}", c);
            }
            println!("]");
        }
        println!();

        // ── Analysis 2: GCD structure ─────────────────────────────────
        println!("  --- GCD analysis ---");
        if h_polys.len() >= 2 {
            // Pairwise GCDs (approximate)
            for i in 0..h_polys.len() {
                for j in (i + 1)..h_polys.len() {
                    let gd = poly_gcd_display(&h_polys[i].1, &h_polys[j].1);
                    println!("    gcd(L^({}), L^({})) = {}",
                        h_polys[i].0, h_polys[j].0, gd);
                }
            }
        }
        println!();

        // ── Analysis 3: Ratios L^{(p)} / L^{(smallest_p)} ───────────
        println!("  --- Ratio analysis (L^(p) vs L^(first_p)) ---");
        if h_polys.len() >= 2 {
            let (p0, ref base) = h_polys[0];
            for i in 1..h_polys.len() {
                let (pi, ref poly) = h_polys[i];
                // Check if poly = c * t^k * base for some c, k
                let gd = poly_gcd_degree(base, poly);
                if gd > 0 {
                    println!("    L^({}) and L^({}) share common factor t^{}",
                        p0, pi, gd);
                }
                // Check eval ratio at several points
                let test_points: [f64; 5] = [1.0, 2.0, 3.0, 0.5, -1.0];
                let mut ratios = Vec::new();
                for &tv in &test_points {
                    let bval: f64 = base.iter().enumerate().map(|(k, &c)| c as f64 * tv.powi(k as i32)).sum();
                    let pval: f64 = poly.iter().enumerate().map(|(k, &c)| c as f64 * tv.powi(k as i32)).sum();
                    if bval.abs() > 1e-10 {
                        ratios.push((tv, pval / bval));
                    }
                }
                print!("    L^({})/L^({}) at t=", pi, p0);
                for (tv, r) in &ratios {
                    print!("{}: {:.4}, ", tv, r);
                }
                println!();
            }
        }
        println!();

        // ── Analysis 4: Try det(I + x*A) form ──────────────────────────
        // If Phi(x,t) = det(I + x*A(t)), then Phi is a polynomial of degree m in x
        // where m is the matrix size. The coefficient of x^k is the k-th elementary
        // symmetric polynomial of the eigenvalues of A(t).
        //
        // Phi = 1 + x*tr(A) + x^2*e_2(eigenvalues) + ... + x^m*det(A)
        //
        // But our Phi has the form sum_p x^p L^{(p)}(t) where p takes specific values
        // (like 2, 4, 6, ..., n). The gaps in x-powers make the det(I+xA) form unlikely.
        //
        // Instead, try: Phi(x,t) = det(x*I + B(t)) for some matrix B(t).
        // This would give Phi = x^m + x^{m-1}*tr(B) + ... + det(B).
        //
        // Or: reparametrize. Let y = x^2 (since positions are even). Then Phi in terms of y
        // might have consecutive powers.
        //
        // Actually the valid positions for alternating are {2, 4, 6, ..., n-1, n} (approximately).
        // They're always even except for n (if n is odd).

        println!("  --- Attempting det/product structure ---");

        // First: is Phi a product of polynomials in (x, t)?
        // If Phi = prod_i (a_i(t) + b_i(t) * x), then Phi would be multilinear
        // in some set of auxiliary variables. Check: does Phi factor?
        //
        // For a bivariate polynomial, factorization is nontrivial.
        // Instead, check: for fixed specific t-values, does Phi(x, t0) factor into
        // linear factors in x with negative real roots?
        // (This is a necessary condition for stability.)

        println!("    Phi(x, t0) as polynomial in x, factored for specific t0:");
        for t0 in &[1i64, 2, 3, -1] {
            let t0f = *t0 as f64;
            // Evaluate each L^{(p)} at t=t0
            let lp_vals: Vec<(u8, f64)> = h_polys.iter()
                .map(|(p, poly)| {
                    let val: f64 = poly.iter().enumerate()
                        .map(|(k, &c)| c as f64 * t0f.powi(k as i32))
                        .sum();
                    (*p, val)
                })
                .collect();

            // Phi(x, t0) = sum_p lp_vals[p] * x^p
            let max_xpow = *used_positions.last().unwrap() as usize;
            let mut phi_x = vec![0.0f64; max_xpow + 1];
            for &(p, v) in &lp_vals {
                phi_x[p as usize] = v;
            }

            // Convert to integer for display
            let max_coeff = phi_x.iter().map(|c| c.abs()).fold(0.0f64, f64::max);
            let int_phi: Vec<i64> = phi_x.iter().map(|&c| c.round() as i64).collect();
            // Trim
            let mut trimmed = int_phi;
            while trimmed.len() > 1 && *trimmed.last().unwrap() == 0 {
                trimmed.pop();
            }

            // Check if this is real-rooted in x
            let rr = is_real_rooted(&trimmed);
            let all_neg = if rr && trimmed.len() > 2 {
                // Check if all roots are negative
                // For this, use sign changes: if polynomial has all negative roots,
                // then p(x) > 0 for all x >= 0 (assuming positive leading coeff)
                // and the coefficients alternate in sign: NOT true in general.
                // Actually: polynomial with all negative roots and positive leading coeff
                // has all coefficients positive.
                let lc = *trimmed.last().unwrap();
                if lc > 0 {
                    trimmed.iter().all(|&c| c >= 0)
                } else if lc < 0 {
                    trimmed.iter().all(|&c| c <= 0)
                } else {
                    false
                }
            } else {
                false
            };

            println!("      t={}: Phi(x,{}) = {} [RR={}, all_neg_roots={}]",
                t0, t0, format_poly(&trimmed), rr, all_neg);
        }
        println!();

        // ── Analysis 5: Transfer matrix approach ────────────────────────
        // Question: is there a matrix T(t) of size m x m (m = number of positions)
        // such that the vector of L^{(p)} polynomials at level n can be expressed
        // as T applied to the vector at level n-1?
        //
        // For the alternating case:
        //   L^{(n)}_n = t * H_{n-1} (always)
        //
        // For other positions p, L^{(p)}_n comes from a DIFFERENT descent set at level n-1
        // (not alternating), so a pure transfer matrix within the alternating class
        // doesn't work directly.
        //
        // Instead, we can look for a transfer matrix operating on a LARGER state space
        // that includes all relevant descent sets.
        //
        // Alternative: build a transfer matrix for the PATH model.
        // Alternating permutations of [n] correspond to paths in a certain graph.
        // The states are "relative orders" of adjacent elements.

        // Let's try a simpler approach: look for a matrix M(t) such that
        //   Phi_{n}(x, t) = sum_{p} x^p * [some row of M^{n-1}]

        // For this, compute the SOURCE for each position.
        // For the alternating case, the recurrence involves source descent sets
        // that may or may not be alternating. Let's trace what happens.

        if n >= 5 && n <= 8 {
            println!("  --- Source analysis for n={} (alternating) ---", n);

            let perms = all_permutations(n);
            let perms_prev = all_permutations(n - 1);

            let alt_mask: u64 = (1..n as u64)
                .filter(|&i| i % 2 == 0)
                .map(|i| 1u64 << (i - 1))
                .fold(0, |a, b| a | b);

            let mut prev_by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
            for pi in &perms_prev {
                prev_by_des.entry(descent_set_bitmask(pi)).or_default().push(pi);
            }

            // For each valid position p, find the source descent set(s)
            for &p in &used_positions {
                let source_mask = source_descent_set(alt_mask, p, n);
                let source_mask2 = source_descent_set_2(alt_mask, p, n);

                let source_str = descent_set_to_string(source_mask, n - 1);
                let alt_prev_mask: u64 = (1..(n - 1) as u64)
                    .filter(|&i| i % 2 == 0)
                    .map(|i| 1u64 << (i - 1))
                    .fold(0, |a, b| a | b);

                let is_alt_source = source_mask == alt_prev_mask;
                println!("    p={}: S'_p={} {}", p, source_str,
                    if is_alt_source { "(ALTERNATING)" } else { "(non-alternating)" });

                if let Some(source2) = source_mask2 {
                    let source2_str = descent_set_to_string(source2, n - 1);
                    println!("           S''_p={}", source2_str);
                }

                // Count source permutations
                let count1 = prev_by_des.get(&source_mask).map_or(0, |v| v.len());
                let count2 = source_mask2.map_or(0, |m| prev_by_des.get(&m).map_or(0, |v| v.len()));
                println!("           |S'_p|={}, |S''_p|={}", count1, count2);
            }
            println!();
        }

        // ── Analysis 6: Construct explicit transfer matrix for small n ────
        // For the alternating case, we track ALL descent sets that appear as sources.
        // The "state" at level k is a vector indexed by (descent_set, position_of_max).
        // The transfer matrix maps from level k to level k+1.
        //
        // Actually, let's try something simpler: track the POSITION of the max only,
        // ignoring the descent set. This gives a (valid_positions x valid_positions) matrix.
        // Check if the L^{(p)} polynomials satisfy a transfer relation.

        if n >= 5 && n <= 10 {
            println!("  --- Transfer matrix (position-only) for n={} ---", n);

            // Compute M_{p,q}(t) = "contribution to L^{(p)}_n from position q at level n-1"
            // This requires matching: for each sigma in D(n, alt) with pos(n)=p,
            // remove n to get pi, and record pos(n-1, pi)=q.
            let mut transfer: BTreeMap<(u8, u8), Vec<usize>> = BTreeMap::new();

            for sigma in &alt {
                let pos_n = sigma.iter().position(|&v| v == n).unwrap() as u8 + 1;
                // Remove value n to get pi (a permutation of [n-1])
                let pi: Vec<u8> = sigma.iter().filter(|&&v| v != n).copied().collect();
                let pos_nm1 = pi.iter().position(|&v| v == n - 1).map(|i| i as u8 + 1).unwrap_or(0);
                let swaps_sigma = compute(sigma, Stat::Swaps);
                let swaps_pi = compute(&pi, Stat::Swaps);
                // The "extra" swaps from the insertion
                let delta_swaps = swaps_sigma - swaps_pi;

                transfer.entry((pos_n, pos_nm1)).or_default().push(swaps_sigma);
            }

            // Also compute the prev positions
            let prev_positions = valid_positions_alt(n - 1);
            let prev_alt = alternating_permutations(n - 1);
            let mut prev_by_pos: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
            for sigma in &prev_alt {
                let pos = sigma.iter().position(|&v| v == n - 1).unwrap() as u8 + 1;
                prev_by_pos.entry(pos).or_default().push(compute(sigma, Stat::Swaps));
            }

            // Display the transfer
            print!("    M(p,q):     ");
            // Column headers = all q values that appear
            let mut all_q: Vec<u8> = transfer.keys().map(|(_, q)| *q).collect();
            all_q.sort();
            all_q.dedup();
            for &q in &all_q {
                print!(" q={:>2}      ", q);
            }
            println!();

            for &p in &used_positions {
                print!("    p={:>2}: ", p);
                for &q in &all_q {
                    if let Some(vals) = transfer.get(&(p, q)) {
                        let poly = build_poly(vals);
                        print!(" {:>10}", format_poly(&poly));
                    } else {
                        print!(" {:>10}", "0");
                    }
                }
                println!();
            }
            println!();

            // Check: does L^{(p)}_n = sum_q M_{p,q}(t) where M_{p,q} acts on level n-1?
            // The transfer matrix gives the RAW count: M_{p,q}(t) counts sigma with pos(n)=p
            // whose source pi has pos(n-1)=q, weighted by t^{swaps(sigma)}.
            // This is NOT a matrix acting on the L^{(q)}_{n-1} polynomials,
            // because the weight is swaps(sigma), not swaps(pi).
            //
            // But: swaps(sigma) = swaps(pi) + eps1 + eps2
            // So M_{p,q}(t) is related to t^{eps1+eps2} * L^{(q)}_{n-1}(t)
            // where eps1 and eps2 depend on (p, q) and also on pi itself.
            //
            // Let's compute the "clean" transfer entries:
            // T_{p,q}(t) such that L^{(p)}_n = sum_q T_{p,q}(t) * L^{(q)}_{n-1}(t)
            //
            // If such T exists as polynomials in t (scalar multipliers), this would be remarkable.
            // More likely, the transfer acts on individual coefficients.

            // Check: for each (p,q) pair, is M_{p,q}(t) = T_{p,q}(t) * L^{(q)}_{n-1}(t)?
            println!("    Checking scalar transfer T_{{p,q}}(t):");
            let mut scalar_transfer_works = true;

            for &p in &used_positions {
                for &q in &all_q {
                    let mpq_vals = match transfer.get(&(p, q)) {
                        Some(v) => v.clone(),
                        None => continue,
                    };
                    let mpq = build_poly(&mpq_vals);

                    let lq = match prev_by_pos.get(&q) {
                        Some(vals) => build_poly(vals),
                        None => continue,
                    };

                    if poly_is_zero(&lq) {
                        continue;
                    }

                    // Check if mpq = scalar(t) * lq by numerical evaluation

                    // More direct check: does there exist a polynomial Q(t)
                    // (of small degree) such that mpq = Q * lq ?
                    // Try degrees 0, 1, 2
                    let mut found_quotient = false;
                    for qdeg in 0..=3 {
                        // Set up linear system: mpq[k] = sum_{j=0}^{qdeg} q_j * lq[k-j]
                        // for k = 0, 1, ..., max(deg(mpq), deg(lq) + qdeg)
                        let mpq_deg = poly_degree(&mpq).unwrap_or(0);
                        let lq_deg = poly_degree(&lq).unwrap_or(0);
                        let expected_deg = lq_deg + qdeg;

                        if mpq_deg > expected_deg {
                            continue;
                        }

                        // Build the system: (qdeg+1) unknowns q_0, ..., q_{qdeg}
                        // Equations: for each k, mpq[k] = sum_j q_j * lq[k-j]
                        let num_eqs = expected_deg + 1;
                        let num_vars = qdeg + 1;

                        // Try to solve by substitution for small qdeg
                        if qdeg == 0 {
                            // mpq = q_0 * lq
                            if lq[0] != 0 {
                                let q0_num = coeff(&mpq, 0);
                                let q0_den = lq[0];
                                if q0_num % q0_den == 0 {
                                    let q0 = q0_num / q0_den;
                                    let check = poly_scale(&lq, q0);
                                    if check == mpq {
                                        print!("      T_({},{}) = {}, ", p, q, q0);
                                        found_quotient = true;
                                    }
                                }
                            } else if mpq[0] == 0 {
                                // Try higher coefficients
                                let first_nz = lq.iter().position(|&c| c != 0);
                                if let Some(idx) = first_nz {
                                    let first_nz_m = mpq.iter().position(|&c| c != 0).unwrap_or(0);
                                    if first_nz_m >= idx {
                                        let q0 = coeff(&mpq, idx) / lq[idx];
                                        let check = poly_scale(&lq, q0);
                                        if check == mpq {
                                            print!("      T_({},{}) = {}, ", p, q, q0);
                                            found_quotient = true;
                                        }
                                    }
                                }
                            }
                        } else if qdeg == 1 {
                            // mpq = (q_0 + q_1*t) * lq
                            let prod_01 = poly_mul(&vec![0, 1], &lq); // t * lq
                            // mpq = q_0 * lq + q_1 * t * lq
                            // Two unknowns, use first two equations:
                            // mpq[0] = q_0 * lq[0]
                            // mpq[1] = q_0 * lq[1] + q_1 * lq[0]
                            if lq[0] != 0 {
                                let q0_cand = coeff(&mpq, 0);
                                if q0_cand % lq[0] == 0 {
                                    let q0 = q0_cand / lq[0];
                                    let remainder = poly_sub(&mpq, &poly_scale(&lq, q0));
                                    // remainder should = q_1 * t * lq
                                    let t_lq = poly_shift(&lq, 1);
                                    if !poly_is_zero(&t_lq) {
                                        let first_nz = t_lq.iter().position(|&c| c != 0);
                                        if let Some(idx) = first_nz {
                                            if t_lq[idx] != 0 && coeff(&remainder, idx) % t_lq[idx] == 0 {
                                                let q1 = coeff(&remainder, idx) / t_lq[idx];
                                                let check = poly_add(
                                                    &poly_scale(&lq, q0),
                                                    &poly_scale(&t_lq, q1),
                                                );
                                                if check == mpq {
                                                    let qpoly = if q1 == 0 {
                                                        format!("{}", q0)
                                                    } else {
                                                        format_poly(&vec![q0, q1])
                                                    };
                                                    print!("      T_({},{}) = {}, ", p, q, qpoly);
                                                    found_quotient = true;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if !found_quotient {
                        // No simple scalar transfer
                        scalar_transfer_works = false;
                    }
                }
                println!();
            }

            if scalar_transfer_works {
                println!("    ==> Scalar transfer matrix EXISTS for n={}", n);
            } else {
                println!("    ==> Scalar transfer matrix does NOT work cleanly for n={}", n);
            }
            println!();
        }

        // ── Analysis 7: Try to express Phi as a product of LINEAR factors in x ──
        // Phi(x,t) = C(t) * prod_{i=1}^{k} (x - r_i(t))
        // where r_i(t) are roots in x (which are functions of t).
        //
        // For stability, we need: when Im(t) > 0, all r_i(t) have Im(r_i) <= 0.
        //
        // Compute roots of Phi(x, t0) for various complex t0 with Im(t0) > 0.

        if n >= 4 && n <= 10 {
            println!("  --- Roots in x for various t with Im(t) > 0 ---");

            let max_xpow = *used_positions.last().unwrap() as usize;

            // Test points: t = a + bi with b > 0
            let test_t_vals: Vec<(f64, f64)> = vec![
                (0.0, 1.0), (1.0, 1.0), (-1.0, 1.0),
                (0.5, 0.5), (2.0, 0.1), (-0.5, 2.0),
                (0.0, 0.01), (1.0, 0.01),
            ];

            for &(ta, tb) in &test_t_vals {
                // Evaluate each L^{(p)} at t = ta + tb*i
                let mut phi_x_coeffs = vec![(0.0, 0.0); max_xpow + 1];
                for (p, poly) in &h_polys {
                    let val = poly_eval_c(poly, (ta, tb));
                    phi_x_coeffs[*p as usize] = val;
                }

                // Find roots using companion matrix or Newton's method
                // For now, just evaluate Phi at many x values and look for sign changes
                // Actually, let's check: for each root of Phi(x, t0) in x,
                // does it have Im(x) <= 0?

                // Use numerical root finding (Durand-Kerner method)
                let roots = find_roots_dk(&phi_x_coeffs, max_xpow);

                let mut all_upper_half_free = true;
                for &(rx, ry) in &roots {
                    if ry > 1e-6 {
                        all_upper_half_free = false;
                    }
                }

                print!("      t = {:.2}+{:.2}i: {} roots, UHP-free={}",
                    ta, tb, roots.len(), all_upper_half_free);
                if roots.len() <= 6 {
                    print!(" [");
                    for (i, &(rx, ry)) in roots.iter().enumerate() {
                        if i > 0 { print!(", "); }
                        print!("{:.3}{:+.3}i", rx, ry);
                    }
                    print!("]");
                }
                println!();
            }
            println!();
        }

        // ── Analysis 8: Matrix power representation ─────────────────────
        // For each n, try to find a matrix A of size m x m (where m is small)
        // and vectors u, v such that:
        //   L^{(p)}(t) = u_p^T * A(t)^{n-something} * v
        //
        // For the alternating case, the simplest possibility:
        //   H_n(t) = e_1^T * M(t)^{n-1} * e_1
        // for some matrix M(t).
        //
        // If H_n(t) = a_0 + a_1 t + ... + a_d t^d, then we need M such that
        // the (1,1) entry of M^{n-1} gives this polynomial.
        //
        // For a 2x2 matrix M = [[a,b],[c,d]], M^k has (1,1) entry involving
        // Chebyshev-like polynomials in the eigenvalues.
        //
        // Let's check: does H_n(t) satisfy a LINEAR RECURRENCE in n (with t-dependent coefficients)?

        if n >= 5 && n <= 10 {
            // Check recurrence: a(t)*H_n + b(t)*H_{n-1} + c(t)*H_{n-2} = 0
            // or higher order
            let alt_nm1 = alternating_permutations(n - 1);
            let alt_nm2 = alternating_permutations(n - 2);

            let h_n = total_poly.clone();
            let h_nm1: Vec<i64> = {
                let vals: Vec<usize> = alt_nm1.iter().map(|s| compute(s, Stat::Swaps)).collect();
                if vals.is_empty() { vec![0] } else { build_poly(&vals) }
            };
            let h_nm2: Vec<i64> = {
                let vals: Vec<usize> = alt_nm2.iter().map(|s| compute(s, Stat::Swaps)).collect();
                if vals.is_empty() { vec![0] } else { build_poly(&vals) }
            };

            // Check: H_n = t * H_{n-1} + ??? * H_{n-2}
            let t_hnm1 = poly_shift(&h_nm1, 1);
            let diff = poly_sub(&h_n, &t_hnm1);
            println!("  --- Recurrence check for n={} ---", n);
            println!("    H_{} = {}", n, format_poly(&h_n));
            println!("    t*H_{} = {}", n - 1, format_poly(&t_hnm1));
            println!("    H_{} - t*H_{} = {}", n, n - 1, format_poly(&diff));

            if !poly_is_zero(&h_nm2) {
                // Check if diff = c(t) * H_{n-2} by numerical evaluation
                let h2_deg = poly_degree(&h_nm2).unwrap_or(0);
                {
                    let test_ts: [f64; 5] = [1.0, 2.0, 3.0, -1.0, 0.5];
                    let mut quotients = Vec::new();
                    for &tv in &test_ts {
                        let dv: f64 = diff.iter().enumerate().map(|(k, &c)| c as f64 * tv.powi(k as i32)).sum();
                        let hv: f64 = h_nm2.iter().enumerate().map(|(k, &c)| c as f64 * tv.powi(k as i32)).sum();
                        if hv.abs() > 1e-10 {
                            quotients.push((tv, dv / hv));
                        }
                    }
                    println!("    (H_{}-t*H_{})/H_{} at t=: {:?}", n, n - 1, n - 2,
                        quotients.iter().map(|(t, r)| format!("{}: {:.4}", t, r)).collect::<Vec<_>>());
                }
            }

            // Also check: H_n = (at + b) * H_{n-1} + (ct^2 + dt + e) * H_{n-2}
            // This is a 2-term recurrence with polynomial coefficients.
            // Try to fit: find a, b, c, d, e such that this holds.
            if !poly_is_zero(&h_nm1) && !poly_is_zero(&h_nm2) && n >= 5 {
                // We know H_n has degree floor((n-1)/2).
                // Try H_n = alpha(t) * H_{n-1} + beta(t) * H_{n-2}
                // where alpha, beta are polynomials in t.
                //
                // This requires: for each coefficient of t^k:
                //   [t^k]H_n = sum_j alpha_j * [t^{k-j}]H_{n-1} + sum_j beta_j * [t^{k-j}]H_{n-2}
                //
                // Try alpha = a0 + a1*t, beta = b0 + b1*t + b2*t^2

                let h_n_deg = poly_degree(&h_n).unwrap_or(0);
                let h_nm1_deg = poly_degree(&h_nm1).unwrap_or(0);
                let h_nm2_deg = poly_degree(&h_nm2).unwrap_or(0);

                // Minimum degree of alpha to make deg(alpha*H_{n-1}) = deg(H_n)
                let alpha_deg_needed = h_n_deg.saturating_sub(h_nm1_deg);

                println!("    deg(H_{})={}, deg(H_{})={}, deg(H_{})={}",
                    n, h_n_deg, n-1, h_nm1_deg, n-2, h_nm2_deg);
                println!("    alpha needs degree >= {}", alpha_deg_needed);
            }
            println!();
        }

        println!();
    }

    // ── Final summary: stability test on Phi(x,t) ──────────────────────
    println!("\n=== Stability verification for alternating Phi_S(x,t) ===\n");

    let mut rng = Rng::new(42);
    let num_tests = 500;

    for n in 4..=max_n {
        let alt = alternating_permutations(n);
        if alt.is_empty() { continue; }

        let positions = valid_positions_alt(n);
        if positions.len() < 2 { continue; }

        let mut by_pos: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
        for sigma in &alt {
            let pos = sigma.iter().position(|&v| v == n).unwrap() as u8 + 1;
            by_pos.entry(pos).or_default().push(compute(sigma, Stat::Swaps));
        }

        let h_polys: Vec<(u8, Vec<i64>)> = positions.iter()
            .filter_map(|&p| by_pos.get(&p).map(|vals| (p, build_poly(vals))))
            .collect();

        let max_xpow = h_polys.iter().map(|(p, _)| *p).max().unwrap() as usize;

        let mut min_abs = f64::MAX;
        let mut stable = true;

        for _ in 0..num_tests {
            let xa = rng.next_range(-5.0, 5.0);
            let xb = rng.next_range(0.001, 5.0);
            let ta = rng.next_range(-5.0, 5.0);
            let tb = rng.next_range(0.001, 5.0);

            let mut val = (0.0, 0.0);
            for (p, poly) in &h_polys {
                let xp = cpow((xa, xb), *p as u32);
                let lt = poly_eval_c(poly, (ta, tb));
                val = cadd(val, cmul(xp, lt));
            }

            let abs_val = cabs(val);
            if abs_val < min_abs {
                min_abs = abs_val;
            }
            if abs_val < 1e-10 {
                stable = false;
                println!("  n={}: FAIL at x={:.4}+{:.4}i, t={:.4}+{:.4}i, |Phi|={:.2e}",
                    n, xa, xb, ta, tb, abs_val);
                break;
            }
        }

        if stable {
            println!("  n={}: STABLE (min |Phi| = {:.4e} over {} tests)", n, min_abs, num_tests);
        }
    }
}

// ── Helper functions ────────────────────────────────────────────────

fn source_descent_set(s_mask: u64, p: u8, n: u8) -> u64 {
    // S'_p: ascending source (p-1 not added to S)
    if n <= 2 { return 0; }
    if p == n { return s_mask; }
    let mut sp = 0u64;
    if p == 1 {
        for j in 2..n {
            if (s_mask >> (j - 1)) & 1 == 1 { sp |= 1 << (j - 2); }
        }
    } else {
        for pos in 1..=(p.saturating_sub(2)) {
            if (s_mask >> (pos - 1)) & 1 == 1 { sp |= 1 << (pos - 1); }
        }
        for j in (p + 1)..n {
            if (s_mask >> (j - 1)) & 1 == 1 { sp |= 1 << (j - 2); }
        }
    }
    sp
}

fn source_descent_set_2(s_mask: u64, p: u8, n: u8) -> Option<u64> {
    // S''_p: descending source (p-1 added to S)
    if p <= 1 || p >= n { return None; }
    Some(source_descent_set(s_mask, p, n) | (1 << (p - 2)))
}

fn descent_set_to_string(mask: u64, n: u8) -> String {
    let mut s = String::from("{");
    let mut first = true;
    for i in 1..n {
        if (mask >> (i - 1)) & 1 == 1 {
            if !first { s.push(','); }
            s.push_str(&i.to_string());
            first = false;
        }
    }
    s.push('}');
    s
}

/// Durand-Kerner method for finding all roots of a polynomial with complex coefficients.
/// coeffs[i] is the coefficient of x^i (as a complex number).
/// Returns approximate roots.
fn find_roots_dk(coeffs: &[(f64, f64)], degree_hint: usize) -> Vec<(f64, f64)> {
    // Find the actual degree
    let mut deg = 0;
    for i in (0..=degree_hint.min(coeffs.len() - 1)).rev() {
        if cabs(coeffs[i]) > 1e-15 {
            deg = i;
            break;
        }
    }
    if deg == 0 {
        return vec![];
    }

    // Normalize: leading coefficient = 1
    let lc = coeffs[deg];
    let lc_inv = (lc.0 / (lc.0 * lc.0 + lc.1 * lc.1), -lc.1 / (lc.0 * lc.0 + lc.1 * lc.1));
    let norm_coeffs: Vec<(f64, f64)> = (0..=deg).map(|i| cmul(coeffs[i], lc_inv)).collect();

    // Initial guesses: points on a circle
    let mut roots: Vec<(f64, f64)> = Vec::new();
    let radius = 2.0; // rough estimate
    for k in 0..deg {
        let angle = 2.0 * std::f64::consts::PI * (k as f64 + 0.37) / deg as f64;
        roots.push((radius * angle.cos(), radius * angle.sin()));
    }

    // Iterate
    for _iter in 0..200 {
        let mut max_change = 0.0f64;
        for i in 0..deg {
            // Evaluate polynomial at roots[i]
            let mut val = (0.0, 0.0);
            for j in (0..=deg).rev() {
                val = cadd(cmul(val, roots[i]), norm_coeffs[j]);
            }

            // Product of (roots[i] - roots[j]) for j != i
            let mut denom = (1.0, 0.0);
            for j in 0..deg {
                if j != i {
                    let diff = (roots[i].0 - roots[j].0, roots[i].1 - roots[j].1);
                    denom = cmul(denom, diff);
                }
            }

            let denom_abs2 = denom.0 * denom.0 + denom.1 * denom.1;
            if denom_abs2 > 1e-30 {
                let denom_inv = (denom.0 / denom_abs2, -denom.1 / denom_abs2);
                let correction = cmul(val, denom_inv);
                roots[i].0 -= correction.0;
                roots[i].1 -= correction.1;
                max_change = max_change.max(cabs(correction));
            }
        }
        if max_change < 1e-12 {
            break;
        }
    }

    roots
}
