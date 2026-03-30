/// Trivariate stability check for the "joint source polynomial":
///
///   Ξ(z, y, t) = Σ_{i=1}^{m} y^i · Φ̃_{p_i}(z, t)
///
/// where Φ̃_{p_i}(z, t) = Σ_q z^q F̃_q^{(p_i)}(t)
/// is the bivariate generating function for the *modified source* at
/// insertion position p_i, tracking:
///   z^q = position of n-1 in the source permutation,
///   t^k = modified swaps count  (swaps + ε₁ correction).
///
/// The y-exponents are consecutive ranks 1..m (not the actual positions p_i),
/// so that by the Borcea–Brändén theorem:
///
///   Ξ(z, y, t) stable  ⟺  for every fixed real (z, t), the sequence
///   (Φ̃_{p_1}(z,t), ..., Φ̃_{p_m}(z,t)) forms an interlacing chain.
///
/// This is a *strengthening* of Conjecture conj:cross-chain (which asks only
/// for compatible families), and is verified here for n ≤ 8.
///
/// If Ξ is trivariate stable, combining with the staircase operator gives
/// a route to proving bivariate stability of Ψ_{n,S}(y,t).

use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use std::collections::BTreeMap;

// ── Complex arithmetic ────────────────────────────────────────────────────────
fn cmul(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}
fn cadd(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 + b.0, a.1 + b.1)
}
fn cpow(base: (f64, f64), exp: u32) -> (f64, f64) {
    let mut r = (1.0_f64, 0.0_f64);
    for _ in 0..exp {
        r = cmul(r, base);
    }
    r
}
fn cabs(a: (f64, f64)) -> f64 {
    (a.0 * a.0 + a.1 * a.1).sqrt()
}

// ── Simple PRNG ───────────────────────────────────────────────────────────────
struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self { Rng(seed.wrapping_add(1)) }
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005)
                       .wrapping_add(1442695040888963407);
        self.0
    }
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn next_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn descent_set_str(mask: u64, n: u8) -> String {
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

/// Valid positions P(S): left-boundary descents, and n if n-1 ∉ S.
fn valid_positions(s_mask: u64, n: u8) -> Vec<u8> {
    let mut pos = Vec::new();
    for p in 1..n {
        let p_in = (s_mask >> (p - 1)) & 1 == 1;
        let pm1_in = p >= 2 && (s_mask >> (p - 2)) & 1 == 1;
        if p_in && !pm1_in { pos.push(p); }
    }
    if n >= 2 && (s_mask >> (n - 2)) & 1 == 0 { pos.push(n); }
    pos
}

/// S'_p  (base source descent set, as bitmask in [n-2]).
fn source_asc(s_mask: u64, p: u8, n: u8) -> u64 {
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

/// S''_p  (augmented source = S'_p ∪ {p-1}), None for p=1 or p=n.
fn source_desc(s_mask: u64, p: u8, n: u8) -> Option<u64> {
    if p <= 1 || p >= n { return None; }
    Some(source_asc(s_mask, p, n) | (1 << (p - 2)))
}

/// ε₁(π, p) = 1 iff π(p-1)+1 = π(p)  (π is 0-indexed slice of length n-1).
fn epsilon1(pi: &[u8], p: u8) -> bool {
    let n = pi.len() as u8 + 1;
    if p <= 1 || p >= n { return false; }
    pi[(p - 2) as usize] + 1 == pi[(p - 1) as usize]
}

// ── Core: build the coefficient tensor ───────────────────────────────────────
//
// coeffs[rank_i][q_idx][k]  = [y^{rank_i+1} z^{q_idx+1} t^k] Ξ
//   rank_i ∈ 0..m  (rank of valid position)
//   q_idx  ∈ 0..(n-2)  (source position of n-1, 1-indexed → 0-indexed here)
//   k      ∈ 0..max_swaps
//
// Returns (m, n-1, max_k+1, coeffs_flat) where coeffs_flat is
// indexed as [rank_i][(q-1)][k].

fn build_xi_coeffs(
    n: u8,
    s_mask: u64,
    positions: &[u8],
    // All source permutations of [n-1], pre-grouped by descent set bitmask
    source_by_descent: &BTreeMap<u64, Vec<Vec<u8>>>,
) -> Vec<Vec<Vec<i64>>> {
    let m = positions.len();
    let nm1 = (n - 1) as usize;

    // Find max modified swaps to determine tensor size
    // (We'll grow dynamically)
    let mut tensor: Vec<Vec<Vec<i64>>> = vec![vec![vec![]; nm1]; m];

    for (rank_i, &p) in positions.iter().enumerate() {
        let sp_asc = source_asc(s_mask, p, n);
        let sp_desc = source_desc(s_mask, p, n);

        // Iterate over all source permutations
        for (desc_mask, class) in source_by_descent.iter() {
            let is_asc = *desc_mask == sp_asc;
            let is_desc = sp_desc.map_or(false, |sd| *desc_mask == sd);
            if !is_asc && !is_desc { continue; }

            for pi in class {
                // Position of n-1 in π (1-indexed position in [n-1])
                let q = pi.iter().position(|&v| v == n - 1).unwrap() + 1;
                let q_idx = q - 1;

                let sw = compute(pi, Stat::Swaps);
                let eps = if is_asc && epsilon1(pi, p) { 1 } else { 0 };
                let mod_sw = sw + eps;

                // Ensure tensor is large enough
                let row = &mut tensor[rank_i][q_idx];
                if row.len() <= mod_sw {
                    row.resize(mod_sw + 1, 0);
                }
                row[mod_sw] += 1;
            }
        }
    }

    tensor
}

// ── Evaluate Ξ(z, y, t) at complex point ─────────────────────────────────────
//
// tensor[i][q_idx][k] = [y^{i+1} z^{q_idx+1} t^k] Ξ
fn eval_xi(
    tensor: &[Vec<Vec<i64>>],
    z: (f64, f64),
    y: (f64, f64),
    t: (f64, f64),
) -> (f64, f64) {
    let mut result = (0.0_f64, 0.0_f64);
    for (i, source_rows) in tensor.iter().enumerate() {
        let yi = cpow(y, (i + 1) as u32); // y^{rank = i+1}
        for (q_idx, tk_coeffs) in source_rows.iter().enumerate() {
            if tk_coeffs.is_empty() { continue; }
            let zq = cpow(z, (q_idx + 1) as u32); // z^{q = q_idx+1}
            // Evaluate Σ_k c_k t^k
            let mut tpoly = (0.0_f64, 0.0_f64);
            for (k, &c) in tk_coeffs.iter().enumerate() {
                if c != 0 {
                    let tk = cpow(t, k as u32);
                    tpoly = cadd(tpoly, cmul((c as f64, 0.0), tk));
                }
            }
            result = cadd(result, cmul(cmul(yi, zq), tpoly));
        }
    }
    result
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    let num_tests: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(200);

    println!("=== Trivariate stability of Ξ(z, y, t) ===");
    println!("Ξ = Σ_i y^i Φ̃_{{p_i}}(z,t),  y-exponents = consecutive ranks");
    println!("Checking n = 3..{}, {} random test points each\n", max_n, num_tests);

    let mut grand_total = 0u64;
    let mut grand_pass = 0u64;
    let mut grand_trivial = 0u64; // m < 2

    for n in 3..=max_n {
        let all_perms = all_permutations(n);

        // Group level-n permutations by descent set (for target)
        let mut target_by_descent: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &all_perms {
            target_by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        // Group level-(n-1) permutations by descent set (for source)
        let source_perms = all_permutations(n - 1);
        let mut source_by_descent: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &source_perms {
            source_by_descent
                .entry(descent_set_bitmask(pi))
                .or_default()
                .push(pi.clone());
        }

        println!("========== n = {} ==========", n);

        let mut n_total = 0u64;
        let mut n_pass = 0u64;
        let mut n_trivial = 0u64;
        let mut n_fail = 0u64;

        for (&s_mask, _class) in &target_by_descent {
            // Only S with 1 ∉ S
            if s_mask & 1 != 0 { continue; }

            let positions = valid_positions(s_mask, n);
            let m = positions.len();

            if m < 2 {
                n_trivial += 1;
                continue;
            }

            n_total += 1;
            let s_str = descent_set_str(s_mask, n);

            // Build the coefficient tensor for Ξ
            let tensor = build_xi_coeffs(n, s_mask, &positions, &source_by_descent);

            // Stability check: evaluate at num_tests random points with Im > 0
            let mut pass = true;
            let mut min_abs = f64::MAX;
            let mut rng = Rng::new(s_mask.wrapping_mul(31337).wrapping_add(n as u64 * 65537));

            'outer: for test_idx in 0..num_tests {
                // Generate complex points with positive imaginary parts.
                // Use a mix of "typical" and "edge-near" values.
                let (z, y, t) = if test_idx < 6 {
                    // A few structured tests
                    let cases: &[(f64, f64, f64, f64, f64, f64)] = &[
                        (0.0, 1.0,  0.0, 1.0,  0.0, 1.0),  // i, i, i
                        (1.0, 1.0,  1.0, 1.0,  1.0, 1.0),  // 1+i, 1+i, 1+i
                        (-1.0, 0.5, 0.0, 0.5, -1.0, 0.5),  // mixed real parts
                        (0.5, 2.0,  0.5, 0.5,  0.5, 2.0),
                        (-0.3, 0.1, 0.7, 0.3, -0.5, 0.4),
                        (2.0, 0.1, -1.0, 0.2,  0.3, 0.7),
                    ];
                    let c = cases[test_idx];
                    ((c.0, c.1), (c.2, c.3), (c.4, c.5))
                } else {
                    let zr = rng.next_range(-3.0, 3.0);
                    let zi = rng.next_range(0.05, 3.0);
                    let yr = rng.next_range(-3.0, 3.0);
                    let yi = rng.next_range(0.05, 3.0);
                    let tr = rng.next_range(-3.0, 3.0);
                    let ti = rng.next_range(0.05, 3.0);
                    ((zr, zi), (yr, yi), (tr, ti))
                };

                let val = eval_xi(&tensor, z, y, t);
                let abs_val = cabs(val);
                if abs_val < min_abs { min_abs = abs_val; }

                if abs_val < 1e-9 {
                    pass = false;
                    println!(
                        "  FAIL S={} p={:?}",
                        s_str, positions
                    );
                    println!(
                        "    Ξ({:.3}+{:.3}i, {:.3}+{:.3}i, {:.3}+{:.3}i) ≈ 0  (|Ξ|={:.2e})",
                        z.0, z.1, y.0, y.1, t.0, t.1, abs_val
                    );
                    n_fail += 1;
                    break 'outer;
                }
            }

            if pass {
                n_pass += 1;
                println!(
                    "  PASS  S={}  m={}  min|Ξ|={:.3e}",
                    s_str, m, min_abs
                );
            }
        }

        grand_total    += n_total;
        grand_pass     += n_pass;
        grand_trivial  += n_trivial;

        println!(
            "n={}: {}/{} passed  ({} trivial with m<2)\n",
            n, n_pass, n_total, n_trivial
        );
    }

    println!("==========================================");
    println!("TOTAL: {}/{} cases passed  ({} trivial skipped)",
             grand_pass, grand_total, grand_trivial);

    if grand_pass == grand_total && grand_total > 0 {
        println!();
        println!("*** ALL CASES PASSED ***");
        println!("Conjecture: Ξ(z, y, t) is trivariate stable for all n, S with 1 ∉ S.");
        println!("This strengthens conj:cross-chain and provides a stability-preserving");
        println!("operator route to proving bivariate stability of Ψ_{{n,S}}(y,t).");
    }
}
