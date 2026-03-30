/// Trivariate stability check for Υ(x, z, t):
///
///   Υ(x, z, t) = Σ_{σ ∈ D(n,S)} x^{σ⁻¹(n)} z^{σ⁻¹(n-1)} t^{swaps(σ)}
///
/// This tracks: x = position of n, z = position of n-1, t = swaps weight.
///
/// Key properties:
///   - z=1 specialization gives Φ_{n,S}(x,t) (bivariate stability conj:stable)
///   - Trivariate stability ⟹ bivariate stability of Φ via real specialization
///
/// Operator formulation:
///   Υ_{n,S}(x,z,t) = Σ_{p ∈ P(S)} x^p · M_p[Q_p](z,t)
///
/// where Q_p(z,t) = Σ_q z^q F̃_q^{(p)}(t) is the source bivariate poly and
///   M_p[Q](z,t) = z·Q(z,t) + (t-z)·tr_{p-2}(Q)(z,t) + (1-z)·z^{p-1}·[z^{p-1}]Q(z,t)
///
/// The operator acts on monomials z^q t^k as:
///   q ≤ p-2  (Zone L): z^q t^k  →  z^q t^{k+1}   (multiply by t)
///   q = p-1  (Zone M): z^{p-1} t^k  →  z^{p-1} t^k  (unchanged)
///   q ≥ p    (Zone R): z^q t^k  →  z^{q+1} t^k    (shift z-degree by 1)
///
/// Two tests are performed:
///   1. Direct: Build Υ directly from target permutations, check trivariate stability.
///   2. Operator: Build via M_p applied to source polynomials, verify it matches Υ,
///      then check if individual M_p[Q_p] are bivariate stable.
///   3. Operator generality: Check if M_p maps *arbitrary* stable polynomials to stable.

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
fn csub(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    (a.0 - b.0, a.1 - b.1)
}
fn cpow(base: (f64, f64), exp: u32) -> (f64, f64) {
    let mut r = (1.0_f64, 0.0_f64);
    for _ in 0..exp { r = cmul(r, base); }
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

fn source_desc(s_mask: u64, p: u8, n: u8) -> Option<u64> {
    if p <= 1 || p >= n { return None; }
    Some(source_asc(s_mask, p, n) | (1 << (p - 2)))
}

fn epsilon1(pi: &[u8], p: u8) -> bool {
    let n = pi.len() as u8 + 1;
    if p <= 1 || p >= n { return false; }
    pi[(p - 2) as usize] + 1 == pi[(p - 1) as usize]
}

// ── Coefficient types ─────────────────────────────────────────────────────────
//
// A bivariate polynomial in (z,t) is represented as a 2D array:
//   coeffs[q_idx][k] = coefficient of z^{q_idx+1} t^k
// (q is 1-indexed: position of element, so q ∈ {1,...,n-1})
//
// For Υ(x,z,t), we add an outer x-dimension:
//   upsilon[p_idx][q_idx][k] = [x^p z^q t^k] Υ
// where p ∈ {1,...,n} (1-indexed position of n)

// ── Build Υ directly from target permutations ─────────────────────────────────

fn build_upsilon_direct(
    n: u8,
    s_mask: u64,
    all_perms: &[Vec<u8>],
) -> Vec<Vec<Vec<i64>>> {
    // upsilon[p-1][q-1][k] = [x^p z^q t^k] Υ
    let n_usize = n as usize;
    let mut upsilon: Vec<Vec<Vec<i64>>> = vec![vec![vec![]; n_usize]; n_usize];

    for pi in all_perms {
        if descent_set_bitmask(pi) != s_mask { continue; }
        let p = pi.iter().position(|&v| v == n).unwrap() + 1;  // pos of n (1-indexed)
        let q = pi.iter().position(|&v| v == n - 1).unwrap() + 1;  // pos of n-1 (1-indexed)
        let sw = compute(pi, Stat::Swaps);

        let row = &mut upsilon[p - 1][q - 1];
        if row.len() <= sw { row.resize(sw + 1, 0); }
        row[sw] += 1;
    }

    upsilon
}

// ── Build Υ via operator M_p from source polynomials ──────────────────────────
//
// Q_p(z,t) = Σ_q z^q F̃_q^{(p)}(t) where F̃_q^{(p)} uses modified swaps.
// M_p[Q](z,t): zone L → multiply by t, zone M → unchanged, zone R → shift z by 1.
// Result: [x^p coefficient of Υ] = M_p[Q_p](z,t).

fn build_q_source(
    n: u8,
    p: u8,
    s_mask: u64,
    source_by_descent: &BTreeMap<u64, Vec<Vec<u8>>>,
) -> Vec<Vec<i64>> {
    // q_poly[q-1][k] = coeff of z^q t^k in Q_p(z,t)
    let nm1 = (n - 1) as usize;
    let mut q_poly: Vec<Vec<i64>> = vec![vec![]; nm1];

    let sp_asc = source_asc(s_mask, p, n);
    let sp_desc = source_desc(s_mask, p, n);

    for (desc_mask, class) in source_by_descent.iter() {
        let is_asc = *desc_mask == sp_asc;
        let is_desc = sp_desc.map_or(false, |sd| *desc_mask == sd);
        if !is_asc && !is_desc { continue; }

        for pi in class {
            let q = pi.iter().position(|&v| v == n - 1).unwrap() + 1;
            let sw = compute(pi, Stat::Swaps);
            let eps = if is_asc && epsilon1(pi, p) { 1 } else { 0 };
            let mod_sw = sw + eps;

            let row = &mut q_poly[q - 1];
            if row.len() <= mod_sw { row.resize(mod_sw + 1, 0); }
            row[mod_sw] += 1;
        }
    }

    q_poly
}

fn apply_mp(q_poly: &[Vec<i64>], p: u8) -> Vec<Vec<i64>> {
    // Apply M_p: zone L (q ≤ p-2) → multiply t, zone M (q=p-1) → unchanged,
    // zone R (q ≥ p) → shift z-index by 1 (q → q+1).
    let n_minus_1 = q_poly.len();
    // Output has one more z-slot (from zone R shift)
    let mut out: Vec<Vec<i64>> = vec![vec![]; n_minus_1 + 1];

    for (q_idx, t_coeffs) in q_poly.iter().enumerate() {
        let q = q_idx + 1;  // q is 1-indexed
        if t_coeffs.is_empty() { continue; }

        if p >= 2 && q <= (p - 2) as usize {
            // Zone L: multiply by t → shift k by 1 (prepend 0)
            let out_q_idx = q_idx;
            let row = &mut out[out_q_idx];
            let needed = t_coeffs.len() + 1;
            if row.len() < needed { row.resize(needed, 0); }
            for (k, &c) in t_coeffs.iter().enumerate() {
                row[k + 1] += c;
            }
        } else if q == (p as usize).saturating_sub(1) {
            // Zone M: unchanged (q = p-1, out_q_idx = q_idx)
            let out_q_idx = q_idx;
            let row = &mut out[out_q_idx];
            if row.len() < t_coeffs.len() { row.resize(t_coeffs.len(), 0); }
            for (k, &c) in t_coeffs.iter().enumerate() {
                row[k] += c;
            }
        } else {
            // Zone R: shift z-index by 1 (q → q+1, so out_q_idx = q_idx + 1)
            let out_q_idx = q_idx + 1;
            if out_q_idx >= out.len() { out.resize(out_q_idx + 1, vec![]); }
            let row = &mut out[out_q_idx];
            if row.len() < t_coeffs.len() { row.resize(t_coeffs.len(), 0); }
            for (k, &c) in t_coeffs.iter().enumerate() {
                row[k] += c;
            }
        }
    }

    out
}

fn build_upsilon_operator(
    n: u8,
    s_mask: u64,
    positions: &[u8],
    source_by_descent: &BTreeMap<u64, Vec<Vec<u8>>>,
) -> Vec<Vec<Vec<i64>>> {
    let n_usize = n as usize;
    let mut upsilon: Vec<Vec<Vec<i64>>> = vec![vec![]; n_usize];

    for &p in positions {
        let q_src = build_q_source(n, p, s_mask, source_by_descent);
        let mp_result = apply_mp(&q_src, p);
        upsilon[p as usize - 1] = mp_result;
    }

    upsilon
}

// ── Evaluate Υ(x,z,t) at complex point ───────────────────────────────────────

fn eval_upsilon(
    upsilon: &[Vec<Vec<i64>>],
    x: (f64, f64),
    z: (f64, f64),
    t: (f64, f64),
) -> (f64, f64) {
    let mut result = (0.0_f64, 0.0_f64);
    for (p_idx, z_rows) in upsilon.iter().enumerate() {
        if z_rows.is_empty() { continue; }
        let xp = cpow(x, (p_idx + 1) as u32);
        for (q_idx, t_coeffs) in z_rows.iter().enumerate() {
            if t_coeffs.is_empty() { continue; }
            let zq = cpow(z, (q_idx + 1) as u32);
            let mut tpoly = (0.0_f64, 0.0_f64);
            for (k, &c) in t_coeffs.iter().enumerate() {
                if c != 0 {
                    tpoly = cadd(tpoly, cmul((c as f64, 0.0), cpow(t, k as u32)));
                }
            }
            result = cadd(result, cmul(cmul(xp, zq), tpoly));
        }
    }
    result
}

// Evaluate bivariate poly (z_rows[q_idx][k] = coeff of z^{q+1} t^k)
fn eval_bivariate(
    z_rows: &[Vec<i64>],
    z: (f64, f64),
    t: (f64, f64),
) -> (f64, f64) {
    let mut result = (0.0_f64, 0.0_f64);
    for (q_idx, t_coeffs) in z_rows.iter().enumerate() {
        if t_coeffs.is_empty() { continue; }
        let zq = cpow(z, (q_idx + 1) as u32);
        let mut tpoly = (0.0_f64, 0.0_f64);
        for (k, &c) in t_coeffs.iter().enumerate() {
            if c != 0 {
                tpoly = cadd(tpoly, cmul((c as f64, 0.0), cpow(t, k as u32)));
            }
        }
        result = cadd(result, cmul(zq, tpoly));
    }
    result
}

// ── Upsilon equality check ────────────────────────────────────────────────────

fn upsilon_equal(a: &[Vec<Vec<i64>>], b: &[Vec<Vec<i64>>]) -> bool {
    let max_p = a.len().max(b.len());
    for p_idx in 0..max_p {
        let a_rows = if p_idx < a.len() { &a[p_idx] } else { &[] as &[Vec<i64>] };
        let b_rows = if p_idx < b.len() { &b[p_idx] } else { &[] as &[Vec<i64>] };
        let max_q = a_rows.len().max(b_rows.len());
        for q_idx in 0..max_q {
            let a_tk = if q_idx < a_rows.len() { &a_rows[q_idx] } else { &[] as &[i64] };
            let b_tk = if q_idx < b_rows.len() { &b_rows[q_idx] } else { &[] as &[i64] };
            let max_k = a_tk.len().max(b_tk.len());
            for k in 0..max_k {
                let av = if k < a_tk.len() { a_tk[k] } else { 0 };
                let bv = if k < b_tk.len() { b_tk[k] } else { 0 };
                if av != bv { return false; }
            }
        }
    }
    true
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);

    let num_tests: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(300);

    println!("=== Trivariate stability of Υ(x,z,t) ===");
    println!("Υ = Σ_{{σ∈D(n,S)}} x^{{pos(n)}} z^{{pos(n-1)}} t^{{swaps(σ)}}");
    println!("Checking n=3..{}, {} random test points each\n", max_n, num_tests);

    let mut grand_total  = 0u64;
    let mut grand_pass   = 0u64;
    let mut grand_match  = 0u64;  // operator == direct

    for n in 3..=max_n {
        let all_perms = all_permutations(n);

        // Group level-n perms by descent set
        let mut by_descent: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &all_perms {
            by_descent.entry(descent_set_bitmask(pi)).or_default().push(pi.clone());
        }

        // Group level-(n-1) perms by descent set (for source)
        let source_perms = all_permutations(n - 1);
        let mut source_by_descent: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &source_perms {
            source_by_descent.entry(descent_set_bitmask(pi)).or_default().push(pi.clone());
        }

        println!("========== n = {} ==========", n);

        let mut n_total = 0u64;
        let mut n_pass  = 0u64;
        let mut n_match = 0u64;

        for (&s_mask, _class) in &by_descent {
            if s_mask & 1 != 0 { continue; }  // only S with 1 ∉ S

            let positions = valid_positions(s_mask, n);
            if positions.len() < 1 { continue; }

            n_total += 1;
            let s_str = descent_set_str(s_mask, n);

            // Build Υ both ways
            let ups_direct   = build_upsilon_direct(n, s_mask, &all_perms);
            let ups_operator = build_upsilon_operator(n, s_mask, &positions, &source_by_descent);

            // Verify they match
            let matched = upsilon_equal(&ups_direct, &ups_operator);
            if matched { n_match += 1; }

            // Stability check: evaluate at num_tests random points with Im > 0
            let mut pass = true;
            let mut min_abs = f64::MAX;
            let mut rng = Rng::new(s_mask.wrapping_mul(31337).wrapping_add(n as u64 * 65537));

            'outer: for test_idx in 0..num_tests {
                let (x, z, t) = if test_idx < 6 {
                    let cases: &[(f64,f64, f64,f64, f64,f64)] = &[
                        (0.0, 1.0,  0.0, 1.0,  0.0, 1.0),
                        (1.0, 1.0,  1.0, 1.0,  1.0, 1.0),
                        (-1.0, 0.5, 0.0, 0.5, -1.0, 0.5),
                        (0.5, 2.0,  0.5, 0.5,  0.5, 2.0),
                        (-0.3, 0.1, 0.7, 0.3, -0.5, 0.4),
                        (2.0, 0.1, -1.0, 0.2,  0.3, 0.7),
                    ];
                    let c = cases[test_idx];
                    ((c.0, c.1), (c.2, c.3), (c.4, c.5))
                } else {
                    let xr = rng.next_range(-3.0, 3.0);
                    let xi = rng.next_range(0.05, 3.0);
                    let zr = rng.next_range(-3.0, 3.0);
                    let zi = rng.next_range(0.05, 3.0);
                    let tr = rng.next_range(-3.0, 3.0);
                    let ti = rng.next_range(0.05, 3.0);
                    ((xr, xi), (zr, zi), (tr, ti))
                };

                let val = eval_upsilon(&ups_direct, x, z, t);
                let abs_val = cabs(val);
                if abs_val < min_abs { min_abs = abs_val; }

                if abs_val < 1e-9 {
                    pass = false;
                    println!("  FAIL S={} pos={:?}", s_str, positions);
                    println!("    Υ({:.3}+{:.3}i, {:.3}+{:.3}i, {:.3}+{:.3}i) ≈ 0  (|Υ|={:.2e})",
                        x.0, x.1, z.0, z.1, t.0, t.1, abs_val);
                    break 'outer;
                }
            }

            if pass { n_pass += 1; }

            let match_str = if matched { "✓" } else { "MISMATCH" };
            println!("  S={}  pos={:?}  {}  {}  min|Υ|={:.3e}",
                s_str, positions,
                if pass { "PASS" } else { "FAIL" },
                match_str,
                min_abs);
        }

        grand_total += n_total;
        grand_pass  += n_pass;
        grand_match += n_match;

        println!("n={}: {}/{} passed, {}/{} operator==direct\n",
            n, n_pass, n_total, n_match, n_total);
    }

    println!("==========================================");
    println!("TOTAL: {}/{} trivariate stable, {}/{} operator matches",
             grand_pass, grand_total, grand_match, grand_total);

    if grand_pass == grand_total && grand_total > 0 {
        println!();
        println!("*** ALL CASES PASSED ***");
        println!("Conjecture: Υ(x,z,t) is trivariate stable for all n, S with 1 ∉ S.");
        println!("This implies bivariate stability of Φ_{{n,S}}(x,t) via z↦1.");
    }

    // ── Separate: test if M_p maps stable polynomials to stable ──────────────
    println!();
    println!("=== Operator stability test ===");
    println!("Does M_p map bivariate stable polynomials to bivariate stable?");
    println!("Test: evaluate M_p[Q] at random Im>0 points for various n, S, p.");

    let mut op_total = 0u64;
    let mut op_pass  = 0u64;

    for n in 3..=max_n {
        let source_perms = all_permutations(n - 1);
        let mut source_by_descent: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &source_perms {
            source_by_descent.entry(descent_set_bitmask(pi)).or_default().push(pi.clone());
        }

        let all_perms = all_permutations(n);
        let mut by_descent: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        for pi in &all_perms {
            by_descent.entry(descent_set_bitmask(pi)).or_default().push(pi.clone());
        }

        for (&s_mask, _class) in &by_descent {
            if s_mask & 1 != 0 { continue; }
            let positions = valid_positions(s_mask, n);

            for &p in &positions {
                let q_src = build_q_source(n, p, s_mask, &source_by_descent);
                let mp_q  = apply_mp(&q_src, p);

                op_total += 1;
                let mut pass = true;
                let mut rng = Rng::new(
                    s_mask.wrapping_mul(99991)
                         .wrapping_add(p as u64 * 7)
                         .wrapping_add(n as u64 * 1009)
                );

                'op_outer: for test_idx in 0..num_tests {
                    let (z, t) = if test_idx < 4 {
                        let cases: &[(f64,f64,f64,f64)] = &[
                            (0.0, 1.0, 0.0, 1.0),
                            (1.0, 1.0, 1.0, 1.0),
                            (0.5, 0.5, 0.5, 0.5),
                            (-1.0, 0.3, 0.5, 0.7),
                        ];
                        let c = cases[test_idx];
                        ((c.0, c.1), (c.2, c.3))
                    } else {
                        let zr = rng.next_range(-3.0, 3.0);
                        let zi = rng.next_range(0.05, 3.0);
                        let tr = rng.next_range(-3.0, 3.0);
                        let ti = rng.next_range(0.05, 3.0);
                        ((zr, zi), (tr, ti))
                    };

                    let val = eval_bivariate(&mp_q, z, t);
                    if cabs(val) < 1e-9 {
                        pass = false;
                        let s_str = descent_set_str(s_mask, n);
                        println!("  FAIL M_p[Q_p]: n={} S={} p={}", n, s_str, p);
                        println!("    M_{}[Q]({:.3}+{:.3}i, {:.3}+{:.3}i) ≈ 0",
                            p, z.0, z.1, t.0, t.1);
                        break 'op_outer;
                    }
                }

                if pass { op_pass += 1; }
            }
        }
    }

    println!("M_p[Q_p] bivariate stable: {}/{}", op_pass, op_total);

    if op_pass == op_total {
        println!();
        println!("All M_p[Q_p] appear bivariate stable!");
        println!("This supports: M_p preserves bivariate stability (at least on source polys).");
    }
}
