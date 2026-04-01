/// Does the operator M_p preserve bivariate stability on *arbitrary* stable polynomials?
///
/// M_p acts on z^q t^k as (0-indexed q):
///   Zone L (q ≤ p-2):  z^q t^k → z^q t^{k+1}    (multiply by t)
///   Zone M (q = p-1):  z^{p-1} t^k → z^{p-1} t^k  (unchanged)
///   Zone R (q ≥ p):    z^q t^k → z^{q+1} t^k      (shift z-degree)
///
/// We generate random stable polynomials as products of linear factors
///   (a·z + b·t + c)  with  a, b, c > 0.
/// Such products are stable: Im(a·z + b·t + c) = a·Im(z) + b·Im(t) > 0
/// whenever Im(z), Im(t) > 0, so no zero can occur.
///
/// For each (p, degree, random Q), we check: does M_p[Q] vanish at any
/// randomly sampled point with Im(z), Im(t) > 0?
///
/// We also test the COMPOSITION:  applying M_{p_1} then M_{p_2}, to see
/// whether iterated application remains stable (relevant for the inductive
/// step in proving conj:upsilon-stable).

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
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn next_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

// ── Bivariate polynomial: coeffs[q][k] = [z^q t^k] P, 0-indexed ──────────────
type BivPoly = Vec<Vec<f64>>;

fn biv_zero() -> BivPoly {
    vec![]
}

fn biv_get(p: &BivPoly, q: usize, k: usize) -> f64 {
    if q < p.len() && k < p[q].len() {
        p[q][k]
    } else {
        0.0
    }
}

fn biv_set(p: &mut BivPoly, q: usize, k: usize, v: f64) {
    if q >= p.len() {
        p.resize(q + 1, vec![]);
    }
    if k >= p[q].len() {
        p[q].resize(k + 1, 0.0);
    }
    p[q][k] = v;
}

fn biv_add_to(p: &mut BivPoly, q: usize, k: usize, v: f64) {
    if q >= p.len() {
        p.resize(q + 1, vec![]);
    }
    if k >= p[q].len() {
        p[q].resize(k + 1, 0.0);
    }
    p[q][k] += v;
}

/// Multiply bivariate polynomial by linear factor (a·z + b·t + c), all reals.
fn biv_mul_linear(p: &BivPoly, a: f64, b: f64, c: f64) -> BivPoly {
    let mut out = biv_zero();
    let qmax = p.len();
    for q in 0..qmax {
        let kmax = if q < p.len() { p[q].len() } else { 0 };
        for k in 0..kmax {
            let v = biv_get(p, q, k);
            if v == 0.0 {
                continue;
            }
            // ·c term
            biv_add_to(&mut out, q, k, c * v);
            // ·a·z term (q → q+1)
            biv_add_to(&mut out, q + 1, k, a * v);
            // ·b·t term (k → k+1)
            biv_add_to(&mut out, q, k + 1, b * v);
        }
    }
    out
}

/// Build a random stable polynomial as a product of `deg` linear factors.
/// Factors: (a·z + b·t + c) with a, b, c ∈ (0.1, 3.0).
fn random_stable_poly(rng: &mut Rng, deg: usize) -> BivPoly {
    // Start with the constant 1
    let mut p: BivPoly = vec![vec![1.0]];
    for _ in 0..deg {
        let a = rng.next_range(0.1, 3.0);
        let b = rng.next_range(0.1, 3.0);
        let c = rng.next_range(0.1, 3.0);
        p = biv_mul_linear(&p, a, b, c);
    }
    p
}

/// Apply M_p to a bivariate polynomial (0-indexed z-degree).
fn apply_mp_float(q_poly: &BivPoly, p: usize) -> BivPoly {
    let mut out = biv_zero();
    for (q, t_coeffs) in q_poly.iter().enumerate() {
        for (k, &c) in t_coeffs.iter().enumerate() {
            if c == 0.0 {
                continue;
            }
            if p >= 2 && q <= p - 2 {
                // Zone L: multiply by t → k → k+1, same q
                biv_add_to(&mut out, q, k + 1, c);
            } else if p >= 1 && q == p - 1 {
                // Zone M: unchanged
                biv_add_to(&mut out, q, k, c);
            } else {
                // Zone R (q >= p): shift z → q+1, same k
                biv_add_to(&mut out, q + 1, k, c);
            }
        }
    }
    out
}

/// Evaluate bivariate polynomial at complex (z, t).
fn eval_biv_complex(p: &BivPoly, z: (f64, f64), t: (f64, f64)) -> (f64, f64) {
    let mut result = (0.0_f64, 0.0_f64);
    for (q, t_coeffs) in p.iter().enumerate() {
        if t_coeffs.is_empty() {
            continue;
        }
        let zq = cpow(z, q as u32);
        for (k, &c) in t_coeffs.iter().enumerate() {
            if c == 0.0 {
                continue;
            }
            let tk = cpow(t, k as u32);
            result = cadd(result, cmul((c, 0.0), cmul(zq, tk)));
        }
    }
    result
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    let num_polys: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(2000);
    let num_eval: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(50);
    let max_deg: usize = std::env::args()
        .nth(3)
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);

    println!("=== M_p stability-preservation test ===");
    println!("Random stable polys: products of (az+bt+c) with a,b,c>0");
    println!(
        "{} polynomials × {} eval points, degree 1..{}\n",
        num_polys, num_eval, max_deg
    );

    let mut total_tests = 0u64;
    let mut total_pass = 0u64;
    let mut min_abs_seen = f64::MAX;

    // Test p = 1..max_p for polynomials of various degrees
    for p in 1..=8usize {
        let mut p_tests = 0u64;
        let mut p_pass = 0u64;
        let mut p_min = f64::MAX;

        let mut rng = Rng::new(p as u64 * 1000003 + 42);

        for poly_idx in 0..num_polys {
            let deg = 1 + poly_idx % max_deg; // cycle through degrees 1..max_deg
            let q_poly = random_stable_poly(&mut rng, deg);
            let mp_q = apply_mp_float(&q_poly, p);

            let mut pass = true;
            let mut min_abs = f64::MAX;

            for eval_idx in 0..num_eval {
                // Sample (z, t) with Im > 0
                let (z, t) = if eval_idx < 4 {
                    let cases = [
                        (0.0_f64, 1.0, 0.0, 1.0),
                        (1.0, 1.0, 1.0, 1.0),
                        (0.5, 0.5, 0.5, 0.5),
                        (-1.0, 0.3, 0.5, 0.7),
                    ];
                    let c = cases[eval_idx];
                    ((c.0, c.1), (c.2, c.3))
                } else {
                    let zr = rng.next_range(-5.0, 5.0);
                    let zi = rng.next_range(0.01, 5.0);
                    let tr = rng.next_range(-5.0, 5.0);
                    let ti = rng.next_range(0.01, 5.0);
                    ((zr, zi), (tr, ti))
                };

                let val = eval_biv_complex(&mp_q, z, t);
                let abs_val = cabs(val);
                if abs_val < min_abs {
                    min_abs = abs_val;
                }
                if abs_val < 1e-10 {
                    pass = false;
                    println!("  FAIL p={} deg={} poly#{}: |M_p[Q]|={:.2e} at ({:.3}+{:.3}i, {:.3}+{:.3}i)",
                        p, deg, poly_idx, abs_val, z.0, z.1, t.0, t.1);
                    // Print the polynomial for diagnosis
                    print!("    Q =");
                    for (q, t_c) in q_poly.iter().enumerate() {
                        for (k, &c) in t_c.iter().enumerate() {
                            if c.abs() > 1e-12 {
                                print!(" + {:.3}·z^{}t^{}", c, q, k);
                            }
                        }
                    }
                    println!();
                    break;
                }
            }

            if pass {
                p_pass += 1;
                if min_abs < p_min {
                    p_min = min_abs;
                }
            }
            p_tests += 1;
        }

        total_tests += p_tests;
        total_pass += p_pass;
        if p_min < min_abs_seen {
            min_abs_seen = p_min;
        }

        println!(
            "p={}: {}/{} passed  min|M_p[Q]|={:.3e}",
            p, p_pass, p_tests, p_min
        );
    }

    println!(
        "\nTOTAL: {}/{} passed  (global min|M_p[Q]|={:.3e})",
        total_pass, total_tests, min_abs_seen
    );

    if total_pass == total_tests {
        println!();
        println!("*** ALL PASSED ***");
        println!("M_p appears to preserve bivariate stability on random stable polynomials.");
        println!("Conjecture: M_p is a stability-preserving operator for all p ≥ 1.");
    }

    // ── Composition test ──────────────────────────────────────────────────────
    println!();
    println!("=== Composition test: M_{{p2}} o M_{{p1}} ===");
    println!("Check that composing two M_p operators preserves stability.");

    let mut comp_total = 0u64;
    let mut comp_pass = 0u64;

    for p1 in 1..=6usize {
        for p2 in 1..=6usize {
            let mut rng = Rng::new(p1 as u64 * 7 + p2 as u64 * 31 + 999983);
            let mut pass_count = 0u32;

            for poly_idx in 0..200 {
                let deg = 1 + poly_idx % max_deg;
                let q_poly = random_stable_poly(&mut rng, deg);
                let after_p1 = apply_mp_float(&q_poly, p1);
                let after_p2 = apply_mp_float(&after_p1, p2);

                let mut ok = true;
                for _ in 0..20 {
                    let zr = rng.next_range(-5.0, 5.0);
                    let zi = rng.next_range(0.01, 5.0);
                    let tr = rng.next_range(-5.0, 5.0);
                    let ti = rng.next_range(0.01, 5.0);
                    if cabs(eval_biv_complex(&after_p2, (zr, zi), (tr, ti))) < 1e-10 {
                        ok = false;
                        println!("  FAIL M_{} ∘ M_{}: poly#{} deg={}", p2, p1, poly_idx, deg);
                        break;
                    }
                }

                if ok {
                    pass_count += 1;
                }
                comp_total += 1;
            }
            comp_pass += pass_count as u64;
        }
    }

    println!("Compositions: {}/{} passed", comp_pass, comp_total);
    if comp_pass == comp_total {
        println!("All compositions M_{{p2}} ∘ M_{{p1}} appear stability-preserving.");
    }
}
