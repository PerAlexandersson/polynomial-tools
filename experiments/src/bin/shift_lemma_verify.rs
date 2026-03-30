/// Verify shift lemma hypotheses for the (t-1) correction in the transfer matrix.
///
/// The transfer matrix decomposes as:
///   M(t) = M_staircase(t) + (t-1) · M_correction(t)
///
/// Applied to the position-refined input vector v = (L^{(q)}_{n-1,S'}(t)):
///   L^{(p)}_{n,S}(t) = f_p(t) + (t-1) · h_p(t)
///
/// where f_p = M_staircase · v and h_p = M_correction · v.
///
/// The shift lemma requires: h_p ≪ f_p and f_p(0) ≥ h_p(0).
///
/// We verify this for ALL descent sets S ⊆ [n-1] and all valid positions p.
///
/// Additionally, we check whether the full output vector
/// (f_1 + (t-1)h_1, ..., f_n + (t-1)h_n) forms an interlacing sequence.
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{
    check_weak_interlacing, format_poly, is_real_rooted,
};
use std::collections::BTreeMap;

fn build_poly(perms: &[&Vec<u8>]) -> Vec<i64> {
    if perms.is_empty() {
        return vec![0];
    }
    let max_s = perms.iter().map(|s| compute(s, Stat::Swaps)).max().unwrap_or(0);
    let mut coeffs = vec![0i64; max_s + 1];
    for s in perms {
        coeffs[compute(s, Stat::Swaps)] += 1;
    }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
        coeffs.pop();
    }
    coeffs
}

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut c = vec![0i64; len];
    for (i, &v) in a.iter().enumerate() { c[i] += v; }
    for (i, &v) in b.iter().enumerate() { c[i] += v; }
    c
}

/// Multiply polynomial by t.
fn poly_mul_t(a: &[i64]) -> Vec<i64> {
    if a.is_empty() || (a.len() == 1 && a[0] == 0) {
        return vec![0];
    }
    let mut c = vec![0i64; a.len() + 1];
    for (i, &v) in a.iter().enumerate() { c[i + 1] = v; }
    c
}

/// Multiply polynomial by (t-1): result = t*a - a.
fn poly_mul_t_minus_1(a: &[i64]) -> Vec<i64> {
    let ta = poly_mul_t(a);
    let len = ta.len().max(a.len());
    let mut c = vec![0i64; len];
    for (i, &v) in ta.iter().enumerate() { c[i] += v; }
    for (i, &v) in a.iter().enumerate() { c[i] -= v; }
    while c.len() > 1 && *c.last().unwrap() == 0 { c.pop(); }
    c
}

fn weakly_interlace(f: &[i64], g: &[i64]) -> bool {
    if f.len() <= 1 || g.len() <= 1 {
        return true;
    }
    let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
    check_weak_interlacing(small, large) == Some(true)
}

fn descent_set_str(s: u64, n: u8) -> String {
    let positions: Vec<String> = (0..n - 1)
        .filter(|&i| s & (1 << i) != 0)
        .map(|i| (i + 1).to_string())
        .collect();
    if positions.is_empty() { "∅".to_string() }
    else { format!("{{{}}}", positions.join(",")) }
}

/// Compute ε₁ for inserting n at position p in π.
fn epsilon1(pi: &[u8], p: usize, n: u8) -> usize {
    if p <= 1 || p >= n as usize { return 0; }
    if pi[p - 2] + 1 == pi[p - 1] { 1 } else { 0 }
}

/// Compute ε₂ for inserting n at position p, given pos of n-1 in π is q.
fn epsilon2(q: usize, p: usize) -> usize {
    if q + 2 <= p { 1 } else { 0 }
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    println!("═══ Shift lemma verification for transfer matrix ═══\n");

    let mut total_checks = 0u64;
    let mut shift_ok = 0u64;
    let mut shift_fail = 0u64;
    let mut interlace_ok = 0u64;
    let mut interlace_fail = 0u64;

    for n in 4..=max_n {
        let all_nm1 = all_permutations(n - 1);
        let all_n = all_permutations(n);

        // Group target by descent set
        let mut target_by_des: BTreeMap<u64, Vec<&Vec<u8>>> = BTreeMap::new();
        for s in &all_n {
            target_by_des.entry(descent_set_bitmask(s)).or_default().push(s);
        }

        let mut n_shift_ok = 0u32;
        let mut n_shift_fail = 0u32;
        let mut n_interlace_ok = 0u32;
        let mut n_interlace_fail = 0u32;
        let mut n_total = 0u32;
        let mut fail_examples: Vec<String> = Vec::new();

        for (&target_ds, _target_perms) in &target_by_des {
            // Determine valid insertion positions
            let valid_pos: Vec<usize> = (1..=n as usize)
                .filter(|&p| {
                    let ok_descent = p >= n as usize || (target_ds & (1 << (p - 1))) != 0;
                    let ok_ascent = p <= 1 || (target_ds & (1 << (p - 2))) == 0;
                    ok_descent && ok_ascent
                })
                .collect();

            if valid_pos.len() < 2 {
                continue; // trivial case
            }

            n_total += 1;

            // For each valid position p, compute f_p (staircase) and h_p (correction).
            // f_p(t) = Σ_q (Σ_{π: pos(n-1)=q, compatible} t^{ε₂} · t^{swaps(π)})
            // h_p(t) = Σ_q (Σ_{π: pos(n-1)=q, compatible, ε₁=1} t^{ε₂} · t^{swaps(π)})
            //
            // Then L^{(p)} = f_p + (t-1) · h_p.

            let mut f_polys: Vec<(usize, Vec<i64>)> = Vec::new();
            let mut h_polys: Vec<(usize, Vec<i64>)> = Vec::new();
            let mut l_polys: Vec<(usize, Vec<i64>)> = Vec::new();

            for &p in &valid_pos {
                let mut f_p = vec![0i64; 1];
                let mut h_p = vec![0i64; 1];

                for pi in &all_nm1 {
                    // Insert n at position p
                    let mut sigma = Vec::with_capacity(n as usize);
                    for i in 0..p - 1 { sigma.push(pi[i]); }
                    sigma.push(n);
                    for i in p - 1..pi.len() { sigma.push(pi[i]); }

                    if descent_set_bitmask(&sigma) != target_ds {
                        continue;
                    }

                    let sw = compute(pi, Stat::Swaps);
                    let q = pi.iter().position(|&v| v == n - 1).unwrap() + 1;
                    let e1 = epsilon1(pi, p, n);
                    let e2 = epsilon2(q, p);

                    // f_p gets t^{ε₂} · t^{sw} for ALL compatible π
                    let f_deg = e2 + sw;
                    if f_deg >= f_p.len() { f_p.resize(f_deg + 1, 0); }
                    f_p[f_deg] += 1;

                    // h_p gets t^{ε₂} · t^{sw} only for π with ε₁=1
                    if e1 == 1 {
                        let h_deg = e2 + sw;
                        if h_deg >= h_p.len() { h_p.resize(h_deg + 1, 0); }
                        h_p[h_deg] += 1;
                    }
                }

                // Trim
                while f_p.len() > 1 && *f_p.last().unwrap() == 0 { f_p.pop(); }
                while h_p.len() > 1 && *h_p.last().unwrap() == 0 { h_p.pop(); }

                // L^{(p)} = f_p + (t-1) · h_p
                let correction = poly_mul_t_minus_1(&h_p);
                let l_p = poly_add(&f_p, &correction);

                f_polys.push((p, f_p));
                h_polys.push((p, h_p));
                l_polys.push((p, l_p));
            }

            // Verify shift lemma for each position p:
            // Need: h_p ≪ f_p (weak interlacing) and f_p(0) ≥ h_p(0)
            let mut all_shift_ok = true;
            for i in 0..f_polys.len() {
                let (p, ref f) = f_polys[i];
                let (_, ref h) = h_polys[i];
                let (_, ref l) = l_polys[i];

                // Check f is real-rooted
                let f_rr = f.len() <= 2 || is_real_rooted(f);
                // Check h is real-rooted
                let h_rr = h.len() <= 2 || is_real_rooted(h);
                // Check l is real-rooted
                let l_rr = l.len() <= 2 || is_real_rooted(l);

                // Check h ≪ f
                let h_int_f = weakly_interlace(h, f);

                // Check f(0) ≥ h(0)
                let f0 = if f.is_empty() { 0 } else { f[0] };
                let h0 = if h.is_empty() { 0 } else { h[0] };
                let f0_ge_h0 = f0 >= h0;

                total_checks += 1;

                if !h_int_f || !f0_ge_h0 || !l_rr {
                    all_shift_ok = false;
                    if fail_examples.len() < 5 {
                        fail_examples.push(format!(
                            "n={} Des={} p={}: h≪f:{} f(0)≥h(0):{} f_rr:{} h_rr:{} L_rr:{}\n      f={}\n      h={}\n      L={}",
                            n, descent_set_str(target_ds, n), p,
                            h_int_f, f0_ge_h0, f_rr, h_rr, l_rr,
                            format_poly(f), format_poly(h), format_poly(l),
                        ));
                    }
                } else {
                    shift_ok += 1;
                }
            }

            if all_shift_ok {
                n_shift_ok += 1;
            } else {
                n_shift_fail += 1;
                shift_fail += f_polys.len() as u64 - shift_ok;
            }

            // Check if the output L^{(p)} form an interlacing sequence
            // (pairwise weakly interlacing)
            let mut pairwise_ok = true;
            for i in 0..l_polys.len() {
                for j in i + 1..l_polys.len() {
                    if !weakly_interlace(&l_polys[i].1, &l_polys[j].1) {
                        pairwise_ok = false;
                    }
                }
            }

            if pairwise_ok {
                n_interlace_ok += 1;
                interlace_ok += 1;
            } else {
                n_interlace_fail += 1;
                interlace_fail += 1;
            }
        }

        println!(
            "n={}: shift lemma: {}/{} ok | output interlacing: {}/{} ok",
            n, n_shift_ok, n_total, n_interlace_ok, n_total,
        );
        for ex in &fail_examples {
            println!("  {}", ex);
        }
    }

    println!("\n═══ Summary ═══");
    println!(
        "Shift lemma: {} checks ok, {} fail",
        shift_ok, shift_fail,
    );
    println!(
        "Output interlacing: {} ok, {} fail",
        interlace_ok, interlace_fail,
    );
}
