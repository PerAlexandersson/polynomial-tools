//! Test whether reversed AA propagates through row-deletion recursion for
//! restricted hit polynomials.
//!
//! For each 312-avoiding board λ with n ≤ max_n and each sub-partition μ ⊆ λ:
//! 1. Compute A_k at the (λ, μ) level by brute force.
//! 2. Build the λ^+ board (prepend row of length m+1) and all compatible μ^+.
//! 3. Compute A_{k'}^+ at the (λ^+, μ^+) level by brute force.
//! 4. Test reversed AA at the λ^+ level.
//! 5. Verify the recursion: A_{k'}^+ = Σ_{j<k'} A_j + t·Σ_{j≥k'} A_j (when μ^+_1=0)
//!    and the general form with hit weights.
//! 6. Examine consecutive differences A_{k'}^+ - A_{k'+1}^+ and check shift-lemma structure.
//!
//! Usage: cargo run --release --bin hit_revaa_prop -- 5

use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly, is_real_rooted};
use std::collections::BTreeSet;

// ── Polynomial helpers ──

fn pt(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
    v
}
fn pz(p: &[i64]) -> bool {
    p.iter().all(|&c| c == 0)
}
fn pa(a: &[i64], b: &[i64]) -> Vec<i64> {
    let l = a.len().max(b.len());
    let mut r = vec![0i64; l];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] += v;
    }
    pt(&r)
}
fn ps(a: &[i64], b: &[i64]) -> Vec<i64> {
    let l = a.len().max(b.len());
    let mut r = vec![0i64; l];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] -= v;
    }
    pt(&r)
}
fn pmt(p: &[i64]) -> Vec<i64> {
    let mut r = vec![0i64; p.len() + 1];
    for (i, &v) in p.iter().enumerate() {
        r[i + 1] = v;
    }
    pt(&r)
}
fn pdeg(p: &[i64]) -> Option<usize> {
    let v = pt(p);
    if pz(&v) {
        None
    } else {
        Some(v.len() - 1)
    }
}
fn pscale(p: &[i64], c: i64) -> Vec<i64> {
    pt(&p.iter().map(|&x| x * c).collect::<Vec<_>>())
}

fn interlaces(f: &[i64], g: &[i64]) -> bool {
    let f = pt(f);
    let g = pt(g);
    if pz(&f) {
        return true;
    }
    if pz(&g) {
        return false;
    }
    match check_weak_interlacing(&f, &g) {
        Some(true) => true,
        Some(false) => false,
        None => {
            let df = pdeg(&f);
            let dg = pdeg(&g);
            match (df, dg) {
                (Some(df), Some(dg)) if df == dg => {
                    let tf = pmt(&f);
                    match check_weak_interlacing(&g, &tf) {
                        Some(b) => b,
                        None => false,
                    }
                }
                _ => false,
            }
        }
    }
}

fn show(name: &str, c: [u64; 2]) {
    if c[0] == 0 {
        println!("  {}: (no data)", name);
    } else if c[1] == 0 {
        println!("  {}: {}/{} ALL PASS", name, c[0], c[0]);
    } else {
        println!("  {}: {}/{} pass ({} FAIL)", name, c[0] - c[1], c[0], c[1]);
    }
}

// ── Combinatorial infrastructure ──

fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> {
    let n = perm.len();
    let mut vis: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut q: BTreeSet<Vec<u8>> = BTreeSet::new();
    q.insert(perm.to_vec());
    while let Some(cur) = q.pop_last() {
        for i in 0..n {
            for j in i + 1..n {
                if cur[i] > cur[j] {
                    let mut c = cur.clone();
                    c.swap(i, j);
                    if !vis.contains(&c) {
                        q.insert(c);
                    }
                }
            }
        }
        vis.insert(cur);
    }
    vis.into_iter().collect()
}

fn board_to_perm(b: &[u8]) -> Vec<u8> {
    let n = b.len();
    let mut p = vec![0u8; n];
    let mut u = vec![false; n + 1];
    for i in 0..n {
        for c in (1..=(b[i] as usize).min(n)).rev() {
            if !u[c] {
                p[i] = c as u8;
                u[c] = true;
                break;
            }
        }
    }
    p
}

fn is_312_avoiding(perm: &[u8]) -> bool {
    let n = perm.len();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                if perm[k] < perm[i] && perm[i] < perm[j] {
                    return false;
                }
            }
        }
    }
    true
}

fn gen_boards(n: usize) -> Vec<Vec<u8>> {
    let mut r = vec![];
    let mut c = vec![];
    gb(n, n, 0, &mut c, &mut r);
    r
}

fn gb(n: usize, mx: usize, d: usize, c: &mut Vec<u8>, r: &mut Vec<Vec<u8>>) {
    if d == n {
        r.push(c.clone());
        return;
    }
    for v in (d + 1).max(if d > 0 { c[d - 1] as usize } else { 1 })..=mx {
        c.push(v as u8);
        gb(n, mx, d + 1, c, r);
        c.pop();
    }
}

/// Generate all sub-partitions mu with mu_i <= lambda_i and mu weakly decreasing.
fn sub_partitions(lambda: &[u8]) -> Vec<Vec<u8>> {
    let mut result = Vec::new();
    let mut mu = vec![0u8; lambda.len()];
    fn gen(lambda: &[u8], mu: &mut Vec<u8>, pos: usize, max_val: u8, result: &mut Vec<Vec<u8>>) {
        if pos == lambda.len() {
            result.push(mu.clone());
            return;
        }
        let upper = lambda[pos].min(max_val);
        for v in 0..=upper {
            mu[pos] = v;
            gen(lambda, mu, pos + 1, v, result);
        }
    }
    gen(lambda, &mut mu, 0, lambda[0], &mut result);
    result
}

/// Count hits: #{i : sigma_i > mu_i}
fn count_hits(sigma: &[u8], mu: &[u8]) -> usize {
    (0..sigma.len())
        .filter(|&i| {
            let mu_i = if i < mu.len() { mu[i] as usize } else { 0 };
            (sigma[i] as usize) > mu_i
        })
        .count()
}

/// Build λ^+ = (n+1, λ_1+1, ..., λ_n+1) from λ.
/// For a 312-avoiding Ferrers board λ = [λ_1, ..., λ_n] with λ_i >= i and weakly increasing,
/// λ^+ = [n+1, λ_1+1, ..., λ_n+1].
fn lambda_plus(lambda: &[u8]) -> Vec<u8> {
    let n = lambda.len();
    let mut lp = Vec::with_capacity(n + 1);
    lp.push((n + 1) as u8);
    for &v in lambda {
        lp.push(v + 1);
    }
    lp
}

/// Generate all sub-partitions mu^+ of lambda^+ that are compatible.
/// mu^+ is a weakly decreasing sequence with mu^+_i <= lambda^+_i.
fn sub_partitions_plus(lambda_plus: &[u8]) -> Vec<Vec<u8>> {
    sub_partitions(lambda_plus)
}

struct FailInfo {
    board: Vec<u8>,
    mu: Vec<u8>,
    kp1: usize,
    kp2: usize,
    f: Vec<i64>,
    g: Vec<i64>,
}

fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(5);
    println!("================================================================");
    println!("  Reversed AA propagation for hit polynomials");
    println!("  312-avoiding Ferrers boards, n <= {}", max_n);
    println!("================================================================\n");

    // Counters: [total, fails]
    let mut c_revaa_base = [0u64; 2]; // revAA at base level (λ, μ)
    let mut c_revaa_plus = [0u64; 2]; // revAA at λ^+ level
    let mut c_revaa_plus_wg = [0u64; 2]; // revAA at λ^+ within-group
    let mut c_revaa_plus_xg = [0u64; 2]; // revAA at λ^+ cross-group

    // Recursion verification
    let mut c_recursion = [0u64; 2]; // recursion formula matches

    // Consecutive difference structure
    let mut c_diff_tm1 = [0u64; 2]; // A_{k'}^+ - A_{k'+1}^+ is (t-1)*h
    let mut c_shift_lemma = [0u64; 2]; // h << A_{k'+1}^+ (shift lemma hypothesis)
    let mut c_shift_result = [0u64; 2]; // A_{k'+1}^+ << A_{k'+1}^+ + (t-1)*h = A_{k'}^+

    // Failure details
    let mut revaa_plus_fails: Vec<FailInfo> = Vec::new();
    let mut shift_fails: Vec<FailInfo> = Vec::new();

    let mut total_pairs = 0u64;
    let mut total_plus_pairs = 0u64;

    // Track some example consecutive differences
    let mut example_diffs: Vec<(Vec<u8>, Vec<u8>, usize, usize, Vec<i64>, Vec<i64>, Vec<i64>)> =
        Vec::new();

    for n in 2..=max_n {
        let boards = gen_boards(n);
        let mut n_boards = 0u64;
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            n_boards += 1;

            let m = board[0] as usize;

            // ── BASE LEVEL: verify revAA at (λ, μ) ──
            let ideal = bruhat_lower_ideal(&perm);
            let subs = sub_partitions(board);

            for mu in &subs {
                total_pairs += 1;

                let mut a_k: Vec<Vec<i64>> = vec![vec![0i64]; m + 1];
                for sigma in &ideal {
                    let k = sigma[0] as usize;
                    if k == 0 || k > m {
                        continue;
                    }
                    let hits = count_hits(sigma, mu);
                    while a_k[k].len() <= hits {
                        a_k[k].push(0);
                    }
                    a_k[k][hits] += 1;
                }
                for k in 1..=m {
                    a_k[k] = pt(&a_k[k]);
                }

                let active: Vec<usize> = (1..=m).filter(|&k| !pz(&a_k[k])).collect();

                // revAA at base level
                for &k in &active {
                    for &kp in &active {
                        if kp >= k {
                            continue;
                        }
                        c_revaa_base[0] += 1;
                        if !interlaces(&a_k[k], &a_k[kp]) {
                            c_revaa_base[1] += 1;
                        }
                    }
                }
            }

            // ── λ^+ LEVEL ──
            let lp = lambda_plus(board);
            let perm_plus = board_to_perm(&lp);
            if !is_312_avoiding(&perm_plus) {
                continue;
            }
            let mp = lp[0] as usize; // = m + 1 for standard boards, but let's be safe

            let ideal_plus = bruhat_lower_ideal(&perm_plus);
            let subs_plus = sub_partitions(&lp);

            for mu_plus in &subs_plus {
                total_plus_pairs += 1;

                // Compute A_{k'}^+ by brute force at (λ^+, μ^+)
                let mut a_kp: Vec<Vec<i64>> = vec![vec![0i64]; mp + 1];
                // Also compute D_{k'}^+ and U_{k'}^+ (descent/ascent at position 1)
                let mut d_kp: Vec<Vec<i64>> = vec![vec![0i64]; mp + 1];
                let mut u_kp: Vec<Vec<i64>> = vec![vec![0i64]; mp + 1];

                let np = lp.len(); // n + 1
                for sigma in &ideal_plus {
                    let k = sigma[0] as usize;
                    if k == 0 || k > mp {
                        continue;
                    }
                    let hits = count_hits(sigma, mu_plus);
                    while a_kp[k].len() <= hits {
                        a_kp[k].push(0);
                    }
                    a_kp[k][hits] += 1;

                    let poly = if np >= 2 && sigma[0] > sigma[1] {
                        &mut d_kp[k]
                    } else {
                        &mut u_kp[k]
                    };
                    while poly.len() <= hits {
                        poly.push(0);
                    }
                    poly[hits] += 1;
                }
                for k in 1..=mp {
                    a_kp[k] = pt(&a_kp[k]);
                    d_kp[k] = pt(&d_kp[k]);
                    u_kp[k] = pt(&u_kp[k]);
                }

                let active_plus: Vec<usize> = (1..=mp).filter(|&k| !pz(&a_kp[k])).collect();

                // ── Test reversed AA at λ^+ ──
                for &k in &active_plus {
                    for &kp in &active_plus {
                        if kp >= k {
                            continue;
                        }
                        c_revaa_plus[0] += 1;
                        if !interlaces(&a_kp[k], &a_kp[kp]) {
                            c_revaa_plus[1] += 1;
                            if revaa_plus_fails.len() < 10 {
                                revaa_plus_fails.push(FailInfo {
                                    board: lp.clone(),
                                    mu: mu_plus.clone(),
                                    kp1: k,
                                    kp2: kp,
                                    f: a_kp[k].clone(),
                                    g: a_kp[kp].clone(),
                                });
                            }
                        }

                        // Determine groups: k' > l' uses col l', k' < l' uses col l'-1
                        // Here l' is the last entry... but we're not refining by l' here.
                        // For this test: within-group = both in upper or both in lower
                        // relative to the "new row boundary".
                        // Actually, for the standard recursion, the group depends on the
                        // relationship to the last entry. Since we're testing ALL pairs,
                        // let's just count all + separate by position relative to some threshold.
                    }
                }

                // ── Verify recursion structure ──
                // At λ^+ level, removing the first row: the remaining σ[1..] is in
                // the Bruhat ideal of λ (with shifted values). The hit contribution
                // from position 1 is t^{[k' > μ^+_1]}.
                //
                // So A_{k'}^+(λ^+, μ^+) = t^{[k'>μ^+_1]} · Σ_j (contribution from
                // removing k' from positions 2..n+1).
                //
                // More precisely: for σ with σ_1 = k', the remaining entries σ[1..]
                // form a permutation of {1,...,n+1}\{k'} in Bruhat order.
                // After standardizing, this relates to the λ-level Bruhat ideal.
                //
                // The hit weight: position i>0 contributes t^{[σ_i > μ^+_i]}.
                //
                // Let's verify by decomposing A_{k'}^+ into the formula:
                // For mu^+_1 = 0: position 1 always hits (k' > 0), so extra factor t.
                //   Remaining: σ[1..] in {1..n+1}\{k'}, standardized to {1..n}.
                //   Remove column k' from the board: λ' with (λ')_i = λ_i - [k' ≤ λ_i+1]...
                //   Actually the recursion is more subtle. Let's just verify the formula
                //   A_{k'}^+ = Σ_{j<k'} A_j^{(mu_rest)} + t·Σ_{j≥k'} A_j^{(mu_rest)}
                //   only when mu^+_1 = 0 (the ascent/peak case).

                // For now: test the structural relationship via consecutive differences.

                // ── Consecutive differences ──
                for i in 0..active_plus.len() {
                    if i + 1 >= active_plus.len() {
                        continue;
                    }
                    let k1 = active_plus[i]; // larger k'
                    let k2 = active_plus[i + 1]; // smaller k'
                                                 // Wait: active_plus is sorted ascending, so active_plus[i] < active_plus[i+1]
                                                 // reversed AA means larger interlaces smaller.
                                                 // consecutive: k' and k'+1
                    if k2 != k1 + 1 {
                        continue;
                    }

                    // diff = A_{k1}^+ - A_{k2}^+  (k1 < k2, so this is "smaller minus larger")
                    // For reversed AA we want A_{k2} << A_{k1} (k2 > k1),
                    // i.e., A_{k1} is on the RIGHT.
                    // The shift-lemma form: A_{k1} = A_{k2} + (t-1)*h?
                    // Actually reversed: A_{k2} << A_{k1}, and k2 > k1.
                    // Let's compute diff = A_{k2} - A_{k1} and check if = (t-1)*h.
                    let diff = ps(&a_kp[k2], &a_kp[k1]);

                    // Check if diff = (t-1)*h, i.e., diff has constant term 0 and
                    // (diff + h) / t = h, so diff = t*h - h, meaning diff[0] = -h[0],
                    // diff[i] = h[i-1] - h[i] for i >= 1.
                    // Equivalently: diff factors through (t-1).
                    // (t-1) divides f iff f(1) = 0.
                    let val_at_1: i64 = diff.iter().sum();
                    c_diff_tm1[0] += 1;
                    if val_at_1 == 0 && !pz(&diff) {
                        // diff = (t-1)*h. Recover h:
                        // diff = -h[0] + (h[0]-h[1])t + (h[1]-h[2])t^2 + ...
                        // So h[0] = -diff[0], h[1] = h[0] - diff[1], etc.
                        let mut h = Vec::new();
                        let mut running: i64 = 0;
                        for &c in &diff {
                            running -= c;
                            h.push(running);
                        }
                        // Remove trailing zeros and the last entry (which should be 0)
                        let h = pt(&h);

                        // Verify: (t-1)*h should equal diff
                        // (t-1)*h = t*h - h
                        let th = pmt(&h);
                        let check = ps(&th, &h);
                        let check = pt(&check);
                        let diff_t = pt(&diff);
                        let ok = check == diff_t;
                        if !ok {
                            c_diff_tm1[1] += 1;
                        } else {
                            // Check shift-lemma: h << A_{k1}?
                            c_shift_lemma[0] += 1;
                            if !pz(&h) && !interlaces(&h, &a_kp[k1]) {
                                c_shift_lemma[1] += 1;
                                if shift_fails.len() < 10 {
                                    shift_fails.push(FailInfo {
                                        board: lp.clone(),
                                        mu: mu_plus.clone(),
                                        kp1: k1,
                                        kp2: k2,
                                        f: h.clone(),
                                        g: a_kp[k1].clone(),
                                    });
                                }
                            }

                            // Also check: A_{k1} << A_{k1} + (t-1)*h = A_{k2}
                            // i.e., does A_{k1} << A_{k2}? (This IS the reversed AA condition)
                            c_shift_result[0] += 1;
                            if !interlaces(&a_kp[k2], &a_kp[k1]) {
                                c_shift_result[1] += 1;
                            }

                            // Save example
                            if example_diffs.len() < 15 && !pz(&h) && pdeg(&h).unwrap_or(0) >= 1 {
                                example_diffs.push((
                                    lp.clone(),
                                    mu_plus.clone(),
                                    k1,
                                    k2,
                                    a_kp[k1].clone(),
                                    a_kp[k2].clone(),
                                    h.clone(),
                                ));
                            }
                        }
                    } else if !pz(&diff) {
                        c_diff_tm1[1] += 1;
                    }
                    // If diff is zero, it trivially factors.
                }

                // ── Within-group vs cross-group at λ^+ ──
                // The recursion for A_{k'}^+ at λ^+ with last entry l':
                // For k' > l': uses column l' from λ level ("upper group")
                // For k' < l': uses column l'-1 from λ level ("lower group")
                // We need to refine by last entry to separate groups.

                // Compute A_{k',l'}^+ refined by both first and last entry
                let mut a_kp_lp: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; mp + 1]; mp + 1];
                for sigma in &ideal_plus {
                    let k = sigma[0] as usize;
                    let l = *sigma.last().unwrap() as usize;
                    if k == 0 || k > mp || l == 0 || l > mp {
                        continue;
                    }
                    let hits = count_hits(sigma, mu_plus);
                    while a_kp_lp[k][l].len() <= hits {
                        a_kp_lp[k][l].push(0);
                    }
                    a_kp_lp[k][l][hits] += 1;
                }
                for k in 1..=mp {
                    for l in 1..=mp {
                        a_kp_lp[k][l] = pt(&a_kp_lp[k][l]);
                    }
                }

                // For each last entry l', test within-group and cross-group revAA
                for lp_val in 1..=mp {
                    let mut upper: Vec<usize> = Vec::new(); // k' > lp_val
                    let mut lower: Vec<usize> = Vec::new(); // k' < lp_val
                    for k in 1..=mp {
                        if k == lp_val {
                            continue;
                        }
                        if pz(&a_kp_lp[k][lp_val]) {
                            continue;
                        }
                        if k > lp_val {
                            upper.push(k);
                        } else {
                            lower.push(k);
                        }
                    }

                    // Within-group: both upper or both lower
                    for group in [&upper, &lower] {
                        for &k in group {
                            for &kp in group {
                                if kp >= k {
                                    continue;
                                }
                                c_revaa_plus_wg[0] += 1;
                                if !interlaces(&a_kp_lp[k][lp_val], &a_kp_lp[kp][lp_val]) {
                                    c_revaa_plus_wg[1] += 1;
                                }
                            }
                        }
                    }

                    // Cross-group: one from upper, one from lower
                    for &k in &upper {
                        for &kp in &lower {
                            // k > lp_val > kp, so k > kp: reversed AA says A_k << A_kp
                            c_revaa_plus_xg[0] += 1;
                            if !interlaces(&a_kp_lp[k][lp_val], &a_kp_lp[kp][lp_val]) {
                                c_revaa_plus_xg[1] += 1;
                            }
                        }
                    }
                    // Also: lower k > upper kp? No: upper has k>l', lower has k<l', so upper always > lower.
                }
            }
        }
        println!(
            "n={}: {} valid boards, {} (lambda,mu) pairs, {} (lambda+,mu+) pairs so far",
            n, n_boards, total_pairs, total_plus_pairs
        );
    }

    println!("\n================================================================");
    println!("  BASE LEVEL (lambda, mu)");
    println!("================================================================");
    show("revAA: A_k << A_{k'} (k > k')", c_revaa_base);

    println!("\n================================================================");
    println!("  LAMBDA^+ LEVEL (lambda+, mu+)");
    println!("================================================================");
    show("revAA at lambda+: all pairs", c_revaa_plus);
    show("revAA at lambda+: within-group (same col)", c_revaa_plus_wg);
    show("revAA at lambda+: cross-group (diff col)", c_revaa_plus_xg);

    println!("\n================================================================");
    println!("  CONSECUTIVE DIFFERENCE STRUCTURE");
    println!("================================================================");
    show("A_{k+1}^+ - A_k^+ = (t-1)*h", c_diff_tm1);
    show("  h << A_k^+ (shift-lemma hyp)", c_shift_lemma);
    show("  A_{k+1} << A_k (shift-lemma result)", c_shift_result);

    // Print failure details
    if !revaa_plus_fails.is_empty() {
        println!(
            "\n  --- revAA failures at lambda+ (first {}) ---",
            revaa_plus_fails.len()
        );
        for f in &revaa_plus_fails {
            println!(
                "    lambda+={:?} mu+={:?}: A_{} << A_{} FAILS",
                f.board, f.mu, f.kp1, f.kp2
            );
            println!("      A_{} = {}", f.kp1, format_poly(&f.f));
            println!("      A_{} = {}", f.kp2, format_poly(&f.g));
        }
    }

    if !shift_fails.is_empty() {
        println!(
            "\n  --- Shift-lemma failures (h << A_k) (first {}) ---",
            shift_fails.len()
        );
        for f in &shift_fails {
            println!(
                "    lambda+={:?} mu+={:?}: h << A_{} FAILS (k'={}, k'+1={})",
                f.board, f.mu, f.kp1, f.kp1, f.kp2
            );
            println!("      h = {}", format_poly(&f.f));
            println!("      A_{} = {}", f.kp1, format_poly(&f.g));
        }
    }

    // Print example consecutive differences
    if !example_diffs.is_empty() {
        println!("\n================================================================");
        println!(
            "  EXAMPLE CONSECUTIVE DIFFERENCES (first {}):",
            example_diffs.len()
        );
        println!("================================================================");
        for (lp, mu, k1, k2, ak1, ak2, h) in &example_diffs {
            println!("  lambda+={:?} mu+={:?}, k'={}, k'+1={}", lp, mu, k1, k2);
            println!("    A_{} = {}", k1, format_poly(ak1));
            println!("    A_{} = {}", k2, format_poly(ak2));
            println!(
                "    h = {}  (A_{{{}}} - A_{{{}}} = (t-1)*h)",
                format_poly(h),
                k2,
                k1
            );
            println!("    h << A_{}? {}", k1, interlaces(h, ak1));
            println!();
        }
    }

    println!("\nTotal (lambda,mu) pairs: {}", total_pairs);
    println!("Total (lambda+,mu+) pairs: {}", total_plus_pairs);
    println!("\nDone.");
}
