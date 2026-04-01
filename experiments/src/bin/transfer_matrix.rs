/// Compute the transfer matrix for the insertion step:
/// L^{(p)}_{n,S}(t) = Σ_q M_{p,q}(t) · L^{(q)}_{n-1,S'_p}(t)
///
/// where M_{p,q}(t) encodes the ε₁ and ε₂ contributions.
///
/// For each n and target descent set S, compute the transfer matrix
/// and check if it has staircase structure.
///
/// Also compute: for each (n, S, insertion position p), what fraction
/// of permutations have ε₁=1 vs ε₁=0, and can the ε₁ contribution
/// be separated cleanly?
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly, is_real_rooted};
use std::collections::BTreeMap;

fn build_poly(perms: &[&Vec<u8>]) -> Vec<i64> {
    if perms.is_empty() {
        return vec![0];
    }
    let max_s = perms
        .iter()
        .map(|s| compute(s, Stat::Swaps))
        .max()
        .unwrap_or(0);
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
    for (i, &v) in a.iter().enumerate() {
        c[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        c[i] += v;
    }
    c
}

fn poly_mul_t(a: &[i64]) -> Vec<i64> {
    let mut c = vec![0i64; a.len() + 1];
    for (i, &v) in a.iter().enumerate() {
        c[i + 1] = v;
    }
    c
}

fn descent_set_str(s: u64, n: u8) -> String {
    let positions: Vec<String> = (0..n - 1)
        .filter(|&i| s & (1 << i) != 0)
        .map(|i| (i + 1).to_string())
        .collect();
    if positions.is_empty() {
        "∅".to_string()
    } else {
        format!("{{{}}}", positions.join(","))
    }
}

/// For target descent set S in S_n, compute the source descent set S'
/// when inserting n at position p (1-indexed).
///
/// Returns None if the insertion at position p is incompatible with S.
fn source_descent_set(target_s: u64, n: u8, p: usize) -> Option<u64> {
    // Check compatibility:
    // p < n: position p must be a descent in σ (p ∈ S, since σ(p)=n > σ(p+1))
    // p > 1: position p-1 must be an ascent in σ (p-1 ∉ S, since σ(p-1) < n)

    if p < n as usize {
        // Position p (0-indexed bit: p-1) must be in S
        if target_s & (1 << (p - 1)) == 0 {
            return None;
        }
    }
    if p > 1 {
        // Position p-1 (0-indexed bit: p-2) must NOT be in S
        if target_s & (1 << (p - 2)) != 0 {
            return None;
        }
    }

    // Compute source descent set S' ⊆ [n-2]:
    // Positions 1..p-2 in π: same as S
    // Position p-1 in π: was position p-1 in σ, which is forced ascent. But in π,
    //   position p-1 now compares π(p-1) with π(p) = σ(p+1).
    //   This comparison depends on whether p-1 was originally an ascent/descent.
    //   Actually: π(p-1) = σ(p-1), π(p) = σ(p+1).
    //   In σ: position p is descent (σ(p)=n > σ(p+1)), position p-1 is ascent (σ(p-1) < n).
    //   But σ(p-1) vs σ(p+1)? This depends on whether p in S' (the shifted version).
    //   Position p-1 in π corresponds to: position p+1 in σ is shifted to position p in π.
    //   So Des(π) at position p-1 is: π(p-1) > π(p) iff σ(p-1) > σ(p+1).
    //   This is σ's descent at position p-1... no wait.
    //   Let me be careful. In σ: σ(1), ..., σ(p-1), n, σ(p+1), ..., σ(n).
    //   In π: π(1)=σ(1), ..., π(p-1)=σ(p-1), π(p)=σ(p+1), ..., π(n-1)=σ(n).
    //   So Des(π) ∩ [p-2] = Des(σ) ∩ [p-2] (same values, same positions).
    //   Position p-1 in π: π(p-1) vs π(p) = σ(p-1) vs σ(p+1).
    //   Positions ≥ p in π: position j in π = position j+1 in σ, so
    //     Des(π) at position j ≥ p: π(j) > π(j+1) iff σ(j+1) > σ(j+2),
    //     i.e., j ∈ Des(π) iff j+1 ∈ Des(σ) = S.

    // So: Des(π) = (S ∩ {1,...,p-2}) ∪ X_{p-1} ∪ {j-1 : j ∈ S, j ≥ p+1}
    // where X_{p-1} = {p-1} if σ(p-1) > σ(p+1) — but this depends on the permutation!

    // For p = 1: π(j) = σ(j+1), so Des(π) = {j-1 : j ∈ S, j ≥ 2} = {j-1 : j ∈ S ∩ {2,...,n-1}}
    if p == 1 {
        let mut s_prime = 0u64;
        for j in 1..n as u64 - 1 {
            if target_s & (1 << j) != 0 {
                s_prime |= 1 << (j - 1);
            }
        }
        return Some(s_prime);
    }

    // For p = n: π = σ without last entry, Des(π) = S ∩ [n-2]
    if p == n as usize {
        let mask = (1u64 << (n - 2)) - 1;
        return Some(target_s & mask);
    }

    // For 1 < p < n: position p-1 in π is ambiguous (depends on specific perm).
    // The source descent set at position p-1 is NOT determined by S and p alone.
    // This means the recursion doesn't factor through a single S'.
    // We need to split by whether p-1 ∈ Des(π) or not.
    // Return None to indicate we need both cases.
    None
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    println!("═══ Transfer matrix structure for insertion step ═══\n");

    // For each (n, S), compute the actual transfer matrix numerically.
    // Group permutations by (target descent set, pos of max) and
    // trace back to (source descent set, pos of max in π).

    for n in 4..=max_n.min(9) {
        let all_n = all_permutations(n);
        let all_nm1 = all_permutations(n - 1);

        // Group target perms by (descent set, pos of max)
        let mut target_groups: BTreeMap<(u64, usize), Vec<&Vec<u8>>> = BTreeMap::new();
        for s in &all_n {
            let ds = descent_set_bitmask(s);
            let pos_max = s.iter().position(|&v| v == n).unwrap() + 1;
            target_groups.entry((ds, pos_max)).or_default().push(s);
        }

        // Group source perms by (descent set, pos of max)
        let mut source_groups: BTreeMap<(u64, usize), Vec<&Vec<u8>>> = BTreeMap::new();
        for s in &all_nm1 {
            let ds = descent_set_bitmask(s);
            let pos_max = s.iter().position(|&v| v - 1 == 0).map(|_| 0).unwrap_or(0);
            // pos of max value n-1
            let pos_nm1 = s.iter().position(|&v| v == n - 1).unwrap() + 1;
            source_groups.entry((ds, pos_nm1)).or_default().push(s);
        }

        // For the alternating descent set, compute the transfer matrix
        // between position-refined polynomials.
        let alt_des = {
            let mut d = 0u64;
            for i in (1..n as u64).step_by(2) {
                if i < n as u64 {
                    d |= 1 << i; // descent at even 1-indexed positions (0-indexed bit: i)
                }
            }
            // Wait, alternating is up-down: σ(1)<σ(2)>σ(3)<σ(4)>...
            // Des at position 2, 4, 6, ... (1-indexed)
            // Bit 1, 3, 5, ... (0-indexed)
            let mut d2 = 0u64;
            for i in (2..=n as u64 - 1).step_by(2) {
                d2 |= 1 << (i - 1);
            }
            d2
        };

        // Compute transfer matrix for alternating case
        println!("n={}, alternating Des={}:", n, descent_set_str(alt_des, n));

        // Valid positions for max n in alternating perms: even positions (peaks)
        let valid_pos: Vec<usize> = (1..=n as usize)
            .filter(|&p| {
                // p < n: bit p-1 must be set in alt_des
                // p > 1: bit p-2 must NOT be set in alt_des
                let ok_descent = p >= n as usize || (alt_des & (1 << (p - 1))) != 0;
                let ok_ascent = p <= 1 || (alt_des & (1 << (p - 2))) == 0;
                ok_descent && ok_ascent
            })
            .collect();

        println!("  Valid positions for max: {:?}", valid_pos);

        // For each valid position p, and each source position q of value n-1,
        // compute the "matrix entry" M[p][q] as a polynomial in t.
        // M[p][q](t) = Σ_{π: pos(n-1)=q, inserting n at p gives target} t^{ε₁+ε₂}

        // Actually, let's compute empirically by tracing insertions.
        let mut matrix_entries: BTreeMap<(usize, usize, u64), Vec<i64>> = BTreeMap::new();

        for pi in &all_nm1 {
            let source_ds = descent_set_bitmask(pi);
            let q = pi.iter().position(|&v| v == n - 1).unwrap() + 1; // pos of n-1

            for &p in &valid_pos {
                // Insert n at position p
                let mut sigma = Vec::with_capacity(n as usize);
                for i in 0..p - 1 {
                    sigma.push(pi[i]);
                }
                sigma.push(n);
                for i in p - 1..pi.len() {
                    sigma.push(pi[i]);
                }

                let target_ds = descent_set_bitmask(&sigma);
                if target_ds != alt_des {
                    continue;
                }

                let sw_pi = compute(pi, Stat::Swaps);
                let sw_sigma = compute(&sigma, Stat::Swaps);
                let eps = sw_sigma - sw_pi;

                let entry = matrix_entries
                    .entry((p, q, source_ds))
                    .or_insert(vec![0i64; 3]);
                if eps < entry.len() {
                    entry[eps] += 1;
                } else {
                    entry.resize(eps + 1, 0);
                    entry[eps] += 1;
                }
            }
        }

        // Print the transfer matrix grouped by target position
        // For each source descent set that contributes:
        let source_des_sets: Vec<u64> = matrix_entries
            .keys()
            .map(|&(_, _, ds)| ds)
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect();

        for &sds in &source_des_sets {
            println!("  Source Des={}:", descent_set_str(sds, n - 1));

            // Source positions (q values that appear)
            let source_pos: Vec<usize> = matrix_entries
                .keys()
                .filter(|&&(_, _, ds)| ds == sds)
                .map(|&(_, q, _)| q)
                .collect::<std::collections::BTreeSet<_>>()
                .into_iter()
                .collect();

            print!("    q\\p  ");
            for &p in &valid_pos {
                print!("{:>8}", p);
            }
            println!();

            for &q in &source_pos {
                print!("    q={:<2}: ", q);
                for &p in &valid_pos {
                    match matrix_entries.get(&(p, q, sds)) {
                        Some(entry) => {
                            // Format as polynomial
                            let parts: Vec<String> = entry
                                .iter()
                                .enumerate()
                                .filter(|(_, &c)| c > 0)
                                .map(|(d, &c)| {
                                    if d == 0 {
                                        format!("{}", c)
                                    } else if c == 1 {
                                        format!("t^{}", d)
                                    } else {
                                        format!("{}t^{}", c, d)
                                    }
                                })
                                .collect();
                            let s = if parts.is_empty() {
                                "0".to_string()
                            } else {
                                parts.join("+")
                            };
                            print!("{:>8}", s);
                        }
                        None => print!("{:>8}", "0"),
                    }
                }
                println!();
            }
        }

        println!();
    }
}
