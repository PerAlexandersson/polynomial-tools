/// Fast DP computation of H_n(t) using backward-reachability pruning.
///
/// Only track descent sets that can eventually lead to the alternating
/// descent set at the final step. This dramatically reduces the state space.
use combpoly::statistics::descent_set_bitmask;
use polynomial_tools::real_rootedness::{
    format_poly, gamma_coefficients, is_gamma_positive, is_real_rooted,
};
use std::collections::{HashMap, HashSet};
use std::time::Instant;

/// Compute the alternating descent bitmask for size n.
/// Alternating (up-down): descents at positions 2, 4, 6, ...
/// Bits: 1, 3, 5, ... (0-indexed)
fn alt_des(n: u8) -> u64 {
    let mut d = 0u64;
    for i in (2..=n as u64 - 1).step_by(2) {
        d |= 1 << (i - 1);
    }
    d
}

/// Given target descent set `target` in S_n, compute all valid insertion
/// positions and their corresponding source descent sets in S_{n-1}.
/// Returns Vec<(position_p, source_des)>.
fn source_descent_sets(target: u64, n: u8) -> Vec<(usize, u64)> {
    let mut result = Vec::new();

    for p in 1..=n as usize {
        // Check validity: p ∈ target (or p=n) and p-1 ∉ target
        let ok_descent = p >= n as usize || (target & (1 << (p - 1))) != 0;
        let ok_ascent = p <= 1 || (target & (1 << (p - 2))) == 0;
        if !ok_descent || !ok_ascent {
            continue;
        }

        // Compute source descent set S' ⊆ [n-2]
        // Removing position p from the permutation:
        // For positions j < p-1 in σ: Des(π) inherits Des(σ)
        // Position p-1 in σ is forced ascent (< n), position p is forced descent (n >)
        // After removing position p, the remaining positions shift:
        //   π(j) = σ(j) for j < p
        //   π(j) = σ(j+1) for j ≥ p
        // Des(π) at position j:
        //   j < p-2: same as Des(σ) at j (bit j in target)
        //   j = p-2: π(p-2) vs π(p-1) = σ(p-2) vs σ(p-1) when p>2
        //            σ(p-1) < n and σ(p-2) could be anything — this is NOT determined by target alone
        //   j = p-1: π(p-1) vs π(p) = σ(p-1) vs σ(p+1) — also not determined

        // Actually, the source descent set at position p-1 depends on the specific
        // permutation, not just target and p. So there are MULTIPLE possible source
        // descent sets for each p.

        // For the forward DP, we don't need to enumerate source sets this way.
        // Instead, during the forward step, we'll compute the new descent set
        // from the old one + insertion position.

        // For backward reachability, we need: which old descent sets D' can
        // produce target when we insert n at position p?

        // New Des from old D' and insertion at p:
        // For p=1: Des(σ) = {1} ∪ {j+1 : j ∈ D'}
        // For p=n: Des(σ) = D'
        // For 1<p<n: Des(σ) = (D' ∩ [p-2]) ∪ {p} ∪ {j+1 : j ∈ D', j ≥ p}
        //   BUT position p-1 is forced ascent, and the bit at p-1 in D'
        //   gets absorbed (whether D' had p-1 as descent or not, σ has ascent there).
        //   Also position p-1 in π maps to positions p-1 in σ (which is ascent)
        //   and position p in π maps to position p+1 in σ.

        // Let me think about it differently. Given target Des and p, what D' works?
        //
        // For p=n: D' = target (restricted to [n-2], which target already is since n-1 ∉ target for valid p=n)
        //   Actually target ⊆ [n-1], and D' ⊆ [n-2]. If p=n is valid then n-1 ∉ target.
        //   D' = target ∩ [n-2] = target (since highest bit < n-1).
        //   But wait, we need target bits up to n-2. target has bits 0..n-2 (positions 1..n-1).
        //   D' has bits 0..n-3 (positions 1..n-2). So D' = target & ((1<<(n-2))-1).
        //
        // For p=1: target = {1} ∪ {j+1 : j ∈ D'}, so D' = {j-1 : j ∈ target, j ≥ 2}
        //   D' bit k = target bit (k+1), for k = 0..n-3.
        //
        // For 1<p<n: target restricted:
        //   - bits 0..p-3 (positions 1..p-2): must equal D' bits 0..p-3
        //   - bit p-2 (position p-1): must be 0 in target (forced ascent) — already checked
        //   - bit p-1 (position p): must be 1 in target (forced descent) — already checked
        //   - bits p..n-2 (positions p+1..n-1): target bit j = D' bit (j-1)
        //
        // BUT: position p-1 in D' is free! (It was position p-1 in π, but in σ
        // this becomes a forced ascent regardless of D'.)
        // So D' at bit p-2 can be 0 or 1 — both give the same target.

        // This means for 1<p<n, there are TWO possible source descent sets:
        // one with bit p-2 set and one without.

        if p == n as usize {
            let d_prime = target & ((1u64 << (n - 2)) - 1);
            result.push((p, d_prime));
        } else if p == 1 {
            let mut d_prime = 0u64;
            for j in 1..n as u64 - 1 {
                if target & (1 << j) != 0 {
                    d_prime |= 1 << (j - 1);
                }
            }
            result.push((p, d_prime));
        } else {
            // 1 < p < n: two possible source sets (bit p-2 free)
            for bit_val in 0..=1u64 {
                let mut d_prime = 0u64;
                // bits 0..p-3: from target
                for j in 0..p.saturating_sub(2) {
                    if target & (1 << j) != 0 {
                        d_prime |= 1 << j;
                    }
                }
                // bit p-2: free
                d_prime |= bit_val << (p - 2);
                // bits p-1..n-3: from target bits p..n-2
                for j in p as u64..n as u64 - 1 {
                    if target & (1 << j) != 0 {
                        d_prime |= 1 << (j - 1);
                    }
                }
                result.push((p, d_prime));
            }
        }
    }
    result
}

/// Compute backward-reachable descent sets starting from the alternating
/// descent set at size `max_n`, going back to size 1.
fn compute_relevant_sets(max_n: u8) -> Vec<HashSet<u64>> {
    let mut relevant: Vec<HashSet<u64>> = vec![HashSet::new(); max_n as usize + 1];

    // At size max_n, only the alternating descent set matters
    relevant[max_n as usize].insert(alt_des(max_n));

    // Work backwards
    for n in (2..=max_n).rev() {
        let current_relevant = relevant[n as usize].clone();
        for &target in &current_relevant {
            let sources = source_descent_sets(target, n);
            for (_, src_des) in sources {
                relevant[n as usize - 1].insert(src_des);
            }
        }
    }

    relevant
}

type State = (u64, u16, u8, u64); // (des, swaps, pos_max, adj)

fn insert_step_pruned(
    table: &HashMap<State, i64>,
    old_n: u8,
    relevant_new: &HashSet<u64>,
) -> HashMap<State, i64> {
    let new_n = old_n + 1;
    let mut new_table: HashMap<State, i64> = HashMap::new();

    for (&(des, sw, pos_max, adj), &count) in table {
        if count == 0 {
            continue;
        }
        let pos_max = pos_max as usize;

        for p in 1..=new_n as usize {
            // Compute new descent set
            let new_des: u64 = if p == 1 {
                let mut d = 1u64; // bit 0
                for j in 1..old_n as u64 {
                    if des & (1 << (j - 1)) != 0 {
                        d |= 1 << j;
                    }
                }
                d
            } else if p == new_n as usize {
                des
            } else {
                let mut d = 0u64;
                for j in 0..p.saturating_sub(2) {
                    if des & (1 << j) != 0 {
                        d |= 1 << j;
                    }
                }
                d |= 1 << (p - 1); // forced descent at p
                for j in (p - 1)..old_n as usize - 1 {
                    if des & (1 << j) != 0 {
                        d |= 1 << (j + 1);
                    }
                }
                d
            };

            // Prune: only keep if new_des is relevant
            if !relevant_new.contains(&new_des) {
                continue;
            }

            // ε₁
            let e1: u16 = if p > 1 && p < new_n as usize {
                if p >= 2 && (adj & (1 << (p - 2))) != 0 {
                    1
                } else {
                    0
                }
            } else {
                0
            };

            // ε₂
            let e2: u16 = if pos_max + 2 <= p { 1 } else { 0 };

            let new_sw = sw + e1 + e2;
            let new_pos_max = p as u8;

            // New adjacency bitmask
            let mut new_adj: u64 = 0;
            for j in 0..new_n as usize - 1 {
                let is_adj = if j + 2 <= p.saturating_sub(1) {
                    (adj & (1 << j)) != 0
                } else if j + 2 == p {
                    p >= 2 && pos_max == p - 1
                } else if j + 1 == p {
                    false
                } else if j >= p {
                    j >= 1 && (adj & (1 << (j - 1))) != 0
                } else {
                    false
                };
                if is_adj {
                    new_adj |= 1 << j;
                }
            }

            *new_table
                .entry((new_des, new_sw, new_pos_max, new_adj))
                .or_insert(0) += count;
        }
    }

    new_table
}

fn extract_alternating_poly(table: &HashMap<State, i64>, n: u8) -> Vec<i64> {
    let target = alt_des(n);
    let max_sw = table
        .keys()
        .filter(|&&(des, _, _, _)| des == target)
        .map(|&(_, sw, _, _)| sw as usize)
        .max()
        .unwrap_or(0);
    let mut coeffs = vec![0i64; max_sw + 1];
    for (&(des, sw, _, _), &count) in table {
        if des == target {
            coeffs[sw as usize] += count;
        }
    }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
        coeffs.pop();
    }
    coeffs
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(18);

    println!("Computing relevant descent sets...");
    let t0 = Instant::now();
    let relevant = compute_relevant_sets(max_n);
    println!("Done in {:?}\n", t0.elapsed());

    for n in 1..=max_n {
        println!(
            "  n={:>2}: {} relevant descent sets (of {} total)",
            n,
            relevant[n as usize].len(),
            1u64 << (n.saturating_sub(1))
        );
    }

    println!("\n═══ H_n(t) via pruned DP ═══\n");

    let mut table: HashMap<State, i64> = HashMap::new();
    table.insert((0, 0, 1, 0), 1);

    let h1 = extract_alternating_poly(&table, 1);
    println!("H_{:<2}(t) = {}", 1, format_poly(&h1));

    let mut h_polys: Vec<Vec<i64>> = vec![h1];

    for n in 2..=max_n {
        let t = Instant::now();
        table = insert_step_pruned(&table, n - 1, &relevant[n as usize]);
        let elapsed = t.elapsed();

        let poly = extract_alternating_poly(&table, n);
        let rr = poly.len() <= 2 || is_real_rooted(&poly);
        let gp = is_gamma_positive(&poly);
        let states = table.len();

        let poly_str = format_poly(&poly);
        let poly_display = if poly_str.len() > 90 {
            format!("{}...", &poly_str[..87])
        } else {
            poly_str
        };

        println!(
            "H_{:<2}(t) = {:<92} {:>10.2?} st={:<10} rr:{} γ+:{}",
            n, poly_display, elapsed, states, rr, gp,
        );
        h_polys.push(poly);

        if states > 100_000_000 {
            println!("  (stopping: state table > 100M)");
            break;
        }
    }

    // Print coefficient table
    println!("\n═══ Coefficient table ═══\n");
    for (i, p) in h_polys.iter().enumerate() {
        let n = i + 1;
        let cs: Vec<String> = p.iter().map(|c| c.to_string()).collect();
        println!("n={:>2}: [{}]", n, cs.join(", "));
    }

    // γ-coefficients
    println!("\n═══ γ-coefficients ═══\n");
    for (i, p) in h_polys.iter().enumerate() {
        let n = i + 1;
        if let Some(g) = gamma_coefficients(p) {
            let gs: Vec<String> = g.iter().map(|c| c.to_string()).collect();
            println!("n={:>2}: [{}]", n, gs.join(", "));
        }
    }
}
