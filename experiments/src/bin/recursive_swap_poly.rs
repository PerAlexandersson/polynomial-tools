/// Recursive computation of L_{n,S}(t) via the insertion lemma,
/// and exploration of recurrences when the descent set grows with n.
///
/// Part 1: Recursive computation of H_n(t) (alternating case).
///   Build H_n from H_{n-1} data using the insertion lemma.
///   State: (source descent set, position of max, adjacency info at each valid insertion point)
///   This avoids enumerating all permutations — instead we maintain a DP table.
///
/// Part 2: Growing descent sets.
///   Fix a "base" descent set pattern and grow it as n increases:
///   S_n = S_{n-1} ∪ {n-1} (add descent at the new last position)
///   or S_n = {k, 2k, ...} ∩ [n-1] (k-alternating pattern)
///   Look for polynomial recurrences in n for these families.
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};
use polynomial_tools::real_rootedness::{format_poly, is_real_rooted, check_weak_interlacing};
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};
use std::collections::BTreeMap;
use std::time::Instant;

fn build_poly_from_perms(perms: &[&Vec<u8>]) -> Vec<i64> {
    if perms.is_empty() { return vec![0]; }
    let max_s = perms.iter().map(|s| compute(s, Stat::Swaps)).max().unwrap_or(0);
    let mut coeffs = vec![0i64; max_s + 1];
    for s in perms { coeffs[compute(s, Stat::Swaps)] += 1; }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 { coeffs.pop(); }
    coeffs
}

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut c = vec![0i64; len];
    for (i, &v) in a.iter().enumerate() { c[i] += v; }
    for (i, &v) in b.iter().enumerate() { c[i] += v; }
    c
}

fn descent_set_str(s: u64, n: u8) -> String {
    let positions: Vec<String> = (0..n - 1)
        .filter(|&i| s & (1 << i) != 0)
        .map(|i| (i + 1).to_string())
        .collect();
    if positions.is_empty() { "∅".to_string() }
    else { format!("{{{}}}", positions.join(",")) }
}

/// DP-based recursive computation of L_{n,S}(t) using the insertion lemma.
///
/// State: for each permutation σ, we track:
///   - descent set (bitmask)
///   - swaps count
///   - position of max value (1-indexed)
///   - for each pair of adjacent positions (i, i+1): whether σ(i)+1 = σ(i+1) (adjacency)
///     encoded as a bitmask: bit j set iff σ(j)+1 = σ(j+1), 0-indexed
///
/// The tuple (des_bitmask, swaps, pos_max, adj_bitmask) defines the state.
/// We maintain a frequency table: state → count.
///
/// When inserting value n+1 at position p:
///   - ε₂ = 1 iff pos_max ≤ p-2
///   - ε₁ = 1 iff p>1 and p<n+1 and bit (p-2) of adj_bitmask is set
///     (i.e., σ(p-1)+1 = σ(p) in the pre-insertion permutation)
///   - new descent set: computed from old des + p
///   - new pos_max: p (since n+1 is the new max)... wait, no:
///     We're inserting n+1, so the new max is at position p.
///     But we need to track the STATE of the new permutation for future insertions.
///     After inserting n+1 at position p, the new pos_max = p.
///
/// Adjacency update after inserting n+1 at position p:
///   The adjacency bitmask changes:
///   - For positions before p-1 in the new perm: same as old (0-indexed: j < p-2)
///   - Position p-1 in new perm: new pair is (σ(p-1), n+1). σ(p-1)+1 = n+1 iff σ(p-1) = n.
///     So adj bit at p-2 (0-indexed) in new perm = (old value at pos p-1 == n) = (pos_max == p-1?
///     No, pos_max is position of the old max value n, so σ(pos_max) = n.)
///     adj[p-2] in new = 1 iff pos_max == p-1 (in old perm, 1-indexed)
///   - Position p in new perm: pair is (n+1, σ'(p+1)). n+1 is max, so (n+1)+1 doesn't exist.
///     But we need σ(p)+1 = σ(p+1) in new. σ(p) = n+1, σ(p+1) = old σ(p) (shifted).
///     (n+1)+1 = n+2 which doesn't exist, so this is never an adjacency.
///     Wait, adjacency is σ(j)+1 = σ(j+1). At position p: σ(p)=n+1, σ(p+1)=old_σ(p).
///     n+1+1 = n+2 > n+1 ≥ all values. So adj[p-1] = 0 in new perm.
///   - Position p+1 in new perm: pair (old_σ(p), old_σ(p+1)). This was the pair at (p, p+1) in old.
///     WAIT — the insertion shifts positions ≥ p by 1. In the new perm:
///     new_σ(j) = old_σ(j) for j < p
///     new_σ(p) = n+1
///     new_σ(j) = old_σ(j-1) for j > p
///
///   So adjacencies in new perm at position j (0-indexed pair (j, j+1)):
///     j < p-2: same as old adj[j]  (both values unchanged)
///     j = p-2: new_σ(p-1)+1 = new_σ(p)? → old_σ(p-1)+1 = n+1? → old_σ(p-1) = n → old pos_max = p-1
///     j = p-1: new_σ(p)+1 = new_σ(p+1)? → (n+1)+1 = old_σ(p)? → impossible (n+2 > all values). So 0.
///     j = p: new_σ(p+1)+1 = new_σ(p+2)? → old_σ(p)+1 = old_σ(p+1)? → old adj[p-1]
///     j > p: new_σ(j+1)+1 = new_σ(j+2)? → old_σ(j)+1 = old_σ(j+1)? → old adj[j-1]
///
///   Summary: new_adj[j] for j = 0..n-1 (n positions in new perm of size n+1):
///     j < p-2: old_adj[j]
///     j = p-2: 1 iff old pos_max == p-1 (and p >= 2)
///     j = p-1: 0 (always, since n+1 is max)
///     j >= p: old_adj[j-1]
///
/// This is fully determined by (old_adj, old_pos_max, p).

type State = (u64, usize, usize, u64); // (des_bitmask, swaps, pos_max, adj_bitmask)
type StateTable = BTreeMap<State, i64>;

fn insert_step(table: &StateTable, old_n: u8) -> StateTable {
    let new_n = old_n + 1;
    let mut new_table: StateTable = BTreeMap::new();

    for (&(des, sw, pos_max, adj), &count) in table {
        if count == 0 { continue; }

        // Try each insertion position p = 1..=new_n
        for p in 1..=new_n as usize {
            // Compute ε₁
            let e1 = if p > 1 && p < new_n as usize {
                // Check adjacency at position p-1 in old perm (0-indexed: p-2)
                if p >= 2 && (adj & (1 << (p - 2))) != 0 { 1 } else { 0 }
            } else {
                0
            };

            // Compute ε₂
            let e2 = if pos_max + 2 <= p { 1 } else { 0 };

            // New swaps count
            let new_sw = sw + e1 + e2;

            // New descent set
            let mut new_des: u64 = 0;
            if p == 1 {
                // Position 1: forced descent (n+1 > σ(2))
                new_des |= 1; // bit 0
                // Positions > 1: shifted from old
                for j in 1..old_n as u64 {
                    if des & (1 << (j - 1)) != 0 {
                        new_des |= 1 << j;
                    }
                }
            } else if p == new_n as usize {
                // Inserting at end: Des unchanged
                new_des = des;
            } else {
                // 1 < p < new_n
                // Positions 1..p-2: same as old
                for j in 0..p.saturating_sub(2) {
                    if des & (1 << j) != 0 {
                        new_des |= 1 << j;
                    }
                }
                // Position p-1: forced ascent (bit p-2 NOT set)
                // Position p: forced descent (bit p-1 set)
                new_des |= 1 << (p - 1);
                // Positions > p: shifted from old positions ≥ p
                for j in (p - 1)..old_n as usize - 1 {
                    if des & (1 << j) != 0 {
                        new_des |= 1 << (j + 1);
                    }
                }
            }

            // New position of max = p
            let new_pos_max = p;

            // New adjacency bitmask
            let mut new_adj: u64 = 0;
            for j in 0..new_n as usize - 1 {
                let is_adj = if j + 2 <= p.saturating_sub(1) {
                    // j < p-2: old_adj[j]
                    (adj & (1 << j)) != 0
                } else if j + 2 == p {
                    // j = p-2: 1 iff old pos_max == p-1
                    p >= 2 && pos_max == p - 1
                } else if j + 1 == p {
                    // j = p-1: always 0
                    false
                } else {
                    // j >= p: old_adj[j-1]
                    j >= 1 && (adj & (1 << (j - 1))) != 0
                };
                if is_adj {
                    new_adj |= 1 << j;
                }
            }

            *new_table.entry((new_des, new_sw, new_pos_max, new_adj)).or_insert(0) += count;
        }
    }

    new_table
}

/// Extract L_{n,S}(t) from the state table for a given descent set S.
fn extract_poly(table: &StateTable, target_des: u64) -> Vec<i64> {
    let max_sw = table.keys()
        .filter(|&&(des, _, _, _)| des == target_des)
        .map(|&(_, sw, _, _)| sw)
        .max()
        .unwrap_or(0);
    let mut coeffs = vec![0i64; max_sw + 1];
    for (&(des, sw, _, _), &count) in table {
        if des == target_des {
            coeffs[sw] += count;
        }
    }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 { coeffs.pop(); }
    coeffs
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(16);

    // Part 1: Generate H_n(t) recursively using DP
    println!("═══ Part 1: Recursive H_n(t) via DP ═══\n");

    // Base case: n=1, the single permutation [1].
    // des=0, swaps=0, pos_max=1, adj=0
    let mut table: StateTable = BTreeMap::new();
    table.insert((0, 0, 1, 0), 1);

    let mut h_polys: Vec<Vec<i64>> = Vec::new();

    // n=1 polynomial
    let h1 = extract_poly(&table, 0);
    println!("H_1(t) = {}", format_poly(&h1));
    h_polys.push(h1);

    for n in 2..=max_n {
        let t = Instant::now();
        table = insert_step(&table, n - 1);
        let elapsed = t.elapsed();

        // Alternating descent set for size n
        let alt_des = {
            let mut d = 0u64;
            for i in (2..=n as u64 - 1).step_by(2) {
                d |= 1 << (i - 1);
            }
            if n >= 2 { d } else { 0 }
        };

        let poly = extract_poly(&table, alt_des);
        let rr = poly.len() <= 2 || is_real_rooted(&poly);
        let states = table.len();

        println!(
            "H_{:<2}(t) = {:<70} {:>8.2?} states={:<8} {}",
            n, format_poly(&poly), elapsed, states,
            if rr { "✓ rr" } else { "✗" }
        );
        h_polys.push(poly);

        // If state table is getting huge, warn
        if states > 5_000_000 {
            println!("  (state table too large, stopping)");
            break;
        }
    }

    // Recurrence search on H_n
    println!("\n=== Recurrence search for H_n(t) ===");
    if h_polys.len() >= 6 {
        let opts = AdaptiveSearchOptions {
            max_rec_len: 4,
            max_var_deg: 3,
            max_idx_deg: 3,
            max_diff_deg: 2,
            try_inhomogeneous: true,
            try_denominator: true,
            max_denom_var_deg: 2,
            max_denom_idx_deg: 2,
            min_margin: 2,
            verbose: false,
        };
        match find_recurrence_adaptive(&h_polys, &opts) {
            Some(res) => println!("FOUND: {}", res.recurrence),
            None => println!("(no recurrence found with these parameters)"),
        }
    }

    // Part 2: Growing descent set families
    println!("\n═══ Part 2: Growing descent set families ═══\n");

    // Family A: S_n = {2, 4, ..., 2⌊(n-1)/2⌋} (alternating, already done above)

    // Family B: k-alternating for k=3
    // S_n = {3, 6, 9, ...} ∩ [n-1]
    println!("--- k=3 alternating: S_n = {{3,6,9,...}} ∩ [n-1] ---");
    let mut k3_polys: Vec<Vec<i64>> = Vec::new();
    // Recompute from scratch for this family
    let mut table3: StateTable = BTreeMap::new();
    table3.insert((0, 0, 1, 0), 1);

    for n in 1..=max_n.min(14) {
        if n > 1 {
            table3 = insert_step(&table3, n - 1);
        }
        let k3_des = {
            let mut d = 0u64;
            for i in (3..=n as u64 - 1).step_by(3) {
                d |= 1 << (i - 1);
            }
            d
        };
        let poly = extract_poly(&table3, k3_des);
        if poly.iter().sum::<i64>() > 0 || n <= 3 {
            let rr = poly.len() <= 2 || is_real_rooted(&poly);
            println!("H_{}^(3)(t) = {:<50} {}", n, format_poly(&poly), if rr { "✓" } else { "✗" });
            k3_polys.push(poly);
        }
        if table3.len() > 5_000_000 { break; }
    }

    if k3_polys.len() >= 6 {
        println!("\nRecurrence search for H_n^(3)(t):");
        let opts = AdaptiveSearchOptions::default();
        match find_recurrence_adaptive(&k3_polys, &opts) {
            Some(res) => println!("FOUND: {}", res.recurrence),
            None => println!("(no recurrence found)"),
        }
    }

    // Family C: "full descent at end" — S_n = {n-1}
    // (descent only at the last position)
    println!("\n--- Descent at last position: S_n = {{n-1}} ---");
    let mut last_des_polys: Vec<Vec<i64>> = Vec::new();
    let mut table_c: StateTable = BTreeMap::new();
    table_c.insert((0, 0, 1, 0), 1);
    for n in 1..=max_n.min(16) {
        if n > 1 {
            table_c = insert_step(&table_c, n - 1);
        }
        let des = if n >= 2 { 1u64 << (n - 2) } else { 0 };
        let poly = extract_poly(&table_c, des);
        if poly.iter().sum::<i64>() > 0 || n <= 2 {
            println!("L_{{n,{{n-1}}}}(t), n={}: {}", n, format_poly(&poly));
            last_des_polys.push(poly);
        }
        if table_c.len() > 5_000_000 { break; }
    }
    if last_des_polys.len() >= 6 {
        println!("\nRecurrence search for L_{{n,{{n-1}}}}(t):");
        let opts = AdaptiveSearchOptions::default();
        match find_recurrence_adaptive(&last_des_polys, &opts) {
            Some(res) => println!("FOUND: {}", res.recurrence),
            None => println!("(no recurrence found)"),
        }
    }

    // Family D: "growing from fixed base" — start with S={2}, extend to S={2,n-1}
    // i.e., S_n = {2, n-1} for n ≥ 4
    println!("\n--- Des = {{2, n-1}} ---");
    let mut base_ext_polys: Vec<Vec<i64>> = Vec::new();
    let mut table_d: StateTable = BTreeMap::new();
    table_d.insert((0, 0, 1, 0), 1);
    for n in 1..=max_n.min(16) {
        if n > 1 {
            table_d = insert_step(&table_d, n - 1);
        }
        if n >= 4 {
            let des = (1u64 << 1) | (1u64 << (n - 2)); // bits 1 and n-2
            let poly = extract_poly(&table_d, des);
            if poly.iter().sum::<i64>() > 0 {
                let rr = poly.len() <= 2 || is_real_rooted(&poly);
                println!("L_{{n,{{2,n-1}}}}(t), n={}: {:<50} {}", n, format_poly(&poly), if rr { "✓" } else { "✗" });
                base_ext_polys.push(poly);
            }
        }
        if table_d.len() > 5_000_000 { break; }
    }
    if base_ext_polys.len() >= 6 {
        println!("\nRecurrence search:");
        let opts = AdaptiveSearchOptions::default();
        match find_recurrence_adaptive(&base_ext_polys, &opts) {
            Some(res) => println!("FOUND: {}", res.recurrence),
            None => println!("(no recurrence found)"),
        }
    }

    // Family E: Staircase — S_n = {1, 2, ..., ⌊n/2⌋}
    println!("\n--- Des = {{1,2,...,⌊n/2⌋}} ---");
    let mut stair_polys: Vec<Vec<i64>> = Vec::new();
    let mut table_e: StateTable = BTreeMap::new();
    table_e.insert((0, 0, 1, 0), 1);
    for n in 1..=max_n.min(16) {
        if n > 1 {
            table_e = insert_step(&table_e, n - 1);
        }
        if n >= 3 {
            let k = n / 2;
            let des = (1u64 << k) - 1; // bits 0..k-1
            let poly = extract_poly(&table_e, des);
            if poly.iter().sum::<i64>() > 0 {
                let rr = poly.len() <= 2 || is_real_rooted(&poly);
                println!("L(n={}, Des={{1..{}}}): {:<50} {}", n, k, format_poly(&poly), if rr { "✓" } else { "✗" });
                stair_polys.push(poly);
            }
        }
        if table_e.len() > 5_000_000 { break; }
    }
    if stair_polys.len() >= 6 {
        println!("\nRecurrence search:");
        let opts = AdaptiveSearchOptions::default();
        match find_recurrence_adaptive(&stair_polys, &opts) {
            Some(res) => println!("FOUND: {}", res.recurrence),
            None => println!("(no recurrence found)"),
        }
    }
}
