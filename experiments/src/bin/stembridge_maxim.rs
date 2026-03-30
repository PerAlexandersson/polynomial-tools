/// Search from the permutation with largest imaginary root, using
/// both transpositions and insertions of (n+1).
use combinatoric_core::poset::Poset;
use polynomial_tools::real_rootedness;
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

fn perm_poset(w: &[u8]) -> Poset {
    let n = w.len();
    let mut edges = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if w[i] < w[j] {
                edges.push((i, j));
            }
        }
    }
    Poset::new(n, &edges).to_hasse_diagram()
}

fn avoids_4321(w: &[u8]) -> bool {
    let n = w.len();
    for a in 0..n {
        for b in (a + 1)..n {
            if w[a] <= w[b] { continue; }
            for c in (b + 1)..n {
                if w[b] <= w[c] { continue; }
                for d in (c + 1)..n {
                    if w[c] > w[d] { return false; }
                }
            }
        }
    }
    true
}

fn contains_321(w: &[u8]) -> bool {
    let n = w.len();
    for a in 0..n {
        for b in (a + 1)..n {
            if w[a] <= w[b] { continue; }
            for c in (b + 1)..n {
                if w[b] > w[c] { return true; }
            }
        }
    }
    false
}

/// BFS over arbitrary transpositions, collecting 4321-avoiding perms.
fn bfs_transpositions(start: &[u8], max_dist: usize) -> Vec<(Vec<u8>, String)> {
    let n = start.len();
    let mut visited: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut queue: VecDeque<(Vec<u8>, usize)> = VecDeque::new();
    let mut result = Vec::new();

    visited.insert(start.to_vec());
    queue.push_back((start.to_vec(), 0));

    while let Some((perm, dist)) = queue.pop_front() {
        if avoids_4321(&perm) {
            result.push((perm.clone(), format!("trans_d{}", dist)));
        }
        if dist < max_dist {
            for i in 0..n {
                for j in (i + 1)..n {
                    let mut nb = perm.clone();
                    nb.swap(i, j);
                    if visited.insert(nb.clone()) {
                        queue.push_back((nb, dist + 1));
                    }
                }
            }
        }
    }
    result
}

/// Insert value (n+1) at each position in w, producing n+1 permutations of [n+1].
fn insertions(w: &[u8]) -> Vec<(Vec<u8>, String)> {
    let n = w.len();
    let new_val = (n + 1) as u8;
    let mut result = Vec::new();
    for pos in 0..=n {
        let mut new_perm = Vec::with_capacity(n + 1);
        new_perm.extend_from_slice(&w[..pos]);
        new_perm.push(new_val);
        new_perm.extend_from_slice(&w[pos..]);
        if avoids_4321(&new_perm) {
            result.push((new_perm, format!("ins_{}_at_{}", new_val, pos)));
        }
    }
    result
}

/// Apply insertions to all permutations in a set, dedup.
fn insert_round(perms: &[(Vec<u8>, String)]) -> Vec<(Vec<u8>, String)> {
    let mut seen: BTreeSet<Vec<u8>> = BTreeSet::new();
    for (p, _) in perms {
        seen.insert(p.clone());
    }
    let mut result = Vec::new();
    for (p, label) in perms {
        for (ins, ins_label) in insertions(p) {
            if seen.insert(ins.clone()) {
                result.push((ins, format!("{}+{}", label, ins_label)));
            }
        }
    }
    result
}

fn main() {
    // Round 3: n=21 winner from previous round.
    // Origin: ins_20_at_19 + ins_21_at_14 from n=19 base.
    let mut w0: Vec<u8> = vec![2, 4, 5, 8, 10, 1, 12, 3, 15, 6, 17, 7, 19, 9, 11, 13, 14, 16, 18];
    w0.push(20);       // ins_20_at_19 (append at end)
    w0.insert(14, 21); // ins_21_at_14
    assert!(avoids_4321(&w0));
    let n = w0.len();

    println!("=== Starting from max-|Im| permutation (n={}) ===", n);
    println!("w0 = {:?}\n", w0);

    // Phase 1: Transpositions up to distance 3
    let max_dist: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(2);

    eprintln!("Phase 1: BFS transpositions, max distance {}...", max_dist);
    let trans = bfs_transpositions(&w0, max_dist);
    eprintln!("  4321-avoiding from transpositions: {}", trans.len());

    // Phase 2: Insertions of n+1 into w0 only (keep manageable at large n)
    eprintln!("Phase 2: Insertions of {} into w0...", n + 1);
    let base_for_ins = vec![(w0.clone(), "base".into())];
    let ins = insert_round(&base_for_ins);
    eprintln!("  4321-avoiding from insertions (n={}): {}", n + 1, ins.len());

    // Phase 3: skip for now — n+2 is too slow at this size
    let ins2: Vec<(Vec<u8>, String)> = Vec::new();

    // Combine all, dedup by permutation
    let mut all: BTreeMap<Vec<u8>, String> = BTreeMap::new();
    for (p, l) in trans.into_iter().chain(ins).chain(ins2) {
        all.entry(p).or_insert(l);
    }
    let all_perms: Vec<(Vec<u8>, String)> = all.into_iter().collect();
    eprintln!("Total unique 4321-avoiding permutations: {}", all_perms.len());

    // Test real-rootedness (parallel)
    eprintln!("Testing real-rootedness...");
    let t0 = std::time::Instant::now();

    let not_rr: Vec<(Vec<u8>, String, Vec<i64>)> = all_perms
        .par_iter()
        .filter_map(|(perm, label)| {
            let p = perm_poset(perm);
            let hstar = if perm.len() <= 19 {
                p.p_eulerian_polynomial()
            } else {
                p.order_polytope_hstar()
            };
            if !real_rootedness::is_real_rooted(&hstar) {
                Some((perm.clone(), label.clone(), hstar))
            } else {
                None
            }
        })
        .collect();

    eprintln!("Done in {:.1?}", t0.elapsed());

    // Group by n
    let mut by_n: BTreeMap<usize, (usize, usize)> = BTreeMap::new(); // n -> (tested, failed)
    for (p, _, _) in &not_rr {
        by_n.entry(p.len()).or_insert((0, 0)).1 += 1;
    }
    for (p, _) in &all_perms {
        by_n.entry(p.len()).or_insert((0, 0)).0 += 1;
    }

    println!("=== Results by size ===\n");
    for (n, (tested, failed)) in &by_n {
        println!("n={}: {}/{} not real-rooted", n, failed, tested);
    }
    println!("\nTotal: {}/{} not real-rooted\n", not_rr.len(), all_perms.len());

    // Full table, sorted by n then polynomial
    println!("=== Non-real-rooted examples ===\n");
    println!(
        "{:<4} {:<4} {:<5} {:<60} {:<20}",
        "n", "321?", "deg", "h* coefficients", "origin"
    );
    println!("{}", "-".repeat(110));

    let mut sorted = not_rr.clone();
    sorted.sort_by(|a, b| a.0.len().cmp(&b.0.len()).then(a.2.cmp(&b.2)));

    for (perm, label, pe) in &sorted {
        println!(
            "{:<4} {:<4} {:<5} {:<60} {}",
            perm.len(),
            if contains_321(perm) { "yes" } else { "no" },
            pe.len() - 1,
            format!("{:?}", pe),
            label
        );
    }

    // Output for Mathematica root-finding
    println!("\n=== Mathematica input (paste into NSolve) ===\n");
    println!("polys = {{");
    for (i, (_, _, pe)) in sorted.iter().enumerate() {
        let coeffs: Vec<String> = pe.iter().map(|c| c.to_string()).collect();
        println!(
            "  {{{}}}{}",
            coeffs.join(","),
            if i + 1 < sorted.len() { "," } else { "" }
        );
    }
    println!("}};");
}
