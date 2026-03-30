/// Search for small non-real-rooted h* by deleting values from known examples.
///
/// For each known non-RR permutation, try all single deletions (remove a value,
/// standardize), check real-rootedness. Iterate downward.
use combinatoric_core::poset::Poset;
use polynomial_tools::real_rootedness;
use rayon::prelude::*;
use std::collections::BTreeSet;

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

/// Remove the value at position `pos` and standardize.
fn delete_and_standardize(w: &[u8], pos: usize) -> Vec<u8> {
    let removed = w[pos];
    let mut result: Vec<u8> = Vec::with_capacity(w.len() - 1);
    for (i, &v) in w.iter().enumerate() {
        if i == pos { continue; }
        // Standardize: subtract 1 from values > removed
        if v > removed { result.push(v - 1); } else { result.push(v); }
    }
    result
}

/// Generate all single deletions of w.
fn all_deletions(w: &[u8]) -> Vec<Vec<u8>> {
    (0..w.len()).map(|pos| delete_and_standardize(w, pos)).collect()
}

fn longest_dec_subseq(w: &[u8]) -> usize {
    let mut tails: Vec<u8> = Vec::new();
    for &x in w {
        let rx = 255 - x;
        match tails.binary_search(&rx) {
            Ok(_) => {}
            Err(pos) => {
                if pos == tails.len() { tails.push(rx); } else { tails[pos] = rx; }
            }
        }
    }
    tails.len()
}

fn main() {
    // All 22 non-RR from Stembridge neighborhood (n=17) + the 321-containing one
    let seeds: Vec<Vec<u8>> = vec![
        vec![2,4,6,8,10,1,12,3,15,5,17,7,9,11,13,14,16],
        vec![2,5,6,8,10,1,12,3,15,4,17,7,9,11,13,14,16],
        vec![2,4,5,8,10,1,12,3,15,6,17,7,9,11,13,14,16],
        vec![2,4,7,8,10,1,12,3,15,5,17,6,9,11,13,14,16],
        vec![2,4,6,8,9,1,12,3,15,5,17,7,10,11,13,14,16],
        vec![2,4,6,8,11,1,12,3,15,5,17,7,9,10,13,14,16],
        vec![2,5,7,8,10,1,12,3,15,4,17,6,9,11,13,14,16],
        vec![2,5,6,8,11,1,12,3,15,4,17,7,9,10,13,14,16],
        vec![2,4,5,9,10,1,12,3,15,6,17,7,8,11,13,14,16],
        vec![2,4,5,8,9,1,12,3,15,6,17,7,10,11,13,14,16],
        vec![2,4,5,8,11,1,12,3,15,6,17,7,9,10,13,14,16],
        vec![2,4,5,8,10,1,12,3,15,6,17,7,9,11,13,14,16],  // dup but harmless
        vec![2,4,6,8,10,1,12,3,15,5,17,7,9,11,13,14,16],  // Stembridge dup
        vec![2,4,6,8,10,1,12,3,15,5,17,7,9,11,13,14,16],
        vec![3,4,6,8,10,12,2,1,15,5,17,7,9,11,13,14,16],  // 321-containing
        vec![2,6,7,8,10,1,12,3,15,4,17,5,9,11,13,14,16],
        vec![2,5,7,8,11,1,12,3,15,4,17,6,9,10,13,14,16],
        vec![2,4,5,7,9,1,12,3,15,6,17,8,10,11,13,14,16],
        vec![2,4,5,9,11,1,12,3,15,6,17,7,8,10,13,14,16],
        vec![2,4,5,8,10,1,12,3,15,6,17,7,9,11,13,14,16],
        vec![2,4,6,8,10,1,12,3,15,5,17,7,9,11,13,14,16],
        vec![2,4,6,8,10,1,12,3,15,5,17,7,9,11,13,14,16],
        // Also the n=19 and n=21 winners
        vec![2,4,5,8,10,1,12,3,15,6,17,7,19,9,11,13,14,16,18],
        vec![2,4,5,8,10,1,12,3,15,6,17,7,19,9,21,11,13,14,16,18,20],
    ];

    let min_target = 12;

    // Track all non-RR perms found, by size
    let mut found: BTreeSet<Vec<u8>> = BTreeSet::new();
    for w in &seeds { found.insert(w.clone()); }

    // Collect all seeds by size, process from largest down
    let max_n = seeds.iter().map(|w| w.len()).max().unwrap();
    let mut current_level: Vec<Vec<u8>> = seeds.iter()
        .filter(|w| w.len() == max_n)
        .cloned()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();

    // Also add smaller seeds directly
    let mut extra_by_n: std::collections::BTreeMap<usize, Vec<Vec<u8>>> = std::collections::BTreeMap::new();
    for w in &seeds {
        if w.len() < max_n {
            extra_by_n.entry(w.len()).or_default().push(w.clone());
        }
    }

    println!("=== Searching downward by deletion + standardization ===\n");
    println!("Seeds: {} unique permutations, n={} to {}\n", found.len(), min_target+1, max_n);

    let mut n = max_n;
    while n > min_target && !current_level.is_empty() {
        // Generate all deletions
        let mut candidates: BTreeSet<Vec<u8>> = BTreeSet::new();
        for w in &current_level {
            for del in all_deletions(w) {
                candidates.insert(del);
            }
        }

        let cands: Vec<Vec<u8>> = candidates.into_iter().collect();
        let target_n = n - 1;

        // Test in parallel
        let non_rr: Vec<(Vec<u8>, Vec<i64>)> = cands
            .par_iter()
            .filter_map(|w| {
                let p = perm_poset(w);
                let pe = p.p_eulerian_polynomial();
                if !real_rootedness::is_real_rooted(&pe) {
                    Some((w.clone(), pe))
                } else {
                    None
                }
            })
            .collect();

        println!(
            "n={}: {}/{} deletions are non-real-rooted",
            target_n,
            non_rr.len(),
            cands.len()
        );

        if !non_rr.is_empty() {
            for (w, pe) in &non_rr {
                let lds = longest_dec_subseq(w);
                println!(
                    "  lds={} w={:?}  h*={}",
                    lds, w, real_rootedness::format_poly(pe)
                );
            }
            println!();
        }

        // Prepare next round: merge with any extra seeds at this size
        current_level = non_rr.iter().map(|(w, _)| w.clone()).collect();
        if let Some(extras) = extra_by_n.get(&target_n) {
            for w in extras {
                if !current_level.contains(w) {
                    // Check it's actually non-RR
                    let p = perm_poset(w);
                    let pe = p.p_eulerian_polynomial();
                    if !real_rootedness::is_real_rooted(&pe) {
                        current_level.push(w.clone());
                    }
                }
            }
        }
        for w in &current_level {
            found.insert(w.clone());
        }
        n = target_n;
    }

    // Extra: try transpositions at n=16 from deletions of n=17 non-RR
    let n17_perms: Vec<Vec<u8>> = found.iter().filter(|w| w.len() == 17).cloned().collect();
    let mut n16_cands: BTreeSet<Vec<u8>> = BTreeSet::new();
    for w in &n17_perms {
        for del in all_deletions(w) {
            n16_cands.insert(del.clone());
            // Also try 1 transposition from each deletion
            for i in 0..del.len() {
                for j in (i + 1)..del.len() {
                    let mut t = del.clone();
                    t.swap(i, j);
                    n16_cands.insert(t);
                }
            }
        }
    }
    let n16_vec: Vec<Vec<u8>> = n16_cands.into_iter().collect();
    let n16_rr: Vec<(Vec<u8>, Vec<i64>)> = n16_vec
        .par_iter()
        .filter_map(|w| {
            let p = perm_poset(w);
            let pe = p.p_eulerian_polynomial();
            if !real_rootedness::is_real_rooted(&pe) {
                Some((w.clone(), pe))
            } else {
                None
            }
        })
        .collect();
    println!(
        "\nn=16 extended search (deletions + 1 transposition): {}/{} not real-rooted",
        n16_rr.len(),
        n16_vec.len()
    );
    for (w, pe) in n16_rr.iter().take(5) {
        println!("  w={:?}\n  h*={}\n", w, real_rootedness::format_poly(pe));
    }

    // Summary
    println!("\n=== Summary ===\n");
    let mut by_n: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
    for w in &found {
        *by_n.entry(w.len()).or_insert(0) += 1;
    }
    for (n, count) in &by_n {
        println!("n={}: {} non-real-rooted permutations found", n, count);
    }
    let smallest = found.iter().min_by_key(|w| w.len());
    if let Some(w) = smallest {
        println!("\nSmallest: n={}, w={:?}", w.len(), w);
        let p = perm_poset(w);
        let pe = p.p_eulerian_polynomial();
        println!("  h*={}", real_rootedness::format_poly(&pe));
    }
}
