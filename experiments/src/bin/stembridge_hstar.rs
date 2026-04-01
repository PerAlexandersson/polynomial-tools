/// Explore permutation posets near Stembridge's example for non-real-rooted h*.
///
/// Starting from Stembridge's 321-avoiding permutation, generates all
/// permutations within 1, 2, 3 adjacent transpositions, constructs the
/// permutation poset, and checks real-rootedness of h* = P-Eulerian.
use combinatoric_core::poset::Poset;
use polynomial_tools::real_rootedness;
use rayon::prelude::*;
use std::collections::{BTreeSet, VecDeque};

/// Construct the permutation poset: i <_P j iff i < j and w[i] < w[j].
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
    // Must reduce to Hasse diagram — transitive edges bloat the DP frontier.
    Poset::new(n, &edges).to_hasse_diagram()
}

/// BFS over arbitrary transpositions, keeping only 4321-avoiding permutations.
fn bfs_av4321(start: &[u8], max_dist: usize) -> (Vec<(Vec<u8>, usize)>, Vec<usize>, Vec<usize>) {
    let n = start.len();
    let mut visited: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut queue: VecDeque<(Vec<u8>, usize)> = VecDeque::new();
    let mut result = Vec::new();
    let mut bfs_by_dist = vec![0usize; max_dist + 1];
    let mut av_by_dist = vec![0usize; max_dist + 1];

    visited.insert(start.to_vec());
    queue.push_back((start.to_vec(), 0));

    while let Some((perm, dist)) = queue.pop_front() {
        bfs_by_dist[dist] += 1;
        if avoids_4321(&perm) {
            av_by_dist[dist] += 1;
            result.push((perm.clone(), dist));
        }
        if dist < max_dist {
            for i in 0..n {
                for j in (i + 1)..n {
                    let mut neighbor = perm.clone();
                    neighbor.swap(i, j);
                    if visited.insert(neighbor.clone()) {
                        queue.push_back((neighbor, dist + 1));
                    }
                }
            }
        }
    }
    (result, bfs_by_dist, av_by_dist)
}

fn avoids_4321(w: &[u8]) -> bool {
    let n = w.len();
    for a in 0..n {
        for b in (a + 1)..n {
            if w[a] <= w[b] {
                continue;
            }
            for c in (b + 1)..n {
                if w[b] <= w[c] {
                    continue;
                }
                for d in (c + 1)..n {
                    if w[c] > w[d] {
                        return false;
                    }
                }
            }
        }
    }
    true
}

fn main() {
    let stembridge: Vec<u8> = vec![2, 4, 6, 8, 10, 1, 12, 3, 15, 5, 17, 7, 9, 11, 13, 14, 16];
    let max_dist: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(4);

    eprintln!(
        "BFS from Stembridge's w (n={}), max distance {} (arbitrary transpositions)...",
        stembridge.len(),
        max_dist
    );
    let (av4321, bfs_by_dist, av_by_dist) = bfs_av4321(&stembridge, max_dist);
    let total_bfs: usize = bfs_by_dist.iter().sum();
    let total_av = av4321.len();
    eprintln!("BFS total: {}, 4321-avoiding: {}", total_bfs, total_av);

    println!("=== Counts ===\n");
    for d in 0..=max_dist {
        println!(
            "Distance {}: {} BFS total, {} are 4321-avoiding",
            d, bfs_by_dist[d], av_by_dist[d]
        );
    }
    println!("\nTotal 4321-avoiding to test: {}\n", total_av);

    // Test real-rootedness via Ehrhart h* (frontier DP + reciprocity, parallel)
    eprintln!(
        "Testing real-rootedness of {} permutation posets (via Ehrhart DP, parallel)...",
        total_av
    );
    let t0 = std::time::Instant::now();

    // Use p_eulerian_polynomial (streaming linext + des counting) — faster than
    // Ehrhart DP for these posets since 321-avoiding perms have few extensions.
    let results: Vec<(Vec<u8>, usize, Vec<i64>)> = av4321
        .par_iter()
        .filter_map(|(perm, dist)| {
            let p = perm_poset(perm);
            let pe = p.p_eulerian_polynomial();
            if !real_rootedness::is_real_rooted(&pe) {
                Some((perm.clone(), *dist, pe))
            } else {
                None
            }
        })
        .collect();

    eprintln!("Done in {:.1?}", t0.elapsed());

    let not_rr = results;
    let mut fail_by_dist = vec![0usize; max_dist + 1];
    for (_, dist, _) in &not_rr {
        fail_by_dist[*dist] += 1;
    }

    println!("=== Real-rootedness results ===\n");
    for d in 0..=max_dist {
        println!(
            "Distance {}: {}/{} not real-rooted",
            d, fail_by_dist[d], av_by_dist[d]
        );
    }
    println!("\nTotal: {}/{} not real-rooted\n", not_rr.len(), total_av);

    // Full table of all non-real-rooted examples
    fn contains_321(w: &[u8]) -> bool {
        let n = w.len();
        for a in 0..n {
            for b in (a + 1)..n {
                if w[a] <= w[b] {
                    continue;
                }
                for c in (b + 1)..n {
                    if w[b] > w[c] {
                        return true;
                    }
                }
            }
        }
        false
    }

    fn longest_dec_subseq(w: &[u8]) -> usize {
        // Patience sorting for LDS length
        let mut tails: Vec<u8> = Vec::new();
        for &x in w {
            // For decreasing: track longest decreasing = longest increasing of reverse
            let rx = 255 - x;
            match tails.binary_search(&rx) {
                Ok(_) => {}
                Err(pos) => {
                    if pos == tails.len() {
                        tails.push(rx);
                    } else {
                        tails[pos] = rx;
                    }
                }
            }
        }
        tails.len()
    }

    println!(
        "\n=== Full table of all {} non-real-rooted examples ===\n",
        not_rr.len()
    );
    println!(
        "{:<4} {:<4} {:<5} {:<5} {:<60} {}",
        "dist", "321?", "lds", "deg", "h* polynomial", "log-c"
    );
    println!("{}", "-".repeat(120));

    let mut count_321 = 0;
    for (perm, dist, pe) in &not_rr {
        let has_321 = contains_321(perm);
        if has_321 {
            count_321 += 1;
        }
        let lds = longest_dec_subseq(perm);
        let lc = real_rootedness::is_log_concave(pe);
        println!(
            "{:<4} {:<4} {:<5} {:<5} {:<60} {}",
            dist,
            if has_321 { "yes" } else { "no" },
            lds,
            pe.len() - 1,
            real_rootedness::format_poly(pe),
            if lc { "yes" } else { "no" }
        );
    }
    println!("\nContains 321: {}/{}", count_321, not_rr.len());

    // Print the full permutations for 321-containing ones
    let with_321: Vec<_> = not_rr.iter().filter(|(p, _, _)| contains_321(p)).collect();
    if !with_321.is_empty() {
        println!("\n=== Non-real-rooted examples containing 321 ===\n");
        for (perm, dist, pe) in &with_321 {
            println!("  dist={} w={:?}", dist, perm);
            println!("  h*={}", real_rootedness::format_poly(pe));
            println!();
        }
    }

    // Analysis: how many distinct polynomials? Group by polynomial.
    let mut poly_set: BTreeSet<Vec<i64>> = BTreeSet::new();
    for (_, _, pe) in &not_rr {
        poly_set.insert(pe.clone());
    }
    println!("=== Analysis ===\n");
    println!(
        "Distinct h* polynomials: {} (out of {} non-RR examples)",
        poly_set.len(),
        not_rr.len()
    );

    // Group by degree
    let mut by_deg: std::collections::BTreeMap<usize, Vec<&(Vec<u8>, usize, Vec<i64>)>> =
        std::collections::BTreeMap::new();
    for entry in &not_rr {
        by_deg.entry(entry.2.len() - 1).or_default().push(entry);
    }
    for (deg, entries) in &by_deg {
        println!("\nDegree {}: {} examples", deg, entries.len());
        // Leading coefficient
        let leading: Vec<i64> = entries
            .iter()
            .map(|(_, _, pe)| *pe.last().unwrap())
            .collect();
        let sums: Vec<i64> = entries
            .iter()
            .map(|(_, _, pe)| pe.iter().sum::<i64>())
            .collect();
        println!("  Leading coefficients: {:?}", leading);
        println!("  Sum of coeffs (= #linext): {:?}", sums);
        // Coefficient ranges per position
        let d = *deg;
        print!("  Coeff ranges: ");
        for k in 0..=d {
            let vals: Vec<i64> = entries.iter().map(|(_, _, pe)| pe[k]).collect();
            let lo = vals.iter().min().unwrap();
            let hi = vals.iter().max().unwrap();
            if lo == hi {
                print!("[{}] ", lo);
            } else {
                print!("[{}..{}] ", lo, hi);
            }
        }
        println!();
    }

    // How many positions differ between Stembridge and each neighbor?
    println!("\n--- Diff from Stembridge ---");
    let stem = &not_rr[0].2;
    for (perm, dist, pe) in &not_rr[1..] {
        let diffs: Vec<(usize, i64, i64)> = pe
            .iter()
            .enumerate()
            .zip(stem.iter().chain(std::iter::repeat(&0)))
            .filter(|((_, a), b)| a != b)
            .map(|((i, a), b)| (i, *b, *a))
            .collect();
        println!("  dist={} deg={} Δcoeffs={:?}", dist, pe.len() - 1, diffs);
    }
}
