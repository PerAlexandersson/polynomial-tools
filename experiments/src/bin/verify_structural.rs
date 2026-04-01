/// Verify three structural conjectures on pattern-avoiding permutations
/// for all n up to 8.
///
/// Conjecture 1 (Av_213): For i < n, position i is a right-to-left maximum
///   iff i is a descent.
///
/// Conjecture 2 (Av_231): DesSet(σ) and RlminSet(σ) partition [n].
///
/// Conjecture 3 (Av_321): ExcSet(σ) and RlminSet(σ) partition [n].
use combpoly::permutation::avoiding_permutations;
use combpoly::statistics::{compute_set, SetStat};
use std::collections::BTreeSet;

fn main() {
    let max_n: u8 = 8;

    println!(
        "=== Structural Conjecture Verification (n = 1..{}) ===\n",
        max_n
    );

    let mut all_ok = true;

    // -----------------------------------------------------------
    // Conjecture 1: Av_213 — RlmaxSet ∩ [n-1] == DesSet
    // (for positions i < n, i is a right-to-left max iff i is a descent)
    // -----------------------------------------------------------
    println!("Conjecture 1: Av_213 — DesSet == RlmaxSet ∩ [n-1]");
    for n in 1..=max_n {
        let perms = avoiding_permutations(n, &[vec![2, 1, 3]]);
        let mut fail_count = 0;
        for sigma in &perms {
            let des = compute_set(sigma, SetStat::DesSet);
            let rlmax = compute_set(sigma, SetStat::RlmaxSet);
            // Restrict rlmax to positions < n (i.e., 1-indexed positions in [1..n-1])
            let rlmax_restricted: BTreeSet<usize> =
                rlmax.iter().copied().filter(|&i| i < n as usize).collect();
            if des != rlmax_restricted {
                println!(
                    "  FAIL n={}: σ={:?}, DesSet={:?}, RlmaxSet∩[n-1]={:?}",
                    n, sigma, des, rlmax_restricted
                );
                fail_count += 1;
            }
        }
        if fail_count == 0 {
            println!("  n={}: OK ({} permutations checked)", n, perms.len());
        } else {
            println!("  n={}: {} FAILURES out of {}", n, fail_count, perms.len());
            all_ok = false;
        }
    }
    println!();

    // -----------------------------------------------------------
    // Conjecture 2: Av_231 — DesSet ∪ RlminSet = [n], disjoint
    // -----------------------------------------------------------
    println!("Conjecture 2: Av_231 — DesSet and RlminSet partition [n]");
    for n in 1..=max_n {
        let perms = avoiding_permutations(n, &[vec![2, 3, 1]]);
        let full: BTreeSet<usize> = (1..=n as usize).collect();
        let mut fail_count = 0;
        for sigma in &perms {
            let des = compute_set(sigma, SetStat::DesSet);
            let rlmin = compute_set(sigma, SetStat::RlminSet);
            let union: BTreeSet<usize> = des.union(&rlmin).copied().collect();
            let intersection: BTreeSet<usize> = des.intersection(&rlmin).copied().collect();
            if union != full || !intersection.is_empty() {
                println!(
                    "  FAIL n={}: σ={:?}, DesSet={:?}, RlminSet={:?}, union={:?}, inter={:?}",
                    n, sigma, des, rlmin, union, intersection
                );
                fail_count += 1;
            }
        }
        if fail_count == 0 {
            println!("  n={}: OK ({} permutations checked)", n, perms.len());
        } else {
            println!("  n={}: {} FAILURES out of {}", n, fail_count, perms.len());
            all_ok = false;
        }
    }
    println!();

    // -----------------------------------------------------------
    // Conjecture 3: Av_321 — ExcSet ∪ RlminSet = [n], disjoint
    // -----------------------------------------------------------
    println!("Conjecture 3: Av_321 — ExcSet and RlminSet partition [n]");
    for n in 1..=max_n {
        let perms = avoiding_permutations(n, &[vec![3, 2, 1]]);
        let full: BTreeSet<usize> = (1..=n as usize).collect();
        let mut fail_count = 0;
        for sigma in &perms {
            let exc = compute_set(sigma, SetStat::ExcSet);
            let rlmin = compute_set(sigma, SetStat::RlminSet);
            let union: BTreeSet<usize> = exc.union(&rlmin).copied().collect();
            let intersection: BTreeSet<usize> = exc.intersection(&rlmin).copied().collect();
            if union != full || !intersection.is_empty() {
                println!(
                    "  FAIL n={}: σ={:?}, ExcSet={:?}, RlminSet={:?}, union={:?}, inter={:?}",
                    n, sigma, exc, rlmin, union, intersection
                );
                fail_count += 1;
            }
        }
        if fail_count == 0 {
            println!("  n={}: OK ({} permutations checked)", n, perms.len());
        } else {
            println!("  n={}: {} FAILURES out of {}", n, fail_count, perms.len());
            all_ok = false;
        }
    }
    println!();

    // -----------------------------------------------------------
    // Summary
    // -----------------------------------------------------------
    if all_ok {
        println!("ALL THREE CONJECTURES VERIFIED for n = 1..{}.", max_n);
    } else {
        println!("SOME CONJECTURES FAILED — see details above.");
        std::process::exit(1);
    }
}
