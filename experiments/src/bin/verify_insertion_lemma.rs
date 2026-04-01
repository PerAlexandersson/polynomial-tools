/// Verify the insertion lemma (Lemma 4.1 in FixedDescentSwaps.tex).
///
/// For π ∈ S_{n-1}, inserting the value n at position p gives σ ∈ S_n with:
///   swaps(σ) = swaps(π) + ε₁ + ε₂
///
/// where:
///   ε₁ = 1 iff 1 < p < n and π(p-1)+1 = π(p)  (insertion breaks an ascending adjacency)
///   ε₂ = 1 iff π⁻¹(n-1) ≤ p-2                  (n-1 is left of insertion, not adjacent)
///
/// Also verifies the descent set transformation:
/// For 1 < p < n:
///   Des(σ) = (Des(π) ∩ [p-2]) ∪ {p} ∪ {j+1 : j ∈ Des(π), j ≥ p}
///   (position p-1 is forced ascent, position p is forced descent)
/// For p = 1: Des(σ) = {1} ∪ {j+1 : j ∈ Des(π)}
/// For p = n: Des(σ) = Des(π)  (no new descent; n at the end)
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, descent_set_bitmask, Stat};

/// Insert value n at 1-indexed position p into π ∈ S_{n-1}.
/// Returns σ ∈ S_n.
fn insert(pi: &[u8], p: usize) -> Vec<u8> {
    let n = pi.len() as u8 + 1;
    let mut sigma = Vec::with_capacity(n as usize);
    for i in 0..p - 1 {
        sigma.push(pi[i]);
    }
    sigma.push(n);
    for i in p - 1..pi.len() {
        sigma.push(pi[i]);
    }
    sigma
}

/// Compute ε₁: did we break an ascending adjacency at the insertion point?
/// ε₁ = 1 iff 1 < p < n and π(p-1) + 1 = π(p)
/// (using 1-indexed positions in π)
fn epsilon1(pi: &[u8], p: usize) -> usize {
    let n = pi.len() + 1;
    if p <= 1 || p >= n {
        return 0;
    }
    // 1-indexed: π(p-1) is pi[p-2], π(p) is pi[p-1]
    if pi[p - 2] + 1 == pi[p - 1] {
        1
    } else {
        0
    }
}

/// Compute ε₂: does n-1 appear left of insertion and not adjacent?
/// ε₂ = 1 iff π⁻¹(n-1) ≤ p-2
/// (π⁻¹(n-1) is 1-indexed position of value n-1 in π)
fn epsilon2(pi: &[u8], p: usize) -> usize {
    let n_minus_1 = pi.len() as u8; // = n-1
                                    // Find 1-indexed position of n-1 in π
    let pos_nm1 = pi.iter().position(|&v| v == n_minus_1).unwrap() + 1; // 1-indexed
    if pos_nm1 + 2 <= p {
        // pos_nm1 ≤ p - 2
        1
    } else {
        0
    }
}

/// Predict the descent set bitmask of σ from π and insertion position p.
fn predicted_descent_bitmask(pi: &[u8], p: usize) -> u64 {
    let n = pi.len() + 1;
    let des_pi = descent_set_bitmask(pi);

    if p == n {
        // Inserting at end: Des(σ) = Des(π)
        // (n is last, σ(n) = n, no new descent since n is max)
        // But wait: σ(n-1) = π(n-1) < n = σ(n), so position n-1 is ascent.
        // The rest is unchanged.
        return des_pi;
    }

    let mut des_sigma: u64 = 0;

    if p == 1 {
        // Position 1 is forced descent (σ(1)=n > σ(2)=π(1))
        des_sigma |= 1; // bit 0 = position 1 descent
                        // Positions > 1: j ∈ Des(σ) iff j-1 ∈ Des(π)
        for j in 1..pi.len() {
            if des_pi & (1 << (j - 1)) != 0 {
                des_sigma |= 1 << j;
            }
        }
    } else {
        // 1 < p < n
        // Positions 1..p-2: same as Des(π)
        for j in 0..p.saturating_sub(2) {
            if des_pi & (1 << j) != 0 {
                des_sigma |= 1 << j;
            }
        }
        // Position p-1: forced ascent (σ(p-1) = π(p-1) < n = σ(p))
        // (bit p-2 is NOT set, regardless of π)
        // Position p: forced descent (σ(p) = n > σ(p+1) = π(p))
        des_sigma |= 1 << (p - 1); // bit p-1 = position p
                                   // Positions > p: j ∈ Des(σ) iff j-1 ∈ Des(π), for j > p
        for j in p..pi.len() {
            if des_pi & (1 << (j - 1)) != 0 {
                des_sigma |= 1 << j;
            }
        }
    }

    des_sigma
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    println!("=== Insertion Lemma Verification ===\n");

    let mut total_checks = 0u64;
    let mut swap_ok = 0u64;
    let mut swap_fail = 0u64;
    let mut des_ok = 0u64;
    let mut des_fail = 0u64;

    for n in 2..=max_n {
        let perms = all_permutations(n - 1);
        let mut n_checks = 0u64;
        let mut n_swap_fail = 0u64;
        let mut n_des_fail = 0u64;

        for pi in &perms {
            let sw_pi = compute(pi, Stat::Swaps);

            for p in 1..=n as usize {
                let sigma = insert(pi, p);
                let sw_sigma = compute(&sigma, Stat::Swaps);
                let e1 = epsilon1(pi, p);
                let e2 = epsilon2(pi, p);
                let predicted_swaps = sw_pi + e1 + e2;

                n_checks += 1;

                if sw_sigma != predicted_swaps {
                    n_swap_fail += 1;
                    if n_swap_fail <= 5 {
                        println!(
                            "  SWAP FAIL: n={} π={:?} p={} swaps(π)={} ε₁={} ε₂={} predicted={} actual={}",
                            n, pi, p, sw_pi, e1, e2, predicted_swaps, sw_sigma
                        );
                    }
                }

                // Verify descent set prediction
                let actual_des = descent_set_bitmask(&sigma);
                let predicted_des = predicted_descent_bitmask(pi, p);
                if actual_des != predicted_des {
                    n_des_fail += 1;
                    if n_des_fail <= 5 {
                        println!(
                            "  DES FAIL: n={} π={:?} p={} actual={:b} predicted={:b}",
                            n, pi, p, actual_des, predicted_des
                        );
                    }
                }
            }
        }

        total_checks += n_checks;
        swap_ok += n_checks - n_swap_fail;
        swap_fail += n_swap_fail;
        des_ok += n_checks - n_des_fail;
        des_fail += n_des_fail;

        println!(
            "n={:>2}: {:>8} insertions, swaps: {} ok / {} fail, des: {} ok / {} fail",
            n,
            n_checks,
            n_checks - n_swap_fail,
            n_swap_fail,
            n_checks - n_des_fail,
            n_des_fail,
        );
    }

    println!("\n=== Summary ===");
    println!("Total: {} insertions checked", total_checks,);
    println!("Swaps lemma: {} ok, {} failures", swap_ok, swap_fail);
    println!("Descent set: {} ok, {} failures", des_ok, des_fail);

    if swap_fail == 0 && des_fail == 0 {
        println!(
            "\n✓ Both the insertion lemma and descent set prediction verified for all n ≤ {}.",
            max_n
        );
    }
}
