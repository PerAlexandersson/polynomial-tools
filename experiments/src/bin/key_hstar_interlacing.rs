use combinatoric_core::key_polynomial::key_ehrhart_polynomial;
use combinatoric_core::partition::Partition;
use combpoly::permutation::all_permutations;
use polynomial_tools::real_rootedness::{
    check_weak_interlacing, ehrhart_to_hstar, format_poly, is_real_rooted,
};
/// Explore h*-vector interlacing for key polynomials along Bruhat order.
///
/// For a fixed partition λ and permutations σ in S_n, compute the Ehrhart
/// polynomial of the key polytope GT(λ, σ), extract the h*-vector, and check
/// whether h*-vectors interlace along covering relations in Bruhat order.
///
/// Question: if σ ⋗ τ in Bruhat order, does h*(λ,τ) ≪ h*(λ,σ)?
use std::collections::BTreeMap;

/// Inversion count (= Bruhat length).
fn inv_count(perm: &[u8]) -> usize {
    let n = perm.len();
    let mut c = 0;
    for i in 0..n {
        for j in i + 1..n {
            if perm[i] > perm[j] {
                c += 1;
            }
        }
    }
    c
}

/// Bruhat covering relations: σ ⋗ τ iff σ = τ·(a,b) with a < b, σ[a] > σ[b],
/// and no position c with a < c < b and σ[b] < σ[c] < σ[a].
fn bruhat_covers(sigma: &[u8]) -> Vec<Vec<u8>> {
    let n = sigma.len();
    let mut covers = Vec::new();
    for a in 0..n {
        for b in a + 1..n {
            if sigma[a] > sigma[b] {
                // Check no c in (a,b) with sigma[b] < sigma[c] < sigma[a]
                let lo = sigma[b];
                let hi = sigma[a];
                let blocked = (a + 1..b).any(|c| sigma[c] > lo && sigma[c] < hi);
                if !blocked {
                    let mut tau = sigma.to_vec();
                    tau.swap(a, b);
                    covers.push(tau);
                }
            }
        }
    }
    covers
}

fn perm_str(perm: &[u8]) -> String {
    perm.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join("")
}

fn main() {
    let test_cases: Vec<(Vec<u32>, usize)> = vec![
        (vec![2, 1], 3),
        (vec![3, 1], 3),
        (vec![2, 2], 3),
        (vec![3, 2, 1], 4),
        (vec![2, 1], 4),
        (vec![4, 2, 1], 4),
    ];

    for (parts, n) in &test_cases {
        let lambda = Partition::new(parts.clone());
        let n = *n;
        println!("========================================");
        println!("λ = {:?}, n = {}", parts, n);
        println!("========================================");

        let perms = all_permutations(n as u8);

        // Compute h*-vectors for all permutations
        let mut hstar_map: BTreeMap<Vec<u8>, Vec<i64>> = BTreeMap::new();
        let mut ehrhart_display: BTreeMap<Vec<u8>, String> = BTreeMap::new();

        for perm in &perms {
            let sigma: Vec<usize> = perm.iter().map(|&x| x as usize).collect();
            let ep = key_ehrhart_polynomial(&lambda, &sigma, None);
            let hstar = ehrhart_to_hstar(&ep.coeffs);
            ehrhart_display.insert(perm.clone(), ep.display());
            hstar_map.insert(perm.clone(), hstar);
        }

        // Print all h*-vectors, grouped by Bruhat length
        let mut by_length: BTreeMap<usize, Vec<&Vec<u8>>> = BTreeMap::new();
        for perm in &perms {
            by_length.entry(inv_count(perm)).or_default().push(perm);
        }

        println!("\nh*-vectors by Bruhat length:");
        for (len, perms_at_len) in &by_length {
            for perm in perms_at_len {
                let h = &hstar_map[*perm];
                let rr = if h.len() <= 1 {
                    true
                } else {
                    is_real_rooted(h)
                };
                println!(
                    "  ℓ={} σ={}: h* = {}  {}  [Ehrhart: {}]",
                    len,
                    perm_str(perm),
                    format_poly(h),
                    if rr { "RR" } else { "NOT-RR" },
                    ehrhart_display[*perm],
                );
            }
        }

        // Check interlacing along Bruhat covers
        println!("\nBruhat cover interlacing (σ ⋗ τ, checking h*(τ) ≪ h*(σ)):");
        let mut total = 0;
        let mut pass = 0;
        let mut fail_details = Vec::new();

        for sigma in &perms {
            for tau in &bruhat_covers(sigma) {
                let h_sigma = &hstar_map[sigma];
                let h_tau = &hstar_map[tau];
                total += 1;

                // Skip trivial cases (both constant)
                if h_sigma.len() <= 1 && h_tau.len() <= 1 {
                    pass += 1;
                    continue;
                }

                let result = check_weak_interlacing(h_tau, h_sigma);
                match result {
                    Some(true) => {
                        pass += 1;
                    }
                    _ => {
                        fail_details.push((sigma.clone(), tau.clone()));
                        println!(
                            "  FAIL: {} ⋗ {}:  h*({}) = {}  vs  h*({}) = {}",
                            perm_str(sigma),
                            perm_str(tau),
                            perm_str(sigma),
                            format_poly(h_sigma),
                            perm_str(tau),
                            format_poly(h_tau),
                        );
                    }
                }
            }
        }

        println!(
            "\nResult: {}/{} covers interlace{}",
            pass,
            total,
            if fail_details.is_empty() {
                " — ALL PASS".to_string()
            } else {
                format!(" — {} failures", fail_details.len())
            }
        );

        // Also check: does h*(τ) ≪ h*(σ) for ALL σ > τ (not just covers)?
        println!("\nAll Bruhat-comparable pairs (σ > τ, checking h*(τ) ≪ h*(σ)):");
        let mut all_total = 0;
        let mut all_pass = 0;
        let mut all_fail = 0;

        for sigma in &perms {
            for tau in &perms {
                if sigma == tau {
                    continue;
                }
                // Check σ > τ: τ is in the Bruhat lower ideal of σ
                // Quick check: inv(σ) > inv(τ) is necessary
                if inv_count(sigma) <= inv_count(tau) {
                    continue;
                }
                // Use the tableau criterion or just generate the ideal
                // For small n, just check membership
                // Actually, for efficiency, precompute ideals
            }
        }

        // Precompute: for each perm, its Bruhat lower ideal
        let ideal_map: BTreeMap<Vec<u8>, Vec<Vec<u8>>> = perms
            .iter()
            .map(|p| (p.clone(), combpoly::order::bruhat_lower_ideal(p)))
            .collect();

        for sigma in &perms {
            let ideal = &ideal_map[sigma];
            for tau in ideal {
                if tau == sigma {
                    continue;
                }
                let h_sigma = &hstar_map[sigma];
                let h_tau = &hstar_map[tau];
                all_total += 1;

                if h_sigma.len() <= 1 && h_tau.len() <= 1 {
                    all_pass += 1;
                    continue;
                }

                let result = check_weak_interlacing(h_tau, h_sigma);
                match result {
                    Some(true) => all_pass += 1,
                    _ => {
                        all_fail += 1;
                        if all_fail <= 5 {
                            println!(
                                "  FAIL: {} > {}:  h*({}) = {}  vs  h*({}) = {}",
                                perm_str(sigma),
                                perm_str(tau),
                                perm_str(sigma),
                                format_poly(h_sigma),
                                perm_str(tau),
                                format_poly(h_tau),
                            );
                        }
                    }
                }
            }
        }

        println!(
            "Result: {}/{} comparable pairs interlace{}",
            all_pass,
            all_total,
            if all_fail == 0 {
                " — ALL PASS".to_string()
            } else {
                format!(" — {} failures", all_fail)
            }
        );

        println!();
    }
}
