/// Test whether colMap_k(s_μ) is Schur-positive when every part of μ
/// appears a multiple of k times.
///
/// colMap_k is defined on the monomial basis:
///   colMap_k(m_μ) = m_{(μ'/k)'}  if all parts of μ' are divisible by k,
///                 = 0             otherwise.
///
/// We test k = 2, 3 and |μ| ≤ 20.
use std::collections::BTreeMap;
use std::time::Instant;

use combinatoric_core::Partition;
use sym_poly_sym::{Basis, SymmetricFunction};

/// Apply colMap_k to a symmetric function.
/// Works by converting to monomial basis, then filtering and remapping.
fn col_map_k(f: &SymmetricFunction<i64>, k: u32) -> SymmetricFunction<i64> {
    let mono = f.to_monomial_basis();
    let mut result_terms: BTreeMap<Partition, i64> = BTreeMap::new();

    for (partition, &coeff) in mono.terms() {
        let conj = partition.conjugate_partition();
        // Check if all parts of the conjugate are divisible by k
        if conj.parts().iter().all(|&p| p % k == 0) {
            // Divide each part of the conjugate by k, then conjugate back
            let divided: Vec<u32> = conj.parts().iter().map(|&p| p / k).collect();
            let new_partition = Partition::from_sorted(divided).conjugate_partition();
            *result_terms.entry(new_partition).or_insert(0) += coeff;
        }
    }

    SymmetricFunction::from_terms(Basis::Monomial, result_terms)
}

/// Generate all partitions of n where every distinct part appears
/// with multiplicity divisible by k.
fn partitions_with_k_multiplicities(n: u32, k: u32) -> Vec<Partition> {
    Partition::all_of_size(n)
        .into_iter()
        .filter(|p| {
            let parts = p.parts();
            if parts.is_empty() {
                return true;
            }
            let mut i = 0;
            while i < parts.len() {
                let val = parts[i];
                let mut count = 0u32;
                while i < parts.len() && parts[i] == val {
                    count += 1;
                    i += 1;
                }
                if count % k != 0 {
                    return false;
                }
            }
            true
        })
        .collect()
}

/// Generate partitions of the form λ^k (each row repeated k times).
fn partitions_lambda_k(n: u32, k: u32) -> Vec<Partition> {
    if n % k != 0 {
        return vec![];
    }
    let m = n / k;
    Partition::all_of_size(m)
        .into_iter()
        .map(|lam| {
            let mut parts = Vec::new();
            for &p in lam.parts() {
                for _ in 0..k {
                    parts.push(p);
                }
            }
            Partition::from_sorted(parts)
        })
        .collect()
}

fn test_family(k: u32, max_size: u32, label: &str, gen: fn(u32, u32) -> Vec<Partition>) {
    println!("========================================");
    println!("  k = {},  |μ| ≤ {},  family: {}", k, max_size, label);
    println!("========================================");

    let mut total_tested = 0u64;
    let mut total_positive = 0u64;
    let mut failures: Vec<(Partition, String)> = Vec::new();

    for n in 1..=max_size {
        let parts = gen(n, k);
        if parts.is_empty() {
            continue;
        }

        let t0 = Instant::now();
        let count = parts.len();
        let mut ok_this_n = 0u64;

        for mu in &parts {
            let s_mu = SymmetricFunction::<i64>::schur_symmetric(mu.clone());
            let result = col_map_k(&s_mu, k);

            if result.is_zero() {
                continue;
            }

            let schur_result = result.to_schur_basis();
            total_tested += 1;

            if schur_result.positive_coefficients() {
                total_positive += 1;
                ok_this_n += 1;
            } else {
                // Extract just the negative terms for a compact summary
                let neg: Vec<String> = schur_result
                    .terms()
                    .iter()
                    .filter(|(_, &c)| c < 0)
                    .map(|(p, c)| format!("{}*s[{}]", c, p))
                    .collect();
                let summary = neg.join(", ");
                failures.push((mu.clone(), summary.clone()));
                println!("  FAIL: μ = {},  negative terms: {}", mu, summary);
            }
        }

        let elapsed = t0.elapsed();
        println!(
            "  |μ| = {:>2}: {:>4} partitions, {:>4} nonzero, all s-positive  {:>8.2?}",
            n, count, ok_this_n, elapsed
        );
    }

    println!("  ----------------------------------------");
    if failures.is_empty() {
        println!("  All {} nonzero cases are Schur-positive! ✓", total_tested);
    } else {
        println!(
            "  {}/{} Schur-positive, {} failures ✗",
            total_positive,
            total_tested,
            failures.len()
        );
        for (mu, neg) in &failures {
            println!("    μ = {}: {}", mu, neg);
        }
    }
    println!();
}

fn main() {
    let max_size: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(20);

    // Test 1: μ = λ^k (each row repeated k times)
    for k in [2u32, 3] {
        test_family(k, max_size, "λ^k (rows repeated)", partitions_lambda_k);
    }

    // Test 2: every part multiplicity divisible by k
    for k in [2u32, 3] {
        test_family(k, max_size, "mult(i) ≡ 0 mod k", |n, k| {
            partitions_with_k_multiplicities(n, k)
        });
    }
}
