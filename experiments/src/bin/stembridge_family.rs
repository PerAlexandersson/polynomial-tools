// Check permutations similar to Stembridge's: [2,4,6,...,2k,1,2k+2,3,...,n]
// These are 321-avoiding with a specific "zigzag" structure.
use combpoly::order;
use combpoly::permutation::contains_pattern;
use polynomial_tools::real_rootedness as polynomial;
use combpoly::statistics::Stat;

fn stembridge_type(n: u8, k: u8) -> Vec<u8> {
    // First k positions: 2,4,6,...,2k (evens)
    // Then: 1, 2k+2, 3, 2k+4, 5, ...  (interleaving remaining odds and evens)
    let mut perm = Vec::new();
    for i in 1..=k { perm.push(2*i); }
    // remaining values sorted
    let used: std::collections::HashSet<u8> = perm.iter().copied().collect();
    let remaining: Vec<u8> = (1..=n).filter(|x| !used.contains(x)).collect();
    // interleave: odd-indexed remaining go ascending
    perm.extend_from_slice(&remaining);
    perm
}

fn main() {
    // Try the exact Stembridge permutation and variations
    let examples: Vec<(String, Vec<u8>)> = vec![
        ("Stembridge n=17".into(), vec![2,4,6,8,10,1,12,3,15,5,17,7,9,11,13,14,16]),
    ];

    for (label, w) in &examples {
        let n = w.len();
        println!("=== {} ===", label);
        println!("w = {:?}", w);
        println!("avoids 321: {}", !contains_pattern(w, &[3,2,1]));

        let coeffs = order::weak_ideal_poly(w, Stat::Des);
        let sz: i64 = coeffs.iter().sum();
        println!("ideal={}, deg={}", sz, coeffs.len()-1);
        println!("poly: {}", polynomial::format_poly(&coeffs));
        println!("real-rooted: {}", polynomial::is_real_rooted(&coeffs));
        println!("log-concave: {}", polynomial::is_log_concave(&coeffs));
        println!("ULC: {}", polynomial::is_ultra_log_concave(&coeffs));
        println!();
    }

    // Now try truncated/extended versions to find smallest counterexample
    println!("=== Searching for smallest counterexample in Av_321 ===");
    for n in 10..=17u8 {
        // Try many 321-avoiding permutations with large inversions
        let perms = combpoly::permutation::avoiding_permutations(n, &[vec![3,2,1]]);
        let mut worst = None;
        let mut worst_inv = 0;
        for w in &perms {
            let inv: usize = (0..w.len()).flat_map(|i| (i+1..w.len()).filter(move |&j| w[i] > w[j])).count();
            if inv > worst_inv {
                worst_inv = inv;
                worst = Some(w.clone());
            }
        }
        // Check permutations with most inversions (most likely to fail RR)
        let mut fail_count = 0;
        // Check top 1% by inversions
        let mut inv_perms: Vec<(usize, &Vec<u8>)> = perms.iter()
            .map(|w| {
                let inv: usize = (0..w.len()).flat_map(|i| (i+1..w.len()).filter(move |&j| w[i] > w[j])).count();
                (inv, w)
            }).collect();
        inv_perms.sort_by(|a,b| b.0.cmp(&a.0));
        let top = std::cmp::max(1, inv_perms.len() / 100);
        for (inv, w) in &inv_perms[..top] {
            let coeffs = order::weak_ideal_poly(w, Stat::Des);
            if !polynomial::is_real_rooted(&coeffs) {
                fail_count += 1;
                if fail_count == 1 {
                    println!("n={}: FAIL! w={:?}, inv={}, coeffs={:?}", n, w, inv, coeffs);
                }
            }
        }
        if fail_count == 0 {
            println!("n={}: top {}  perms (by inv) all pass RR (max inv={})", n, top, worst_inv);
        } else {
            println!("n={}: {} failures in top {} perms", n, fail_count, top);
        }
    }
}
