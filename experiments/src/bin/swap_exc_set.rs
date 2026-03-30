/// Long swaps on permutations with a fixed excedance set.
///
/// For S ⊆ [n], let E(n, S) = {σ ∈ S_n : Exc(σ) = S}
/// where Exc(σ) = {i : σ(i) > i}.
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, Stat};
use polynomial_tools::real_rootedness::{format_poly, is_real_rooted};

fn excedance_set(sigma: &[u8]) -> u32 {
    let mut s = 0u32;
    for (i, &v) in sigma.iter().enumerate() {
        if v > (i as u8 + 1) {
            s |= 1 << i;
        }
    }
    s
}

fn set_str(s: u32, n: u8) -> String {
    let positions: Vec<String> = (0..n)
        .filter(|&i| s & (1 << i) != 0)
        .map(|i| (i + 1).to_string())
        .collect();
    if positions.is_empty() {
        "∅".to_string()
    } else {
        format!("{{{}}}", positions.join(","))
    }
}

fn build_poly(perms: &[&Vec<u8>]) -> Vec<i64> {
    if perms.is_empty() { return vec![0]; }
    let max_s = perms.iter().map(|s| compute(s, Stat::Swaps)).max().unwrap_or(0);
    let mut coeffs = vec![0i64; max_s + 1];
    for s in perms { coeffs[compute(s, Stat::Swaps)] += 1; }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 { coeffs.pop(); }
    coeffs
}

fn main() {
    let max_n: u8 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(9);

    for n in 3..=max_n {
        let all = all_permutations(n);

        let mut by_exc: std::collections::BTreeMap<u32, Vec<&Vec<u8>>> =
            std::collections::BTreeMap::new();
        for s in &all {
            by_exc.entry(excedance_set(s)).or_default().push(s);
        }

        let total_sets = by_exc.len();
        let mut rr_count = 0;
        let mut fail_count = 0;
        let mut failures = Vec::new();

        for (&es, perms) in &by_exc {
            let poly = build_poly(perms);
            let rr = if poly.len() <= 2 { true } else { is_real_rooted(&poly) };
            if rr {
                rr_count += 1;
            } else {
                fail_count += 1;
                failures.push((es, perms.len(), poly));
            }
        }

        println!(
            "n={}: {} excedance sets, {} real-rooted, {} failures",
            n, total_sets, rr_count, fail_count
        );

        for (es, count, poly) in failures.iter().take(8) {
            println!(
                "  FAIL: Exc={} |E|={} {}",
                set_str(*es, n), count, format_poly(poly)
            );
        }
        if fail_count > 8 {
            println!("  ... and {} more failures", fail_count - 8);
        }
    }
}
