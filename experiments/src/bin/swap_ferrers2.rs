use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, Stat};
use polynomial_tools::real_rootedness::{
    check_interlacing_sturm, check_weak_interlacing, format_poly, is_real_rooted,
};
fn perms_on_board(lambda: &[u8]) -> Vec<Vec<u8>> {
    let n = lambda.len() as u8;
    if n == 0 {
        return vec![vec![]];
    }
    all_permutations(n)
        .into_iter()
        .filter(|s| s.iter().enumerate().all(|(i, &v)| v <= lambda[i]))
        .collect()
}
fn build_poly(perms: &[Vec<u8>]) -> Vec<i64> {
    if perms.is_empty() {
        return vec![0];
    }
    let max_s = perms
        .iter()
        .map(|s| compute(s, Stat::Swaps))
        .max()
        .unwrap_or(0);
    let mut coeffs = vec![0i64; max_s + 1];
    for s in perms {
        coeffs[compute(s, Stat::Swaps)] += 1;
    }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
        coeffs.pop();
    }
    coeffs
}
fn check_il(f: &[i64], g: &[i64]) -> &'static str {
    if f.len() <= 1 || g.len() <= 1 {
        return "triv";
    }
    let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
    match check_interlacing_sturm(small, large) {
        Some(true) => "✓",
        Some(false) => match check_weak_interlacing(small, large) {
            Some(true) => "✓w",
            _ => "✗",
        },
        None => "?",
    }
}
fn main() {
    // For n=5, enumerate ALL Ferrers boards λ with λ_i >= i (weakly increasing)
    // and check real-rootedness. Then for containment pairs, check interlacing.
    let n = 5u8;
    let mut boards: Vec<Vec<u8>> = Vec::new();
    // Generate all weakly increasing λ with λ_i >= i and λ_i <= n
    fn gen(n: u8, pos: usize, min_val: u8, current: &mut Vec<u8>, result: &mut Vec<Vec<u8>>) {
        if pos == n as usize {
            result.push(current.clone());
            return;
        }
        let lo = min_val.max(pos as u8 + 1);
        for v in lo..=n {
            current.push(v);
            gen(n, pos + 1, v, current, result);
            current.pop();
        }
    }
    gen(n, 0, 1, &mut Vec::new(), &mut boards);

    println!("=== All Ferrers boards for n={} with λ_i >= i ===\n", n);
    let mut all_rr = true;
    let mut polys: Vec<(Vec<u8>, Vec<i64>)> = Vec::new();
    for b in &boards {
        let perms = perms_on_board(b);
        let poly = build_poly(&perms);
        let rr = if poly.len() <= 2 {
            true
        } else {
            is_real_rooted(&poly)
        };
        if !rr {
            all_rr = false;
        }
        let name: String = b
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "{:<12} |S|={:<5} {:<50} {}",
            name,
            perms.len(),
            format_poly(&poly),
            if rr { "✓" } else { "✗" }
        );
        polys.push((b.clone(), poly));
    }
    println!("\nAll real-rooted: {}", if all_rr { "YES" } else { "NO" });

    // Check containment interlacing: if λ ⊆ μ (entrywise), does L_λ interlace L_μ?
    println!("\n=== Containment interlacing failures ===");
    let mut fail_count = 0;
    let mut check_count = 0;
    for i in 0..polys.len() {
        for j in (i + 1)..polys.len() {
            let (ref la, ref pa) = polys[i];
            let (ref lb, ref pb) = polys[j];
            // Check if la ⊆ lb (entrywise)
            if la.iter().zip(lb.iter()).all(|(&a, &b)| a <= b) {
                if pa.len() > 1 && pb.len() > 1 {
                    check_count += 1;
                    let r = check_il(pa, pb);
                    if r == "✗" {
                        fail_count += 1;
                        let na: String = la
                            .iter()
                            .map(|v| v.to_string())
                            .collect::<Vec<_>>()
                            .join(",");
                        let nb: String = lb
                            .iter()
                            .map(|v| v.to_string())
                            .collect::<Vec<_>>()
                            .join(",");
                        println!("  {} ≪ {}: ✗", na, nb);
                    }
                }
            }
        }
    }
    println!(
        "Checked {} containment pairs, {} failures",
        check_count, fail_count
    );
}
