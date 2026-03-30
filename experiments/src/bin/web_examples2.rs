use combpoly::permutation::all_permutations;
use combpoly::polynomial_builder::build_generating_polynomial;
use combpoly::statistics::Stat;
use polynomial_tools::sequences::*;
use polynomial_tools::real_rootedness::format_poly;

fn trim(p: &[i64]) -> Vec<i64> {
    let end = p.iter().rposition(|&c| c != 0).map_or(0, |i| i + 1);
    if end == 0 { vec![0] } else { p[..end].to_vec() }
}

fn main() {
    println!("=== Descent on derangements d_n(t) ===");
    for n in 2..=14u8 {
        let perms = all_permutations(n);
        let derangements: Vec<Vec<u8>> = perms
            .into_iter()
            .filter(|w| w.iter().enumerate().all(|(i, &v)| v != (i as u8 + 1)))
            .collect();
        let poly = build_generating_polynomial(&derangements, Stat::Des);
        let t = trim(&poly);
        // Print as comma-separated coefficients for easy copy
        let coeffs: Vec<String> = t.iter().map(|c| c.to_string()).collect();
        println!("# d_{}(t)", n);
        println!("{}", coeffs.join(", "));
    }

    println!("\n=== Eulerian A_n(t) ===");
    let ep = eulerian_polynomials(14);
    for (i, p) in ep.iter().enumerate() {
        let t = trim(p);
        let coeffs: Vec<String> = t.iter().map(|c| c.to_string()).collect();
        println!("{}", coeffs.join(", "));
    }

    println!("\n=== Narayana N_n(t) ===");
    let np = narayana_polynomials(14);
    for (i, p) in np.iter().enumerate() {
        let t = trim(p);
        let coeffs: Vec<String> = t.iter().map(|c| c.to_string()).collect();
        println!("{}", coeffs.join(", "));
    }
}
