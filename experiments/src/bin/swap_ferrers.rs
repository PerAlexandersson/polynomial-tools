/// Long swaps on permutations fitting a Ferrers board.
///
/// A permutation σ fits a Ferrers board λ = (λ₁ ≤ λ₂ ≤ ... ≤ λ_n)
/// if σ(i) ≤ λ_i for all i (French convention, weakly increasing λ).
///
/// We compute the long-swaps generating polynomial and check
/// real-rootedness and interlacing across dilations kλ.
use combpoly::permutation::all_permutations;
use combpoly::statistics::{compute, Stat};
use polynomial_tools::real_rootedness::{
    check_interlacing_sturm, check_weak_interlacing, format_poly, is_real_rooted,
};

/// Enumerate all permutations of [n] fitting board λ.
/// λ is given in weakly increasing order (French convention):
/// row i allows columns 1..=λ[i].
fn perms_on_board(lambda: &[u8]) -> Vec<Vec<u8>> {
    let n = lambda.len() as u8;
    if n == 0 {
        return vec![vec![]];
    }
    let all = all_permutations(n);
    all.into_iter()
        .filter(|sigma| sigma.iter().enumerate().all(|(i, &v)| v <= lambda[i]))
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
    // Test various Ferrers boards
    // French convention: λ = (λ₁ ≤ λ₂ ≤ ... ≤ λ_n) weakly increasing
    let boards: Vec<(&str, Vec<u8>)> = vec![
        // Staircase shapes
        ("1,2,3", vec![1, 2, 3]),
        ("1,2,3,4", vec![1, 2, 3, 4]),
        ("1,2,3,4,5", vec![1, 2, 3, 4, 5]),
        ("1,2,3,4,5,6", vec![1, 2, 3, 4, 5, 6]),
        ("1,2,3,4,5,6,7", vec![1, 2, 3, 4, 5, 6, 7]),
        // Near-staircase
        ("2,3,4,5", vec![2, 3, 4, 5]),
        ("2,3,4,5,6", vec![2, 3, 4, 5, 6]),
        ("2,3,4,5,6,7", vec![2, 3, 4, 5, 6, 7]),
        // Rectangular (full S_n)
        ("3,3,3", vec![3, 3, 3]),
        ("4,4,4,4", vec![4, 4, 4, 4]),
        // Other shapes
        ("1,3,3", vec![1, 3, 3]),
        ("1,2,4,4", vec![1, 2, 4, 4]),
        ("2,2,3,4", vec![2, 2, 3, 4]),
        ("1,3,4,4", vec![1, 3, 4, 4]),
        ("2,4,4,4", vec![2, 4, 4, 4]),
        ("1,2,3,5,5", vec![1, 2, 3, 5, 5]),
        ("2,3,5,5,5", vec![2, 3, 5, 5, 5]),
        ("1,2,4,5,5", vec![1, 2, 4, 5, 5]),
        ("3,3,4,5,5", vec![3, 3, 4, 5, 5]),
    ];

    for (name, lambda) in &boards {
        let perms = perms_on_board(lambda);
        let poly = build_poly(&perms);
        let rr = if poly.len() <= 2 {
            true
        } else {
            is_real_rooted(&poly)
        };
        println!(
            "λ={:<15} |perms|={:<6} {:<55} {}",
            name,
            perms.len(),
            format_poly(&poly),
            if rr { "✓rr" } else { "✗rr" }
        );
    }

    // Now try dilation: for a fixed shape λ, consider kλ for k=1,2,3,...
    // and check interlacing of the sequence
    println!("\n=== Dilation interlacing ===\n");

    let base_shapes: Vec<(&str, Vec<u8>)> = vec![
        ("1,2,3", vec![1, 2, 3]),
        ("1,2,3,4", vec![1, 2, 3, 4]),
        ("2,3,4", vec![2, 3, 4]),
        ("1,3,3", vec![1, 3, 3]),
        ("2,2,3,4", vec![2, 2, 3, 4]),
    ];

    for (name, base) in &base_shapes {
        println!("--- λ = {} ---", name);
        let n = base.len() as u8;
        let max_k = if n <= 4 { 5 } else { 3 };

        let mut polys: Vec<Vec<i64>> = Vec::new();
        for k in 1..=max_k {
            let dilated: Vec<u8> = base
                .iter()
                .map(|&v| {
                    let d = v as u32 * k as u32;
                    d.min(n as u32) as u8 // cap at n (can't exceed n columns)
                })
                .collect();
            let perms = perms_on_board(&dilated);
            let poly = build_poly(&perms);
            let rr = if poly.len() <= 2 {
                true
            } else {
                is_real_rooted(&poly)
            };
            let dname: String = dilated
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(",");
            println!(
                "  k={}: λ={:<12} |perms|={:<6} {:<50} {}",
                k,
                dname,
                perms.len(),
                format_poly(&poly),
                if rr { "✓" } else { "✗" }
            );
            polys.push(poly);
        }

        // Check interlacing
        for i in 0..polys.len() - 1 {
            println!(
                "    k={} ≪ k={}: {}",
                i + 1,
                i + 2,
                check_il(&polys[i], &polys[i + 1])
            );
        }
        println!();
    }

    // Also: for staircase λ=(1,2,...,n), permutations fitting are exactly Av(312)!
    // Check this and compare with the earlier Av(312) results.
    println!("=== Staircase = Av(312) check ===\n");
    for n in 3..=8u8 {
        let staircase: Vec<u8> = (1..=n).collect();
        let board_perms = perms_on_board(&staircase);
        let av312 = combpoly::permutation::avoiding_permutations(n, &[vec![3, 1, 2]]);
        let bp = build_poly(&board_perms);
        let ap = build_poly(&av312);
        println!(
            "n={}: staircase={}, Av(312)={}, match={}",
            n,
            board_perms.len(),
            av312.len(),
            bp == ap
        );
    }
}
