/// Fast verification of Catalan search order polynomials.
///
/// Since backtrack_π(S_n) = π⁻¹, we just compute stat on Av_{τ⁻¹}(n).
/// This avoids the expensive backtrack computation and lets us go to n=12+.
use combpoly::permutation::{all_permutations, contains_pattern};
use combpoly::statistics::{self, Stat};
use polynomial_tools::real_rootedness as polynomial;
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};
use std::env;
use std::io::{self, Write};

/// Compute the inverse of a permutation (1-indexed internally)
fn inverse_perm(p: &[u8]) -> Vec<u8> {
    let n = p.len();
    let mut inv = vec![0u8; n];
    for i in 0..n {
        inv[(p[i] - 1) as usize] = (i + 1) as u8;
    }
    inv
}

fn main() {
    let max_n: u8 = env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(11);

    // (search order pattern τ, target stat class = Av_{τ⁻¹})
    // 231⁻¹ = 312, 312⁻¹ = 231, 132⁻¹ = 132, 213⁻¹ = 213, 321⁻¹ = 321
    let cases: Vec<(&str, Vec<u8>, &str)> = vec![
        // Narayana cases (des)
        ("Av_231→des = des on Av_312", vec![3, 1, 2], "des"),
        ("Av_132→des = des on Av_132", vec![1, 3, 2], "des"),
        ("Av_321→des = des on Av_321", vec![3, 2, 1], "des"),
        // Exc cases
        ("Av_231→exc = exc on Av_312", vec![3, 1, 2], "exc"),
        ("Av_312→exc = exc on Av_231", vec![2, 3, 1], "exc"),
        ("Av_321→exc = exc on Av_321", vec![3, 2, 1], "exc"),
        // Peak cases
        (
            "Av_132→peak = peak on Av_132 (ClassA)",
            vec![1, 3, 2],
            "peak",
        ),
        (
            "Av_231→peak = peak on Av_312 (ClassB)",
            vec![3, 1, 2],
            "peak",
        ),
        // Valley
        ("Av_231→valley = valley on Av_312", vec![3, 1, 2], "valley"),
    ];

    let stat_map = |s: &str| -> Stat {
        match s {
            "des" => Stat::Des,
            "exc" => Stat::Exc,
            "peak" => Stat::Peak,
            "valley" => Stat::Valley,
            _ => panic!("unknown stat"),
        }
    };

    for (label, avoid_pat, stat_name) in &cases {
        let stat = stat_map(stat_name);
        println!("{}", "=".repeat(60));
        println!("  {}", label);
        println!("{}", "=".repeat(60));

        let mut polys: Vec<Vec<i64>> = Vec::new();

        for n in 1..=max_n {
            let all = all_permutations(n);
            // Compute stat on Av_{τ⁻¹}(n) = permutations avoiding avoid_pat
            let avoiding: Vec<&Vec<u8>> = all
                .iter()
                .filter(|p| !contains_pattern(p, avoid_pat))
                .collect();

            let mut coeffs: Vec<i64> = Vec::new();
            for sigma in &avoiding {
                let val = statistics::compute(sigma, stat);
                while coeffs.len() <= val {
                    coeffs.push(0);
                }
                coeffs[val] += 1;
            }
            if coeffs.is_empty() {
                coeffs.push(0);
            }

            let rr = polynomial::is_real_rooted(&coeffs);
            let sum: i64 = coeffs.iter().sum();

            print!("  n={:>2}: [", n);
            for (j, &c) in coeffs.iter().enumerate() {
                if j > 0 {
                    print!(", ");
                }
                print!("{}", c);
            }
            println!("]  sum={}  RR={}", sum, if rr { "yes" } else { "NO" });

            polys.push(coeffs);
            eprint!("\r  n={} done   ", n);
            io::stderr().flush().ok();
        }
        eprintln!();

        // Recurrence search: with LHS denominator f(n,t)·P(n) = ...
        println!("\n  Recurrence search (with LHS factor):");
        for (rlabel, opts) in &[
            (
                "2-term, f(n)·P(n)",
                AdaptiveSearchOptions {
                    max_rec_len: 1,
                    max_var_deg: 2,
                    max_idx_deg: 2,
                    max_diff_deg: 0,
                    try_denominator: true,
                    max_denom_var_deg: 1,
                    max_denom_idx_deg: 1,
                    verbose: false,
                    ..Default::default()
                },
            ),
            (
                "3-term, f(n)·P(n)",
                AdaptiveSearchOptions {
                    max_rec_len: 2,
                    max_var_deg: 2,
                    max_idx_deg: 2,
                    max_diff_deg: 0,
                    try_denominator: true,
                    max_denom_var_deg: 1,
                    max_denom_idx_deg: 1,
                    verbose: false,
                    ..Default::default()
                },
            ),
            (
                "3-term, f(n,t)·P(n)",
                AdaptiveSearchOptions {
                    max_rec_len: 2,
                    max_var_deg: 2,
                    max_idx_deg: 2,
                    max_diff_deg: 0,
                    try_denominator: true,
                    max_denom_var_deg: 2,
                    max_denom_idx_deg: 2,
                    verbose: false,
                    ..Default::default()
                },
            ),
            (
                "2-term+d, f(n)·P(n)",
                AdaptiveSearchOptions {
                    max_rec_len: 1,
                    max_var_deg: 2,
                    max_idx_deg: 2,
                    max_diff_deg: 1,
                    try_denominator: true,
                    max_denom_var_deg: 1,
                    max_denom_idx_deg: 1,
                    verbose: false,
                    ..Default::default()
                },
            ),
            (
                "3-term+d, f(n)·P(n)",
                AdaptiveSearchOptions {
                    max_rec_len: 2,
                    max_var_deg: 2,
                    max_idx_deg: 2,
                    max_diff_deg: 1,
                    try_denominator: true,
                    max_denom_var_deg: 1,
                    max_denom_idx_deg: 1,
                    verbose: false,
                    ..Default::default()
                },
            ),
        ] {
            eprint!("\r    searching {} ...         ", rlabel);
            io::stderr().flush().ok();
            match find_recurrence_adaptive(&polys, opts) {
                Some(r) => println!(
                    "    {}: {}  (eqs={}, unk={})",
                    rlabel, r.recurrence, r.num_equations, r.num_unknowns
                ),
                None => println!("    {}: (none found)", rlabel),
            }
        }
        eprintln!();

        println!();
    }
}
