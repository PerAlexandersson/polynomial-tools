/// Search for palindromic backtrack polynomial sequences.
///
/// Systematically checks all (constraint set, statistic) pairs for
/// palindromic (symmetric) coefficient sequences, filtering out
/// known Eulerian and Narayana cases.
use combpoly::permutation::{
    all_permutations, avoiding_permutations, backtrack_image, backtrack_stat_poly,
    backtrack_stat_poly_subset, contains_pattern, is_alternating, is_derangement, is_grassmann,
    is_involution, is_separable, is_simsun, is_vexillary, parse_sequence,
};
use combpoly::statistics::{compute, Stat};
use std::env;

fn is_palindromic(p: &[i64]) -> bool {
    let n = p.len();
    if n <= 1 {
        return true;
    }
    // Trim trailing zeros
    let end = p.iter().rposition(|&x| x != 0).unwrap_or(0);
    let start = p.iter().position(|&x| x != 0).unwrap_or(0);
    let trimmed = &p[start..=end];
    let m = trimmed.len();
    (0..m / 2).all(|i| trimmed[i] == trimmed[m - 1 - i])
}

fn format_poly(p: &[i64]) -> String {
    p.iter()
        .map(|c| c.to_string())
        .collect::<Vec<_>>()
        .join(", ")
}

fn main() {
    let max_n: u8 = env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(8);

    let patterns = ["123", "132", "213", "231", "312", "321"];
    let stats: &[(&str, Stat)] = &[
        ("des", Stat::Des),
        ("exc", Stat::Exc),
        ("peak", Stat::Peak),
        ("lrmax", Stat::Lrmax),
        ("inv", Stat::Inv),
        ("maj", Stat::Maj),
        ("fix", Stat::Fix),
        ("cyc", Stat::Cyc),
    ];

    println!(
        "=== Palindromic search: pattern-avoiding S-sets (n=1..{}) ===\n",
        max_n
    );

    for &(stat_name, stat) in stats {
        for pat_str in &patterns {
            let pat = parse_sequence(pat_str);
            let mut polys: Vec<Vec<i64>> = Vec::new();
            let mut all_palindromic = true;

            for n in 1..=max_n {
                let s_set = avoiding_permutations(n, &[pat.clone()]);
                let coeffs = backtrack_stat_poly(n, &s_set, stat);
                if !is_palindromic(&coeffs) {
                    all_palindromic = false;
                    break;
                }
                polys.push(coeffs);
            }

            if all_palindromic && polys.len() >= 4 {
                // Check it's not trivially constant
                let nontrivial = polys.iter().any(|p| p.len() > 2);
                if nontrivial {
                    println!(
                        "PALINDROMIC: Av_{} + {} (n=1..{})",
                        pat_str, stat_name, max_n
                    );
                    for (i, p) in polys.iter().enumerate() {
                        println!("  n={}: [{}]", i + 1, format_poly(p));
                    }
                    println!();
                }
            }
        }
    }

    // Named constraint sets
    println!(
        "\n=== Palindromic search: named S-sets (n=1..{}) ===\n",
        max_n
    );

    let named_sets: Vec<(&str, Box<dyn Fn(&[u8]) -> bool>)> = vec![
        ("involutions", Box::new(|p: &[u8]| is_involution(p))),
        ("derangements", Box::new(|p: &[u8]| is_derangement(p))),
        ("alternating", Box::new(|p: &[u8]| is_alternating(p))),
        ("vexillary", Box::new(|p: &[u8]| is_vexillary(p))),
        ("grassmann", Box::new(|p: &[u8]| is_grassmann(p))),
        ("separable", Box::new(|p: &[u8]| is_separable(p))),
        ("simsun", Box::new(|p: &[u8]| is_simsun(p))),
    ];

    for &(stat_name, stat) in stats {
        for (set_name, filter_fn) in &named_sets {
            let mut polys: Vec<Vec<i64>> = Vec::new();
            let mut all_palindromic = true;

            for n in 1..=max_n {
                let all = all_permutations(n);
                let s_set: Vec<Vec<u8>> = all.into_iter().filter(|p| filter_fn(p)).collect();
                if s_set.is_empty() {
                    polys.push(vec![0]);
                    continue;
                }
                let coeffs = backtrack_stat_poly(n, &s_set, stat);
                if !is_palindromic(&coeffs) {
                    all_palindromic = false;
                    break;
                }
                polys.push(coeffs);
            }

            if all_palindromic && polys.len() >= 4 {
                let nontrivial = polys.iter().any(|p| p.len() > 2);
                if nontrivial {
                    println!("PALINDROMIC: {} + {} (n=1..{})", set_name, stat_name, max_n);
                    for (i, p) in polys.iter().enumerate() {
                        println!("  n={}: [{}]", i + 1, format_poly(p));
                    }
                    println!();
                }
            }
        }
    }

    // Restricted search orders (Catalan): Av_pi x Av_S + stat
    println!(
        "\n=== Palindromic search: Catalan search orders (n=1..{}) ===\n",
        max_n
    );

    for &(stat_name, stat) in stats {
        for pi_pat_str in &patterns {
            let pi_p = parse_sequence(pi_pat_str);
            for s_pat_str in &patterns {
                let s_p = parse_sequence(s_pat_str);
                let mut polys: Vec<Vec<i64>> = Vec::new();
                let mut all_palindromic = true;

                for n in 1..=max_n {
                    let all = all_permutations(n);
                    let pi_set: Vec<Vec<u8>> = all
                        .iter()
                        .filter(|p| !contains_pattern(p, &pi_p))
                        .cloned()
                        .collect();
                    let s_set: Vec<Vec<u8>> = all
                        .into_iter()
                        .filter(|p| !contains_pattern(p, &s_p))
                        .collect();
                    let coeffs = backtrack_stat_poly_subset(&pi_set, &s_set, stat);
                    if !is_palindromic(&coeffs) {
                        all_palindromic = false;
                        break;
                    }
                    polys.push(coeffs);
                }

                if all_palindromic && polys.len() >= 4 {
                    let nontrivial = polys.iter().any(|p| p.len() > 2);
                    if nontrivial {
                        println!(
                            "PALINDROMIC: Av_{} pi + Av_{} S + {} (n=1..{})",
                            pi_pat_str, s_pat_str, stat_name, max_n
                        );
                        for (i, p) in polys.iter().enumerate() {
                            println!("  n={}: [{}]", i + 1, format_poly(p));
                        }
                        println!();
                    }
                }
            }
        }
    }

    // Bonus: try two-pattern avoidance sets
    println!(
        "\n=== Palindromic search: two-pattern avoidance (n=1..{}) ===\n",
        max_n
    );

    for &(stat_name, stat) in stats {
        for i in 0..patterns.len() {
            for j in i + 1..patterns.len() {
                let pat1 = parse_sequence(patterns[i]);
                let pat2 = parse_sequence(patterns[j]);
                let mut polys: Vec<Vec<i64>> = Vec::new();
                let mut all_palindromic = true;

                for n in 1..=max_n {
                    let s_set = avoiding_permutations(n, &[pat1.clone(), pat2.clone()]);
                    if s_set.is_empty() {
                        polys.push(vec![0]);
                        continue;
                    }
                    let coeffs = backtrack_stat_poly(n, &s_set, stat);
                    if !is_palindromic(&coeffs) {
                        all_palindromic = false;
                        break;
                    }
                    polys.push(coeffs);
                }

                if all_palindromic && polys.len() >= 4 {
                    let nontrivial = polys.iter().any(|p| p.len() > 2);
                    if nontrivial {
                        println!(
                            "PALINDROMIC: Av_{{{},{}}} + {} (n=1..{})",
                            patterns[i], patterns[j], stat_name, max_n
                        );
                        for (i, p) in polys.iter().enumerate() {
                            println!("  n={}: [{}]", i + 1, format_poly(p));
                        }
                        println!();
                    }
                }
            }
        }
    }
}
