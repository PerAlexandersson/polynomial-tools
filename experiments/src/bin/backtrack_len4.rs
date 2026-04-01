/// Systematic study of backtrack polynomials for length-4 pattern avoidance.
///
/// For each τ ∈ S_4 and stat ∈ {des, exc, peak}, computes
///   P_n(t) = Σ_{π ∈ S_n} t^{stat(backtrack_π(Av_τ(n)))}
/// through n = max_n, checks real-rootedness, identifies equivalence classes,
/// and searches for recurrences.
use combpoly::permutation::{all_permutations, backtrack_stat_poly, contains_pattern};
use combpoly::statistics::Stat;
use polynomial_tools::real_rootedness as polynomial;
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};
use std::collections::BTreeMap;
use std::env;
use std::io::{self, Write};

fn pat_label(pat: &[u8]) -> String {
    pat.iter().map(|x| x.to_string()).collect()
}

fn main() {
    let max_n: u8 = env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(8);

    let run_recurrence = env::args().any(|a| a == "-r" || a == "--recurrence");
    let verbose = env::args().any(|a| a == "-v" || a == "--verbose");

    let patterns = all_permutations(4);
    let stats: Vec<(&str, Stat)> =
        vec![("des", Stat::Des), ("exc", Stat::Exc), ("peak", Stat::Peak)];

    // Precompute all permutations and avoidance sets
    eprintln!("Precomputing avoidance sets for n=1..{} ...", max_n);
    let mut all_perms: Vec<Vec<Vec<u8>>> = Vec::new();
    // all_perms[n-1] = all permutations of size n
    for n in 1..=max_n {
        all_perms.push(all_permutations(n));
    }

    for &(stat_name, stat) in &stats {
        println!("\n{}", "=".repeat(72));
        println!("  STAT = {}    (n = 1..{})", stat_name, max_n);
        println!("{}", "=".repeat(72));

        // For each pattern, compute the polynomial sequence
        // Key: polynomial sequence as string -> list of patterns
        let mut equiv_classes: BTreeMap<String, Vec<String>> = BTreeMap::new();
        let mut pat_data: Vec<(String, Vec<Vec<i64>>, Option<usize>)> = Vec::new();

        for pat in &patterns {
            let label = pat_label(pat);
            let mut polys: Vec<Vec<i64>> = Vec::new();

            for n in 1..=max_n {
                let perms = &all_perms[(n - 1) as usize];
                let s_set: Vec<Vec<u8>> = perms
                    .iter()
                    .filter(|p| !contains_pattern(p, pat))
                    .cloned()
                    .collect();

                let coeffs = backtrack_stat_poly(n, &s_set, stat);
                polys.push(coeffs);
            }
            eprint!("\r  Av_{} + {} done          ", label, stat_name);
            io::stderr().flush().ok();

            // Check real-rootedness
            let first_fail = polys
                .iter()
                .enumerate()
                .find(|(_, p)| !polynomial::is_real_rooted(p))
                .map(|(i, _)| i + 1);

            // Build equivalence class key from coefficient arrays
            let key: String = polys
                .iter()
                .map(|p| {
                    p.iter()
                        .map(|c| c.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .collect::<Vec<_>>()
                .join("|");

            equiv_classes.entry(key).or_default().push(label.clone());

            pat_data.push((label, polys, first_fail));
        }
        eprintln!();

        // Print equivalence classes
        println!("\n--- Equivalence classes (patterns giving identical arrays) ---\n");
        let mut class_id = 0;
        let mut printed_pats = std::collections::BTreeSet::new();

        for (_key, pats) in &equiv_classes {
            class_id += 1;
            let rr_status = pat_data
                .iter()
                .find(|(l, _, _)| l == &pats[0])
                .map(|(_, _, f)| match f {
                    None => format!("RR: yes (n≤{})", max_n),
                    Some(n) => format!("RR: FAILS at n={}", n),
                })
                .unwrap_or_default();

            let pats_str = pats.join(" = ");
            println!("  Class {}: Av_{{{}}}    {}", class_id, pats_str, rr_status);
            for p in pats {
                printed_pats.insert(p.clone());
            }
        }

        // Print triangles for each equivalence class representative
        if verbose {
            println!("\n--- Triangles (one per class) ---\n");
            for (_key, pats) in &equiv_classes {
                let rep = &pats[0];
                let (_, polys, _) = pat_data.iter().find(|(l, _, _)| l == rep).unwrap();
                println!("  Av_{} + {}:", rep, stat_name);
                for (i, p) in polys.iter().enumerate() {
                    print!("    n={}: ", i + 1);
                    for (j, &c) in p.iter().enumerate() {
                        if j > 0 {
                            print!(" ");
                        }
                        print!("{}", c);
                    }
                    println!();
                }
                println!();
            }
        }

        // Print real-rootedness summary table
        println!("\n--- Real-rootedness summary ---\n");
        println!("  {:>6}  {:>12}", "Av_τ", "Real-rooted?");
        println!("  {:->6}  {:->12}", "", "");
        for (label, _, first_fail) in &pat_data {
            let status = match first_fail {
                None => format!("yes (n≤{})", max_n),
                Some(n) => format!("FAILS n={}", n),
            };
            println!("  {:>6}  {:>12}", label, status);
        }

        // Recurrence search
        if run_recurrence {
            println!("\n--- Recurrence search (one per equivalence class) ---\n");
            for (_key, pats) in &equiv_classes {
                let rep = &pats[0];
                let (_, polys, _) = pat_data.iter().find(|(l, _, _)| l == rep).unwrap();

                if polys.len() < 4 {
                    continue;
                }

                let search = AdaptiveSearchOptions {
                    max_rec_len: 3,
                    max_var_deg: 3,
                    max_idx_deg: 3,
                    max_diff_deg: 1,
                    verbose: false,
                    ..Default::default()
                };

                let class_label = if pats.len() > 1 {
                    format!("Av_{{{}}}", pats.join("="))
                } else {
                    format!("Av_{}", rep)
                };

                eprint!(
                    "\r  Searching recurrence for {} + {} ...         ",
                    class_label, stat_name
                );
                io::stderr().flush().ok();

                match find_recurrence_adaptive(polys, &search) {
                    Some(result) => {
                        println!("  {} + {}: {}", class_label, stat_name, result.recurrence);
                    }
                    None => {
                        println!("  {} + {}: (none found)", class_label, stat_name);
                    }
                }
            }
            eprintln!();
        }
    }
}
