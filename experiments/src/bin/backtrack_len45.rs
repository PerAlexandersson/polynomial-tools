/// Focused study: length-4 class representatives through n=9,
/// and length-5 equivalence class survey through n=7.
use combpoly::permutation::{all_permutations, backtrack_stat_poly, contains_pattern};
use combpoly::statistics::Stat;
use polynomial_tools::real_rootedness as polynomial;
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};
use std::collections::BTreeMap;
use std::io::{self, Write};

fn pat_label(pat: &[u8]) -> String {
    pat.iter().map(|x| x.to_string()).collect()
}

fn compute_polys(pat: &[u8], stat: Stat, max_n: u8, all_perms: &[Vec<Vec<u8>>]) -> Vec<Vec<i64>> {
    let mut polys = Vec::new();
    for n in 1..=max_n {
        let perms = &all_perms[(n - 1) as usize];
        let s_set: Vec<Vec<u8>> = perms
            .iter()
            .filter(|p| !contains_pattern(p, pat))
            .cloned()
            .collect();
        polys.push(backtrack_stat_poly(n, &s_set, stat));
    }
    polys
}

fn print_triangle(polys: &[Vec<i64>]) {
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
}

fn search_rec(polys: &[Vec<i64>], label: &str) {
    // Try progressively more complex searches
    let searches = [
        (
            "simple",
            AdaptiveSearchOptions {
                max_rec_len: 2,
                max_var_deg: 2,
                max_idx_deg: 2,
                max_diff_deg: 1,
                verbose: false,
                ..Default::default()
            },
        ),
        (
            "medium",
            AdaptiveSearchOptions {
                max_rec_len: 3,
                max_var_deg: 2,
                max_idx_deg: 2,
                max_diff_deg: 1,
                verbose: false,
                ..Default::default()
            },
        ),
        (
            "wide",
            AdaptiveSearchOptions {
                max_rec_len: 3,
                max_var_deg: 3,
                max_idx_deg: 3,
                max_diff_deg: 1,
                verbose: false,
                ..Default::default()
            },
        ),
    ];

    for (level, opts) in &searches {
        eprint!("\r    [{}] searching...          ", level);
        io::stderr().flush().ok();
        if let Some(result) = find_recurrence_adaptive(polys, opts) {
            println!("    Recurrence [{}]: {}", level, result.recurrence);
            return;
        }
    }
    println!("    Recurrence: none found");
}

fn main() {
    // =========================================================
    // Part 1: Length-4 class representatives through n=9
    // =========================================================
    println!("{}", "=".repeat(72));
    println!("  PART 1: Length-4 class representatives (n=1..9)");
    println!("{}", "=".repeat(72));

    let max_n_4: u8 = 9;
    eprintln!("Precomputing permutations for n=1..{} ...", max_n_4);
    let all_perms_9: Vec<Vec<Vec<u8>>> = (1..=max_n_4).map(|n| all_permutations(n)).collect();

    // Des class representatives (from the equivalence classes we found)
    let des_reps: Vec<(&str, Vec<u8>)> = vec![
        ("1324 (Eulerian, 12 patterns)", vec![1, 3, 2, 4]),
        ("1234", vec![1, 2, 3, 4]),
        ("4321", vec![4, 3, 2, 1]),
        ("1342=2341", vec![1, 3, 4, 2]),
        ("1432=2431", vec![1, 4, 3, 2]),
        ("2134=3412", vec![2, 1, 3, 4]),
        ("2143=4312", vec![2, 1, 4, 3]),
        ("3421", vec![3, 4, 2, 1]),
        ("1243", vec![1, 2, 4, 3]),
    ];

    // Peak class representatives
    let peak_reps: Vec<(&str, Vec<u8>)> = vec![
        ("1324 (standard peak, 12 pats)", vec![1, 3, 2, 4]),
        ("1234=2134=4312=4321 (4 pats)", vec![1, 2, 3, 4]),
        ("1342=1432=2341=2431 (4 pats)", vec![1, 3, 4, 2]),
        ("1243=2143=3412=3421 (4 pats)", vec![1, 2, 4, 3]),
    ];

    println!("\n--- DESCENTS ---\n");
    for (label, pat) in &des_reps {
        eprint!("\r  Computing Av_{} + des ...          ", label);
        io::stderr().flush().ok();
        let polys = compute_polys(pat, Stat::Des, max_n_4, &all_perms_9);
        let rr = polys
            .iter()
            .enumerate()
            .find(|(_, p)| !polynomial::is_real_rooted(p))
            .map(|(i, _)| format!("FAILS n={}", i + 1))
            .unwrap_or_else(|| format!("yes (n≤{})", max_n_4));
        println!("\n  Av_{} + des    RR: {}", label, rr);
        print_triangle(&polys);
        search_rec(&polys, label);
    }

    println!("\n--- PEAKS ---\n");
    for (label, pat) in &peak_reps {
        eprint!("\r  Computing Av_{} + peak ...          ", label);
        io::stderr().flush().ok();
        let polys = compute_polys(pat, Stat::Peak, max_n_4, &all_perms_9);
        let rr = polys
            .iter()
            .enumerate()
            .find(|(_, p)| !polynomial::is_real_rooted(p))
            .map(|(i, _)| format!("FAILS n={}", i + 1))
            .unwrap_or_else(|| format!("yes (n≤{})", max_n_4));
        println!("\n  Av_{} + peak    RR: {}", label, rr);
        print_triangle(&polys);
        search_rec(&polys, label);
    }

    // =========================================================
    // Part 2: Length-5 equivalence class survey (n=1..7)
    // =========================================================
    println!("\n{}", "=".repeat(72));
    println!("  PART 2: Length-5 equivalence classes (n=1..7)");
    println!("{}", "=".repeat(72));

    let max_n_5: u8 = 7;
    // Reuse all_perms_9 up to index 6 (n=7)
    let patterns5 = all_permutations(5);
    let stats: Vec<(&str, Stat)> =
        vec![("des", Stat::Des), ("peak", Stat::Peak), ("exc", Stat::Exc)];

    for &(stat_name, stat) in &stats {
        println!(
            "\n--- {} (length-5 patterns, n=1..{}) ---\n",
            stat_name, max_n_5
        );

        let mut equiv_classes: BTreeMap<String, Vec<String>> = BTreeMap::new();
        let mut class_data: BTreeMap<String, (Vec<Vec<i64>>, Option<usize>)> = BTreeMap::new();

        for (pi, pat) in patterns5.iter().enumerate() {
            let label = pat_label(pat);
            let polys = compute_polys(pat, stat, max_n_5, &all_perms_9);
            let first_fail = polys
                .iter()
                .enumerate()
                .find(|(_, p)| !polynomial::is_real_rooted(p))
                .map(|(i, _)| i + 1);

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

            equiv_classes
                .entry(key.clone())
                .or_default()
                .push(label.clone());
            class_data.entry(key).or_insert((polys, first_fail));

            if (pi + 1) % 20 == 0 {
                eprint!(
                    "\r  {}/{} patterns done...          ",
                    pi + 1,
                    patterns5.len()
                );
                io::stderr().flush().ok();
            }
        }
        eprintln!("\r  All {} patterns done.              ", patterns5.len());

        // Print summary
        println!("  {} equivalence classes:\n", equiv_classes.len());

        let mut class_id = 0;
        for (key, pats) in &equiv_classes {
            class_id += 1;
            let (polys, first_fail) = class_data.get(key).unwrap();
            let rr = match first_fail {
                None => format!("RR: yes (n≤{})", max_n_5),
                Some(n) => format!("RR: FAILS n={}", n),
            };
            let size = pats.len();
            // Show first few patterns if class is large
            let pats_display = if pats.len() <= 6 {
                pats.join(" = ")
            } else {
                format!("{} = {} = ... ({} patterns)", pats[0], pats[1], pats.len())
            };
            println!(
                "  Class {:2}: [size {:3}] Av_{{{}}}    {}",
                class_id, size, pats_display, rr
            );
        }

        // For classes with clean recurrences (small classes), search
        println!();
        let mut searched = 0;
        for (key, pats) in &equiv_classes {
            let (polys, _) = class_data.get(key).unwrap();
            if polys.len() < 4 {
                continue;
            }

            let pats_display = if pats.len() <= 4 {
                pats.join("=")
            } else {
                format!("{}=...({} pats)", pats[0], pats.len())
            };

            eprint!(
                "\r  Searching rec for Av_{{{}}} + {} ...          ",
                pats_display, stat_name
            );
            io::stderr().flush().ok();

            let opts = AdaptiveSearchOptions {
                max_rec_len: 2,
                max_var_deg: 2,
                max_idx_deg: 2,
                max_diff_deg: 1,
                verbose: false,
                ..Default::default()
            };
            if let Some(result) = find_recurrence_adaptive(polys, &opts) {
                println!(
                    "  Av_{{{}}} + {}: {}",
                    pats_display, stat_name, result.recurrence
                );
                searched += 1;
            }
        }
        if searched == 0 {
            // Try medium search on the largest class
            let largest = equiv_classes.values().max_by_key(|v| v.len()).unwrap();
            let (polys, _) = class_data
                .get(
                    &equiv_classes
                        .iter()
                        .find(|(_, v)| std::ptr::eq(v.as_slice(), largest.as_slice()))
                        .unwrap()
                        .0
                        .clone(),
                )
                .unwrap();
            let opts = AdaptiveSearchOptions {
                max_rec_len: 2,
                max_var_deg: 2,
                max_idx_deg: 2,
                max_diff_deg: 1,
                verbose: false,
                ..Default::default()
            };
            let pats_display = if largest.len() <= 4 {
                largest.join("=")
            } else {
                format!("{}=...({} pats)", largest[0], largest.len())
            };
            eprint!(
                "\r  Searching rec for largest class Av_{{{}}} + {} ...          ",
                pats_display, stat_name
            );
            io::stderr().flush().ok();
            match find_recurrence_adaptive(polys, &opts) {
                Some(r) => println!(
                    "  Av_{{{}}} + {}: {}",
                    pats_display, stat_name, r.recurrence
                ),
                None => println!("  No simple recurrence found for largest class."),
            }
        }
        eprintln!();
    }
}
