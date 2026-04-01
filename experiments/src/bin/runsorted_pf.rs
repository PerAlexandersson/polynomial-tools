//! Enumerate parking functions under three "run-sorted" restrictions:
//!   1. Strict: min(R_1) < min(R_2) < ... < min(R_k)
//!   2. Weak:   min(R_1) <= min(R_2) <= ... <= min(R_k)
//!   3. Lex:    runs in lexicographic order R_1 <_lex R_2 <_lex ... <_lex R_k
//!
//! where R_1, R_2, ..., R_k are the maximal ascending runs.

use combpoly::parking::all_parking_functions;
use combpoly::statistics;
use combpoly::statistics::Stat;
use polynomial_tools::{format_poly, is_real_rooted};

/// Extract maximal ascending runs from a word.
fn ascending_runs(w: &[u8]) -> Vec<&[u8]> {
    if w.is_empty() {
        return vec![];
    }
    let mut runs = Vec::new();
    let mut start = 0;
    for i in 1..w.len() {
        if w[i] <= w[i - 1] {
            runs.push(&w[start..i]);
            start = i;
        }
    }
    runs.push(&w[start..]);
    runs
}

/// Strictly ascending runs: strict inequality w[i] < w[i+1]
fn strict_ascending_runs(w: &[u8]) -> Vec<&[u8]> {
    if w.is_empty() {
        return vec![];
    }
    let mut runs = Vec::new();
    let mut start = 0;
    for i in 1..w.len() {
        if w[i] < w[i - 1] {
            runs.push(&w[start..i]);
            start = i;
        }
    }
    runs.push(&w[start..]);
    runs
}

/// Check: min of each run strictly increasing
fn is_run_sorted_strict(w: &[u8]) -> bool {
    let runs = ascending_runs(w);
    if runs.len() <= 1 {
        return true;
    }
    for i in 1..runs.len() {
        let min_prev = runs[i - 1].iter().min().unwrap();
        let min_cur = runs[i].iter().min().unwrap();
        if min_cur <= min_prev {
            return false;
        }
    }
    true
}

/// Check: min of each run weakly increasing
fn is_run_sorted_weak(w: &[u8]) -> bool {
    let runs = ascending_runs(w);
    if runs.len() <= 1 {
        return true;
    }
    for i in 1..runs.len() {
        let min_prev = runs[i - 1].iter().min().unwrap();
        let min_cur = runs[i].iter().min().unwrap();
        if min_cur < min_prev {
            return false;
        }
    }
    true
}

/// Check: runs in lexicographic order
fn is_run_sorted_lex(w: &[u8]) -> bool {
    let runs = ascending_runs(w);
    if runs.len() <= 1 {
        return true;
    }
    for i in 1..runs.len() {
        if runs[i] <= runs[i - 1] {
            return false;
        }
    }
    true
}

fn accumulate_poly(coeffs: &mut Vec<i64>, val: usize) {
    if val >= coeffs.len() {
        coeffs.resize(val + 1, 0);
    }
    coeffs[val] += 1;
}

fn poly_info(objects: &[&Vec<u8>], stat: Stat) -> (Vec<i64>, bool) {
    let mut coeffs: Vec<i64> = Vec::new();
    for w in objects {
        accumulate_poly(&mut coeffs, statistics::compute(w, stat));
    }
    if coeffs.is_empty() {
        coeffs.push(0);
    }
    let rr = if coeffs.iter().all(|&c| c == 0) {
        true
    } else {
        is_real_rooted(&coeffs)
    };
    (coeffs, rr)
}

fn main() {
    let max_n: u8 = 9;

    // Also try with non-strict ascending runs (allowing equal consecutive elements)
    // For words (not permutations), a "run" could mean w[i] <= w[i+1] (weakly ascending)
    // or w[i] < w[i+1] (strictly ascending). For words with repeated letters,
    // the standard is weakly ascending (non-decreasing) runs.
    // Let's do both interpretations.

    println!("# Run-sorted parking functions\n");
    println!("Ascending runs use WEAK inequality (non-decreasing): w[i] <= w[i+1]\n");

    // Collect counts
    println!("## Count table\n");
    println!("| n | PF(n) | strict | weak | lex |");
    println!("|---|-------|--------|------|-----|");

    let mut all_strict: Vec<Vec<Vec<u8>>> = Vec::new();
    let mut all_weak: Vec<Vec<Vec<u8>>> = Vec::new();
    let mut all_lex: Vec<Vec<Vec<u8>>> = Vec::new();

    for n in 1..=max_n {
        let pfs = all_parking_functions(n);
        let total = pfs.len();

        let strict: Vec<Vec<u8>> = pfs
            .iter()
            .filter(|w| is_run_sorted_strict(w))
            .cloned()
            .collect();
        let weak: Vec<Vec<u8>> = pfs
            .iter()
            .filter(|w| is_run_sorted_weak(w))
            .cloned()
            .collect();
        let lex: Vec<Vec<u8>> = pfs
            .iter()
            .filter(|w| is_run_sorted_lex(w))
            .cloned()
            .collect();

        println!(
            "| {} | {} | {} | {} | {} |",
            n,
            total,
            strict.len(),
            weak.len(),
            lex.len()
        );

        all_strict.push(strict);
        all_weak.push(weak);
        all_lex.push(lex);
    }

    // Now try with strictly ascending runs
    println!("\n---\n");
    println!("Ascending runs use STRICT inequality: w[i] < w[i+1]\n");
    println!("| n | PF(n) | strict | weak | lex |");
    println!("|---|-------|--------|------|-----|");

    let mut all_strict2: Vec<Vec<Vec<u8>>> = Vec::new();
    let mut all_weak2: Vec<Vec<Vec<u8>>> = Vec::new();
    let mut all_lex2: Vec<Vec<Vec<u8>>> = Vec::new();

    for n in 1..=max_n {
        let pfs = all_parking_functions(n);
        let total = pfs.len();

        let strict: Vec<Vec<u8>> = pfs
            .iter()
            .filter(|w| {
                let runs = strict_ascending_runs(w);
                runs.len() <= 1
                    || (1..runs.len())
                        .all(|i| runs[i].iter().min().unwrap() > runs[i - 1].iter().min().unwrap())
            })
            .cloned()
            .collect();

        let weak: Vec<Vec<u8>> = pfs
            .iter()
            .filter(|w| {
                let runs = strict_ascending_runs(w);
                runs.len() <= 1
                    || (1..runs.len())
                        .all(|i| runs[i].iter().min().unwrap() >= runs[i - 1].iter().min().unwrap())
            })
            .cloned()
            .collect();

        let lex: Vec<Vec<u8>> = pfs
            .iter()
            .filter(|w| {
                let runs = strict_ascending_runs(w);
                runs.len() <= 1 || (1..runs.len()).all(|i| runs[i] > runs[i - 1])
            })
            .cloned()
            .collect();

        println!(
            "| {} | {} | {} | {} | {} |",
            n,
            total,
            strict.len(),
            weak.len(),
            lex.len()
        );

        all_strict2.push(strict);
        all_weak2.push(weak);
        all_lex2.push(lex);
    }

    // Stat polynomials for the non-decreasing run version
    let stats = [Stat::Des, Stat::Exc, Stat::Inv, Stat::Maj, Stat::Peak];

    for (label, collection) in [
        ("strict (non-decr runs)", &all_strict),
        ("weak (non-decr runs)", &all_weak),
        ("lex (non-decr runs)", &all_lex),
        ("strict (strict-asc runs)", &all_strict2),
        ("weak (strict-asc runs)", &all_weak2),
        ("lex (strict-asc runs)", &all_lex2),
    ] {
        println!("\n## Polynomials: {}\n", label);
        for stat in &stats {
            println!("### stat = {:?}\n", stat);
            for (idx, objects) in collection.iter().enumerate() {
                let n = idx + 1;
                let refs: Vec<&Vec<u8>> = objects.iter().collect();
                if refs.is_empty() {
                    println!("n={}: 0 objects", n);
                    continue;
                }
                let (coeffs, rr) = poly_info(&refs, *stat);
                println!(
                    "n={}: {:>8} objects, rr={:5}, poly = {}",
                    n,
                    refs.len(),
                    rr,
                    format_poly(&coeffs)
                );
            }
            println!();
        }
    }
}
