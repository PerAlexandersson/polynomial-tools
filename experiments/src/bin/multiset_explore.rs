//! Explore the multiset rook-Eulerian conjecture:
//!   {R_{λ_1}, R_{λ_1-1}, ..., R_1} is an interlacing sequence.
//!
//! R_j(λ,α;t) = ascent-generating poly for words on board λ with content α
//!              and first entry j.
//!
//! Recursion: R_i(λ+,α+e_i;t) = Σ_{j≤i} R_j(λ,α;t) + t·Σ_{j>i} R_j(λ,α;t)
//!
//! The matrix for this recursion has entries 1 (j≤i) or t (j>i), which is
//! a staircase matrix similar to the rook-Eulerian one.
//! Check if Brändén's Corollary 8.7 applies.

use std::collections::HashMap;

fn exact_interlaces(fc: &[i64], gc: &[i64]) -> bool {
    polynomial_tools::check_weak_interlacing(fc, gc).unwrap_or(false)
}

fn pt(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
    v
}
fn pz(p: &[i64]) -> bool {
    p.iter().all(|&c| c == 0)
}
fn pa(a: &[i64], b: &[i64]) -> Vec<i64> {
    let l = a.len().max(b.len());
    let mut r = vec![0i64; l];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] += v;
    }
    pt(&r)
}
fn pf(p: &[i64]) -> String {
    let p = pt(p);
    if pz(&p) {
        "0".into()
    } else {
        let mut t = vec![];
        for (i, &c) in p.iter().enumerate() {
            if c == 0 {
                continue;
            }
            match (c, i) {
                (c, 0) => t.push(format!("{}", c)),
                (1, 1) => t.push("t".into()),
                (c, 1) => t.push(format!("{}t", c)),
                (1, e) => t.push(format!("t^{}", e)),
                (c, e) => t.push(format!("{}t^{}", c, e)),
            }
        }
        t.join(" + ")
    }
}

/// Generate all words with given content on a board.
fn words_on_board(content: &[usize], board: &[u8]) -> Vec<Vec<u8>> {
    let n = board.len();
    let k = content.len();
    let total: usize = content.iter().sum();
    if total != n {
        return vec![];
    }
    let mut result = Vec::new();
    let mut current = Vec::with_capacity(n);
    let mut remaining = content.to_vec();
    gen_words(0, n, k, board, &mut remaining, &mut current, &mut result);
    result
}

fn gen_words(
    pos: usize,
    n: usize,
    k: usize,
    board: &[u8],
    remaining: &mut Vec<usize>,
    current: &mut Vec<u8>,
    result: &mut Vec<Vec<u8>>,
) {
    if pos == n {
        result.push(current.clone());
        return;
    }
    let hi = (board[pos] as usize).min(k);
    for letter in 1..=hi {
        let j = letter - 1;
        if remaining[j] > 0 {
            remaining[j] -= 1;
            current.push(letter as u8);
            gen_words(pos + 1, n, k, board, remaining, current, result);
            current.pop();
            remaining[j] += 1;
        }
    }
}

fn ascents(w: &[u8]) -> usize {
    (1..w.len()).filter(|&i| w[i - 1] < w[i]).count()
}

fn main() {
    println!("=== Multiset rook-Eulerian conjecture ===\n");
    println!("Testing: {{R_m, R_{{m-1}}, ..., R_1}} is interlacing (reversed).\n");

    let mut total_boards = 0u64;
    let mut total_pairs = 0u64;
    let mut interl_fails = 0u64;
    let mut rr_fails = 0u64;

    // Test various (board, content) pairs
    // Boards: Ferrers boards λ = (λ_1,...,λ_n) with λ_i ≥ 1, λ non-decreasing
    // Contents: α = (α_1,...,α_k) with Σα_i = n, k ≤ λ_1

    // Generate test cases systematically
    for n in 2..=7 {
        // Generate some boards of length n
        let boards = gen_ferrers_boards(n);
        for board in &boards {
            let m = board[0] as usize; // λ_1 = max letter
                                       // Generate some contents
            let contents = gen_contents(n, m);
            for content in &contents {
                total_boards += 1;

                let words = words_on_board(content, board);
                if words.is_empty() {
                    continue;
                }

                // Compute R_j for each first entry j
                let mut r_j: Vec<Vec<i64>> = vec![vec![0]; m + 1]; // 1-indexed
                for w in &words {
                    let j = w[0] as usize;
                    let asc = ascents(w);
                    if j <= m {
                        while r_j[j].len() <= asc {
                            r_j[j].push(0);
                        }
                        r_j[j][asc] += 1;
                    }
                }
                for j in 1..=m {
                    r_j[j] = pt(&r_j[j]);
                }

                // Check total is real-rooted
                let mut total_poly = vec![0i64];
                for j in 1..=m {
                    total_poly = pa(&total_poly, &r_j[j]);
                }
                if !pz(&total_poly) {
                    if !polynomial_tools::is_real_rooted(&total_poly) {
                        rr_fails += 1;
                        println!(
                            "FAIL RR: board={:?}, content={:?}, poly={}",
                            board,
                            content,
                            pf(&total_poly)
                        );
                    }
                }

                // Check reversed interlacing: R_m ≼ R_{m-1} ≼ ... ≼ R_1
                for j in (2..=m).rev() {
                    let l = j;
                    let j2 = j - 1;
                    if pz(&r_j[l]) || pz(&r_j[j2]) {
                        continue;
                    }
                    total_pairs += 1;
                    if !exact_interlaces(&r_j[l], &r_j[j2]) {
                        interl_fails += 1;
                        if interl_fails <= 5 {
                            println!(
                                "FAIL: R_{} ≼ R_{} for board={:?}, content={:?}",
                                l, j2, board, content
                            );
                            println!("  R_{} = {}", l, pf(&r_j[l]));
                            println!("  R_{} = {}", j2, pf(&r_j[j2]));
                        }
                    }
                }
            }
        }
        println!(
            "n={}: boards={}, pairs={}, RR_fails={}, interl_fails={}",
            n, total_boards, total_pairs, rr_fails, interl_fails
        );
    }

    println!("\n=== SUMMARY ===");
    println!("Total (board,content) pairs: {}", total_boards);
    println!("Interlacing pairs tested: {}", total_pairs);
    println!("Real-rootedness failures: {}", rr_fails);
    println!("Interlacing failures: {}", interl_fails);
    if rr_fails == 0 && interl_fails == 0 {
        println!("\nALL CONDITIONS HOLD.");
    }

    // Print the paper's example
    println!("\n=== Paper example: board=22233, content=(2,2,1) ===");
    let board = vec![2u8, 2, 2, 3, 3];
    let content = vec![2usize, 2, 1];
    let words = words_on_board(&content, &board);
    let m = 3;
    let mut r_j: Vec<Vec<i64>> = vec![vec![0]; m + 1];
    for w in &words {
        let j = w[0] as usize;
        let asc = ascents(w);
        while r_j[j].len() <= asc {
            r_j[j].push(0);
        }
        r_j[j][asc] += 1;
    }
    for j in 1..=m {
        r_j[j] = pt(&r_j[j]);
        println!("R_{} = {}", j, pf(&r_j[j]));
    }
    let mut total_poly = vec![0i64];
    for j in 1..=m {
        total_poly = pa(&total_poly, &r_j[j]);
    }
    println!("R = {}", pf(&total_poly));
}

/// Generate Ferrers boards of length n: non-decreasing, λ_i ≥ 1.
fn gen_ferrers_boards(n: usize) -> Vec<Vec<u8>> {
    let max_val = n + 2; // allow some room
    let mut results = vec![];
    let mut current = vec![];
    gfb(n, max_val, 0, &mut current, &mut results);
    results
}

fn gfb(n: usize, max_val: usize, depth: usize, current: &mut Vec<u8>, results: &mut Vec<Vec<u8>>) {
    if depth == n {
        results.push(current.clone());
        return;
    }
    let min_val = 1.max(if depth > 0 {
        current[depth - 1] as usize
    } else {
        1
    });
    for v in min_val..=max_val.min(n + 2) {
        current.push(v as u8);
        gfb(n, max_val, depth + 1, current, results);
        current.pop();
    }
}

/// Generate some contents with sum = n and max letter ≤ max_k.
fn gen_contents(n: usize, max_k: usize) -> Vec<Vec<usize>> {
    let mut results = vec![];
    let k = max_k.min(n); // at most n distinct letters
    gen_comp(n, k, &mut vec![], &mut results);
    results
}

fn gen_comp(
    remaining: usize,
    parts_left: usize,
    current: &mut Vec<usize>,
    results: &mut Vec<Vec<usize>>,
) {
    if parts_left == 0 {
        if remaining == 0 {
            results.push(current.clone());
        }
        return;
    }
    if parts_left == 1 {
        current.push(remaining);
        results.push(current.clone());
        current.pop();
        return;
    }
    // Each part ≥ 0, but at least 1 for the last part? No, content can have 0's.
    // Actually content entries should be ≥ 0 with sum = n.
    // But having α_i = 0 means letter i doesn't appear. Let's allow it for generality.
    // To keep the search space manageable, require each α_i ≥ 1.
    for v in 1..=remaining.saturating_sub(parts_left - 1) {
        current.push(v);
        gen_comp(remaining - v, parts_left - 1, current, results);
        current.pop();
    }
}
