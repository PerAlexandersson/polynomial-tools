/// Compute example polynomial sequences for the polynomial-tools web interface.
use combpoly::permutation::all_permutations;
use combpoly::polynomial_builder::build_generating_polynomial;
use combpoly::statistics::Stat;
use polynomial_tools::real_rootedness::format_poly;
use polynomial_tools::sequences::*;

fn main() {
    println!("=== Eulerian polynomials A_n(t) ===");
    let ep = eulerian_polynomials(10);
    for (i, p) in ep.iter().enumerate() {
        println!("# A_{}(t)", i + 1);
        let trimmed: Vec<i64> = p
            .iter()
            .copied()
            .rev()
            .skip_while(|&c| c == 0)
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect();
        println!("{}", format_poly(&trimmed));
    }

    println!("\n=== Narayana polynomials N_n(t) ===");
    let np = narayana_polynomials(10);
    for (i, p) in np.iter().enumerate() {
        println!("# N_{}(t)", i + 1);
        let trimmed: Vec<i64> = p
            .iter()
            .copied()
            .rev()
            .skip_while(|&c| c == 0)
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect();
        println!("{}", format_poly(&trimmed));
    }

    println!("\n=== Descent on derangements d_n(t) ===");
    for n in 2..=10u8 {
        let perms = all_permutations(n);
        let derangements: Vec<Vec<u8>> = perms
            .into_iter()
            .filter(|w| w.iter().enumerate().all(|(i, &v)| v != (i as u8 + 1)))
            .collect();
        let poly = build_generating_polynomial(&derangements, Stat::Des);
        println!("# d_{}(t)  ({} derangements)", n, derangements.len());
        println!("{}", format_poly(&poly));
    }

    // Fibonacci matching polynomials on 2×n grid
    // f_1 = 1 + t, f_2 = 1 + 3t + t^2, f_n = (1+t)f_{n-1} + t*f_{n-2}
    // (matchings of ladder graph L_n by number of edges)
    println!("\n=== Fibonacci matching polynomials (2×n grid) ===");
    let mut fib: Vec<Vec<i64>> = Vec::new();
    // f_0 = 1 (empty matching)
    fib.push(vec![1]);
    // f_1 = 1 + t (empty or the one edge)
    // Actually for 2x1 grid (single edge): matchings are {} and {e}, so f_1 = 1 + t
    fib.push(vec![1, 1]);
    // For 2×n grid: f_n = (1+2t)*f_{n-1} - t^2*f_{n-2} might not be right.
    // Let me use the transfer matrix method for matchings on P_2 × P_n.
    // Actually the standard recurrence for matchings on 2×n grid:
    // Let M(n,k) = matchings of 2×n grid with k edges.
    // If the rightmost column (col n) has no edges involving it: M(n-1,k)
    // One horizontal edge in row 1 or row 2: 2*M(n-1,k-1)... this gets complex.
    //
    // Let me just enumerate matchings directly for small n using a DP approach.
    // State = which vertices in the current column are already matched to the left.
    // Column has 2 vertices: states are 00, 01, 10, 11 (bitmask).
    // Transition: from state s of col i to state s' of col i+1, with some edges used.

    // Actually, let me just compute these via a simpler recurrence.
    // m(n) = total matchings of 2×n grid.
    // m(0) = 1, m(1) = 2, m(n) = 2*m(n-1) + m(n-2) - wait that's not right either.
    // Let me just enumerate.

    // For the web examples, let me use a known result:
    // The matching polynomial of the ladder graph L_n (= P_2 × P_n) tracks matchings by size.
    // These are related to Fibonacci but I'll compute them via transfer matrix.

    // States for column boundary: subset of {top, bottom} already matched.
    // 4 states: 0b00=neither, 0b01=bottom, 0b10=top, 0b11=both
    // For each column, we can:
    // - place vertical edge (matches both vertices): need both free, next state = 00, +1 edge
    // - place horizontal top edge: need top free, next state has top matched, +1 edge
    // - place horizontal bottom edge: need bottom free, next state has bottom matched, +1 edge
    // - place both horizontal edges: need both free, next state = 11, +2 edges
    // - place nothing: next state inherits

    // This is getting complex. Let me use a simple DP.
    // dp[col][state] = polynomial (Vec<i64>) counting matchings

    let max_grid = 10;
    let mut results: Vec<Vec<i64>> = Vec::new();

    for n in 1..=max_grid {
        // Transfer matrix DP for matchings on 2×n grid
        // dp[mask] = polynomial where mask indicates which vertices in current column
        // are matched by horizontal edges FROM the previous column.
        // mask: bit 0 = top vertex pre-matched, bit 1 = bottom vertex pre-matched

        let mut dp: Vec<Vec<i64>> = vec![vec![]; 4];
        dp[0] = vec![1]; // start: no vertices pre-matched

        for _col in 0..n {
            let mut ndp: Vec<Vec<i64>> = vec![vec![]; 4];

            for inmask in 0..4u8 {
                if dp[inmask as usize].is_empty() {
                    continue;
                }
                let p = &dp[inmask as usize];
                let top_free = inmask & 1 == 0;
                let bot_free = inmask & 2 == 0;

                // Enumerate what we do in this column:
                // 1. Vertical edge (if both free): uses both, output mask depends on horizontal to right
                // 2. Horizontal edges to the right
                // 3. Nothing

                // Actually, let me think of it differently.
                // For each column, the free vertices (not pre-matched from left) can be:
                // - left unmatched
                // - matched by a vertical edge within the column (need both free)
                // - matched by a horizontal edge to the right (sets the corresponding bit in output mask)

                // Let's enumerate all valid local configurations:
                // Case by case based on which vertices are free

                if top_free && bot_free {
                    // Both free. Options:
                    // a) do nothing: both pass through unmatched, out=00
                    poly_add_to(&mut ndp[0], p, 0);
                    // b) vertical edge: both matched, out=00, +1 edge
                    poly_add_to(&mut ndp[0], p, 1);
                    // c) horizontal top to right: top matched, out=01 (top of next col pre-matched)
                    poly_add_to(&mut ndp[1], p, 1);
                    // d) horizontal bot to right: bot matched, out=10
                    poly_add_to(&mut ndp[2], p, 1);
                    // e) both horizontal to right: out=11, +2 edges
                    poly_add_to(&mut ndp[3], p, 2);
                    // f) vertical + (nothing, already handled)
                    // g) horizontal top + horizontal bot to right: same as (e)? No, that's 2 horizontal.
                    // h) vertical edge + horizontal edges: can't, vertical uses both
                }
                if top_free && !bot_free {
                    // Top free, bottom pre-matched.
                    // a) do nothing with top
                    poly_add_to(&mut ndp[0], p, 0);
                    // b) horizontal top to right
                    poly_add_to(&mut ndp[1], p, 1);
                }
                if !top_free && bot_free {
                    // Top pre-matched, bottom free.
                    // a) do nothing with bottom
                    poly_add_to(&mut ndp[0], p, 0);
                    // b) horizontal bot to right
                    poly_add_to(&mut ndp[2], p, 1);
                }
                if !top_free && !bot_free {
                    // Both pre-matched, nothing to do
                    poly_add_to(&mut ndp[0], p, 0);
                }
            }

            dp = ndp;
        }

        // Final: only state 0 (no pending horizontal edges) is valid
        let poly = if dp[0].is_empty() {
            vec![0]
        } else {
            dp[0].clone()
        };
        println!("# F_{}(t)  (matchings on 2×{} grid)", n, n);
        println!("{}", format_poly(&poly));
        results.push(poly);
    }
}

fn poly_add_to(target: &mut Vec<i64>, source: &[i64], shift: usize) {
    let needed = source.len() + shift;
    if target.len() < needed {
        target.resize(needed, 0);
    }
    for (i, &c) in source.iter().enumerate() {
        target[i + shift] += c;
    }
}
