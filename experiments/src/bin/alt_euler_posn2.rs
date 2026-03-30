/// Explore interlacing of H_{n,j}(t) in various orderings.
///
/// Key observation: H_{n,n} = H_{n-1} (total), and there's a
/// palindromic symmetry H_{n,j} <-> t*rev(H_{n,n+2-j}).
/// Try interlacing from outside in, and other orderings.
use combpoly::permutation::alternating_permutations;
use combpoly::statistics::{compute, Stat};
use polynomial_tools::real_rootedness::{
    check_interlacing, check_interlacing_sturm, check_weak_interlacing, format_poly,
};

fn build_poly(perms: &[&Vec<u8>]) -> Vec<i64> {
    if perms.is_empty() { return vec![0]; }
    let max_s = perms.iter().map(|s| compute(s, Stat::Swaps)).max().unwrap_or(0);
    let mut coeffs = vec![0i64; max_s + 1];
    for s in perms { coeffs[compute(s, Stat::Swaps)] += 1; }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 { coeffs.pop(); }
    coeffs
}

fn check_il(f: &[i64], g: &[i64]) -> &'static str {
    if f.len() <= 1 || g.len() <= 1 { return "trivial"; }
    // Try f << g (f smaller degree)
    let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
    match check_interlacing_sturm(small, large) {
        Some(true) => "✓",
        Some(false) => {
            // Try weak interlacing
            match check_weak_interlacing(small, large) {
                Some(true) => "✓ (weak)",
                _ => "✗",
            }
        }
        None => "?",
    }
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(13);

    for n in 3..=max_n {
        let alt = alternating_permutations(n);

        // Group by position of n (0-indexed internally, 1-indexed for display)
        let mut by_pos: Vec<Vec<&Vec<u8>>> = vec![vec![]; n as usize];
        for s in &alt {
            let pos = s.iter().position(|&v| v == n).unwrap();
            by_pos[pos].push(s);
        }

        // Collect active positions (even positions, 1-indexed = peaks)
        let mut active: Vec<(usize, Vec<i64>)> = Vec::new();
        for (pos, group) in by_pos.iter().enumerate() {
            if group.is_empty() { continue; }
            active.push((pos + 1, build_poly(group)));
        }

        println!("=== n={} ({} perms, {} positions) ===", n, alt.len(), active.len());
        for (pos, poly) in &active {
            let count: i64 = poly.iter().sum();
            println!("  pos={:>2}: count={:<6} {}", pos, count, format_poly(poly));
        }

        if active.len() < 2 {
            println!();
            continue;
        }

        // Ordering 1: outside-in (last, first, second-to-last, second, ...)
        let mut outside_in: Vec<usize> = Vec::new();
        let (mut lo, mut hi) = (0usize, active.len() - 1);
        let mut from_hi = true;
        while lo <= hi {
            if from_hi {
                outside_in.push(hi);
                if hi == 0 { break; }
                hi -= 1;
            } else {
                outside_in.push(lo);
                lo += 1;
            }
            from_hi = !from_hi;
        }

        println!("  Outside-in order: {}", outside_in.iter()
            .map(|&i| format!("pos {}", active[i].0)).collect::<Vec<_>>().join(" → "));
        for w in outside_in.windows(2) {
            let (p1, ref f) = active[w[0]];
            let (p2, ref g) = active[w[1]];
            println!("    pos {} ≪ pos {}: {}", p1, p2, check_il(f, g));
        }

        // Ordering 2: last to first (decreasing position)
        println!("  Decreasing order:");
        for i in (0..active.len() - 1).rev() {
            let (p1, ref f) = active[i + 1];
            let (p2, ref g) = active[i];
            println!("    pos {} ≪ pos {}: {}", p1, p2, check_il(f, g));
        }

        // Ordering 3: first to last (increasing position)
        println!("  Increasing order:");
        for i in 0..active.len() - 1 {
            let (p1, ref f) = active[i];
            let (p2, ref g) = active[i + 1];
            println!("    pos {} ≪ pos {}: {}", p1, p2, check_il(f, g));
        }

        // Ordering 4: cumulative sums from the right
        // S_k = H_{n,n} + H_{n,n-2} + ... + H_{n,k}
        println!("  Cumulative sums (from right):");
        let mut cum_polys: Vec<(String, Vec<i64>)> = Vec::new();
        let mut running = vec![0i64; 1];
        for i in (0..active.len()).rev() {
            let (pos, ref p) = active[i];
            // Add p to running
            let new_len = running.len().max(p.len());
            running.resize(new_len, 0);
            for (j, &c) in p.iter().enumerate() {
                running[j] += c;
            }
            while running.len() > 1 && *running.last().unwrap() == 0 {
                running.pop();
            }
            cum_polys.push((format!("Σ pos≥{}", pos), running.clone()));
        }
        cum_polys.reverse();
        for (label, poly) in &cum_polys {
            let count: i64 = poly.iter().sum();
            println!("    {}: count={:<6} {}", label, count, format_poly(poly));
        }
        println!("  Cumulative interlacing:");
        for w in cum_polys.windows(2) {
            let (ref l1, ref f) = w[0];
            let (ref l2, ref g) = w[1];
            let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
            println!("    {} ≪ {}: {}", l1, l2, check_il(small, large));
        }

        println!();
    }
}
