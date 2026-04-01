/// Refine H_n(t) by position of the maximum value n:
/// H_{n,j}(t) = sum_{σ ∈ A_n, σ(j)=n} t^{swaps(σ)}.
///
/// Since alternating perms are up-down, n must appear at an even
/// position (1-indexed), i.e., at a peak.
use combpoly::permutation::alternating_permutations;
use combpoly::statistics::{compute, Stat};
use polynomial_tools::real_rootedness::{check_interlacing_sturm, format_poly, is_real_rooted};
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};

fn build_poly(perms: &[&Vec<u8>]) -> Vec<i64> {
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

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(13);

    // Collect all H_{n,j} polys indexed by position j
    // For the recurrence search, collect sequences for each fixed position j
    let mut polys_by_pos: Vec<Vec<Vec<i64>>> = Vec::new(); // polys_by_pos[j][n_idx]

    for n in 1..=max_n {
        let alt = alternating_permutations(n);

        println!("=== n={} ({} perms) ===", n, alt.len());

        // Group by position of n (0-indexed)
        let mut by_pos: Vec<Vec<&Vec<u8>>> = vec![vec![]; n as usize];
        for s in &alt {
            let pos = s.iter().position(|&v| v == n).unwrap();
            by_pos[pos].push(s);
        }

        let mut active_polys: Vec<(usize, Vec<i64>)> = Vec::new();
        for (pos, group) in by_pos.iter().enumerate() {
            if group.is_empty() {
                continue;
            }
            let poly = build_poly(group);
            let count: i64 = poly.iter().sum();
            let rr = if poly.len() <= 2 {
                true
            } else {
                is_real_rooted(&poly)
            };
            let pos1 = pos + 1; // 1-indexed
            println!(
                "  pos={:>2} ({}): count={:<6} {:<50} {}",
                pos1,
                if pos1 % 2 == 0 { "peak" } else { "    " },
                count,
                format_poly(&poly),
                if rr { "✓ rr" } else { "✗ NOT rr" }
            );
            active_polys.push((pos1, poly));
        }

        // Check interlacing among the position-refined polys
        if active_polys.len() >= 2 {
            println!("  Interlacing (consecutive positions):");
            for i in 0..active_polys.len() - 1 {
                let (p1, ref f) = active_polys[i];
                let (p2, ref g) = active_polys[i + 1];
                if f.len() <= 1 || g.len() <= 1 {
                    continue;
                }
                let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
                match check_interlacing_sturm(small, large) {
                    Some(true) => println!("    pos {} ≪ pos {}: ✓", p1, p2),
                    Some(false) => println!("    pos {} ≪ pos {}: ✗", p1, p2),
                    None => println!("    pos {} ≪ pos {}: ?", p1, p2),
                }
            }
        }

        // Store polys for recurrence search (only even positions = peaks)
        while polys_by_pos.len() < n as usize {
            polys_by_pos.push(Vec::new());
        }
        for (pos1, poly) in &active_polys {
            if *pos1 <= polys_by_pos.len() {
                polys_by_pos[pos1 - 1].push(poly.clone());
            }
        }

        println!();
    }

    // Search for recurrences in individual position sequences
    println!("=== Recurrence search for fixed-position sequences ===\n");
    let opts = AdaptiveSearchOptions {
        max_rec_len: 3,
        max_var_deg: 2,
        max_idx_deg: 2,
        max_diff_deg: 1,
        try_inhomogeneous: false,
        try_denominator: true,
        max_denom_var_deg: 1,
        max_denom_idx_deg: 1,
        min_margin: 2,
        verbose: false,
    };

    for (j, seq) in polys_by_pos.iter().enumerate() {
        if seq.len() < 4 {
            continue;
        }
        let pos1 = j + 1;
        print!("  pos={}: {} terms, ", pos1, seq.len());
        match find_recurrence_adaptive(seq, &opts) {
            Some(res) => {
                let s = format!("{}", res.recurrence);
                if s.len() > 200 {
                    println!("overfitted (coeff too large)");
                } else {
                    println!("FOUND: {}", s);
                }
            }
            None => println!("none"),
        }
    }
}
