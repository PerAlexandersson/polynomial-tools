/// Explore routes to proving real-rootedness and interlacing of H_n(t).
///
/// 1. γ-positivity (palindromic case k=2): implies real-rootedness
/// 2. Position-refined interlacing chain: H_{n,2} ≪ H_{n,4} ≪ ... ≪ H_{n,n}
/// 3. Cumulative sums: t·H_{n-1} ≪ S_{n,n-2} ≪ ... ≪ H_n
/// 4. Hermite-Biehler for non-palindromic k≥3: even/odd part interlacing
/// 5. Common interlacers: do all H_{n,j} interlace a single polynomial?
/// 6. γ-coefficients: look for combinatorial patterns
use combpoly::permutation::{all_permutations, alternating_permutations};
use combpoly::statistics::{compute, Stat};
use polynomial_tools::real_rootedness::{
    check_interlacing, check_weak_interlacing, format_poly, gamma_coefficients, is_gamma_positive,
    is_real_rooted,
};

/// Build swaps polynomial on a set of permutations.
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

/// k-alternating: descent at position j iff k | j.
fn is_k_alternating(sigma: &[u8], k: usize) -> bool {
    for i in 0..sigma.len() - 1 {
        let pos = i + 1; // 1-indexed
        let is_descent = sigma[i] > sigma[i + 1];
        let should_descend = pos % k == 0;
        if is_descent != should_descend {
            return false;
        }
    }
    true
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(12);

    // ========================================================
    // Part 1: γ-positivity of H_n(t) (k=2, palindromic case)
    // ========================================================
    println!("═══ Part 1: γ-positivity of H_n(t) ═══\n");

    let mut h_polys: Vec<Vec<i64>> = Vec::new();

    for n in 1..=max_n {
        let alt = alternating_permutations(n);
        let refs: Vec<&Vec<u8>> = alt.iter().collect();
        let poly = build_poly(&refs);
        let gp = is_gamma_positive(&poly);
        let gamma = gamma_coefficients(&poly).unwrap_or_default();

        println!(
            "H_{:<2}(t) = {:<55} γ-pos: {:<5} γ = {:?}",
            n,
            format_poly(&poly),
            gp,
            gamma,
        );
        h_polys.push(poly);
    }

    // ========================================================
    // Part 2: Position-refined interlacing chain
    // ========================================================
    println!("\n═══ Part 2: Position-refined H_{{n,j}}(t) interlacing ═══\n");

    for n in (4..=max_n.min(14)).step_by(2) {
        let alt = alternating_permutations(n);

        // Group by position of max value n
        let mut by_pos: std::collections::BTreeMap<usize, Vec<&Vec<u8>>> =
            std::collections::BTreeMap::new();
        for s in &alt {
            let pos = s.iter().position(|&v| v == n).unwrap() + 1; // 1-indexed
            by_pos.entry(pos).or_default().push(s);
        }

        let positions: Vec<usize> = by_pos.keys().copied().collect();
        let pos_polys: Vec<(usize, Vec<i64>)> = positions
            .iter()
            .map(|&p| (p, build_poly(by_pos.get(&p).unwrap())))
            .collect();

        println!("n={}:", n);
        for (p, poly) in &pos_polys {
            let gp = is_gamma_positive(poly);
            println!(
                "  H_{{{},{}}}(t) = {:<50} γ-pos: {}",
                n,
                p,
                format_poly(poly),
                gp,
            );
        }

        // Check pairwise interlacing in position order
        let mut chain_ok = true;
        for i in 0..pos_polys.len() - 1 {
            let (p1, f) = &pos_polys[i];
            let (p2, g) = &pos_polys[i + 1];
            let interlaces = if f.len() <= g.len() {
                check_interlacing(f, g)
            } else {
                check_interlacing(g, f)
            };
            let ok = interlaces == Some(true);
            if !ok {
                chain_ok = false;
                println!("  H_{{{},{}}} ≪ H_{{{},{}}}: ✗", n, p1, n, p2);
            }
        }
        if chain_ok {
            let chain_str: String = positions
                .iter()
                .map(|p| format!("H_{{{},{}}}", n, p))
                .collect::<Vec<_>>()
                .join(" ≪ ");
            println!("  ✓ {}", chain_str);
        }

        // Check: do ALL H_{n,j} interlace H_{n,n} = t·H_{n-1}?
        let last_poly = &pos_polys.last().unwrap().1;
        let mut all_interlace_last = true;
        for (p, poly) in &pos_polys[..pos_polys.len() - 1] {
            let int = if poly.len() <= last_poly.len() {
                check_interlacing(poly, last_poly)
            } else {
                check_interlacing(last_poly, poly)
            };
            if int != Some(true) {
                all_interlace_last = false;
                println!(
                    "  H_{{{},{}}} does NOT interlace H_{{{},{}}}",
                    n,
                    p,
                    n,
                    positions.last().unwrap()
                );
            }
        }
        if all_interlace_last {
            println!(
                "  ✓ All H_{{{},j}} interlace H_{{{},{}}} = t·H_{}(t)",
                n,
                n,
                positions.last().unwrap(),
                n - 1
            );
        }

        // Cumulative sums from the right
        let mut cum_polys: Vec<Vec<i64>> = Vec::new();
        let mut running = vec![0i64; pos_polys.last().unwrap().1.len()];
        for (_, poly) in pos_polys.iter().rev() {
            // Add poly to running sum
            if running.len() < poly.len() {
                running.resize(poly.len(), 0);
            }
            for (i, &c) in poly.iter().enumerate() {
                running[i] += c;
            }
            cum_polys.push(running.clone());
        }
        cum_polys.reverse();
        // cum_polys[0] = H_n, cum_polys[last] = H_{n,n}

        let mut cum_chain_ok = true;
        for i in 0..cum_polys.len() - 1 {
            let (f, g) = (&cum_polys[i + 1], &cum_polys[i]);
            let int = if f.len() <= g.len() {
                check_weak_interlacing(f, g)
            } else {
                check_weak_interlacing(g, f)
            };
            if int != Some(true) {
                cum_chain_ok = false;
            }
        }
        if cum_chain_ok {
            println!("  ✓ Cumulative sum chain S_{{n,n}} ≪ ... ≪ S_{{n,2}} = H_n");
        }

        println!();
    }

    // ========================================================
    // Part 3: k=3 (non-palindromic) — Hermite-Biehler
    // ========================================================
    println!("═══ Part 3: k=3 Hermite-Biehler (even/odd part interlacing) ═══\n");

    for n in 4..=max_n.min(11) {
        let all = all_permutations(n);
        let k3: Vec<&Vec<u8>> = all.iter().filter(|s| is_k_alternating(s, 3)).collect();
        if k3.is_empty() {
            continue;
        }
        let poly = build_poly(&k3);
        let rr = is_real_rooted(&poly);

        // Hermite-Biehler: split into even and odd parts
        let mut even_part: Vec<i64> = Vec::new();
        let mut odd_part: Vec<i64> = Vec::new();
        for (i, &c) in poly.iter().enumerate() {
            if i % 2 == 0 {
                even_part.push(c);
            } else {
                odd_part.push(c);
            }
        }
        // Trim
        while even_part.len() > 1 && *even_part.last().unwrap() == 0 {
            even_part.pop();
        }
        while odd_part.len() > 1 && *odd_part.last().unwrap() == 0 {
            odd_part.pop();
        }

        let even_rr = even_part.len() <= 2 || is_real_rooted(&even_part);
        let odd_rr = odd_part.len() <= 2 || is_real_rooted(&odd_part);

        let (smaller, larger) = if even_part.len() <= odd_part.len() {
            (&even_part, &odd_part)
        } else {
            (&odd_part, &even_part)
        };
        let parts_interlace = if smaller.len() <= 1 {
            true
        } else {
            check_interlacing(smaller, larger) == Some(true)
        };

        println!(
            "H_{}^(3)(t) = {:<50} rr:{} even_rr:{} odd_rr:{} parts_interlace:{}",
            n,
            format_poly(&poly),
            rr,
            even_rr,
            odd_rr,
            parts_interlace,
        );
    }

    // ========================================================
    // Part 4: γ-coefficient patterns — look for combinatorics
    // ========================================================
    println!("\n═══ Part 4: γ-coefficient sequences ═══\n");
    println!("Looking for patterns in γ-coefficients of H_n(t):\n");

    for n in 1..=max_n {
        let gamma = gamma_coefficients(&h_polys[n as usize - 1]).unwrap_or_default();
        // γ_0 should always be 1
        // Print the sequence for OEIS lookup
        let gamma_str: String = gamma
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join(", ");
        println!("n={:>2}: γ = [{}]", n, gamma_str);
    }

    println!("\nγ_1 values (for OEIS):");
    let gamma1: Vec<String> = (1..=max_n as usize)
        .map(|n| {
            let g = gamma_coefficients(&h_polys[n - 1]).unwrap_or_default();
            if g.len() > 1 {
                g[1].to_string()
            } else {
                "-".to_string()
            }
        })
        .collect();
    println!("  {}", gamma1.join(", "));

    println!("\nγ_2 values (for OEIS):");
    let gamma2: Vec<String> = (1..=max_n as usize)
        .map(|n| {
            let g = gamma_coefficients(&h_polys[n - 1]).unwrap_or_default();
            if g.len() > 2 {
                g[2].to_string()
            } else {
                "-".to_string()
            }
        })
        .collect();
    println!("  {}", gamma2.join(", "));

    // ========================================================
    // Part 5: Individual H_{n,j} γ-positivity
    // ========================================================
    println!("\n═══ Part 5: Are individual H_{{n,j}}(t) γ-positive? ═══\n");

    for n in (6..=max_n.min(14)).step_by(2) {
        let alt = alternating_permutations(n);
        let mut by_pos: std::collections::BTreeMap<usize, Vec<&Vec<u8>>> =
            std::collections::BTreeMap::new();
        for s in &alt {
            let pos = s.iter().position(|&v| v == n).unwrap() + 1;
            by_pos.entry(pos).or_default().push(s);
        }

        let mut all_gp = true;
        for (&p, perms) in &by_pos {
            let poly = build_poly(perms);
            let gp = is_gamma_positive(&poly);
            if !gp {
                all_gp = false;
                println!(
                    "  n={} j={}: NOT γ-positive: {}  γ={:?}",
                    n,
                    p,
                    format_poly(&poly),
                    gamma_coefficients(&poly).unwrap_or_default(),
                );
            }
        }
        if all_gp {
            println!("n={}: all H_{{n,j}} are γ-positive", n);
        }
    }

    // ========================================================
    // Part 6: Log-concavity of γ-coefficients
    // ========================================================
    println!("\n═══ Part 6: Are γ-coefficients log-concave / unimodal? ═══\n");
    for n in 1..=max_n {
        let gamma = gamma_coefficients(&h_polys[n as usize - 1]).unwrap_or_default();
        if gamma.len() < 3 {
            continue;
        }
        let mut lc = true;
        for i in 1..gamma.len() - 1 {
            if gamma[i] * gamma[i] < gamma[i - 1] * gamma[i + 1] {
                lc = false;
                break;
            }
        }
        let mut unimodal = true;
        let mut increasing = true;
        for i in 1..gamma.len() {
            if gamma[i] < gamma[i - 1] {
                increasing = false;
            } else if !increasing && gamma[i] > gamma[i - 1] {
                unimodal = false;
                break;
            }
        }
        println!(
            "n={:>2}: γ = {:?}  log-concave: {}  unimodal: {}",
            n, gamma, lc, unimodal
        );
    }
}
