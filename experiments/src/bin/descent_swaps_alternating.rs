/// Focus on the ALTERNATING case: S = {2, 4, 6, ...} ∩ [n-1].
///
/// Special properties:
/// - H_{n,n} = t · H_{n-1} (max at last position)
/// - Palindromic symmetry H_{n,j} ↔ H_{n,n+2-j}
/// - Complement σ → σ̄ maps pos(max=n)=j to pos(min=1)=n+1-j
///
/// Can we find a MATRIX RECURRENCE relating (H_{n,j})_j to (H_{n-1,j})_j
/// that satisfies Brändén's condition?
///
/// Also: check the cumulative interlacing chain
/// C_{n,2} ≺ C_{n,4} ≺ ... ≺ C_{n,n} = H_n
/// and see if this can be proved from the recurrence.
///
use combpoly::permutation::alternating_permutations;
use combpoly::statistics::{compute, Stat};
use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly, is_real_rooted};

fn build_poly(vals: &[usize]) -> Vec<i64> {
    if vals.is_empty() {
        return vec![0];
    }
    let max_s = *vals.iter().max().unwrap();
    let mut coeffs = vec![0i64; max_s + 1];
    for &s in vals {
        coeffs[s] += 1;
    }
    while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
        coeffs.pop();
    }
    coeffs
}

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &c) in a.iter().enumerate() {
        r[i] += c;
    }
    for (i, &c) in b.iter().enumerate() {
        r[i] += c;
    }
    while r.len() > 1 && *r.last().unwrap() == 0 {
        r.pop();
    }
    r
}

fn check_compat(f: &[i64], g: &[i64]) -> bool {
    if f.len() <= 1 || g.len() <= 1 {
        return true;
    }
    if !is_real_rooted(f) || !is_real_rooted(g) {
        return false;
    }
    let weights: Vec<(i64, i64)> = vec![
        (1, 1),
        (1, 2),
        (2, 1),
        (1, 3),
        (3, 1),
        (1, 5),
        (5, 1),
        (2, 3),
        (3, 2),
        (1, 7),
        (7, 1),
    ];
    let maxlen = f.len().max(g.len());
    for (a, b) in &weights {
        let mut c = vec![0i64; maxlen];
        for (i, &v) in f.iter().enumerate() {
            c[i] += a * v;
        }
        for (i, &v) in g.iter().enumerate() {
            c[i] += b * v;
        }
        while c.len() > 1 && *c.last().unwrap() == 0 {
            c.pop();
        }
        if c.len() > 1 && !is_real_rooted(&c) {
            return false;
        }
    }
    true
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(14);

    println!("=== Alternating case: detailed analysis ===\n");

    let mut prev_h: Option<Vec<i64>> = None;

    for n in 1..=max_n {
        let alt = alternating_permutations(n);
        if alt.is_empty() {
            continue;
        }

        // Group by position of max (value n)
        let mut by_pos: Vec<Vec<usize>> = vec![vec![]; n as usize];
        for sigma in &alt {
            let pos = sigma.iter().position(|&v| v == n).unwrap();
            by_pos[pos].push(compute(sigma, Stat::Swaps));
        }

        let h_polys: Vec<(u8, Vec<i64>)> = by_pos
            .iter()
            .enumerate()
            .filter(|(_, v)| !v.is_empty())
            .map(|(pos, vals)| ((pos + 1) as u8, build_poly(vals)))
            .collect();

        let h_n: Vec<i64> = {
            let all_vals: Vec<usize> = alt.iter().map(|s| compute(s, Stat::Swaps)).collect();
            build_poly(&all_vals)
        };

        println!(
            "═══ n = {} ({} alt perms, {} positions) ═══",
            n,
            alt.len(),
            h_polys.len()
        );
        println!("H_{}(t) = {}", n, format_poly(&h_n));

        for (j, poly) in &h_polys {
            println!("  H_{},{} = {}", n, j, format_poly(poly));
        }

        // Check: H_{n,n} = t · H_{n-1}?
        if let Some(ref ph) = prev_h {
            let t_prev: Vec<i64> = {
                let mut v = vec![0i64; ph.len() + 1];
                for (i, &c) in ph.iter().enumerate() {
                    v[i + 1] = c;
                }
                v
            };
            let last_poly = &h_polys.last().unwrap().1;
            let last_pos = h_polys.last().unwrap().0;
            if *last_poly == t_prev {
                println!("  H_{},{} = t·H_{} ✓", n, last_pos, n - 1);
            } else {
                println!("  H_{},{} ≠ t·H_{} ✗", n, last_pos, n - 1);
            }
        }

        // Check H_{n-1} interlaces H_n
        if let Some(ref ph) = prev_h {
            if ph.len() > 1 && h_n.len() > 1 {
                let (small, large) = if ph.len() <= h_n.len() {
                    (ph, &h_n)
                } else {
                    (&h_n, ph)
                };
                match check_weak_interlacing(small, large) {
                    Some(true) => println!("  H_{} ≺ H_{} ✓", n - 1, n),
                    Some(false) => println!("  H_{} ≺ H_{} ✗ FAIL", n - 1, n),
                    None => println!("  H_{} ≺ H_{} ?", n - 1, n),
                }
            }
        }

        // Compatible family check
        if h_polys.len() >= 2 {
            let mut all_compat = true;
            for i in 0..h_polys.len() {
                for j in (i + 1)..h_polys.len() {
                    if !check_compat(&h_polys[i].1, &h_polys[j].1) {
                        all_compat = false;
                        println!(
                            "  COMPAT FAIL: H_{},{} vs H_{},{}",
                            n, h_polys[i].0, n, h_polys[j].0
                        );
                    }
                }
            }
            if all_compat {
                println!("  Compatible family ✓");
            }
        }

        // Cumulative chain: C_{n,2} ≺ C_{n,4} ≺ ... ≺ C_{n,n}
        if h_polys.len() >= 2 {
            let mut cum = vec![0i64; 1];
            let mut cum_polys: Vec<(u8, Vec<i64>)> = Vec::new();
            for (j, poly) in &h_polys {
                cum = poly_add(&cum, poly);
                cum_polys.push((*j, cum.clone()));
            }

            let mut chain_ok = true;
            for i in 0..cum_polys.len() - 1 {
                let (j1, ref c1) = cum_polys[i];
                let (j2, ref c2) = cum_polys[i + 1];
                if c1.len() > 1 && c2.len() > 1 {
                    let (s, l) = if c1.len() <= c2.len() {
                        (c1, c2)
                    } else {
                        (c2, c1)
                    };
                    match check_weak_interlacing(s, l) {
                        Some(true) => {}
                        _ => {
                            chain_ok = false;
                            println!("  Chain FAIL: C_{},{} ≺ C_{},{}", n, j1, n, j2);
                        }
                    }
                }
            }
            if chain_ok {
                println!("  Cumulative chain ✓");
            }
        }

        // Transfer matrix: express H_{n,j} in terms of the H_{n-1,k}
        // For p=n (last position): H_{n,n} = t · H_{n-1} (already checked)
        // For other positions: H_{n,j} involves OTHER descent classes
        //
        // Key: what fraction of H_{n,j} comes from the alternating source?
        if n >= 4 && h_polys.len() >= 2 {
            let alt_prev = alternating_permutations(n - 1);
            let total_prev = alt_prev.len();

            // For each valid position j in H_n:
            // Count how many source permutations come from A_{n-1}
            // vs other descent classes
            for (j, _) in &h_polys {
                let p = *j;
                // Count sources from alternating class
                // Source S'_p: those in S_{n-1} that insert at p giving alt descent set
                // For the alternating case, S'_p depends on p.
                // S'_n = {2,4,...,n-2} = alt descent set for n-1
                // S'_p for p < n involves non-alternating descent sets

                // Let's just count: of the target permutations σ with max at j,
                // remove the max to get π. How many π are alternating?
                let count_alt_source: usize = alt
                    .iter()
                    .filter(|sigma| sigma.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                    .filter(|sigma| {
                        // Remove n from position p, get π
                        let mut pi: Vec<u8> = sigma.iter().filter(|&&v| v != n).copied().collect();
                        // Check if π is alternating (up-down)
                        if pi.len() < 2 {
                            return true;
                        }
                        for i in 0..pi.len() - 1 {
                            if i % 2 == 0 && pi[i] >= pi[i + 1] {
                                return false;
                            }
                            if i % 2 == 1 && pi[i] <= pi[i + 1] {
                                return false;
                            }
                        }
                        true
                    })
                    .count();

                let total_at_j: usize = alt
                    .iter()
                    .filter(|sigma| sigma.iter().position(|&v| v == n).unwrap() as u8 + 1 == p)
                    .count();

                if total_at_j > 0 {
                    let frac = count_alt_source as f64 / total_at_j as f64;
                    if p == n {
                        // Should be 100% from alternating source
                        print!(
                            "  p={}: {}/{} from A_{{n-1}} ({:.0}%)",
                            p,
                            count_alt_source,
                            total_at_j,
                            frac * 100.0
                        );
                        if frac == 1.0 {
                            println!(" ✓");
                        } else {
                            println!(" ✗");
                        }
                    } else if frac > 0.0 {
                        println!(
                            "  p={}: {}/{} from A_{{n-1}} ({:.0}%)",
                            p,
                            count_alt_source,
                            total_at_j,
                            frac * 100.0
                        );
                    }
                }
            }
        }

        prev_h = Some(h_n);
        println!();
    }
}
