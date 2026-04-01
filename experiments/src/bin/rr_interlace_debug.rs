/// Debug the position-refined interlacing.
/// The chain H_{n,2} ≪ H_{n,4} ≪ ... ≪ H_{n,n} fails in practice.
/// Try alternative interlacing notions and orderings.
use combpoly::permutation::alternating_permutations;
use combpoly::statistics::{compute, Stat};
use polynomial_tools::real_rootedness::{
    check_interlacing, check_interlacing_sturm, check_weak_interlacing, format_poly, is_real_rooted,
};

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

fn poly_add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut c = vec![0i64; len];
    for (i, &v) in a.iter().enumerate() {
        c[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        c[i] += v;
    }
    c
}

fn interlace_str(f: &[i64], g: &[i64]) -> &'static str {
    if f.len() <= 1 || g.len() <= 1 {
        return "trivial";
    }
    // Try both orderings with both methods
    let bezout_fg = check_interlacing(f, g);
    let bezout_gf = check_interlacing(g, f);
    let weak_fg = check_weak_interlacing(f, g);
    let weak_gf = check_weak_interlacing(g, f);

    if bezout_fg == Some(true) || bezout_gf == Some(true) {
        "strict ✓"
    } else if weak_fg == Some(true) || weak_gf == Some(true) {
        "weak ✓"
    } else {
        "✗"
    }
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(13);

    for n in (4..=max_n).step_by(1) {
        let alt = alternating_permutations(n);

        let mut by_pos: std::collections::BTreeMap<usize, Vec<&Vec<u8>>> =
            std::collections::BTreeMap::new();
        for s in &alt {
            let pos = s.iter().position(|&v| v == n).unwrap() + 1;
            by_pos.entry(pos).or_default().push(s);
        }

        let positions: Vec<usize> = by_pos.keys().copied().collect();
        if positions.len() < 2 {
            continue;
        }
        let pos_polys: Vec<(usize, Vec<i64>)> = positions
            .iter()
            .map(|&p| (p, build_poly(by_pos.get(&p).unwrap())))
            .collect();

        println!("n={}:", n);

        // Check all pairs for interlacing (both orderings)
        println!("  Pairwise interlacing:");
        for i in 0..pos_polys.len() {
            for j in i + 1..pos_polys.len() {
                let (pi, fi) = &pos_polys[i];
                let (pj, fj) = &pos_polys[j];
                let status = interlace_str(fi, fj);
                if status != "✗" {
                    println!("    H_{{{},{}}} ~ H_{{{},{}}}: {}", n, pi, n, pj, status);
                }
            }
        }

        // Check reversed order: H_{n,n} ≪ H_{n,n-2} ≪ ... ≪ H_{n,2}
        print!("  Reverse chain (large→small pos): ");
        let mut rev_ok = true;
        for i in (0..pos_polys.len() - 1).rev() {
            let f = &pos_polys[i + 1].1;
            let g = &pos_polys[i].1;
            let weak = check_weak_interlacing(
                if f.len() <= g.len() { f } else { g },
                if f.len() <= g.len() { g } else { f },
            );
            if weak != Some(true) {
                rev_ok = false;
                break;
            }
        }
        println!("{}", if rev_ok { "✓" } else { "✗" });

        // Check cumulative sums from the LEFT: C_k = H_{n,2} + ... + H_{n,k}
        print!("  Cumulative L→R chain (C_2 ≪ C_4 ≪ ... ≪ H_n): ");
        let mut cum_lr = vec![0i64; 1];
        let mut cum_lr_polys = Vec::new();
        for (_, poly) in &pos_polys {
            cum_lr = poly_add(&cum_lr, poly);
            cum_lr_polys.push(cum_lr.clone());
        }
        let mut lr_ok = true;
        for i in 0..cum_lr_polys.len() - 1 {
            let f = &cum_lr_polys[i];
            let g = &cum_lr_polys[i + 1];
            let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
            if small.len() <= 1 {
                continue;
            }
            let weak = check_weak_interlacing(small, large);
            if weak != Some(true) {
                lr_ok = false;
                break;
            }
        }
        println!("{}", if lr_ok { "✓" } else { "✗" });

        // Check cumulative sums from the RIGHT: S_k = H_{n,k} + H_{n,k+2} + ... + H_{n,n}
        print!("  Cumulative R→L chain (S_n ≪ S_{{n-2}} ≪ ... ≪ H_n): ");
        let mut cum_rl = vec![0i64; 1];
        let mut cum_rl_polys = Vec::new();
        for (_, poly) in pos_polys.iter().rev() {
            cum_rl = poly_add(&cum_rl, poly);
            cum_rl_polys.push(cum_rl.clone());
        }
        cum_rl_polys.reverse();
        let mut rl_ok = true;
        for i in 0..cum_rl_polys.len() - 1 {
            // S_{n,k+2} ≪ S_{n,k}
            let f = &cum_rl_polys[i + 1]; // smaller cumulative
            let g = &cum_rl_polys[i]; // larger cumulative
            let (small, large) = if f.len() <= g.len() { (f, g) } else { (g, f) };
            if small.len() <= 1 {
                continue;
            }
            let weak = check_weak_interlacing(small, large);
            if weak != Some(true) {
                rl_ok = false;
                break;
            }
        }
        println!("{}", if rl_ok { "✓" } else { "✗" });

        // Check: does each H_{n,j} individually interlace H_n?
        let h_n = &cum_lr_polys.last().unwrap();
        print!("  All H_{{n,j}} interlace H_n: ");
        let mut all_int_hn = true;
        for (_, poly) in &pos_polys {
            let (small, large) = if poly.len() <= h_n.len() {
                (poly.as_slice(), h_n.as_slice())
            } else {
                (h_n.as_slice(), poly.as_slice())
            };
            if small.len() <= 1 {
                continue;
            }
            let weak = check_weak_interlacing(small, large);
            if weak != Some(true) {
                all_int_hn = false;
                break;
            }
        }
        println!("{}", if all_int_hn { "✓" } else { "✗" });

        // Check: H_{n-1} interlaces H_n
        if n >= 5 {
            let alt_prev = alternating_permutations(n - 1);
            let refs_prev: Vec<&Vec<u8>> = alt_prev.iter().collect();
            let h_prev = build_poly(&refs_prev);
            let (small, large) = if h_prev.len() <= h_n.len() {
                (h_prev.as_slice(), h_n.as_slice())
            } else {
                (h_n.as_slice(), h_prev.as_slice())
            };
            let int = check_weak_interlacing(small, large);
            println!(
                "  H_{} ≪ H_{}: {}",
                n - 1,
                n,
                if int == Some(true) { "✓" } else { "✗" }
            );
        }

        println!();
    }
}
