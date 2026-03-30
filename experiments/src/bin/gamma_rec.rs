/// Search for recurrences in the γ-coefficient sequences of H_n(t).
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};

fn main() {
    // γ-coefficients from H_n_data.md, n=1..20
    let gamma: Vec<Vec<i64>> = vec![
        vec![1],
        vec![1],
        vec![1],
        vec![1, 1],
        vec![1, 4],
        vec![1, 10, 5],
        vec![1, 21, 36],
        vec![1, 40, 159, 45],
        vec![1, 72, 556, 528],
        vec![1, 125, 1697, 3612, 665],
        vec![1, 212, 4747, 18920, 11440],
        vec![1, 354, 12516, 84288, 111757, 14457],
        vec![1, 585, 31656, 337072, 817920, 342528],
        vec![1, 960, 77703, 1249182, 5000345, 4523742, 433741],
        vec![1, 1568, 186589, 4378841, 27029020, 44064496, 13534016],
        vec![1, 2553, 440819, 14721672, 133665931, 353768720, 232203115, 17160421],
        vec![1, 4148, 1028820, 47938784, 618365140, 2482366544, 2903667648, 681920256],
        vec![1, 6730, 2379195, 152285588, 2717591717, 15775809859, 29618072105, 14743462160, 865407905],
        vec![1, 10909, 5463913, 474452060, 11472119711, 92961984688, 261641726752, 230170645248, 42664785664],
        vec![1, 17672, 12482147, 1455626336, 46901500041, 516343128806, 2076390817606, 2906867224568, 1134808215161, 54179057649],
    ];

    // Extract individual γ_k sequences (fixed k, varying n)
    let max_k = 9;
    for k in 0..=max_k {
        let seq: Vec<i64> = gamma.iter()
            .filter_map(|g| if g.len() > k { Some(g[k]) } else { None })
            .collect();
        if seq.len() < 6 { continue; }

        println!("=== γ_{} sequence ({} terms): {:?}... ===", k, seq.len(),
            &seq[..seq.len().min(10)]);

        // Try as a plain integer sequence recurrence
        // Convert to degree-0 polynomials
        let poly_seq: Vec<Vec<i64>> = seq.iter().map(|&v| vec![v]).collect();

        for (desc, opts) in &[
            ("len=2 id=2", AdaptiveSearchOptions {
                max_rec_len: 2, max_var_deg: 0, max_idx_deg: 2, max_diff_deg: 0,
                try_inhomogeneous: false, try_denominator: false,
                max_denom_var_deg: 0, max_denom_idx_deg: 0, min_margin: 2, verbose: false }),
            ("len=3 id=2", AdaptiveSearchOptions {
                max_rec_len: 3, max_var_deg: 0, max_idx_deg: 2, max_diff_deg: 0,
                try_inhomogeneous: false, try_denominator: false,
                max_denom_var_deg: 0, max_denom_idx_deg: 0, min_margin: 2, verbose: false }),
            ("len=2 id=3", AdaptiveSearchOptions {
                max_rec_len: 2, max_var_deg: 0, max_idx_deg: 3, max_diff_deg: 0,
                try_inhomogeneous: false, try_denominator: false,
                max_denom_var_deg: 0, max_denom_idx_deg: 0, min_margin: 2, verbose: false }),
            ("len=3 id=3", AdaptiveSearchOptions {
                max_rec_len: 3, max_var_deg: 0, max_idx_deg: 3, max_diff_deg: 0,
                try_inhomogeneous: false, try_denominator: false,
                max_denom_var_deg: 0, max_denom_idx_deg: 0, min_margin: 2, verbose: false }),
            ("len=4 id=2", AdaptiveSearchOptions {
                max_rec_len: 4, max_var_deg: 0, max_idx_deg: 2, max_diff_deg: 0,
                try_inhomogeneous: false, try_denominator: false,
                max_denom_var_deg: 0, max_denom_idx_deg: 0, min_margin: 2, verbose: false }),
            ("len=2 id=4", AdaptiveSearchOptions {
                max_rec_len: 2, max_var_deg: 0, max_idx_deg: 4, max_diff_deg: 0,
                try_inhomogeneous: false, try_denominator: false,
                max_denom_var_deg: 0, max_denom_idx_deg: 0, min_margin: 2, verbose: false }),
        ] {
            let t = std::time::Instant::now();
            match find_recurrence_adaptive(&poly_seq, opts) {
                Some(res) => {
                    let s = format!("{}", res.recurrence);
                    if s.len() > 300 {
                        println!("  {}: overfitted ({} chars) {:?}", desc, s.len(), t.elapsed());
                    } else {
                        println!("  {}: FOUND {:?}\n    {}", desc, t.elapsed(), s);
                    }
                }
                None => println!("  {}: none {:?}", desc, t.elapsed()),
            }
        }
        println!();
    }

    // Also try the γ-vectors as polynomial sequences (treating γ_j as coeff of t^j)
    println!("=== γ-vectors as polynomial sequence ===\n");
    // Skip the first few trivial ones
    let gamma_from4: Vec<Vec<i64>> = gamma[3..].to_vec(); // n=4..20

    for (desc, opts) in &[
        ("len=2 vd=1 id=1 d=0", AdaptiveSearchOptions {
            max_rec_len: 2, max_var_deg: 1, max_idx_deg: 1, max_diff_deg: 0,
            try_inhomogeneous: false, try_denominator: false,
            max_denom_var_deg: 0, max_denom_idx_deg: 0, min_margin: 2, verbose: false }),
        ("len=2 vd=2 id=2 d=0", AdaptiveSearchOptions {
            max_rec_len: 2, max_var_deg: 2, max_idx_deg: 2, max_diff_deg: 0,
            try_inhomogeneous: false, try_denominator: false,
            max_denom_var_deg: 0, max_denom_idx_deg: 0, min_margin: 2, verbose: false }),
        ("len=3 vd=1 id=1 d=0", AdaptiveSearchOptions {
            max_rec_len: 3, max_var_deg: 1, max_idx_deg: 1, max_diff_deg: 0,
            try_inhomogeneous: false, try_denominator: false,
            max_denom_var_deg: 0, max_denom_idx_deg: 0, min_margin: 2, verbose: false }),
        ("len=2 vd=2 id=2 d=1", AdaptiveSearchOptions {
            max_rec_len: 2, max_var_deg: 2, max_idx_deg: 2, max_diff_deg: 1,
            try_inhomogeneous: false, try_denominator: false,
            max_denom_var_deg: 0, max_denom_idx_deg: 0, min_margin: 2, verbose: false }),
        ("len=3 vd=2 id=2 d=0", AdaptiveSearchOptions {
            max_rec_len: 3, max_var_deg: 2, max_idx_deg: 2, max_diff_deg: 0,
            try_inhomogeneous: false, try_denominator: false,
            max_denom_var_deg: 0, max_denom_idx_deg: 0, min_margin: 2, verbose: false }),
        ("len=2 vd=3 id=3 d=0", AdaptiveSearchOptions {
            max_rec_len: 2, max_var_deg: 3, max_idx_deg: 3, max_diff_deg: 0,
            try_inhomogeneous: false, try_denominator: false,
            max_denom_var_deg: 0, max_denom_idx_deg: 0, min_margin: 2, verbose: false }),
        ("len=3 vd=1 id=2 d=1", AdaptiveSearchOptions {
            max_rec_len: 3, max_var_deg: 1, max_idx_deg: 2, max_diff_deg: 1,
            try_inhomogeneous: false, try_denominator: false,
            max_denom_var_deg: 0, max_denom_idx_deg: 0, min_margin: 2, verbose: false }),
        ("len=2 vd=2 id=3 d=1", AdaptiveSearchOptions {
            max_rec_len: 2, max_var_deg: 2, max_idx_deg: 3, max_diff_deg: 1,
            try_inhomogeneous: false, try_denominator: false,
            max_denom_var_deg: 0, max_denom_idx_deg: 0, min_margin: 2, verbose: false }),
    ] {
        let t = std::time::Instant::now();
        match find_recurrence_adaptive(&gamma_from4, opts) {
            Some(res) => {
                let s = format!("{}", res.recurrence);
                if s.len() > 400 {
                    println!("{}: overfitted ({} chars) {:?}", desc, s.len(), t.elapsed());
                } else {
                    println!("{}: FOUND {:?}\n  {}\n", desc, t.elapsed(), s);
                }
            }
            None => println!("{}: none {:?}", desc, t.elapsed()),
        }
    }
}
