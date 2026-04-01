/// Search for recurrence of Family III (Av(231) long swaps).
use combpoly::permutation::avoiding_permutations;
use combpoly::statistics::{compute, Stat};
use polynomial_tools::recurrence::{
    find_polynomial_recurrence, find_recurrence_adaptive, AdaptiveSearchOptions, RecurrenceOptions,
};

fn build_poly(perms: &[Vec<u8>]) -> Vec<i64> {
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
    let max_n = 14u8;
    let pattern = vec![2u8, 3, 1]; // Av(231)

    eprintln!("Computing Family III (Av(231)) to n={}...", max_n);
    let polys: Vec<Vec<i64>> = (1..=max_n)
        .map(|n| build_poly(&avoiding_permutations(n, &[pattern.clone()])))
        .collect();
    eprintln!("Done.\n");

    // Try specific configs
    let configs: Vec<(&str, RecurrenceOptions)> = vec![
        (
            "d=2 t1 n1 D0",
            RecurrenceOptions {
                rec_len: 2,
                var_deg: 1,
                idx_deg: 1,
                diff_deg: 0,
                homogeneous: true,
                denom_var_deg: 0,
                denom_idx_deg: 0,
            },
        ),
        (
            "d=2 t1 n2 D0",
            RecurrenceOptions {
                rec_len: 2,
                var_deg: 1,
                idx_deg: 2,
                diff_deg: 0,
                homogeneous: true,
                denom_var_deg: 0,
                denom_idx_deg: 0,
            },
        ),
        (
            "d=2 t2 n2 D0",
            RecurrenceOptions {
                rec_len: 2,
                var_deg: 2,
                idx_deg: 2,
                diff_deg: 0,
                homogeneous: true,
                denom_var_deg: 0,
                denom_idx_deg: 0,
            },
        ),
        (
            "d=2 t1 n1 D1",
            RecurrenceOptions {
                rec_len: 2,
                var_deg: 1,
                idx_deg: 1,
                diff_deg: 1,
                homogeneous: true,
                denom_var_deg: 0,
                denom_idx_deg: 0,
            },
        ),
        (
            "d=2 t2 n1 D1",
            RecurrenceOptions {
                rec_len: 2,
                var_deg: 2,
                idx_deg: 1,
                diff_deg: 1,
                homogeneous: true,
                denom_var_deg: 0,
                denom_idx_deg: 0,
            },
        ),
        (
            "d=2 t1 n2 D1",
            RecurrenceOptions {
                rec_len: 2,
                var_deg: 1,
                idx_deg: 2,
                diff_deg: 1,
                homogeneous: true,
                denom_var_deg: 0,
                denom_idx_deg: 0,
            },
        ),
        // With denominator
        (
            "d=2 t1 n1 D1 dn1",
            RecurrenceOptions {
                rec_len: 2,
                var_deg: 1,
                idx_deg: 1,
                diff_deg: 1,
                homogeneous: true,
                denom_var_deg: 0,
                denom_idx_deg: 1,
            },
        ),
        (
            "d=2 t2 n1 D1 dn1",
            RecurrenceOptions {
                rec_len: 2,
                var_deg: 2,
                idx_deg: 1,
                diff_deg: 1,
                homogeneous: true,
                denom_var_deg: 0,
                denom_idx_deg: 1,
            },
        ),
        (
            "d=2 t1 n2 D1 dn1",
            RecurrenceOptions {
                rec_len: 2,
                var_deg: 1,
                idx_deg: 2,
                diff_deg: 1,
                homogeneous: true,
                denom_var_deg: 0,
                denom_idx_deg: 1,
            },
        ),
        (
            "d=2 t2 n2 D1 dn1",
            RecurrenceOptions {
                rec_len: 2,
                var_deg: 2,
                idx_deg: 2,
                diff_deg: 1,
                homogeneous: true,
                denom_var_deg: 0,
                denom_idx_deg: 1,
            },
        ),
    ];

    for (name, opts) in &configs {
        print!("{:<25} ", name);
        match find_polynomial_recurrence(&polys, opts) {
            Some(rec) => {
                let s = format!("{}", rec);
                if s.len() > 300 {
                    println!("overfitted ({} chars)", s.len());
                } else {
                    println!("FOUND: {}", s);
                }
            }
            None => println!("none"),
        }
    }

    // Adaptive with decent margin
    println!("\n=== Adaptive (margin=2) ===");
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
        verbose: true,
    };
    match find_recurrence_adaptive(&polys, &opts) {
        Some(res) => {
            let s = format!("{}", res.recurrence);
            if s.len() > 300 {
                println!("overfitted ({} chars)", s.len());
            } else {
                println!("FOUND: {}", s);
            }
            println!(
                "Params: rec_len={} var_deg={} idx_deg={} diff_deg={} denom=({},{})",
                res.opts.rec_len,
                res.opts.var_deg,
                res.opts.idx_deg,
                res.opts.diff_deg,
                res.opts.denom_var_deg,
                res.opts.denom_idx_deg
            );
        }
        None => println!("(none)"),
    }
}
