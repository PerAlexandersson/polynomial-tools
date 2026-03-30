/// Search for recurrences of H_n(t) with various parameter settings.
use combpoly::permutation::alternating_permutations;
use combpoly::statistics::{compute, Stat};
use polynomial_tools::recurrence::{
    find_polynomial_recurrence, find_recurrence_adaptive,
    AdaptiveSearchOptions, RecurrenceOptions,
};

fn compute_polys(max_n: u8) -> Vec<Vec<i64>> {
    let mut polys = Vec::new();
    for n in 1..=max_n {
        let alt = alternating_permutations(n);
        let max_s = alt.iter().map(|s| compute(s, Stat::Swaps)).max().unwrap_or(0);
        let mut coeffs = vec![0i64; max_s + 1];
        for s in &alt {
            coeffs[compute(s, Stat::Swaps)] += 1;
        }
        polys.push(coeffs);
    }
    polys
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(13);

    eprintln!("Computing H_n for n=1..{}", max_n);
    let polys = compute_polys(max_n);
    eprintln!("Done. Searching for recurrences...\n");

    // Try specific small parameter combinations manually
    let configs: Vec<(&str, RecurrenceOptions)> = vec![
        ("depth=2, t^1, n^1, D=0", RecurrenceOptions {
            rec_len: 2, var_deg: 1, idx_deg: 1, diff_deg: 0,
            homogeneous: true, denom_var_deg: 0, denom_idx_deg: 0,
        }),
        ("depth=2, t^1, n^2, D=0", RecurrenceOptions {
            rec_len: 2, var_deg: 1, idx_deg: 2, diff_deg: 0,
            homogeneous: true, denom_var_deg: 0, denom_idx_deg: 0,
        }),
        ("depth=2, t^2, n^2, D=0", RecurrenceOptions {
            rec_len: 2, var_deg: 2, idx_deg: 2, diff_deg: 0,
            homogeneous: true, denom_var_deg: 0, denom_idx_deg: 0,
        }),
        ("depth=2, t^1, n^1, D=1", RecurrenceOptions {
            rec_len: 2, var_deg: 1, idx_deg: 1, diff_deg: 1,
            homogeneous: true, denom_var_deg: 0, denom_idx_deg: 0,
        }),
        ("depth=2, t^2, n^2, D=1", RecurrenceOptions {
            rec_len: 2, var_deg: 2, idx_deg: 2, diff_deg: 1,
            homogeneous: true, denom_var_deg: 0, denom_idx_deg: 0,
        }),
        ("depth=3, t^1, n^1, D=0", RecurrenceOptions {
            rec_len: 3, var_deg: 1, idx_deg: 1, diff_deg: 0,
            homogeneous: true, denom_var_deg: 0, denom_idx_deg: 0,
        }),
        ("depth=3, t^2, n^2, D=0", RecurrenceOptions {
            rec_len: 3, var_deg: 2, idx_deg: 2, diff_deg: 0,
            homogeneous: true, denom_var_deg: 0, denom_idx_deg: 0,
        }),
        ("depth=3, t^1, n^1, D=1", RecurrenceOptions {
            rec_len: 3, var_deg: 1, idx_deg: 1, diff_deg: 1,
            homogeneous: true, denom_var_deg: 0, denom_idx_deg: 0,
        }),
        ("depth=3, t^2, n^1, D=1", RecurrenceOptions {
            rec_len: 3, var_deg: 2, idx_deg: 1, diff_deg: 1,
            homogeneous: true, denom_var_deg: 0, denom_idx_deg: 0,
        }),
        ("depth=3, t^2, n^2, D=1", RecurrenceOptions {
            rec_len: 3, var_deg: 2, idx_deg: 2, diff_deg: 1,
            homogeneous: true, denom_var_deg: 0, denom_idx_deg: 0,
        }),
        // With denominator (rational in n)
        ("depth=2, t^1, n^1, D=1, denom n^1", RecurrenceOptions {
            rec_len: 2, var_deg: 1, idx_deg: 1, diff_deg: 1,
            homogeneous: true, denom_var_deg: 0, denom_idx_deg: 1,
        }),
        ("depth=2, t^2, n^1, D=1, denom n^1", RecurrenceOptions {
            rec_len: 2, var_deg: 2, idx_deg: 1, diff_deg: 1,
            homogeneous: true, denom_var_deg: 0, denom_idx_deg: 1,
        }),
        ("depth=2, t^2, n^2, D=1, denom n^1", RecurrenceOptions {
            rec_len: 2, var_deg: 2, idx_deg: 2, diff_deg: 1,
            homogeneous: true, denom_var_deg: 0, denom_idx_deg: 1,
        }),
        ("depth=2, t^1, n^2, D=1, denom n^1", RecurrenceOptions {
            rec_len: 2, var_deg: 1, idx_deg: 2, diff_deg: 1,
            homogeneous: true, denom_var_deg: 0, denom_idx_deg: 1,
        }),
    ];

    for (name, opts) in &configs {
        let unknowns = (opts.rec_len + 1) * (opts.var_deg + 1) * (opts.idx_deg + 1)
            * (opts.diff_deg + 1);
        print!("{:<45} unknowns~{:<4} ", name, unknowns);
        match find_polynomial_recurrence(&polys, opts) {
            Some(rec) => {
                // Check if coefficients are small (real recurrence) vs huge (overfitting)
                let max_coeff: usize = format!("{}", rec)
                    .len();
                if max_coeff > 200 {
                    println!("FOUND but coefficients too large (overfitting)");
                } else {
                    println!("FOUND: {}", rec);
                }
            }
            None => println!("none"),
        }
    }

    // Also try the adaptive search with high margin
    println!("\n=== Adaptive search with min_margin=3 ===");
    let opts = AdaptiveSearchOptions {
        max_rec_len: 4,
        max_var_deg: 2,
        max_idx_deg: 2,
        max_diff_deg: 1,
        try_inhomogeneous: false,
        try_denominator: true,
        max_denom_var_deg: 1,
        max_denom_idx_deg: 1,
        min_margin: 3,
        verbose: true,
    };
    match find_recurrence_adaptive(&polys, &opts) {
        Some(res) => {
            println!("Found: {}", res.recurrence);
            println!("Params: rec_len={} var_deg={} idx_deg={} diff_deg={}",
                res.opts.rec_len, res.opts.var_deg, res.opts.idx_deg, res.opts.diff_deg);
        }
        None => println!("(no recurrence found)"),
    }
}
