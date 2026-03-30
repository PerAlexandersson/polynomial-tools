/// Search for recurrences in the peak-nesting polynomial sequences.
///
/// Computes pnest_n(t) and pnest_n^c(t) for n = 0..max_n, then searches for recurrences.
use combpoly::catalan;
use polynomial_tools::real_rootedness as polynomial;
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};
use std::env;
use std::io::{self, Write};

fn main() {
    let max_n: usize = env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(14);

    eprintln!("Computing peak-nesting polynomials for n = 0..{}", max_n);

    let mut pnest_all: Vec<Vec<i64>> = Vec::new();
    let mut pnest_conn: Vec<Vec<i64>> = Vec::new();

    for n in 0..=max_n {
        let (poly, poly_c) = catalan::peak_nesting_polys(n);
        eprintln!(
            "  n={:2}  pnest(t) = {:30}  pnest^c(t) = {}",
            n,
            polynomial::format_poly(&poly),
            polynomial::format_poly(&poly_c),
        );
        io::stderr().flush().ok();
        pnest_all.push(poly);
        pnest_conn.push(poly_c);
    }

    let rational = env::args().any(|a| a == "--rational" || a == "-r");

    if !rational {
        // Phase 1: polynomial recurrence (no denominator)
        let search = AdaptiveSearchOptions {
            max_rec_len: 5,
            max_var_deg: 3,
            max_idx_deg: 3,
            max_diff_deg: 2,
            verbose: true,
            ..Default::default()
        };

        for (name, polys) in [
            ("pnest_n(t) — all Dyck paths", &pnest_all),
            ("pnest_n^c(t) — connected", &pnest_conn),
        ] {
            eprintln!();
            println!("=== {} (n=0..{}) ===", name, max_n);
            match find_recurrence_adaptive(polys, &search) {
                Some(result) => {
                    eprintln!(
                        "Found after {} candidates (unknowns={}, equations={}, \
                         rec_len={}, var_deg={}, idx_deg={}, diff_deg={})",
                        result.candidates_tried,
                        result.num_unknowns,
                        result.num_equations,
                        result.opts.rec_len,
                        result.opts.var_deg,
                        result.opts.idx_deg,
                        result.opts.diff_deg,
                    );
                    println!("{}", result.recurrence);
                }
                None => println!("No recurrence found."),
            }
        }
    } else {
        // Phase 2: rational recurrence (with LHS denominator)
        let search_rat = AdaptiveSearchOptions {
            max_rec_len: 4,
            max_var_deg: 3,
            max_idx_deg: 3,
            max_diff_deg: 2,
            try_denominator: true,
            max_denom_var_deg: 2,
            max_denom_idx_deg: 2,
            verbose: true,
            ..Default::default()
        };

        for (name, polys) in [
            ("pnest_n(t) — all, RATIONAL", &pnest_all),
            ("pnest_n^c(t) — connected, RATIONAL", &pnest_conn),
        ] {
            eprintln!();
            println!("=== {} (n=0..{}) ===", name, max_n);
            match find_recurrence_adaptive(polys, &search_rat) {
                Some(result) => {
                    eprintln!(
                        "Found after {} candidates (unknowns={}, equations={}, \
                         rec_len={}, var_deg={}, idx_deg={}, diff_deg={}, \
                         denom_var={}, denom_idx={})",
                        result.candidates_tried,
                        result.num_unknowns,
                        result.num_equations,
                        result.opts.rec_len,
                        result.opts.var_deg,
                        result.opts.idx_deg,
                        result.opts.diff_deg,
                        result.opts.denom_var_deg,
                        result.opts.denom_idx_deg,
                    );
                    println!("{}", result.recurrence);
                }
                None => println!("No recurrence found."),
            }
        }
    }
}
