/// Benchmark: compute chromatic symmetric functions for all area sequences
/// of length n, expanded in the elementary basis.
use combinatoric_core::graph::Graph;
use combpoly::catalan::all_area_sequences;
use std::time::Instant;
use sym_poly_sym::chromatic_symmetric;

fn main() {
    for n in 5..=9 {
        let areas = all_area_sequences(n);
        let count = areas.len();
        println!("n={}: {} area sequences", n, count);

        let t = Instant::now();
        let mut e_pos_count = 0usize;
        let mut total = 0usize;

        for area in &areas {
            let g = Graph::unit_interval(area);
            let x = chromatic_symmetric::<i64>(&g);
            let e = x.to_elementary_basis();
            total += 1;

            // Check e-positivity (all coefficients non-negative)
            let all_nonneg = e.terms().values().all(|&c| c >= 0);
            if all_nonneg {
                e_pos_count += 1;
            }
        }

        let elapsed = t.elapsed();
        println!(
            "  computed {} chromatic symm funcs in e-basis: {:.2?}",
            total, elapsed
        );
        println!(
            "  e-positive: {}/{} ({:.1}%)",
            e_pos_count,
            total,
            100.0 * e_pos_count as f64 / total as f64
        );
        println!("  avg per graph: {:.2?}", elapsed / total as u32);
        println!();
    }
}
// Note: extended to n=8,9 with the multiset permutation approach
