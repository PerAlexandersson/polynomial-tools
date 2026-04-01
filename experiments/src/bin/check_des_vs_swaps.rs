// Check if swaps = des on alternating permutations
use combpoly::permutation::alternating_permutations;
use combpoly::statistics::{compute, Stat};

fn main() {
    for n in 3..=8u8 {
        let alt = alternating_permutations(n);
        let mut same = true;
        for s in &alt {
            let sw = compute(s, Stat::Swaps);
            let des = compute(s, Stat::Des);
            if sw != des {
                same = false;
                break;
            }
        }
        // Also compare generating polynomials
        let max_sw = alt
            .iter()
            .map(|s| compute(s, Stat::Swaps))
            .max()
            .unwrap_or(0);
        let max_des = alt.iter().map(|s| compute(s, Stat::Des)).max().unwrap_or(0);
        let mut sw_poly = vec![0i64; max_sw + 1];
        let mut des_poly = vec![0i64; max_des + 1];
        for s in &alt {
            sw_poly[compute(s, Stat::Swaps)] += 1;
            des_poly[compute(s, Stat::Des)] += 1;
        }
        println!(
            "n={}: swaps==des? {}  sw_poly={:?}  des_poly={:?}",
            n, same, sw_poly, des_poly
        );
    }
}
