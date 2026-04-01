use combinatoric_core::key_polynomial::*;
use combinatoric_core::partition::Partition;

fn main() {
    let lambda = Partition::new(vec![2, 1]);
    let patterns = gt_patterns_with_equalities_n(&lambda, 3, &[]);
    println!("GT patterns for (2,1,0), n=3: {} total", patterns.len());
    for p in &patterns {
        let w = p.weight();
        println!("  rows: {:?}  weight: {:?}", p.rows, w);
    }
}
