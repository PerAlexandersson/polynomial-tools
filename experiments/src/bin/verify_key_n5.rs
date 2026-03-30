use combinatoric_core::partition::Partition;
use combinatoric_core::key_polynomial::key_polynomial_weights;
use std::time::Instant;

fn main() {
    let lambda = Partition::new(vec![3, 2, 1]);
    let perms: Vec<Vec<usize>> = vec![
        vec![1,2,3,4,5],
        vec![2,1,3,4,5],
        vec![1,3,2,4,5],
        vec![3,2,1,4,5],
        vec![2,1,4,3,5],
        vec![5,4,3,2,1],
        vec![3,1,4,2,5],
        vec![2,4,1,5,3],
    ];

    println!("{:<16} {:>10} {:>12} {:>8}", "sigma", "monomials", "totalCoeff", "time");
    println!("{}", "-".repeat(50));

    for sigma in &perms {
        let t = Instant::now();
        let w = key_polynomial_weights(&lambda, sigma);
        let elapsed = t.elapsed();
        let monomials = w.len();
        let total_coeff: u64 = w.values().sum();
        println!(
            "{:<16} {:>10} {:>12} {:>6.1?}",
            format!("{:?}", sigma), monomials, total_coeff, elapsed
        );
    }
}
