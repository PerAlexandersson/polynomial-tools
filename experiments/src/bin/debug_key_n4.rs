use combinatoric_core::key_polynomial::*;
use combinatoric_core::partition::Partition;

fn main() {
    let lambda = Partition::new(vec![2, 1]);
    let sigma: Vec<usize> = vec![2, 4, 3, 1];
    let n = 4;

    let w0: Vec<usize> = (1..=n).rev().collect();
    let w0_sigma: Vec<usize> = w0.iter().map(|&i| sigma[i - 1]).collect();
    println!("w0 = {:?}", w0);
    println!("sigma = {:?}", sigma);
    println!("w0_sigma = {:?}", w0_sigma);

    let faces = reduced_kogan_faces(n, &w0_sigma);
    println!(
        "\nReduced Kogan faces of type {:?}: {} faces",
        w0_sigma,
        faces.len()
    );
    for (idx, eqs) in faces.iter().enumerate() {
        let pats = gt_patterns_with_equalities_n(&lambda, n, eqs);
        println!("  Face {}: eqs={:?}, {} patterns", idx, eqs, pats.len());
        for p in &pats {
            println!("    rows: {:?} weight: {:?}", p.rows, p.weight());
        }
    }

    let weights = key_polynomial_weights(&lambda, &sigma);
    println!("\nTotal: {} GT patterns", weights.values().sum::<u64>());
    for (w, c) in &weights {
        println!("  {:?}: {}", w, c);
    }
}
