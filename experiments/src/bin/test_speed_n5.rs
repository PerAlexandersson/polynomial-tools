use combinatoric_core::key_polynomial::key_ehrhart_polynomial;
use combinatoric_core::partition::Partition;
use std::time::Instant;

fn main() {
    let lambda = Partition::new(vec![4, 3, 2, 1]);
    // Try w_0 = [5,4,3,2,1] (hardest case)
    let t = Instant::now();
    let ep = key_ehrhart_polynomial(&lambda, &[5, 4, 3, 2, 1], None);
    println!("w0: {} ({:.2?})", ep.display(), t.elapsed());

    // Try a mid-length perm
    let t = Instant::now();
    let ep = key_ehrhart_polynomial(&lambda, &[3, 1, 5, 2, 4], None);
    println!("31524: {} ({:.2?})", ep.display(), t.elapsed());

    // Try identity
    let t = Instant::now();
    let ep = key_ehrhart_polynomial(&lambda, &[1, 2, 3, 4, 5], None);
    println!("id: {} ({:.2?})", ep.display(), t.elapsed());
}
