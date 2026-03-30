//! Verify directed interlacing works correctly after the fix.
use polynomial_tools::real_rootedness::check_weak_interlacing;

fn main() {
    // f = t+3 (root -3), g = t+1 (root -1).
    // f ≪ g: -3 ≤ -1 ✓ (f on left)
    // g ≪ f: -1 ≤ -3 ✗ (g NOT on left)
    let f = vec![3, 1]; // t + 3
    let g = vec![1, 1]; // t + 1
    println!("f = t+3, g = t+1");
    println!("  f ≪ g (should be true):  {:?}", check_weak_interlacing(&f, &g));
    println!("  g ≪ f (should be false): {:?}", check_weak_interlacing(&g, &f));

    // Same degree: f = (t+4)(t+2) = t²+6t+8, g = (t+3)(t+1) = t²+4t+3
    // Roots of f: -4, -2. Roots of g: -3, -1.
    // f ≪ g: -4 ≤ -3 ≤ -2 ≤ -1 ✓
    // g ≪ f: -3 ≤ -4? ✗
    let f2 = vec![8, 6, 1]; // t²+6t+8
    let g2 = vec![3, 4, 1]; // t²+4t+3
    println!("\nf = t²+6t+8 (roots -4,-2), g = t²+4t+3 (roots -3,-1)");
    println!("  f ≪ g (should be true):  {:?}", check_weak_interlacing(&f2, &g2));
    println!("  g ≪ f (should be false): {:?}", check_weak_interlacing(&g2, &f2));

    // deg(f) = deg(g) + 1: f = (t+3)(t+1) = t²+4t+3, g = t+2
    // f ≪ g: -3 ≤ -2 ≤ -1 ✓
    // g ≪ f: should be None (deg(g) < deg(f), wrong for g ≪ f)
    let f3 = vec![3, 4, 1]; // t²+4t+3
    let g3 = vec![2, 1];    // t+2
    println!("\nf = t²+4t+3 (roots -3,-1), g = t+2 (root -2)");
    println!("  f ≪ g (should be None, deg(f) > deg(g)): {:?}", check_weak_interlacing(&f3, &g3));
    println!("  g ≪ f (should be true):  {:?}", check_weak_interlacing(&g3, &f3));

    // Eulerian: A_3 = 1+4t+t², reversed with t·A_3 = t+4t²+t³
    // 1+4t+t² ≪ t+4t²+t³? Roots: ≈-3.73,-0.27 vs 0,-3.73,-0.27
    // Pattern: -3.73 ≤ -3.73 ≤ -0.27 ≤ -0.27 ≤ 0 ✓
    let a3 = vec![1, 4, 1];    // 1+4t+t²
    let ta3 = vec![0, 1, 4, 1]; // t+4t²+t³
    println!("\nA_3 = 1+4t+t², t·A_3 = t+4t²+t³");
    println!("  A_3 ≪ t·A_3 (should be true): {:?}", check_weak_interlacing(&a3, &ta3));
    println!("  t·A_3 ≪ A_3 (should be None): {:?}", check_weak_interlacing(&ta3, &a3));
}
