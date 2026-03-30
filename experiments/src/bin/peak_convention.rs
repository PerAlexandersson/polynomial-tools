//! Quick convention check: does exact_interlaces(f, g) match the paper's f ≪ g?
//!
//! Paper: f ≪ g means roots of f are more negative than roots of g:
//!   ... ≤ a_2 ≤ b_2 ≤ a_1 ≤ b_1 ≤ 0  (a = f roots, b = g roots, desc index)
//! In ascending arrays: f[0] ≤ g[0] ≤ f[1] ≤ g[1] ≤ ... (same degree)
//!
//! Test with known examples:
//!   f = 1+t (root -1), g = t (root 0).
//!   Paper: f ≪ g iff -1 ≤ 0. TRUE.
//!   Paper: g ≪ f iff 0 ≤ -1. FALSE.

fn exact_interlaces(fc: &[i64], gc: &[i64]) -> bool {
    polynomial_tools::check_weak_interlacing(fc, gc).unwrap_or(false)
}

fn main() {
    // f = 1+t (root -1), g = t (root 0). Same degree 1.
    // Paper: f ≪ g iff root(f) ≤ root(g), i.e., -1 ≤ 0. TRUE.
    let f = vec![1, 1]; // 1 + t
    let g = vec![0, 1]; // t

    println!("f = 1+t (root -1), g = t (root 0)");
    println!("exact_interlaces(f, g) = {}", exact_interlaces(&f, &g));
    println!("exact_interlaces(g, f) = {}", exact_interlaces(&g, &f));
    println!("Paper: f ≪ g should be TRUE (since -1 ≤ 0)");
    println!("Paper: g ≪ f should be FALSE (since 0 > -1)");

    println!();

    // f = 2+t (root -2), g = 1+t (root -1). Same degree 1.
    // Paper: f ≪ g iff -2 ≤ -1. TRUE.
    let f2 = vec![2, 1]; // 2 + t
    let g2 = vec![1, 1]; // 1 + t
    println!("f = 2+t (root -2), g = 1+t (root -1)");
    println!("exact_interlaces(f, g) = {}", exact_interlaces(&f2, &g2));
    println!("exact_interlaces(g, f) = {}", exact_interlaces(&g2, &f2));
    println!("Paper: f ≪ g should be TRUE (since -2 ≤ -1)");

    println!();

    // f = 1 (constant), g = 1+t. deg(f)+1 = deg(g).
    // Paper: f ≪ g iff g is RR with non-pos roots and f has pos leading coeff. TRUE.
    let f3 = vec![1];
    let g3 = vec![1, 1];
    println!("f = 1 (constant), g = 1+t (root -1)");
    println!("exact_interlaces(f, g) = {}", exact_interlaces(&f3, &g3));
    println!("exact_interlaces(g, f) = {}", exact_interlaces(&g3, &f3));
    println!("Paper: f ≪ g should be TRUE (const ≪ linear with neg root)");

    println!();

    // Same degree 2: f = (t+1)(t+4) = 4+5t+t^2, g = (t+2)(t+3) = 6+5t+t^2
    // Roots: f has -4, -1; g has -3, -2.
    // Paper: f ≪ g iff -4 ≤ -3 ≤ -1... wait need full interleaving.
    // Ascending: f[-4, -1], g[-3, -2].
    // Paper: f[0] ≤ g[0] ≤ f[1] ≤ g[1]: -4 ≤ -3 ≤ -1 ≤ -2? NO (-1 > -2).
    // So f ≪ g should be FALSE.
    // Paper: g ≪ f: g[0] ≤ f[0]? -3 ≤ -4? NO. So g ≪ f also FALSE.
    // Neither interlaces!
    let f4 = vec![4, 5, 1]; // (t+1)(t+4)
    let g4 = vec![6, 5, 1]; // (t+2)(t+3)
    println!("f = (t+1)(t+4) roots [-4,-1], g = (t+2)(t+3) roots [-3,-2]");
    println!("exact_interlaces(f, g) = {}", exact_interlaces(&f4, &g4));
    println!("exact_interlaces(g, f) = {}", exact_interlaces(&g4, &f4));
    println!("Paper: NEITHER f≪g nor g≪f (roots don't interleave)");

    println!();

    // f = (t+1)(t+3) = 3+4t+t^2, g = (t+2)(t+4) = 8+6t+t^2
    // Roots: f[-3,-1], g[-4,-2]. Ascending same.
    // Paper f ≪ g: f[0]≤g[0]? -3≤-4? NO. FALSE.
    // Paper g ≪ f: g[0]≤f[0]? -4≤-3 ✓. f[0]≤g[1]? -3≤-2 ✓. g[1]≤f[1]? -2≤-1 ✓. TRUE.
    let f5 = vec![3, 4, 1]; // (t+1)(t+3)
    let g5 = vec![8, 6, 1]; // (t+2)(t+4)
    println!("f = (t+1)(t+3) roots [-3,-1], g = (t+2)(t+4) roots [-4,-2]");
    println!("exact_interlaces(f, g) = {}", exact_interlaces(&f5, &g5));
    println!("exact_interlaces(g, f) = {}", exact_interlaces(&g5, &f5));
    println!("Paper: f≪g FALSE, g≪f TRUE");
}
