//! Test stability of F(s,t) = perm(M(s)) for staircase matrices.
//! M(s) is M_{μ'} with one entry replaced by s.
//! 
//! Stability means: F(s,t) ≠ 0 whenever Im(s) > 0 and Im(t) > 0.
//! 
//! Test by evaluating F at many points (s,t) in the upper half-plane
//! and checking if it's ever zero (or very close to zero).
//!
//! Also test: for RANDOM staircase matrices (not just from Ferrers boards),
//! does the permanent remain nonzero in the upper half-plane?
//!
//! And test the COUNTEREXAMPLE structure: does perm([[1,s],[t,1]]) = 1+st
//! fail, while staircase perm([[s,t],[1,1]]) = s+t succeed?

use std::collections::BTreeSet;

#[derive(Clone, Copy, Debug)]
struct C64 { re: f64, im: f64 }
impl C64 {
    fn new(re: f64, im: f64) -> Self { Self { re, im } }
    fn mul(self, other: Self) -> Self {
        C64::new(self.re * other.re - self.im * other.im,
                 self.re * other.im + self.im * other.re)
    }
    fn add(self, other: Self) -> Self {
        C64::new(self.re + other.re, self.im + other.im)
    }
    fn norm2(self) -> f64 { self.re * self.re + self.im * self.im }
}

/// Evaluate permanent of an n×n complex matrix (brute force, n ≤ 8).
fn perm_complex(mat: &[Vec<C64>]) -> C64 {
    let n = mat.len();
    if n == 0 { return C64::new(1.0, 0.0); }
    let mut result = C64::new(0.0, 0.0);
    // Generate all permutations
    let mut indices: Vec<usize> = (0..n).collect();
    loop {
        let mut prod = C64::new(1.0, 0.0);
        for i in 0..n { prod = prod.mul(mat[i][indices[i]]); }
        result = result.add(prod);
        if !next_perm(&mut indices) { break; }
    }
    result
}

fn next_perm(a: &mut [usize]) -> bool {
    let n = a.len();
    if n <= 1 { return false; }
    let mut i = n - 2;
    while a[i] >= a[i + 1] { if i == 0 { return false; } i -= 1; }
    let mut j = n - 1;
    while a[j] <= a[i] { j -= 1; }
    a.swap(i, j);
    a[i + 1..].reverse();
    true
}

fn main() {
    println!("=== Stability tests for staircase permanent F(s,t) ===\n");
    
    // Test points in the upper half-plane
    let test_points: Vec<C64> = vec![
        C64::new(0.0, 1.0),    // i
        C64::new(1.0, 1.0),    // 1+i
        C64::new(-1.0, 1.0),   // -1+i
        C64::new(0.0, 0.1),    // 0.1i
        C64::new(0.0, 10.0),   // 10i
        C64::new(-0.5, 0.5),
        C64::new(2.0, 0.3),
        C64::new(-3.0, 0.01),
        C64::new(0.5, 2.0),
        C64::new(-0.1, 0.1),
    ];
    
    // Test 1: Known failure - non-staircase [[1,s],[t,1]]
    println!("--- Test 1: Non-staircase [[1,s],[t,1]] = 1 + st ---");
    let mut found_zero = false;
    for &s in &test_points {
        for &t in &test_points {
            let mat = vec![
                vec![C64::new(1.0, 0.0), s],
                vec![t, C64::new(1.0, 0.0)],
            ];
            let p = perm_complex(&mat);
            if p.norm2() < 1e-20 {
                println!("  ZERO at s={:.2}+{:.2}i, t={:.2}+{:.2}i: F={:.6}+{:.6}i",
                    s.re, s.im, t.re, t.im, p.re, p.im);
                found_zero = true;
            }
        }
    }
    if !found_zero { println!("  No zeros found (but we know s=t=i gives zero)"); }
    // Exact check at s=t=i
    let s = C64::new(0.0, 1.0); let t = C64::new(0.0, 1.0);
    let mat = vec![vec![C64::new(1.0,0.0), s], vec![t, C64::new(1.0,0.0)]];
    let p = perm_complex(&mat);
    println!("  At s=t=i: F = {:.6} + {:.6}i (norm²={:.2e})", p.re, p.im, p.norm2());
    
    // Test 2: Staircase [[s,t],[1,1]] = s + t (should be stable)
    println!("\n--- Test 2: Staircase [[s,t],[1,1]] = s + t ---");
    let mut min_norm = f64::MAX;
    for &s in &test_points { for &t in &test_points {
        let mat = vec![vec![s, t], vec![C64::new(1.0,0.0), C64::new(1.0,0.0)]];
        let p = perm_complex(&mat);
        if p.norm2() < min_norm { min_norm = p.norm2(); }
    }}
    println!("  Min |F|² = {:.2e} (stable if >> 0)", min_norm);
    
    // Test 3: Staircase matrices from actual Ferrers boards
    println!("\n--- Test 3: Staircase hit matrices M(s) from Ferrers boards ---");
    // For each board and sub-partition, build M(s) and test
    let mut total = 0u64;
    let mut min_norm_all = f64::MAX;
    let mut worst_case = String::new();
    
    for n in 2..=5usize {
        let boards = gen_boards(n);
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) { continue; }
            let m = *board.iter().max().unwrap() as usize;
            let subs = sub_partitions(board);
            for mu in &subs {
                for j in 0..n {
                    if mu[j] >= board[j] { continue; }
                    let mut mu_prime = mu.clone();
                    mu_prime[j] += 1;
                    if j > 0 && mu_prime[j] > mu_prime[j-1] { continue; }
                    let col_j = mu[j] as usize + 1;
                    // Build M(s): matrix for mu_prime with entry (j, col_j) = s
                    // Test at several (s,t) in upper half-plane
                    for &s_val in &test_points {
                        for &t_val in &test_points {
                            let mut mat = Vec::new();
                            for i in 0..n {
                                let mut row = Vec::new();
                                for k in 1..=m {
                                    if k > board[i] as usize {
                                        row.push(C64::new(0.0, 0.0));
                                    } else if i == j && k == col_j {
                                        row.push(s_val); // the s entry
                                    } else if k <= mu_prime[i] as usize {
                                        row.push(C64::new(1.0, 0.0));
                                    } else {
                                        row.push(t_val);
                                    }
                                }
                                mat.push(row);
                            }
                            let p = perm_complex(&mat);
                            total += 1;
                            if p.norm2() < min_norm_all {
                                min_norm_all = p.norm2();
                                worst_case = format!("board={:?} mu'={:?} j={} s={:.1}+{:.1}i t={:.1}+{:.1}i F={:.4}+{:.4}i",
                                    board, mu_prime, j, s_val.re, s_val.im, t_val.re, t_val.im, p.re, p.im);
                            }
                        }
                    }
                }
            }
        }
    }
    println!("  {} evaluations tested", total);
    println!("  Min |F|² = {:.2e}", min_norm_all);
    println!("  Worst case: {}", worst_case);
    if min_norm_all > 1e-10 {
        println!("  CONCLUSION: No zeros found in upper half-plane → consistent with stability <<<");
    } else {
        println!("  WARNING: Found near-zero → stability may fail!");
    }
    
    // Test 4: Random staircase matrices (not from Ferrers boards)
    println!("\n--- Test 4: Random staircase matrices ---");
    let mut rand_tests = 0u64;
    let mut rand_min_norm = f64::MAX;
    // Generate random staircase: boundaries b_1 >= b_2 >= ... >= b_n, board widths w_i >= b_i
    for n in 2..=5 {
        for _ in 0..100 {
            let m = n + 2; // width
            let mut boundaries: Vec<usize> = (0..n).map(|i| (m - i) % (m + 1)).collect();
            boundaries.sort_unstable_by(|a, b| b.cmp(a)); // decreasing
            let j = 0; // always replace entry in first row
            let col_s = boundaries[j]; // boundary position
            if col_s == 0 || col_s > m { continue; }
            
            for &s_val in &test_points {
                for &t_val in &test_points {
                    let mut mat = Vec::new();
                    for i in 0..n {
                        let mut row = Vec::new();
                        for k in 0..m {
                            if i == j && k + 1 == col_s {
                                row.push(s_val);
                            } else if k < boundaries[i] {
                                row.push(C64::new(1.0, 0.0));
                            } else {
                                row.push(t_val);
                            }
                        }
                        mat.push(row);
                    }
                    let p = perm_complex(&mat);
                    rand_tests += 1;
                    if p.norm2() < rand_min_norm { rand_min_norm = p.norm2(); }
                }
            }
        }
    }
    println!("  {} evaluations", rand_tests);
    println!("  Min |F|² = {:.2e}", rand_min_norm);
    if rand_min_norm > 1e-10 {
        println!("  CONCLUSION: Consistent with stability <<<");
    }
}

fn board_to_perm(b: &[u8]) -> Vec<u8> { let n = b.len(); let mut p = vec![0u8; n]; let mut u = vec![false; n+1]; for i in 0..n { for c in (1..=(b[i] as usize).min(n)).rev() { if !u[c] { p[i] = c as u8; u[c] = true; break; } } } p }
fn is_312_avoiding(perm: &[u8]) -> bool { let n = perm.len(); for i in 0..n { for j in i+1..n { for k in j+1..n { if perm[k] < perm[i] && perm[i] < perm[j] { return false; } } } } true }
fn gen_boards(n: usize) -> Vec<Vec<u8>> { let mut r = vec![]; let mut c = vec![]; gb(n, n, 0, &mut c, &mut r); r }
fn gb(n: usize, mx: usize, d: usize, c: &mut Vec<u8>, r: &mut Vec<Vec<u8>>) { if d == n { r.push(c.clone()); return; } for v in (d+1).max(if d > 0 { c[d-1] as usize } else { 1 })..=mx { c.push(v as u8); gb(n, mx, d+1, c, r); c.pop(); } }
fn sub_partitions(lambda: &[u8]) -> Vec<Vec<u8>> { let n = lambda.len(); let mut result = Vec::new(); let mut mu = vec![0u8; n];
    fn gen(lam: &[u8], mu: &mut Vec<u8>, pos: usize, mx: u8, res: &mut Vec<Vec<u8>>) { if pos == lam.len() { res.push(mu.clone()); return; } let u = lam[pos].min(mx); for v in 0..=u { mu[pos] = v; gen(lam, mu, pos+1, v, res); } }
    gen(lambda, &mut mu, 0, lambda[0], &mut result); result }
