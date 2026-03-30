//! Explore restricted hit polynomials H_μ^λ(t) for sub-partitions μ ⊆ λ.
//! M_{i,j} = t if μ_i < j ≤ λ_i, 1 if j ≤ μ_i, 0 if j > λ_i.
use polynomial_tools::real_rootedness::{is_real_rooted, format_poly};

fn pt(p: &[i64]) -> Vec<i64> { let mut v = p.to_vec(); while v.len() > 1 && *v.last().unwrap() == 0 { v.pop(); } v }
fn pz(p: &[i64]) -> bool { p.iter().all(|&c| c == 0) }
fn pa(a: &[i64], b: &[i64]) -> Vec<i64> { let l = a.len().max(b.len()); let mut r = vec![0i64; l]; for (i, &v) in a.iter().enumerate() { r[i] += v; } for (i, &v) in b.iter().enumerate() { r[i] += v; } pt(&r) }
fn pmul(a: &[i64], b: &[i64]) -> Vec<i64> { if a.is_empty() || b.is_empty() { return vec![0]; } let mut r = vec![0i64; a.len() + b.len() - 1]; for (i, &av) in a.iter().enumerate() { for (j, &bv) in b.iter().enumerate() { r[i+j] += av * bv; } } pt(&r) }

fn permanent(mat: &[Vec<Vec<i64>>]) -> Vec<i64> {
    let n = mat.len();
    if n == 0 { return vec![1]; }
    let m = mat[0].len();
    let mut result = vec![0i64];
    for j in 0..m {
        if pz(&mat[0][j]) { continue; }
        let mut sub = Vec::new();
        for i in 1..n {
            let mut row = Vec::new();
            for jj in 0..m { if jj == j { continue; } row.push(mat[i][jj].clone()); }
            sub.push(row);
        }
        result = pa(&result, &pmul(&mat[0][j], &permanent(&sub)));
    }
    pt(&result)
}

fn hit_poly(lambda: &[u8], mu: &[u8]) -> Vec<i64> {
    let n = lambda.len();
    let m = *lambda.iter().max().unwrap_or(&0) as usize;
    let mut mat = Vec::new();
    for i in 0..n {
        let mu_i = if i < mu.len() { mu[i] as usize } else { 0 };
        let lam_i = lambda[i] as usize;
        let mut row = Vec::new();
        for j in 0..m {
            let j1 = j + 1; // 1-indexed
            if j1 > lam_i { row.push(vec![0i64]); }
            else if j1 <= mu_i { row.push(vec![1i64]); }
            else { row.push(vec![0, 1]); } // t
        }
        mat.push(row);
    }
    permanent(&mat)
}

/// Generate all sub-partitions mu with mu_i <= lambda_i and mu weakly decreasing.
fn sub_partitions(lambda: &[u8]) -> Vec<Vec<u8>> {
    let n = lambda.len();
    let mut result = Vec::new();
    let mut mu = vec![0u8; n];
    fn gen(lambda: &[u8], mu: &mut Vec<u8>, pos: usize, max_val: u8, result: &mut Vec<Vec<u8>>) {
        if pos == lambda.len() { result.push(mu.clone()); return; }
        let upper = lambda[pos].min(max_val);
        for v in 0..=upper { mu[pos] = v; gen(lambda, mu, pos + 1, v, result); }
    }
    gen(lambda, &mut mu, 0, lambda[0], &mut result);
    result
}

fn main() {
    let max_n: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(5);
    
    println!("=== Restricted hit polynomials H_μ^λ(t) ===\n");
    
    let mut total = 0u64;
    let mut rr_count = 0u64;
    let mut non_rr = Vec::new();
    
    // Test all 312-avoiding Ferrers boards
    for n in 2..=max_n {
        // Full boards [n,...,n] with all sub-partitions
        let lambda: Vec<u8> = vec![n as u8; n];
        let subs = sub_partitions(&lambda);
        for mu in &subs {
            let h = hit_poly(&lambda, mu);
            total += 1;
            if is_real_rooted(&pt(&h)) || pt(&h).len() <= 2 {
                rr_count += 1;
            } else {
                non_rr.push((lambda.clone(), mu.clone(), pt(&h)));
            }
        }
    }
    
    // Also test non-square Ferrers boards
    let test_boards: Vec<Vec<u8>> = vec![
        vec![3,3,2], vec![4,3,2], vec![4,4,3], vec![4,4,3,2],
        vec![5,4,3], vec![5,5,4,3], vec![5,5,3,2],
    ];
    for lambda in &test_boards {
        let subs = sub_partitions(lambda);
        for mu in &subs {
            let h = hit_poly(lambda, mu);
            total += 1;
            if is_real_rooted(&pt(&h)) || pt(&h).len() <= 2 {
                rr_count += 1;
            } else {
                non_rr.push((lambda.clone(), mu.clone(), pt(&h)));
            }
        }
    }
    
    println!("Total (lambda, mu) pairs tested: {}", total);
    println!("Real-rooted: {}/{}", rr_count, total);
    if !non_rr.is_empty() {
        println!("\nNon-real-rooted examples (first 10):");
        for (lam, mu, h) in non_rr.iter().take(10) {
            println!("  λ={:?} μ={:?}: H = {}", lam, mu, format_poly(&h));
        }
    } else {
        println!("\nALL REAL-ROOTED!");
    }
    
    // Print some examples
    println!("\n=== Selected examples ===");
    for &(ref lam, ref mu_desc) in &[
        (vec![4u8,4,4,4], "excedance (μ=1,2,3,4)"),
        (vec![4u8,4,4,4], "shifted exc s=1 (μ=2,3,4,4)"),
        (vec![4u8,4,4,4], "staircase (μ=0,1,2,3)"),
        (vec![4u8,4,4,4], "half (μ=2,2,2,2)"),
    ] {
        let mu: Vec<u8> = match mu_desc {
            s if s.contains("excedance") => (1..=lam.len() as u8).collect(),
            s if s.contains("s=1") => (2..=lam.len() as u8 + 1).map(|x| x.min(lam[0])).collect(),
            s if s.contains("staircase") => (0..lam.len() as u8).collect(),
            s if s.contains("half") => vec![2; lam.len()],
            _ => vec![0; lam.len()],
        };
        let h = hit_poly(&lam, &mu);
        println!("  λ={:?} μ={:?} ({}): H = {}", lam, mu, mu_desc, format_poly(&pt(&h)));
    }
}
