use std::collections::BTreeSet;
fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> {
    let n = perm.len(); let mut vis: BTreeSet<Vec<u8>> = BTreeSet::new(); let mut q: BTreeSet<Vec<u8>> = BTreeSet::new();
    q.insert(perm.to_vec()); while let Some(cur) = q.pop_last() { for i in 0..n { for j in i+1..n {
        if cur[i] > cur[j] { let mut c = cur.clone(); c.swap(i, j); if !vis.contains(&c) { q.insert(c); } }
    }} vis.insert(cur); } vis.into_iter().collect()
}
fn main() {
    // Board (2,2,3), perm = [2,1,3] (312-avoiding ✓)
    let perm = vec![2u8, 2, 3]; // board, not perm!
    // board_to_perm for (2,2,3): row 0: max avail ≤ 2 = 2, row 1: max avail ≤ 2 = 1, row 2: max avail ≤ 3 = 3
    let p = vec![2u8, 1, 3];
    let ideal = bruhat_lower_ideal(&p);
    println!("Board (2,2,3), perm [2,1,3], ideal size: {}", ideal.len());
    
    let mu_prime = vec![2u8, 2, 2]; // μ' = (2,2,2) after cover
    let row_j = 2usize; // 0-indexed
    
    for k in 1..=3 {
        let mut poly = vec![0i64; 4];
        for sigma in &ideal {
            if sigma[row_j] as usize != k { continue; }
            let mut hits = 0;
            for i in 0..3 {
                if i == row_j { continue; }
                if sigma[i] as usize > mu_prime[i] as usize { hits += 1; }
            }
            poly[hits] += 1;
        }
        // Trim trailing zeros
        while poly.len() > 1 && *poly.last().unwrap() == 0 { poly.pop(); }
        println!("  C_{{2,{}}} (sigma[2]={}) = {:?}  deg={}", k, k, poly, poly.len()-1);
    }
}
