//! Check: what is the ACTUAL recursion for shifted excedance on S_n?
//! F^(s)_j(S_n) = weight(j,s) * perm(submatrix).
//! The submatrix has rows 2..n (relabeled 1..n-1) with columns {1..n}\{j}.
//! Entry at new row i', column j': t if j' > i'+1+s, else 1.
//! This is NOT F^(s+1)(S_{n-1}) because columns are {1..n}\{j}, not {1..n-1}.
//! But maybe we can RELABEL columns? Removing column j and relabeling
//! j+1..n as j..n-1 changes the staircase in a j-dependent way.
//!
//! Key idea: maybe refine by first AND last entry, or by the "column gap".
use polynomial_tools::real_rootedness::format_poly;

fn pt(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
    v
}
fn pz(p: &[i64]) -> bool {
    p.iter().all(|&c| c == 0)
}
fn pa(a: &[i64], b: &[i64]) -> Vec<i64> {
    let l = a.len().max(b.len());
    let mut r = vec![0i64; l];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] += v;
    }
    pt(&r)
}
fn pmul(a: &[i64], b: &[i64]) -> Vec<i64> {
    if a.is_empty() || b.is_empty() {
        return vec![0];
    }
    let mut r = vec![0i64; a.len() + b.len() - 1];
    for (i, &av) in a.iter().enumerate() {
        for (j, &bv) in b.iter().enumerate() {
            r[i + j] += av * bv;
        }
    }
    pt(&r)
}

/// Compute shifted excedance polynomial on S_n by brute force.
/// shift s: exc at position i means sigma_i > i + s.
fn shifted_exc_bruteforce(perm_prefix: &[u8], remaining: &[u8], s: i32, n: usize) -> Vec<i64> {
    if remaining.is_empty() {
        let exc = (0..n)
            .filter(|&i| perm_prefix[i] as i32 > (i as i32 + 1) + s)
            .count();
        let mut p = vec![0i64; exc + 1];
        p[exc] = 1;
        return p;
    }
    let mut result = vec![0i64];
    for (idx, &val) in remaining.iter().enumerate() {
        let mut next_perm = perm_prefix.to_vec();
        next_perm.push(val);
        let mut next_rem: Vec<u8> = remaining.to_vec();
        next_rem.remove(idx);
        result = pa(
            &result,
            &shifted_exc_bruteforce(&next_perm, &next_rem, s, n),
        );
    }
    result
}

/// Compute F^(s)_{first=j, last=l}(S_n) refined by first AND last entry.
fn shifted_exc_refined(n: usize, s: i32) -> Vec<Vec<Vec<i64>>> {
    // result[j][l] = polynomial for first entry j, last entry l
    let mut result = vec![vec![vec![0i64]; n + 1]; n + 1];
    let vals: Vec<u8> = (1..=n as u8).collect();
    // Generate all permutations
    fn gen(prefix: &[u8], remaining: &[u8], n: usize, s: i32, result: &mut Vec<Vec<Vec<i64>>>) {
        if remaining.is_empty() {
            let j = prefix[0] as usize;
            let l = *prefix.last().unwrap() as usize;
            let exc = (0..n)
                .filter(|&i| prefix[i] as i32 > (i as i32 + 1) + s)
                .count();
            while result[j][l].len() <= exc {
                result[j][l].push(0);
            }
            result[j][l][exc] += 1;
            return;
        }
        for (idx, &val) in remaining.iter().enumerate() {
            let mut np = prefix.to_vec();
            np.push(val);
            let mut nr: Vec<u8> = remaining.to_vec();
            nr.remove(idx);
            gen(&np, &nr, n, s, result);
        }
    }
    gen(&[], &vals, n, s, &mut result);
    result
}

fn main() {
    // Print F^(s)_{j,l}(S_n) for small n and look for recurrence patterns
    for n in 3..=5 {
        println!("=== S_{} ===", n);
        for s in 0..n as i32 {
            println!("  s={}:", s);
            let refined = shifted_exc_refined(n, s);
            // Print compact table
            for j in 1..=n {
                for l in 1..=n {
                    if j == l {
                        continue;
                    }
                    let p = pt(&refined[j][l]);
                    if !pz(&p) {
                        print!("    F({},{}) = {}  ", j, l, format_poly(&p));
                    }
                }
            }
            println!();
            // Check: Σ_l F(j,l) and Σ_j F(j,l)
            let mut by_first: Vec<Vec<i64>> = vec![vec![0i64]; n + 1];
            let mut by_last: Vec<Vec<i64>> = vec![vec![0i64]; n + 1];
            for j in 1..=n {
                for l in 1..=n {
                    if j != l {
                        by_first[j] = pa(&by_first[j], &refined[j][l]);
                        by_last[l] = pa(&by_last[l], &refined[j][l]);
                    }
                }
            }
            print!("    Σ_l: ");
            for j in 1..=n {
                if !pz(&by_first[j]) {
                    print!("F({},*) = {}  ", j, format_poly(&pt(&by_first[j])));
                }
            }
            println!();
            print!("    Σ_j: ");
            for l in 1..=n {
                if !pz(&by_last[l]) {
                    print!("F(*,{}) = {}  ", l, format_poly(&pt(&by_last[l])));
                }
            }
            println!();
        }
        println!();
    }
}
