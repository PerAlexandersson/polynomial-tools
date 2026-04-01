/// Explore palindromicity and gamma-positivity of S_{L(T)}(t) for trees T.
///
/// For each tree, compute S_{L(T)}(t), strip the leading factor of t^k,
/// and check whether the resulting polynomial is palindromic.
/// When palindromic, compute the gamma-vector and check non-negativity.
use combinatoric_core::graph::Graph;
use num_bigint::BigInt;
use num_traits::{ToPrimitive, Zero};
use polynomial_tools::real_rootedness::{
    format_poly, gamma_coefficients, is_gamma_positive, is_palindromic,
};
use std::collections::BTreeSet;

/// Generate tree from Prüfer sequence.
fn prufer_to_tree(seq: &[usize], n: usize) -> Graph {
    let mut degree = vec![1usize; n];
    for &v in seq {
        degree[v] += 1;
    }
    let mut edges = Vec::new();
    let mut seq_iter = seq.iter();
    for _ in 0..n - 2 {
        let &v = seq_iter.next().unwrap();
        for u in 0..n {
            if degree[u] == 1 {
                edges.push((u, v));
                degree[u] -= 1;
                degree[v] -= 1;
                break;
            }
        }
    }
    let remaining: Vec<usize> = (0..n).filter(|&v| degree[v] == 1).collect();
    edges.push((remaining[0], remaining[1]));
    Graph::new(n, &edges)
}

fn degree_seq(g: &Graph) -> Vec<usize> {
    let mut ds: Vec<usize> = (0..g.num_vertices()).map(|v| g.degree(v)).collect();
    ds.sort();
    ds
}

/// Strip leading zeros and convert BigInt -> i64 (dividing out GCD).
/// Returns None if coefficients overflow i64 even after GCD reduction.
fn to_i64_stripped(coeffs: &[BigInt]) -> Option<(usize, Vec<i64>)> {
    let first_nz = coeffs
        .iter()
        .position(|c| !c.is_zero())
        .unwrap_or(coeffs.len());
    let tail = &coeffs[first_nz..];
    if tail.is_empty() {
        return Some((first_nz, vec![]));
    }
    // GCD reduction
    let mut g = BigInt::from(0);
    for c in tail {
        if !c.is_zero() {
            g = num::integer::gcd(
                if g < BigInt::from(0) { -&g } else { g.clone() },
                if c < &BigInt::from(0) { -c } else { c.clone() },
            );
        }
    }
    if g.is_zero() {
        return Some((first_nz, vec![0; tail.len()]));
    }
    let reduced: Option<Vec<i64>> = tail.iter().map(|c| (c / &g).to_i64()).collect();
    reduced.map(|r| (first_nz, r))
}

/// Build a spider Sp(k, m): k legs of length m from a central vertex.
fn build_spider(k: usize, m: usize) -> Graph {
    let n = 1 + k * m;
    let mut edges = Vec::new();
    for leg in 0..k {
        let base = 1 + leg * m;
        edges.push((0, base));
        for j in 1..m {
            edges.push((base + j - 1, base + j));
        }
    }
    Graph::new(n, &edges)
}

/// Build a caterpillar Cat(a_1, ..., a_k).
fn build_caterpillar(a: &[usize]) -> Graph {
    let k = a.len();
    let total: usize = k + a.iter().sum::<usize>();
    let mut edges = Vec::new();
    for i in 0..k - 1 {
        edges.push((i, i + 1));
    }
    let mut next = k;
    for i in 0..k {
        for _ in 0..a[i] {
            edges.push((i, next));
            next += 1;
        }
    }
    Graph::new(total, &edges)
}

fn check_tree(
    name: &str,
    g: &Graph,
    palindromic_list: &mut Vec<String>,
    not_gamma: &mut Vec<String>,
) {
    let coeffs = g.sink_polynomial_tree_bigint();
    let (min_pow, stripped) = match to_i64_stripped(&coeffs) {
        Some(v) => v,
        None => {
            eprintln!("  SKIP {} (coeffs overflow i64)", name);
            return;
        }
    };
    if stripped.len() <= 1 {
        return; // trivial
    }
    if is_palindromic(&stripped) {
        let gp = is_gamma_positive(&stripped);
        let gamma = gamma_coefficients(&stripped).unwrap();
        let tag = if gp { "gamma+" } else { "NOT gamma+" };
        println!(
            "PALINDROMIC  {:<40} t^{} * ({})  gamma={:?}  [{}]",
            name,
            min_pow,
            format_poly(&stripped),
            gamma,
            tag,
        );
        palindromic_list.push(name.to_string());
        if !gp {
            not_gamma.push(name.to_string());
        }
    }
}

fn main() {
    println!("=== Palindromicity & gamma-positivity of S_{{L(T)}}(t) for trees ===\n");

    let mut palindromic_list = Vec::new();
    let mut not_gamma = Vec::new();
    let mut total_trees = 0usize;

    // All trees up to n vertices (Prüfer enum is O(n^{n-2}), cap at 9)
    let max_n = 9;
    for n in 2..=max_n {
        let mut seen = BTreeSet::new();
        if n == 2 {
            total_trees += 1;
            check_tree(
                "P_2",
                &Graph::path(2),
                &mut palindromic_list,
                &mut not_gamma,
            );
            continue;
        }
        let seq_len = n - 2;
        let total = (n as u64).pow(seq_len as u32);
        for code in 0..total {
            let mut seq = vec![0usize; seq_len];
            let mut c = code;
            for i in (0..seq_len).rev() {
                seq[i] = (c % n as u64) as usize;
                c /= n as u64;
            }
            let g = prufer_to_tree(&seq, n);
            let ds = degree_seq(&g);
            if !seen.insert(ds.clone()) {
                continue;
            }
            total_trees += 1;
            let name = format!("n={} ds={:?}", n, ds);
            check_tree(&name, &g, &mut palindromic_list, &mut not_gamma);
        }
        eprintln!("n={} done ({} trees so far)", n, total_trees);
    }

    // Named families: Spiders Sp(k, m)
    println!("\n--- Spiders Sp(k,m) ---");
    for k in 2..=10 {
        for m in 1..=8 {
            let g = build_spider(k, m);
            let name = format!("Sp({},{})", k, m);
            check_tree(&name, &g, &mut palindromic_list, &mut not_gamma);
        }
    }

    // Double stars S(a,b)
    println!("\n--- Double stars S(a,b) ---");
    for a in 1..=8 {
        for b in a..=8 {
            let g = Graph::complete_bipartite(1, a + b + 0); // hack: build manually
                                                             // Build S(a,b): edge {0,1}, a leaves at 0, b leaves at 1
            let n = 2 + a + b;
            let mut edges = vec![(0usize, 1usize)];
            for i in 0..a {
                edges.push((0, 2 + i));
            }
            for i in 0..b {
                edges.push((1, 2 + a + i));
            }
            let g = Graph::new(n, &edges);
            let name = format!("S({},{})", a, b);
            check_tree(&name, &g, &mut palindromic_list, &mut not_gamma);
        }
    }

    // Symmetric caterpillars Cat(a,...,a)
    println!("\n--- Uniform caterpillars Cat(a x k) ---");
    for k in 2..=12 {
        for a in 0..=5 {
            let avec = vec![a; k];
            let g = build_caterpillar(&avec);
            let name = format!("Cat({}x{})", a, k);
            check_tree(&name, &g, &mut palindromic_list, &mut not_gamma);
        }
    }

    // Palindromic caterpillars Cat(a_1,...,a_k) with a_i = a_{k+1-i}
    println!("\n--- Palindromic-spine caterpillars ---");
    for k in 3..=8 {
        let half = (k + 1) / 2;
        for code in 0..(5usize.pow(half as u32)) {
            let mut avec = vec![0usize; k];
            let mut c = code;
            for i in 0..half {
                avec[i] = c % 5;
                avec[k - 1 - i] = avec[i];
                c /= 5;
            }
            let g = build_caterpillar(&avec);
            let name = format!("Cat({:?})", avec);
            check_tree(&name, &g, &mut palindromic_list, &mut not_gamma);
        }
    }

    // Summary
    println!("\n=== SUMMARY ===");
    println!("Total trees examined: {}", total_trees);
    println!("Palindromic (stripped): {}", palindromic_list.len());
    if not_gamma.is_empty() {
        println!("All palindromic cases are gamma-positive!");
    } else {
        println!("NOT gamma-positive:");
        for name in &not_gamma {
            println!("  {}", name);
        }
    }
}
