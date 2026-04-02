/// Cross-check: sink_polynomial_tree() vs acyclic_sink_polynomial() on trees.
/// Also benchmarks the tree method on large caterpillars.
use combinatoric_core::graph::Graph;
use polynomial_tools::real_rootedness::is_real_rooted;
use std::time::Instant;

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

fn main() {
    println!("=== Cross-check: tree method vs 2^m brute-force ===\n");

    let mut all_ok = true;

    // Test all trees up to n=8 via Prüfer sequences
    for n in 2..=8 {
        let mut seen = std::collections::BTreeSet::new();
        let mut count = 0;
        let mut mismatches = 0;

        if n == 2 {
            let g = Graph::path(2);
            let tree = g.line_graph().acyclic_sink_polynomial();
            let fast_tree = g.sink_polynomial_tree();
            count = 1;
            if tree != fast_tree {
                mismatches += 1;
                println!("  MISMATCH P_2: fast={:?} tree={:?}", tree, fast_tree);
            }
        } else {
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
                let mut ds: Vec<usize> = (0..n).map(|v| g.degree(v)).collect();
                ds.sort();
                if !seen.insert(ds) {
                    continue;
                }
                count += 1;

                let fast = g.line_graph().acyclic_sink_polynomial();
                let tree = g.sink_polynomial_tree();
                if fast != tree {
                    mismatches += 1;
                    all_ok = false;
                    let degs: Vec<usize> = (0..n).map(|v| g.degree(v)).collect();
                    println!(
                        "  MISMATCH n={} degs={:?}: fast={:?} tree={:?}",
                        n, degs, fast, tree
                    );
                }
            }
        }

        let status = if mismatches == 0 { "OK" } else { "FAIL" };
        println!("n={}: {} trees tested [{}]", n, count, status);
    }

    if all_ok {
        println!("\nAll cross-checks passed!\n");
    } else {
        println!("\nSOME MISMATCHES FOUND!\n");
        return;
    }

    // Benchmark: large caterpillars using tree method
    println!("=== Tree method on large caterpillars ===\n");
    let big_cats: Vec<(&str, Vec<usize>)> = vec![
        ("Cat(5,5,5,5,5)", vec![5, 5, 5, 5, 5]),
        ("Cat(10,10,10)", vec![10, 10, 10]),
        ("Cat(0,20,0)", vec![0, 20, 0]),
        ("Cat(3,3,3,3,3,3,3)", vec![3, 3, 3, 3, 3, 3, 3]),
        ("Cat(1,1,1,1,1,1,1,1,1,1,1,1)", vec![1; 12]),
        ("Cat(0,0,...,0) = P_20", vec![0; 20]),
        ("Cat(0,0,...,0) = P_50", vec![0; 50]),
        ("Star K_{1,20}", vec![20]),
    ];

    for (name, a) in &big_cats {
        // Handle single-vertex spine (star)
        let g = if a.len() == 1 {
            Graph::complete_bipartite(1, a[0])
        } else {
            build_caterpillar(a)
        };

        let t0 = Instant::now();
        let poly = g.sink_polynomial_tree();
        let ms = t0.elapsed().as_micros();
        let rr = is_real_rooted(&poly);

        println!(
            "{:<35} |V|={:<4} poly={:?}  [{}, {}μs]",
            name,
            g.num_vertices(),
            poly,
            if rr { "RR" } else { "NOT RR" },
            ms
        );
    }
}
