/// Experiment: Compute the sink polynomial of L(T) for all trees T
/// and check real-rootedness.
///
/// Conjecture: S_{L(G)}(t) is always real-rooted.
///
/// For a tree T, the number of acyclic orientations of L(T)
/// equals the product of the degrees of the vertices of T.
/// We filter to keep this product below a threshold.
use combinatoric_core::graph::Graph;
use num_bigint::BigInt;
use num_traits::{Signed, ToPrimitive, Zero};
use polynomial_tools::real_rootedness::is_real_rooted;

/// Check real-rootedness of a BigInt polynomial.
/// Divides out the GCD of coefficients to reduce to i64 range.
fn is_real_rooted_bigint(coeffs: &[BigInt]) -> bool {
    if coeffs.is_empty() {
        return true;
    }
    // Compute GCD of all nonzero coefficients
    let mut g = BigInt::from(0);
    for c in coeffs {
        if !c.is_zero() {
            g = num::integer::gcd(g.abs(), c.abs());
        }
    }
    if g.is_zero() {
        return true;
    }
    // Divide out GCD and convert to i64
    let reduced: Vec<i64> = coeffs
        .iter()
        .map(|c| {
            let r = c / &g;
            r.to_i64()
                .expect("coefficient too large even after GCD reduction")
        })
        .collect();
    is_real_rooted(&reduced)
}

fn check_and_print_tree(name: &str, g: &Graph) -> bool {
    let coeffs = g.sink_polynomial_tree_bigint();
    let rr = is_real_rooted_bigint(&coeffs);
    let status = if rr { "RR" } else { "NOT RR <<<" };
    let degs: Vec<usize> = (0..g.num_vertices()).map(|v| g.degree(v)).collect();
    println!(
        "{:<30} degs={:<20} S_{{L(T)}}(t) = {:?}  [{}]",
        name,
        format!("{:?}", degs),
        coeffs,
        status
    );
    rr
}

fn check_and_print_fast(name: &str, g: &Graph) -> bool {
    let lg = g.line_graph();
    let coeffs = lg.acyclic_sink_polynomial();
    let rr = is_real_rooted(&coeffs);
    let status = if rr { "RR" } else { "NOT RR <<<" };
    println!(
        "{:<30} |V(L)|={:<3} |E(L)|={:<3} S(t) = {:?}  [{}]",
        name,
        lg.num_vertices(),
        lg.num_edges(),
        coeffs,
        status
    );
    rr
}

/// Build a caterpillar Cat(a_1, ..., a_k):
/// Path v_0 - v_1 - ... - v_{k-1}, with a_i leaves attached to v_i.
fn build_caterpillar(a: &[usize]) -> Graph {
    let k = a.len();
    let total: usize = k + a.iter().sum::<usize>();
    let mut edges = Vec::new();
    // Path edges
    for i in 0..k - 1 {
        edges.push((i, i + 1));
    }
    // Leaf edges
    let mut next_vertex = k;
    for i in 0..k {
        for _ in 0..a[i] {
            edges.push((i, next_vertex));
            next_vertex += 1;
        }
    }
    assert_eq!(next_vertex, total);
    Graph::new(total, &edges)
}

fn caterpillar_name(a: &[usize]) -> String {
    let parts: Vec<String> = a.iter().map(|x| x.to_string()).collect();
    format!("Cat({})", parts.join(","))
}

fn main() {
    let mut all_rr = true;
    let mut warnings = Vec::new();

    // --- Paths (tree method, no size limit) ---
    println!("=== Paths (tree method) ===");
    for n in 2..=30 {
        let g = Graph::path(n);
        let rr = check_and_print_tree(&format!("P_{}", n), &g);
        if !rr {
            warnings.push(format!("P_{}", n));
        }
        all_rr &= rr;
    }
    println!();

    // --- Cycles (2^m method) ---
    println!("=== Cycles (brute-force) ===");
    for n in 3..=12 {
        let g = Graph::cycle(n);
        let rr = check_and_print_fast(&format!("C_{}", n), &g);
        if !rr {
            warnings.push(format!("C_{}", n));
        }
        all_rr &= rr;
    }
    println!();

    // --- Caterpillar graphs (tree method, no size limit on L(T)) ---
    let max_leaf = 6;
    for k in 2..=6 {
        println!(
            "=== Caterpillars with k={} spine vertices, a_i <= {} ===",
            k, max_leaf
        );
        let mut count = 0u64;

        let total_combos = (max_leaf + 1_usize).pow(k as u32);
        for code in 0..total_combos {
            let mut a = vec![0usize; k];
            let mut c = code;
            for i in (0..k).rev() {
                a[i] = c % (max_leaf + 1);
                c /= max_leaf + 1;
            }

            let g = build_caterpillar(&a);
            let name = caterpillar_name(&a);
            count += 1;
            let rr = check_and_print_tree(&name, &g);
            if !rr {
                warnings.push(name.clone());
                println!("  WARNING: {} is NOT real-rooted!", name);
            }
            all_rr &= rr;
        }

        println!("  ({} caterpillars tested)", count);
        println!();
    }

    // --- Summary ---
    if all_rr {
        println!("ALL tested graphs have real-rooted S_{{L(G)}}(t)!");
    } else {
        println!("COUNTEREXAMPLES FOUND:");
        for w in &warnings {
            println!("  - {}", w);
        }
    }
}
