/// Experiment: Sink polynomial over ALL orientations (not just acyclic)
/// for line graphs L(G).
///
/// S^{all}_{L(G)}(t) = sum over all 2^m orientations of L(G) of t^{#sinks}
///
/// Question: Is this polynomial always real-rooted when G is a line graph?
use combinatoric_core::graph::Graph;
use polynomial_tools::real_rootedness::is_real_rooted;

fn check(name: &str, g: &Graph) {
    let lg = g.line_graph();
    let m = lg.num_edges();
    if m > 20 {
        println!("{:<25} |E(L)|={:<3} -- skipped (too many edges)", name, m);
        return;
    }
    let coeffs = lg.sink_polynomial_all_orientations();
    let rr = is_real_rooted(&coeffs);
    let status = if rr { "RR" } else { "NOT RR <<<" };

    // Also compute the acyclic version for comparison
    let acyclic = lg.sink_polynomial_fast();

    println!(
        "{:<25} |V(L)|={:<3} |E(L)|={:<3}  all={:?}  [{}]",
        name,
        lg.num_vertices(),
        m,
        coeffs,
        status,
    );
    println!("{:<25} {:>10}          acyc={:?}", "", "", acyclic,);
}

fn main() {
    println!("=== Sink poly over ALL orientations of L(G) ===\n");

    // Paths: L(P_n) = P_{n-1}, a tree, so all orientations = acyclic
    println!("--- Paths (should match acyclic, since L(P_n) is a tree) ---");
    for n in 3..=6 {
        check(&format!("P_{}", n), &Graph::path(n));
    }

    // Cycles: L(C_n) = C_n
    println!("\n--- Cycles ---");
    for n in 3..=10 {
        check(&format!("C_{}", n), &Graph::cycle(n));
    }

    // Complete graphs: L(K_n)
    println!("\n--- Complete graphs ---");
    for n in 3..=5 {
        check(&format!("K_{}", n), &Graph::complete(n));
    }

    // Complete bipartite
    println!("\n--- Complete bipartite ---");
    for a in 2..=4 {
        for b in a..=4 {
            if a * b <= 12 {
                check(
                    &format!("K_{{{},{}}}", a, b),
                    &Graph::complete_bipartite(a, b),
                );
            }
        }
    }

    // Stars: L(K_{1,n}) = K_n (complete, so acyclic = all for tournaments? No, K_n can have cyclic orientations)
    println!("\n--- Stars K_{{1,n}} ---");
    for n in 2..=6 {
        check(&format!("K_{{1,{}}}", n), &Graph::complete_bipartite(1, n));
    }

    // Petersen
    println!("\n--- Petersen ---");
    check("Petersen", &Graph::petersen());

    println!("\nDone.");
}
