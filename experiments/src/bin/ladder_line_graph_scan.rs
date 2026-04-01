use std::env;

use combinatoric_core::chromatic::chromatic_symmetric;
use combinatoric_core::graph::Graph;

fn main() {
    let mut max_n = 5usize;

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--max-n" => {
                let value = args.next().expect("missing value after --max-n");
                max_n = value.parse().expect("--max-n must be an integer");
            }
            _ => {
                eprintln!(
                    "usage: cargo run -p experiments --bin ladder_line_graph_scan -- [--max-n N]"
                );
                std::process::exit(1);
            }
        }
    }

    for n in 2..=max_n {
        let ladder = ladder_graph(n);
        let line_graph = ladder.line_graph();
        let x = chromatic_symmetric::<i64>(&line_graph);
        let e = x.to_elementary_basis();
        let sum: i64 = e.terms().values().copied().sum();

        println!(
            "Ladder n={}: G has {} vertices / {} edges; L(G) has {} vertices / {} edges",
            n,
            ladder.num_vertices(),
            ladder.num_edges(),
            line_graph.num_vertices(),
            line_graph.num_edges()
        );
        println!("  X_(L(G)) in e-basis: {}", e);
        println!("  sum of e-coefficients: {}", sum);

        if e.positive_coefficients() {
            println!("  e-positive");
        } else {
            println!("  not e-positive");
            println!("  negative e-terms:");
            for (partition, coeff) in e.terms() {
                if *coeff < 0 {
                    println!("    coeff of e[{}] = {}", partition, coeff);
                }
            }
            return;
        }
    }
}

fn ladder_graph(n: usize) -> Graph {
    assert!(n >= 2, "ladder graph requires n >= 2");

    let top = |i: usize| i;
    let bottom = |i: usize| n + i;

    let mut edges = Vec::new();

    for i in 0..n {
        edges.push((top(i), bottom(i)));
    }

    for i in 0..(n - 1) {
        edges.push((top(i), top(i + 1)));
        edges.push((bottom(i), bottom(i + 1)));
    }

    Graph::new(2 * n, &edges)
}
