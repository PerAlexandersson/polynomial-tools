use std::env;

use combinatoric_core::graph::Graph;
use sym_poly_sym::chromatic_symmetric;

fn main() {
    let mut max_m = 4usize;
    let mut max_n = 6usize;

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--max-m" => {
                let value = args.next().expect("missing value after --max-m");
                max_m = value.parse().expect("--max-m must be an integer");
            }
            "--max-n" => {
                let value = args.next().expect("missing value after --max-n");
                max_n = value.parse().expect("--max-n must be an integer");
            }
            _ => {
                eprintln!(
                    "usage: cargo run -p experiments --bin rectangle_line_graph_scan -- [--max-m M] [--max-n N]"
                );
                std::process::exit(1);
            }
        }
    }

    for m in 1..=max_m {
        for n in m..=max_n {
            let g = Graph::complete_bipartite(m, n);
            let lg = g.line_graph();
            let x = chromatic_symmetric::<i64>(&lg);
            let e = x.to_elementary_basis();

            println!(
                "K_{{{}, {}}}: L(G) has {} vertices and {} edges",
                m,
                n,
                lg.num_vertices(),
                lg.num_edges()
            );

            if e.positive_coefficients() {
                println!("  e-positive");
            } else {
                println!("  not e-positive");
                println!("  X_(L(K_{{{}, {}}})) in e-basis: {}", m, n, e);
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

    println!("no e-negativity found for rectangles K_(m,n) in the scanned range");
}
