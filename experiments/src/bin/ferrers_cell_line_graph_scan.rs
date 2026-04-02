use std::collections::BTreeMap;
use std::env;

use combinatoric_core::graph::Graph;
use combinatoric_core::partition::Partition;
use sym_poly_sym::chromatic_symmetric;

fn main() {
    let mut max_cells = 8u32;
    let mut print_all = false;

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--max-cells" => {
                let value = args.next().expect("missing value after --max-cells");
                max_cells = value.parse().expect("--max-cells must be an integer");
            }
            "--all" => print_all = true,
            _ => {
                eprintln!(
                    "usage: cargo run -p experiments --bin ferrers_cell_line_graph_scan -- [--max-cells N] [--all]"
                );
                std::process::exit(1);
            }
        }
    }

    for cells in 1..=max_cells {
        let shapes = Partition::all_of_size(cells);
        println!("cells = {}: {} Ferrers shapes", cells, shapes.len());

        let mut e_positive = 0usize;
        let mut first_failure = None;

        for shape in &shapes {
            let g = ferrers_cell_graph(shape);
            let lg = g.line_graph();
            let x = chromatic_symmetric::<i64>(&lg);
            let e = x.to_elementary_basis();

            if print_all {
                let sum: i64 = e.terms().values().copied().sum();
                println!(
                    "  shape {}: |V(L)|={}, |E(L)|={}, sum={}, e-positive={}",
                    shape,
                    lg.num_vertices(),
                    lg.num_edges(),
                    sum,
                    e.positive_coefficients()
                );
            }

            if e.positive_coefficients() {
                e_positive += 1;
                continue;
            }

            first_failure = Some((shape.clone(), g, lg, x, e));
            break;
        }

        if let Some((shape, g, lg, x, e)) = first_failure {
            println!("  first non-e-positive shape: {}", shape);
            println!("  cell graph edges: {:?}", g.edges());
            println!("  line graph edges: {:?}", lg.edges());
            println!("  X_(L(G)) in e-basis: {}", e);
            println!("  X_(L(G)) in s-basis: {}", x.to_schur_basis());
            println!("  negative e-terms:");
            for (partition, coeff) in e.terms() {
                if *coeff < 0 {
                    println!("    coeff of e[{}] = {}", partition, coeff);
                }
            }
            return;
        }

        println!("  all {} shapes are e-positive", e_positive);
    }

    println!(
        "no e-negativity found for Ferrers cell-adjacency graphs up to {} cells",
        max_cells
    );
}

fn ferrers_cell_graph(shape: &Partition) -> Graph {
    let boxes = shape.diagram_boxes();
    let mut index = BTreeMap::new();
    for (i, &cell) in boxes.iter().enumerate() {
        index.insert(cell, i);
    }

    let mut edges = Vec::new();
    for &(r, c) in &boxes {
        if let Some(&j) = index.get(&(r + 1, c)) {
            edges.push((index[&(r, c)], j));
        }
        if let Some(&j) = index.get(&(r, c + 1)) {
            edges.push((index[&(r, c)], j));
        }
    }

    Graph::new(boxes.len(), &edges)
}
