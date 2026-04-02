use std::env;
use std::path::PathBuf;

use combinatoric_core::graph::Graph;
use sym_poly_sym::chromatic_symmetric;

fn main() {
    let mut graph_dir = PathBuf::from("/home/paxinum/Dropbox/mathematica-packages");
    let mut max_n = 12usize;
    let mut caterpillars_only = false;

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--max-n" => {
                let value = args.next().expect("missing value after --max-n");
                max_n = value.parse().expect("--max-n must be an integer");
            }
            "--graph-dir" => {
                let value = args.next().expect("missing value after --graph-dir");
                graph_dir = PathBuf::from(value);
            }
            "--caterpillars" => caterpillars_only = true,
            _ => {
                eprintln!(
                    "usage: cargo run -p experiments --bin tree_line_graph_chromatic_scan -- [--max-n N] [--graph-dir DIR] [--caterpillars]"
                );
                std::process::exit(1);
            }
        }
    }

    for n in 1..=max_n {
        let mut trees = trees_on_n(n, &graph_dir);
        if caterpillars_only {
            trees.retain(is_caterpillar);
        }
        if trees.is_empty() {
            let family = if caterpillars_only {
                "caterpillars"
            } else {
                "trees"
            };
            println!("n = {}: no {}", n, family);
            continue;
        }

        let family = if caterpillars_only {
            "caterpillars"
        } else {
            "trees"
        };
        println!("n = {}: {} {}", n, trees.len(), family);

        let mut e_positive = 0usize;
        let mut first_failure = None;

        for (idx, tree) in trees.iter().enumerate() {
            let line_graph = tree.line_graph();
            let x = chromatic_symmetric::<i64>(&line_graph);
            let e = x.to_elementary_basis();

            if e.positive_coefficients() {
                e_positive += 1;
                continue;
            }

            first_failure = Some((idx + 1, tree.clone(), line_graph, x, e));
            break;
        }

        if let Some((index, tree, line_graph, x, e)) = first_failure {
            println!(
                "  first non-e-positive line graph comes from tree #{}",
                index
            );
            println!("  tree edges: {:?}", tree.edges());
            println!("  line graph edges: {:?}", line_graph.edges());
            println!("  X_(L(T)) in e-basis: {}", e);
            println!("  X_(L(T)) in s-basis: {}", x.to_schur_basis());
            println!("  negative e-terms:");
            for (partition, coeff) in e.terms() {
                if *coeff < 0 {
                    println!("    coeff of e[{}] = {}", partition, coeff);
                }
            }
            return;
        }

        println!("  all {} line graphs are e-positive", e_positive);
    }

    println!(
        "no e-negativity found for chromatic symmetric functions of line graphs of {} up to n = {}",
        if caterpillars_only {
            "caterpillars"
        } else {
            "trees"
        },
        max_n
    );
}

fn trees_on_n(n: usize, graph_dir: &PathBuf) -> Vec<Graph> {
    match n {
        1 => vec![Graph::new(1, &[])],
        2 => vec![Graph::path(2)],
        3 => vec![Graph::path(3)],
        4 => vec![Graph::path(4), Graph::new(4, &[(0, 1), (0, 2), (0, 3)])],
        _ => {
            let path = graph_dir.join(format!("trees{}.g6", n));
            if !path.exists() {
                return vec![];
            }
            Graph::all_from_graph6_file(&path)
                .unwrap_or_else(|err| panic!("failed to read {}: {}", path.display(), err))
        }
    }
}

fn is_caterpillar(tree: &Graph) -> bool {
    if tree.num_vertices() <= 2 {
        return true;
    }

    let core: Vec<_> = (0..tree.num_vertices())
        .filter(|&v| tree.degree(v) != 1)
        .collect();
    if core.len() <= 1 {
        return true;
    }

    let spine = tree.induced_subgraph(&core);
    spine.is_connected() && (0..spine.num_vertices()).all(|v| spine.degree(v) <= 2)
}
