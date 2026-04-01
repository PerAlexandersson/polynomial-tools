use std::collections::BTreeMap;
use std::env;
use std::path::PathBuf;

use combinatoric_core::chromatic::{chromatic_symmetric, first_bad_edge_symmetric};
use combinatoric_core::graph::Graph;
use combinatoric_core::symmetric_function::{Basis, SymmetricFunction};

fn main() {
    let mut graph_dir = PathBuf::from("/home/paxinum/Dropbox/mathematica-packages");
    let mut max_n = 12usize;

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
            _ => {
                eprintln!(
                    "usage: cargo run -p experiments --bin ordered_bad_edge_tree_scan -- [--max-n N] [--graph-dir DIR]"
                );
                std::process::exit(1);
            }
        }
    }

    verify_recurrence_and_formula();

    for n in 1..=max_n {
        let trees = trees_on_n(n, &graph_dir);
        if trees.is_empty() {
            println!("n = {}: missing tree data", n);
            continue;
        }

        println!("n = {}: {} trees", n, trees.len());

        let mut schur_positive = 0usize;
        let mut first_failure = None;

        for (idx, raw_tree) in trees.iter().enumerate() {
            let tree = sorted_graph(raw_tree);
            let f = first_bad_edge_recurrence(&tree);
            let schur = f.to_schur_basis();

            if schur.positive_coefficients() {
                schur_positive += 1;
                continue;
            }

            first_failure = Some((idx + 1, tree, f, schur));
            break;
        }

        if let Some((index, tree, f, schur)) = first_failure {
            println!(
                "  first Schur-negative tree is #{} with edges {:?}",
                index,
                tree.edges()
            );
            println!("  f_(T, Ω) in e-basis: {}", f.to_elementary_basis());
            println!("  f_(T, Ω) in s-basis: {}", schur);
            println!("  negative Schur terms:");
            for (partition, coeff) in schur.terms() {
                if *coeff < 0 {
                    println!("    coeff of s[{}] = {}", partition, coeff);
                }
            }
            return;
        }

        println!("  all {} trees are Schur-positive", schur_positive);
    }

    println!("no Schur-negativity found for trees up to n = {}", max_n);
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

fn verify_recurrence_and_formula() {
    let samples = [
        Graph::path(3),
        Graph::path(4),
        Graph::new(5, &[(0, 1), (0, 2), (0, 3), (0, 4)]),
    ];

    for raw in samples {
        let tree = sorted_graph(&raw);
        let order = tree.edges().to_vec();
        let direct = first_bad_edge_symmetric::<i64>(&tree, &order);
        let recursive = first_bad_edge_recurrence(&tree);
        assert_eq!(
            direct,
            recursive,
            "recurrence mismatch on {:?}",
            tree.edges()
        );
    }
}

fn first_bad_edge_recurrence(tree: &Graph) -> SymmetricFunction<i64> {
    let mut memo = BTreeMap::new();
    first_bad_edge_recurrence_memo(tree, &mut memo)
}

fn first_bad_edge_recurrence_memo(
    graph: &Graph,
    memo: &mut BTreeMap<String, SymmetricFunction<i64>>,
) -> SymmetricFunction<i64> {
    let graph = sorted_graph(graph);
    let key = graph_key(&graph);
    if let Some(value) = memo.get(&key) {
        return value.clone();
    }

    let value = match graph.edges().last().copied() {
        None => SymmetricFunction::zero(Basis::Monomial),
        Some((u, v)) => {
            let contracted = sorted_graph(&graph.contract_edge(u, v));
            let deleted = sorted_graph(&graph.delete_edge(u, v));
            chromatic_symmetric::<i64>(&contracted) + first_bad_edge_recurrence_memo(&deleted, memo)
        }
    };

    memo.insert(key, value.clone());
    value
}

fn sorted_graph(graph: &Graph) -> Graph {
    let mut edges = graph.edges().to_vec();
    edges.sort_unstable();
    Graph::new(graph.num_vertices(), &edges)
}

fn graph_key(graph: &Graph) -> String {
    format!("{}:{:?}", graph.num_vertices(), graph.edges())
}
