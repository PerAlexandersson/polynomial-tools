use std::collections::BTreeMap;
use std::env;
use std::path::PathBuf;

use combinatoric_core::graph::Graph;
use sym_poly_sym::{chromatic_symmetric, first_bad_edge_symmetric, Basis, SymmetricFunction};

#[derive(Clone, Copy, PartialEq, Eq)]
enum Family {
    Connected,
    Trees,
}

fn main() {
    let mut graph_dir = PathBuf::from("/home/paxinum/Dropbox/mathematica-packages");
    let mut max_n = 10usize;
    let mut family = Family::Connected;

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
            "--trees" => family = Family::Trees,
            "--connected" => family = Family::Connected,
            _ => {
                eprintln!(
                    "usage: cargo run -p experiments --bin ordered_bad_edge_e_scan -- [--connected|--trees] [--max-n N] [--graph-dir DIR]"
                );
                std::process::exit(1);
            }
        }
    }

    verify_recurrence_and_formula();

    for n in family_start(family)..=max_n {
        let graphs = match family {
            Family::Connected => connected_graphs_on_n(n, &graph_dir),
            Family::Trees => trees_on_n(n, &graph_dir),
        };

        if graphs.is_empty() {
            println!("n = {}: missing data", n);
            continue;
        }

        let family_name = match family {
            Family::Connected => "connected graphs",
            Family::Trees => "trees",
        };
        println!("n = {}: {} {}", n, graphs.len(), family_name);

        let mut e_positive = 0usize;
        let mut first_failure = None;

        for (idx, raw_graph) in graphs.iter().enumerate() {
            let graph = sorted_graph(raw_graph);
            let f = first_bad_edge_recurrence(&graph);
            let e = f.to_elementary_basis();

            if e.positive_coefficients() {
                e_positive += 1;
                continue;
            }

            first_failure = Some((idx + 1, graph, f, e));
            break;
        }

        if let Some((index, graph, f, e)) = first_failure {
            println!(
                "  first non-e-positive object is #{} with edges {:?}",
                index,
                graph.edges()
            );
            println!("  f_(G, Ω) in e-basis: {}", e);
            println!("  f_(G, Ω) in s-basis: {}", f.to_schur_basis());
            println!("  negative e-terms:");
            for (partition, coeff) in e.terms() {
                if *coeff < 0 {
                    println!("    coeff of e[{}] = {}", partition, coeff);
                }
            }
            return;
        }

        println!("  all {} objects are e-positive", e_positive);
    }

    let family_name = match family {
        Family::Connected => "connected graphs",
        Family::Trees => "trees",
    };
    println!(
        "no e-negativity found for {} up to n = {}",
        family_name, max_n
    );
}

fn family_start(family: Family) -> usize {
    match family {
        Family::Connected => 2,
        Family::Trees => 1,
    }
}

fn connected_graphs_on_n(n: usize, graph_dir: &PathBuf) -> Vec<Graph> {
    let path = graph_dir.join(format!("graph{}c.g6", n));
    if !path.exists() {
        return vec![];
    }
    Graph::all_from_graph6_file(&path)
        .unwrap_or_else(|err| panic!("failed to read {}: {}", path.display(), err))
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
        Graph::cycle(4),
        Graph::new(5, &[(0, 1), (0, 3), (0, 4), (1, 2)]),
    ];

    for raw in samples {
        let graph = sorted_graph(&raw);
        let order = graph.edges().to_vec();
        let direct = first_bad_edge_symmetric::<i64>(&graph, &order);
        let recursive = first_bad_edge_recurrence(&graph);
        assert_eq!(
            direct,
            recursive,
            "recurrence mismatch on {:?}",
            graph.edges()
        );
    }
}

fn first_bad_edge_recurrence(graph: &Graph) -> SymmetricFunction<i64> {
    let mut memo = BTreeMap::new();
    first_bad_edge_recurrence_memo(graph, &mut memo)
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
