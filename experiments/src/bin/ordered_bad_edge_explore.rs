use std::collections::BTreeMap;

use combinatoric_core::graph::Graph;
use combinatoric_core::partition::Partition;
use sym_poly_sym::{chromatic_symmetric, first_bad_edge_symmetric, Basis, SymmetricFunction};

#[derive(Clone)]
struct OrderedCase {
    label: String,
    graph: Graph,
    orders: Vec<(&'static str, Vec<(usize, usize)>)>,
}

fn main() {
    let cases = vec![
        OrderedCase {
            label: "P3".to_string(),
            graph: Graph::path(3),
            orders: vec![("natural", Graph::path(3).edges().to_vec())],
        },
        OrderedCase {
            label: "paw".to_string(),
            graph: Graph::new(4, &[(0, 1), (0, 2), (1, 2), (2, 3)]),
            orders: vec![
                ("triangle-first", vec![(0, 1), (0, 2), (1, 2), (2, 3)]),
                ("leaf-first", vec![(2, 3), (0, 1), (0, 2), (1, 2)]),
            ],
        },
        OrderedCase {
            label: "diamond".to_string(),
            graph: Graph::new(4, &[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]),
            orders: vec![
                ("natural", vec![(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]),
                ("outer-first", vec![(1, 3), (2, 3), (0, 1), (0, 2), (1, 2)]),
            ],
        },
    ];

    for case in &cases {
        run_case(case);
    }

    if let Some(case) = find_order_sensitive_case(5) {
        run_case(&case);
    } else {
        println!("no order-sensitive connected example found up to 5 vertices");
    }
}

fn run_case(case: &OrderedCase) {
    println!("============================================================");
    println!(
        "{} on {} vertices with edges {}",
        case.label,
        case.graph.num_vertices(),
        format_edges(case.graph.edges())
    );

    let mut computed = Vec::new();
    for (order_label, order) in &case.orders {
        let report = analyze_order(&case.graph, order_label, order);
        computed.push(report);
    }

    if computed.len() >= 2 {
        println!();
        println!("order comparison:");
        for idx in 1..computed.len() {
            let lhs = &computed[0];
            let rhs = &computed[idx];
            let diff = lhs.function.clone() - rhs.function.clone();
            if diff.is_zero() {
                println!(
                    "  {} and {} give the same symmetric function",
                    lhs.label, rhs.label
                );
            } else {
                println!(
                    "  {} - {} = {}",
                    lhs.label,
                    rhs.label,
                    diff.to_monomial_basis()
                );
            }
        }
    }

    println!();
}

fn find_order_sensitive_case(max_vertices: usize) -> Option<OrderedCase> {
    for n in 4..=max_vertices {
        let all_edges = Graph::complete(n).edges().to_vec();
        let total_graphs = 1usize << all_edges.len();

        for mask in 0..total_graphs {
            let edges: Vec<_> = all_edges
                .iter()
                .enumerate()
                .filter_map(|(idx, &edge)| ((mask >> idx) & 1 == 1).then_some(edge))
                .collect();

            if edges.len() < 2 {
                continue;
            }

            let graph = Graph::new(n, &edges);
            if !graph.is_connected() {
                continue;
            }

            let natural = graph.edges().to_vec();
            let reverse: Vec<_> = natural.iter().rev().copied().collect();
            if natural == reverse {
                continue;
            }

            let f_natural = first_bad_edge_symmetric::<i64>(&graph, &natural);
            let f_reverse = first_bad_edge_symmetric::<i64>(&graph, &reverse);
            if f_natural != f_reverse {
                return Some(OrderedCase {
                    label: format!(
                        "first order-sensitive connected graph found up to {} vertices",
                        n
                    ),
                    graph,
                    orders: vec![("natural", natural), ("reverse", reverse)],
                });
            }
        }
    }

    None
}

struct OrderReport {
    label: String,
    function: SymmetricFunction<i64>,
}

fn analyze_order(g: &Graph, order_label: &str, order: &[(usize, usize)]) -> OrderReport {
    println!();
    println!("order {}: {}", order_label, format_edges(order));

    let via_formula = first_bad_edge_symmetric::<i64>(g, order).to_monomial_basis();
    let via_direct = direct_first_bad_edge(g, order);
    assert_eq!(
        via_formula, via_direct,
        "direct enumeration disagrees for order {}",
        order_label
    );

    println!("  f_(G, Ω) in m-basis: {}", via_formula);
    println!(
        "  f_(G, Ω) in e-basis: {}",
        via_formula.to_elementary_basis()
    );

    println!("  decomposition Σ_e X_(G_e):");
    for (idx, &(u, v)) in order.iter().enumerate() {
        let prefix_edges: Vec<_> = order[..idx].iter().copied().map(normalize_edge).collect();
        let prefix_graph = Graph::new(g.num_vertices(), &prefix_edges);
        let ge = prefix_graph.contract_edge(u, v);
        let xge = chromatic_symmetric::<i64>(&ge).to_monomial_basis();
        println!(
            "    e{} = {}: prefix edges {}, G_e has {} vertices / edges {}, X_(G_e) = {}",
            idx + 1,
            format_edge((u, v)),
            format_edges(&prefix_edges),
            ge.num_vertices(),
            format_edges(ge.edges()),
            xge
        );
    }

    println!("  direct coefficient check by partition:");
    for lambda in Partition::all_of_size(g.num_vertices().saturating_sub(1) as u32) {
        let coeff = via_formula.coefficient(&lambda);
        if coeff != 0 {
            println!("    coeff of m[{}] = {}", lambda, coeff);
        }
    }

    println!("  trivial specializations f(1^k):");
    for k in 1..=4 {
        let spec = via_formula.trivial_specialization(k);
        let total = (k as i64).pow(g.num_vertices() as u32);
        let proper = g.chromatic_polynomial_eval(k as usize) as i64;
        let nonproper = total - proper;
        println!(
            "    k = {}: f(1^k) = {}, non-proper {}^{} - χ_G({}) = {}",
            k,
            spec,
            k,
            g.num_vertices(),
            k,
            nonproper
        );
        assert_eq!(spec, nonproper);
    }

    OrderReport {
        label: order_label.to_string(),
        function: via_formula,
    }
}

fn direct_first_bad_edge(g: &Graph, order: &[(usize, usize)]) -> SymmetricFunction<i64> {
    let degree = g.num_vertices().saturating_sub(1) as u32;
    let mut terms = BTreeMap::new();

    for lambda in Partition::all_of_size(degree) {
        let coeff = direct_coefficient(g, order, &lambda);
        if coeff != 0 {
            terms.insert(lambda, coeff);
        }
    }

    SymmetricFunction::from_terms(Basis::Monomial, terms)
}

fn direct_coefficient(g: &Graph, order: &[(usize, usize)], lambda: &Partition) -> i64 {
    let target: Vec<usize> = lambda.parts().iter().map(|&part| part as usize).collect();
    let mut coeff = 0;

    for chosen_color in 0..target.len() {
        let mut word = Vec::with_capacity(g.num_vertices());
        for (color, &count) in target.iter().enumerate() {
            let frequency = count + usize::from(color == chosen_color);
            for _ in 0..frequency {
                word.push(color);
            }
        }
        word.sort_unstable();

        loop {
            if let Some((u, _v)) = first_bad_edge(&word, order) {
                if word[u] == chosen_color {
                    coeff += 1;
                }
            }
            if !next_multiset_perm(&mut word) {
                break;
            }
        }
    }

    coeff
}

fn first_bad_edge(coloring: &[usize], order: &[(usize, usize)]) -> Option<(usize, usize)> {
    order
        .iter()
        .copied()
        .map(normalize_edge)
        .find(|&(u, v)| coloring[u] == coloring[v])
}

fn next_multiset_perm(perm: &mut [usize]) -> bool {
    let n = perm.len();
    if n <= 1 {
        return false;
    }

    let mut i = n - 2;
    loop {
        if perm[i] < perm[i + 1] {
            break;
        }
        if i == 0 {
            return false;
        }
        i -= 1;
    }

    let mut j = n - 1;
    while perm[j] <= perm[i] {
        j -= 1;
    }
    perm.swap(i, j);
    perm[i + 1..].reverse();
    true
}

fn normalize_edge((u, v): (usize, usize)) -> (usize, usize) {
    if u < v {
        (u, v)
    } else {
        (v, u)
    }
}

fn format_edge(edge: (usize, usize)) -> String {
    let (u, v) = normalize_edge(edge);
    format!("({},{})", u, v)
}

fn format_edges(edges: &[(usize, usize)]) -> String {
    if edges.is_empty() {
        return "[]".to_string();
    }
    let pieces: Vec<_> = edges.iter().copied().map(format_edge).collect();
    format!("[{}]", pieces.join(", "))
}
