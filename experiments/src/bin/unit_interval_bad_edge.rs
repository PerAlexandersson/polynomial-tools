use std::env;

use combinatoric_core::graph::Graph;
use sym_poly_sym::first_bad_edge_symmetric;

fn main() {
    let mut args = env::args().skip(1);
    let Some(area_arg) = args.next() else {
        eprintln!(
            "usage: cargo run -p experiments --bin unit_interval_bad_edge -- <area-sequence>"
        );
        eprintln!("example: cargo run -p experiments --bin unit_interval_bad_edge -- 012331");
        std::process::exit(1);
    };

    let area = parse_area_sequence(&area_arg);
    let mut order = Graph::unit_interval(&area).edges().to_vec();
    order.sort_unstable();
    let graph = Graph::new(area.len(), &order);
    let f = first_bad_edge_symmetric::<i64>(&graph, &order);

    println!("area sequence: {}", area_arg);
    println!("vertices: {}", graph.num_vertices());
    println!("lex edge order Ω: {:?}", order);
    println!("f_(G, Ω) in m-basis: {}", f.to_monomial_basis());
    println!("f_(G, Ω) in e-basis: {}", f.to_elementary_basis());
}

fn parse_area_sequence(s: &str) -> Vec<u8> {
    assert!(!s.is_empty(), "area sequence must be nonempty");
    s.chars()
        .map(|ch| {
            ch.to_digit(10)
                .expect("area sequence must use decimal digits") as u8
        })
        .collect()
}
