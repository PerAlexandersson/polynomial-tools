use combinatoric_core::graph::Graph;
use sym_poly_core::Partition;
use sym_poly_sym::{Basis, SymmetricFunction};

type Sym = SymmetricFunction<i64>;

fn partition(parts: &[u32]) -> Partition {
    Partition::new(parts.to_vec())
}

fn set_partitions(n: usize) -> Vec<Vec<Vec<usize>>> {
    let mut partitions = vec![Vec::<Vec<usize>>::new()];
    for vertex in 0..n {
        let mut next = Vec::new();
        for blocks in partitions {
            for block_index in 0..blocks.len() {
                let mut extended = blocks.clone();
                extended[block_index].push(vertex);
                next.push(extended);
            }
            let mut with_singleton = blocks;
            with_singleton.push(vec![vertex]);
            next.push(with_singleton);
        }
        partitions = next;
    }
    partitions
}

fn bond_partitions(graph: &Graph) -> Vec<Vec<Vec<usize>>> {
    set_partitions(graph.num_vertices())
        .into_iter()
        .filter(|blocks| {
            blocks
                .iter()
                .all(|block| graph.induced_subgraph(block).is_connected())
        })
        .collect()
}

fn one() -> Sym {
    Sym::complete_h_symmetric(Partition::empty())
}

/// Compute M_G from the recurrence in González D'León--Wachs, Theorem 5.7.
fn mobius_symmetric(graph: &Graph) -> Sym {
    assert!(graph.is_connected());
    if graph.num_vertices() == 1 {
        return one();
    }

    let mut sum = Sym::zero(Basis::CompleteH);
    for blocks in bond_partitions(graph)
        .into_iter()
        .filter(|blocks| blocks.len() > 1)
    {
        let h = Sym::complete_h_symmetric(partition(&[(blocks.len() - 1) as u32]));
        let product = blocks.iter().fold(one(), |acc, block| {
            acc.multiply(&mobius_symmetric(&graph.induced_subgraph(block)))
        });
        sum = sum + h.multiply(&product);
    }
    -sum
}

/// Compute Psi_G = sum_{pi in Pi_G} M_{G|pi}.
fn chromatic_mobius_symmetric(graph: &Graph) -> Sym {
    bond_partitions(graph)
        .into_iter()
        .fold(Sym::zero(Basis::CompleteH), |sum, blocks| {
            let product = blocks.iter().fold(one(), |acc, block| {
                acc.multiply(&mobius_symmetric(&graph.induced_subgraph(block)))
            });
            sum + product
        })
}

fn assert_coefficient(function: &Sym, parts: &[u32], expected: i64) {
    assert_eq!(function.coefficient(&partition(parts)), expected);
}

fn main() {
    let path = Graph::path(3);
    let complete = Graph::complete(3);

    let m_path = mobius_symmetric(&path).to_monomial_basis();
    let m_complete = mobius_symmetric(&complete).to_monomial_basis();
    assert_coefficient(&m_path, &[2], 1);
    assert_coefficient(&m_path, &[1, 1], 3);
    assert_coefficient(&m_complete, &[2], 2);
    assert_coefficient(&m_complete, &[1, 1], 5);

    let psi_path = chromatic_mobius_symmetric(&path).to_elementary_basis();
    let psi_complete = chromatic_mobius_symmetric(&complete).to_elementary_basis();
    for (parts, expected) in [
        (&[][..], 1),
        (&[1][..], -2),
        (&[2][..], 1),
        (&[1, 1][..], 1),
    ] {
        assert_coefficient(&psi_path, parts, expected);
    }
    for (parts, expected) in [
        (&[][..], 1),
        (&[1][..], -3),
        (&[2][..], 1),
        (&[1, 1][..], 2),
    ] {
        assert_coefficient(&psi_complete, parts, expected);
    }

    println!("M_P3 in the m-basis: {m_path}");
    println!("M_K3 in the m-basis: {m_complete}");
    println!("Psi_P3 in the e-basis: {psi_path}");
    println!("Psi_K3 in the e-basis: {psi_complete}");
}
