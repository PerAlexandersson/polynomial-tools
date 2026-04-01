/// Debug chromatic polynomial DC on P_2
use combinatoric_core::graph::Graph;

fn main() {
    // P_2: single edge 0-1
    let g = Graph::path(2);
    println!(
        "G = P_2: {} vertices, edges: {:?}",
        g.num_vertices(),
        g.edges()
    );

    // Delete edge (0,1): empty graph on 2 vertices -> chi = t^2 = [0, 0, 1]
    let g_del = g.delete_edge(0, 1);
    println!(
        "G-e: {} vertices, edges: {:?}",
        g_del.num_vertices(),
        g_del.edges()
    );
    let chi_del = g_del.chromatic_polynomial();
    println!("chi(G-e) = {:?}", chi_del);

    // Contract edge (0,1): single vertex, no edges -> chi = t = [0, 1]
    let g_con = g.contract_edge(0, 1);
    println!(
        "G/e: {} vertices, edges: {:?}",
        g_con.num_vertices(),
        g_con.edges()
    );
    let chi_con = g_con.chromatic_polynomial();
    println!("chi(G/e) = {:?}", chi_con);

    // chi(P_2) = chi(G-e) - chi(G/e) = t^2 - t = [0, -1, 1]
    let chi = g.chromatic_polynomial();
    println!("chi(P_2) = {:?}", chi);
    // Expected: t(t-1) = t^2 - t = [0, -1, 1]
}
