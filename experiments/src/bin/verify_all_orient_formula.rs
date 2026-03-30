/// Verify: S^{all}_{L(G)}(t) = 2^{|E(L(G))|} * mu_G(c_1(t-1), ..., c_m(t-1))
/// where c_i = 1/2^{deg_{L(G)}(e_i)} and mu_G is the matching polynomial.
///
/// We check this by:
/// 1. Computing S^{all}_{L(G)}(t) by brute force (all 2^m orientations of L(G))
/// 2. Computing it via the independence-set / matching formula
/// 3. Comparing

use combinatoric_core::graph::Graph;
use std::collections::BTreeSet;

/// Enumerate all matchings of G. Returns list of matchings,
/// each matching is a set of edge indices.
fn all_matchings(g: &Graph) -> Vec<Vec<usize>> {
    let edges = g.edges();
    let m = edges.len();
    let mut result = Vec::new();
    for mask in 0u64..(1u64 << m) {
        let mut vertices_used = BTreeSet::new();
        let mut matching = Vec::new();
        let mut is_matching = true;
        for i in 0..m {
            if mask & (1u64 << i) != 0 {
                let (u, v) = edges[i];
                if vertices_used.contains(&u) || vertices_used.contains(&v) {
                    is_matching = false;
                    break;
                }
                vertices_used.insert(u);
                vertices_used.insert(v);
                matching.push(i);
            }
        }
        if is_matching {
            result.push(matching);
        }
    }
    result
}

/// Compute S^{all}_G(t) by brute force: iterate over all 2^m orientations
fn brute_force_all_orient(g: &Graph) -> Vec<i64> {
    g.sink_polynomial_all_orientations()
}

/// Compute S^{all}_{L(G)}(t) via the matching formula.
/// S = 2^{|E(L)|} * sum_{M in matchings(G)} prod_{e in M} (t-1)/2^{deg_L(e)}
///
/// We compute this as a polynomial in t by expanding (t-1)^k for each matching of size k,
/// weighted appropriately.
fn matching_formula(g: &Graph) -> Vec<i64> {
    let lg = g.line_graph();
    let m_lg = lg.num_edges(); // |E(L(G))|
    let n_lg = lg.num_vertices(); // |V(L(G))| = |E(G)|

    // Degree of each vertex in L(G)
    let lg_degrees: Vec<usize> = (0..n_lg).map(|v| lg.degree(v)).collect();

    // All matchings of G (= independent sets of L(G))
    let matchings = all_matchings(g);

    // We compute using rational arithmetic to avoid floating point.
    // S = 2^{m_lg} * sum_M prod_{e in M} (t-1) / 2^{deg_L(e)}
    // = sum_M 2^{m_lg - sum_{e in M} deg_L(e)} * (t-1)^{|M|}
    //
    // Since 2^{m_lg - sum deg} might not be integer if sum deg > m_lg,
    // let's compute the LCM denominator approach.
    // Actually, each deg_L(e) <= n_lg, and sum deg <= m_lg always?
    // Not necessarily. Let's use exact integer arithmetic.
    //
    // Rewrite: S = sum_M 2^{m_lg - D_M} * (t-1)^{|M|}
    // where D_M = sum_{e in M} deg_L(e).
    // We need 2^{m_lg - D_M} to be an integer, i.e., m_lg >= D_M.
    //
    // Is this always true? For an independent set I of L(G),
    // sum_{v in I} deg(v) <= |E(L(G))| since the edges incident to I
    // are disjoint (I is independent) and are a subset of all edges.
    // Yes! Since I is independent, no edge of L(G) is incident to two vertices of I,
    // so the edges incident to I are counted exactly once each.
    // Thus sum deg_L(v) for v in I <= |E(L(G))|. ✓

    let max_matching_size = matchings.iter().map(|m| m.len()).max().unwrap_or(0);
    let max_deg = n_lg + 1; // safe upper bound for polynomial degree
    let mut result = vec![0i64; max_deg];

    for matching in &matchings {
        let k = matching.len();
        let d_m: usize = matching.iter().map(|&e| lg_degrees[e]).sum();
        let power_of_2 = 1i64 << (m_lg - d_m);

        // (t-1)^k expanded: sum_{j=0}^{k} C(k,j) * t^j * (-1)^{k-j}
        let mut binom = vec![0i64; k + 1];
        binom[0] = 1;
        for i in 0..k {
            let mut new_binom = vec![0i64; k + 1];
            for j in 0..=i {
                // multiply by (t - 1): new[j] += binom[j] * (-1), new[j+1] += binom[j] * 1
                new_binom[j] -= binom[j];
                new_binom[j + 1] += binom[j];
            }
            binom = new_binom;
        }

        for (j, &coeff) in binom.iter().enumerate() {
            result[j] += power_of_2 * coeff;
        }
    }

    // Strip trailing zeros
    while result.last() == Some(&0) {
        result.pop();
    }
    result
}

fn check(name: &str, g: &Graph) {
    let lg = g.line_graph();
    let m_lg = lg.num_edges();
    if m_lg > 20 {
        println!("{:<20} |E(L)|={} -- skipped", name, m_lg);
        return;
    }

    let brute = brute_force_all_orient(&lg);
    let formula = matching_formula(g);

    let ok = brute == formula;
    let status = if ok { "✓" } else { "MISMATCH !!!" };

    println!(
        "{:<20} brute={:?}  formula={:?}  [{}]",
        name, brute, formula, status
    );
}

fn main() {
    println!("=== Verifying all-orientations matching formula ===\n");

    for n in 3..=7 {
        check(&format!("P_{}", n), &Graph::path(n));
    }
    for n in 3..=8 {
        check(&format!("C_{}", n), &Graph::cycle(n));
    }
    for n in 3..=5 {
        check(&format!("K_{}", n), &Graph::complete(n));
    }
    for a in 1..=4 {
        for b in a..=4 {
            check(&format!("K_{{{},{}}}", a, b), &Graph::complete_bipartite(a, b));
        }
    }
    check("Petersen", &Graph::petersen());

    println!("\nDone.");
}
