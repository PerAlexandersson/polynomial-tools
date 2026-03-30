/// Explore interlacing relations in the tree sink polynomial recursion.
///
/// For each tree T, compute A_v(t) and B_v(t) at each vertex, and check:
/// 1. Are A_v, B_v, A_v+B_v all real-rooted?
/// 2. Does A_v interlace B_v? (or vice versa)
/// 3. Does A_v interlace A_v+B_v?
/// 4. Does t·A_v + B_v have any special interlacing with A_v+B_v?
/// 5. Under leaf deletion: does S_{L(T')} interlace S_{L(T)}?

use combinatoric_core::graph::Graph;
use num_bigint::BigInt;
use num_traits::{Signed, ToPrimitive, Zero};
use polynomial_tools::real_rootedness::{
    check_interlacing, check_weak_interlacing, is_real_rooted,
};

/// Convert BigInt poly to i64 by dividing out GCD. Returns None if too large.
fn to_i64_poly(p: &[BigInt]) -> Option<Vec<i64>> {
    if p.is_empty() {
        return Some(vec![]);
    }
    let mut g = BigInt::from(0);
    for c in p {
        if !c.is_zero() {
            g = num::integer::gcd(g.abs(), c.abs());
        }
    }
    if g.is_zero() {
        return Some(vec![0; p.len()]);
    }
    p.iter()
        .map(|c| (c / &g).to_i64())
        .collect()
}

fn is_rr_big(p: &[BigInt]) -> bool {
    match to_i64_poly(p) {
        Some(v) => is_real_rooted(&v),
        None => {
            eprintln!("  WARNING: coefficients too large for real-rootedness check");
            true // assume true
        }
    }
}

fn interlaces_big(f: &[BigInt], g: &[BigInt]) -> Option<bool> {
    let f64 = to_i64_poly(f)?;
    let g64 = to_i64_poly(g)?;
    check_interlacing(&f64, &g64)
}

fn weak_interlaces_big(f: &[BigInt], g: &[BigInt]) -> Option<bool> {
    let f64 = to_i64_poly(f)?;
    let g64 = to_i64_poly(g)?;
    check_weak_interlacing(&f64, &g64)
}

fn bpoly_add(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let len = a.len().max(b.len());
    let mut r = vec![BigInt::from(0); len];
    for (i, c) in a.iter().enumerate() { r[i] += c; }
    for (i, c) in b.iter().enumerate() { r[i] += c; }
    while r.last().map_or(false, |c| c.is_zero()) { r.pop(); }
    r
}

fn bpoly_mul_t(p: &[BigInt]) -> Vec<BigInt> {
    if p.is_empty() { return vec![]; }
    let mut r = vec![BigInt::from(0); p.len() + 1];
    for (i, c) in p.iter().enumerate() { r[i + 1] = c.clone(); }
    r
}

/// Generate tree from Prüfer sequence.
fn prufer_to_tree(seq: &[usize], n: usize) -> Graph {
    let mut degree = vec![1usize; n];
    for &v in seq { degree[v] += 1; }
    let mut edges = Vec::new();
    let mut seq_iter = seq.iter();
    for _ in 0..n - 2 {
        let &v = seq_iter.next().unwrap();
        for u in 0..n {
            if degree[u] == 1 {
                edges.push((u, v));
                degree[u] -= 1;
                degree[v] -= 1;
                break;
            }
        }
    }
    let remaining: Vec<usize> = (0..n).filter(|&v| degree[v] == 1).collect();
    edges.push((remaining[0], remaining[1]));
    Graph::new(n, &edges)
}

fn main() {
    println!("=== Interlacing exploration for tree sink polynomials ===\n");

    // Counters for relations
    let mut total_vertices = 0u64;
    let mut a_rr = 0u64;
    let mut b_rr = 0u64;
    let mut ab_rr = 0u64;
    let mut tab_rr = 0u64; // t*A + B
    let mut a_interlaces_b = 0u64;        // A ≪ B (strict)
    let mut a_weak_interlaces_b = 0u64;   // A ≪ B (weak)
    let mut a_interlaces_ab = 0u64;       // A ≪ A+B
    let mut a_weak_interlaces_ab = 0u64;
    let mut b_weak_interlaces_tab = 0u64; // B ≪ tA+B
    let mut a_weak_interlaces_tab = 0u64; // A ≪ tA+B
    let mut ab_weak_interlaces_tab = 0u64; // A+B ≪ tA+B

    // Also check leaf deletion interlacing
    let mut leaf_del_count = 0u64;
    let mut leaf_del_interlace = 0u64;
    let mut leaf_del_weak_interlace = 0u64;

    for n in 2..=9 {
        let mut seen = std::collections::BTreeSet::new();

        if n == 2 {
            process_tree(&Graph::path(2),
                &mut total_vertices, &mut a_rr, &mut b_rr, &mut ab_rr, &mut tab_rr,
                &mut a_interlaces_b, &mut a_weak_interlaces_b,
                &mut a_interlaces_ab, &mut a_weak_interlaces_ab,
                &mut b_weak_interlaces_tab, &mut a_weak_interlaces_tab,
                &mut ab_weak_interlaces_tab,
                &mut leaf_del_count, &mut leaf_del_interlace, &mut leaf_del_weak_interlace,
            );
            continue;
        }

        let seq_len = n - 2;
        let total = (n as u64).pow(seq_len as u32);
        for code in 0..total {
            let mut seq = vec![0usize; seq_len];
            let mut c = code;
            for i in (0..seq_len).rev() {
                seq[i] = (c % n as u64) as usize;
                c /= n as u64;
            }
            let g = prufer_to_tree(&seq, n);
            let mut ds: Vec<usize> = (0..n).map(|v| g.degree(v)).collect();
            ds.sort();
            if !seen.insert(ds) { continue; }

            process_tree(&g,
                &mut total_vertices, &mut a_rr, &mut b_rr, &mut ab_rr, &mut tab_rr,
                &mut a_interlaces_b, &mut a_weak_interlaces_b,
                &mut a_interlaces_ab, &mut a_weak_interlaces_ab,
                &mut b_weak_interlaces_tab, &mut a_weak_interlaces_tab,
                &mut ab_weak_interlaces_tab,
                &mut leaf_del_count, &mut leaf_del_interlace, &mut leaf_del_weak_interlace,
            );
        }

        println!("n={}: done", n);
    }

    println!("\n=== Summary ===");
    println!("Total non-root vertices examined: {}", total_vertices);
    println!();
    println!("Real-rootedness:");
    println!("  A_v real-rooted:     {}/{}", a_rr, total_vertices);
    println!("  B_v real-rooted:     {}/{}", b_rr, total_vertices);
    println!("  A_v+B_v real-rooted: {}/{}", ab_rr, total_vertices);
    println!("  t*A_v+B_v real-rooted: {}/{}", tab_rr, total_vertices);
    println!();
    println!("Interlacing (at non-root non-leaf vertices):");
    println!("  A_v ≪ B_v (strict):         {}/{}", a_interlaces_b, total_vertices);
    println!("  A_v ≪ B_v (weak):           {}/{}", a_weak_interlaces_b, total_vertices);
    println!("  A_v ≪ A_v+B_v (strict):     {}/{}", a_interlaces_ab, total_vertices);
    println!("  A_v ≪ A_v+B_v (weak):       {}/{}", a_weak_interlaces_ab, total_vertices);
    println!("  B_v ≪ t*A_v+B_v (weak):     {}/{}", b_weak_interlaces_tab, total_vertices);
    println!("  A_v ≪ t*A_v+B_v (weak):     {}/{}", a_weak_interlaces_tab, total_vertices);
    println!("  A_v+B_v ≪ t*A_v+B_v (weak): {}/{}", ab_weak_interlaces_tab, total_vertices);
    println!();
    println!("Leaf deletion interlacing:");
    println!("  S_{{L(T')}} ≪ S_{{L(T)}} (strict): {}/{}", leaf_del_interlace, leaf_del_count);
    println!("  S_{{L(T')}} ≪ S_{{L(T)}} (weak):   {}/{}", leaf_del_weak_interlace, leaf_del_count);
}

#[allow(clippy::too_many_arguments)]
fn process_tree(
    g: &Graph,
    total_vertices: &mut u64,
    a_rr: &mut u64, b_rr: &mut u64, ab_rr: &mut u64, tab_rr: &mut u64,
    a_interlaces_b: &mut u64, a_weak_interlaces_b: &mut u64,
    a_interlaces_ab: &mut u64, a_weak_interlaces_ab: &mut u64,
    b_weak_interlaces_tab: &mut u64, a_weak_interlaces_tab: &mut u64,
    ab_weak_interlaces_tab: &mut u64,
    leaf_del_count: &mut u64,
    leaf_del_interlace: &mut u64,
    leaf_del_weak_interlace: &mut u64,
) {
    let n = g.num_vertices();
    let (order, parent, a_poly, b_poly) = g.tree_ab_polynomials();
    let root = order[0];

    // Check per-vertex properties (non-root, non-leaf)
    for &v in &order {
        if v == root { continue; }
        // Skip leaves (A=1, B=0, trivial)
        if a_poly[v] == vec![BigInt::from(1)] && b_poly[v].is_empty() {
            continue;
        }

        *total_vertices += 1;

        let av = &a_poly[v];
        let bv = &b_poly[v];
        let abv = bpoly_add(av, bv);
        let tav = bpoly_mul_t(av);
        let tabv = bpoly_add(&tav, bv);

        if is_rr_big(av) { *a_rr += 1; }
        if is_rr_big(bv) { *b_rr += 1; }
        if is_rr_big(&abv) { *ab_rr += 1; }
        if is_rr_big(&tabv) { *tab_rr += 1; }

        if let Some(true) = interlaces_big(av, bv) { *a_interlaces_b += 1; }
        if let Some(true) = weak_interlaces_big(av, bv) { *a_weak_interlaces_b += 1; }
        if let Some(true) = interlaces_big(av, &abv) { *a_interlaces_ab += 1; }
        if let Some(true) = weak_interlaces_big(av, &abv) { *a_weak_interlaces_ab += 1; }
        if let Some(true) = weak_interlaces_big(bv, &tabv) { *b_weak_interlaces_tab += 1; }
        if let Some(true) = weak_interlaces_big(av, &tabv) { *a_weak_interlaces_tab += 1; }
        if let Some(true) = weak_interlaces_big(&abv, &tabv) { *ab_weak_interlaces_tab += 1; }
    }

    // Leaf deletion interlacing: for each leaf, delete it and compare sink polynomials
    if n >= 3 {
        let s_full = g.sink_polynomial_tree_bigint();
        for v in 0..n {
            if g.degree(v) != 1 { continue; }
            // Build T' = T - v
            let verts: Vec<usize> = (0..n).filter(|&u| u != v).collect();
            let t_prime = g.induced_subgraph(&verts);
            let s_prime = t_prime.sink_polynomial_tree_bigint();

            *leaf_del_count += 1;
            if let Some(true) = interlaces_big(&s_prime, &s_full) {
                *leaf_del_interlace += 1;
            }
            if let Some(true) = weak_interlaces_big(&s_prime, &s_full) {
                *leaf_del_weak_interlace += 1;
            }
        }
    }
}
