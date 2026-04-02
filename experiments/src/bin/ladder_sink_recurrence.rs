use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::env;
use std::time::Instant;

use combinatoric_core::graph::Graph;
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};

type Poly = Vec<u128>;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct FrontierState {
    size: u8,
    out_bits: u8,
    reach_rows: [u8; 3],
}

fn ladder_graph(n: usize) -> Graph {
    assert!(n >= 2, "ladder graph requires n >= 2");

    let top = |i: usize| i;
    let bottom = |i: usize| n + i;

    let mut edges = Vec::new();

    for i in 0..n {
        edges.push((top(i), bottom(i)));
    }

    for i in 0..(n - 1) {
        edges.push((top(i), top(i + 1)));
        edges.push((bottom(i), bottom(i + 1)));
    }

    Graph::new(2 * n, &edges)
}

fn format_poly_u128(coeffs: &[u128]) -> String {
    let mut terms = Vec::new();
    for (k, &c) in coeffs.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let term = match k {
            0 => format!("{c}"),
            1 => {
                if c == 1 {
                    "t".to_string()
                } else {
                    format!("{c}t")
                }
            }
            _ => {
                if c == 1 {
                    format!("t^{k}")
                } else {
                    format!("{c}t^{k}")
                }
            }
        };
        terms.push(term);
    }
    if terms.is_empty() {
        "0".to_string()
    } else {
        terms.join(" + ")
    }
}

fn rows_get(rows: &[u8; 6], from: usize, to: usize) -> bool {
    (rows[from] >> to) & 1 == 1
}

fn rows_set(rows: &mut [u8; 6], from: usize, to: usize) {
    rows[from] |= 1 << to;
}

fn add_arc(rows: &mut [u8; 6], size: usize, from: usize, to: usize) -> bool {
    if from == to || rows_get(rows, to, from) {
        return false;
    }
    if rows_get(rows, from, to) {
        return true;
    }
    rows_set(rows, from, to);

    let mut pred_mask = 1u8 << from;
    for v in 0..size {
        if rows_get(rows, v, from) {
            pred_mask |= 1 << v;
        }
    }

    let succ_mask = rows[to] | (1u8 << to);
    for v in 0..size {
        if (pred_mask >> v) & 1 == 1 {
            rows[v] |= succ_mask;
        }
    }
    true
}

fn poly_add_assign(dst: &mut Poly, src: &[u128]) {
    if dst.len() < src.len() {
        dst.resize(src.len(), 0);
    }
    for (i, &c) in src.iter().enumerate() {
        dst[i] += c;
    }
    while dst.last() == Some(&0) {
        dst.pop();
    }
}

fn poly_mul(a: &[u128], b: &[u128]) -> Poly {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![0; a.len() + b.len() - 1];
    for (i, &ca) in a.iter().enumerate() {
        for (j, &cb) in b.iter().enumerate() {
            r[i + j] += ca * cb;
        }
    }
    while r.last() == Some(&0) {
        r.pop();
    }
    r
}

fn shift_poly_add(dst: &mut Poly, src: &[u128], shift: usize) {
    if dst.len() < src.len() + shift {
        dst.resize(src.len() + shift, 0);
    }
    for (i, &c) in src.iter().enumerate() {
        dst[i + shift] += c;
    }
    while dst.last() == Some(&0) {
        dst.pop();
    }
}

fn restricted_state(
    rows: &[u8; 6],
    out_bits: u8,
    old_size: usize,
    new_size: usize,
) -> FrontierState {
    let mut reach_rows = [0u8; 3];
    let mut new_out = 0u8;
    for i in 0..new_size {
        let src = old_size + i;
        if (out_bits >> src) & 1 == 1 {
            new_out |= 1 << i;
        }
        for j in 0..new_size {
            let dst = old_size + j;
            if rows_get(rows, src, dst) {
                reach_rows[i] |= 1 << j;
            }
        }
    }
    FrontierState {
        size: new_size as u8,
        out_bits: new_out,
        reach_rows,
    }
}

fn frontier_transition(
    dp: &BTreeMap<FrontierState, Poly>,
    old_size: usize,
    new_size: usize,
    step_edges: &[(usize, usize)],
) -> BTreeMap<FrontierState, Poly> {
    let union_size = old_size + new_size;
    let mut next = BTreeMap::new();

    for (state, poly) in dp {
        assert_eq!(state.size as usize, old_size);

        for mask in 0u32..(1u32 << step_edges.len()) {
            let mut rows = [0u8; 6];
            let mut out_bits = state.out_bits;

            for i in 0..old_size {
                rows[i] = state.reach_rows[i];
            }

            let mut ok = true;
            for (bit, &(u, v)) in step_edges.iter().enumerate() {
                let (from, to) = if (mask >> bit) & 1 == 0 {
                    (u, v)
                } else {
                    (v, u)
                };
                if !add_arc(&mut rows, union_size, from, to) {
                    ok = false;
                    break;
                }
                out_bits |= 1 << from;
            }
            if !ok {
                continue;
            }

            let forgotten_sinks = (0..old_size)
                .filter(|&v| ((out_bits >> v) & 1) == 0)
                .count();
            let key = restricted_state(&rows, out_bits, old_size, new_size);
            let entry = next.entry(key).or_insert_with(Vec::new);
            shift_poly_add(entry, poly, forgotten_sinks);
        }
    }

    next
}

fn frontier_finalize(dp: &BTreeMap<FrontierState, Poly>) -> Poly {
    let mut total = Vec::new();
    for (state, poly) in dp {
        let sinks = (0..state.size as usize)
            .filter(|&v| ((state.out_bits >> v) & 1) == 0)
            .count();
        shift_poly_add(&mut total, poly, sinks);
    }
    total
}

fn ladder_line_graph_sink_frontier(n: usize) -> Poly {
    assert!(n >= 2, "ladder graph requires n >= 2");

    let init = FrontierState {
        size: 1,
        out_bits: 0,
        reach_rows: [0, 0, 0],
    };
    let mut dp = BTreeMap::new();
    dp.insert(init, vec![1]);

    let first_step = [(0usize, 2usize), (0, 3), (1, 2), (1, 3)];
    dp = frontier_transition(&dp, 1, 3, &first_step);

    let generic_step = [(0usize, 4usize), (0, 5), (3, 4), (3, 5), (1, 4), (2, 5)];
    for _ in 2..n {
        dp = frontier_transition(&dp, 3, 3, &generic_step);
    }

    frontier_finalize(&dp)
}

fn first_step_distribution() -> BTreeMap<FrontierState, Poly> {
    let init = FrontierState {
        size: 1,
        out_bits: 0,
        reach_rows: [0, 0, 0],
    };
    let mut dp = BTreeMap::new();
    dp.insert(init, vec![1]);
    let first_step = [(0usize, 2usize), (0, 3), (1, 2), (1, 3)];
    frontier_transition(&dp, 1, 3, &first_step)
}

fn build_detailed_transfer() -> (Vec<FrontierState>, Vec<Poly>, Vec<Vec<Poly>>, Vec<Poly>) {
    let generic_step = [(0usize, 4usize), (0, 5), (3, 4), (3, 5), (1, 4), (2, 5)];
    let start_map = first_step_distribution();

    let mut seen = BTreeSet::new();
    let mut queue = VecDeque::new();
    for &state in start_map.keys() {
        seen.insert(state);
        queue.push_back(state);
    }

    while let Some(state) = queue.pop_front() {
        let mut single = BTreeMap::new();
        single.insert(state, vec![1]);
        let next = frontier_transition(&single, 3, 3, &generic_step);
        for &next_state in next.keys() {
            if seen.insert(next_state) {
                queue.push_back(next_state);
            }
        }
    }

    let states: Vec<FrontierState> = seen.into_iter().collect();
    let mut index = BTreeMap::new();
    for (i, &state) in states.iter().enumerate() {
        index.insert(state, i);
    }

    let mut start = vec![Vec::new(); states.len()];
    for (state, poly) in &start_map {
        start[index[state]] = poly.clone();
    }

    let mut matrix = vec![vec![Vec::new(); states.len()]; states.len()];
    for &state in &states {
        let mut single = BTreeMap::new();
        single.insert(state, vec![1]);
        let next = frontier_transition(&single, 3, 3, &generic_step);
        for (next_state, poly) in next {
            matrix[index[&state]][index[&next_state]] = poly;
        }
    }

    let mut term = vec![Vec::new(); states.len()];
    for &state in &states {
        let mut single = BTreeMap::new();
        single.insert(state, vec![1]);
        term[index[&state]] = frontier_finalize(&single);
    }

    (states, start, matrix, term)
}

fn apply_transfer(vec: &[Poly], matrix: &[Vec<Poly>]) -> Vec<Poly> {
    let mut next = vec![Vec::new(); vec.len()];
    for (i, row) in matrix.iter().enumerate() {
        for (j, entry) in row.iter().enumerate() {
            if vec[i].is_empty() || entry.is_empty() {
                continue;
            }
            let term = poly_mul(&vec[i], entry);
            poly_add_assign(&mut next[j], &term);
        }
    }
    next
}

fn ladder_line_graph_sink_transfer(
    n: usize,
    start: &[Poly],
    matrix: &[Vec<Poly>],
    term: &[Poly],
) -> Poly {
    assert!(n >= 2, "ladder graph requires n >= 2");
    let mut vec = start.to_vec();
    for _ in 0..(n - 2) {
        vec = apply_transfer(&vec, matrix);
    }
    let mut total = Vec::new();
    for i in 0..vec.len() {
        if vec[i].is_empty() || term[i].is_empty() {
            continue;
        }
        let contribution = poly_mul(&vec[i], &term[i]);
        poly_add_assign(&mut total, &contribution);
    }
    total
}

fn relation_signature(state: FrontierState) -> String {
    let r_a = (state.reach_rows[0] >> 1) & 1 == 1;
    let a_r = (state.reach_rows[1] >> 0) & 1 == 1;
    let r_b = (state.reach_rows[0] >> 2) & 1 == 1;
    let b_r = (state.reach_rows[2] >> 0) & 1 == 1;
    let a_b = (state.reach_rows[1] >> 2) & 1 == 1;
    let b_a = (state.reach_rows[2] >> 1) & 1 == 1;

    let rel = match (r_a, a_r, r_b, b_r, a_b, b_a) {
        (true, false, true, false, false, false) => "r<a,b",
        (true, false, true, false, true, false) => "r<a<b",
        (true, false, true, false, false, true) => "r<b<a",
        (false, true, true, false, true, false) => "a<r<b",
        (true, false, false, true, false, true) => "b<r<a",
        (false, true, false, true, true, false) => "a<b<r",
        (false, true, false, true, false, true) => "b<a<r",
        (false, true, false, true, false, false) => "a,b<r",
        _ => "other",
    };

    format!(
        "{} ; out={:03b} ; rows=[{:03b},{:03b},{:03b}]",
        rel, state.out_bits, state.reach_rows[0], state.reach_rows[1], state.reach_rows[2]
    )
}

fn print_transfer_data(
    states: &[FrontierState],
    start: &[Poly],
    matrix: &[Vec<Poly>],
    term: &[Poly],
) {
    println!("=== Frontier transfer data ===");
    println!("states: {}", states.len());
    for (i, state) in states.iter().enumerate() {
        println!("state {i}: {}", relation_signature(*state));
    }

    println!("\ninitial:");
    for (i, poly) in start.iter().enumerate() {
        if !poly.is_empty() {
            println!("  state {i} -> {}", format_poly_u128(poly));
        }
    }

    println!("\nterminal:");
    for (i, poly) in term.iter().enumerate() {
        println!("  state {i} -> {}", format_poly_u128(poly));
    }

    println!("\nnonzero transfer entries:");
    for i in 0..states.len() {
        for j in 0..states.len() {
            if !matrix[i][j].is_empty() {
                println!("  {i} -> {j} : {}", format_poly_u128(&matrix[i][j]));
            }
        }
    }
}

fn as_i64_polys(polys: &[Poly]) -> Vec<Vec<i64>> {
    let mut result = Vec::new();
    for poly in polys {
        let mut row = Vec::with_capacity(poly.len());
        for &c in poly {
            if c > i64::MAX as u128 {
                return result;
            }
            row.push(c as i64);
        }
        result.push(row);
    }
    result
}

fn run_search(name: &str, polys: &[Vec<i64>], opts: AdaptiveSearchOptions) {
    let t0 = Instant::now();
    print!("{name}: ");
    match find_recurrence_adaptive(polys, &opts) {
        Some(res) => {
            println!(
                "FOUND in {:?} [{} unknowns, {} equations, {} candidates]",
                t0.elapsed(),
                res.num_unknowns,
                res.num_equations,
                res.candidates_tried
            );
            println!("  {}", res.recurrence);
        }
        None => {
            println!("none in {:?}", t0.elapsed());
        }
    }
}

fn main() {
    let mut max_n = 12usize;
    let mut verify_up_to = 12usize;
    let mut print_matrix = true;

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--max-n" => {
                let value = args.next().expect("missing value after --max-n");
                max_n = value.parse().expect("--max-n must be an integer");
            }
            "--verify-up-to" => {
                let value = args.next().expect("missing value after --verify-up-to");
                verify_up_to = value.parse().expect("--verify-up-to must be an integer");
            }
            "--no-matrix" => {
                print_matrix = false;
            }
            _ => {
                eprintln!(
                    "usage: ladder_sink_recurrence [--max-n N] [--verify-up-to N] [--no-matrix]"
                );
                std::process::exit(1);
            }
        }
    }

    let (states, start, matrix, term) = build_detailed_transfer();
    if print_matrix {
        print_transfer_data(&states, &start, &matrix, &term);
        println!();
    }

    let mut polys = Vec::new();

    println!("=== Acyclic sink polynomials for L(Ladder_n) via frontier transfer ===\n");
    for n in 2..=max_n {
        let t0 = Instant::now();
        let coeffs = ladder_line_graph_sink_transfer(n, &start, &matrix, &term);
        let elapsed = t0.elapsed();

        print!("n={n}: coeffs={coeffs:?} [{elapsed:?}]");
        if n <= verify_up_to {
            let frontier = ladder_line_graph_sink_frontier(n);
            if frontier != coeffs {
                println!("\n  frontier mismatch: {frontier:?}");
                std::process::exit(1);
            }
            if n <= 5 {
                let ladder = ladder_graph(n);
                let brute = ladder.line_graph().acyclic_sink_polynomial();
                let brute_u128: Vec<u128> = brute.iter().map(|&c| c as u128).collect();
                if brute_u128 != coeffs {
                    println!("\n  brute force mismatch: {brute:?}");
                    std::process::exit(1);
                }
            }
            print!(" verified");
        }
        println!();
        println!("  {}", format_poly_u128(&coeffs));
        polys.push(coeffs);
    }

    println!("\n=== CSV for polytool ===");
    for p in &polys {
        let row = p
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join(", ");
        println!("{row}");
    }

    let i64_polys = as_i64_polys(&polys);
    if i64_polys.len() < 4 {
        println!("\n=== Recurrence search ===");
        println!("not enough i64-convertible terms for a useful search");
        return;
    }

    println!(
        "\n=== Recurrence search on first {} i64-convertible terms ===",
        i64_polys.len()
    );
    run_search(
        "constant in n, homogeneous",
        &i64_polys,
        AdaptiveSearchOptions {
            min_rec_len: 1,
            max_rec_len: 6,
            min_var_deg: 0,
            max_var_deg: 4,
            min_idx_deg: 0,
            max_idx_deg: 0,
            min_diff_deg: 0,
            max_diff_deg: 0,
            min_margin: 3,
            ..AdaptiveSearchOptions::default()
        },
    );
}
