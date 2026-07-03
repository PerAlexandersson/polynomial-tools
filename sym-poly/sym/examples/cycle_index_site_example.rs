use std::collections::BTreeMap;

use num_rational::Ratio;
use sym_poly_core::Partition;
use sym_poly_sym::{Basis, SymmetricFunction};

type Q = Ratio<i64>;

fn p(parts: &[u32]) -> Partition {
    Partition::new(parts.to_vec())
}

fn q(n: i64, d: i64) -> Q {
    Ratio::new(n, d)
}

fn permutations(n: usize) -> Vec<Vec<usize>> {
    fn backtrack(current: &mut Vec<usize>, unused: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if unused.is_empty() {
            out.push(current.clone());
            return;
        }
        for idx in 0..unused.len() {
            let value = unused.remove(idx);
            current.push(value);
            backtrack(current, unused, out);
            current.pop();
            unused.insert(idx, value);
        }
    }

    let mut out = Vec::new();
    let mut current = Vec::new();
    let mut unused = (0..n).collect();
    backtrack(&mut current, &mut unused, &mut out);
    out
}

fn cycle_type(permutation: &[usize]) -> Vec<u32> {
    let mut visited = vec![false; permutation.len()];
    let mut parts = Vec::new();
    for start in 0..permutation.len() {
        if visited[start] {
            continue;
        }
        let mut length = 0;
        let mut current = start;
        loop {
            visited[current] = true;
            length += 1;
            current = permutation[current];
            if current == start {
                break;
            }
        }
        parts.push(length);
    }
    parts.sort_unstable_by(|a, b| b.cmp(a));
    parts
}

fn cycle_index_symmetric_group(n: usize) -> SymmetricFunction<Q> {
    let perms = permutations(n);
    let mut counts = BTreeMap::<Partition, i64>::new();
    for permutation in &perms {
        *counts
            .entry(Partition::new(cycle_type(permutation)))
            .or_insert(0) += 1;
    }

    let denominator = perms.len() as i64;
    let terms = counts
        .into_iter()
        .map(|(partition, count)| (partition, q(count, denominator)))
        .collect();
    SymmetricFunction::from_terms(Basis::PowerSum, terms)
}

fn main() {
    let cycle_index = cycle_index_symmetric_group(3);
    assert_eq!(cycle_index.coefficient(&p(&[1, 1, 1])), q(1, 6));
    assert_eq!(cycle_index.coefficient(&p(&[2, 1])), q(1, 2));
    assert_eq!(cycle_index.coefficient(&p(&[3])), q(1, 3));

    let schur = cycle_index.to_schur_basis();
    assert_eq!(schur.coefficient(&p(&[3])), q(1, 1));
    assert_eq!(schur.coefficient(&p(&[2, 1])), q(0, 1));
    assert_eq!(schur.coefficient(&p(&[1, 1, 1])), q(0, 1));

    println!("Z(S_3) = 1/6 p_111 + 1/2 p_21 + 1/3 p_3 = s_3");
}
