use combinatoric_core::key_polynomial::key_ehrhart_polynomial;
use combinatoric_core::partition::Partition;
use combpoly::permutation::all_permutations;
use polynomial_tools::ehrhart_to_hstar;
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
struct ChainData {
    labels: Vec<(u8, u8)>,
    weight: u64,
}

#[derive(Clone, Debug)]
struct MatchRecord {
    sigma: Vec<u8>,
    target: Vec<i64>,
    chain_poly: Vec<u64>,
}

fn inv_count(perm: &[u8]) -> usize {
    let mut c = 0;
    for i in 0..perm.len() {
        for j in i + 1..perm.len() {
            if perm[i] > perm[j] {
                c += 1;
            }
        }
    }
    c
}

fn bruhat_covers_down_with_labels(sigma: &[u8]) -> Vec<(Vec<u8>, (u8, u8))> {
    let n = sigma.len();
    let mut covers = Vec::new();
    for a in 0..n {
        for b in a + 1..n {
            if sigma[a] > sigma[b] {
                let lo = sigma[b];
                let hi = sigma[a];
                let blocked = (a + 1..b).any(|c| sigma[c] > lo && sigma[c] < hi);
                if !blocked {
                    let mut tau = sigma.to_vec();
                    tau.swap(a, b);
                    covers.push((tau, ((a + 1) as u8, (b + 1) as u8)));
                }
            }
        }
    }
    covers
}

fn enumerate_chains_down(
    sigma: &[u8],
    id: &[u8],
    suffix_labels: &mut Vec<(u8, u8)>,
    suffix_weight: u64,
    out: &mut Vec<ChainData>,
) {
    if sigma == id {
        let mut labels = suffix_labels.clone();
        labels.reverse();
        out.push(ChainData {
            labels,
            weight: suffix_weight,
        });
        return;
    }

    for (tau, label) in bruhat_covers_down_with_labels(sigma) {
        suffix_labels.push(label);
        enumerate_chains_down(
            &tau,
            id,
            suffix_labels,
            suffix_weight * (label.1 - label.0) as u64,
            out,
        );
        suffix_labels.pop();
    }
}

fn all_saturated_chains(sigma: &[u8]) -> Vec<ChainData> {
    let id: Vec<u8> = (1..=sigma.len() as u8).collect();
    let mut out = Vec::new();
    let mut suffix_labels = Vec::new();
    enumerate_chains_down(sigma, &id, &mut suffix_labels, 1, &mut out);
    out
}

fn trim_trailing_zeros(mut v: Vec<i64>) -> Vec<i64> {
    while v.len() > 1 && v.last() == Some(&0) {
        v.pop();
    }
    v
}

fn next_permutation<T: Ord>(a: &mut [T]) -> bool {
    if a.len() < 2 {
        return false;
    }
    let mut i = a.len() - 2;
    loop {
        if a[i] < a[i + 1] {
            break;
        }
        if i == 0 {
            return false;
        }
        i -= 1;
    }
    let mut j = a.len() - 1;
    while a[j] <= a[i] {
        j -= 1;
    }
    a.swap(i, j);
    a[i + 1..].reverse();
    true
}

fn weighted_chain_descent_poly(
    chains: &[ChainData],
    label_ranks: &BTreeMap<(u8, u8), usize>,
) -> Vec<u64> {
    let mut coeffs = Vec::<u64>::new();
    for chain in chains {
        let mut desc = 0usize;
        for w in chain.labels.windows(2) {
            let left = label_ranks[&w[0]];
            let right = label_ranks[&w[1]];
            if left > right {
                desc += 1;
            }
        }
        if coeffs.len() <= desc {
            coeffs.resize(desc + 1, 0);
        }
        coeffs[desc] += chain.weight;
    }
    if coeffs.is_empty() {
        coeffs.push(1);
    }
    coeffs
}

fn perm_str_u8(perm: &[u8]) -> String {
    perm.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join("")
}

fn perm_str_usize(perm: &[(u8, u8)]) -> String {
    perm.iter()
        .map(|(i, j)| format!("({},{})", i, j))
        .collect::<Vec<_>>()
        .join(" < ")
}

fn main() {
    let lambda = Partition::new(vec![3, 2, 1]);
    let n = 4u8;
    let perms = all_permutations(n);

    let mut target_hstar = BTreeMap::<Vec<u8>, Vec<i64>>::new();
    let mut chains_by_sigma = BTreeMap::<Vec<u8>, Vec<ChainData>>::new();

    println!("Computing key h* and saturated Bruhat chains for lambda=(3,2,1,0), n=4...");
    for perm in &perms {
        let sigma_usize: Vec<usize> = perm.iter().map(|&x| x as usize).collect();
        let ep = key_ehrhart_polynomial(&lambda, &sigma_usize, None);
        let hstar = trim_trailing_zeros(ehrhart_to_hstar(&ep.coeffs));
        let chains = all_saturated_chains(perm);
        target_hstar.insert(perm.clone(), hstar);
        chains_by_sigma.insert(perm.clone(), chains);
    }

    let mut labels = vec![(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)];
    let mut exact_orders = Vec::<Vec<(u8, u8)>>::new();
    let mut best_score = 0usize;
    let mut best_examples = Vec::<(Vec<(u8, u8)>, Vec<MatchRecord>)>::new();
    let mut tested = 0usize;

    loop {
        tested += 1;
        let label_ranks: BTreeMap<(u8, u8), usize> = labels
            .iter()
            .enumerate()
            .map(|(i, &lab)| (lab, i))
            .collect();

        let mut matches = Vec::<MatchRecord>::new();
        for perm in &perms {
            let chain_poly = weighted_chain_descent_poly(&chains_by_sigma[perm], &label_ranks);
            let target = &target_hstar[perm];
            if chain_poly.iter().map(|&x| x as i64).collect::<Vec<_>>() == *target {
                matches.push(MatchRecord {
                    sigma: perm.clone(),
                    target: target.clone(),
                    chain_poly,
                });
            }
        }

        if matches.len() == perms.len() {
            exact_orders.push(labels.clone());
        } else {
            if matches.len() > best_score {
                best_score = matches.len();
                best_examples.clear();
            }
            if matches.len() == best_score && best_examples.len() < 5 {
                best_examples.push((labels.clone(), matches));
            }
        }

        if !next_permutation(&mut labels) {
            break;
        }
    }

    println!("Tested {} total orders on transposition labels.", tested);
    println!();

    if exact_orders.is_empty() {
        println!(
            "No reflection order on the six labels gives an exact match for all 24 permutations."
        );
        println!("Best score: {} exact matches out of 24.", best_score);
        println!();
        println!("Best orders found:");
        for (order, matches) in &best_examples {
            println!("  order: {}", perm_str_usize(order));
            let mut by_length = BTreeMap::<usize, Vec<String>>::new();
            for rec in matches {
                by_length
                    .entry(inv_count(&rec.sigma))
                    .or_default()
                    .push(perm_str_u8(&rec.sigma));
            }
            for (len, sigmas) in by_length {
                println!("    inv {}: {}", len, sigmas.join(", "));
            }
        }
    } else {
        println!(
            "Found {} exact reflection orders matching the trimmed key h* for all 24 permutations:",
            exact_orders.len()
        );
        for order in exact_orders {
            println!("  {}", perm_str_usize(&order));
        }
    }

    println!();
    println!("Selected comparisons for the lex order (1,2)<(1,3)<(1,4)<(2,3)<(2,4)<(3,4):");
    let lex_order = vec![(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)];
    let lex_ranks: BTreeMap<(u8, u8), usize> = lex_order
        .iter()
        .enumerate()
        .map(|(i, &lab)| (lab, i))
        .collect();
    for perm in &perms {
        let chain_poly = weighted_chain_descent_poly(&chains_by_sigma[perm], &lex_ranks);
        let target = &target_hstar[perm];
        if chain_poly.iter().map(|&x| x as i64).collect::<Vec<_>>() == *target {
            println!(
                "  {}: MATCH  {}",
                perm_str_u8(perm),
                format!("{:?}", target)
            );
        }
    }
}
