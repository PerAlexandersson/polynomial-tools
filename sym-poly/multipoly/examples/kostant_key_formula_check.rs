use std::collections::BTreeMap;

use sym_poly_multipoly::key_polynomial::key_polynomial;

type Weight = Vec<i32>;
type Root = Vec<i32>;

fn dominant_rearrangement(alpha: &[u32]) -> Vec<i32> {
    let mut lambda: Vec<u32> = alpha.to_vec();
    lambda.sort_unstable_by(|a, b| b.cmp(a));
    lambda.into_iter().map(|x| x as i32).collect()
}

fn sorting_reduced_word(alpha: &[u32]) -> Vec<usize> {
    let mut perm = alpha.to_vec();
    let mut word = Vec::new();
    let n = perm.len();
    if n < 2 {
        return word;
    }
    loop {
        let mut swapped = false;
        for i in 0..n - 1 {
            if perm[i] < perm[i + 1] {
                perm.swap(i, i + 1);
                word.push(i);
                swapped = true;
            }
        }
        if !swapped {
            break;
        }
    }
    word
}

fn simple_root(n: usize, i: usize) -> Root {
    let mut root = vec![0; n];
    root[i] = 1;
    root[i + 1] = -1;
    root
}

fn apply_simple(v: &mut Weight, i: usize) {
    v.swap(i, i + 1);
}

fn sub(a: &Weight, b: &Weight) -> Weight {
    a.iter().zip(b).map(|(x, y)| x - y).collect()
}

fn canonical_positive_root(beta: &Root) -> (Root, bool) {
    let plus = beta.iter().position(|&x| x == 1).unwrap();
    let minus = beta.iter().position(|&x| x == -1).unwrap();
    if plus < minus {
        (beta.clone(), true)
    } else {
        let mut gamma = beta.clone();
        for x in &mut gamma {
            *x = -*x;
        }
        (gamma, false)
    }
}

fn root_height(v: &Weight) -> i32 {
    -v.iter()
        .enumerate()
        .map(|(i, x)| (i as i32) * x)
        .sum::<i32>()
}

fn kostant_partition(roots: &[Root], target: &Weight) -> i64 {
    if target.iter().sum::<i32>() != 0 {
        return 0;
    }
    let bound = root_height(target);
    if bound < 0 {
        return 0;
    }
    fn rec(
        idx: usize,
        roots: &[Root],
        target: Weight,
        bound: i32,
        memo: &mut BTreeMap<(usize, Weight), i64>,
    ) -> i64 {
        if idx == roots.len() {
            return if target.iter().all(|&x| x == 0) { 1 } else { 0 };
        }
        if let Some(value) = memo.get(&(idx, target.clone())) {
            return *value;
        }
        let mut total = 0;
        for m in 0..=bound {
            let mut next = target.clone();
            for (x, r) in next.iter_mut().zip(&roots[idx]) {
                *x -= m * r;
            }
            if root_height(&next) <= bound - m {
                total += rec(idx + 1, roots, next, bound - m, memo);
            }
        }
        memo.insert((idx, target), total);
        total
    }
    rec(0, roots, target.clone(), bound, &mut BTreeMap::new())
}

fn key_coeff_kostant(alpha: &[u32], mu: &[u32]) -> i64 {
    let n = alpha.len();
    let lambda = dominant_rearrangement(alpha);
    let mut word = sorting_reduced_word(alpha);
    word.reverse();
    let mu_i: Weight = mu.iter().map(|&x| x as i32).collect();
    let mut total = 0;

    for (sign, base, roots) in demazure_kostant_terms(n, &lambda, &word) {
        let target = sub(&base, &mu_i);
        total += sign * kostant_partition(&roots, &target);
    }

    total
}

fn canonicalize_term(sign: &mut i64, base: &mut Weight, roots: &mut Vec<Root>) {
    for root in roots {
        let (gamma, positive) = canonical_positive_root(root);
        if positive {
            *root = gamma;
        } else {
            *sign = -*sign;
            for (x, g) in base.iter_mut().zip(&gamma) {
                *x -= *g;
            }
            *root = gamma;
        }
    }
}

fn demazure_kostant_terms(
    n: usize,
    lambda: &Weight,
    application_word: &[usize],
) -> Vec<(i64, Weight, Vec<Root>)> {
    let mut terms = vec![(1i64, lambda.clone(), Vec::<Root>::new())];
    for &i in application_word {
        let alpha_i = simple_root(n, i);
        let mut next_terms = Vec::new();
        for (sign, base, roots) in terms {
            // Identity branch: F/(1-e^{-alpha_i}).
            let mut roots0 = roots.clone();
            roots0.push(alpha_i.clone());
            next_terms.push((sign, base.clone(), roots0));

            // Reflection branch: - e^{-alpha_i} s_i(F)/(1-e^{-alpha_i}).
            let mut sign1 = -sign;
            let mut base1 = base.clone();
            apply_simple(&mut base1, i);
            for (x, a) in base1.iter_mut().zip(&alpha_i) {
                *x -= *a;
            }
            let mut roots1 = roots.clone();
            for root in &mut roots1 {
                apply_simple(root, i);
            }
            roots1.push(alpha_i.clone());
            canonicalize_term(&mut sign1, &mut base1, &mut roots1);
            next_terms.push((sign1, base1, roots1));
        }
        terms = next_terms;
    }
    terms
}

fn compositions_of_degree(n: usize, degree: u32) -> Vec<Vec<u32>> {
    fn rec(n: usize, degree: u32, prefix: &mut Vec<u32>, out: &mut Vec<Vec<u32>>) {
        if prefix.len() == n - 1 {
            let mut item = prefix.clone();
            item.push(degree);
            out.push(item);
            return;
        }
        for x in 0..=degree {
            prefix.push(x);
            rec(n, degree - x, prefix, out);
            prefix.pop();
        }
    }
    let mut out = Vec::new();
    rec(n, degree, &mut Vec::new(), &mut out);
    out
}

fn check_alpha(alpha: &[u32]) {
    let key = key_polynomial::<i64>(alpha);
    let degree: u32 = alpha.iter().sum();
    let mut kostant_terms = BTreeMap::new();
    for mu in compositions_of_degree(alpha.len(), degree) {
        let coeff = key_coeff_kostant(alpha, &mu);
        if coeff != 0 {
            kostant_terms.insert(mu, coeff);
        }
    }
    assert_eq!(
        key.terms(),
        &kostant_terms,
        "Kostant formula mismatch for alpha={alpha:?}"
    );
    println!("alpha={alpha:?}: verified {} terms", key.terms().len());
    for (mu, coeff) in key.terms() {
        println!(
            "  x^{mu:?}: key={coeff}, kostant={}",
            key_coeff_kostant(alpha, mu)
        );
    }
}

fn main() {
    for alpha in [
        vec![0, 2],
        vec![1, 0, 2],
        vec![0, 1, 2],
        vec![2, 0, 1],
        vec![1, 2, 0],
    ] {
        check_alpha(&alpha);
    }
}
