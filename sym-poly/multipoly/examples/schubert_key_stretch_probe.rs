use std::collections::BTreeMap;

use combinatoric_core::next_permutation;
use sym_poly_multipoly::{
    diagram_weight, format_diagram, is_yamanouchi, key_polynomial, kohnert_diagrams, rothe_diagram,
    schubert_polynomial, schubert_to_key, MultiPoly,
};

fn lehmer_code(perm: &[usize]) -> Vec<usize> {
    let n = perm.len();
    let mut code = vec![0; n];
    for i in 0..n {
        for j in i + 1..n {
            if perm[j] < perm[i] {
                code[i] += 1;
            }
        }
    }
    code
}

fn from_lehmer_code(code: &[usize]) -> Vec<usize> {
    let mut available: Vec<usize> = (1..=code.len()).collect();
    let mut perm = Vec::with_capacity(code.len());
    for &c in code {
        perm.push(available.remove(c));
    }
    perm
}

fn stretch_perm(perm: &[usize], dilation: usize) -> Vec<usize> {
    let code = lehmer_code(perm);
    let needed_len = code
        .iter()
        .enumerate()
        .filter(|(_, c)| **c > 0)
        .map(|(i, c)| i + 1 + dilation * c)
        .max()
        .unwrap_or(perm.len())
        .max(perm.len());
    let mut stretched_code = vec![0; needed_len];
    for (i, c) in code.iter().enumerate() {
        stretched_code[i] = dilation * c;
    }
    from_lehmer_code(&stretched_code)
}

fn key_support_summary(perm: &[usize]) -> (usize, bool, BTreeMap<Vec<u32>, i64>) {
    let expansion = schubert_to_key::<i64>(perm);
    let positive = expansion.terms().values().all(|c| *c > 0);
    let terms = expansion
        .terms()
        .iter()
        .map(|(alpha, coeff)| (alpha.parts().to_vec(), *coeff))
        .collect::<BTreeMap<_, _>>();
    (terms.len(), positive, terms)
}

fn key_expansion_by_sparse_triangular_reduction(poly: &MultiPoly<i64>) -> BTreeMap<Vec<u32>, i64> {
    let mut residual = poly.clone();
    let mut expansion = BTreeMap::new();

    while !residual.is_zero() {
        let (alpha, coeff) = residual
            .terms()
            .iter()
            .next()
            .map(|(alpha, coeff)| (alpha.clone(), *coeff))
            .expect("nonzero residual has a leading monomial");
        let key = key_polynomial::<i64>(&alpha);
        assert_eq!(
            key.coefficient(&alpha),
            1,
            "key leading coefficient is not 1 for {alpha:?}"
        );
        residual = residual - key.scale(&coeff);
        *expansion.entry(alpha).or_insert(0) += coeff;
    }

    expansion.retain(|_, coeff| *coeff != 0);
    expansion
}

fn print_one(perm: &[usize], max_dilation: usize, show_terms: bool) {
    println!("u={perm:?}, code={:?}", lehmer_code(perm));
    for dilation in 1..=max_dilation {
        let stretched = stretch_perm(perm, dilation);
        let code = lehmer_code(&stretched);
        let (count, positive, terms) = key_support_summary(&stretched);
        println!(
            "  N={dilation}: N*u={stretched:?}, code={code:?}, key_terms={count}, positive={positive}"
        );
        if show_terms {
            for (alpha, coeff) in terms {
                println!("    {coeff:>2} * kappa_{alpha:?}");
            }
        }
    }
}

fn check_21354_formula(max_dilation: usize) {
    let perm = vec![2, 1, 3, 5, 4];
    for dilation in 1..=max_dilation {
        let stretched = stretch_perm(&perm, dilation);
        let num_vars = stretched.len();
        let schubert = schubert_polynomial::<i64>(&stretched);
        let mut rhs = MultiPoly::zero(num_vars);
        for j in 0..=dilation {
            let mut alpha = vec![0; num_vars];
            alpha[0] = (dilation + j) as u32;
            alpha[3] = (dilation - j) as u32;
            rhs = rhs + key_polynomial::<i64>(&alpha);
        }
        println!(
            "N={dilation}: terms in formula={}, monomials={}, formula_ok={}",
            dilation + 1,
            schubert.terms().len(),
            schubert == rhs
        );
    }
}

fn check_21543_formula(max_dilation: usize) {
    let perm = vec![2, 1, 5, 4, 3];
    for dilation in 1..=max_dilation {
        let stretched = stretch_perm(&perm, dilation);
        let num_vars = stretched.len();
        let schubert = schubert_polynomial::<i64>(&stretched);
        let n = dilation as u32;
        let mut rhs = MultiPoly::zero(num_vars);
        let mut alphas = Vec::new();

        let mut alpha = vec![0; num_vars];
        alpha[0] = n;
        alpha[2] = 2 * n;
        alpha[3] = n;
        alphas.push(alpha);

        let mut alpha = vec![0; num_vars];
        alpha[0] = n + 1;
        alpha[2] = 2 * n;
        alpha[3] = n - 1;
        alphas.push(alpha);

        let mut alpha = vec![0; num_vars];
        alpha[0] = 2 * n + 1;
        alpha[2] = n - 1;
        alpha[3] = n;
        alphas.push(alpha);

        for alpha in &alphas {
            rhs = rhs + key_polynomial::<i64>(alpha);
        }
        println!(
            "N={dilation}: alphas={:?}, schubert_terms={}, formula_ok={}",
            alphas
                .iter()
                .map(|alpha| trim_u32(alpha))
                .collect::<Vec<_>>(),
            schubert.terms().len(),
            schubert == rhs
        );
        assert_eq!(schubert, rhs);
    }
}

fn check_s5_hard_key_formulas(max_dilation: usize) {
    let cases = [
        ("24153", vec![2, 4, 1, 5, 3], 1usize),
        ("25143", vec![2, 5, 1, 4, 3], 1usize),
        ("32154", vec![3, 2, 1, 5, 4], 2usize),
        ("42153", vec![4, 2, 1, 5, 3], 1usize),
        ("52143", vec![5, 2, 1, 4, 3], 1usize),
        ("32514", vec![3, 2, 5, 1, 4], 1usize),
        ("32541", vec![3, 2, 5, 4, 1], 1usize),
        ("31542", vec![3, 1, 5, 4, 2], 1usize),
        ("21543", vec![2, 1, 5, 4, 3], 1usize),
    ];

    for (name, perm, first_dilation) in cases {
        println!("{name}: u={perm:?}, code={:?}", lehmer_code(&perm));
        for dilation in first_dilation..=max_dilation {
            let stretched = stretch_perm(&perm, dilation);
            let schubert = schubert_polynomial::<i64>(&stretched);
            let alphas = s5_hard_key_alphas(name, dilation, stretched.len());
            let mut rhs = MultiPoly::zero(stretched.len());
            for alpha in &alphas {
                rhs = rhs + key_polynomial::<i64>(alpha);
            }
            println!(
                "  N={dilation}: terms={}, alphas={:?}, formula_ok={}",
                alphas.len(),
                alphas
                    .iter()
                    .map(|alpha| trim_u32(alpha))
                    .collect::<Vec<_>>(),
                schubert == rhs
            );
            assert_eq!(schubert, rhs);
        }
    }
}

#[derive(Clone, Debug)]
struct AffineKeyTerm {
    slope: Vec<i32>,
    intercept: Vec<i32>,
    coeff: i64,
}

#[derive(Clone, Debug)]
struct AffineKeyExpansion {
    first_dilation: usize,
    terms: Vec<AffineKeyTerm>,
}

#[derive(Clone, Debug)]
struct RankedKeyExample {
    perm: Vec<usize>,
    code: Vec<usize>,
    expansion: AffineKeyExpansion,
}

fn compare_yamanouchi_to_key(perm: &[usize], dilation: usize, max_diagrams: usize) {
    let stretched = stretch_perm(perm, dilation);
    let initial = rothe_diagram(&stretched);
    let diagrams = match kohnert_diagrams(&initial, max_diagrams) {
        Ok(diagrams) => diagrams,
        Err(err) => {
            println!("{err}");
            return;
        }
    };

    let mut yam_weights = BTreeMap::<Vec<u32>, usize>::new();
    let mut yam_diagrams = Vec::new();
    for diagram in &diagrams {
        if is_yamanouchi(&initial, diagram) {
            let weight = diagram_weight(diagram);
            *yam_weights.entry(weight.clone()).or_insert(0) += 1;
            yam_diagrams.push((weight, diagram.clone()));
        }
    }

    let key_terms = key_terms_for_stretch(perm, dilation)
        .into_iter()
        .collect::<BTreeMap<_, _>>();
    println!(
        "u={perm:?}, N={dilation}, N*u={stretched:?}, cells={}, KD={}, Yam={}",
        initial.len(),
        diagrams.len(),
        yam_diagrams.len()
    );
    println!("Yamanouchi weights:");
    for (weight, count) in &yam_weights {
        println!("  {count:>2} * {weight:?}");
    }
    println!("key expansion weights:");
    for (weight, coeff) in &key_terms {
        println!("  {coeff:>2} * {weight:?}");
    }
    let key_as_counts = key_terms
        .iter()
        .filter_map(|(weight, coeff)| usize::try_from(*coeff).ok().map(|c| (weight.clone(), c)))
        .collect::<BTreeMap<_, _>>();
    println!("matches_key_support={}", yam_weights == key_as_counts);
    if yam_diagrams.len() <= 12 {
        println!("Yamanouchi diagrams:");
        for (weight, diagram) in yam_diagrams {
            println!("  wt={weight:?}: {}", format_diagram(&diagram));
        }
    }
}

fn s5_hard_key_alphas(name: &str, dilation: usize, num_vars: usize) -> Vec<Vec<u32>> {
    let n = dilation as u32;
    match name {
        "24153" => vec![
            composition(num_vars, &[(1, n), (2, 2 * n), (4, n)]),
            composition(num_vars, &[(1, n + 1), (2, 2 * n), (4, n - 1)]),
        ],
        "25143" => vec![
            composition(num_vars, &[(1, n), (2, 3 * n), (4, n)]),
            composition(num_vars, &[(1, n + 1), (2, 3 * n), (4, n - 1)]),
        ],
        "32154" => vec![
            composition(num_vars, &[(1, 2 * n), (2, n), (4, n)]),
            composition(num_vars, &[(1, 2 * n), (2, n + 1), (4, n - 1)]),
        ],
        "42153" => vec![
            composition(num_vars, &[(1, 3 * n), (2, n), (4, n)]),
            composition(num_vars, &[(1, 3 * n), (2, n + 1), (4, n - 1)]),
        ],
        "52143" => vec![
            composition(num_vars, &[(1, 4 * n), (2, n), (4, n)]),
            composition(num_vars, &[(1, 4 * n), (2, n + 1), (4, n - 1)]),
        ],
        "32514" => vec![
            composition(num_vars, &[(1, 2 * n), (2, n), (3, 2 * n)]),
            composition(num_vars, &[(1, 2 * n + 1), (2, n), (3, 2 * n - 1)]),
        ],
        "32541" => vec![
            composition(num_vars, &[(1, 2 * n), (2, n), (3, 2 * n), (4, n)]),
            composition(num_vars, &[(1, 2 * n + 1), (2, n), (3, 2 * n - 1), (4, n)]),
        ],
        "31542" => vec![
            composition(num_vars, &[(1, 2 * n), (3, 2 * n), (4, n)]),
            composition(num_vars, &[(1, 2 * n + 1), (3, 2 * n - 1), (4, n)]),
        ],
        "21543" => vec![
            composition(num_vars, &[(1, n), (3, 2 * n), (4, n)]),
            composition(num_vars, &[(1, n + 1), (3, 2 * n), (4, n - 1)]),
            composition(num_vars, &[(1, 2 * n + 1), (3, n - 1), (4, n)]),
        ],
        _ => panic!("unknown hard S5 case {name}"),
    }
}

fn key_terms_for_stretch(perm: &[usize], dilation: usize) -> Vec<(Vec<u32>, i64)> {
    let stretched = stretch_perm(perm, dilation);
    let schubert = schubert_polynomial::<i64>(&stretched);
    key_expansion_by_sparse_triangular_reduction(&schubert)
        .into_iter()
        .map(|(alpha, coeff)| (trim_u32(&alpha), coeff))
        .collect()
}

fn pad_i32(values: &[u32], len: usize) -> Vec<i32> {
    let mut result = vec![0; len];
    for (idx, value) in values.iter().enumerate() {
        result[idx] = *value as i32;
    }
    result
}

fn trim_i32(values: &[i32]) -> Vec<i32> {
    let mut result = values.to_vec();
    while result.last() == Some(&0) {
        result.pop();
    }
    result
}

fn diff_i32(lhs: &[i32], rhs: &[i32]) -> Vec<i32> {
    lhs.iter().zip(rhs.iter()).map(|(a, b)| a - b).collect()
}

fn find_affine_key_expansion(
    expansions: &[Vec<(Vec<u32>, i64)>],
    max_dilation: usize,
) -> Option<AffineKeyExpansion> {
    for first_dilation in 1..max_dilation {
        let slice = &expansions[first_dilation - 1..];
        let term_count = slice.first()?.len();
        if slice.iter().any(|terms| terms.len() != term_count) {
            continue;
        }

        let mut affine_terms = Vec::new();
        let mut works = true;
        for term_idx in 0..term_count {
            let coeff = slice[0][term_idx].1;
            if slice.iter().any(|terms| terms[term_idx].1 != coeff) {
                works = false;
                break;
            }

            let max_len = slice
                .iter()
                .map(|terms| terms[term_idx].0.len())
                .max()
                .unwrap_or(0);
            let padded = slice
                .iter()
                .map(|terms| pad_i32(&terms[term_idx].0, max_len))
                .collect::<Vec<_>>();
            if padded.len() < 2 {
                works = false;
                break;
            }

            let slope = diff_i32(&padded[1], &padded[0]);
            if padded
                .windows(2)
                .any(|pair| diff_i32(&pair[1], &pair[0]) != slope)
            {
                works = false;
                break;
            }

            let intercept = padded[0]
                .iter()
                .zip(slope.iter())
                .map(|(value, delta)| value - (first_dilation as i32) * delta)
                .collect::<Vec<_>>();
            affine_terms.push(AffineKeyTerm {
                slope: trim_i32(&slope),
                intercept: trim_i32(&intercept),
                coeff,
            });
        }

        if works {
            return Some(AffineKeyExpansion {
                first_dilation,
                terms: affine_terms,
            });
        }
    }

    None
}

fn format_vec_i32(values: &[i32]) -> String {
    if values.is_empty() {
        return "()".to_string();
    }
    let entries = values
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(",");
    format!("({entries})")
}

fn format_affine_expansion(expansion: &AffineKeyExpansion) -> String {
    expansion
        .terms
        .iter()
        .map(|term| {
            format!(
                "{}*kappa_{{N{}+{}}}",
                term.coeff,
                format_vec_i32(&term.slope),
                format_vec_i32(&term.intercept)
            )
        })
        .collect::<Vec<_>>()
        .join(" + ")
}

fn evaluate_affine(term: &AffineKeyTerm, dilation: usize) -> Vec<i32> {
    let len = term.slope.len().max(term.intercept.len());
    let mut result = Vec::with_capacity(len);
    for idx in 0..len {
        let slope = term.slope.get(idx).copied().unwrap_or(0);
        let intercept = term.intercept.get(idx).copied().unwrap_or(0);
        result.push((dilation as i32) * slope + intercept);
    }
    trim_i32(&result)
}

fn shape_and_rank_perm(alpha: &[i32]) -> (Vec<i32>, Vec<usize>) {
    let mut indexed = alpha
        .iter()
        .copied()
        .enumerate()
        .collect::<Vec<(usize, i32)>>();
    indexed.sort_by(|(left_idx, left), (right_idx, right)| {
        right.cmp(left).then_with(|| left_idx.cmp(right_idx))
    });

    let shape = trim_i32(&indexed.iter().map(|(_, value)| *value).collect::<Vec<_>>());
    let mut rank_perm = vec![0; alpha.len()];
    for (rank, (idx, _)) in indexed.iter().enumerate() {
        rank_perm[*idx] = rank + 1;
    }
    (shape, rank_perm)
}

fn format_rank_perm(values: &[usize]) -> String {
    if values.is_empty() {
        return "()".to_string();
    }
    let entries = values
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(",");
    format!("({entries})")
}

fn profile_key_stretch(n: usize, max_dilation: usize, show_formulas: bool) {
    let mut distribution = BTreeMap::<(usize, usize), usize>::new();
    let mut failures = Vec::new();
    let mut interesting = Vec::new();

    for perm in all_perms(n) {
        let expansions = (1..=max_dilation)
            .map(|dilation| key_terms_for_stretch(&perm, dilation))
            .collect::<Vec<_>>();
        let Some(affine) = find_affine_key_expansion(&expansions, max_dilation) else {
            failures.push(perm);
            continue;
        };

        let term_count = affine.terms.len();
        *distribution
            .entry((affine.first_dilation, term_count))
            .or_insert(0) += 1;

        if show_formulas && (term_count > 1 || affine.first_dilation > 1) {
            interesting.push((perm, affine));
        }
    }

    println!("S_{n}, checked N=1..{max_dilation}");
    println!("distribution by (first stable N, number of key terms):");
    for ((first_dilation, term_count), count) in distribution {
        println!("  ({first_dilation}, {term_count}) -> {count}");
    }

    if failures.is_empty() {
        println!("all permutations fit an affine key-support model in the checked range");
    } else {
        println!("failures:");
        for perm in failures {
            println!("  {perm:?}, code={:?}", lehmer_code(&perm));
        }
    }

    if show_formulas {
        println!("nontrivial affine formulas:");
        for (perm, affine) in interesting {
            println!(
                "  {:?}, code={:?}, first_N={}: {}",
                perm,
                lehmer_code(&perm),
                affine.first_dilation,
                format_affine_expansion(&affine)
            );
        }
    }
}

fn top_key_stretch_examples(n: usize, max_dilation: usize, count: usize) {
    let mut examples = Vec::new();

    for perm in all_perms(n) {
        let expansions = (1..=max_dilation)
            .map(|dilation| key_terms_for_stretch(&perm, dilation))
            .collect::<Vec<_>>();
        let Some(expansion) = find_affine_key_expansion(&expansions, max_dilation) else {
            println!(
                "no affine model found in checked range for {:?}, code={:?}",
                perm,
                lehmer_code(&perm)
            );
            continue;
        };
        examples.push(RankedKeyExample {
            code: lehmer_code(&perm),
            perm,
            expansion,
        });
    }

    examples.sort_by(|left, right| {
        right
            .expansion
            .terms
            .len()
            .cmp(&left.expansion.terms.len())
            .then_with(|| {
                right
                    .expansion
                    .first_dilation
                    .cmp(&left.expansion.first_dilation)
            })
            .then_with(|| {
                right
                    .code
                    .iter()
                    .sum::<usize>()
                    .cmp(&left.code.iter().sum())
            })
            .then_with(|| left.perm.cmp(&right.perm))
    });

    println!("Top {count} examples in S_{n}, checked N=1..{max_dilation}");
    for (idx, example) in examples.into_iter().take(count).enumerate() {
        let max_coeff = example
            .expansion
            .terms
            .iter()
            .map(|term| term.coeff.abs())
            .max()
            .unwrap_or(0);
        println!(
            "{}. u={:?}, code={:?}, first_N={}, key_terms={}, max_coeff={}",
            idx + 1,
            example.perm,
            example.code,
            example.expansion.first_dilation,
            example.expansion.terms.len(),
            max_coeff
        );
        for term in &example.expansion.terms {
            let alpha = evaluate_affine(term, example.expansion.first_dilation);
            let (shape, rank_perm) = shape_and_rank_perm(&alpha);
            println!(
                "   {} * kappa_{{N{}+{}}}; at N0 alpha={}, shape={}, rank={}",
                term.coeff,
                format_vec_i32(&term.slope),
                format_vec_i32(&term.intercept),
                format_vec_i32(&alpha),
                format_vec_i32(&shape),
                format_rank_perm(&rank_perm)
            );
        }
    }
}

fn composition(num_vars: usize, entries: &[(usize, u32)]) -> Vec<u32> {
    let mut alpha = vec![0; num_vars];
    for &(idx, value) in entries {
        alpha[idx - 1] = value;
    }
    alpha
}

fn print_one_triangular(perm: &[usize], max_dilation: usize) {
    println!("u={perm:?}, code={:?}", lehmer_code(perm));
    for dilation in 1..=max_dilation {
        let stretched = stretch_perm(perm, dilation);
        let code = lehmer_code(&stretched);
        let schubert = schubert_polynomial::<i64>(&stretched);
        let terms = key_expansion_by_sparse_triangular_reduction(&schubert);
        println!(
            "  N={dilation}: N*u={stretched:?}, code={code:?}, key_terms={}",
            terms.len()
        );
        for (alpha, coeff) in terms {
            println!("    {coeff:>2} * kappa_{alpha:?}");
        }
    }
}

fn all_perms(n: usize) -> Vec<Vec<usize>> {
    let mut perm: Vec<usize> = (1..=n).collect();
    let mut result = Vec::new();
    loop {
        result.push(perm.clone());
        if !next_permutation(&mut perm) {
            break;
        }
    }
    result
}

fn trim_u32(values: &[u32]) -> Vec<u32> {
    let mut result = values.to_vec();
    while result.last() == Some(&0) {
        result.pop();
    }
    result
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("--check-21354") {
        let max_dilation = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(5);
        check_21354_formula(max_dilation);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--check-21543") {
        let max_dilation = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(6);
        check_21543_formula(max_dilation);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--check-s5-hard") {
        let max_dilation = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(8);
        check_s5_hard_key_formulas(max_dilation);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--profile") {
        let n = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(5);
        let max_dilation = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(5);
        let show_formulas = args.get(4).map(String::as_str) == Some("--show");
        profile_key_stretch(n, max_dilation, show_formulas);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--top") {
        let n = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(6);
        let max_dilation = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(5);
        let count = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(10);
        top_key_stretch_examples(n, max_dilation, count);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--yam") {
        let Some(perm_arg) = args.get(2) else {
            panic!("usage: --yam <comma-separated permutation> [dilation] [max diagrams]");
        };
        let perm = perm_arg
            .split(',')
            .map(|s| {
                s.parse::<usize>()
                    .expect("permutation entries are integers")
            })
            .collect::<Vec<_>>();
        let dilation = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1);
        let max_diagrams = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(200_000);
        compare_yamanouchi_to_key(&perm, dilation, max_diagrams);
        return;
    }

    if args.get(1).map(String::as_str) == Some("--triangular") {
        let Some(perm_arg) = args.get(2) else {
            panic!("usage: --triangular <comma-separated permutation> [max dilation]");
        };
        let perm = perm_arg
            .split(',')
            .map(|s| {
                s.parse::<usize>()
                    .expect("permutation entries are integers")
            })
            .collect::<Vec<_>>();
        let max_dilation = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(3);
        print_one_triangular(&perm, max_dilation);
        return;
    }

    if args.len() > 1 {
        let perm = args[1]
            .split(',')
            .map(|s| {
                s.parse::<usize>()
                    .expect("permutation entries are integers")
            })
            .collect::<Vec<_>>();
        let max_dilation = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3);
        print_one(&perm, max_dilation, true);
        return;
    }

    for n in 1..=4 {
        println!("S_{n}");
        for perm in all_perms(n) {
            let max_dilation = if n <= 3 { 2 } else { 1 };
            print_one(&perm, max_dilation, false);
        }
    }
}
