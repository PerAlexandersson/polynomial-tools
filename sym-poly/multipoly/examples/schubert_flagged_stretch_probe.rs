use std::collections::BTreeMap;

use sym_poly_multipoly::{schubert_polynomial, MultiPoly};

type Expansion = BTreeMap<(Vec<u32>, Vec<usize>), i64>;

#[derive(Clone, Copy)]
enum ThreeTermCorrection {
    AddRightStripSubtractRectangle,
    SubtractRightStrip,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("--scan") {
        let n = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4);
        let max_dilation = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(5);
        scan_permutations(n, max_dilation);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--profile") {
        let n = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(5);
        let max_dilation = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(6);
        profile_permutations(n, max_dilation);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--rank") {
        let n = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4);
        let max_dilation = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(5);
        rank_scan(n, max_dilation);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--word") {
        let Some(perm_arg) = args.get(2) else {
            panic!("usage: --word <comma-separated permutation> [max dilation]");
        };
        let perm = parse_perm(perm_arg);
        let max_dilation = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(4);
        print_downward_words(&perm, max_dilation);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--transition-step") {
        let Some(perm_arg) = args.get(2) else {
            panic!("usage: --transition-step <comma-separated permutation> [dilation]");
        };
        let perm = parse_perm(perm_arg);
        let dilation = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(2);
        print_transition_step(&perm, dilation);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--collapse-21534") {
        let max_dilation = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(6);
        check_21534_collapse(max_dilation);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--collapse-21543") {
        let max_dilation = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(6);
        check_21543_collapse(max_dilation);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--collapse-2143") {
        let max_dilation = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(8);
        check_2143_collapse(max_dilation);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--collapse-s5-hard") {
        let max_dilation = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(8);
        check_s5_hard_collapses(max_dilation);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--bjs") {
        let n = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(5);
        let max_dilation = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(4);
        check_bjs_321_avoiding(n, max_dilation);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--s5-blocks") {
        let max_dilation = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(5);
        let show_diff = args.get(3).is_some_and(|s| s == "diff");
        check_s5_three_term_blocks(max_dilation, show_diff);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--search-skew") {
        let Some(perm_arg) = args.get(2) else {
            panic!("usage: --search-skew <comma-separated permutation> [dilation] [max rows]");
        };
        let perm = parse_perm(perm_arg);
        let dilation = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(2);
        let max_rows = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(5);
        search_shifted_skew(&perm, dilation, max_rows);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--search-product-skew") {
        let Some(perm_arg) = args.get(2) else {
            panic!(
                "usage: --search-product-skew <comma-separated permutation> [dilation] [max rows] [max flag]"
            );
        };
        let perm = parse_perm(perm_arg);
        let dilation = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(2);
        let max_rows = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(4);
        let max_flag = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(4);
        search_product_of_two_shifted_skews(&perm, dilation, max_rows, max_flag);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--search-21543-blocks") {
        let dilation = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(2);
        let max_rows = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(5);
        search_21543_blocks(dilation, max_rows);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--residual") {
        let Some(perm_arg) = args.get(2) else {
            panic!("usage: --residual <comma-separated permutation> [dilation]");
        };
        let perm = parse_perm(perm_arg);
        let dilation = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(2);
        print_shifted_residual(&perm, dilation);
        return;
    }

    let examples = if args.len() > 1 {
        let perm = parse_perm(&args[1]);
        let max_n = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3);
        vec![(perm, max_n)]
    } else {
        vec![
            (vec![2, 1, 4, 3], 5),
            (vec![2, 1, 3, 5, 4], 5),
            (vec![2, 4, 1, 3], 4),
            (vec![3, 1, 4, 2], 4),
        ]
    };

    for (perm, max_n) in examples {
        println!("u={perm:?}, code={:?}", lehmer_code(&perm));
        for dilation in 1..=max_n {
            let stretched = stretch_perm(&perm, dilation);
            let expansion = transition_flagged_expansion(&stretched);
            let verified = verify_transition_expansion(&stretched, &expansion);
            println!(
                "  N={dilation}: N*u={stretched:?}, code={:?}, leaves={}, verified={verified}",
                lehmer_code(&stretched),
                expansion.len()
            );

            for ((shift, leaf), coeff) in &expansion {
                let leaf_poly = schubert_polynomial::<i64>(leaf);
                let flagged = find_vexillary_flagged_schur(leaf, &leaf_poly);
                let shift_trimmed = trim_u32(shift);
                match flagged {
                    Some((shape, flags)) => println!(
                        "    {coeff:+} * x^{shift_trimmed:?} * s_{shape:?}({flags:?}) leaf={leaf:?}"
                    ),
                    None => println!(
                        "    {coeff:+} * x^{shift_trimmed:?} * S_{leaf:?} (no flagged match found)"
                    ),
                }
            }
        }
    }
}

fn check_21534_collapse(max_dilation: usize) {
    let perm = vec![2, 1, 5, 3, 4];
    println!(
        "finite product collapse for u={perm:?}, code={:?}",
        lehmer_code(&perm)
    );
    for dilation in 1..=max_dilation {
        let stretched = stretch_perm(&perm, dilation);
        let num_vars = stretched.len();
        let schubert = schubert_polynomial::<i64>(&stretched);
        let collapsed = collapse_21534_formula(dilation, num_vars);
        println!(
            "  N={dilation}: schubert_terms={}, collapsed_terms={}, match={}",
            schubert.terms().len(),
            collapsed.terms().len(),
            schubert == collapsed
        );
        assert_eq!(schubert, collapsed);
    }
}

fn check_2143_collapse(max_dilation: usize) {
    let perm = vec![2, 1, 4, 3];
    println!(
        "eventual skew-block collapse for u={perm:?}, code={:?}",
        lehmer_code(&perm)
    );
    for dilation in 1..=max_dilation {
        let stretched = stretch_perm(&perm, dilation);
        let num_vars = stretched.len();
        let schubert = schubert_polynomial::<i64>(&stretched);
        let n = dilation as u32;
        let collapsed = monomial(num_vars, &[(1, n)]) * complete_h(num_vars, n, 2, 3)
            + monomial(num_vars, &[(1, n + 1)]) * complete_h(num_vars, n - 1, 2, 3);
        println!(
            "  N={dilation}: schubert_terms={}, collapsed_terms={}, match={}",
            schubert.terms().len(),
            collapsed.terms().len(),
            schubert == collapsed
        );
        assert_eq!(schubert, collapsed);
    }
}

fn check_21543_collapse(max_dilation: usize) {
    let perm = vec![2, 1, 5, 4, 3];
    println!(
        "two-block residual collapse for u={perm:?}, code={:?}",
        lehmer_code(&perm)
    );
    for dilation in 1..=max_dilation {
        let stretched = stretch_perm(&perm, dilation);
        let schubert = schubert_polynomial::<i64>(&stretched);
        let collapsed = collapse_21543_formula(dilation, stretched.len());
        let residual = divide_by_monomial(
            &schubert,
            &monomial_shift(stretched.len(), &[(1, dilation as u32)]),
        );
        let support_block = block_21543_support(dilation, stretched.len());
        let overlap_block = block_21543_overlap(dilation, stretched.len());
        println!(
            "  N={dilation}: schubert_terms={}, collapsed_terms={}, support_terms={}, overlap_terms={}, match={}",
            schubert.terms().len(),
            collapsed.terms().len(),
            support_block.terms().len(),
            overlap_block.terms().len(),
            schubert == collapsed
        );
        assert_eq!(residual, support_block + overlap_block);
        assert_eq!(schubert, collapsed);
    }
}

fn collapse_21543_formula(dilation: usize, num_vars: usize) -> MultiPoly<i64> {
    let n = dilation as u32;
    monomial(num_vars, &[(1, n)])
        * (block_21543_support(dilation, num_vars) + block_21543_overlap(dilation, num_vars))
}

fn block_21543_support(dilation: usize, num_vars: usize) -> MultiPoly<i64> {
    let n = dilation as u32;
    support_poly(num_vars, 4, 3 * n, |exp| {
        exp[0] <= n + 1
            && exp[1] <= 2 * n
            && exp[2] <= 2 * n
            && exp[3] <= n
            && exp[0] + exp[1] <= 2 * n + 1
            && exp[0] + exp[2] <= 2 * n + 1
    })
}

fn block_21543_overlap(dilation: usize, num_vars: usize) -> MultiPoly<i64> {
    let n = dilation as u32;
    support_poly(num_vars, 4, 3 * n, |exp| {
        1 <= exp[0]
            && exp[0] <= n + 1
            && 1 <= exp[1]
            && exp[0] + exp[1] <= 2 * n
            && exp[0] + exp[2] <= 2 * n
            && exp[3] <= n - 1
    })
}

fn support_poly(
    num_vars: usize,
    active_vars: usize,
    degree: u32,
    predicate: impl Fn(&[u32]) -> bool,
) -> MultiPoly<i64> {
    assert!(active_vars <= num_vars);
    let mut terms = BTreeMap::new();
    let mut exp = vec![0; num_vars];
    support_poly_inner(degree, 0, active_vars, &mut exp, &predicate, &mut terms);
    MultiPoly::from_terms(num_vars, terms)
}

fn support_poly_inner(
    remaining: u32,
    var: usize,
    active_vars: usize,
    exp: &mut [u32],
    predicate: &impl Fn(&[u32]) -> bool,
    terms: &mut BTreeMap<Vec<u32>, i64>,
) {
    if var + 1 == active_vars {
        exp[var] = remaining;
        if predicate(exp) {
            terms.insert(exp.to_vec(), 1);
        }
        exp[var] = 0;
        return;
    }

    for value in 0..=remaining {
        exp[var] = value;
        support_poly_inner(
            remaining - value,
            var + 1,
            active_vars,
            exp,
            predicate,
            terms,
        );
    }
    exp[var] = 0;
}

fn check_s5_hard_collapses(max_dilation: usize) {
    let cases: [(Vec<usize>, fn(usize, usize) -> MultiPoly<i64>, &str); 3] = [
        (vec![3, 2, 5, 1, 4], collapse_32514_formula, "32514 block"),
        (
            vec![3, 2, 5, 4, 1],
            collapse_32541_formula,
            "x4^N times the 32514 block",
        ),
        (
            vec![3, 1, 5, 4, 2],
            collapse_31542_formula,
            "summed 31542 chain",
        ),
    ];

    for (perm, formula, label) in cases {
        println!("{label} for u={perm:?}, code={:?}", lehmer_code(&perm));
        for dilation in 1..=max_dilation {
            let stretched = stretch_perm(&perm, dilation);
            let schubert = schubert_polynomial::<i64>(&stretched);
            let collapsed = formula(dilation, stretched.len());
            println!(
                "  N={dilation}: schubert_terms={}, collapsed_terms={}, match={}",
                schubert.terms().len(),
                collapsed.terms().len(),
                schubert == collapsed
            );
            assert_eq!(schubert, collapsed);
        }
    }
}

fn collapse_32514_formula(dilation: usize, num_vars: usize) -> MultiPoly<i64> {
    let n = dilation as u32;
    let base = monomial(num_vars, &[(1, 2 * n), (2, n), (3, n)]);
    let residual = complete_h(num_vars, n, 2, 3)
        + monomial(num_vars, &[(1, 1)]) * complete_h(num_vars, n - 1, 2, 3);
    base * residual
}

fn collapse_32541_formula(dilation: usize, num_vars: usize) -> MultiPoly<i64> {
    let n = dilation as u32;
    monomial(num_vars, &[(4, n)]) * collapse_32514_formula(dilation, num_vars)
}

fn collapse_31542_formula(dilation: usize, num_vars: usize) -> MultiPoly<i64> {
    let n = dilation as u32;
    let block_sum = complete_h(num_vars, 3 * n - 1, 2, 4)
        - monomial(num_vars, &[(4, n + 1)]) * complete_h(num_vars, 2 * n - 2, 2, 4)
        - monomial(num_vars, &[(2, 2 * n)]) * complete_h(num_vars, n - 1, 2, 4)
        - monomial(num_vars, &[(3, 2 * n)]) * complete_h(num_vars, n - 1, 2, 4);
    let residual = monomial(num_vars, &[(2, 2 * n)]) * complete_h(num_vars, n, 3, 4)
        + (monomial(num_vars, &[(1, 1)]) + monomial(num_vars, &[(3, 1)])) * block_sum;
    monomial(num_vars, &[(1, 2 * n)]) * residual
}

fn check_bjs_321_avoiding(n: usize, max_dilation: usize) {
    println!("BJS flagged skew check for 321-avoiding permutations in S_{n}");
    for perm in all_perms(n) {
        if !avoids_pattern(&perm, &[3, 2, 1]) {
            continue;
        }

        let mut skipped = Vec::new();
        for dilation in 2..=max_dilation {
            let stretched = stretch_perm(&perm, dilation);
            if !avoids_pattern(&stretched, &[3, 2, 1]) {
                skipped.push(dilation);
                continue;
            }
            let Some((lambda, mu, flags)) = bjs_shape_and_flags(&stretched) else {
                skipped.push(dilation);
                continue;
            };
            let schubert = schubert_polynomial::<i64>(&stretched);
            let skew = flagged_skew_schur(&lambda, &mu, &flags, stretched.len());
            if schubert != skew {
                println!(
                    "  candidate formula mismatch for u={perm:?}, N={dilation}, lambda={lambda:?}, mu={mu:?}, flags={flags:?}"
                );
            }
        }

        let stretched = stretch_perm(&perm, max_dilation);
        if !skipped.is_empty() {
            println!(
                "  u={perm:?}, code={:?}, stretched not 321-avoiding for N={skipped:?}",
                lehmer_code(&perm)
            );
        } else if let Some((lambda, mu, flags)) = bjs_shape_and_flags(&stretched) {
            println!(
                "  u={perm:?}, code={:?}, at N={max_dilation}: lambda={lambda:?}, mu={mu:?}, flags={flags:?}",
                lehmer_code(&perm)
            );
        }
    }
}

fn check_s5_three_term_blocks(max_dilation: usize, show_diff: bool) {
    let cases = [
        (
            vec![2, 4, 1, 5, 3],
            2usize,
            vec![2, 2, 4],
            ThreeTermCorrection::AddRightStripSubtractRectangle,
        ),
        (
            vec![2, 5, 1, 4, 3],
            3usize,
            vec![2, 2, 4],
            ThreeTermCorrection::AddRightStripSubtractRectangle,
        ),
        (
            vec![3, 2, 1, 5, 4],
            2usize,
            vec![1, 2, 4],
            ThreeTermCorrection::SubtractRightStrip,
        ),
        (
            vec![4, 2, 1, 5, 3],
            3usize,
            vec![1, 2, 4],
            ThreeTermCorrection::SubtractRightStrip,
        ),
        (
            vec![5, 2, 1, 4, 3],
            4usize,
            vec![1, 2, 4],
            ThreeTermCorrection::SubtractRightStrip,
        ),
    ];

    for (perm, first_slope, flags, correction_kind) in cases {
        println!(
            "three-term skew block candidate for u={perm:?}, code={:?}",
            lehmer_code(&perm)
        );
        for dilation in 2..=max_dilation {
            let n = dilation;
            let stretched = stretch_perm(&perm, dilation);
            let schubert = schubert_polynomial::<i64>(&stretched);
            let candidate = flagged_skew_schur(
                &[first_slope * n, n + 1, n],
                &[0, 1, 0],
                &flags,
                stretched.len(),
            );
            let corrected = candidate.clone()
                + s5_three_term_correction(first_slope, dilation, stretched.len(), correction_kind);
            println!(
                "  N={dilation}: schubert_terms={}, candidate_terms={}, match={}, corrected_terms={}, corrected_match={}",
                schubert.terms().len(),
                candidate.terms().len(),
                schubert == candidate,
                corrected.terms().len(),
                schubert == corrected
            );
            if show_diff && schubert != candidate {
                let diff = schubert - candidate;
                print_signed_support("    schubert-candidate", &diff, 120);
            }
        }
    }
}

fn s5_three_term_correction(
    first_slope: usize,
    dilation: usize,
    num_vars: usize,
    correction_kind: ThreeTermCorrection,
) -> MultiPoly<i64> {
    let n = dilation as u32;
    match correction_kind {
        ThreeTermCorrection::AddRightStripSubtractRectangle => {
            let right_strip = monomial(num_vars, &[(1, n), (2, first_slope as u32 * n)])
                * complete_h(num_vars, n, 3, 4);
            let rectangle = monomial(num_vars, &[(1, n + 2), (2, n)])
                * complete_h(num_vars, (first_slope as u32 - 1) * n - 1, 1, 2)
                * complete_h(num_vars, n - 1, 3, 4);
            right_strip - rectangle
        }
        ThreeTermCorrection::SubtractRightStrip => {
            let right_strip = monomial(num_vars, &[(1, first_slope as u32 * n + 1), (2, n)])
                * complete_h(num_vars, n - 1, 3, 4);
            -right_strip
        }
    }
}

fn print_signed_support(label: &str, poly: &MultiPoly<i64>, max_terms: usize) {
    let positive = poly.terms().values().filter(|&&coeff| coeff > 0).count();
    let negative = poly.terms().values().filter(|&&coeff| coeff < 0).count();
    let mut by_coeff = BTreeMap::new();
    for &coeff in poly.terms().values() {
        *by_coeff.entry(coeff).or_insert(0usize) += 1;
    }
    println!(
        "{label}: terms={}, positive={}, negative={}, coeffs={by_coeff:?}",
        poly.terms().len(),
        positive,
        negative
    );
    for (idx, (exp, coeff)) in poly.terms().iter().enumerate() {
        if idx == max_terms {
            println!("      ...");
            break;
        }
        println!("      {coeff:+} x^{:?}", trim_u32(exp));
    }
}

fn search_shifted_skew(perm: &[usize], dilation: usize, max_rows: usize) {
    let stretched = stretch_perm(perm, dilation);
    let poly = schubert_polynomial::<i64>(&stretched);
    let shift = componentwise_min_exponent(&poly);
    let residual = divide_by_monomial(&poly, &shift);
    let degree = residual.total_degree().unwrap_or(0) as usize;
    println!(
        "search shifted skew for u={perm:?}, N={dilation}, stretched={stretched:?}, shift={:?}, residual_degree={degree}, residual_terms={}",
        trim_u32(&shift),
        residual.terms().len()
    );

    let mut matches = 0usize;
    for rows in 1..=max_rows.min(degree.max(1)) {
        for lambda in partitions_with_len_at_most(degree, rows) {
            let mut padded_lambda = lambda.clone();
            padded_lambda.resize(rows, 0);
            for mu in subpartitions(&padded_lambda) {
                if skew_size(&padded_lambda, &mu) != degree {
                    continue;
                }
                for flags in weak_flags(rows, stretched.len()) {
                    let skew = flagged_skew_schur(&padded_lambda, &mu, &flags, stretched.len());
                    if skew == residual {
                        println!("  match lambda={padded_lambda:?}, mu={mu:?}, flags={flags:?}");
                        matches += 1;
                        if matches >= 20 {
                            return;
                        }
                    }
                }
            }
        }
    }

    if matches == 0 {
        println!("  no single flagged skew match found");
    }
}

fn print_shifted_residual(perm: &[usize], dilation: usize) {
    let stretched = stretch_perm(perm, dilation);
    let poly = schubert_polynomial::<i64>(&stretched);
    let shift = componentwise_min_exponent(&poly);
    let residual = divide_by_monomial(&poly, &shift);
    println!(
        "shifted residual for u={perm:?}, N={dilation}, stretched={stretched:?}, shift={:?}",
        trim_u32(&shift)
    );
    print_signed_support("  residual", &residual, 200);
}

fn search_product_of_two_shifted_skews(
    perm: &[usize],
    dilation: usize,
    max_rows: usize,
    max_flag: usize,
) {
    let stretched = stretch_perm(perm, dilation);
    let poly = schubert_polynomial::<i64>(&stretched);
    let shift = componentwise_min_exponent(&poly);
    let residual = divide_by_monomial(&poly, &shift);
    let degree = residual.total_degree().unwrap_or(0) as usize;
    println!(
        "search product of two shifted skews for u={perm:?}, N={dilation}, shift={:?}, residual_degree={degree}, residual_terms={}",
        trim_u32(&shift),
        residual.terms().len()
    );
    search_product_of_two_skews_for_poly(&residual, degree, max_rows, max_flag);
}

fn search_product_of_two_skews_for_poly(
    poly: &MultiPoly<i64>,
    degree: usize,
    max_rows: usize,
    max_flag: usize,
) {
    let candidates = (0..=degree)
        .map(|d| skew_candidates(d, max_rows, max_flag, poly.num_vars()))
        .collect::<Vec<_>>();
    for (d, values) in candidates.iter().enumerate() {
        println!("  degree {d}: unique candidates={}", values.len());
    }

    let mut matches = 0usize;
    for d in 0..=degree {
        for (left_poly, left_desc) in &candidates[d] {
            for (right_poly, right_desc) in &candidates[degree - d] {
                if left_poly.clone() * right_poly.clone() == *poly {
                    println!("  match: {left_desc} * {right_desc}");
                    matches += 1;
                    if matches >= 20 {
                        return;
                    }
                }
            }
        }
    }

    if matches == 0 {
        println!("  no product of two flagged skew matches found");
    }
}

fn search_21543_blocks(dilation: usize, max_rows: usize) {
    let num_vars = stretch_perm(&[2, 1, 5, 4, 3], dilation).len();
    let blocks = [
        (
            "support",
            block_21543_support(dilation, num_vars),
            3 * dilation,
        ),
        (
            "overlap",
            block_21543_overlap(dilation, num_vars),
            3 * dilation,
        ),
    ];

    for (name, block, degree) in blocks {
        println!(
            "search {name} block for 21543, N={dilation}: degree={degree}, terms={}",
            block.terms().len()
        );
        println!("  single flagged skew search:");
        search_skew_for_poly(&block, degree, max_rows);
        println!("  product of two flagged skew search:");
        search_product_of_two_skews_for_poly(&block, degree, max_rows, 4);
    }
}

fn search_skew_for_poly(poly: &MultiPoly<i64>, degree: usize, max_rows: usize) {
    let mut matches = 0usize;
    for rows in 1..=max_rows.min(degree.max(1)) {
        for lambda in partitions_with_len_at_most(degree, rows) {
            let mut padded_lambda = lambda.clone();
            padded_lambda.resize(rows, 0);
            for mu in subpartitions(&padded_lambda) {
                if skew_size(&padded_lambda, &mu) != degree {
                    continue;
                }
                for flags in weak_flags(rows, poly.num_vars()) {
                    let skew = flagged_skew_schur(&padded_lambda, &mu, &flags, poly.num_vars());
                    if skew == *poly {
                        println!("  match lambda={padded_lambda:?}, mu={mu:?}, flags={flags:?}");
                        matches += 1;
                        if matches >= 20 {
                            return;
                        }
                    }
                }
            }
        }
    }

    if matches == 0 {
        println!("  no single flagged skew match found");
    }
}

fn skew_candidates(
    degree: usize,
    max_rows: usize,
    max_flag: usize,
    num_vars: usize,
) -> Vec<(MultiPoly<i64>, String)> {
    if degree == 0 {
        return vec![(MultiPoly::constant(num_vars, 1), "1".to_string())];
    }

    let mut by_signature: BTreeMap<Vec<(Vec<u32>, i64)>, (MultiPoly<i64>, String)> =
        BTreeMap::new();
    for rows in 1..=max_rows.min(degree) {
        for lambda in partitions_with_len_at_most(degree, rows) {
            let mut padded_lambda = lambda.clone();
            padded_lambda.resize(rows, 0);
            for mu in subpartitions(&padded_lambda) {
                if skew_size(&padded_lambda, &mu) != degree {
                    continue;
                }
                for flags in weak_flags(rows, max_flag) {
                    let skew = flagged_skew_schur(&padded_lambda, &mu, &flags, num_vars);
                    if skew.is_zero() {
                        continue;
                    }
                    let desc = format!("s_{padded_lambda:?}/{mu:?}({flags:?})");
                    by_signature
                        .entry(poly_signature(&skew))
                        .or_insert((skew, desc));
                }
            }
        }
    }
    by_signature.into_values().collect()
}

fn componentwise_min_exponent(poly: &MultiPoly<i64>) -> Vec<u32> {
    let mut result = vec![u32::MAX; poly.num_vars()];
    for exp in poly.terms().keys() {
        for (idx, &power) in exp.iter().enumerate() {
            result[idx] = result[idx].min(power);
        }
    }
    for entry in &mut result {
        if *entry == u32::MAX {
            *entry = 0;
        }
    }
    result
}

fn divide_by_monomial(poly: &MultiPoly<i64>, shift: &[u32]) -> MultiPoly<i64> {
    let terms = poly
        .terms()
        .iter()
        .map(|(exp, coeff)| {
            let residual = exp
                .iter()
                .zip(shift)
                .map(|(&value, &offset)| value - offset)
                .collect::<Vec<_>>();
            (residual, *coeff)
        })
        .collect();
    MultiPoly::from_terms(poly.num_vars(), terms)
}

fn partitions_with_len_at_most(total: usize, max_len: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut current = Vec::new();
    partitions_inner(total, total, max_len, &mut current, &mut result);
    result
}

fn partitions_inner(
    remaining: usize,
    max_part: usize,
    max_len: usize,
    current: &mut Vec<usize>,
    result: &mut Vec<Vec<usize>>,
) {
    if remaining == 0 {
        result.push(current.clone());
        return;
    }
    if current.len() == max_len {
        return;
    }
    for part in (1..=remaining.min(max_part)).rev() {
        current.push(part);
        partitions_inner(remaining - part, part, max_len, current, result);
        current.pop();
    }
}

fn subpartitions(lambda: &[usize]) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut current = vec![0; lambda.len()];
    subpartitions_inner(
        lambda,
        0,
        lambda.first().copied().unwrap_or(0),
        &mut current,
        &mut result,
    );
    result
}

fn subpartitions_inner(
    lambda: &[usize],
    row: usize,
    previous: usize,
    current: &mut [usize],
    result: &mut Vec<Vec<usize>>,
) {
    if row == lambda.len() {
        result.push(current.to_vec());
        return;
    }
    for value in (0..=lambda[row].min(previous)).rev() {
        current[row] = value;
        subpartitions_inner(lambda, row + 1, value, current, result);
    }
}

fn skew_size(lambda: &[usize], mu: &[usize]) -> usize {
    lambda.iter().sum::<usize>() - mu.iter().sum::<usize>()
}

fn weak_flags(rows: usize, max_flag: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut current = Vec::new();
    weak_flags_inner(rows, max_flag, 1, &mut current, &mut result);
    result
}

fn weak_flags_inner(
    rows: usize,
    max_flag: usize,
    min_flag: usize,
    current: &mut Vec<usize>,
    result: &mut Vec<Vec<usize>>,
) {
    if current.len() == rows {
        result.push(current.clone());
        return;
    }
    for flag in min_flag..=max_flag {
        current.push(flag);
        weak_flags_inner(rows, max_flag, flag, current, result);
        current.pop();
    }
}

fn bjs_shape_and_flags(perm: &[usize]) -> Option<(Vec<usize>, Vec<usize>, Vec<usize>)> {
    if !avoids_pattern(perm, &[3, 2, 1]) {
        return None;
    }

    let excedances = (1..=perm.len())
        .filter(|&idx| perm[idx - 1] > idx)
        .collect::<Vec<_>>();
    if excedances.is_empty() {
        return Some((vec![], vec![], vec![]));
    }

    let e = excedances
        .iter()
        .map(|&idx| perm[idx - 1] - 1)
        .collect::<Vec<_>>();
    let t = excedances.len();
    let mut lambda = Vec::with_capacity(t);
    let mut mu = Vec::with_capacity(t);
    for row in 0..t {
        lambda.push(e[t - 1 - row].checked_sub(t - row - 1)?);
        mu.push(excedances[t - 1 - row].checked_sub(t - row)?);
    }
    Some((lambda, mu, excedances))
}

fn collapse_21534_formula(dilation: usize, num_vars: usize) -> MultiPoly<i64> {
    let n = dilation as u32;
    let x1n_x2n = monomial(num_vars, &[(1, n), (2, n)]);
    let x1n = monomial(num_vars, &[(1, n)]);
    let x1n_x2 = monomial(num_vars, &[(1, n), (2, 1)]);

    let h_n_13 = complete_h(num_vars, n, 1, 3);
    let h_np1_13 = complete_h(num_vars, n + 1, 1, 3);
    let h_nm1_23 = complete_h(num_vars, n.saturating_sub(1), 2, 3);

    // This is
    // x1^N x2^N h_N(x1,x2,x3)
    //   + x1^N (h_{N+1}(x1,x2,x3)-x2 h_N(x1,x2,x3)) h_{N-1}(x2,x3).
    x1n_x2n * h_n_13.clone() + x1n * h_np1_13 * h_nm1_23.clone() - x1n_x2 * h_n_13 * h_nm1_23
}

fn print_downward_words(perm: &[usize], max_dilation: usize) {
    println!("u={perm:?}, code={:?}", lehmer_code(perm));
    for dilation in 1..=max_dilation {
        let stretched = stretch_perm(perm, dilation);
        let word = top_operator_word(&stretched);
        println!(
            "  N={dilation}: len={}, degree={}, word_len={}, runs={:?}, word={:?}",
            stretched.len(),
            inversion_count(&stretched),
            word.len(),
            descending_runs(&word),
            word
        );
    }
}

fn print_transition_step(perm: &[usize], dilation: usize) {
    let stretched = stretch_perm(perm, dilation);
    let (r, s) = transition_indices(&stretched);
    let mut v = stretched.clone();
    v.swap(r, s);
    let target_len = inversion_count(&stretched);
    println!(
        "u={perm:?}, N={dilation}, stretched={stretched:?}, code={:?}, r={}, s={}, v={v:?}, v_code={:?}",
        lehmer_code(&stretched),
        r + 1,
        s + 1,
        lehmer_code(&v)
    );
    println!("  x_part: x_{} * S_v", r + 1);
    for j in 0..r {
        let mut child = v.clone();
        child.swap(j, r);
        if inversion_count(&child) == target_len {
            println!(
                "  child j={}: perm={child:?}, code={:?}, vexillary={}, avoids321={}",
                j + 1,
                lehmer_code(&child),
                is_vexillary(&child),
                avoids_pattern(&child, &[3, 2, 1])
            );
        }
    }
}

fn descending_runs(word: &[usize]) -> Vec<(usize, usize)> {
    if word.is_empty() {
        return vec![];
    }

    let mut runs = Vec::new();
    let mut start = word[0];
    let mut previous = word[0];
    for &entry in &word[1..] {
        if entry + 1 == previous {
            previous = entry;
        } else {
            runs.push((start, previous));
            start = entry;
            previous = entry;
        }
    }
    runs.push((start, previous));
    runs
}

fn scan_permutations(n: usize, max_dilation: usize) {
    let mut total = 0usize;
    let mut stable_from_one = 0usize;
    let mut stable_from_two = 0usize;
    let mut max_leaves = 0usize;
    let mut changing = Vec::new();

    for perm in all_perms(n) {
        total += 1;
        let counts = (1..=max_dilation)
            .map(|dilation| {
                let stretched = stretch_perm(&perm, dilation);
                let expansion = transition_flagged_expansion(&stretched);
                assert!(
                    verify_transition_expansion(&stretched, &expansion),
                    "transition expansion failed for {perm:?}, N={dilation}"
                );
                expansion.len()
            })
            .collect::<Vec<_>>();

        max_leaves = max_leaves.max(*counts.iter().max().unwrap_or(&0));
        if counts.windows(2).all(|w| w[0] == w[1]) {
            stable_from_one += 1;
        } else if counts.len() >= 2 && counts[1..].windows(2).all(|w| w[0] == w[1]) {
            stable_from_two += 1;
        } else {
            changing.push((perm, counts));
        }
    }

    println!(
        "scan S_{n}, N=1..{max_dilation}: total={total}, stable_from_1={stable_from_one}, stable_from_2={}, still_changing={}, max_leaves={max_leaves}",
        stable_from_one + stable_from_two,
        changing.len()
    );
    for (perm, counts) in changing {
        println!(
            "  changing u={perm:?}, code={:?}, leaf_counts={counts:?}",
            lehmer_code(&perm)
        );
    }
}

fn profile_permutations(n: usize, max_dilation: usize) {
    let mut by_counts: BTreeMap<Vec<usize>, usize> = BTreeMap::new();
    let mut by_profile: BTreeMap<(bool, bool, usize, usize, usize, usize), usize> = BTreeMap::new();
    let mut changing = Vec::new();
    let mut nonvexillary = Vec::new();

    println!(
        "perm,code,leaves,threshold,degree,vexillary,avoids321,descents,essential,code_support"
    );
    for perm in all_perms(n) {
        let code = lehmer_code(&perm);
        let counts = (1..=max_dilation)
            .map(|dilation| transition_flagged_expansion(&stretch_perm(&perm, dilation)).len())
            .collect::<Vec<_>>();
        let (threshold, degree) = eventual_polynomial_profile(&counts).unwrap_or((1, max_dilation));
        let vexillary = is_vexillary(&perm);
        let avoids321 = avoids_pattern(&perm, &[3, 2, 1]);
        let count2143 = pattern_count(&perm, &[2, 1, 4, 3]);
        let count321 = pattern_count(&perm, &[3, 2, 1]);
        let descents = descent_count(&perm);
        let essential = essential_set_size(&perm);
        let code_support = code.iter().filter(|&&entry| entry > 0).count();

        *by_counts.entry(counts.clone()).or_insert(0) += 1;
        *by_profile
            .entry((vexillary, avoids321, descents, essential, threshold, degree))
            .or_insert(0) += 1;

        if threshold > 1 || degree > 0 {
            changing.push((
                perm.clone(),
                code.clone(),
                counts.clone(),
                threshold,
                degree,
            ));
        }
        if !vexillary {
            nonvexillary.push((
                perm.clone(),
                code.clone(),
                counts.clone(),
                threshold,
                degree,
                avoids321,
                count2143,
                count321,
                descents,
                essential,
            ));
        }

        println!(
            "{perm:?},{code:?},{counts:?},{threshold},{degree},{vexillary},{avoids321},{descents},{essential},{code_support}"
        );
    }

    println!("\ncount sequence distribution:");
    for (counts, multiplicity) in by_counts {
        println!("  {counts:?}: {multiplicity}");
    }

    println!(
        "\nprofile distribution (vexillary, avoids321, descents, essential, threshold, degree):"
    );
    for (profile, multiplicity) in by_profile {
        println!("  {profile:?}: {multiplicity}");
    }

    println!("\nnonstable-from-N=1 raw transition leaves:");
    for (perm, code, counts, threshold, degree) in changing {
        println!(
            "  u={perm:?}, code={code:?}, leaves={counts:?}, threshold={threshold}, degree={degree}"
        );
    }

    println!("\nnonvexillary cases:");
    for (
        perm,
        code,
        counts,
        threshold,
        degree,
        avoids321,
        count2143,
        count321,
        descents,
        essential,
    ) in nonvexillary
    {
        println!(
            "  u={perm:?}, code={code:?}, leaves={counts:?}, threshold={threshold}, degree={degree}, avoids321={avoids321}, #2143={count2143}, #321={count321}, descents={descents}, essential={essential}"
        );
    }
}

fn rank_scan(n: usize, max_dilation: usize) {
    println!("rank scan S_{n}; columns are monomials, ranks computed modulo two primes");
    for dilation in 1..=max_dilation {
        let mut schuberts = Vec::new();
        let mut transition_terms: BTreeMap<Vec<(Vec<u32>, i64)>, MultiPoly<i64>> = BTreeMap::new();
        let mut total_transition_leaves = 0usize;

        for perm in all_perms(n) {
            let stretched = stretch_perm(&perm, dilation);
            schuberts.push(schubert_polynomial::<i64>(&stretched));

            let expansion = transition_flagged_expansion(&stretched);
            total_transition_leaves += expansion.len();
            for ((shift, leaf), _coeff) in expansion {
                let shifted =
                    MultiPoly::x_power(stretched.len(), shift) * schubert_polynomial::<i64>(&leaf);
                transition_terms.insert(poly_signature(&shifted), shifted);
            }
        }

        let schubert_rank = modular_rank_pair(&schuberts);
        let transition_polys = transition_terms.values().cloned().collect::<Vec<_>>();
        let transition_rank = modular_rank_pair(&transition_polys);
        println!(
            "  N={dilation}: schuberts={}, schubert_rank={:?}, raw_leaves={}, unique_transition_terms={}, transition_rank={:?}",
            schuberts.len(),
            schubert_rank,
            total_transition_leaves,
            transition_polys.len(),
            transition_rank
        );
    }
}

fn poly_signature(poly: &MultiPoly<i64>) -> Vec<(Vec<u32>, i64)> {
    poly.terms()
        .iter()
        .map(|(exp, coeff)| (trim_u32(exp), *coeff))
        .collect()
}

fn modular_rank_pair(polys: &[MultiPoly<i64>]) -> (usize, usize) {
    (
        modular_rank(polys, 1_000_000_007),
        modular_rank(polys, 1_000_000_009),
    )
}

fn modular_rank(polys: &[MultiPoly<i64>], prime: i64) -> usize {
    let mut monomial_index = BTreeMap::new();
    for poly in polys {
        for exp in poly.terms().keys() {
            let next = monomial_index.len();
            monomial_index.entry(trim_u32(exp)).or_insert(next);
        }
    }

    let width = monomial_index.len();
    let mut matrix = polys
        .iter()
        .map(|poly| {
            let mut row = vec![0i64; width];
            for (exp, coeff) in poly.terms() {
                let col = monomial_index[&trim_u32(exp)];
                row[col] = coeff.rem_euclid(prime);
            }
            row
        })
        .collect::<Vec<_>>();

    let mut rank = 0usize;
    for col in 0..width {
        let Some(pivot) = (rank..matrix.len()).find(|&row| matrix[row][col] != 0) else {
            continue;
        };
        matrix.swap(rank, pivot);
        let inv = mod_inverse(matrix[rank][col], prime);
        for entry in &mut matrix[rank][col..] {
            *entry = (*entry * inv).rem_euclid(prime);
        }

        for row in 0..matrix.len() {
            if row == rank || matrix[row][col] == 0 {
                continue;
            }
            let factor = matrix[row][col];
            for j in col..width {
                matrix[row][j] = (matrix[row][j] - factor * matrix[rank][j]).rem_euclid(prime);
            }
        }
        rank += 1;
        if rank == matrix.len() {
            break;
        }
    }
    rank
}

fn mod_inverse(value: i64, prime: i64) -> i64 {
    let mut base = value.rem_euclid(prime);
    let mut exp = prime - 2;
    let mut result = 1i64;
    while exp > 0 {
        if exp % 2 == 1 {
            result = (result * base).rem_euclid(prime);
        }
        base = (base * base).rem_euclid(prime);
        exp /= 2;
    }
    result
}

fn eventual_polynomial_profile(values: &[usize]) -> Option<(usize, usize)> {
    let mut best = None;
    for start in 0..values.len() {
        if values.len() - start < 3 {
            continue;
        }
        if let Some(degree) = polynomial_degree_from_values(&values[start..]) {
            match best {
                None => best = Some((start + 1, degree)),
                Some((best_start, best_degree)) => {
                    if degree < best_degree || (degree == best_degree && start + 1 < best_start) {
                        best = Some((start + 1, degree));
                    }
                }
            }
        }
    }
    best
}

fn polynomial_degree_from_values(values: &[usize]) -> Option<usize> {
    if values.is_empty() {
        return None;
    }

    let mut differences = values.iter().map(|&value| value as i64).collect::<Vec<_>>();
    for degree in 0..values.len() {
        if differences.windows(2).all(|window| window[0] == window[1]) {
            return Some(degree);
        }
        differences = differences
            .windows(2)
            .map(|window| window[1] - window[0])
            .collect();
    }
    None
}

fn transition_flagged_expansion(perm: &[usize]) -> Expansion {
    let mut memo = BTreeMap::new();
    transition_flagged_expansion_inner(perm, &mut memo)
}

fn transition_flagged_expansion_inner(
    perm: &[usize],
    memo: &mut BTreeMap<Vec<usize>, Expansion>,
) -> Expansion {
    if let Some(result) = memo.get(perm) {
        return result.clone();
    }

    let n = perm.len();
    if is_vexillary(perm) {
        let mut result = BTreeMap::new();
        result.insert((vec![0u32; n], perm.to_vec()), 1);
        memo.insert(perm.to_vec(), result.clone());
        return result;
    }

    let (r, s) = transition_indices(perm);
    let mut v = perm.to_vec();
    v.swap(r, s);
    let target_len = inversion_count(perm);

    let mut result = BTreeMap::new();

    let x_part = transition_flagged_expansion_inner(&v, memo);
    for ((mut shift, leaf), coeff) in x_part {
        shift[r] += 1;
        *result.entry((shift, leaf)).or_insert(0) += coeff;
    }

    for j in 0..r {
        let mut child = v.clone();
        child.swap(j, r);
        if inversion_count(&child) == target_len {
            let child_expansion = transition_flagged_expansion_inner(&child, memo);
            for (term, coeff) in child_expansion {
                *result.entry(term).or_insert(0) += coeff;
            }
        }
    }

    result.retain(|_, coeff| *coeff != 0);
    memo.insert(perm.to_vec(), result.clone());
    result
}

fn transition_indices(perm: &[usize]) -> (usize, usize) {
    let r = (0..perm.len() - 1)
        .rev()
        .find(|&idx| perm[idx] > perm[idx + 1])
        .expect("non-vexillary permutations have a descent");
    let s = (r + 1..perm.len())
        .rev()
        .find(|&idx| perm[idx] < perm[r])
        .expect("a descent gives a transition partner");
    (r, s)
}

fn verify_transition_expansion(perm: &[usize], expansion: &Expansion) -> bool {
    let n = perm.len();
    let mut rhs = MultiPoly::zero(n);
    for ((shift, leaf), coeff) in expansion {
        let shifted = MultiPoly::x_power(n, shift.clone()) * schubert_polynomial::<i64>(leaf);
        rhs = rhs + shifted.scale(coeff);
    }
    rhs == schubert_polynomial::<i64>(perm)
}

fn is_vexillary(perm: &[usize]) -> bool {
    let n = perm.len();
    for a in 0..n {
        for b in a + 1..n {
            for c in b + 1..n {
                for d in c + 1..n {
                    if perm[b] < perm[a] && perm[a] < perm[d] && perm[d] < perm[c] {
                        return false;
                    }
                }
            }
        }
    }
    true
}

fn avoids_pattern(perm: &[usize], pattern: &[usize]) -> bool {
    let k = pattern.len();
    let mut indices = Vec::with_capacity(k);
    !contains_pattern_inner(perm, pattern, 0, &mut indices)
}

fn pattern_count(perm: &[usize], pattern: &[usize]) -> usize {
    let mut indices = Vec::with_capacity(pattern.len());
    pattern_count_inner(perm, pattern, 0, &mut indices)
}

fn pattern_count_inner(
    perm: &[usize],
    pattern: &[usize],
    next_index: usize,
    indices: &mut Vec<usize>,
) -> usize {
    if indices.len() == pattern.len() {
        return usize::from(has_relative_order(perm, pattern, indices));
    }

    let mut result = 0;
    for index in next_index..perm.len() {
        indices.push(index);
        result += pattern_count_inner(perm, pattern, index + 1, indices);
        indices.pop();
    }
    result
}

fn contains_pattern_inner(
    perm: &[usize],
    pattern: &[usize],
    next_index: usize,
    indices: &mut Vec<usize>,
) -> bool {
    if indices.len() == pattern.len() {
        return has_relative_order(perm, pattern, indices);
    }

    for index in next_index..perm.len() {
        indices.push(index);
        if contains_pattern_inner(perm, pattern, index + 1, indices) {
            return true;
        }
        indices.pop();
    }
    false
}

fn has_relative_order(perm: &[usize], pattern: &[usize], indices: &[usize]) -> bool {
    for i in 0..pattern.len() {
        for j in i + 1..pattern.len() {
            if (perm[indices[i]] < perm[indices[j]]) != (pattern[i] < pattern[j]) {
                return false;
            }
        }
    }
    true
}

fn descent_count(perm: &[usize]) -> usize {
    perm.windows(2)
        .filter(|window| window[0] > window[1])
        .count()
}

fn essential_set_size(perm: &[usize]) -> usize {
    let n = perm.len();
    let inverse = inverse_perm(perm);
    let mut rothe = vec![vec![false; n + 1]; n + 1];
    for i in 1..=n {
        for j in 1..perm[i - 1] {
            if inverse[j - 1] > i {
                rothe[i][j] = true;
            }
        }
    }

    let mut result = 0;
    for i in 1..=n {
        for j in 1..=n {
            if rothe[i][j] && !rothe.get(i + 1).is_some_and(|row| row[j]) && !rothe[i][j + 1] {
                result += 1;
            }
        }
    }
    result
}

fn find_vexillary_flagged_schur(
    leaf: &[usize],
    poly: &MultiPoly<i64>,
) -> Option<(Vec<usize>, Vec<usize>)> {
    let num_vars = poly.num_vars();
    let mut shape = lehmer_code(leaf)
        .into_iter()
        .filter(|&part| part > 0)
        .collect::<Vec<_>>();
    shape.sort_by(|a, b| b.cmp(a));

    if shape.is_empty() {
        return (*poly == MultiPoly::constant(num_vars, 1)).then_some((shape, vec![]));
    }

    let mut flags = Vec::with_capacity(shape.len());
    if search_flags(poly, &shape, num_vars, 0, 1, &mut flags) {
        return Some((shape, flags));
    }
    None
}

fn search_flags(
    poly: &MultiPoly<i64>,
    shape: &[usize],
    num_vars: usize,
    row: usize,
    min_flag: usize,
    flags: &mut Vec<usize>,
) -> bool {
    if row == shape.len() {
        return flagged_schur(shape, flags, num_vars) == *poly;
    }

    let lower = min_flag.max(row + 1);
    for flag in lower..=num_vars {
        flags.push(flag);
        if search_flags(poly, shape, num_vars, row + 1, flag, flags) {
            return true;
        }
        flags.pop();
    }
    false
}

fn flagged_schur(lambda: &[usize], flags: &[usize], num_vars: usize) -> MultiPoly<i64> {
    assert_eq!(lambda.len(), flags.len());
    assert!(flags.iter().all(|&b| b <= num_vars));

    let mut terms = BTreeMap::new();
    let mut filling = lambda
        .iter()
        .map(|&row_len| vec![0usize; row_len])
        .collect::<Vec<_>>();
    enumerate_cell(lambda, flags, num_vars, &mut filling, 0, 0, &mut terms);
    MultiPoly::from_terms(num_vars, terms)
}

fn flagged_skew_schur(
    lambda: &[usize],
    mu: &[usize],
    flags: &[usize],
    num_vars: usize,
) -> MultiPoly<i64> {
    assert_eq!(lambda.len(), mu.len());
    assert_eq!(lambda.len(), flags.len());
    if lambda.is_empty() {
        return MultiPoly::constant(num_vars, 1);
    }

    let max_width = *lambda.iter().max().unwrap();
    let mut filling = vec![vec![0usize; max_width]; lambda.len()];
    let mut cells = Vec::new();
    for row in 0..lambda.len() {
        assert!(mu[row] <= lambda[row]);
        for col in mu[row]..lambda[row] {
            cells.push((row, col));
        }
    }

    let mut terms = BTreeMap::new();
    enumerate_skew_cell(
        lambda,
        mu,
        flags,
        num_vars,
        &cells,
        &mut filling,
        0,
        &mut terms,
    );
    MultiPoly::from_terms(num_vars, terms)
}

fn enumerate_skew_cell(
    lambda: &[usize],
    mu: &[usize],
    flags: &[usize],
    num_vars: usize,
    cells: &[(usize, usize)],
    filling: &mut [Vec<usize>],
    cell_index: usize,
    terms: &mut BTreeMap<Vec<u32>, i64>,
) {
    if cell_index == cells.len() {
        let mut exp = vec![0u32; num_vars];
        for &(row, col) in cells {
            exp[filling[row][col] - 1] += 1;
        }
        *terms.entry(exp).or_insert(0) += 1;
        return;
    }

    let (row, col) = cells[cell_index];
    let row_min = if col == mu[row] {
        1
    } else {
        filling[row][col - 1]
    };
    let col_min = if row == 0 || col < mu[row - 1] || lambda[row - 1] <= col {
        1
    } else {
        filling[row - 1][col] + 1
    };
    let min_entry = row_min.max(col_min);

    for entry in min_entry..=flags[row] {
        filling[row][col] = entry;
        enumerate_skew_cell(
            lambda,
            mu,
            flags,
            num_vars,
            cells,
            filling,
            cell_index + 1,
            terms,
        );
    }
    filling[row][col] = 0;
}

fn complete_h(num_vars: usize, degree: u32, first_var: usize, last_var: usize) -> MultiPoly<i64> {
    let mut result = BTreeMap::new();
    let mut exp = vec![0u32; num_vars];
    complete_h_inner(degree, first_var - 1, last_var - 1, &mut exp, &mut result);
    MultiPoly::from_terms(num_vars, result)
}

fn complete_h_inner(
    remaining: u32,
    var: usize,
    last_var: usize,
    exp: &mut [u32],
    result: &mut BTreeMap<Vec<u32>, i64>,
) {
    if var == last_var {
        exp[var] = remaining;
        result.insert(exp.to_vec(), 1);
        exp[var] = 0;
        return;
    }

    for power in 0..=remaining {
        exp[var] = power;
        complete_h_inner(remaining - power, var + 1, last_var, exp, result);
    }
    exp[var] = 0;
}

fn monomial(num_vars: usize, factors: &[(usize, u32)]) -> MultiPoly<i64> {
    let mut exp = vec![0u32; num_vars];
    for &(var, power) in factors {
        exp[var - 1] += power;
    }
    MultiPoly::x_power(num_vars, exp)
}

fn monomial_shift(num_vars: usize, factors: &[(usize, u32)]) -> Vec<u32> {
    let mut exp = vec![0u32; num_vars];
    for &(var, power) in factors {
        exp[var - 1] += power;
    }
    exp
}

fn enumerate_cell(
    lambda: &[usize],
    flags: &[usize],
    num_vars: usize,
    filling: &mut [Vec<usize>],
    row: usize,
    col: usize,
    terms: &mut BTreeMap<Vec<u32>, i64>,
) {
    if row == lambda.len() {
        let mut exp = vec![0u32; num_vars];
        for row_values in filling.iter() {
            for &entry in row_values {
                exp[entry - 1] += 1;
            }
        }
        *terms.entry(exp).or_insert(0) += 1;
        return;
    }

    if col == lambda[row] {
        enumerate_cell(lambda, flags, num_vars, filling, row + 1, 0, terms);
        return;
    }

    let row_min = if col == 0 { 1 } else { filling[row][col - 1] };
    let col_min = if row == 0 || col >= lambda[row - 1] {
        1
    } else {
        filling[row - 1][col] + 1
    };
    let min_entry = row_min.max(col_min);

    for entry in min_entry..=flags[row] {
        filling[row][col] = entry;
        enumerate_cell(lambda, flags, num_vars, filling, row, col + 1, terms);
    }
    filling[row][col] = 0;
}

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

fn from_lehmer_code(code: &[usize]) -> Vec<usize> {
    let mut available: Vec<usize> = (1..=code.len()).collect();
    let mut perm = Vec::with_capacity(code.len());
    for &c in code {
        perm.push(available.remove(c));
    }
    perm
}

fn inversion_count(perm: &[usize]) -> usize {
    let mut result = 0;
    for i in 0..perm.len() {
        for j in i + 1..perm.len() {
            if perm[i] > perm[j] {
                result += 1;
            }
        }
    }
    result
}

fn top_operator_word(perm: &[usize]) -> Vec<usize> {
    let w0 = (1..=perm.len()).rev().collect::<Vec<_>>();
    let inv = inverse_perm(perm);
    reduced_word(&compose_perm(&inv, &w0))
}

fn reduced_word(perm: &[usize]) -> Vec<usize> {
    let mut pi = perm.to_vec();
    let mut word = Vec::new();
    loop {
        let mut found = false;
        for i in 0..pi.len().saturating_sub(1) {
            if pi[i] > pi[i + 1] {
                word.push(i + 1);
                pi.swap(i, i + 1);
                found = true;
                break;
            }
        }
        if !found {
            break;
        }
    }
    word
}

fn inverse_perm(perm: &[usize]) -> Vec<usize> {
    let mut inverse = vec![0; perm.len()];
    for (idx, &value) in perm.iter().enumerate() {
        inverse[value - 1] = idx + 1;
    }
    inverse
}

fn compose_perm(a: &[usize], b: &[usize]) -> Vec<usize> {
    assert_eq!(a.len(), b.len());
    b.iter().map(|&value| a[value - 1]).collect()
}

fn parse_perm(text: &str) -> Vec<usize> {
    text.split(',')
        .map(|s| {
            s.parse::<usize>()
                .expect("permutation entries are integers")
        })
        .collect()
}

fn all_perms(n: usize) -> Vec<Vec<usize>> {
    let mut perm = (1..=n).collect::<Vec<_>>();
    let mut result = Vec::new();
    loop {
        result.push(perm.clone());
        if !next_permutation(&mut perm) {
            break;
        }
    }
    result
}

fn next_permutation<T: Ord>(perm: &mut [T]) -> bool {
    if perm.len() <= 1 {
        return false;
    }

    let mut i = perm.len() - 2;
    while perm[i] >= perm[i + 1] {
        if i == 0 {
            return false;
        }
        i -= 1;
    }

    let mut j = perm.len() - 1;
    while perm[j] <= perm[i] {
        j -= 1;
    }

    perm.swap(i, j);
    perm[i + 1..].reverse();
    true
}

fn trim_u32(values: &[u32]) -> Vec<u32> {
    let mut result = values.to_vec();
    while result.last() == Some(&0) {
        result.pop();
    }
    result
}
