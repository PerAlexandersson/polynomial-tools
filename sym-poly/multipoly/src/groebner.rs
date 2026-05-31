//! Buchberger algorithm for small exact quotient computations.
//!
//! This is still a classical Buchberger implementation, not an F4/F5 engine,
//! but it includes the first criteria and caches needed by quotient-module
//! computations.

use std::collections::{BTreeMap, HashSet};

use sym_poly_core::{Field, Ring};

use crate::division::{
    matrix_normal_forms_with_leading_terms, multiply_by_monomial, normal_form,
    normal_form_with_leading_terms,
};
use crate::monomial_order::{leading_term, monomial_divides, LeadingTerm, MonomialOrder};
use crate::MultiPoly;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GroebnerOptions {
    /// Buchberger product criterion: if leading monomials are relatively
    /// prime, the corresponding S-polynomial reduces to zero.
    pub skip_relatively_prime_leading_terms: bool,
    /// Buchberger chain criterion using already handled S-pairs.
    pub skip_chain_criterion_pairs: bool,
    /// Reduce all currently minimal lcm-degree S-polynomials in one matrix
    /// batch. Disable this to compare against scalar Buchberger reduction.
    pub batch_reduce_same_lcm_degree: bool,
}

impl Default for GroebnerOptions {
    fn default() -> Self {
        Self {
            skip_relatively_prime_leading_terms: true,
            skip_chain_criterion_pairs: true,
            batch_reduce_same_lcm_degree: true,
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BuchbergerStats {
    pub input_generators: usize,
    pub monic_generators: usize,
    pub pairs_processed: usize,
    pub pairs_skipped_relatively_prime: usize,
    pub pairs_skipped_chain_criterion: usize,
    pub matrix_reduction_batches: usize,
    pub matrix_reduction_polynomials: usize,
    pub matrix_zero_rechecks: usize,
    pub matrix_zero_rechecks_nonzero: usize,
    pub batch_fallbacks: usize,
    pub zero_remainders: usize,
    pub nonzero_remainders: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroebnerComputation<C: Ring> {
    pub basis: Vec<MultiPoly<C>>,
    pub stats: BuchbergerStats,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GroebnerBasis<C: Ring> {
    pub num_vars: usize,
    pub order: MonomialOrder,
    pub generators: Vec<MultiPoly<C>>,
    pub leading_terms: Vec<LeadingTerm<C>>,
    pub stats: BuchbergerStats,
}

impl<C: Field> GroebnerBasis<C> {
    pub fn new(generators: Vec<MultiPoly<C>>, order: MonomialOrder) -> Self {
        Self::with_options(generators, order, GroebnerOptions::default())
    }

    pub fn with_options(
        generators: Vec<MultiPoly<C>>,
        order: MonomialOrder,
        options: GroebnerOptions,
    ) -> Self {
        let input_num_vars = generators.first().map(MultiPoly::num_vars).unwrap_or(0);
        let computation = reduced_groebner_basis_with_stats(&generators, order, options);
        let leading_terms = computation
            .basis
            .iter()
            .map(|polynomial| leading_term(polynomial, order).expect("nonzero polynomial"))
            .collect();
        let num_vars = computation
            .basis
            .first()
            .map(MultiPoly::num_vars)
            .unwrap_or(input_num_vars);
        Self {
            num_vars,
            order,
            generators: computation.basis,
            leading_terms,
            stats: computation.stats,
        }
    }

    pub fn normal_form(&self, polynomial: &MultiPoly<C>) -> MultiPoly<C> {
        normal_form_with_leading_terms(
            polynomial,
            &self.generators,
            &self.leading_terms,
            self.order,
        )
    }
}

pub fn s_polynomial<C: Field>(
    f: &MultiPoly<C>,
    g: &MultiPoly<C>,
    order: MonomialOrder,
) -> MultiPoly<C> {
    assert_eq!(
        f.num_vars(),
        g.num_vars(),
        "polynomials have different rings"
    );
    let lt_f = leading_term(f, order).expect("first polynomial is zero");
    let lt_g = leading_term(g, order).expect("second polynomial is zero");
    let lcm = monomial_lcm(&lt_f.exponents, &lt_g.exponents);
    let exp_f = monomial_difference(&lcm, &lt_f.exponents);
    let exp_g = monomial_difference(&lcm, &lt_g.exponents);

    multiply_by_monomial(f, &exp_f, C::one() / lt_f.coefficient)
        - multiply_by_monomial(g, &exp_g, C::one() / lt_g.coefficient)
}

fn s_polynomial_with_leading_terms<C: Field>(
    f: &MultiPoly<C>,
    g: &MultiPoly<C>,
    lt_f: &LeadingTerm<C>,
    lt_g: &LeadingTerm<C>,
) -> MultiPoly<C> {
    assert_eq!(
        f.num_vars(),
        g.num_vars(),
        "polynomials have different rings"
    );
    let lcm = monomial_lcm(&lt_f.exponents, &lt_g.exponents);
    let exp_f = monomial_difference(&lcm, &lt_f.exponents);
    let exp_g = monomial_difference(&lcm, &lt_g.exponents);

    multiply_by_monomial(f, &exp_f, C::one() / lt_f.coefficient.clone())
        - multiply_by_monomial(g, &exp_g, C::one() / lt_g.coefficient.clone())
}

pub fn buchberger_basis<C: Field>(
    generators: &[MultiPoly<C>],
    order: MonomialOrder,
) -> Vec<MultiPoly<C>> {
    buchberger_basis_with_options(generators, order, GroebnerOptions::default())
}

pub fn buchberger_basis_with_options<C: Field>(
    generators: &[MultiPoly<C>],
    order: MonomialOrder,
    options: GroebnerOptions,
) -> Vec<MultiPoly<C>> {
    buchberger_basis_with_stats(generators, order, options).basis
}

pub fn buchberger_basis_with_stats<C: Field>(
    generators: &[MultiPoly<C>],
    order: MonomialOrder,
    options: GroebnerOptions,
) -> GroebnerComputation<C> {
    if generators.is_empty() {
        return GroebnerComputation {
            basis: Vec::new(),
            stats: BuchbergerStats::default(),
        };
    }
    let num_vars = generators[0].num_vars();
    assert!(
        generators
            .iter()
            .all(|polynomial| polynomial.num_vars() == num_vars),
        "all generators must have the same number of variables"
    );

    let mut basis: Vec<_> = generators
        .iter()
        .filter(|polynomial| !polynomial.is_zero())
        .map(|polynomial| make_monic(polynomial, order))
        .collect();
    let mut basis_lts = basis
        .iter()
        .map(|polynomial| leading_term(polynomial, order).expect("nonzero polynomial"))
        .collect::<Vec<_>>();
    let mut pairs = CriticalPairQueue::new(&basis_lts);
    let mut handled_pairs = HashSet::new();
    let mut stats = BuchbergerStats {
        input_generators: generators.len(),
        monic_generators: basis.len(),
        ..BuchbergerStats::default()
    };

    while !pairs.is_empty() {
        let pair_batch = pairs.take_next_batch(options.batch_reduce_same_lcm_degree);
        let mut active_pairs = Vec::new();
        let mut s_polynomials = Vec::new();

        for (i, j) in pair_batch {
            stats.pairs_processed += 1;
            let pair = canonical_pair(i, j);
            if options.skip_relatively_prime_leading_terms
                && leading_monomials_are_relatively_prime(&basis_lts[i], &basis_lts[j])
            {
                stats.pairs_skipped_relatively_prime += 1;
                handled_pairs.insert(pair);
                continue;
            }
            if options.skip_chain_criterion_pairs
                && chain_criterion_applies(i, j, &basis_lts, &handled_pairs)
            {
                stats.pairs_skipped_chain_criterion += 1;
                handled_pairs.insert(pair);
                continue;
            }

            active_pairs.push(pair);
            s_polynomials.push(s_polynomial_with_leading_terms(
                &basis[i],
                &basis[j],
                &basis_lts[i],
                &basis_lts[j],
            ));
        }

        if s_polynomials.is_empty() {
            continue;
        }

        let remainders = if s_polynomials.len() == 1 || !options.batch_reduce_same_lcm_degree {
            s_polynomials
                .iter()
                .map(|s| normal_form_with_leading_terms(s, &basis, &basis_lts, order))
                .collect()
        } else {
            stats.matrix_reduction_batches += 1;
            stats.matrix_reduction_polynomials += s_polynomials.len();
            let mut remainders =
                matrix_normal_forms_with_leading_terms(&s_polynomials, &basis, &basis_lts, order);
            for (s, remainder) in s_polynomials.iter().zip(remainders.iter_mut()) {
                if !remainder.is_zero() {
                    continue;
                }
                stats.matrix_zero_rechecks += 1;
                let scalar_remainder = normal_form_with_leading_terms(s, &basis, &basis_lts, order);
                if !scalar_remainder.is_zero() {
                    stats.matrix_zero_rechecks_nonzero += 1;
                    *remainder = scalar_remainder;
                }
            }
            remainders
        };

        for (pair, remainder) in active_pairs.into_iter().zip(remainders) {
            handled_pairs.insert(pair);
            if remainder.is_zero() {
                stats.zero_remainders += 1;
                continue;
            }
            stats.nonzero_remainders += 1;
            let new_index = basis.len();
            let new_polynomial = make_monic(&remainder, order);
            let new_lt = leading_term(&new_polynomial, order).expect("nonzero polynomial");
            basis.push(new_polynomial);
            basis_lts.push(new_lt);
            for i in 0..new_index {
                pairs.insert(i, new_index, &basis_lts);
            }
        }
    }

    if options.batch_reduce_same_lcm_degree && !is_groebner_basis(&basis, order) {
        let mut fallback_options = options;
        fallback_options.batch_reduce_same_lcm_degree = false;
        let mut fallback = buchberger_basis_with_stats(generators, order, fallback_options);
        fallback.stats.batch_fallbacks += 1;
        return fallback;
    }

    GroebnerComputation { basis, stats }
}

pub fn reduced_groebner_basis<C: Field>(
    generators: &[MultiPoly<C>],
    order: MonomialOrder,
) -> Vec<MultiPoly<C>> {
    reduced_groebner_basis_with_options(generators, order, GroebnerOptions::default())
}

pub fn reduced_groebner_basis_with_options<C: Field>(
    generators: &[MultiPoly<C>],
    order: MonomialOrder,
    options: GroebnerOptions,
) -> Vec<MultiPoly<C>> {
    reduced_groebner_basis_with_stats(generators, order, options).basis
}

pub fn reduced_groebner_basis_with_stats<C: Field>(
    generators: &[MultiPoly<C>],
    order: MonomialOrder,
    options: GroebnerOptions,
) -> GroebnerComputation<C> {
    let computation = buchberger_basis_with_stats(generators, order, options);
    let mut basis = computation.basis;
    remove_leading_monomial_multiples(&mut basis, order);
    let mut reduced = Vec::new();

    for i in 0..basis.len() {
        let others: Vec<_> = basis
            .iter()
            .enumerate()
            .filter_map(|(j, polynomial)| (i != j).then_some(polynomial.clone()))
            .collect();
        let remainder = normal_form(&basis[i], &others, order);
        if !remainder.is_zero() {
            reduced.push(make_monic(&remainder, order));
        }
    }

    reduced.sort_by(|a, b| {
        let a_lt = leading_term(a, order).expect("nonzero polynomial");
        let b_lt = leading_term(b, order).expect("nonzero polynomial");
        order.compare(&b_lt.exponents, &a_lt.exponents)
    });
    if options.batch_reduce_same_lcm_degree
        && (!is_groebner_basis(&reduced, order)
            || !generators_reduce_to_zero(generators, &reduced, order))
    {
        let mut fallback_options = options;
        fallback_options.batch_reduce_same_lcm_degree = false;
        let mut fallback = reduced_groebner_basis_with_stats(generators, order, fallback_options);
        fallback.stats.batch_fallbacks += 1;
        return fallback;
    }
    GroebnerComputation {
        basis: reduced,
        stats: computation.stats,
    }
}

pub fn is_groebner_basis<C: Field>(basis: &[MultiPoly<C>], order: MonomialOrder) -> bool {
    for (i, j) in all_pairs(basis.len()) {
        let s = s_polynomial(&basis[i], &basis[j], order);
        if !normal_form(&s, basis, order).is_zero() {
            return false;
        }
    }
    true
}

fn generators_reduce_to_zero<C: Field>(
    generators: &[MultiPoly<C>],
    basis: &[MultiPoly<C>],
    order: MonomialOrder,
) -> bool {
    generators
        .iter()
        .all(|generator| normal_form(generator, basis, order).is_zero())
}

pub fn make_monic<C: Field>(polynomial: &MultiPoly<C>, order: MonomialOrder) -> MultiPoly<C> {
    let lt = leading_term(polynomial, order).expect("zero polynomial has no monic form");
    polynomial.scale(&(C::one() / lt.coefficient))
}

fn monomial_lcm(a: &[u32], b: &[u32]) -> Vec<u32> {
    assert_eq!(a.len(), b.len(), "monomials have different lengths");
    a.iter().zip(b.iter()).map(|(&x, &y)| x.max(y)).collect()
}

fn leading_monomials_are_relatively_prime<C>(lt_f: &LeadingTerm<C>, lt_g: &LeadingTerm<C>) -> bool {
    lt_f.exponents
        .iter()
        .zip(lt_g.exponents.iter())
        .all(|(&a, &b)| a == 0 || b == 0)
}

fn monomial_difference(a: &[u32], b: &[u32]) -> Vec<u32> {
    assert_eq!(a.len(), b.len(), "monomials have different lengths");
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| {
            assert!(x >= y, "monomial difference would have negative exponent");
            x - y
        })
        .collect()
}

fn all_pairs(n: usize) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            pairs.push((i, j));
        }
    }
    pairs
}

fn canonical_pair(i: usize, j: usize) -> (usize, usize) {
    if i < j {
        (i, j)
    } else {
        (j, i)
    }
}

fn pair_lcm_total_degree<C>(i: usize, j: usize, leading_terms: &[LeadingTerm<C>]) -> u32 {
    leading_terms[i]
        .exponents
        .iter()
        .zip(leading_terms[j].exponents.iter())
        .map(|(&a, &b)| a.max(b))
        .sum()
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct CriticalPairQueue {
    by_lcm_degree: BTreeMap<u32, Vec<(usize, usize)>>,
    len: usize,
}

impl CriticalPairQueue {
    fn new<C>(leading_terms: &[LeadingTerm<C>]) -> Self {
        let mut queue = Self::default();
        for i in 0..leading_terms.len() {
            for j in i + 1..leading_terms.len() {
                queue.insert(i, j, leading_terms);
            }
        }
        queue
    }

    fn is_empty(&self) -> bool {
        self.len == 0
    }

    fn insert<C>(&mut self, i: usize, j: usize, leading_terms: &[LeadingTerm<C>]) {
        let pair = canonical_pair(i, j);
        let degree = pair_lcm_total_degree(pair.0, pair.1, leading_terms);
        self.by_lcm_degree.entry(degree).or_default().push(pair);
        self.len += 1;
    }

    fn take_next_batch(&mut self, batch_same_lcm_degree: bool) -> Vec<(usize, usize)> {
        assert!(self.len > 0, "critical pair queue is empty");
        if batch_same_lcm_degree {
            return self.take_smallest_degree_bucket();
        }
        vec![self.take_best_pair()]
    }

    fn take_smallest_degree_bucket(&mut self) -> Vec<(usize, usize)> {
        let degree = *self
            .by_lcm_degree
            .keys()
            .next()
            .expect("critical pair queue is empty");
        let mut bucket = self
            .by_lcm_degree
            .remove(&degree)
            .expect("degree bucket disappeared");
        self.len -= bucket.len();
        bucket.sort_by_key(|&(i, j)| (i.max(j), i.min(j)));
        bucket
    }

    fn take_best_pair(&mut self) -> (usize, usize) {
        let degree = *self
            .by_lcm_degree
            .keys()
            .next()
            .expect("critical pair queue is empty");
        let bucket = self
            .by_lcm_degree
            .get_mut(&degree)
            .expect("degree bucket disappeared");
        let best_index = bucket
            .iter()
            .enumerate()
            .min_by_key(|entry| {
                let &(i, j) = entry.1;
                (i.max(j), i.min(j))
            })
            .map(|(index, _)| index)
            .expect("degree bucket is empty");
        let pair = bucket.swap_remove(best_index);
        if bucket.is_empty() {
            self.by_lcm_degree.remove(&degree);
        }
        self.len -= 1;
        pair
    }
}

fn chain_criterion_applies<C>(
    i: usize,
    j: usize,
    leading_terms: &[LeadingTerm<C>],
    handled_pairs: &HashSet<(usize, usize)>,
) -> bool {
    let lcm_ij = monomial_lcm(&leading_terms[i].exponents, &leading_terms[j].exponents);
    (0..leading_terms.len()).any(|k| {
        if k == i || k == j {
            return false;
        }
        if !handled_pairs.contains(&canonical_pair(i, k))
            || !handled_pairs.contains(&canonical_pair(k, j))
        {
            return false;
        }
        let lcm_ik = monomial_lcm(&leading_terms[i].exponents, &leading_terms[k].exponents);
        let lcm_kj = monomial_lcm(&leading_terms[k].exponents, &leading_terms[j].exponents);
        monomial_divides(&lcm_ik, &lcm_ij) && monomial_divides(&lcm_kj, &lcm_ij)
    })
}

fn remove_leading_monomial_multiples<C: Ring>(basis: &mut Vec<MultiPoly<C>>, order: MonomialOrder) {
    let leading_terms: Vec<_> = basis
        .iter()
        .map(|polynomial| leading_term(polynomial, order).expect("nonzero polynomial"))
        .collect();
    let keep: Vec<_> = leading_terms
        .iter()
        .enumerate()
        .map(|(i, lt_i)| {
            !leading_terms
                .iter()
                .enumerate()
                .any(|(j, lt_j)| i != j && monomial_divides(&lt_j.exponents, &lt_i.exponents))
        })
        .collect();

    let old_basis = std::mem::take(basis);
    *basis = old_basis
        .into_iter()
        .enumerate()
        .filter_map(|(i, polynomial)| keep[i].then_some(polynomial))
        .collect();
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Ratio;
    use sym_poly_core::{qt_rational_monomial, QtRationalFunction, Ring};

    type Q = Ratio<i64>;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn mono(exponents: &[u32], coefficient: i64) -> MultiPoly<Q> {
        MultiPoly::monomial(2, exponents.to_vec(), q(coefficient))
    }

    fn mono_n(num_vars: usize, exponents: &[u32], coefficient: i64) -> MultiPoly<Q> {
        MultiPoly::monomial(num_vars, exponents.to_vec(), q(coefficient))
    }

    fn constant(value: i64) -> MultiPoly<Q> {
        MultiPoly::constant(2, q(value))
    }

    #[test]
    fn test_s_polynomial_cancels_leading_terms() {
        let f = mono(&[1, 0], 1) + mono(&[0, 1], 1);
        let g = mono(&[1, 0], 1) + constant(-1);
        let s = s_polynomial(&f, &g, MonomialOrder::Lex);

        assert_eq!(
            leading_term(&s, MonomialOrder::Lex).unwrap().exponents,
            vec![0, 1]
        );
    }

    #[test]
    fn test_buchberger_for_already_groebner_basis() {
        let g1 = mono(&[1, 0], 1) + mono(&[0, 1], 1);
        let g2 = mono(&[0, 2], 1) + constant(-1);
        let basis = buchberger_basis(&[g1, g2], MonomialOrder::Lex);

        assert!(is_groebner_basis(&basis, MonomialOrder::Lex));
    }

    #[test]
    fn test_groebner_basis_preserves_variable_count_for_zero_ideal() {
        let zero = MultiPoly::<Q>::zero(3);
        let gb = GroebnerBasis::new(vec![zero], MonomialOrder::Lex);
        let polynomial = mono_n(3, &[1, 0, 0], 1) + mono_n(3, &[0, 1, 0], 2);

        assert_eq!(gb.num_vars, 3);
        assert!(gb.generators.is_empty());
        assert!(gb.leading_terms.is_empty());
        assert_eq!(gb.normal_form(&polynomial), polynomial);
    }

    #[test]
    fn test_buchberger_supports_standard_monomial_orders() {
        let e1 = mono_n(3, &[1, 0, 0], 1) + mono_n(3, &[0, 1, 0], 1) + mono_n(3, &[0, 0, 1], 1);
        let e2 = mono_n(3, &[1, 1, 0], 1) + mono_n(3, &[1, 0, 1], 1) + mono_n(3, &[0, 1, 1], 1);
        let e3 = mono_n(3, &[1, 1, 1], 1);
        let generators = vec![e1, e2, e3];

        for order in MonomialOrder::STANDARD_ORDERS {
            let basis = reduced_groebner_basis(&generators, order);
            assert!(
                is_groebner_basis(&basis, order),
                "{order} did not produce a Groebner basis"
            );
            assert!(
                generators_reduce_to_zero(&generators, &basis, order),
                "{order} basis does not generate the original ideal"
            );
        }
    }

    #[test]
    fn test_product_criterion_skips_relatively_prime_leading_terms() {
        let g1 = mono(&[1, 0], 1) + constant(1);
        let g2 = mono(&[0, 1], 1) + constant(-1);
        let computation =
            buchberger_basis_with_stats(&[g1, g2], MonomialOrder::Lex, GroebnerOptions::default());

        assert_eq!(computation.stats.pairs_processed, 1);
        assert_eq!(computation.stats.pairs_skipped_relatively_prime, 1);
        assert_eq!(computation.stats.zero_remainders, 0);
        assert!(is_groebner_basis(&computation.basis, MonomialOrder::Lex));
    }

    #[test]
    fn test_product_criterion_can_be_disabled() {
        let g1 = mono(&[1, 0], 1) + constant(1);
        let g2 = mono(&[0, 1], 1) + constant(-1);
        let computation = buchberger_basis_with_stats(
            &[g1, g2],
            MonomialOrder::Lex,
            GroebnerOptions {
                skip_relatively_prime_leading_terms: false,
                ..GroebnerOptions::default()
            },
        );

        assert_eq!(computation.stats.pairs_processed, 1);
        assert_eq!(computation.stats.pairs_skipped_relatively_prime, 0);
        assert_eq!(computation.stats.zero_remainders, 1);
        assert!(is_groebner_basis(&computation.basis, MonomialOrder::Lex));
    }

    #[test]
    fn test_pair_lcm_total_degree_without_lcm_allocation() {
        let leading_terms = vec![
            LeadingTerm {
                exponents: vec![2, 0, 1],
                coefficient: q(1),
            },
            LeadingTerm {
                exponents: vec![1, 3, 0],
                coefficient: q(1),
            },
        ];

        assert_eq!(pair_lcm_total_degree(0, 1, &leading_terms), 6);
    }

    #[test]
    fn test_chain_criterion_skips_redundant_pair() {
        let g1 = mono(&[2, 0], 1);
        let g2 = mono(&[1, 1], 1);
        let g3 = mono(&[0, 2], 1);
        let computation = buchberger_basis_with_stats(
            &[g1, g2, g3],
            MonomialOrder::Lex,
            GroebnerOptions {
                skip_relatively_prime_leading_terms: false,
                skip_chain_criterion_pairs: true,
                ..GroebnerOptions::default()
            },
        );

        assert_eq!(computation.stats.pairs_processed, 3);
        assert_eq!(computation.stats.pairs_skipped_relatively_prime, 0);
        assert_eq!(computation.stats.pairs_skipped_chain_criterion, 1);
        assert_eq!(computation.stats.zero_remainders, 2);
        assert!(is_groebner_basis(&computation.basis, MonomialOrder::Lex));
    }

    #[test]
    fn test_batch_buchberger_matches_scalar_buchberger() {
        let f = mono(&[2, 0], 1) - mono(&[0, 1], 1);
        let g = mono(&[1, 1], 1) - constant(1);
        let h = mono(&[0, 3], 1) - constant(1);
        let generators = vec![f, g, h];
        let scalar = reduced_groebner_basis_with_options(
            &generators,
            MonomialOrder::Lex,
            GroebnerOptions {
                batch_reduce_same_lcm_degree: false,
                ..GroebnerOptions::default()
            },
        );
        let batched = reduced_groebner_basis_with_options(
            &generators,
            MonomialOrder::Lex,
            GroebnerOptions::default(),
        );

        assert_eq!(batched, scalar);
        assert!(is_groebner_basis(&batched, MonomialOrder::Lex));
    }

    #[test]
    fn test_batch_buchberger_records_matrix_reduction_stats() {
        let f = mono(&[2, 0], 1) - mono(&[0, 1], 1);
        let g = mono(&[1, 1], 1) - constant(1);
        let h = mono(&[0, 2], 1) - mono(&[1, 0], 1);
        let computation = buchberger_basis_with_stats(
            &[f, g, h],
            MonomialOrder::Lex,
            GroebnerOptions {
                skip_relatively_prime_leading_terms: false,
                skip_chain_criterion_pairs: false,
                batch_reduce_same_lcm_degree: true,
            },
        );

        assert!(computation.stats.matrix_reduction_batches > 0);
        assert!(computation.stats.matrix_reduction_polynomials >= 2);
        assert!(is_groebner_basis(&computation.basis, MonomialOrder::Lex));
    }

    #[test]
    fn test_reduced_groebner_basis_preserves_artin_s3_pure_power() {
        let e1 = mono_n(3, &[1, 0, 0], 1) + mono_n(3, &[0, 1, 0], 1) + mono_n(3, &[0, 0, 1], 1);
        let e2 = mono_n(3, &[1, 1, 0], 1) + mono_n(3, &[1, 0, 1], 1) + mono_n(3, &[0, 1, 1], 1);
        let e3 = mono_n(3, &[1, 1, 1], 1);
        let generators = vec![e1, e2, e3];
        let gb = GroebnerBasis::new(generators.clone(), MonomialOrder::Lex);

        assert!(gb
            .leading_terms
            .iter()
            .any(|leading_term| leading_term.exponents == vec![0, 0, 3]));
        assert!(generators_reduce_to_zero(
            &generators,
            &gb.generators,
            MonomialOrder::Lex
        ));
    }

    #[test]
    fn test_reduced_groebner_basis_for_two_linear_relations() {
        let g1 = mono(&[1, 0], 1) + mono(&[0, 1], 1);
        let g2 = mono(&[0, 1], 1) + constant(-1);
        let basis = reduced_groebner_basis(&[g1, g2], MonomialOrder::Lex);

        assert_eq!(
            basis,
            vec![
                mono(&[1, 0], 1) + constant(1),
                mono(&[0, 1], 1) + constant(-1),
            ]
        );
        assert!(is_groebner_basis(&basis, MonomialOrder::Lex));
    }

    #[test]
    fn test_groebner_basis_normal_form() {
        let g1 = mono(&[1, 0], 1) + mono(&[0, 1], 1);
        let g2 = mono(&[0, 1], 1) + constant(-1);
        let gb = GroebnerBasis::new(vec![g1, g2], MonomialOrder::Lex);
        let f = mono(&[2, 0], 1) + mono(&[0, 2], 1);

        assert_eq!(gb.normal_form(&f), constant(2));
    }

    #[test]
    fn test_groebner_basis_over_qt_rational_functions() {
        type K = QtRationalFunction<Q>;

        let parameter = qt_rational_monomial::<Q>(1, 0);
        let generator =
            MultiPoly::monomial(1, vec![1], parameter.clone()) - MultiPoly::constant(1, K::one());
        let gb = GroebnerBasis::new(vec![generator], MonomialOrder::Lex);
        let normal = gb.normal_form(&MultiPoly::var(1, 0));

        assert_eq!(normal, MultiPoly::constant(1, K::one() / parameter.clone()));
        assert_eq!(gb.stats.pairs_processed, 0);
    }

    #[test]
    fn test_buchberger_adds_missing_s_polynomial() {
        let f = mono(&[2, 0], 1) - mono(&[0, 1], 1);
        let g = mono(&[1, 1], 1) - constant(1);
        let basis = buchberger_basis(&[f, g], MonomialOrder::Lex);

        assert!(basis.len() > 2);
        assert!(is_groebner_basis(&basis, MonomialOrder::Lex));
    }
}
