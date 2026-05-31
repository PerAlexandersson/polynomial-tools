use crate::{CheckedRange, EvaluationDraft};
use anyhow::{Context, Result};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::panic::{catch_unwind, UnwindSafe};

type ComputeFn = fn(usize) -> Result<Vec<BigInt>>;

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct PolynomialFamilyInfo {
    pub id: String,
    pub label: String,
    pub symbol: String,
    pub source: String,
    pub min_n: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct ComputedPolynomial {
    pub family_id: String,
    pub label: String,
    pub symbol: String,
    pub n: usize,
    pub coefficients: Vec<String>,
    pub degree: usize,
    pub polynomial: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct FamilyCheckItem {
    pub n: usize,
    pub real_rooted: bool,
    pub degree: usize,
    pub polynomial: String,
    pub coefficients: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct CheckFamilyRealRootednessReport {
    pub family_id: String,
    pub n_min: usize,
    pub n_max: usize,
    pub all_real_rooted: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub first_failure_n: Option<usize>,
    pub items: Vec<FamilyCheckItem>,
}

#[derive(Clone)]
struct FamilyEntry {
    id: &'static str,
    label: &'static str,
    symbol: &'static str,
    source: &'static str,
    min_n: usize,
    compute: ComputeFn,
}

#[derive(Clone)]
pub struct PolynomialFamilyRegistry {
    families: BTreeMap<&'static str, FamilyEntry>,
}

impl PolynomialFamilyRegistry {
    pub fn new() -> Self {
        Self {
            families: BTreeMap::new(),
        }
    }

    pub fn register(
        &mut self,
        id: &'static str,
        label: &'static str,
        symbol: &'static str,
        source: &'static str,
        min_n: usize,
        compute: ComputeFn,
    ) {
        self.families.insert(
            id,
            FamilyEntry {
                id,
                label,
                symbol,
                source,
                min_n,
                compute,
            },
        );
    }

    pub fn list(&self) -> Vec<PolynomialFamilyInfo> {
        self.families.values().map(FamilyEntry::info).collect()
    }

    pub fn contains(&self, family_id: &str) -> bool {
        self.families.contains_key(family_id)
    }

    pub fn compute(&self, family_id: &str, n: usize) -> Result<ComputedPolynomial> {
        let entry = self
            .families
            .get(family_id)
            .with_context(|| format!("unknown polynomial family id '{family_id}'"))?;
        if n < entry.min_n {
            anyhow::bail!(
                "family '{}' is defined in this registry only for n >= {}",
                entry.id,
                entry.min_n
            );
        }
        let coefficients = (entry.compute)(n)
            .with_context(|| format!("failed to compute family '{}' at n={n}", entry.id))?;
        Ok(computed_polynomial(entry, n, coefficients))
    }

    pub fn check_real_rooted(
        &self,
        family_id: &str,
        n_min: usize,
        n_max: usize,
    ) -> Result<CheckFamilyRealRootednessReport> {
        if n_min > n_max {
            anyhow::bail!("expected n_min <= n_max");
        }
        let mut items = Vec::new();
        let mut first_failure_n = None;
        for n in n_min..=n_max {
            let computed = self.compute(family_id, n)?;
            let bigint_coeffs = computed
                .coefficients
                .iter()
                .map(|coefficient| coefficient.parse::<BigInt>())
                .collect::<Result<Vec<_>, _>>()
                .with_context(|| "internal coefficient serialization failure")?;
            let real_rooted = polynomial_tools::is_real_rooted_bigint_coeffs(&bigint_coeffs);
            if !real_rooted && first_failure_n.is_none() {
                first_failure_n = Some(n);
            }
            items.push(FamilyCheckItem {
                n,
                real_rooted,
                degree: computed.degree,
                polynomial: computed.polynomial,
                coefficients: computed.coefficients,
            });
        }
        Ok(CheckFamilyRealRootednessReport {
            family_id: family_id.to_string(),
            n_min,
            n_max,
            all_real_rooted: first_failure_n.is_none(),
            first_failure_n,
            items,
        })
    }
}

impl Default for PolynomialFamilyRegistry {
    fn default() -> Self {
        default_family_registry()
    }
}

impl FamilyEntry {
    fn info(&self) -> PolynomialFamilyInfo {
        PolynomialFamilyInfo {
            id: self.id.to_string(),
            label: self.label.to_string(),
            symbol: self.symbol.to_string(),
            source: self.source.to_string(),
            min_n: self.min_n,
        }
    }
}

pub fn default_family_registry() -> PolynomialFamilyRegistry {
    let mut registry = PolynomialFamilyRegistry::new();
    registry.register(
        "eulerian_polynomial",
        "Eulerian polynomial",
        "A_n(t)",
        "polynomial-tools::sequences",
        1,
        eulerian_polynomial,
    );
    registry.register(
        "narayana_polynomial",
        "Narayana polynomial",
        "N_n(t)",
        "polynomial-tools::sequences",
        1,
        narayana_polynomial,
    );
    registry.register(
        "type_b_eulerian_polynomial",
        "Type B Eulerian polynomial",
        "B_n(t)",
        "polynomial-tools::sequences",
        0,
        type_b_eulerian_polynomial,
    );
    registry.register(
        "chebyshev_t_polynomial",
        "Chebyshev polynomial of the first kind",
        "T_n(t)",
        "polynomial-tools::sequences",
        0,
        chebyshev_t_polynomial,
    );
    registry.register(
        "chebyshev_u_polynomial",
        "Chebyshev polynomial of the second kind",
        "U_n(t)",
        "polynomial-tools::sequences",
        0,
        chebyshev_u_polynomial,
    );
    registry.register(
        "hermite_polynomial",
        "Probabilist Hermite polynomial",
        "He_n(t)",
        "polynomial-tools::sequences",
        0,
        hermite_polynomial,
    );
    registry.register(
        "derangement_descent_polynomial",
        "Derangement descent polynomial",
        "D_n(t)",
        "polynomial-lab::families",
        1,
        derangement_descent_polynomial,
    );
    registry.register(
        "normalized_derangement_descent_polynomial",
        "Normalized derangement descent polynomial",
        "u_n(t)",
        "polynomial-lab::families",
        2,
        normalized_derangement_descent_polynomial,
    );
    registry.register(
        "reciprocal_eulerian_derivative_polynomial",
        "Reciprocal Eulerian derivative polynomial",
        "V_n(t)",
        "polynomial-lab::families",
        2,
        reciprocal_eulerian_derivative_polynomial,
    );
    registry
}

pub fn real_rooted_evidence_id(
    relation: &str,
    first_failure_n: Option<usize>,
    n_min: usize,
    n_max: usize,
) -> String {
    match first_failure_n {
        Some(n) => format!("{relation}_first_failure_n_{n}"),
        None => format!("{relation}_n_{n_min}_{n_max}"),
    }
}

pub fn real_rooted_evaluation_draft(
    id: String,
    relation_id: String,
    family_id: &str,
    report: &CheckFamilyRealRootednessReport,
) -> Result<EvaluationDraft> {
    let checked_range = Some(CheckedRange {
        n_min: i64::try_from(report.n_min)?,
        n_max: i64::try_from(report.n_max)?,
    });
    let mut extra = BTreeMap::new();
    extra.insert(
        "family_id".to_string(),
        Value::String(family_id.to_string()),
    );
    extra.insert("item_count".to_string(), json!(report.items.len()));
    extra.insert(
        "all_real_rooted".to_string(),
        Value::Bool(report.all_real_rooted),
    );

    let (status, first_failure, failure_reason) = if let Some(n) = report.first_failure_n {
        extra.insert("first_failure_n".to_string(), json!(n));
        (
            "counterexample_found".to_string(),
            Some(json!({ "n": n })),
            Some(format!("{family_id} is not real-rooted at n={n}")),
        )
    } else {
        ("holds_for_checked_domain".to_string(), None, None)
    };

    Ok(EvaluationDraft {
        id,
        relation_id,
        status,
        method: Some(
            "polynomial-lab family registry + polynomial-tools real-rootedness".to_string(),
        ),
        notes: Some(format!(
            "Checked real-rootedness of {family_id} for n={}..{}.",
            report.n_min, report.n_max
        )),
        checked_range,
        first_failure,
        failure_reason,
        timeout_seconds: None,
        extra,
    })
}

fn computed_polynomial(
    entry: &FamilyEntry,
    n: usize,
    coefficients: Vec<BigInt>,
) -> ComputedPolynomial {
    let coefficients = trim_bigint(coefficients);
    ComputedPolynomial {
        family_id: entry.id.to_string(),
        label: entry.label.to_string(),
        symbol: entry.symbol.to_string(),
        n,
        degree: degree_bigint(&coefficients),
        polynomial: format_bigint_poly(&coefficients),
        coefficients: coefficients.iter().map(ToString::to_string).collect(),
    }
}

fn eulerian_polynomial(n: usize) -> Result<Vec<BigInt>> {
    indexed_from_one(n, polynomial_tools::sequences::eulerian_polynomials)
}

fn narayana_polynomial(n: usize) -> Result<Vec<BigInt>> {
    indexed_from_one(n, polynomial_tools::sequences::narayana_polynomials)
}

fn type_b_eulerian_polynomial(n: usize) -> Result<Vec<BigInt>> {
    indexed_from_zero(n, polynomial_tools::sequences::type_b_eulerian_polynomials)
}

fn chebyshev_t_polynomial(n: usize) -> Result<Vec<BigInt>> {
    indexed_from_zero(n, polynomial_tools::sequences::chebyshev_polynomials_t)
}

fn chebyshev_u_polynomial(n: usize) -> Result<Vec<BigInt>> {
    indexed_from_zero(n, polynomial_tools::sequences::chebyshev_polynomials_u)
}

fn hermite_polynomial(n: usize) -> Result<Vec<BigInt>> {
    indexed_from_zero(n, polynomial_tools::sequences::hermite_polynomials)
}

fn indexed_from_one<F>(n: usize, sequence: F) -> Result<Vec<BigInt>>
where
    F: FnOnce(usize) -> Vec<Vec<i64>> + UnwindSafe,
{
    let polys = catch_sequence_overflow(|| sequence(n))?;
    polys
        .get(n - 1)
        .map(|coefficients| i64_coeffs_to_bigint(coefficients))
        .with_context(|| format!("sequence did not return polynomial for n={n}"))
}

fn indexed_from_zero<F>(n: usize, sequence: F) -> Result<Vec<BigInt>>
where
    F: FnOnce(usize) -> Vec<Vec<i64>> + UnwindSafe,
{
    let polys = catch_sequence_overflow(|| sequence(n))?;
    polys
        .get(n)
        .map(|coefficients| i64_coeffs_to_bigint(coefficients))
        .with_context(|| format!("sequence did not return polynomial for n={n}"))
}

fn catch_sequence_overflow<F>(f: F) -> Result<Vec<Vec<i64>>>
where
    F: FnOnce() -> Vec<Vec<i64>> + UnwindSafe,
{
    catch_unwind(f).map_err(|_| anyhow::anyhow!("standard sequence overflowed i64 coefficients"))
}

fn i64_coeffs_to_bigint(coefficients: &[i64]) -> Vec<BigInt> {
    coefficients.iter().map(|&c| BigInt::from(c)).collect()
}

fn derangement_descent_polynomial(n: usize) -> Result<Vec<BigInt>> {
    Ok(derangement_descent_polynomials(n)
        .pop()
        .unwrap_or_else(|| vec![BigInt::zero()]))
}

fn normalized_derangement_descent_polynomial(n: usize) -> Result<Vec<BigInt>> {
    Ok(div_t(&derangement_descent_polynomial(n)?))
}

fn reciprocal_eulerian_derivative_polynomial(n: usize) -> Result<Vec<BigInt>> {
    let a_n = eulerian_polynomial(n)?;
    Ok(reciprocal(&derivative(&a_n)))
}

fn derangement_descent_polynomials(max_n: usize) -> Vec<Vec<BigInt>> {
    if max_n == 0 {
        return vec![vec![BigInt::one()]];
    }
    let mut d = vec![vec![BigInt::zero()]];
    if max_n >= 2 {
        d.push(vec![BigInt::zero(), BigInt::one()]);
    }
    for n in 3..=max_n {
        let prev = &d[n - 2];
        let mut curr = vec![BigInt::zero(); n];
        for k in 0..=n - 2 {
            let a = if k < prev.len() {
                BigInt::from(k + 1) * &prev[k]
            } else {
                BigInt::zero()
            };
            let b = if k >= 1 {
                BigInt::from(n - k) * &prev[k - 1]
            } else {
                BigInt::zero()
            };
            curr[k] = a + b;
        }
        curr[n - 1] = if n % 2 == 0 {
            BigInt::one()
        } else {
            BigInt::zero()
        };
        d.push(trim_bigint(curr));
    }
    d
}

fn derivative(p: &[BigInt]) -> Vec<BigInt> {
    let mut out = Vec::with_capacity(p.len().saturating_sub(1));
    for (i, c) in p.iter().enumerate().skip(1) {
        out.push(BigInt::from(i) * c);
    }
    trim_bigint(out)
}

fn reciprocal(p: &[BigInt]) -> Vec<BigInt> {
    let mut out = p.to_vec();
    out.reverse();
    trim_bigint(out)
}

fn div_t(p: &[BigInt]) -> Vec<BigInt> {
    if p.len() <= 1 {
        vec![BigInt::zero()]
    } else {
        trim_bigint(p[1..].to_vec())
    }
}

fn trim_bigint(mut p: Vec<BigInt>) -> Vec<BigInt> {
    while p.len() > 1 && p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
    if p.is_empty() {
        vec![BigInt::zero()]
    } else {
        p
    }
}

fn degree_bigint(p: &[BigInt]) -> usize {
    p.iter().rposition(|c| !c.is_zero()).unwrap_or(0)
}

fn format_bigint_poly(coefficients: &[BigInt]) -> String {
    let mut terms = Vec::new();
    for (i, coefficient) in coefficients.iter().enumerate() {
        if coefficient.is_zero() {
            continue;
        }
        let abs = if coefficient < &BigInt::zero() {
            -coefficient
        } else {
            coefficient.clone()
        };
        let body = match (abs == BigInt::one(), i) {
            (_, 0) => abs.to_string(),
            (true, 1) => "t".to_string(),
            (false, 1) => format!("{abs}t"),
            (true, _) => format!("t^{i}"),
            (false, _) => format!("{abs}t^{i}"),
        };
        terms.push((coefficient < &BigInt::zero(), body));
    }
    if terms.is_empty() {
        return "0".to_string();
    }
    let mut out = String::new();
    for (index, (negative, body)) in terms.into_iter().enumerate() {
        match (index, negative) {
            (0, true) => out.push('-'),
            (0, false) => {}
            (_, true) => out.push_str(" - "),
            (_, false) => out.push_str(" + "),
        }
        out.push_str(&body);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_contains_standard_and_project_families() {
        let registry = default_family_registry();
        assert!(registry.contains("eulerian_polynomial"));
        assert!(registry.contains("derangement_descent_polynomial"));
        assert!(registry.contains("normalized_derangement_descent_polynomial"));
    }

    #[test]
    fn computes_known_values() {
        let registry = default_family_registry();
        assert_eq!(
            registry
                .compute("eulerian_polynomial", 4)
                .expect("eulerian")
                .coefficients,
            vec!["1", "11", "11", "1"]
        );
        assert_eq!(
            registry
                .compute("derangement_descent_polynomial", 4)
                .expect("derangement")
                .coefficients,
            vec!["0", "4", "4", "1"]
        );
        assert_eq!(
            registry
                .compute("normalized_derangement_descent_polynomial", 4)
                .expect("normalized derangement")
                .coefficients,
            vec!["4", "4", "1"]
        );
        assert_eq!(
            registry
                .compute("reciprocal_eulerian_derivative_polynomial", 4)
                .expect("reciprocal derivative")
                .coefficients,
            vec!["3", "22", "11"]
        );
    }

    #[test]
    fn checks_derangement_real_rootedness_range() {
        let registry = default_family_registry();
        let report = registry
            .check_real_rooted("derangement_descent_polynomial", 2, 8)
            .expect("check real-rootedness");
        assert!(report.all_real_rooted);
        assert_eq!(report.items.len(), 7);
    }
}
