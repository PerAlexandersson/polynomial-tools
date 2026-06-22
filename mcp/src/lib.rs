use polynomial_tools::recurrence::{
    find_recurrence_adaptive, find_recurrence_adaptive_rational, AdaptiveSearchOptions,
};
use polynomial_tools::sequences::{
    chebyshev_polynomials_t, chebyshev_polynomials_u, eulerian_polynomials, hermite_polynomials,
    narayana_polynomials, type_b_eulerian_polynomials,
};
use polynomial_tools::*;
use rmcp::{
    handler::server::{router::tool::ToolRouter, wrapper::Parameters},
    model::{Implementation, ServerCapabilities, ServerInfo},
    tool, tool_handler, tool_router, ErrorData as McpError, Json, ServerHandler,
};
use schemars::{JsonSchema, Schema, SchemaGenerator};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::borrow::Cow;

#[derive(Debug, Clone)]
pub struct PolynomialToolsServer {
    tool_router: ToolRouter<Self>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct PolynomialInput {
    pub coefficients: Option<Vec<i64>>,
    pub expression: Option<String>,
}

impl JsonSchema for PolynomialInput {
    fn schema_name() -> Cow<'static, str> {
        "PolynomialInput".into()
    }

    fn json_schema(generator: &mut SchemaGenerator) -> Schema {
        Schema::try_from(json!({
            "type": "object",
            "additionalProperties": false,
            "properties": {
                "coefficients": <Vec<i64>>::json_schema(generator),
                "expression": String::json_schema(generator)
            },
            "oneOf": [
                { "required": ["coefficients"] },
                { "required": ["expression"] }
            ]
        }))
        .expect("valid PolynomialInput schema")
    }
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PolynomialBatchInput {
    pub polynomials: Option<Vec<PolynomialInput>>,
    pub text: Option<String>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FamilyCheckOptions {
    pub require_real_rooted: Option<bool>,
    pub require_simple_roots: Option<bool>,
    pub require_palindromic: Option<bool>,
    pub require_gamma_positive: Option<bool>,
    pub require_log_concave: Option<bool>,
    pub require_ultra_log_concave: Option<bool>,
    pub check_consecutive_interlacing: Option<bool>,
    pub require_consecutive_weak_interlacing: Option<bool>,
    pub find_recurrence: Option<bool>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LaceCheckRequest {
    pub block_rows: Option<usize>,
    pub block_cols: Option<usize>,
    pub max_minor_size: Option<usize>,
    pub include_matrix: Option<bool>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CheckPolynomialFamilyRequest {
    pub polynomials: Option<Vec<PolynomialInput>>,
    pub text: Option<String>,
    pub sequence: Option<SequenceKind>,
    pub max_n: Option<usize>,
    pub options: Option<FamilyCheckOptions>,
    pub recurrence_options: Option<RecurrenceSearchOptionsInput>,
    pub lace: Option<LaceCheckRequest>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct InterlacingPairRequest {
    pub p: PolynomialInput,
    pub q: PolynomialInput,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceSearchOptionsInput {
    pub skip_prefix: Option<usize>,
    pub require_all_offsets: Option<bool>,
    pub min_rec_len: Option<usize>,
    pub max_rec_len: Option<usize>,
    pub min_var_deg: Option<usize>,
    pub max_var_deg: Option<usize>,
    pub min_idx_deg: Option<usize>,
    pub max_idx_deg: Option<usize>,
    pub min_diff_deg: Option<usize>,
    pub max_diff_deg: Option<usize>,
    pub try_inhomogeneous: Option<bool>,
    pub min_inhomo_var_deg: Option<usize>,
    pub max_inhomo_var_deg: Option<usize>,
    pub min_inhomo_idx_deg: Option<usize>,
    pub max_inhomo_idx_deg: Option<usize>,
    pub try_denominator: Option<bool>,
    pub try_alternating_sign: Option<bool>,
    pub max_denom_var_deg: Option<usize>,
    pub max_denom_idx_deg: Option<usize>,
    pub min_margin: Option<usize>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(untagged)]
pub enum RationalCoefficientInput {
    Integer(i64),
    Text(String),
}

impl From<i64> for RationalCoefficientInput {
    fn from(value: i64) -> Self {
        Self::Integer(value)
    }
}

impl JsonSchema for RationalCoefficientInput {
    fn schema_name() -> Cow<'static, str> {
        "RationalCoefficientInput".into()
    }

    fn json_schema(_generator: &mut SchemaGenerator) -> Schema {
        Schema::try_from(json!({
            "anyOf": [
                { "type": "integer" },
                {
                    "type": "string",
                    "description": "Exact rational coefficient, e.g. \"42\", \"-17\", or \"3/7\"."
                }
            ]
        }))
        .expect("valid RationalCoefficientInput schema")
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct FindRecurrenceRequest {
    pub polynomials: Option<Vec<PolynomialInput>>,
    pub coefficients: Option<Vec<Vec<RationalCoefficientInput>>>,
    pub expressions: Option<Vec<String>>,
    pub text: Option<String>,
    pub options: Option<RecurrenceSearchOptionsInput>,
}

impl JsonSchema for FindRecurrenceRequest {
    fn schema_name() -> Cow<'static, str> {
        "FindRecurrenceRequest".into()
    }

    fn json_schema(generator: &mut SchemaGenerator) -> Schema {
        Schema::try_from(json!({
            "type": "object",
            "additionalProperties": false,
            "properties": {
                "polynomials": <Vec<PolynomialInput>>::json_schema(generator),
                "coefficients": <Vec<Vec<RationalCoefficientInput>>>::json_schema(generator),
                "expressions": <Vec<String>>::json_schema(generator),
                "text": String::json_schema(generator),
                "options": {
                    "anyOf": [
                        generator.subschema_for::<RecurrenceSearchOptionsInput>(),
                        { "type": "null" }
                    ]
                }
            },
            "oneOf": [
                { "required": ["polynomials"] },
                { "required": ["coefficients"] },
                { "required": ["expressions"] },
                { "required": ["text"] }
            ]
        }))
        .expect("valid FindRecurrenceRequest schema")
    }
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EhrhartHstarRequest {
    pub mode: EhrhartHstarMode,
    pub hstar: Option<Vec<i64>>,
    pub ehrhart_coefficients: Option<Vec<String>>,
    pub numerator_coefficients: Option<Vec<i64>>,
    pub denominator: Option<i64>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum EhrhartHstarMode {
    HstarToEhrhart,
    EhrhartToHstar,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenerateSequenceRequest {
    pub sequence: SequenceKind,
    pub max_n: usize,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum SequenceKind {
    Eulerian,
    Narayana,
    TypeBEulerian,
    ChebyshevT,
    ChebyshevU,
    Hermite,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct NormalizedPolynomial {
    pub polynomial: String,
    pub coefficients: Vec<i64>,
    pub degree: usize,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct ParsePolynomialItem {
    pub index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub polynomial: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coefficients: Option<Vec<i64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub degree: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct ParsePolynomialsResponse {
    pub items: Vec<ParsePolynomialItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct PolynomialPropertiesItem {
    pub index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub polynomial: Option<NormalizedPolynomial>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub real_rooted: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub simple_roots: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub palindromic: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gamma_positive: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gamma_coefficients: Option<Vec<i64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub log_concave: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ultra_log_concave: Option<bool>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct PolynomialPropertiesResponse {
    pub items: Vec<PolynomialPropertiesItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct InterlacingPairResult {
    pub p: NormalizedPolynomial,
    pub q: NormalizedPolynomial,
    pub strict: Option<bool>,
    pub weak: Option<bool>,
    pub status: String,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct InterlacingPairItem {
    pub pair_index: usize,
    pub left_index: usize,
    pub right_index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub result: Option<InterlacingPairResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct InterlacingSequenceResponse {
    pub items: Vec<ParsePolynomialItem>,
    pub pairs: Vec<InterlacingPairItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct RealRootsItem {
    pub index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub polynomial: Option<NormalizedPolynomial>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub real_rooted: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub roots: Option<Vec<String>>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct RealRootsResponse {
    pub items: Vec<RealRootsItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct FindRecurrenceResponse {
    pub found: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub recurrence: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub latex: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mathematica: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sage: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unknowns: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub equations: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub candidates_tried: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub parse_errors: Vec<ParsePolynomialItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct ResultantResponse {
    pub p: NormalizedPolynomial,
    pub q: NormalizedPolynomial,
    pub resultant: String,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct DiscriminantItem {
    pub index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub polynomial: Option<NormalizedPolynomial>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub discriminant: Option<String>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct DiscriminantResponse {
    pub items: Vec<DiscriminantItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct EhrhartHstarResponse {
    pub mode: EhrhartHstarMode,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hstar: Option<Vec<i64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ehrhart_coefficients: Option<Vec<String>>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct DisplayedPolynomial {
    pub polynomial: String,
    pub coefficients: Vec<i64>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct MagicBasisReport {
    pub coordinates: Vec<String>,
    pub left_partial_sums: Vec<String>,
    pub right_partial_sums: Vec<String>,
    pub partial_sum_checks: Vec<bool>,
    pub all_nonnegative: bool,
    pub left_leq_right: bool,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct DecompositionReport {
    #[serde(flatten)]
    pub polynomial: NormalizedPolynomial,
    pub reciprocal: DisplayedPolynomial,
    pub a: DisplayedPolynomial,
    pub b: DisplayedPolynomial,
    pub a_real_rooted: bool,
    pub b_real_rooted: bool,
    pub b_interlaces_a: Option<bool>,
    pub reciprocal_interlaces_input: Option<bool>,
    pub alternatingly_increasing: bool,
    pub f_polynomial: DisplayedPolynomial,
    pub r_transform_of_f: DisplayedPolynomial,
    pub r_a: DisplayedPolynomial,
    pub r_b: DisplayedPolynomial,
    pub r_interlaces_f: Option<bool>,
    pub magic: MagicBasisReport,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct DecompositionItem {
    pub index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub report: Option<DecompositionReport>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct DecompositionResponse {
    pub items: Vec<DecompositionItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct GenerateSequenceResponse {
    pub sequence: SequenceKind,
    pub max_n: usize,
    pub polynomials: Vec<NormalizedPolynomial>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct FamilyPolynomialReport {
    pub index: usize,
    #[serde(flatten)]
    pub polynomial: NormalizedPolynomial,
    pub real_rooted: bool,
    pub simple_roots: bool,
    pub palindromic: bool,
    pub gamma_positive: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gamma_coefficients: Option<Vec<i64>>,
    pub log_concave: bool,
    pub ultra_log_concave: bool,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct LaceCheckResponse {
    pub block_rows: usize,
    pub block_cols: usize,
    pub max_minor_size: usize,
    pub rows: usize,
    pub columns: usize,
    pub tnn: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub matrix: Option<Vec<Vec<i64>>>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct CheckPolynomialFamilyResponse {
    pub source: String,
    pub item_count: usize,
    pub all_required_checks_passed: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub first_failure: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub parse_errors: Vec<ParsePolynomialItem>,
    pub items: Vec<FamilyPolynomialReport>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub consecutive_pairs: Vec<InterlacingPairItem>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lace: Option<LaceCheckResponse>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub recurrence: Option<FindRecurrenceResponse>,
    pub markdown: String,
}

type ParsedBatch = Vec<Result<NormalizedPolynomial, String>>;
type RationalParsedBatch = Vec<Result<Vec<BigRational>, String>>;
type FamilySourceParse = (
    String,
    Result<Vec<NormalizedPolynomial>, Vec<ParsePolynomialItem>>,
);

fn invalid_params(message: impl Into<String>) -> McpError {
    McpError::invalid_params(message.into(), None)
}

fn parse_polynomial_input(input: &PolynomialInput) -> Result<Vec<i64>, String> {
    match (&input.coefficients, &input.expression) {
        (Some(coefficients), None) => Ok(coefficients.clone()),
        (None, Some(expression)) => parse_polynomial(expression),
        (Some(_), Some(_)) => {
            Err("expected exactly one of `coefficients` or `expression`, got both".to_string())
        }
        (None, None) => Err("expected exactly one of `coefficients` or `expression`".to_string()),
    }
}

fn parse_required_polynomial(input: &PolynomialInput, name: &str) -> Result<Vec<i64>, McpError> {
    parse_polynomial_input(input).map_err(|e| invalid_params(format!("{name}: {e}")))
}

fn normalize_coefficients(mut coefficients: Vec<i64>) -> Vec<i64> {
    while coefficients.len() > 1 && coefficients.last() == Some(&0) {
        coefficients.pop();
    }
    if coefficients.is_empty() {
        coefficients.push(0);
    }
    coefficients
}

fn normalize_rational_coefficients(mut coefficients: Vec<BigRational>) -> Vec<BigRational> {
    let zero = parse_rational("0").expect("0 parses as BigRational");
    while coefficients.len() > 1 && coefficients.last() == Some(&zero) {
        coefficients.pop();
    }
    if coefficients.is_empty() {
        coefficients.push(zero);
    }
    coefficients
}

fn degree(coefficients: &[i64]) -> usize {
    coefficients.iter().rposition(|&c| c != 0).unwrap_or(0)
}

fn normalize_polynomial(coefficients: Vec<i64>) -> NormalizedPolynomial {
    let coefficients = normalize_coefficients(coefficients);
    NormalizedPolynomial {
        polynomial: format_poly(&coefficients),
        degree: degree(&coefficients),
        coefficients,
    }
}

fn parse_batch(input: &PolynomialBatchInput) -> Result<ParsedBatch, McpError> {
    match (&input.polynomials, &input.text) {
        (Some(polynomials), None) => Ok(polynomials
            .iter()
            .map(|p| parse_polynomial_input(p).map(normalize_polynomial))
            .collect()),
        (None, Some(text)) => Ok(parse_polynomials(text)
            .into_iter()
            .map(|r| r.map(normalize_polynomial))
            .collect()),
        (Some(_), Some(_)) => Err(invalid_params(
            "expected exactly one of `polynomials` or `text`, got both",
        )),
        (None, None) => Err(invalid_params(
            "expected exactly one of `polynomials` or `text`",
        )),
    }
}

fn i64_coefficients_to_rational(coefficients: &[i64]) -> Vec<BigRational> {
    coefficients
        .iter()
        .map(|coeff| parse_rational(&coeff.to_string()).expect("i64 parses as BigRational"))
        .collect()
}

fn parse_rational_coefficient(input: &RationalCoefficientInput) -> Result<BigRational, String> {
    match input {
        RationalCoefficientInput::Integer(value) => parse_rational(&value.to_string()),
        RationalCoefficientInput::Text(value) => parse_rational(value),
    }
}

fn parse_rational_coefficients(
    coefficients: &[RationalCoefficientInput],
) -> Result<Vec<BigRational>, String> {
    coefficients
        .iter()
        .enumerate()
        .map(|(index, coefficient)| {
            parse_rational_coefficient(coefficient)
                .map_err(|error| format!("coefficient {index}: {error}"))
        })
        .collect::<Result<Vec<_>, _>>()
        .map(normalize_rational_coefficients)
}

fn parse_recurrence_batch_rational(
    input: &FindRecurrenceRequest,
) -> Result<RationalParsedBatch, McpError> {
    let source_count = [
        input.polynomials.is_some(),
        input.coefficients.is_some(),
        input.expressions.is_some(),
        input.text.is_some(),
    ]
    .into_iter()
    .filter(|present| *present)
    .count();

    if source_count != 1 {
        return Err(invalid_params(
            "expected exactly one of `polynomials`, `coefficients`, `expressions`, or `text`",
        ));
    }

    if let Some(polynomials) = &input.polynomials {
        let batch = parse_batch(&PolynomialBatchInput {
            polynomials: Some(polynomials.clone()),
            text: None,
        })?;
        return Ok(batch
            .into_iter()
            .map(|item| {
                item.map(|polynomial| i64_coefficients_to_rational(&polynomial.coefficients))
            })
            .collect());
    }

    if let Some(coefficients) = &input.coefficients {
        return Ok(coefficients
            .iter()
            .map(|coefficients| parse_rational_coefficients(coefficients))
            .collect());
    }

    if let Some(expressions) = &input.expressions {
        let polynomials = expressions
            .iter()
            .cloned()
            .map(|expression| PolynomialInput {
                coefficients: None,
                expression: Some(expression),
            })
            .collect();
        let batch = parse_batch(&PolynomialBatchInput {
            polynomials: Some(polynomials),
            text: None,
        })?;
        return Ok(batch
            .into_iter()
            .map(|item| {
                item.map(|polynomial| i64_coefficients_to_rational(&polynomial.coefficients))
            })
            .collect());
    }

    let batch = parse_batch(&PolynomialBatchInput {
        polynomials: None,
        text: input.text.clone(),
    })?;
    Ok(batch
        .into_iter()
        .map(|item| item.map(|polynomial| i64_coefficients_to_rational(&polynomial.coefficients)))
        .collect())
}

fn parse_items(batch: &ParsedBatch) -> Vec<ParsePolynomialItem> {
    batch
        .iter()
        .enumerate()
        .map(|(index, item)| match item {
            Ok(polynomial) => ParsePolynomialItem {
                index,
                ok: true,
                polynomial: Some(polynomial.polynomial.clone()),
                coefficients: Some(polynomial.coefficients.clone()),
                degree: Some(polynomial.degree),
                error: None,
            },
            Err(error) => ParsePolynomialItem {
                index,
                ok: false,
                polynomial: None,
                coefficients: None,
                degree: None,
                error: Some(error.clone()),
            },
        })
        .collect()
}

fn collect_polynomials_or_errors(
    batch: ParsedBatch,
) -> Result<Vec<NormalizedPolynomial>, Vec<ParsePolynomialItem>> {
    let errors: Vec<ParsePolynomialItem> = batch
        .iter()
        .enumerate()
        .filter_map(|(index, item)| match item {
            Ok(_) => None,
            Err(error) => Some(ParsePolynomialItem {
                index,
                ok: false,
                polynomial: None,
                coefficients: None,
                degree: None,
                error: Some(error.clone()),
            }),
        })
        .collect();
    if errors.is_empty() {
        Ok(batch.into_iter().map(Result::unwrap).collect())
    } else {
        Err(errors)
    }
}

fn collect_rational_polynomials_or_errors(
    batch: RationalParsedBatch,
) -> Result<Vec<Vec<BigRational>>, Vec<ParsePolynomialItem>> {
    let errors: Vec<ParsePolynomialItem> = batch
        .iter()
        .enumerate()
        .filter_map(|(index, item)| match item {
            Ok(_) => None,
            Err(error) => Some(ParsePolynomialItem {
                index,
                ok: false,
                polynomial: None,
                coefficients: None,
                degree: None,
                error: Some(error.clone()),
            }),
        })
        .collect();
    if errors.is_empty() {
        Ok(batch.into_iter().map(Result::unwrap).collect())
    } else {
        Err(errors)
    }
}

fn interlacing_result(p: NormalizedPolynomial, q: NormalizedPolynomial) -> InterlacingPairResult {
    let strict = check_interlacing(&p.coefficients, &q.coefficients);
    let weak = check_weak_interlacing(&p.coefficients, &q.coefficients);
    let status = match (strict, weak) {
        (Some(true), _) => "strictly_interlace",
        (_, Some(true)) => "weakly_interlace",
        (Some(false), Some(false)) => "do_not_interlace",
        (Some(false), None) => "not_real_rooted_or_incompatible",
        (None, Some(false)) => "not_real_rooted_or_incompatible",
        (None, None) => "not_real_rooted_or_incompatible",
    }
    .to_string();
    InterlacingPairResult {
        p,
        q,
        strict,
        weak,
        status,
    }
}

fn family_report(index: usize, polynomial: NormalizedPolynomial) -> FamilyPolynomialReport {
    let coefficients = &polynomial.coefficients;
    FamilyPolynomialReport {
        index,
        real_rooted: is_real_rooted(coefficients),
        simple_roots: has_simple_roots(coefficients),
        palindromic: is_palindromic_ignoring_initial_zeros(coefficients),
        gamma_positive: is_gamma_positive_ignoring_initial_zeros(coefficients),
        gamma_coefficients: gamma_coefficients_ignoring_initial_zeros(coefficients),
        log_concave: is_log_concave(coefficients),
        ultra_log_concave: is_ultra_log_concave(coefficients),
        polynomial,
    }
}

fn family_source_and_polynomials(
    input: &CheckPolynomialFamilyRequest,
) -> Result<FamilySourceParse, McpError> {
    let explicit_sources = usize::from(input.polynomials.is_some())
        + usize::from(input.text.is_some())
        + usize::from(input.sequence.is_some());
    if explicit_sources != 1 {
        return Err(invalid_params(
            "expected exactly one source: `polynomials`, `text`, or `sequence`",
        ));
    }

    if let Some(sequence) = input.sequence.clone() {
        let max_n = input
            .max_n
            .ok_or_else(|| invalid_params("`max_n` is required with `sequence`"))?;
        let polynomials = generated_sequence_polynomials(&sequence, max_n)
            .into_iter()
            .map(normalize_polynomial)
            .collect();
        return Ok((
            format!("sequence:{sequence:?}, max_n={max_n}"),
            Ok(polynomials),
        ));
    }

    if input.max_n.is_some() {
        return Err(invalid_params(
            "`max_n` is only valid when `sequence` is provided",
        ));
    }

    let batch = parse_batch(&PolynomialBatchInput {
        polynomials: input.polynomials.clone(),
        text: input.text.clone(),
    })?;
    let source = if input.text.is_some() {
        "text".to_string()
    } else {
        "polynomials".to_string()
    };
    Ok((source, collect_polynomials_or_errors(batch)))
}

fn generated_sequence_polynomials(sequence: &SequenceKind, max_n: usize) -> Vec<Vec<i64>> {
    match sequence {
        SequenceKind::Eulerian => eulerian_polynomials(max_n),
        SequenceKind::Narayana => narayana_polynomials(max_n),
        SequenceKind::TypeBEulerian => type_b_eulerian_polynomials(max_n),
        SequenceKind::ChebyshevT => chebyshev_polynomials_t(max_n),
        SequenceKind::ChebyshevU => chebyshev_polynomials_u(max_n),
        SequenceKind::Hermite => hermite_polynomials(max_n),
    }
}

fn check_family_lace(
    polynomials: &[NormalizedPolynomial],
    request: &LaceCheckRequest,
) -> LaceCheckResponse {
    let coefficients: Vec<Vec<i64>> = polynomials.iter().map(|p| p.coefficients.clone()).collect();
    let max_degree = polynomials.iter().map(|p| p.degree).max().unwrap_or(0);
    let block_rows = request
        .block_rows
        .unwrap_or_else(|| polynomials.len().max(1));
    let block_cols = request
        .block_cols
        .unwrap_or_else(|| block_rows + max_degree + 1);

    match lace_matrix_sequence_i64(&coefficients, block_rows, block_cols) {
        Ok(matrix) => {
            let rows = matrix.len();
            let columns = matrix.first().map(Vec::len).unwrap_or(0);
            let max_minor_size = request
                .max_minor_size
                .unwrap_or_else(|| rows.min(columns).min(4));
            let check = check_lace_sequence_total_nonnegative_i64(
                &coefficients,
                block_rows,
                block_cols,
                max_minor_size,
            );
            LaceCheckResponse {
                block_rows,
                block_cols,
                max_minor_size,
                rows,
                columns,
                tnn: check.is_ok(),
                error: check.err(),
                matrix: request.include_matrix.unwrap_or(false).then_some(matrix),
            }
        }
        Err(error) => LaceCheckResponse {
            block_rows,
            block_cols,
            max_minor_size: request.max_minor_size.unwrap_or(0),
            rows: 0,
            columns: 0,
            tnn: false,
            error: Some(error.to_string()),
            matrix: None,
        },
    }
}

fn family_recurrence_response(
    coefficients: &[Vec<i64>],
    options: Option<RecurrenceSearchOptionsInput>,
) -> FindRecurrenceResponse {
    if coefficients.len() < 3 {
        return FindRecurrenceResponse {
            found: false,
            recurrence: None,
            latex: None,
            mathematica: None,
            sage: None,
            unknowns: None,
            equations: None,
            candidates_tried: None,
            error: Some("need at least 3 polynomials".to_string()),
            parse_errors: Vec::new(),
        };
    }

    let search = apply_recurrence_options(options);
    match find_recurrence_adaptive(coefficients, &search) {
        Some(result) => FindRecurrenceResponse {
            found: true,
            recurrence: Some(format!("{}", result.recurrence)),
            latex: Some(result.recurrence.to_latex()),
            mathematica: Some(result.recurrence.to_mathematica_definition(coefficients)),
            sage: Some(result.recurrence.to_sage_definition(coefficients)),
            unknowns: Some(result.num_unknowns),
            equations: Some(result.num_equations),
            candidates_tried: Some(result.candidates_tried),
            error: None,
            parse_errors: Vec::new(),
        },
        None => FindRecurrenceResponse {
            found: false,
            recurrence: None,
            latex: None,
            mathematica: None,
            sage: None,
            unknowns: None,
            equations: None,
            candidates_tried: None,
            error: Some("no recurrence found within the search bounds".to_string()),
            parse_errors: Vec::new(),
        },
    }
}

fn first_family_failure(
    items: &[FamilyPolynomialReport],
    pairs: &[InterlacingPairItem],
    lace: Option<&LaceCheckResponse>,
    recurrence: Option<&FindRecurrenceResponse>,
    options: Option<&FamilyCheckOptions>,
) -> Option<String> {
    let require_real_rooted = options.and_then(|o| o.require_real_rooted).unwrap_or(true);
    let require_simple_roots = options
        .and_then(|o| o.require_simple_roots)
        .unwrap_or(false);
    let require_palindromic = options.and_then(|o| o.require_palindromic).unwrap_or(false);
    let require_gamma_positive = options
        .and_then(|o| o.require_gamma_positive)
        .unwrap_or(false);
    let require_log_concave = options.and_then(|o| o.require_log_concave).unwrap_or(false);
    let require_ultra_log_concave = options
        .and_then(|o| o.require_ultra_log_concave)
        .unwrap_or(false);
    let require_interlacing = options
        .and_then(|o| o.require_consecutive_weak_interlacing)
        .unwrap_or(false);

    for item in items {
        if require_real_rooted && !item.real_rooted {
            return Some(format!("polynomial {} is not real-rooted", item.index));
        }
        if require_simple_roots && !item.simple_roots {
            return Some(format!(
                "polynomial {} does not have simple roots",
                item.index
            ));
        }
        if require_palindromic && !item.palindromic {
            return Some(format!("polynomial {} is not palindromic", item.index));
        }
        if require_gamma_positive && !item.gamma_positive {
            return Some(format!("polynomial {} is not gamma-positive", item.index));
        }
        if require_log_concave && !item.log_concave {
            return Some(format!("polynomial {} is not log-concave", item.index));
        }
        if require_ultra_log_concave && !item.ultra_log_concave {
            return Some(format!(
                "polynomial {} is not ultra-log-concave",
                item.index
            ));
        }
    }

    if require_interlacing {
        for pair in pairs {
            if pair.result.as_ref().and_then(|r| r.weak) != Some(true) {
                return Some(format!(
                    "consecutive pair {}-{} is not weakly interlacing",
                    pair.left_index, pair.right_index
                ));
            }
        }
    }

    if let Some(lace) = lace {
        if !lace.tnn {
            return Some(format!(
                "finite Lace truncation is not TNN up to minors of size {}",
                lace.max_minor_size
            ));
        }
    }

    if let Some(recurrence) = recurrence {
        if !recurrence.found {
            return Some("recurrence search did not find a recurrence".to_string());
        }
    }

    None
}

fn yes_no(value: bool) -> &'static str {
    if value {
        "yes"
    } else {
        "no"
    }
}

fn opt_yes_no(value: Option<bool>) -> &'static str {
    match value {
        Some(true) => "yes",
        Some(false) => "no",
        None => "n/a",
    }
}

fn family_markdown(
    source: &str,
    all_required_checks_passed: bool,
    first_failure: Option<&str>,
    items: &[FamilyPolynomialReport],
    pairs: &[InterlacingPairItem],
    lace: Option<&LaceCheckResponse>,
    recurrence: Option<&FindRecurrenceResponse>,
) -> String {
    let mut out = String::new();
    out.push_str("# Polynomial family check\n\n");
    out.push_str(&format!("- Source: `{source}`\n"));
    out.push_str(&format!("- Polynomials: {}\n", items.len()));
    out.push_str(&format!(
        "- Required checks: {}\n",
        if all_required_checks_passed {
            "passed"
        } else {
            "failed"
        }
    ));
    if let Some(first_failure) = first_failure {
        out.push_str(&format!("- First failure: {first_failure}\n"));
    }

    out.push_str("\n## Polynomial properties\n\n");
    out.push_str("| i | polynomial | degree | RR | gamma+ | LC | ULC | pal |\n");
    out.push_str("|---:|---|---:|:---:|:---:|:---:|:---:|:---:|\n");
    for item in items {
        out.push_str(&format!(
            "| {} | `{}` | {} | {} | {} | {} | {} | {} |\n",
            item.index,
            item.polynomial.polynomial,
            item.polynomial.degree,
            yes_no(item.real_rooted),
            yes_no(item.gamma_positive),
            yes_no(item.log_concave),
            yes_no(item.ultra_log_concave),
            yes_no(item.palindromic)
        ));
    }

    if !pairs.is_empty() {
        out.push_str("\n## Consecutive interlacing\n\n");
        out.push_str("| pair | weak | strict | status |\n");
        out.push_str("|---:|:---:|:---:|---|\n");
        for pair in pairs {
            if let Some(result) = &pair.result {
                out.push_str(&format!(
                    "| {}-{} | {} | {} | `{}` |\n",
                    pair.left_index,
                    pair.right_index,
                    opt_yes_no(result.weak),
                    opt_yes_no(result.strict),
                    result.status
                ));
            } else {
                out.push_str(&format!(
                    "| {}-{} | n/a | n/a | `{}` |\n",
                    pair.left_index,
                    pair.right_index,
                    pair.error.as_deref().unwrap_or("parse error")
                ));
            }
        }
    }

    if let Some(lace) = lace {
        out.push_str("\n## Finite Lace check\n\n");
        out.push_str(&format!(
            "- Truncation: block rows `{}`, block cols `{}`; matrix size `{} x {}`.\n",
            lace.block_rows, lace.block_cols, lace.rows, lace.columns
        ));
        out.push_str(&format!(
            "- Checked minors up to size `{}`: {}.\n",
            lace.max_minor_size,
            if lace.tnn { "TNN" } else { "not TNN" }
        ));
        if let Some(error) = &lace.error {
            out.push_str(&format!("- Failure: `{error}`\n"));
        }
    }

    if let Some(recurrence) = recurrence {
        out.push_str("\n## Recurrence search\n\n");
        if let Some(found) = &recurrence.recurrence {
            out.push_str(&format!("- Found: `{found}`\n"));
        } else if let Some(error) = &recurrence.error {
            out.push_str(&format!("- Not found: {error}\n"));
        } else {
            out.push_str("- Not found.\n");
        }
    }

    out
}

fn displayed(coefficients: &[i64]) -> DisplayedPolynomial {
    DisplayedPolynomial {
        polynomial: format_poly(coefficients),
        coefficients: coefficients.to_vec(),
    }
}

fn format_rational<T: ToString>(rational: T) -> String {
    rational.to_string()
}

fn parse_rational(input: &str) -> Result<BigRational, String> {
    input
        .parse::<BigRational>()
        .map_err(|e| format!("invalid rational `{input}`: {e}"))
}

fn apply_recurrence_options(input: Option<RecurrenceSearchOptionsInput>) -> AdaptiveSearchOptions {
    let mut options = AdaptiveSearchOptions::default();
    let Some(input) = input else {
        options.verbose = false;
        return options;
    };

    if let Some(value) = input.skip_prefix {
        options.skip_prefix = value;
    }
    if let Some(value) = input.require_all_offsets {
        options.require_all_offsets = value;
    }
    if let Some(value) = input.min_rec_len {
        options.min_rec_len = value;
    }
    if let Some(value) = input.max_rec_len {
        options.max_rec_len = value;
    }
    if let Some(value) = input.min_var_deg {
        options.min_var_deg = value;
    }
    if let Some(value) = input.max_var_deg {
        options.max_var_deg = value;
    }
    if let Some(value) = input.min_idx_deg {
        options.min_idx_deg = value;
    }
    if let Some(value) = input.max_idx_deg {
        options.max_idx_deg = value;
    }
    if let Some(value) = input.min_diff_deg {
        options.min_diff_deg = value;
    }
    if let Some(value) = input.max_diff_deg {
        options.max_diff_deg = value;
    }
    if let Some(value) = input.try_inhomogeneous {
        options.try_inhomogeneous = value;
    }
    if let Some(value) = input.min_inhomo_var_deg {
        options.try_inhomogeneous = true;
        options.min_inhomo_var_deg = value;
    }
    if let Some(value) = input.max_inhomo_var_deg {
        options.try_inhomogeneous = true;
        options.max_inhomo_var_deg = value;
    }
    if let Some(value) = input.min_inhomo_idx_deg {
        options.try_inhomogeneous = true;
        options.min_inhomo_idx_deg = value;
    }
    if let Some(value) = input.max_inhomo_idx_deg {
        options.try_inhomogeneous = true;
        options.max_inhomo_idx_deg = value;
    }
    if let Some(value) = input.try_denominator {
        options.try_denominator = value;
    }
    if let Some(value) = input.try_alternating_sign {
        options.try_alternating_sign = value;
    }
    if let Some(value) = input.max_denom_var_deg {
        options.try_denominator = true;
        options.max_denom_var_deg = value;
    }
    if let Some(value) = input.max_denom_idx_deg {
        options.try_denominator = true;
        options.max_denom_idx_deg = value;
    }
    if let Some(value) = input.min_margin {
        options.min_margin = value;
    }
    options.verbose = false;
    options
}

#[tool_router(router = tool_router)]
impl PolynomialToolsServer {
    pub fn new() -> Self {
        Self {
            tool_router: Self::tool_router(),
        }
    }

    #[tool(description = "Parse and format one or more dense univariate polynomials.")]
    pub fn parse_polynomials(
        &self,
        Parameters(input): Parameters<PolynomialBatchInput>,
    ) -> Result<Json<ParsePolynomialsResponse>, McpError> {
        let batch = parse_batch(&input)?;
        Ok(Json(ParsePolynomialsResponse {
            items: parse_items(&batch),
        }))
    }

    #[tool(
        description = "Check a polynomial family in one exact batch: properties, consecutive interlacing, optional finite Lace TNN, and optional recurrence search."
    )]
    pub fn check_polynomial_family(
        &self,
        Parameters(input): Parameters<CheckPolynomialFamilyRequest>,
    ) -> Result<Json<CheckPolynomialFamilyResponse>, McpError> {
        let (source, polynomials) = family_source_and_polynomials(&input)?;
        let polynomials = match polynomials {
            Ok(polynomials) => polynomials,
            Err(parse_errors) => {
                let first_failure = Some("one or more polynomials failed to parse".to_string());
                let mut markdown = String::new();
                markdown.push_str("# Polynomial family check\n\n");
                markdown.push_str(&format!("- Source: `{source}`\n"));
                markdown.push_str("- Required checks: failed\n");
                markdown.push_str("- First failure: one or more polynomials failed to parse\n\n");
                markdown.push_str("## Parse errors\n\n");
                for item in &parse_errors {
                    markdown.push_str(&format!(
                        "- polynomial {}: {}\n",
                        item.index,
                        item.error.as_deref().unwrap_or("parse error")
                    ));
                }
                return Ok(Json(CheckPolynomialFamilyResponse {
                    source,
                    item_count: 0,
                    all_required_checks_passed: false,
                    first_failure,
                    parse_errors,
                    items: Vec::new(),
                    consecutive_pairs: Vec::new(),
                    lace: None,
                    recurrence: None,
                    markdown,
                }));
            }
        };

        let options = input.options.as_ref();
        let check_consecutive_interlacing = options
            .and_then(|o| o.check_consecutive_interlacing)
            .unwrap_or(true);
        let find_recurrence = options.and_then(|o| o.find_recurrence).unwrap_or(false);

        let items: Vec<_> = polynomials
            .iter()
            .cloned()
            .enumerate()
            .map(|(index, polynomial)| family_report(index, polynomial))
            .collect();

        let consecutive_pairs = if check_consecutive_interlacing {
            polynomials
                .windows(2)
                .enumerate()
                .map(|(pair_index, pair)| InterlacingPairItem {
                    pair_index,
                    left_index: pair_index,
                    right_index: pair_index + 1,
                    ok: true,
                    result: Some(interlacing_result(pair[0].clone(), pair[1].clone())),
                    error: None,
                })
                .collect()
        } else {
            Vec::new()
        };

        let lace = input
            .lace
            .as_ref()
            .map(|request| check_family_lace(&polynomials, request));

        let coefficients: Vec<Vec<i64>> =
            polynomials.iter().map(|p| p.coefficients.clone()).collect();
        let recurrence = find_recurrence
            .then(|| family_recurrence_response(&coefficients, input.recurrence_options.clone()));

        let first_failure = first_family_failure(
            &items,
            &consecutive_pairs,
            lace.as_ref(),
            recurrence.as_ref(),
            options,
        );
        let all_required_checks_passed = first_failure.is_none();
        let markdown = family_markdown(
            &source,
            all_required_checks_passed,
            first_failure.as_deref(),
            &items,
            &consecutive_pairs,
            lace.as_ref(),
            recurrence.as_ref(),
        );

        Ok(Json(CheckPolynomialFamilyResponse {
            source,
            item_count: items.len(),
            all_required_checks_passed,
            first_failure,
            parse_errors: Vec::new(),
            items,
            consecutive_pairs,
            lace,
            recurrence,
            markdown,
        }))
    }

    #[tool(
        description = "Check real-rootedness, gamma-positivity, log-concavity, and related polynomial properties."
    )]
    pub fn polynomial_properties(
        &self,
        Parameters(input): Parameters<PolynomialBatchInput>,
    ) -> Result<Json<PolynomialPropertiesResponse>, McpError> {
        let batch = parse_batch(&input)?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| match item {
                Ok(polynomial) => {
                    let coefficients = polynomial.coefficients.clone();
                    PolynomialPropertiesItem {
                        index,
                        ok: true,
                        error: None,
                        polynomial: Some(polynomial),
                        real_rooted: Some(is_real_rooted(&coefficients)),
                        simple_roots: Some(has_simple_roots(&coefficients)),
                        palindromic: Some(is_palindromic_ignoring_initial_zeros(&coefficients)),
                        gamma_positive: Some(is_gamma_positive_ignoring_initial_zeros(
                            &coefficients,
                        )),
                        gamma_coefficients: gamma_coefficients_ignoring_initial_zeros(
                            &coefficients,
                        ),
                        log_concave: Some(is_log_concave(&coefficients)),
                        ultra_log_concave: Some(is_ultra_log_concave(&coefficients)),
                    }
                }
                Err(error) => PolynomialPropertiesItem {
                    index,
                    ok: false,
                    error: Some(error),
                    polynomial: None,
                    real_rooted: None,
                    simple_roots: None,
                    palindromic: None,
                    gamma_positive: None,
                    gamma_coefficients: None,
                    log_concave: None,
                    ultra_log_concave: None,
                },
            })
            .collect();
        Ok(Json(PolynomialPropertiesResponse { items }))
    }

    #[tool(description = "Check strict and weak interlacing for a pair of polynomials.")]
    pub fn check_interlacing_pair(
        &self,
        Parameters(input): Parameters<InterlacingPairRequest>,
    ) -> Result<Json<InterlacingPairResult>, McpError> {
        let p = normalize_polynomial(parse_required_polynomial(&input.p, "p")?);
        let q = normalize_polynomial(parse_required_polynomial(&input.q, "q")?);
        Ok(Json(interlacing_result(p, q)))
    }

    #[tool(description = "Check strict and weak interlacing for consecutive polynomial pairs.")]
    pub fn check_interlacing_sequence(
        &self,
        Parameters(input): Parameters<PolynomialBatchInput>,
    ) -> Result<Json<InterlacingSequenceResponse>, McpError> {
        let batch = parse_batch(&input)?;
        let items = parse_items(&batch);
        let pairs = batch
            .windows(2)
            .enumerate()
            .map(|(pair_index, pair)| {
                let left_index = pair_index;
                let right_index = pair_index + 1;
                match (&pair[0], &pair[1]) {
                    (Ok(p), Ok(q)) => InterlacingPairItem {
                        pair_index,
                        left_index,
                        right_index,
                        ok: true,
                        result: Some(interlacing_result(p.clone(), q.clone())),
                        error: None,
                    },
                    (left, right) => {
                        let mut messages = Vec::new();
                        if let Err(error) = left {
                            messages.push(format!("polynomial {left_index}: {error}"));
                        }
                        if let Err(error) = right {
                            messages.push(format!("polynomial {right_index}: {error}"));
                        }
                        InterlacingPairItem {
                            pair_index,
                            left_index,
                            right_index,
                            ok: false,
                            result: None,
                            error: Some(messages.join("; ")),
                        }
                    }
                }
            })
            .collect();
        Ok(Json(InterlacingSequenceResponse { items, pairs }))
    }

    #[tool(
        description = "Return rational midpoint representatives for isolated real roots, or real_rooted=false."
    )]
    pub fn real_roots(
        &self,
        Parameters(input): Parameters<PolynomialBatchInput>,
    ) -> Result<Json<RealRootsResponse>, McpError> {
        let batch = parse_batch(&input)?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| match item {
                Ok(polynomial) => match polynomial_tools::real_roots(&polynomial.coefficients) {
                    Some(roots) => RealRootsItem {
                        index,
                        ok: true,
                        error: None,
                        polynomial: Some(polynomial),
                        real_rooted: Some(true),
                        roots: Some(roots.into_iter().map(format_rational).collect()),
                    },
                    None => RealRootsItem {
                        index,
                        ok: true,
                        error: None,
                        polynomial: Some(polynomial),
                        real_rooted: Some(false),
                        roots: None,
                    },
                },
                Err(error) => RealRootsItem {
                    index,
                    ok: false,
                    error: Some(error),
                    polynomial: None,
                    real_rooted: None,
                    roots: None,
                },
            })
            .collect();
        Ok(Json(RealRootsResponse { items }))
    }

    #[tool(
        description = "Search adaptively for a polynomial recurrence. Provide exactly one of `polynomials`, `coefficients`, `expressions`, or `text`."
    )]
    pub fn find_recurrence(
        &self,
        Parameters(input): Parameters<FindRecurrenceRequest>,
    ) -> Result<Json<FindRecurrenceResponse>, McpError> {
        let batch = parse_recurrence_batch_rational(&input)?;
        let polynomials = match collect_rational_polynomials_or_errors(batch) {
            Ok(polynomials) => polynomials,
            Err(parse_errors) => {
                return Ok(Json(FindRecurrenceResponse {
                    found: false,
                    recurrence: None,
                    latex: None,
                    mathematica: None,
                    sage: None,
                    unknowns: None,
                    equations: None,
                    candidates_tried: None,
                    error: Some("one or more polynomials failed to parse".to_string()),
                    parse_errors,
                }));
            }
        };
        if polynomials.len() < 3 {
            return Ok(Json(FindRecurrenceResponse {
                found: false,
                recurrence: None,
                latex: None,
                mathematica: None,
                sage: None,
                unknowns: None,
                equations: None,
                candidates_tried: None,
                error: Some("need at least 3 polynomials".to_string()),
                parse_errors: Vec::new(),
            }));
        }
        let search = apply_recurrence_options(input.options);
        match find_recurrence_adaptive_rational(&polynomials, &search) {
            Some(result) => Ok(Json(FindRecurrenceResponse {
                found: true,
                recurrence: Some(format!("{}", result.recurrence)),
                latex: Some(result.recurrence.to_latex()),
                mathematica: Some(
                    result
                        .recurrence
                        .to_mathematica_definition_rational(&polynomials),
                ),
                sage: Some(result.recurrence.to_sage_definition_rational(&polynomials)),
                unknowns: Some(result.num_unknowns),
                equations: Some(result.num_equations),
                candidates_tried: Some(result.candidates_tried),
                error: None,
                parse_errors: Vec::new(),
            })),
            None => Ok(Json(FindRecurrenceResponse {
                found: false,
                recurrence: None,
                latex: None,
                mathematica: None,
                sage: None,
                unknowns: None,
                equations: None,
                candidates_tried: None,
                error: Some("no recurrence found within the search bounds".to_string()),
                parse_errors: Vec::new(),
            })),
        }
    }

    #[tool(description = "Compute the exact resultant of two polynomials.")]
    pub fn resultant(
        &self,
        Parameters(input): Parameters<InterlacingPairRequest>,
    ) -> Result<Json<ResultantResponse>, McpError> {
        let p = normalize_polynomial(parse_required_polynomial(&input.p, "p")?);
        let q = normalize_polynomial(parse_required_polynomial(&input.q, "q")?);
        let resultant = polynomial_tools::resultant(&p.coefficients, &q.coefficients).to_string();
        Ok(Json(ResultantResponse { p, q, resultant }))
    }

    #[tool(description = "Compute exact discriminants for one or more polynomials.")]
    pub fn discriminant(
        &self,
        Parameters(input): Parameters<PolynomialBatchInput>,
    ) -> Result<Json<DiscriminantResponse>, McpError> {
        let batch = parse_batch(&input)?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| match item {
                Ok(polynomial) => DiscriminantItem {
                    index,
                    ok: true,
                    error: None,
                    discriminant: Some(
                        polynomial_tools::discriminant(&polynomial.coefficients).to_string(),
                    ),
                    polynomial: Some(polynomial),
                },
                Err(error) => DiscriminantItem {
                    index,
                    ok: false,
                    error: Some(error),
                    polynomial: None,
                    discriminant: None,
                },
            })
            .collect();
        Ok(Json(DiscriminantResponse { items }))
    }

    #[tool(
        description = "Convert h*-vectors to Ehrhart polynomials or Ehrhart polynomials to h*-vectors."
    )]
    pub fn ehrhart_hstar(
        &self,
        Parameters(input): Parameters<EhrhartHstarRequest>,
    ) -> Result<Json<EhrhartHstarResponse>, McpError> {
        match input.mode {
            EhrhartHstarMode::HstarToEhrhart => {
                let hstar = input
                    .hstar
                    .ok_or_else(|| invalid_params("`hstar` is required for hstar_to_ehrhart"))?;
                let ehrhart = hstar_to_ehrhart(&hstar)
                    .into_iter()
                    .map(format_rational)
                    .collect();
                Ok(Json(EhrhartHstarResponse {
                    mode: EhrhartHstarMode::HstarToEhrhart,
                    hstar: Some(hstar),
                    ehrhart_coefficients: Some(ehrhart),
                }))
            }
            EhrhartHstarMode::EhrhartToHstar => {
                let hstar = match (
                    input.ehrhart_coefficients,
                    input.numerator_coefficients,
                    input.denominator,
                ) {
                    (Some(coefficients), None, None) => {
                        let coefficients: Result<Vec<_>, _> =
                            coefficients.iter().map(|c| parse_rational(c)).collect();
                        ehrhart_to_hstar(&coefficients.map_err(invalid_params)?)
                    }
                    (None, Some(numerator_coefficients), Some(denominator)) => {
                        if denominator == 0 {
                            return Err(invalid_params("`denominator` must be nonzero"));
                        }
                        ehrhart_to_hstar_with_denom(&numerator_coefficients, denominator)
                    }
                    _ => {
                        return Err(invalid_params(
                            "for ehrhart_to_hstar, provide either `ehrhart_coefficients` or both `numerator_coefficients` and `denominator`",
                        ));
                    }
                };
                Ok(Json(EhrhartHstarResponse {
                    mode: EhrhartHstarMode::EhrhartToHstar,
                    hstar: Some(hstar),
                    ehrhart_coefficients: None,
                }))
            }
        }
    }

    #[tool(
        description = "Analyze symmetric decomposition, R-transform, interlacing checks, and magic-basis coordinates."
    )]
    pub fn analyze_decomposition(
        &self,
        Parameters(input): Parameters<PolynomialBatchInput>,
    ) -> Result<Json<DecompositionResponse>, McpError> {
        let batch = parse_batch(&input)?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| match item {
                Ok(polynomial) => {
                    match analyze_symmetric_decomposition_i64(&polynomial.coefficients) {
                        Ok(analysis) => {
                            let partial_sum_checks = analysis
                                .magic
                                .left_partial_sums
                                .iter()
                                .zip(analysis.magic.right_partial_sums.iter())
                                .map(|(left, right)| left <= right)
                                .collect();
                            DecompositionItem {
                                index,
                                ok: true,
                                report: Some(DecompositionReport {
                                    polynomial,
                                    reciprocal: displayed(&analysis.reciprocal),
                                    a: displayed(&analysis.a),
                                    b: displayed(&analysis.b),
                                    a_real_rooted: analysis.a_real_rooted,
                                    b_real_rooted: analysis.b_real_rooted,
                                    b_interlaces_a: analysis.b_interlaces_a,
                                    reciprocal_interlaces_input: analysis
                                        .reciprocal_interlaces_input,
                                    alternatingly_increasing: analysis.alternatingly_increasing,
                                    f_polynomial: displayed(&analysis.f_polynomial),
                                    r_transform_of_f: displayed(&analysis.r_transform_of_f),
                                    r_a: displayed(&analysis.r_a),
                                    r_b: displayed(&analysis.r_b),
                                    r_interlaces_f: analysis.r_interlaces_f,
                                    magic: MagicBasisReport {
                                        coordinates: analysis
                                            .magic
                                            .coordinates
                                            .iter()
                                            .map(ToString::to_string)
                                            .collect(),
                                        left_partial_sums: analysis
                                            .magic
                                            .left_partial_sums
                                            .iter()
                                            .map(ToString::to_string)
                                            .collect(),
                                        right_partial_sums: analysis
                                            .magic
                                            .right_partial_sums
                                            .iter()
                                            .map(ToString::to_string)
                                            .collect(),
                                        partial_sum_checks,
                                        all_nonnegative: analysis.magic.all_nonnegative,
                                        left_leq_right: analysis.magic.left_leq_right,
                                    },
                                }),
                                error: None,
                            }
                        }
                        Err(error) => DecompositionItem {
                            index,
                            ok: false,
                            report: None,
                            error: Some(error.to_string()),
                        },
                    }
                }
                Err(error) => DecompositionItem {
                    index,
                    ok: false,
                    report: None,
                    error: Some(error),
                },
            })
            .collect();
        Ok(Json(DecompositionResponse { items }))
    }

    #[tool(description = "Generate standard polynomial sequences.")]
    pub fn generate_sequence(
        &self,
        Parameters(input): Parameters<GenerateSequenceRequest>,
    ) -> Result<Json<GenerateSequenceResponse>, McpError> {
        let polynomials = generated_sequence_polynomials(&input.sequence, input.max_n)
            .into_iter()
            .map(normalize_polynomial)
            .collect();

        Ok(Json(GenerateSequenceResponse {
            sequence: input.sequence,
            max_n: input.max_n,
            polynomials,
        }))
    }
}

#[tool_handler(router = self.tool_router)]
impl ServerHandler for PolynomialToolsServer {
    fn get_info(&self) -> ServerInfo {
        ServerInfo {
            protocol_version: Default::default(),
            capabilities: ServerCapabilities::builder().enable_tools().build(),
            server_info: Implementation {
                name: "polynomial-tools".to_string(),
                title: Some("polynomial-tools".to_string()),
                version: env!("CARGO_PKG_VERSION").to_string(),
                description: Some(
                    "Exact univariate polynomial tools for combinatorial research.".to_string(),
                ),
                icons: None,
                website_url: None,
            },
            instructions: Some(
                "Exact univariate polynomial tools for combinatorial research.".to_string(),
            ),
        }
    }
}

impl Default for PolynomialToolsServer {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn coeffs(values: &[i64]) -> PolynomialInput {
        PolynomialInput {
            coefficients: Some(values.to_vec()),
            expression: None,
        }
    }

    fn expr(value: &str) -> PolynomialInput {
        PolynomialInput {
            coefficients: None,
            expression: Some(value.to_string()),
        }
    }

    #[test]
    fn parses_coefficients_and_expression() {
        assert_eq!(
            normalize_polynomial(parse_polynomial_input(&coeffs(&[1, 2, 3, 0])).unwrap()),
            NormalizedPolynomial {
                polynomial: "1 + 2t + 3t^2".to_string(),
                coefficients: vec![1, 2, 3],
                degree: 2,
            }
        );
        assert_eq!(
            normalize_polynomial(parse_polynomial_input(&expr("3t^2 + 2t + 1")).unwrap())
                .coefficients,
            vec![1, 2, 3]
        );
    }

    #[test]
    fn preserves_batch_parse_errors() {
        let batch = parse_batch(&PolynomialBatchInput {
            polynomials: None,
            text: Some("1,2,1\nbad + x\n".to_string()),
        })
        .unwrap();
        let items = parse_items(&batch);
        assert_eq!(items.len(), 2);
        assert!(items[0].ok);
        assert!(!items[1].ok);
    }

    #[test]
    fn checks_eulerian_properties() {
        let server = PolynomialToolsServer::new();
        let Json(response) = server
            .polynomial_properties(Parameters(PolynomialBatchInput {
                polynomials: Some(vec![coeffs(&[1, 11, 11, 1])]),
                text: None,
            }))
            .unwrap();
        let item = &response.items[0];
        assert_eq!(item.real_rooted, Some(true));
        assert_eq!(item.gamma_positive, Some(true));
        assert_eq!(item.gamma_coefficients, Some(vec![1, 8]));
    }

    #[test]
    fn checks_strict_and_weak_interlacing() {
        let server = PolynomialToolsServer::new();
        let Json(strict) = server
            .check_interlacing_pair(Parameters(InterlacingPairRequest {
                p: coeffs(&[-2, 1]),
                q: coeffs(&[3, -4, 1]),
            }))
            .unwrap();
        assert_eq!(strict.strict, Some(true));

        let Json(weak) = server
            .check_interlacing_pair(Parameters(InterlacingPairRequest {
                p: coeffs(&[-1, 1]),
                q: coeffs(&[2, -3, 1]),
            }))
            .unwrap();
        assert_eq!(weak.weak, Some(true));
    }

    #[test]
    fn computes_resultant_and_discriminant() {
        let server = PolynomialToolsServer::new();
        let Json(resultant) = server
            .resultant(Parameters(InterlacingPairRequest {
                p: coeffs(&[2, -3, 1]),
                q: coeffs(&[-3, 1]),
            }))
            .unwrap();
        assert_eq!(resultant.resultant, "2");

        let Json(discriminant) = server
            .discriminant(Parameters(PolynomialBatchInput {
                polynomials: Some(vec![coeffs(&[-1, 0, 1])]),
                text: None,
            }))
            .unwrap();
        assert_eq!(discriminant.items[0].discriminant.as_deref(), Some("4"));
    }

    #[test]
    fn roundtrips_hstar_and_ehrhart() {
        let server = PolynomialToolsServer::new();
        let Json(to_ehrhart) = server
            .ehrhart_hstar(Parameters(EhrhartHstarRequest {
                mode: EhrhartHstarMode::HstarToEhrhart,
                hstar: Some(vec![1, 1, 0]),
                ehrhart_coefficients: None,
                numerator_coefficients: None,
                denominator: None,
            }))
            .unwrap();
        assert_eq!(
            to_ehrhart.ehrhart_coefficients.as_ref().unwrap(),
            &vec!["1".to_string(), "2".to_string(), "1".to_string()]
        );

        let Json(to_hstar) = server
            .ehrhart_hstar(Parameters(EhrhartHstarRequest {
                mode: EhrhartHstarMode::EhrhartToHstar,
                hstar: None,
                ehrhart_coefficients: Some(vec!["1".to_string(), "2".to_string(), "1".to_string()]),
                numerator_coefficients: None,
                denominator: None,
            }))
            .unwrap();
        assert_eq!(to_hstar.hstar, Some(vec![1, 1, 0]));
    }

    #[test]
    fn recurrence_search_examples() {
        let server = PolynomialToolsServer::new();
        let Json(geometric) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: Some(vec![
                    coeffs(&[1]),
                    coeffs(&[2]),
                    coeffs(&[4]),
                    coeffs(&[8]),
                    coeffs(&[16]),
                ]),
                coefficients: None,
                expressions: None,
                text: None,
                options: None,
            }))
            .unwrap();
        assert_eq!(geometric.recurrence.as_deref(), Some("P(n) = 2 P(n-1)"));

        let Json(fibonacci) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: None,
                coefficients: Some(vec![
                    vec![1.into()],
                    vec![1.into()],
                    vec![2.into()],
                    vec![3.into()],
                    vec![5.into()],
                    vec![8.into()],
                ]),
                expressions: None,
                text: None,
                options: None,
            }))
            .unwrap();
        assert_eq!(
            fibonacci.recurrence.as_deref(),
            Some("P(n) = P(n-1) + P(n-2)")
        );

        let Json(expression_geometric) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: None,
                coefficients: None,
                expressions: Some(vec![
                    "1".to_string(),
                    "2".to_string(),
                    "4".to_string(),
                    "8".to_string(),
                    "16".to_string(),
                ]),
                text: None,
                options: None,
            }))
            .unwrap();
        assert_eq!(
            expression_geometric.recurrence.as_deref(),
            Some("P(n) = 2 P(n-1)")
        );

        let Json(rational_geometric) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: None,
                coefficients: Some(vec![
                    vec![RationalCoefficientInput::Text("1/2".to_string())],
                    vec![RationalCoefficientInput::Text("1/4".to_string())],
                    vec![RationalCoefficientInput::Text("1/8".to_string())],
                    vec![RationalCoefficientInput::Text("1/16".to_string())],
                    vec![RationalCoefficientInput::Text("1/32".to_string())],
                ]),
                expressions: None,
                text: None,
                options: None,
            }))
            .unwrap();
        assert_eq!(
            rational_geometric.recurrence.as_deref(),
            Some("P(n) = 1/2 P(n-1)")
        );

        let large = "1267650600228229401496703205376";
        let double_large = "2535301200456458802993406410752";
        let Json(large_geometric) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: None,
                coefficients: Some(vec![
                    vec![RationalCoefficientInput::Text(large.to_string())],
                    vec![RationalCoefficientInput::Text(double_large.to_string())],
                    vec![RationalCoefficientInput::Text(
                        "5070602400912917605986812821504".to_string(),
                    )],
                    vec![RationalCoefficientInput::Text(
                        "10141204801825835211973625643008".to_string(),
                    )],
                    vec![RationalCoefficientInput::Text(
                        "20282409603651670423947251286016".to_string(),
                    )],
                ]),
                expressions: None,
                text: None,
                options: None,
            }))
            .unwrap();
        assert_eq!(
            large_geometric.recurrence.as_deref(),
            Some("P(n) = 2 P(n-1)")
        );
        assert!(large_geometric
            .sage
            .as_deref()
            .is_some_and(|sage| sage.contains(large)));

        let Json(alternating) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: None,
                coefficients: Some(vec![
                    vec![1.into()],
                    vec![1.into()],
                    vec![(-1).into()],
                    vec![(-1).into()],
                    vec![1.into()],
                    vec![1.into()],
                    vec![(-1).into()],
                    vec![(-1).into()],
                ]),
                expressions: None,
                text: None,
                options: Some(
                    serde_json::from_value(json!({
                        "try_alternating_sign": true,
                        "max_rec_len": 1,
                        "max_var_deg": 0,
                        "max_idx_deg": 0,
                        "max_diff_deg": 0
                    }))
                    .unwrap(),
                ),
            }))
            .unwrap();
        assert_eq!(
            alternating.recurrence.as_deref(),
            Some("P(n) = (-1)^n P(n-1)")
        );
        assert!(alternating
            .sage
            .as_deref()
            .is_some_and(|sage| sage.contains("(-1)**n*P(n - 1)")));

        let eulerian = vec![
            coeffs(&[1]),
            coeffs(&[1]),
            coeffs(&[1, 1]),
            coeffs(&[1, 4, 1]),
            coeffs(&[1, 11, 11, 1]),
            coeffs(&[1, 26, 66, 26, 1]),
        ];
        let Json(eulerian_result) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: Some(eulerian),
                coefficients: None,
                expressions: None,
                text: None,
                options: None,
            }))
            .unwrap();
        assert!(eulerian_result.found);

        let conflict = parse_recurrence_batch_rational(&FindRecurrenceRequest {
            polynomials: Some(vec![coeffs(&[1]), coeffs(&[2]), coeffs(&[4])]),
            coefficients: Some(vec![vec![1.into()], vec![2.into()], vec![4.into()]]),
            expressions: None,
            text: None,
            options: None,
        })
        .unwrap_err();
        assert!(format!("{conflict:?}").contains("expected exactly one"));
    }

    #[test]
    fn family_check_reports_properties_and_recurrence() {
        let server = PolynomialToolsServer::new();
        let Json(response) = server
            .check_polynomial_family(Parameters(CheckPolynomialFamilyRequest {
                polynomials: Some(vec![
                    coeffs(&[1]),
                    coeffs(&[2]),
                    coeffs(&[4]),
                    coeffs(&[8]),
                    coeffs(&[16]),
                ]),
                text: None,
                sequence: None,
                max_n: None,
                options: Some(FamilyCheckOptions {
                    require_real_rooted: Some(true),
                    require_simple_roots: None,
                    require_palindromic: None,
                    require_gamma_positive: None,
                    require_log_concave: None,
                    require_ultra_log_concave: None,
                    check_consecutive_interlacing: Some(false),
                    require_consecutive_weak_interlacing: None,
                    find_recurrence: Some(true),
                }),
                recurrence_options: None,
                lace: None,
            }))
            .unwrap();

        assert!(response.all_required_checks_passed);
        assert_eq!(response.item_count, 5);
        assert_eq!(
            response
                .recurrence
                .as_ref()
                .and_then(|r| r.recurrence.as_deref()),
            Some("P(n) = 2 P(n-1)")
        );
        assert!(response.markdown.contains("Polynomial properties"));
    }

    #[test]
    fn family_check_can_require_finite_lace_tnn() {
        let server = PolynomialToolsServer::new();
        let Json(response) = server
            .check_polynomial_family(Parameters(CheckPolynomialFamilyRequest {
                polynomials: Some(vec![
                    coeffs(&[2, 1]),
                    coeffs(&[8, 6, 1]),
                    coeffs(&[3, 4, 1]),
                ]),
                text: None,
                sequence: None,
                max_n: None,
                options: None,
                recurrence_options: None,
                lace: Some(LaceCheckRequest {
                    block_rows: Some(1),
                    block_cols: Some(3),
                    max_minor_size: Some(3),
                    include_matrix: Some(true),
                }),
            }))
            .unwrap();

        assert!(!response.all_required_checks_passed);
        assert_eq!(
            response.first_failure.as_deref(),
            Some("finite Lace truncation is not TNN up to minors of size 3")
        );
        assert_eq!(response.lace.as_ref().map(|l| l.tnn), Some(false));
        assert!(response
            .lace
            .as_ref()
            .and_then(|l| l.error.as_deref())
            .is_some_and(|error| error.contains("det = -1")));
    }

    #[test]
    fn serializes_bigint_and_bigrational_as_strings() {
        let server = PolynomialToolsServer::new();
        let Json(resultant) = server
            .resultant(Parameters(InterlacingPairRequest {
                p: coeffs(&[1, 0, 1]),
                q: coeffs(&[-1, 1]),
            }))
            .unwrap();
        let value = serde_json::to_value(resultant).unwrap();
        assert_eq!(value["resultant"], "2");

        let Json(ehrhart) = server
            .ehrhart_hstar(Parameters(EhrhartHstarRequest {
                mode: EhrhartHstarMode::HstarToEhrhart,
                hstar: Some(vec![1, 0, 0]),
                ehrhart_coefficients: None,
                numerator_coefficients: None,
                denominator: None,
            }))
            .unwrap();
        let value = serde_json::to_value(ehrhart).unwrap();
        assert_eq!(value["ehrhart_coefficients"][1], "3/2");
    }
}
