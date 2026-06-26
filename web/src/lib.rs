use polynomial_tools::recurrence::BigRational as RecurrenceBigRational;
use polynomial_tools::recurrence::*;
use polynomial_tools::*;
use serde::Serialize;
use wasm_bindgen::prelude::*;

// ---------------------------------------------------------------------------
// Result types (serialized to JSON)
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct PolyProps {
    polynomial: String,
    degree: usize,
    coefficients: Vec<i64>,
    real_rooted: bool,
    palindromic: bool,
    gamma_positive: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    gamma_coefficients: Option<Vec<i64>>,
    log_concave: bool,
    ultra_log_concave: bool,
}

#[derive(Serialize)]
struct InterlacingResult {
    p: String,
    q: String,
    strict: bool,
    weak: bool,
    status: String,
}

#[derive(Serialize)]
struct RecurrenceResult {
    found: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    recurrence: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    latex: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    mathematica: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    sage: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    python: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    recurrence_json: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    unknowns: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    weighted_unknowns: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    equations: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    fit_polynomials: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    verification_polynomials: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    candidates_tried: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<String>,
}

#[derive(Serialize)]
#[serde(untagged)]
enum ParsedOrError {
    Ok {
        polynomial: String,
        coefficients: Vec<i64>,
    },
    Err {
        error: String,
    },
}

#[derive(Serialize)]
struct DisplayedPoly {
    polynomial: String,
    coefficients: Vec<i64>,
}

#[derive(Serialize)]
struct MagicBasisReport {
    coordinates: Vec<String>,
    left_partial_sums: Vec<String>,
    right_partial_sums: Vec<String>,
    partial_sum_checks: Vec<bool>,
    all_nonnegative: bool,
    left_leq_right: bool,
}

#[derive(Serialize)]
struct DecompositionResult {
    polynomial: String,
    degree: usize,
    coefficients: Vec<i64>,
    reciprocal: DisplayedPoly,
    a: DisplayedPoly,
    b: DisplayedPoly,
    a_real_rooted: bool,
    b_real_rooted: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    b_interlaces_a: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    reciprocal_interlaces_input: Option<bool>,
    alternatingly_increasing: bool,
    f_polynomial: DisplayedPoly,
    r_transform_of_f: DisplayedPoly,
    r_a: DisplayedPoly,
    r_b: DisplayedPoly,
    #[serde(skip_serializing_if = "Option::is_none")]
    r_interlaces_f: Option<bool>,
    magic: MagicBasisReport,
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn strip_trailing_zeros(coeffs: &[i64]) -> &[i64] {
    let end = coeffs.iter().rposition(|&c| c != 0).map_or(0, |i| i + 1);
    if end == 0 {
        &[0i64; 0]
    } else {
        &coeffs[..end]
    }
}

fn integer_polys_to_rational(polys: &[Vec<i64>]) -> Vec<Vec<RecurrenceBigRational>> {
    polys
        .iter()
        .map(|row| {
            row.iter()
                .map(|coeff| {
                    parse_rational_coeff(&coeff.to_string())
                        .expect("integer coefficients parse as rationals")
                })
                .collect()
        })
        .collect()
}

fn recurrence_json_output(
    result: &AdaptiveSearchResult,
    polys: &[Vec<i64>],
    search: &AdaptiveSearchOptions,
) -> String {
    let rational_polys = integer_polys_to_rational(polys);
    let searched_polys = rational_polys.get(search.skip_prefix..).unwrap_or(&[]);
    let initial_count = result.recurrence.max_offset().min(searched_polys.len());
    let recurrence_json = RecurrenceJson::from_recurrence_rational(
        &result.recurrence,
        1,
        &searched_polys[..initial_count],
        Some(RecurrenceJsonSearch {
            recurrence_text: result.recurrence.to_string(),
            source_rows: polys.len(),
            skip_prefix: search.skip_prefix,
            unknowns: result.num_unknowns,
            weighted_unknowns: result.weighted_unknowns,
            equations: result.num_equations,
            fit_polynomials: result.fit_polynomials,
            verification_polynomials: result.verification_polynomials,
            candidates_tried: result.candidates_tried,
            options: RecurrenceOptionsJson::from(&result.opts),
        }),
    );
    serde_json::to_string_pretty(&recurrence_json).expect("serialize recurrence JSON")
}

fn parse_input(input: &str) -> Vec<Result<Vec<i64>, String>> {
    parse_polynomials(input)
}

fn parse_ok(input: &str) -> (Vec<Vec<i64>>, Vec<String>) {
    let mut polys = Vec::new();
    let mut errors = Vec::new();
    for r in parse_input(input) {
        match r {
            Ok(p) => polys.push(p),
            Err(e) => errors.push(e),
        }
    }
    (polys, errors)
}

fn display_poly(coeffs: &[i64]) -> DisplayedPoly {
    DisplayedPoly {
        polynomial: format_poly(coeffs),
        coefficients: coeffs.to_vec(),
    }
}

// ---------------------------------------------------------------------------
// WASM exports
// ---------------------------------------------------------------------------

#[wasm_bindgen]
pub fn check_properties(input: &str) -> String {
    let mut results: Vec<serde_json::Value> = Vec::new();

    for r in parse_input(input) {
        match r {
            Ok(coeffs) => {
                let c = strip_trailing_zeros(&coeffs);
                let deg = c.iter().rposition(|&x| x != 0).unwrap_or(0);
                let props = PolyProps {
                    polynomial: format_poly(c),
                    degree: deg,
                    coefficients: c.to_vec(),
                    real_rooted: is_real_rooted(c),
                    palindromic: is_palindromic_ignoring_initial_zeros(c),
                    gamma_positive: is_gamma_positive_ignoring_initial_zeros(c),
                    gamma_coefficients: gamma_coefficients_ignoring_initial_zeros(c),
                    log_concave: is_log_concave(c),
                    ultra_log_concave: is_ultra_log_concave(c),
                };
                results.push(serde_json::to_value(&props).unwrap());
            }
            Err(e) => {
                results.push(serde_json::json!({"error": e}));
            }
        }
    }

    serde_json::to_string(&results).unwrap()
}

#[wasm_bindgen]
pub fn check_interlacing_pairs(input: &str) -> String {
    let (polys, _) = parse_ok(input);
    let mut results: Vec<InterlacingResult> = Vec::new();

    for pair in polys.windows(2) {
        let p = strip_trailing_zeros(&pair[0]);
        let q = strip_trailing_zeros(&pair[1]);
        let strict = check_interlacing(p, q) == Some(true);
        let weak = check_weak_interlacing(p, q) == Some(true);
        let status = if strict {
            "strictly interlace".to_string()
        } else if weak {
            "weakly interlace (shared roots)".to_string()
        } else {
            "do NOT interlace".to_string()
        };
        results.push(InterlacingResult {
            p: format_poly(p),
            q: format_poly(q),
            strict,
            weak,
            status,
        });
    }

    serde_json::to_string(&results).unwrap()
}

#[wasm_bindgen]
pub fn compute_resultant(input: &str) -> String {
    let (polys, _) = parse_ok(input);
    if polys.len() < 2 {
        return serde_json::json!({"error": "need exactly two polynomials"}).to_string();
    }
    let r = resultant(&polys[0], &polys[1]);
    serde_json::json!({
        "p": format_poly(&polys[0]),
        "q": format_poly(&polys[1]),
        "resultant": r.to_string()
    })
    .to_string()
}

#[wasm_bindgen]
pub fn compute_discriminant(input: &str) -> String {
    let mut results: Vec<serde_json::Value> = Vec::new();
    for r in parse_input(input) {
        match r {
            Ok(coeffs) => {
                let c = strip_trailing_zeros(&coeffs);
                let d = discriminant(c);
                results.push(serde_json::json!({
                    "polynomial": format_poly(c),
                    "discriminant": d.to_string()
                }));
            }
            Err(e) => {
                results.push(serde_json::json!({"error": e}));
            }
        }
    }
    serde_json::to_string(&results).unwrap()
}

#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn find_recurrence(
    input: &str,
    max_rec_len: u32,
    max_var_deg: u32,
    max_idx_deg: u32,
    max_diff_deg: u32,
    try_inhomogeneous: bool,
    try_denominator: bool,
    try_alternating_sign: bool,
) -> String {
    let (polys, _) = parse_ok(input);
    if polys.len() < 3 {
        return serde_json::to_string(&RecurrenceResult {
            found: false,
            recurrence: None,
            latex: None,
            mathematica: None,
            sage: None,
            python: None,
            recurrence_json: None,
            unknowns: None,
            weighted_unknowns: None,
            equations: None,
            fit_polynomials: None,
            verification_polynomials: None,
            candidates_tried: None,
            error: Some("need at least 3 polynomials".to_string()),
        })
        .unwrap();
    }

    let search = AdaptiveSearchOptions {
        max_rec_len: max_rec_len as usize,
        max_var_deg: max_var_deg as usize,
        max_idx_deg: max_idx_deg as usize,
        max_diff_deg: max_diff_deg as usize,
        try_inhomogeneous,
        try_denominator,
        try_alternating_sign,
        ..Default::default()
    };

    match find_recurrence_adaptive(&polys, &search) {
        Some(res) => serde_json::to_string(&RecurrenceResult {
            found: true,
            recurrence: Some(format!("{}", res.recurrence)),
            latex: Some(res.recurrence.to_latex()),
            mathematica: Some(res.recurrence.to_mathematica_definition(&polys)),
            sage: Some(res.recurrence.to_sage_definition(&polys)),
            python: Some(res.recurrence.to_python_definition(&polys)),
            recurrence_json: Some(recurrence_json_output(&res, &polys, &search)),
            unknowns: Some(res.num_unknowns),
            weighted_unknowns: Some(res.weighted_unknowns),
            equations: Some(res.num_equations),
            fit_polynomials: Some(res.fit_polynomials),
            verification_polynomials: Some(res.verification_polynomials),
            candidates_tried: Some(res.candidates_tried),
            error: None,
        })
        .unwrap(),
        None => serde_json::to_string(&RecurrenceResult {
            found: false,
            recurrence: None,
            latex: None,
            mathematica: None,
            sage: None,
            python: None,
            recurrence_json: None,
            unknowns: None,
            weighted_unknowns: None,
            equations: None,
            fit_polynomials: None,
            verification_polynomials: None,
            candidates_tried: None,
            error: Some("no recurrence found within search bounds".to_string()),
        })
        .unwrap(),
    }
}

#[wasm_bindgen]
pub fn analyze_decompositions(input: &str) -> String {
    let mut results: Vec<serde_json::Value> = Vec::new();

    for r in parse_input(input) {
        match r {
            Ok(coeffs) => {
                let c = strip_trailing_zeros(&coeffs);
                match analyze_symmetric_decomposition_i64(c) {
                    Ok(analysis) => {
                        let report = DecompositionResult {
                            polynomial: format_poly(c),
                            degree: analysis.degree,
                            coefficients: c.to_vec(),
                            reciprocal: display_poly(&analysis.reciprocal),
                            a: display_poly(&analysis.a),
                            b: display_poly(&analysis.b),
                            a_real_rooted: analysis.a_real_rooted,
                            b_real_rooted: analysis.b_real_rooted,
                            b_interlaces_a: analysis.b_interlaces_a,
                            reciprocal_interlaces_input: analysis.reciprocal_interlaces_input,
                            alternatingly_increasing: analysis.alternatingly_increasing,
                            f_polynomial: display_poly(&analysis.f_polynomial),
                            r_transform_of_f: display_poly(&analysis.r_transform_of_f),
                            r_a: display_poly(&analysis.r_a),
                            r_b: display_poly(&analysis.r_b),
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
                                partial_sum_checks: analysis
                                    .magic
                                    .left_partial_sums
                                    .iter()
                                    .zip(analysis.magic.right_partial_sums.iter())
                                    .map(|(left, right)| left <= right)
                                    .collect(),
                                all_nonnegative: analysis.magic.all_nonnegative,
                                left_leq_right: analysis.magic.left_leq_right,
                            },
                        };
                        results.push(serde_json::to_value(&report).unwrap());
                    }
                    Err(e) => {
                        results.push(serde_json::json!({"error": e.to_string()}));
                    }
                }
            }
            Err(e) => {
                results.push(serde_json::json!({"error": e}));
            }
        }
    }

    serde_json::to_string(&results).unwrap()
}

/// Check interlacing between two polynomials given as JSON arrays of i64.
/// Returns JSON: {"strict": bool, "weak": bool}
#[wasm_bindgen]
pub fn check_interlacing_pair(p_json: &str, q_json: &str) -> String {
    let p: Vec<i64> = serde_json::from_str(p_json).unwrap_or_default();
    let q: Vec<i64> = serde_json::from_str(q_json).unwrap_or_default();
    let strict = check_interlacing(&p, &q) == Some(true);
    let weak = if strict {
        true
    } else {
        check_weak_interlacing(&p, &q) == Some(true)
    };
    serde_json::json!({"strict": strict, "weak": weak}).to_string()
}

#[wasm_bindgen]
pub fn parse_and_format(input: &str) -> String {
    let mut results: Vec<ParsedOrError> = Vec::new();
    for r in parse_input(input) {
        match r {
            Ok(coeffs) => {
                let c = strip_trailing_zeros(&coeffs);
                results.push(ParsedOrError::Ok {
                    polynomial: format_poly(c),
                    coefficients: c.to_vec(),
                });
            }
            Err(e) => {
                results.push(ParsedOrError::Err { error: e });
            }
        }
    }
    serde_json::to_string(&results).unwrap()
}
