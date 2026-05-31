use polynomial_lab::{
    default_family_registry, default_lab_root, interlacing_evaluation_draft,
    interlacing_evidence_id, real_rooted_evaluation_draft, real_rooted_evidence_id,
    CheckFamilyInterlacingReport, CheckFamilyRealRootednessReport, CheckedRange,
    ComputedPolynomial, ConjectureDraft, EvaluationDraft, EvaluationFilter, GeneratedFile,
    GoalDraft, ImplicationDraft, InterlacingMode, LabRecord, LabStore, PolynomialFamilyInfo,
    ProjectDraft, ProjectOverview, ProjectReport, TraceGoalSupport, ValidationMode,
    ValidationReport, WrittenEvaluation, WrittenRecord,
};
use rmcp::{
    handler::server::{router::tool::ToolRouter, wrapper::Parameters},
    model::{Implementation, ServerCapabilities, ServerInfo},
    tool, tool_handler, tool_router, ErrorData as McpError, Json, ServerHandler,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct PolynomialLabServer {
    root: PathBuf,
    tool_router: ToolRouter<Self>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EmptyRequest {}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ValidateRequest {
    pub strict: Option<bool>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ProjectRequest {
    pub project_id: String,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OptionalProjectRequest {
    pub project_id: Option<String>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EvaluationSearchRequest {
    pub project_id: Option<String>,
    pub status: Option<String>,
    pub relation_id: Option<String>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TraceGoalRequest {
    pub project_id: String,
    pub goal_id: String,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CreateProjectRequest {
    pub project_id: String,
    pub label: Option<String>,
    pub status: Option<String>,
    pub description: Option<String>,
    pub main_objects: Option<Vec<String>>,
    pub main_goals: Option<Vec<String>>,
    pub source_notes: Option<Vec<String>>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AppendGoalRequest {
    pub project_id: String,
    pub id: String,
    pub statement: String,
    pub label: Option<String>,
    pub status: Option<String>,
    pub objects: Option<Vec<String>>,
    pub motivation: Option<String>,
    pub current_best_route: Option<String>,
    pub depends_on: Option<Vec<String>>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AppendConjectureRequest {
    pub project_id: String,
    pub id: String,
    pub statement: String,
    pub label: Option<String>,
    pub status: Option<String>,
    pub relation: Option<String>,
    pub left: Option<String>,
    pub right: Option<String>,
    pub index_condition: Option<String>,
    pub depends_on: Option<Vec<String>>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AppendImplicationRequest {
    pub project_id: String,
    pub id: String,
    pub from: Vec<String>,
    pub to: String,
    pub label: Option<String>,
    pub status: Option<String>,
    pub explanation: Option<String>,
    pub proof_tags: Option<Vec<String>>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ComputeFamilyRequest {
    pub family_id: String,
    pub n: usize,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CheckFamilyRealRootedRequest {
    pub family_id: String,
    pub n_min: usize,
    pub n_max: usize,
    pub project_id: Option<String>,
    pub relation_id: Option<String>,
    pub id: Option<String>,
    pub append: Option<bool>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CheckFamilyInterlacingRequest {
    pub left_family_id: String,
    pub right_family_id: String,
    pub n_min: usize,
    pub n_max: usize,
    pub mode: Option<InterlacingMode>,
    pub project_id: Option<String>,
    pub relation_id: Option<String>,
    pub id: Option<String>,
    pub append: Option<bool>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GeneratedProjectRequest {
    pub project_id: String,
    pub output: Option<String>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AppendEvaluationRequest {
    pub project_id: String,
    pub id: String,
    pub relation_id: String,
    pub status: String,
    pub method: Option<String>,
    pub notes: Option<String>,
    pub checked_range: Option<CheckedRange>,
    pub extra: Option<BTreeMap<String, Value>>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AppendCounterexampleRequest {
    pub project_id: String,
    pub id: String,
    pub relation_id: String,
    pub method: Option<String>,
    pub notes: Option<String>,
    pub n: Option<i64>,
    pub first_failure: Option<Value>,
    pub failure_reason: Option<String>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AppendTimeoutRequest {
    pub project_id: String,
    pub id: String,
    pub relation_id: String,
    pub seconds: u64,
    pub method: Option<String>,
    pub notes: Option<String>,
    pub checked_range: Option<CheckedRange>,
}

#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct ProjectListResponse {
    pub projects: Vec<ProjectOverview>,
}

#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct RecordListResponse {
    pub records: Vec<LabRecord>,
}

#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct FamilyListResponse {
    pub families: Vec<PolynomialFamilyInfo>,
}

#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct CheckFamilyRealRootedResponse {
    pub report: CheckFamilyRealRootednessReport,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub written_evidence: Option<WrittenEvaluation>,
}

#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct CheckFamilyInterlacingResponse {
    pub report: CheckFamilyInterlacingReport,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub written_evidence: Option<WrittenEvaluation>,
}

#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct MarkdownResponse {
    pub project_id: String,
    pub markdown: String,
}

#[derive(Debug, Clone, Serialize, JsonSchema)]
pub struct HtmlResponse {
    pub project_id: String,
    pub html: String,
}

#[tool_router(router = tool_router)]
impl PolynomialLabServer {
    pub fn new() -> Self {
        Self::with_root(default_lab_root())
    }

    pub fn with_root(root: impl Into<PathBuf>) -> Self {
        Self {
            root: root.into(),
            tool_router: Self::tool_router(),
        }
    }

    #[tool(description = "Validate the polynomial interlacing lab index.")]
    pub fn validate_lab(
        &self,
        Parameters(input): Parameters<ValidateRequest>,
    ) -> Result<Json<ValidationReport>, McpError> {
        let mode = if input.strict.unwrap_or(false) {
            ValidationMode::Strict
        } else {
            ValidationMode::Tolerant
        };
        Ok(Json(self.load_store()?.validate_with_mode(mode)))
    }

    #[tool(description = "List polynomial interlacing lab projects.")]
    pub fn list_projects(
        &self,
        Parameters(_input): Parameters<EmptyRequest>,
    ) -> Result<Json<ProjectListResponse>, McpError> {
        Ok(Json(ProjectListResponse {
            projects: self.load_store()?.project_overviews(),
        }))
    }

    #[tool(description = "Get one polynomial interlacing lab project with related records.")]
    pub fn get_project(
        &self,
        Parameters(input): Parameters<ProjectRequest>,
    ) -> Result<Json<ProjectReport>, McpError> {
        Ok(Json(self.load_store()?.project_report(&input.project_id)))
    }

    #[tool(description = "List goal records, optionally restricted to one project.")]
    pub fn list_goals(
        &self,
        Parameters(input): Parameters<OptionalProjectRequest>,
    ) -> Result<Json<RecordListResponse>, McpError> {
        Ok(Json(RecordListResponse {
            records: self.load_store()?.goals(input.project_id.as_deref()),
        }))
    }

    #[tool(
        description = "List evidence records with optional project, status, and relation filters."
    )]
    pub fn list_evaluations(
        &self,
        Parameters(input): Parameters<EvaluationSearchRequest>,
    ) -> Result<Json<RecordListResponse>, McpError> {
        Ok(Json(RecordListResponse {
            records: self.load_store()?.evaluations(&EvaluationFilter {
                project_id: input.project_id,
                status: input.status,
                relation_id: input.relation_id,
            }),
        }))
    }

    #[tool(description = "Trace indexed implications and evidence that currently support a goal.")]
    pub fn trace_goal_support(
        &self,
        Parameters(input): Parameters<TraceGoalRequest>,
    ) -> Result<Json<TraceGoalSupport>, McpError> {
        Ok(Json(
            self.load_store()?
                .trace_goal_support(&input.project_id, &input.goal_id),
        ))
    }

    #[tool(description = "Create a project skeleton with standard lab subdirectories.")]
    pub fn create_project(
        &self,
        Parameters(input): Parameters<CreateProjectRequest>,
    ) -> Result<Json<WrittenRecord>, McpError> {
        self.load_store()?
            .create_project(ProjectDraft {
                id: input.project_id,
                label: input.label,
                status: Some(input.status.unwrap_or_else(|| "active".to_string())),
                description: input.description,
                main_objects: input.main_objects.unwrap_or_default(),
                main_goals: input.main_goals.unwrap_or_default(),
                source_notes: input.source_notes.unwrap_or_default(),
            })
            .map(Json)
            .map_err(internal_error)
    }

    #[tool(description = "Append a goal TOML record to a project.")]
    pub fn append_goal(
        &self,
        Parameters(input): Parameters<AppendGoalRequest>,
    ) -> Result<Json<WrittenRecord>, McpError> {
        self.load_store()?
            .append_goal(
                &input.project_id,
                GoalDraft {
                    id: input.id,
                    label: input.label,
                    statement: input.statement,
                    status: input.status.unwrap_or_else(|| "open".to_string()),
                    objects: input.objects.unwrap_or_default(),
                    motivation: input.motivation,
                    current_best_route: input.current_best_route,
                    depends_on: input.depends_on.unwrap_or_default(),
                },
            )
            .map(Json)
            .map_err(internal_error)
    }

    #[tool(description = "Append a conjectural relation TOML record to a project.")]
    pub fn append_conjecture(
        &self,
        Parameters(input): Parameters<AppendConjectureRequest>,
    ) -> Result<Json<WrittenRecord>, McpError> {
        self.load_store()?
            .append_conjecture(
                &input.project_id,
                ConjectureDraft {
                    id: input.id,
                    label: input.label,
                    statement: input.statement,
                    status: input.status.unwrap_or_else(|| "plausible".to_string()),
                    relation: input.relation,
                    left: input.left,
                    right: input.right,
                    index_condition: input.index_condition,
                    depends_on: input.depends_on.unwrap_or_default(),
                },
            )
            .map(Json)
            .map_err(internal_error)
    }

    #[tool(description = "Append an implication/dependency TOML record to a project.")]
    pub fn append_implication(
        &self,
        Parameters(input): Parameters<AppendImplicationRequest>,
    ) -> Result<Json<WrittenRecord>, McpError> {
        self.load_store()?
            .append_implication(
                &input.project_id,
                ImplicationDraft {
                    id: input.id,
                    label: input.label,
                    status: input.status.unwrap_or_else(|| "plausible".to_string()),
                    from: input.from,
                    to: input.to,
                    explanation: input.explanation,
                    proof_tags: input.proof_tags.unwrap_or_default(),
                },
            )
            .map(Json)
            .map_err(internal_error)
    }

    #[tool(description = "List proof-rule records from the interlacing toolbox.")]
    pub fn list_proof_rules(
        &self,
        Parameters(_input): Parameters<EmptyRequest>,
    ) -> Result<Json<RecordListResponse>, McpError> {
        Ok(Json(RecordListResponse {
            records: self.load_store()?.proof_rules(),
        }))
    }

    #[tool(
        description = "List search recipes for discovering recurrences and interlacing relations."
    )]
    pub fn list_search_recipes(
        &self,
        Parameters(_input): Parameters<EmptyRequest>,
    ) -> Result<Json<RecordListResponse>, McpError> {
        Ok(Json(RecordListResponse {
            records: self.load_store()?.search_recipes(),
        }))
    }

    #[tool(description = "List registered computable polynomial families.")]
    pub fn list_families(
        &self,
        Parameters(_input): Parameters<EmptyRequest>,
    ) -> Result<Json<FamilyListResponse>, McpError> {
        Ok(Json(FamilyListResponse {
            families: default_family_registry().list(),
        }))
    }

    #[tool(description = "Compute one registered polynomial family at a given n.")]
    pub fn compute_family(
        &self,
        Parameters(input): Parameters<ComputeFamilyRequest>,
    ) -> Result<Json<ComputedPolynomial>, McpError> {
        default_family_registry()
            .compute(&input.family_id, input.n)
            .map(Json)
            .map_err(internal_error)
    }

    #[tool(
        description = "Check real-rootedness of a registered polynomial family over a range, optionally appending evidence."
    )]
    pub fn check_family_real_rooted(
        &self,
        Parameters(input): Parameters<CheckFamilyRealRootedRequest>,
    ) -> Result<Json<CheckFamilyRealRootedResponse>, McpError> {
        let registry = default_family_registry();
        let report = registry
            .check_real_rooted(&input.family_id, input.n_min, input.n_max)
            .map_err(internal_error)?;
        let written_evidence = if input.append.unwrap_or(false) {
            let project_id = input
                .project_id
                .as_deref()
                .ok_or_else(|| invalid_params("`append` requires `project_id`"))?;
            let relation_id = input
                .relation_id
                .as_deref()
                .ok_or_else(|| invalid_params("`append` requires `relation_id`"))?;
            let id = input.id.unwrap_or_else(|| {
                real_rooted_evidence_id(
                    relation_id,
                    report.first_failure_n,
                    input.n_min,
                    input.n_max,
                )
            });
            let draft = real_rooted_evaluation_draft(
                id,
                relation_id.to_string(),
                &input.family_id,
                &report,
            )
            .map_err(internal_error)?;
            Some(
                self.load_store()?
                    .append_evaluation(project_id, draft)
                    .map_err(internal_error)?,
            )
        } else {
            None
        };
        Ok(Json(CheckFamilyRealRootedResponse {
            report,
            written_evidence,
        }))
    }

    #[tool(
        description = "Check directed strict or weak interlacing between two registered polynomial families over a range, optionally appending evidence."
    )]
    pub fn check_family_interlacing(
        &self,
        Parameters(input): Parameters<CheckFamilyInterlacingRequest>,
    ) -> Result<Json<CheckFamilyInterlacingResponse>, McpError> {
        let mode = input.mode.unwrap_or(InterlacingMode::Weak);
        let registry = default_family_registry();
        let report = registry
            .check_interlacing(
                &input.left_family_id,
                &input.right_family_id,
                input.n_min,
                input.n_max,
                mode,
            )
            .map_err(internal_error)?;
        let written_evidence = if input.append.unwrap_or(false) {
            let project_id = input
                .project_id
                .as_deref()
                .ok_or_else(|| invalid_params("`append` requires `project_id`"))?;
            let relation_id = input
                .relation_id
                .as_deref()
                .ok_or_else(|| invalid_params("`append` requires `relation_id`"))?;
            let id = input.id.unwrap_or_else(|| {
                interlacing_evidence_id(
                    relation_id,
                    mode,
                    report.first_failure_n,
                    input.n_min,
                    input.n_max,
                )
            });
            let draft = interlacing_evaluation_draft(id, relation_id.to_string(), &report)
                .map_err(internal_error)?;
            Some(
                self.load_store()?
                    .append_evaluation(project_id, draft)
                    .map_err(internal_error)?,
            )
        } else {
            None
        };
        Ok(Json(CheckFamilyInterlacingResponse {
            report,
            written_evidence,
        }))
    }

    #[tool(description = "Render a project summary as Markdown.")]
    pub fn render_project_markdown(
        &self,
        Parameters(input): Parameters<ProjectRequest>,
    ) -> Result<Json<MarkdownResponse>, McpError> {
        Ok(Json(MarkdownResponse {
            markdown: self
                .load_store()?
                .render_project_markdown(&input.project_id),
            project_id: input.project_id,
        }))
    }

    #[tool(description = "Render a project summary as static HTML.")]
    pub fn render_project_html(
        &self,
        Parameters(input): Parameters<ProjectRequest>,
    ) -> Result<Json<HtmlResponse>, McpError> {
        Ok(Json(HtmlResponse {
            html: self.load_store()?.render_project_html(&input.project_id),
            project_id: input.project_id,
        }))
    }

    #[tool(description = "Write a generated Markdown project summary file.")]
    pub fn write_project_markdown(
        &self,
        Parameters(input): Parameters<GeneratedProjectRequest>,
    ) -> Result<Json<GeneratedFile>, McpError> {
        let output = input.output.as_deref().map(PathBuf::from);
        self.load_store()?
            .write_project_markdown(&input.project_id, output.as_deref())
            .map(Json)
            .map_err(internal_error)
    }

    #[tool(description = "Write a generated HTML project summary file.")]
    pub fn write_project_html(
        &self,
        Parameters(input): Parameters<GeneratedProjectRequest>,
    ) -> Result<Json<GeneratedFile>, McpError> {
        let output = input.output.as_deref().map(PathBuf::from);
        self.load_store()?
            .write_project_html(&input.project_id, output.as_deref())
            .map(Json)
            .map_err(internal_error)
    }

    #[tool(description = "Append a machine-readable evaluation JSON record.")]
    pub fn append_evaluation(
        &self,
        Parameters(input): Parameters<AppendEvaluationRequest>,
    ) -> Result<Json<WrittenEvaluation>, McpError> {
        self.load_store()?
            .append_evaluation(
                &input.project_id,
                EvaluationDraft {
                    id: input.id,
                    relation_id: input.relation_id,
                    status: input.status,
                    method: input.method,
                    notes: input.notes,
                    checked_range: input.checked_range,
                    first_failure: None,
                    failure_reason: None,
                    timeout_seconds: None,
                    extra: input.extra.unwrap_or_default(),
                },
            )
            .map(Json)
            .map_err(internal_error)
    }

    #[tool(description = "Append a counterexample evaluation JSON record.")]
    pub fn append_counterexample(
        &self,
        Parameters(input): Parameters<AppendCounterexampleRequest>,
    ) -> Result<Json<WrittenEvaluation>, McpError> {
        let first_failure = match (input.n, input.first_failure) {
            (Some(n), None) => serde_json::json!({ "n": n }),
            (None, Some(value)) => value,
            (Some(_), Some(_)) => {
                return Err(invalid_params(
                    "use either `n` or `first_failure`, not both",
                ))
            }
            (None, None) => {
                return Err(invalid_params(
                    "counterexample records require `n` or `first_failure`",
                ))
            }
        };
        self.load_store()?
            .append_evaluation(
                &input.project_id,
                EvaluationDraft {
                    id: input.id,
                    relation_id: input.relation_id,
                    status: "counterexample_found".to_string(),
                    method: input.method,
                    notes: input.notes,
                    checked_range: None,
                    first_failure: Some(first_failure),
                    failure_reason: input.failure_reason,
                    timeout_seconds: None,
                    extra: BTreeMap::new(),
                },
            )
            .map(Json)
            .map_err(internal_error)
    }

    #[tool(description = "Append a timeout evaluation JSON record.")]
    pub fn append_timeout(
        &self,
        Parameters(input): Parameters<AppendTimeoutRequest>,
    ) -> Result<Json<WrittenEvaluation>, McpError> {
        self.load_store()?
            .append_evaluation(
                &input.project_id,
                EvaluationDraft {
                    id: input.id,
                    relation_id: input.relation_id,
                    status: "timeout".to_string(),
                    method: input.method,
                    notes: input.notes,
                    checked_range: input.checked_range,
                    first_failure: None,
                    failure_reason: None,
                    timeout_seconds: Some(input.seconds),
                    extra: BTreeMap::new(),
                },
            )
            .map(Json)
            .map_err(internal_error)
    }

    fn load_store(&self) -> Result<LabStore, McpError> {
        LabStore::load(&self.root).map_err(|error| {
            McpError::internal_error(
                format!(
                    "failed to load polynomial lab at {}: {error}",
                    self.root.display()
                ),
                None,
            )
        })
    }
}

fn invalid_params(message: impl Into<String>) -> McpError {
    McpError::invalid_params(message.into(), None)
}

fn internal_error(error: anyhow::Error) -> McpError {
    McpError::internal_error(error.to_string(), None)
}

#[tool_handler(router = self.tool_router)]
impl ServerHandler for PolynomialLabServer {
    fn get_info(&self) -> ServerInfo {
        ServerInfo {
            protocol_version: Default::default(),
            capabilities: ServerCapabilities::builder().enable_tools().build(),
            server_info: Implementation {
                name: "polynomial-lab".to_string(),
                title: Some("polynomial-lab".to_string()),
                version: env!("CARGO_PKG_VERSION").to_string(),
                description: Some(
                    "Structured project data tools for polynomial interlacing labs.".to_string(),
                ),
                icons: None,
                website_url: None,
            },
            instructions: Some(
                "Use these tools to recover project goals, evidence, proof routes, and search recipes for polynomial real-rootedness and interlacing projects."
                    .to_string(),
            ),
        }
    }
}

impl Default for PolynomialLabServer {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture_root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../tests/fixtures/minimal_lab")
    }

    fn server() -> PolynomialLabServer {
        PolynomialLabServer::with_root(fixture_root())
    }

    #[test]
    fn lists_demo_project() {
        let Json(response) = server()
            .list_projects(Parameters(EmptyRequest {}))
            .expect("list projects should work");
        assert!(response
            .projects
            .iter()
            .any(|project| project.id == "demo_project"));
    }

    #[test]
    fn traces_goal_support() {
        let Json(trace) = server()
            .trace_goal_support(Parameters(TraceGoalRequest {
                project_id: "demo_project".to_string(),
                goal_id: "demo_real_rootedness_goal".to_string(),
            }))
            .expect("trace should work");
        assert_eq!(trace.incoming_implications.len(), 1);
        assert!(trace.incoming_implications[0]
            .prerequisites
            .iter()
            .any(|item| item.id == "demo_interlaces_envelope"));
    }

    #[test]
    fn filters_evidence() {
        let Json(response) = server()
            .list_evaluations(Parameters(EvaluationSearchRequest {
                project_id: Some("demo_project".to_string()),
                status: Some("holds_for_checked_domain".to_string()),
                relation_id: None,
            }))
            .expect("list evaluations should work");
        assert_eq!(response.records.len(), 1);
    }

    #[test]
    fn computes_registered_family() {
        let Json(response) = server()
            .compute_family(Parameters(ComputeFamilyRequest {
                family_id: "derangement_descent_polynomial".to_string(),
                n: 4,
            }))
            .expect("compute family should work");
        assert_eq!(response.coefficients, vec!["0", "4", "4", "1"]);
    }

    #[test]
    fn checks_registered_family_interlacing() {
        let Json(response) = server()
            .check_family_interlacing(Parameters(CheckFamilyInterlacingRequest {
                left_family_id: "normalized_derangement_descent_polynomial".to_string(),
                right_family_id: "reciprocal_eulerian_derivative_polynomial".to_string(),
                n_min: 5,
                n_max: 7,
                mode: Some(InterlacingMode::Weak),
                project_id: None,
                relation_id: None,
                id: None,
                append: None,
            }))
            .expect("check interlacing should work");
        assert!(response.report.all_interlacing);
        assert!(response.written_evidence.is_none());
    }
}
