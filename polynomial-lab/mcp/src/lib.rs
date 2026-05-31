use polynomial_lab::{
    default_lab_root, CheckedRange, EvaluationDraft, EvaluationFilter, GeneratedFile, LabRecord,
    LabStore, ProjectOverview, ProjectReport, TraceGoalSupport, ValidationMode, ValidationReport,
    WrittenEvaluation,
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
}
