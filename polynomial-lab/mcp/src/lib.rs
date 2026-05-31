use polynomial_lab::{
    default_lab_root, EvaluationFilter, LabRecord, LabStore, ProjectOverview, ProjectReport,
    TraceGoalSupport, ValidationReport,
};
use rmcp::{
    handler::server::{router::tool::ToolRouter, wrapper::Parameters},
    model::{Implementation, ServerCapabilities, ServerInfo},
    tool, tool_handler, tool_router, ErrorData as McpError, Json, ServerHandler,
};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
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
        Parameters(_input): Parameters<EmptyRequest>,
    ) -> Result<Json<ValidationReport>, McpError> {
        Ok(Json(self.load_store()?.validate()))
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
    use polynomial_lab::DEFAULT_LAB_ROOT;

    fn server() -> PolynomialLabServer {
        PolynomialLabServer::with_root(DEFAULT_LAB_ROOT)
    }

    #[test]
    fn lists_derangement_project() {
        let Json(response) = server()
            .list_projects(Parameters(EmptyRequest {}))
            .expect("list projects should work");
        assert!(response
            .projects
            .iter()
            .any(|project| project.id == "derangement_descents"));
    }

    #[test]
    fn traces_goal_support() {
        let Json(trace) = server()
            .trace_goal_support(Parameters(TraceGoalRequest {
                project_id: "derangement_descents".to_string(),
                goal_id: "derangement_descent_real_rootedness".to_string(),
            }))
            .expect("trace should work");
        assert_eq!(trace.incoming_implications.len(), 1);
        assert!(trace.incoming_implications[0]
            .prerequisites
            .iter()
            .any(|item| item.id
                == "normalized_derangement_descent_interlaces_reciprocal_eulerian_derivative"));
    }

    #[test]
    fn filters_evidence() {
        let Json(response) = server()
            .list_evaluations(Parameters(EvaluationSearchRequest {
                project_id: Some("derangement_descents".to_string()),
                status: Some("verified_range".to_string()),
                relation_id: None,
            }))
            .expect("list evaluations should work");
        assert_eq!(response.records.len(), 1);
    }
}
