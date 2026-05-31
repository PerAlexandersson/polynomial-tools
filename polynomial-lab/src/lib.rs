use anyhow::{Context, Result};
use chrono::{SecondsFormat, Utc};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};

pub mod families;

pub use families::{
    default_family_registry, interlacing_evaluation_draft, interlacing_evidence_id,
    real_rooted_evaluation_draft, real_rooted_evidence_id, CheckFamilyInterlacingReport,
    CheckFamilyRealRootednessReport, ComputedPolynomial, FamilyCheckItem,
    FamilyInterlacingCheckItem, InterlacingMode, PolynomialFamilyInfo, PolynomialFamilyRegistry,
};

pub const DEFAULT_LAB_ROOT: &str = "/workspace/projects/polynomial-interlacing-lab";

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct LabRecord {
    pub id: String,
    pub kind: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub project_id: Option<String>,
    pub path: String,
    pub data: Value,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct ProjectOverview {
    pub id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub status: Option<String>,
    pub path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct ProjectReport {
    pub project_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub project: Option<LabRecord>,
    pub definitions: Vec<LabRecord>,
    pub goals: Vec<LabRecord>,
    pub conjectures: Vec<LabRecord>,
    pub implications: Vec<LabRecord>,
    pub evaluations: Vec<LabRecord>,
    pub cache_records: Vec<LabRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct ValidationReport {
    pub root: String,
    pub ok: bool,
    pub counts: BTreeMap<String, usize>,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct TraceGoalSupport {
    pub project_id: String,
    pub goal_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub goal: Option<LabRecord>,
    pub incoming_implications: Vec<ImplicationTrace>,
    pub direct_evaluations: Vec<LabRecord>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct ImplicationTrace {
    pub implication: LabRecord,
    pub prerequisites: Vec<TracePrerequisite>,
    pub evaluations: Vec<LabRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct TracePrerequisite {
    pub id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub record: Option<LabRecord>,
    pub evaluations: Vec<LabRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq, Default)]
pub struct EvaluationFilter {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub project_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub status: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub relation_id: Option<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, JsonSchema, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum ValidationMode {
    #[default]
    Tolerant,
    Strict,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct CheckedRange {
    pub n_min: i64,
    pub n_max: i64,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct EvaluationDraft {
    pub id: String,
    pub relation_id: String,
    pub status: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub method: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub notes: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub checked_range: Option<CheckedRange>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub first_failure: Option<Value>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub failure_reason: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub timeout_seconds: Option<u64>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub extra: BTreeMap<String, Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct WrittenEvaluation {
    pub path: String,
    pub record: LabRecord,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq)]
pub struct WrittenRecord {
    pub path: String,
    pub record: LabRecord,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct ProjectDraft {
    pub id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub status: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub main_objects: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub main_goals: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source_notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct GoalDraft {
    pub id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    pub statement: String,
    pub status: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub objects: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub motivation: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub current_best_route: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub depends_on: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct ConjectureDraft {
    pub id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    pub statement: String,
    pub status: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub relation: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub left: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub right: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub index_condition: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub depends_on: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct ImplicationDraft {
    pub id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    pub status: String,
    pub from: Vec<String>,
    pub to: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub explanation: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub proof_tags: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema, PartialEq, Eq)]
pub struct GeneratedFile {
    pub path: String,
}

#[derive(Debug, Clone)]
pub struct LabStore {
    root: PathBuf,
    records: Vec<LabRecord>,
}

impl LabStore {
    pub fn load(root: impl AsRef<Path>) -> Result<Self> {
        let root = root.as_ref().to_path_buf();
        let mut loader = LabLoader::new(root);
        loader.load()?;
        Ok(Self {
            root: loader.root,
            records: loader.records,
        })
    }

    pub fn load_default() -> Result<Self> {
        Self::load(default_lab_root())
    }

    pub fn root(&self) -> &Path {
        &self.root
    }

    pub fn records(&self) -> &[LabRecord] {
        &self.records
    }

    pub fn record(&self, id: &str) -> Option<&LabRecord> {
        self.records.iter().find(|record| record.id == id)
    }

    pub fn records_by_kind(&self, kind: &str) -> Vec<LabRecord> {
        self.records
            .iter()
            .filter(|record| record.kind == kind)
            .cloned()
            .collect()
    }

    pub fn project_overviews(&self) -> Vec<ProjectOverview> {
        self.records
            .iter()
            .filter(|record| record.kind == "project")
            .map(|record| ProjectOverview {
                id: record.id.clone(),
                label: record.label.clone(),
                status: value_string(&record.data, "status"),
                path: record.path.clone(),
            })
            .collect()
    }

    pub fn project_report(&self, project_id: &str) -> ProjectReport {
        ProjectReport {
            project_id: project_id.to_string(),
            project: self
                .records
                .iter()
                .find(|record| record.kind == "project" && record.id == project_id)
                .cloned(),
            definitions: self.records_for_project(project_id, |kind| {
                matches!(
                    kind,
                    "family" | "operator" | "relation" | "refinement" | "object"
                )
            }),
            goals: self.records_for_project(project_id, |kind| kind == "goal"),
            conjectures: self.records_for_project(project_id, |kind| kind == "conjecture"),
            implications: self.records_for_project(project_id, |kind| kind == "implication"),
            evaluations: self.records_for_project(project_id, |kind| kind == "evaluation"),
            cache_records: self.records_for_project(project_id, |kind| kind == "cache_record"),
        }
    }

    pub fn goals(&self, project_id: Option<&str>) -> Vec<LabRecord> {
        self.records
            .iter()
            .filter(|record| record.kind == "goal")
            .filter(|record| project_matches(record, project_id))
            .cloned()
            .collect()
    }

    pub fn proof_rules(&self) -> Vec<LabRecord> {
        self.records_by_kind("proof_rule")
    }

    pub fn search_recipes(&self) -> Vec<LabRecord> {
        self.records_by_kind("search_recipe")
    }

    pub fn evaluations(&self, filter: &EvaluationFilter) -> Vec<LabRecord> {
        self.records
            .iter()
            .filter(|record| record.kind == "evaluation")
            .filter(|record| project_matches(record, filter.project_id.as_deref()))
            .filter(|record| {
                filter.status.as_deref().is_none_or(|status| {
                    value_string(&record.data, "status").as_deref() == Some(status)
                })
            })
            .filter(|record| {
                filter.relation_id.as_deref().is_none_or(|relation_id| {
                    value_string(&record.data, "relation_id").as_deref() == Some(relation_id)
                })
            })
            .cloned()
            .collect()
    }

    pub fn trace_goal_support(&self, project_id: &str, goal_id: &str) -> TraceGoalSupport {
        let mut warnings = Vec::new();
        let goal = self.record(goal_id).cloned();
        if goal.is_none() {
            warnings.push(format!("goal id '{goal_id}' is not indexed"));
        }

        let incoming_implications = self
            .records
            .iter()
            .filter(|record| record.kind == "implication")
            .filter(|record| project_matches(record, Some(project_id)))
            .filter(|record| value_string(&record.data, "to").as_deref() == Some(goal_id))
            .map(|record| {
                let prerequisites = value_string_array(&record.data, "from")
                    .into_iter()
                    .map(|id| TracePrerequisite {
                        evaluations: self.evaluations(&EvaluationFilter {
                            project_id: Some(project_id.to_string()),
                            relation_id: Some(id.clone()),
                            ..Default::default()
                        }),
                        record: self.record(&id).cloned(),
                        id,
                    })
                    .collect();
                ImplicationTrace {
                    implication: record.clone(),
                    prerequisites,
                    evaluations: self.evaluations(&EvaluationFilter {
                        project_id: Some(project_id.to_string()),
                        relation_id: Some(record.id.clone()),
                        ..Default::default()
                    }),
                }
            })
            .collect();

        TraceGoalSupport {
            project_id: project_id.to_string(),
            goal_id: goal_id.to_string(),
            goal,
            incoming_implications,
            direct_evaluations: self.evaluations(&EvaluationFilter {
                project_id: Some(project_id.to_string()),
                relation_id: Some(goal_id.to_string()),
                ..Default::default()
            }),
            warnings,
        }
    }

    pub fn append_evaluation(
        &self,
        project_id: &str,
        draft: EvaluationDraft,
    ) -> Result<WrittenEvaluation> {
        self.require_project(project_id)?;
        validate_record_id(&draft.id)?;
        validate_record_id(&draft.relation_id)?;
        validate_record_id(project_id)?;

        let evidence_dir = self.root.join("projects").join(project_id).join("evidence");
        fs::create_dir_all(&evidence_dir)
            .with_context(|| format!("failed to create {}", evidence_dir.display()))?;
        let path = evidence_dir.join(format!("{}.json", draft.id));
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .with_context(|| format!("failed to create evidence file {}", path.display()))?;

        let value = draft_to_evaluation_json(project_id, &draft)?;
        let text = serde_json::to_string_pretty(&value)?;
        writeln!(file, "{text}")
            .with_context(|| format!("failed to write evidence file {}", path.display()))?;

        let record = LabRecord {
            id: draft.id,
            kind: "evaluation".to_string(),
            label: None,
            project_id: Some(project_id.to_string()),
            path: relative_path(&self.root, &path),
            data: value,
        };
        Ok(WrittenEvaluation {
            path: record.path.clone(),
            record,
        })
    }

    pub fn create_project(&self, draft: ProjectDraft) -> Result<WrittenRecord> {
        validate_record_id(&draft.id)?;
        self.ensure_new_record_id(&draft.id)?;
        for id in draft.main_objects.iter().chain(draft.main_goals.iter()) {
            validate_record_id(id)?;
        }

        let project_dir = self.root.join("projects").join(&draft.id);
        let path = project_dir.join("project.toml");
        let value = serde_json::to_value(&draft)?;
        let written = self.write_toml_record(
            &path,
            draft.id.clone(),
            "project",
            draft.label.clone(),
            Some(draft.id.clone()),
            value,
        )?;
        for directory in [
            "definitions",
            "goals",
            "conjectures",
            "implications",
            "evidence",
            "cache",
            "generated",
        ] {
            fs::create_dir_all(project_dir.join(directory)).with_context(|| {
                format!(
                    "failed to create project subdirectory {}",
                    project_dir.join(directory).display()
                )
            })?;
        }
        Ok(written)
    }

    pub fn append_goal(&self, project_id: &str, draft: GoalDraft) -> Result<WrittenRecord> {
        self.require_project(project_id)?;
        validate_record_id(&draft.id)?;
        self.ensure_new_record_id(&draft.id)?;
        validate_record_id(project_id)?;
        for id in draft.objects.iter().chain(draft.depends_on.iter()) {
            validate_record_id(id)?;
        }
        validate_known_status(&draft.status, GOAL_STATUSES, "goal")?;

        let mut value = typed_record_value("goal", &draft)?;
        ensure_project_id(&mut value, project_id);
        let path = self
            .root
            .join("projects")
            .join(project_id)
            .join("goals")
            .join(format!("{}.toml", draft.id));
        self.write_toml_record(
            &path,
            draft.id,
            "goal",
            draft.label,
            Some(project_id.to_string()),
            value,
        )
    }

    pub fn append_conjecture(
        &self,
        project_id: &str,
        draft: ConjectureDraft,
    ) -> Result<WrittenRecord> {
        self.require_project(project_id)?;
        validate_record_id(&draft.id)?;
        self.ensure_new_record_id(&draft.id)?;
        validate_record_id(project_id)?;
        for id in draft
            .left
            .iter()
            .chain(draft.right.iter())
            .chain(draft.depends_on.iter())
        {
            validate_record_id(id)?;
        }
        validate_known_status(&draft.status, CONJECTURE_STATUSES, "conjecture")?;

        let mut value = typed_record_value("conjecture", &draft)?;
        ensure_project_id(&mut value, project_id);
        let path = self
            .root
            .join("projects")
            .join(project_id)
            .join("conjectures")
            .join(format!("{}.toml", draft.id));
        self.write_toml_record(
            &path,
            draft.id,
            "conjecture",
            draft.label,
            Some(project_id.to_string()),
            value,
        )
    }

    pub fn append_implication(
        &self,
        project_id: &str,
        draft: ImplicationDraft,
    ) -> Result<WrittenRecord> {
        self.require_project(project_id)?;
        validate_record_id(&draft.id)?;
        self.ensure_new_record_id(&draft.id)?;
        validate_record_id(project_id)?;
        if draft.from.is_empty() {
            anyhow::bail!("implication '{}' needs at least one prerequisite", draft.id);
        }
        for id in draft.from.iter().chain(std::iter::once(&draft.to)) {
            validate_record_id(id)?;
        }
        validate_known_status(&draft.status, IMPLICATION_STATUSES, "implication")?;

        let mut value = typed_record_value("implication", &draft)?;
        ensure_project_id(&mut value, project_id);
        let path = self
            .root
            .join("projects")
            .join(project_id)
            .join("implications")
            .join(format!("{}.toml", draft.id));
        self.write_toml_record(
            &path,
            draft.id,
            "implication",
            draft.label,
            Some(project_id.to_string()),
            value,
        )
    }

    pub fn render_project_html(&self, project_id: &str) -> String {
        let report = self.project_report(project_id);
        let title = self
            .project_report(project_id)
            .project
            .and_then(|record| record.label)
            .unwrap_or_else(|| project_id.to_string());
        let mut html = String::new();
        html.push_str("<!doctype html>\n<html lang=\"en\">\n<head>\n");
        html.push_str("  <meta charset=\"utf-8\">\n");
        html.push_str(&format!("  <title>{}</title>\n", escape_html(&title)));
        html.push_str(&format!("  <style>{}</style>\n", HTML_STYLE));
        html.push_str("</head>\n<body>\n<main>\n");
        html.push_str(&format!("<h1>{}</h1>\n", escape_html(&title)));

        if let Some(project) = &report.project {
            html.push_str("<dl class=\"metadata\">\n");
            html.push_str(&format!(
                "  <dt>id</dt><dd><code>{}</code></dd>\n",
                escape_html(&project.id)
            ));
            if let Some(status) = value_string(&project.data, "status") {
                html.push_str(&format!(
                    "  <dt>status</dt><dd><code>{}</code></dd>\n",
                    escape_html(&status)
                ));
            }
            html.push_str(&format!(
                "  <dt>source</dt><dd><code>{}</code></dd>\n",
                escape_html(&project.path)
            ));
            html.push_str("</dl>\n");
            if let Some(description) = value_string(&project.data, "description") {
                html.push_str(&format!("<p>{}</p>\n", escape_html(description.trim())));
            }
        } else {
            html.push_str("<p>Project record not found.</p>\n");
        }

        render_html_records_section(&mut html, "Goals", &report.goals);
        render_html_records_section(&mut html, "Definitions", &report.definitions);
        render_html_records_section(&mut html, "Conjectures", &report.conjectures);
        render_html_implications_section(&mut html, &report.implications);
        render_html_evaluations_section(&mut html, &report.evaluations);
        html.push_str("</main>\n</body>\n</html>\n");
        html
    }

    pub fn write_project_markdown(
        &self,
        project_id: &str,
        output: Option<&Path>,
    ) -> Result<GeneratedFile> {
        let path = self.generated_output_path(project_id, output, "project-summary.md")?;
        write_generated_file(&path, &self.render_project_markdown(project_id))?;
        Ok(GeneratedFile {
            path: relative_path(&self.root, &path),
        })
    }

    pub fn write_project_html(
        &self,
        project_id: &str,
        output: Option<&Path>,
    ) -> Result<GeneratedFile> {
        let path = self.generated_output_path(project_id, output, "project-summary.html")?;
        write_generated_file(&path, &self.render_project_html(project_id))?;
        Ok(GeneratedFile {
            path: relative_path(&self.root, &path),
        })
    }

    fn require_project(&self, project_id: &str) -> Result<()> {
        if self
            .records
            .iter()
            .any(|record| record.kind == "project" && record.id == project_id)
        {
            Ok(())
        } else {
            anyhow::bail!("unknown project id '{project_id}'");
        }
    }

    fn ensure_new_record_id(&self, id: &str) -> Result<()> {
        if self.record(id).is_some() {
            anyhow::bail!("record id '{id}' already exists");
        }
        Ok(())
    }

    fn write_toml_record(
        &self,
        path: &Path,
        id: String,
        kind: &str,
        label: Option<String>,
        project_id: Option<String>,
        data: Value,
    ) -> Result<WrittenRecord> {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)
                .with_context(|| format!("failed to create {}", parent.display()))?;
        }
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)
            .with_context(|| format!("failed to create TOML record {}", path.display()))?;
        let text = toml::to_string_pretty(&data)
            .with_context(|| format!("failed to serialize TOML record {}", path.display()))?;
        writeln!(file, "{text}")
            .with_context(|| format!("failed to write TOML record {}", path.display()))?;
        let record = LabRecord {
            id,
            kind: kind.to_string(),
            label,
            project_id,
            path: relative_path(&self.root, path),
            data,
        };
        Ok(WrittenRecord {
            path: record.path.clone(),
            record,
        })
    }

    fn generated_output_path(
        &self,
        project_id: &str,
        output: Option<&Path>,
        default_filename: &str,
    ) -> Result<PathBuf> {
        self.require_project(project_id)?;
        let path = output.map_or_else(
            || {
                self.root
                    .join("projects")
                    .join(project_id)
                    .join("generated")
                    .join(default_filename)
            },
            |path| {
                if path.is_absolute() {
                    path.to_path_buf()
                } else {
                    self.root.join(path)
                }
            },
        );
        Ok(path)
    }

    pub fn validate(&self) -> ValidationReport {
        self.validate_with_mode(ValidationMode::Tolerant)
    }

    pub fn validate_with_mode(&self, mode: ValidationMode) -> ValidationReport {
        let mut warnings = Vec::new();
        let mut errors = Vec::new();
        let mut counts = BTreeMap::new();
        for record in &self.records {
            *counts.entry(record.kind.clone()).or_insert(0) += 1;
        }

        let mut by_id: BTreeMap<&str, Vec<&LabRecord>> = BTreeMap::new();
        for record in &self.records {
            by_id.entry(&record.id).or_default().push(record);
        }
        for (id, records) in &by_id {
            if records.len() > 1 {
                let paths = records
                    .iter()
                    .map(|record| record.path.as_str())
                    .collect::<Vec<_>>()
                    .join(", ");
                add_validation_issue(
                    mode,
                    &mut warnings,
                    &mut errors,
                    format!("duplicate id '{id}' appears in {paths}"),
                );
            }
        }

        for project in self
            .records
            .iter()
            .filter(|record| record.kind == "project")
        {
            for goal_id in value_string_array(&project.data, "main_goals") {
                if self.record(&goal_id).is_none() {
                    add_validation_issue(
                        mode,
                        &mut warnings,
                        &mut errors,
                        format!(
                            "project '{}' references missing main goal '{}'",
                            project.id, goal_id
                        ),
                    );
                }
            }
            for object_id in value_string_array(&project.data, "main_objects") {
                if self.record(&object_id).is_none() {
                    add_validation_issue(
                        mode,
                        &mut warnings,
                        &mut errors,
                        format!(
                            "project '{}' references missing main object '{}'",
                            project.id, object_id
                        ),
                    );
                }
            }
            for source_note in value_string_array(&project.data, "source_notes") {
                if !Path::new(&source_note).exists() {
                    add_validation_issue(
                        mode,
                        &mut warnings,
                        &mut errors,
                        format!(
                            "project '{}' source note does not exist: {}",
                            project.id, source_note
                        ),
                    );
                }
            }
        }

        for record in &self.records {
            validate_record_status(mode, record, &mut warnings, &mut errors);

            for dependency_id in value_string_array(&record.data, "depends_on") {
                if self.record(&dependency_id).is_none() {
                    add_validation_issue(
                        mode,
                        &mut warnings,
                        &mut errors,
                        format!(
                            "record '{}' depends on missing id '{}'",
                            record.id, dependency_id
                        ),
                    );
                }
            }

            if record.kind == "goal" {
                for object_id in value_string_array(&record.data, "objects") {
                    if self.record(&object_id).is_none() {
                        add_validation_issue(
                            mode,
                            &mut warnings,
                            &mut errors,
                            format!(
                                "goal '{}' references missing object '{}'",
                                record.id, object_id
                            ),
                        );
                    }
                }
            }

            if record.kind == "conjecture" {
                for key in ["left", "right"] {
                    if let Some(id) = value_string(&record.data, key) {
                        if self.record(&id).is_none() {
                            add_validation_issue(
                                mode,
                                &mut warnings,
                                &mut errors,
                                format!(
                                    "conjecture '{}' references missing {} '{}'",
                                    record.id, key, id
                                ),
                            );
                        }
                    }
                }
            }

            if record.kind == "implication" {
                for source_id in value_string_array(&record.data, "from") {
                    if self.record(&source_id).is_none() {
                        add_validation_issue(
                            mode,
                            &mut warnings,
                            &mut errors,
                            format!(
                                "implication '{}' has missing prerequisite '{}'",
                                record.id, source_id
                            ),
                        );
                    }
                }
                if let Some(target_id) = value_string(&record.data, "to") {
                    if self.record(&target_id).is_none() {
                        add_validation_issue(
                            mode,
                            &mut warnings,
                            &mut errors,
                            format!(
                                "implication '{}' has missing target '{}'",
                                record.id, target_id
                            ),
                        );
                    }
                } else {
                    errors.push(format!("implication '{}' is missing field 'to'", record.id));
                }
            }

            if record.kind == "evaluation" {
                if let Some(relation_id) = value_string(&record.data, "relation_id") {
                    if self.record(&relation_id).is_none() {
                        add_validation_issue(
                            mode,
                            &mut warnings,
                            &mut errors,
                            format!(
                                "evaluation '{}' references missing relation '{}'",
                                record.id, relation_id
                            ),
                        );
                    }
                } else {
                    add_validation_issue(
                        mode,
                        &mut warnings,
                        &mut errors,
                        format!("evaluation '{}' has no relation_id", record.id),
                    );
                }
            }
        }

        ValidationReport {
            root: self.root.display().to_string(),
            ok: errors.is_empty(),
            counts,
            warnings,
            errors,
        }
    }

    pub fn render_project_markdown(&self, project_id: &str) -> String {
        let report = self.project_report(project_id);
        let mut out = String::new();
        let title = report
            .project
            .as_ref()
            .and_then(|record| record.label.as_deref())
            .unwrap_or(project_id);
        push_line(&mut out, &format!("# {title}"));
        push_line(&mut out, "");

        if let Some(project) = &report.project {
            push_line(&mut out, &format!("- id: `{}`", project.id));
            if let Some(status) = value_string(&project.data, "status") {
                push_line(&mut out, &format!("- status: `{status}`"));
            }
            push_line(&mut out, &format!("- source: `{}`", project.path));
            if let Some(description) = value_string(&project.data, "description") {
                push_line(&mut out, "");
                push_line(&mut out, description.trim());
            }
        } else {
            push_line(&mut out, "Project record not found.");
        }

        render_records_section(&mut out, "Goals", &report.goals);
        render_records_section(&mut out, "Definitions", &report.definitions);
        render_records_section(&mut out, "Conjectures", &report.conjectures);
        render_implications_section(&mut out, &report.implications);
        render_evaluations_section(&mut out, &report.evaluations);

        out
    }

    fn records_for_project<F>(&self, project_id: &str, kind_predicate: F) -> Vec<LabRecord>
    where
        F: Fn(&str) -> bool,
    {
        self.records
            .iter()
            .filter(|record| project_matches(record, Some(project_id)))
            .filter(|record| kind_predicate(&record.kind))
            .cloned()
            .collect()
    }
}

pub fn default_lab_root() -> PathBuf {
    env::var_os("POLY_LAB_ROOT")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(DEFAULT_LAB_ROOT))
}

const EVALUATION_STATUSES: &[&str] = &[
    "holds_for_checked_domain",
    "does_not_hold",
    "counterexample_found",
    "incomplete_scan",
    "timeout",
    "method_not_applicable",
    "proof_strategy_failed",
    "obsolete_or_superseded",
    "proved",
    "verified_range",
];

const GOAL_STATUSES: &[&str] = &[
    "open",
    "partially_proved",
    "proved",
    "false",
    "retired",
    "superseded",
];

const IMPLICATION_STATUSES: &[&str] = &[
    "proved",
    "standard",
    "plausible",
    "open",
    "blocked",
    "false",
    "superseded",
];

const CONJECTURE_STATUSES: &[&str] = &[
    "open",
    "plausible",
    "proved",
    "false",
    "superseded",
    "retired",
];

const HTML_STYLE: &str = "\
body { margin: 0; font-family: system-ui, sans-serif; line-height: 1.5; color: #1f2933; background: #f7f8fa; }
main { max-width: 1100px; margin: 0 auto; padding: 32px 24px; }
h1 { margin: 0 0 16px; font-size: 2rem; }
h2 { margin: 32px 0 12px; font-size: 1.25rem; }
.metadata { display: grid; grid-template-columns: max-content 1fr; gap: 4px 12px; margin: 16px 0; }
.metadata dt { font-weight: 700; }
section { background: #fff; border: 1px solid #d9dee7; border-radius: 6px; padding: 18px 20px; margin: 18px 0; }
ul { padding-left: 22px; }
li { margin: 10px 0; }
code { background: #eef2f7; padding: 1px 4px; border-radius: 4px; }
.detail { margin: 2px 0 0; color: #485465; }
";

fn add_validation_issue(
    mode: ValidationMode,
    warnings: &mut Vec<String>,
    errors: &mut Vec<String>,
    message: String,
) {
    match mode {
        ValidationMode::Tolerant => warnings.push(message),
        ValidationMode::Strict => errors.push(message),
    }
}

fn validate_record_status(
    mode: ValidationMode,
    record: &LabRecord,
    warnings: &mut Vec<String>,
    errors: &mut Vec<String>,
) {
    let Some(status) = value_string(&record.data, "status") else {
        return;
    };
    let allowed = match record.kind.as_str() {
        "evaluation" => Some(EVALUATION_STATUSES),
        "goal" => Some(GOAL_STATUSES),
        "implication" => Some(IMPLICATION_STATUSES),
        "conjecture" => Some(CONJECTURE_STATUSES),
        _ => None,
    };
    if let Some(allowed) = allowed {
        if !allowed.contains(&status.as_str()) {
            add_validation_issue(
                mode,
                warnings,
                errors,
                format!(
                    "record '{}' has unknown {} status '{}'",
                    record.id, record.kind, status
                ),
            );
        }
    }
}

fn validate_record_id(id: &str) -> Result<()> {
    if id.is_empty() {
        anyhow::bail!("record id must not be empty");
    }
    if !id
        .bytes()
        .all(|byte| byte.is_ascii_alphanumeric() || byte == b'_' || byte == b'-')
    {
        anyhow::bail!(
            "record id '{id}' contains unsupported characters; use ASCII letters, digits, '_' or '-'"
        );
    }
    Ok(())
}

fn validate_known_status(status: &str, allowed: &[&str], record_kind: &str) -> Result<()> {
    if allowed.contains(&status) {
        Ok(())
    } else {
        anyhow::bail!("unknown {record_kind} status '{status}'");
    }
}

fn typed_record_value<T: Serialize>(record_type: &str, draft: &T) -> Result<Value> {
    let mut value = serde_json::to_value(draft)?;
    let object = value
        .as_object_mut()
        .with_context(|| "record draft did not serialize to a TOML table")?;
    object.insert("type".to_string(), Value::String(record_type.to_string()));
    Ok(value)
}

fn draft_to_evaluation_json(project_id: &str, draft: &EvaluationDraft) -> Result<Value> {
    if !EVALUATION_STATUSES.contains(&draft.status.as_str()) {
        anyhow::bail!("unknown evaluation status '{}'", draft.status);
    }
    if let Some(range) = &draft.checked_range {
        if range.n_min > range.n_max {
            anyhow::bail!("checked_range must satisfy n_min <= n_max");
        }
    }

    let mut object = serde_json::Map::new();
    object.insert("id".to_string(), Value::String(draft.id.clone()));
    object.insert(
        "project_id".to_string(),
        Value::String(project_id.to_string()),
    );
    object.insert(
        "relation_id".to_string(),
        Value::String(draft.relation_id.clone()),
    );
    object.insert("status".to_string(), Value::String(draft.status.clone()));
    object.insert(
        "created_utc".to_string(),
        Value::String(Utc::now().to_rfc3339_opts(SecondsFormat::Secs, true)),
    );
    object.insert(
        "software".to_string(),
        json!({
            "language": "rust",
            "crate": "polynomial-lab"
        }),
    );
    if let Some(method) = &draft.method {
        object.insert("method".to_string(), Value::String(method.clone()));
    }
    if let Some(notes) = &draft.notes {
        object.insert("notes".to_string(), Value::String(notes.clone()));
    }
    if let Some(range) = &draft.checked_range {
        object.insert(
            "checked_range".to_string(),
            json!({
                "n_min": range.n_min,
                "n_max": range.n_max
            }),
        );
    }
    if let Some(first_failure) = &draft.first_failure {
        object.insert("first_failure".to_string(), first_failure.clone());
    }
    if let Some(failure_reason) = &draft.failure_reason {
        object.insert(
            "failure_reason".to_string(),
            Value::String(failure_reason.clone()),
        );
    }
    if let Some(timeout_seconds) = draft.timeout_seconds {
        object.insert("timeout_seconds".to_string(), json!(timeout_seconds));
    }
    for (key, extra_value) in &draft.extra {
        object.insert(key.clone(), extra_value.clone());
    }
    Ok(Value::Object(object))
}

fn write_generated_file(path: &Path, content: &str) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    fs::write(path, content).with_context(|| format!("failed to write {}", path.display()))
}

fn escape_html(text: &str) -> String {
    let mut escaped = String::with_capacity(text.len());
    for ch in text.chars() {
        match ch {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            '"' => escaped.push_str("&quot;"),
            '\'' => escaped.push_str("&#39;"),
            _ => escaped.push(ch),
        }
    }
    escaped
}

fn project_matches(record: &LabRecord, project_id: Option<&str>) -> bool {
    project_id.is_none_or(|project_id| {
        record.project_id.as_deref() == Some(project_id)
            || (record.kind == "project" && record.id == project_id)
    })
}

fn value_string(value: &Value, key: &str) -> Option<String> {
    value.get(key)?.as_str().map(ToString::to_string)
}

fn value_string_array(value: &Value, key: &str) -> Vec<String> {
    value
        .get(key)
        .and_then(Value::as_array)
        .map(|items| {
            items
                .iter()
                .filter_map(Value::as_str)
                .map(ToString::to_string)
                .collect()
        })
        .unwrap_or_default()
}

struct LabLoader {
    root: PathBuf,
    records: Vec<LabRecord>,
}

impl LabLoader {
    fn new(root: PathBuf) -> Self {
        Self {
            root,
            records: Vec::new(),
        }
    }

    fn load(&mut self) -> Result<()> {
        if !self.root.exists() {
            anyhow::bail!("lab root does not exist: {}", self.root.display());
        }
        self.load_projects()?;
        self.load_top_level_array_records("proof_rules", "rules", "proof_rule")?;
        self.load_top_level_array_records("search_recipes", "recipes", "search_recipe")?;
        self.load_top_level_json_records("cache", "cache_record")?;
        self.records.sort_by(|a, b| {
            a.kind
                .cmp(&b.kind)
                .then_with(|| a.project_id.cmp(&b.project_id))
                .then_with(|| a.id.cmp(&b.id))
                .then_with(|| a.path.cmp(&b.path))
        });
        Ok(())
    }

    fn load_projects(&mut self) -> Result<()> {
        let projects_dir = self.root.join("projects");
        if !projects_dir.exists() {
            return Ok(());
        }

        let mut project_dirs = child_dirs(&projects_dir)?;
        project_dirs.sort();
        for project_dir in project_dirs {
            let project_file = project_dir.join("project.toml");
            if !project_file.exists() {
                continue;
            }
            let project_value = load_toml_json(&project_file)?;
            let project_id = required_id(&project_value, &project_file)?;
            self.push_record(
                project_id.clone(),
                "project".to_string(),
                value_string(&project_value, "label"),
                Some(project_id.clone()),
                &project_file,
                project_value,
            );

            self.load_project_definition_records(&project_dir, &project_id)?;
            self.load_project_toml_directory(&project_dir, &project_id, "goals", "goal")?;
            self.load_project_toml_directory(
                &project_dir,
                &project_id,
                "conjectures",
                "conjecture",
            )?;
            self.load_project_toml_directory(
                &project_dir,
                &project_id,
                "implications",
                "implication",
            )?;
            self.load_project_json_directory(&project_dir, &project_id, "evidence", "evaluation")?;
            self.load_project_json_directory(&project_dir, &project_id, "cache", "cache_record")?;
        }
        Ok(())
    }

    fn load_project_definition_records(
        &mut self,
        project_dir: &Path,
        project_id: &str,
    ) -> Result<()> {
        let definitions_dir = project_dir.join("definitions");
        for path in files_with_extension(&definitions_dir, "toml")? {
            let value = load_toml_json(&path)?;
            let local_project_id =
                value_string(&value, "project_id").unwrap_or_else(|| project_id.to_string());
            for (array_key, item) in table_array_items(&value) {
                let Some(id) = value_string(item, "id") else {
                    continue;
                };
                let mut data = item.clone();
                ensure_project_id(&mut data, &local_project_id);
                let kind = value_string(&data, "type")
                    .unwrap_or_else(|| singular_record_kind(&array_key).to_string());
                self.push_record(
                    id,
                    kind,
                    value_string(&data, "label"),
                    Some(local_project_id.clone()),
                    &path,
                    data,
                );
            }
        }
        Ok(())
    }

    fn load_project_toml_directory(
        &mut self,
        project_dir: &Path,
        project_id: &str,
        directory_name: &str,
        default_kind: &str,
    ) -> Result<()> {
        let directory = project_dir.join(directory_name);
        for path in files_with_extension(&directory, "toml")? {
            let mut value = load_toml_json(&path)?;
            let id = required_id(&value, &path)?;
            ensure_project_id(&mut value, project_id);
            let kind = value_string(&value, "type").unwrap_or_else(|| default_kind.to_string());
            self.push_record(
                id,
                kind,
                value_string(&value, "label"),
                Some(project_id.to_string()),
                &path,
                value,
            );
        }
        Ok(())
    }

    fn load_project_json_directory(
        &mut self,
        project_dir: &Path,
        project_id: &str,
        directory_name: &str,
        default_kind: &str,
    ) -> Result<()> {
        let directory = project_dir.join(directory_name);
        for path in files_with_extension(&directory, "json")? {
            let mut value = load_json(&path)?;
            let id = required_id(&value, &path)?;
            ensure_project_id(&mut value, project_id);
            let kind = value_string(&value, "type").unwrap_or_else(|| default_kind.to_string());
            self.push_record(
                id,
                kind,
                value_string(&value, "label"),
                Some(project_id.to_string()),
                &path,
                value,
            );
        }
        Ok(())
    }

    fn load_top_level_array_records(
        &mut self,
        directory_name: &str,
        array_key: &str,
        kind: &str,
    ) -> Result<()> {
        let directory = self.root.join(directory_name);
        for path in files_with_extension(&directory, "toml")? {
            let value = load_toml_json(&path)?;
            if let Some(items) = value.get(array_key).and_then(Value::as_array) {
                for item in items {
                    let Some(id) = value_string(item, "id") else {
                        continue;
                    };
                    self.push_record(
                        id,
                        kind.to_string(),
                        value_string(item, "label"),
                        None,
                        &path,
                        item.clone(),
                    );
                }
            }
        }
        Ok(())
    }

    fn load_top_level_json_records(
        &mut self,
        directory_name: &str,
        default_kind: &str,
    ) -> Result<()> {
        let directory = self.root.join(directory_name);
        for path in files_with_extension(&directory, "json")? {
            let value = load_json(&path)?;
            let id = required_id(&value, &path)?;
            let kind = value_string(&value, "type").unwrap_or_else(|| default_kind.to_string());
            self.push_record(
                id,
                kind,
                value_string(&value, "label"),
                value_string(&value, "project_id"),
                &path,
                value,
            );
        }
        Ok(())
    }

    fn push_record(
        &mut self,
        id: String,
        kind: String,
        label: Option<String>,
        project_id: Option<String>,
        path: &Path,
        data: Value,
    ) {
        self.records.push(LabRecord {
            id,
            kind: normalize_kind(&kind).to_string(),
            label,
            project_id,
            path: relative_path(&self.root, path),
            data,
        });
    }
}

fn load_toml_json(path: &Path) -> Result<Value> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("failed to read TOML file {}", path.display()))?;
    let value: toml::Value = toml::from_str(&text)
        .with_context(|| format!("failed to parse TOML file {}", path.display()))?;
    serde_json::to_value(value)
        .with_context(|| format!("failed to convert TOML file {} to JSON", path.display()))
}

fn load_json(path: &Path) -> Result<Value> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("failed to read JSON file {}", path.display()))?;
    serde_json::from_str(&text)
        .with_context(|| format!("failed to parse JSON file {}", path.display()))
}

fn required_id(value: &Value, path: &Path) -> Result<String> {
    value_string(value, "id")
        .with_context(|| format!("record in {} is missing string id", path.display()))
}

fn ensure_project_id(value: &mut Value, project_id: &str) {
    if let Some(object) = value.as_object_mut() {
        object
            .entry("project_id".to_string())
            .or_insert_with(|| Value::String(project_id.to_string()));
    }
}

fn table_array_items(value: &Value) -> Vec<(String, &Value)> {
    let Some(object) = value.as_object() else {
        return Vec::new();
    };
    let mut items = Vec::new();
    for (key, value) in object {
        if let Some(array) = value.as_array() {
            for item in array {
                if item.is_object() {
                    items.push((key.clone(), item));
                }
            }
        }
    }
    items
}

fn singular_record_kind(array_key: &str) -> &str {
    match array_key {
        "families" => "family",
        "operators" => "operator",
        "relations" => "relation",
        "refinements" => "refinement",
        "objects" => "object",
        "recipes" => "search_recipe",
        "rules" => "proof_rule",
        _ => array_key.strip_suffix('s').unwrap_or(array_key),
    }
}

fn normalize_kind(kind: &str) -> &str {
    match kind {
        "polynomial_family" => "family",
        "linear_operator" => "operator",
        "proof-rule" => "proof_rule",
        "search-recipe" => "search_recipe",
        other => other,
    }
}

fn child_dirs(root: &Path) -> Result<Vec<PathBuf>> {
    if !root.exists() {
        return Ok(Vec::new());
    }
    let mut dirs = Vec::new();
    for entry in fs::read_dir(root).with_context(|| format!("failed to read {}", root.display()))? {
        let entry = entry.with_context(|| format!("failed to read entry in {}", root.display()))?;
        let path = entry.path();
        if path.is_dir() {
            dirs.push(path);
        }
    }
    Ok(dirs)
}

fn files_with_extension(root: &Path, extension: &str) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    collect_files_with_extension(root, extension, &mut files)?;
    files.sort();
    Ok(files)
}

fn collect_files_with_extension(
    root: &Path,
    extension: &str,
    files: &mut Vec<PathBuf>,
) -> Result<()> {
    if !root.exists() {
        return Ok(());
    }
    for entry in fs::read_dir(root).with_context(|| format!("failed to read {}", root.display()))? {
        let entry = entry.with_context(|| format!("failed to read entry in {}", root.display()))?;
        let path = entry.path();
        if path.is_dir() {
            collect_files_with_extension(&path, extension, files)?;
        } else if path.extension().and_then(|ext| ext.to_str()) == Some(extension) {
            files.push(path);
        }
    }
    Ok(())
}

fn relative_path(root: &Path, path: &Path) -> String {
    path.strip_prefix(root)
        .unwrap_or(path)
        .to_string_lossy()
        .replace('\\', "/")
}

fn render_records_section(out: &mut String, title: &str, records: &[LabRecord]) {
    push_line(out, "");
    push_line(out, &format!("## {title}"));
    if records.is_empty() {
        push_line(out, "");
        push_line(out, "No records.");
        return;
    }
    push_line(out, "");
    for record in records {
        let label = record.label.as_deref().unwrap_or(&record.id);
        push_line(
            out,
            &format!("- `{}` ({}) - {}", record.id, record.kind, label),
        );
        if let Some(statement) = value_string(&record.data, "statement") {
            push_line(out, &format!("  Statement: {}", statement.trim()));
        } else if let Some(definition) = value_string(&record.data, "definition") {
            push_line(out, &format!("  Definition: {}", definition.trim()));
        }
        if let Some(status) = value_string(&record.data, "status") {
            push_line(out, &format!("  Status: `{status}`"));
        }
    }
}

fn render_implications_section(out: &mut String, records: &[LabRecord]) {
    push_line(out, "");
    push_line(out, "## Implications");
    if records.is_empty() {
        push_line(out, "");
        push_line(out, "No records.");
        return;
    }
    push_line(out, "");
    for record in records {
        let sources = value_string_array(&record.data, "from").join(", ");
        let target = value_string(&record.data, "to").unwrap_or_else(|| "?".to_string());
        let label = record.label.as_deref().unwrap_or(&record.id);
        push_line(
            out,
            &format!(
                "- `{}` - {}: [{}] => `{}`",
                record.id, label, sources, target
            ),
        );
        if let Some(status) = value_string(&record.data, "status") {
            push_line(out, &format!("  Status: `{status}`"));
        }
    }
}

fn render_evaluations_section(out: &mut String, records: &[LabRecord]) {
    push_line(out, "");
    push_line(out, "## Evidence");
    if records.is_empty() {
        push_line(out, "");
        push_line(out, "No records.");
        return;
    }
    push_line(out, "");
    for record in records {
        let relation = value_string(&record.data, "relation_id").unwrap_or_else(|| "?".to_string());
        let status = value_string(&record.data, "status").unwrap_or_else(|| "?".to_string());
        push_line(
            out,
            &format!("- `{}` for `{}`: `{}`", record.id, relation, status),
        );
        if let Some(checked_range) = record.data.get("checked_range") {
            push_line(out, &format!("  Checked range: {checked_range}"));
        }
        if let Some(method) = value_string(&record.data, "method") {
            push_line(out, &format!("  Method: `{method}`"));
        }
    }
}

fn render_html_records_section(out: &mut String, title: &str, records: &[LabRecord]) {
    out.push_str(&format!("<section>\n<h2>{}</h2>\n", escape_html(title)));
    if records.is_empty() {
        out.push_str("<p>No records.</p>\n</section>\n");
        return;
    }
    out.push_str("<ul>\n");
    for record in records {
        let label = record.label.as_deref().unwrap_or(&record.id);
        out.push_str(&format!(
            "  <li><code>{}</code> ({}) - {}",
            escape_html(&record.id),
            escape_html(&record.kind),
            escape_html(label)
        ));
        if let Some(statement) = value_string(&record.data, "statement") {
            out.push_str(&format!(
                "<p class=\"detail\">Statement: {}</p>",
                escape_html(statement.trim())
            ));
        } else if let Some(definition) = value_string(&record.data, "definition") {
            out.push_str(&format!(
                "<p class=\"detail\">Definition: {}</p>",
                escape_html(definition.trim())
            ));
        }
        if let Some(status) = value_string(&record.data, "status") {
            out.push_str(&format!(
                "<p class=\"detail\">Status: <code>{}</code></p>",
                escape_html(&status)
            ));
        }
        out.push_str("</li>\n");
    }
    out.push_str("</ul>\n</section>\n");
}

fn render_html_implications_section(out: &mut String, records: &[LabRecord]) {
    out.push_str("<section>\n<h2>Implications</h2>\n");
    if records.is_empty() {
        out.push_str("<p>No records.</p>\n</section>\n");
        return;
    }
    out.push_str("<ul>\n");
    for record in records {
        let sources = value_string_array(&record.data, "from")
            .into_iter()
            .map(|source| format!("<code>{}</code>", escape_html(&source)))
            .collect::<Vec<_>>()
            .join(", ");
        let target = value_string(&record.data, "to").unwrap_or_else(|| "?".to_string());
        let label = record.label.as_deref().unwrap_or(&record.id);
        out.push_str(&format!(
            "  <li><code>{}</code> - {}<p class=\"detail\">[{}] =&gt; \
             <code>{}</code></p>",
            escape_html(&record.id),
            escape_html(label),
            sources,
            escape_html(&target)
        ));
        if let Some(status) = value_string(&record.data, "status") {
            out.push_str(&format!(
                "<p class=\"detail\">Status: <code>{}</code></p>",
                escape_html(&status)
            ));
        }
        out.push_str("</li>\n");
    }
    out.push_str("</ul>\n</section>\n");
}

fn render_html_evaluations_section(out: &mut String, records: &[LabRecord]) {
    out.push_str("<section>\n<h2>Evidence</h2>\n");
    if records.is_empty() {
        out.push_str("<p>No records.</p>\n</section>\n");
        return;
    }
    out.push_str("<ul>\n");
    for record in records {
        let relation = value_string(&record.data, "relation_id").unwrap_or_else(|| "?".to_string());
        let status = value_string(&record.data, "status").unwrap_or_else(|| "?".to_string());
        out.push_str(&format!(
            "  <li><code>{}</code> for <code>{}</code>: <code>{}</code>",
            escape_html(&record.id),
            escape_html(&relation),
            escape_html(&status)
        ));
        if let Some(checked_range) = record.data.get("checked_range") {
            out.push_str(&format!(
                "<p class=\"detail\">Checked range: <code>{}</code></p>",
                escape_html(&checked_range.to_string())
            ));
        }
        if let Some(method) = value_string(&record.data, "method") {
            out.push_str(&format!(
                "<p class=\"detail\">Method: <code>{}</code></p>",
                escape_html(&method)
            ));
        }
        out.push_str("</li>\n");
    }
    out.push_str("</ul>\n</section>\n");
}

fn push_line(out: &mut String, line: &str) {
    out.push_str(line);
    out.push('\n');
}

pub fn format_validation_report(report: &ValidationReport) -> String {
    let mut out = String::new();
    push_line(&mut out, &format!("root: {}", report.root));
    push_line(&mut out, &format!("ok: {}", report.ok));
    push_line(&mut out, "counts:");
    for (kind, count) in &report.counts {
        push_line(&mut out, &format!("  {kind}: {count}"));
    }
    if !report.warnings.is_empty() {
        push_line(&mut out, "warnings:");
        for warning in &report.warnings {
            push_line(&mut out, &format!("  - {warning}"));
        }
    }
    if !report.errors.is_empty() {
        push_line(&mut out, "errors:");
        for error in &report.errors {
            push_line(&mut out, &format!("  - {error}"));
        }
    }
    out
}

pub fn format_project_overviews(projects: &[ProjectOverview]) -> String {
    projects
        .iter()
        .map(|project| {
            let label = project.label.as_deref().unwrap_or("");
            let status = project.status.as_deref().unwrap_or("");
            format!("{}\t{}\t{}\t{}", project.id, status, label, project.path)
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

pub fn format_records(records: &[LabRecord]) -> String {
    records
        .iter()
        .map(|record| {
            let label = record.label.as_deref().unwrap_or("");
            let project = record.project_id.as_deref().unwrap_or("");
            format!(
                "{}\t{}\t{}\t{}\t{}",
                record.id, record.kind, project, label, record.path
            )
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

pub fn format_trace(trace: &TraceGoalSupport) -> String {
    let mut out = String::new();
    push_line(
        &mut out,
        &format!("goal: {} / {}", trace.project_id, trace.goal_id),
    );
    if let Some(goal) = &trace.goal {
        let label = goal.label.as_deref().unwrap_or(&goal.id);
        push_line(&mut out, &format!("goal record: `{}` - {}", goal.id, label));
    }
    if !trace.direct_evaluations.is_empty() {
        push_line(&mut out, "direct evaluations:");
        for evaluation in &trace.direct_evaluations {
            push_line(&mut out, &format!("  - `{}`", evaluation.id));
        }
    }
    push_line(&mut out, "incoming implications:");
    if trace.incoming_implications.is_empty() {
        push_line(&mut out, "  none");
    }
    for implication in &trace.incoming_implications {
        push_line(
            &mut out,
            &format!(
                "  - `{}`: {}",
                implication.implication.id,
                implication
                    .implication
                    .label
                    .as_deref()
                    .unwrap_or(&implication.implication.id)
            ),
        );
        for prerequisite in &implication.prerequisites {
            let status = if prerequisite.record.is_some() {
                "indexed"
            } else {
                "missing"
            };
            push_line(
                &mut out,
                &format!("    * `{}` ({})", prerequisite.id, status),
            );
            for evaluation in &prerequisite.evaluations {
                let evidence_status =
                    value_string(&evaluation.data, "status").unwrap_or_else(|| "?".to_string());
                push_line(
                    &mut out,
                    &format!("      evidence `{}`: {}", evaluation.id, evidence_status),
                );
            }
        }
    }
    if !trace.warnings.is_empty() {
        push_line(&mut out, "warnings:");
        for warning in &trace.warnings {
            push_line(&mut out, &format!("  - {warning}"));
        }
    }
    out
}

pub fn known_ids(records: &[LabRecord]) -> BTreeSet<String> {
    records.iter().map(|record| record.id.clone()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::{Path, PathBuf};
    use std::time::{SystemTime, UNIX_EPOCH};

    fn fixture_root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/minimal_lab")
    }

    fn fixture_store() -> LabStore {
        LabStore::load(fixture_root()).expect("fixture lab root should load")
    }

    fn writable_fixture_root() -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system time")
            .as_nanos();
        let root = env::temp_dir().join(format!("polynomial-lab-fixture-{nonce}"));
        copy_directory(&fixture_root(), &root).expect("copy fixture");
        root
    }

    fn copy_directory(source: &Path, destination: &Path) -> Result<()> {
        fs::create_dir_all(destination)?;
        for entry in fs::read_dir(source)? {
            let entry = entry?;
            let source_path = entry.path();
            let destination_path = destination.join(entry.file_name());
            if source_path.is_dir() {
                copy_directory(&source_path, &destination_path)?;
            } else {
                fs::copy(&source_path, &destination_path)?;
            }
        }
        Ok(())
    }

    #[test]
    fn loads_demo_project() {
        let store = fixture_store();
        let projects = store.project_overviews();
        assert!(projects.iter().any(|project| project.id == "demo_project"));

        let report = store.project_report("demo_project");
        assert!(report.project.is_some());
        assert!(report
            .definitions
            .iter()
            .any(|record| record.id == "demo_polynomial_family"));
        assert!(report
            .goals
            .iter()
            .any(|record| record.id == "demo_real_rootedness_goal"));
    }

    #[test]
    fn validation_keeps_dangling_research_links_as_warnings() {
        let store = fixture_store();
        let report = store.validate();
        assert!(report.ok);
        assert!(report.counts.get("goal").copied().unwrap_or_default() >= 1);
        assert!(report
            .warnings
            .iter()
            .any(|warning| { warning.contains("demo_envelope_real_rooted") }));
    }

    #[test]
    fn strict_validation_promotes_dangling_links_to_errors() {
        let store = fixture_store();
        let report = store.validate_with_mode(ValidationMode::Strict);
        assert!(!report.ok);
        assert!(report
            .errors
            .iter()
            .any(|error| error.contains("demo_envelope_real_rooted")));
    }

    #[test]
    fn traces_goal_to_interlacing_evidence() {
        let store = fixture_store();
        let trace = store.trace_goal_support("demo_project", "demo_real_rootedness_goal");
        assert_eq!(trace.goal_id, "demo_real_rootedness_goal");
        assert_eq!(trace.incoming_implications.len(), 1);
        let prereq = &trace.incoming_implications[0].prerequisites[0];
        assert_eq!(prereq.id, "demo_interlaces_envelope");
        assert!(prereq
            .evaluations
            .iter()
            .any(|record| value_string(&record.data, "status").as_deref()
                == Some("holds_for_checked_domain")));
    }

    #[test]
    fn filters_evaluations_by_status_and_relation() {
        let store = fixture_store();
        let evaluations = store.evaluations(&EvaluationFilter {
            project_id: Some("demo_project".to_string()),
            status: Some("holds_for_checked_domain".to_string()),
            relation_id: Some("demo_interlaces_envelope".to_string()),
        });
        assert_eq!(evaluations.len(), 1);
    }

    #[test]
    fn renders_project_markdown() {
        let store = fixture_store();
        let markdown = store.render_project_markdown("demo_project");
        assert!(markdown.contains("# Demo interlacing project"));
        assert!(markdown.contains("## Goals"));
        assert!(markdown.contains("demo_real_rootedness_goal"));
    }

    #[test]
    fn appends_evaluation_without_overwriting_existing_files() {
        let root = writable_fixture_root();
        let store = LabStore::load(&root).expect("writable fixture should load");
        let written = store
            .append_evaluation(
                "demo_project",
                EvaluationDraft {
                    id: "demo_new_timeout".to_string(),
                    relation_id: "demo_interlaces_envelope".to_string(),
                    status: "timeout".to_string(),
                    method: Some("fixture_timeout".to_string()),
                    notes: None,
                    checked_range: Some(CheckedRange { n_min: 5, n_max: 8 }),
                    first_failure: None,
                    failure_reason: None,
                    timeout_seconds: Some(60),
                    extra: BTreeMap::new(),
                },
            )
            .expect("append evaluation");
        assert_eq!(
            written.path,
            "projects/demo_project/evidence/demo_new_timeout.json"
        );

        let reloaded = LabStore::load(&root).expect("reloaded fixture should load");
        assert_eq!(
            reloaded
                .evaluations(&EvaluationFilter {
                    project_id: Some("demo_project".to_string()),
                    status: Some("timeout".to_string()),
                    relation_id: Some("demo_interlaces_envelope".to_string()),
                })
                .len(),
            1
        );
        assert!(store
            .append_evaluation(
                "demo_project",
                EvaluationDraft {
                    id: "demo_new_timeout".to_string(),
                    relation_id: "demo_interlaces_envelope".to_string(),
                    status: "timeout".to_string(),
                    method: None,
                    notes: None,
                    checked_range: None,
                    first_failure: None,
                    failure_reason: None,
                    timeout_seconds: Some(60),
                    extra: BTreeMap::new(),
                },
            )
            .is_err());
        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn creates_project_and_appends_research_records() {
        let root = writable_fixture_root();
        let store = LabStore::load(&root).expect("writable fixture should load");
        let project = store
            .create_project(ProjectDraft {
                id: "new_project".to_string(),
                label: Some("New project".to_string()),
                status: Some("active".to_string()),
                description: Some("A writable fixture project.".to_string()),
                main_objects: Vec::new(),
                main_goals: Vec::new(),
                source_notes: Vec::new(),
            })
            .expect("create project");
        assert_eq!(project.path, "projects/new_project/project.toml");

        let reloaded = LabStore::load(&root).expect("project should reload");
        let goal = reloaded
            .append_goal(
                "new_project",
                GoalDraft {
                    id: "new_real_rootedness_goal".to_string(),
                    label: Some("New real-rootedness goal".to_string()),
                    statement: "F_n(t) is real-rooted.".to_string(),
                    status: "open".to_string(),
                    objects: vec!["future_family".to_string()],
                    motivation: None,
                    current_best_route: Some("Find an interlacing refinement.".to_string()),
                    depends_on: Vec::new(),
                },
            )
            .expect("append goal");
        assert_eq!(
            goal.path,
            "projects/new_project/goals/new_real_rootedness_goal.toml"
        );

        let reloaded = LabStore::load(&root).expect("goal should reload");
        let conjecture = reloaded
            .append_conjecture(
                "new_project",
                ConjectureDraft {
                    id: "new_interlacing_relation".to_string(),
                    label: Some("New interlacing relation".to_string()),
                    statement: "F_n(t) interlaces G_n(t).".to_string(),
                    status: "plausible".to_string(),
                    relation: Some("weak_interlaces".to_string()),
                    left: Some("future_family".to_string()),
                    right: Some("future_envelope".to_string()),
                    index_condition: Some("n >= 1".to_string()),
                    depends_on: Vec::new(),
                },
            )
            .expect("append conjecture");
        assert_eq!(
            conjecture.path,
            "projects/new_project/conjectures/new_interlacing_relation.toml"
        );

        let reloaded = LabStore::load(&root).expect("conjecture should reload");
        let implication = reloaded
            .append_implication(
                "new_project",
                ImplicationDraft {
                    id: "new_interlacing_implies_goal".to_string(),
                    label: Some("New implication".to_string()),
                    status: "plausible".to_string(),
                    from: vec!["new_interlacing_relation".to_string()],
                    to: "new_real_rootedness_goal".to_string(),
                    explanation: Some("Interlacing implies real-rootedness.".to_string()),
                    proof_tags: vec!["interlacing_implies_real_rootedness".to_string()],
                },
            )
            .expect("append implication");
        assert_eq!(
            implication.path,
            "projects/new_project/implications/new_interlacing_implies_goal.toml"
        );

        let reloaded = LabStore::load(&root).expect("records should reload");
        let trace = reloaded.trace_goal_support("new_project", "new_real_rootedness_goal");
        assert_eq!(trace.incoming_implications.len(), 1);
        assert_eq!(
            trace.incoming_implications[0].prerequisites[0].id,
            "new_interlacing_relation"
        );
        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn writes_generated_markdown_and_html() {
        let root = writable_fixture_root();
        let store = LabStore::load(&root).expect("writable fixture should load");
        let markdown = store
            .write_project_markdown("demo_project", None)
            .expect("write markdown");
        let html = store
            .write_project_html("demo_project", None)
            .expect("write html");
        assert_eq!(
            markdown.path,
            "projects/demo_project/generated/project-summary.md"
        );
        assert_eq!(
            html.path,
            "projects/demo_project/generated/project-summary.html"
        );
        assert!(root.join(&markdown.path).exists());
        assert!(root.join(&html.path).exists());
        let _ = fs::remove_dir_all(root);
    }
}
