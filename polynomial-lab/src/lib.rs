use anyhow::{Context, Result};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

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

    pub fn validate(&self) -> ValidationReport {
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
                warnings.push(format!("duplicate id '{id}' appears in {paths}"));
            }
        }

        for project in self
            .records
            .iter()
            .filter(|record| record.kind == "project")
        {
            for goal_id in value_string_array(&project.data, "main_goals") {
                if self.record(&goal_id).is_none() {
                    warnings.push(format!(
                        "project '{}' references missing main goal '{}'",
                        project.id, goal_id
                    ));
                }
            }
            for object_id in value_string_array(&project.data, "main_objects") {
                if self.record(&object_id).is_none() {
                    warnings.push(format!(
                        "project '{}' references missing main object '{}'",
                        project.id, object_id
                    ));
                }
            }
            for source_note in value_string_array(&project.data, "source_notes") {
                if !Path::new(&source_note).exists() {
                    warnings.push(format!(
                        "project '{}' source note does not exist: {}",
                        project.id, source_note
                    ));
                }
            }
        }

        for record in &self.records {
            for dependency_id in value_string_array(&record.data, "depends_on") {
                if self.record(&dependency_id).is_none() {
                    warnings.push(format!(
                        "record '{}' depends on missing id '{}'",
                        record.id, dependency_id
                    ));
                }
            }

            if record.kind == "goal" {
                for object_id in value_string_array(&record.data, "objects") {
                    if self.record(&object_id).is_none() {
                        warnings.push(format!(
                            "goal '{}' references missing object '{}'",
                            record.id, object_id
                        ));
                    }
                }
            }

            if record.kind == "conjecture" {
                for key in ["left", "right"] {
                    if let Some(id) = value_string(&record.data, key) {
                        if self.record(&id).is_none() {
                            warnings.push(format!(
                                "conjecture '{}' references missing {} '{}'",
                                record.id, key, id
                            ));
                        }
                    }
                }
            }

            if record.kind == "implication" {
                for source_id in value_string_array(&record.data, "from") {
                    if self.record(&source_id).is_none() {
                        warnings.push(format!(
                            "implication '{}' has missing prerequisite '{}'",
                            record.id, source_id
                        ));
                    }
                }
                if let Some(target_id) = value_string(&record.data, "to") {
                    if self.record(&target_id).is_none() {
                        warnings.push(format!(
                            "implication '{}' has missing target '{}'",
                            record.id, target_id
                        ));
                    }
                } else {
                    errors.push(format!("implication '{}' is missing field 'to'", record.id));
                }
            }

            if record.kind == "evaluation" {
                if let Some(relation_id) = value_string(&record.data, "relation_id") {
                    if self.record(&relation_id).is_none() {
                        warnings.push(format!(
                            "evaluation '{}' references missing relation '{}'",
                            record.id, relation_id
                        ));
                    }
                } else {
                    warnings.push(format!("evaluation '{}' has no relation_id", record.id));
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

    fn fixture_store() -> LabStore {
        LabStore::load(DEFAULT_LAB_ROOT).expect("fixture lab root should load")
    }

    #[test]
    fn loads_derangement_descents_project() {
        let store = fixture_store();
        let projects = store.project_overviews();
        assert!(projects
            .iter()
            .any(|project| project.id == "derangement_descents"));

        let report = store.project_report("derangement_descents");
        assert!(report.project.is_some());
        assert!(report
            .definitions
            .iter()
            .any(|record| record.id == "normalized_derangement_descent_polynomial"));
        assert!(report
            .goals
            .iter()
            .any(|record| record.id == "derangement_descent_real_rootedness"));
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
            .any(|warning| { warning.contains("reciprocal_eulerian_derivative_real_rooted") }));
    }

    #[test]
    fn traces_goal_to_interlacing_evidence() {
        let store = fixture_store();
        let trace = store.trace_goal_support(
            "derangement_descents",
            "derangement_descent_real_rootedness",
        );
        assert_eq!(trace.goal_id, "derangement_descent_real_rootedness");
        assert_eq!(trace.incoming_implications.len(), 1);
        let prereq = &trace.incoming_implications[0].prerequisites[0];
        assert_eq!(
            prereq.id,
            "normalized_derangement_descent_interlaces_reciprocal_eulerian_derivative"
        );
        assert!(prereq
            .evaluations
            .iter()
            .any(
                |record| value_string(&record.data, "status").as_deref() == Some("verified_range")
            ));
    }

    #[test]
    fn filters_evaluations_by_status_and_relation() {
        let store = fixture_store();
        let evaluations = store.evaluations(&EvaluationFilter {
            project_id: Some("derangement_descents".to_string()),
            status: Some("verified_range".to_string()),
            relation_id: Some(
                "normalized_derangement_descent_interlaces_reciprocal_eulerian_derivative"
                    .to_string(),
            ),
        });
        assert_eq!(evaluations.len(), 1);
    }

    #[test]
    fn renders_project_markdown() {
        let store = fixture_store();
        let markdown = store.render_project_markdown("derangement_descents");
        assert!(markdown.contains("# Derangement descent real-rootedness"));
        assert!(markdown.contains("## Goals"));
        assert!(markdown.contains("derangement_descent_real_rootedness"));
    }
}
