use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use polynomial_lab::{
    default_family_registry, default_lab_root, format_project_overviews, format_records,
    format_trace, format_validation_report, interlacing_evaluation_draft,
    interlacing_evidence_id_with_offsets, real_rooted_evaluation_draft, real_rooted_evidence_id,
    CheckedRange, ComputedRefinementDraft, ConjectureDraft, EvaluationDraft, EvaluationFilter,
    FamilyIndexOffsets, GoalDraft, ImplicationDraft, InterlacingMode, LabStore, ProjectDraft,
    ValidationMode, WrittenRecord,
};
use serde::Serialize;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::str::FromStr;

#[derive(Debug, Parser)]
#[command(name = "poly-lab")]
#[command(about = "Inspect polynomial interlacing lab project data")]
struct Cli {
    #[arg(long, env = "POLY_LAB_ROOT")]
    root: Option<PathBuf>,

    #[arg(long)]
    json: bool,

    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    InitProject {
        project_id: String,
        #[arg(long)]
        label: Option<String>,
        #[arg(long, default_value = "active")]
        status: String,
        #[arg(long)]
        description: Option<String>,
        #[arg(long = "main-object")]
        main_objects: Vec<String>,
        #[arg(long = "main-goal")]
        main_goals: Vec<String>,
        #[arg(long = "source-note")]
        source_notes: Vec<String>,
    },
    Validate {
        #[arg(long)]
        strict: bool,
    },
    ListProjects,
    GetProject {
        project_id: String,
    },
    ListGoals {
        #[arg(long)]
        project: Option<String>,
    },
    ListEvaluations {
        #[arg(long)]
        project: Option<String>,
        #[arg(long)]
        status: Option<String>,
        #[arg(long)]
        relation: Option<String>,
    },
    ListExperiments {
        #[arg(long)]
        project: Option<String>,
    },
    ListProofRules,
    ListSearchRecipes,
    ListFamilies,
    ComputeFamily {
        family_id: String,
        #[arg(long)]
        n: usize,
    },
    CheckFamilyRealRooted {
        family_id: String,
        #[arg(long)]
        n_min: usize,
        #[arg(long)]
        n_max: usize,
        #[arg(long)]
        project: Option<String>,
        #[arg(long)]
        relation: Option<String>,
        #[arg(long)]
        id: Option<String>,
        #[arg(long)]
        append: bool,
    },
    CheckFamilyInterlacing {
        #[arg(long)]
        left: String,
        #[arg(long)]
        right: String,
        #[arg(long, default_value_t = 0, allow_hyphen_values = true)]
        left_offset: isize,
        #[arg(long, default_value_t = 0, allow_hyphen_values = true)]
        right_offset: isize,
        #[arg(long)]
        n_min: usize,
        #[arg(long)]
        n_max: usize,
        #[arg(long, default_value = "weak")]
        mode: String,
        #[arg(long)]
        project: Option<String>,
        #[arg(long)]
        relation: Option<String>,
        #[arg(long)]
        id: Option<String>,
        #[arg(long)]
        append: bool,
    },
    TraceGoal {
        project_id: String,
        goal_id: String,
    },
    RenderMarkdown {
        project_id: String,
        #[arg(long)]
        write: bool,
        #[arg(long)]
        output: Option<PathBuf>,
    },
    RenderHtml {
        project_id: String,
        #[arg(long)]
        write: bool,
        #[arg(long)]
        output: Option<PathBuf>,
    },
    AppendEvaluation {
        project_id: String,
        id: String,
        #[arg(long)]
        relation: String,
        #[arg(long)]
        status: String,
        #[arg(long)]
        method: Option<String>,
        #[arg(long)]
        notes: Option<String>,
        #[arg(long)]
        n_min: Option<i64>,
        #[arg(long)]
        n_max: Option<i64>,
    },
    AppendGoal {
        project_id: String,
        id: String,
        #[arg(long)]
        statement: String,
        #[arg(long)]
        label: Option<String>,
        #[arg(long, default_value = "open")]
        status: String,
        #[arg(long = "object")]
        objects: Vec<String>,
        #[arg(long)]
        motivation: Option<String>,
        #[arg(long)]
        current_best_route: Option<String>,
        #[arg(long = "depends-on")]
        depends_on: Vec<String>,
    },
    AppendConjecture {
        project_id: String,
        id: String,
        #[arg(long)]
        statement: String,
        #[arg(long)]
        label: Option<String>,
        #[arg(long, default_value = "plausible")]
        status: String,
        #[arg(long)]
        relation: Option<String>,
        #[arg(long)]
        left: Option<String>,
        #[arg(long)]
        right: Option<String>,
        #[arg(long)]
        index_condition: Option<String>,
        #[arg(long = "depends-on")]
        depends_on: Vec<String>,
    },
    AppendImplication {
        project_id: String,
        id: String,
        #[arg(long = "from", required = true)]
        from: Vec<String>,
        #[arg(long)]
        to: String,
        #[arg(long)]
        label: Option<String>,
        #[arg(long, default_value = "plausible")]
        status: String,
        #[arg(long)]
        explanation: Option<String>,
        #[arg(long = "proof-tag")]
        proof_tags: Vec<String>,
    },
    AppendComputedRefinement {
        project_id: String,
        id: String,
        #[arg(long)]
        producer: String,
        #[arg(long)]
        output_kind: String,
        #[arg(long = "index", required = true)]
        indices: Vec<String>,
        #[arg(long)]
        label: Option<String>,
        #[arg(long, default_value = "planned")]
        status: String,
        #[arg(long)]
        description: Option<String>,
        #[arg(long)]
        command: Option<String>,
        #[arg(long = "source-file")]
        source_files: Vec<String>,
        #[arg(long = "depends-on")]
        depends_on: Vec<String>,
    },
    AppendCounterexample {
        project_id: String,
        id: String,
        #[arg(long)]
        relation: String,
        #[arg(long)]
        method: Option<String>,
        #[arg(long)]
        notes: Option<String>,
        #[arg(long)]
        n: Option<i64>,
        #[arg(long)]
        first_failure_json: Option<String>,
        #[arg(long)]
        failure_reason: Option<String>,
    },
    AppendTimeout {
        project_id: String,
        id: String,
        #[arg(long)]
        relation: String,
        #[arg(long)]
        seconds: u64,
        #[arg(long)]
        method: Option<String>,
        #[arg(long)]
        notes: Option<String>,
        #[arg(long)]
        n_min: Option<i64>,
        #[arg(long)]
        n_max: Option<i64>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let root = cli.root.unwrap_or_else(default_lab_root);
    let store = LabStore::load(&root)
        .with_context(|| format!("failed to load polynomial lab at {}", root.display()))?;

    match cli.command {
        Command::InitProject {
            project_id,
            label,
            status,
            description,
            main_objects,
            main_goals,
            source_notes,
        } => {
            let written = store.create_project(ProjectDraft {
                id: project_id,
                label,
                status: Some(status),
                description,
                main_objects,
                main_goals,
                source_notes,
            })?;
            print_written_record(cli.json, &written)?;
        }
        Command::Validate { strict } => {
            let mode = if strict {
                ValidationMode::Strict
            } else {
                ValidationMode::Tolerant
            };
            let report = store.validate_with_mode(mode);
            if cli.json {
                print_json(&report)?;
            } else {
                print!("{}", format_validation_report(&report));
            }
            if !report.ok {
                std::process::exit(1);
            }
        }
        Command::ListProjects => {
            let projects = store.project_overviews();
            if cli.json {
                print_json(&projects)?;
            } else {
                print!("{}", format_project_overviews(&projects));
            }
        }
        Command::GetProject { project_id } => {
            let report = store.project_report(&project_id);
            if cli.json {
                print_json(&report)?;
            } else {
                print!("{}", store.render_project_markdown(&project_id));
            }
        }
        Command::ListGoals { project } => {
            let goals = store.goals(project.as_deref());
            if cli.json {
                print_json(&goals)?;
            } else {
                print!("{}", format_records(&goals));
            }
        }
        Command::ListEvaluations {
            project,
            status,
            relation,
        } => {
            let evaluations = store.evaluations(&EvaluationFilter {
                project_id: project,
                status,
                relation_id: relation,
            });
            if cli.json {
                print_json(&evaluations)?;
            } else {
                print!("{}", format_records(&evaluations));
            }
        }
        Command::ListExperiments { project } => {
            let experiments = store.experiments(project.as_deref());
            if cli.json {
                print_json(&experiments)?;
            } else {
                print!("{}", format_records(&experiments));
            }
        }
        Command::ListProofRules => {
            let rules = store.proof_rules();
            if cli.json {
                print_json(&rules)?;
            } else {
                print!("{}", format_records(&rules));
            }
        }
        Command::ListSearchRecipes => {
            let recipes = store.search_recipes();
            if cli.json {
                print_json(&recipes)?;
            } else {
                print!("{}", format_records(&recipes));
            }
        }
        Command::ListFamilies => {
            let registry = default_family_registry();
            let families = registry.list();
            if cli.json {
                print_json(&families)?;
            } else {
                for family in families {
                    println!(
                        "{}\t{}\t{}\tn >= {}\t{}",
                        family.id, family.symbol, family.label, family.min_n, family.source
                    );
                }
            }
        }
        Command::ComputeFamily { family_id, n } => {
            let registry = default_family_registry();
            let computed = registry.compute(&family_id, n)?;
            if cli.json {
                print_json(&computed)?;
            } else {
                println!("{} at n={}", computed.family_id, computed.n);
                println!("coefficients: [{}]", computed.coefficients.join(", "));
                println!("polynomial: {}", computed.polynomial);
            }
        }
        Command::CheckFamilyRealRooted {
            family_id,
            n_min,
            n_max,
            project,
            relation,
            id,
            append,
        } => {
            let registry = default_family_registry();
            let report = registry.check_real_rooted(&family_id, n_min, n_max)?;
            let written = if append {
                let project = project
                    .as_deref()
                    .ok_or_else(|| anyhow::anyhow!("--append requires --project"))?;
                let relation = relation
                    .as_deref()
                    .ok_or_else(|| anyhow::anyhow!("--append requires --relation"))?;
                let evidence_id = id.unwrap_or_else(|| {
                    real_rooted_evidence_id(relation, report.first_failure_n, n_min, n_max)
                });
                let draft = real_rooted_evaluation_draft(
                    evidence_id,
                    relation.to_string(),
                    &family_id,
                    &report,
                )?;
                Some(store.append_evaluation(project, draft)?)
            } else {
                None
            };

            if cli.json {
                print_json(&json!({
                    "report": report,
                    "written_evidence": written,
                }))?;
            } else {
                println!(
                    "{} real-rooted for n={}..{}: {}",
                    report.family_id, report.n_min, report.n_max, report.all_real_rooted
                );
                if let Some(first_failure_n) = report.first_failure_n {
                    println!("first failure: n={first_failure_n}");
                }
                if let Some(written) = written {
                    println!("wrote {}", written.path);
                }
            }
        }
        Command::CheckFamilyInterlacing {
            left,
            right,
            left_offset,
            right_offset,
            n_min,
            n_max,
            mode,
            project,
            relation,
            id,
            append,
        } => {
            let mode = InterlacingMode::from_str(&mode)?;
            let registry = default_family_registry();
            let report = registry.check_interlacing_with_offsets(
                &left,
                &right,
                n_min,
                n_max,
                FamilyIndexOffsets {
                    left: left_offset,
                    right: right_offset,
                },
                mode,
            )?;
            let written = if append {
                let project = project
                    .as_deref()
                    .ok_or_else(|| anyhow::anyhow!("--append requires --project"))?;
                let relation = relation
                    .as_deref()
                    .ok_or_else(|| anyhow::anyhow!("--append requires --relation"))?;
                let evidence_id = id.unwrap_or_else(|| {
                    interlacing_evidence_id_with_offsets(
                        relation,
                        mode,
                        report.first_failure_n,
                        n_min,
                        n_max,
                        left_offset,
                        right_offset,
                    )
                });
                let draft =
                    interlacing_evaluation_draft(evidence_id, relation.to_string(), &report)?;
                Some(store.append_evaluation(project, draft)?)
            } else {
                None
            };

            if cli.json {
                print_json(&json!({
                    "report": report,
                    "written_evidence": written,
                }))?;
            } else {
                println!(
                    "{} {}-interlaces {} for n={}..{}: {}",
                    cli_indexed_family(&report.left_family_id, report.left_offset),
                    report.mode,
                    cli_indexed_family(&report.right_family_id, report.right_offset),
                    report.n_min,
                    report.n_max,
                    report.all_interlacing
                );
                if let Some(first_failure_n) = report.first_failure_n {
                    println!("first failure: n={first_failure_n}");
                }
                if let Some(written) = written {
                    println!("wrote {}", written.path);
                }
            }
        }
        Command::TraceGoal {
            project_id,
            goal_id,
        } => {
            let trace = store.trace_goal_support(&project_id, &goal_id);
            if cli.json {
                print_json(&trace)?;
            } else {
                print!("{}", format_trace(&trace));
            }
        }
        Command::RenderMarkdown {
            project_id,
            write,
            output,
        } => {
            if write || output.is_some() {
                let generated = store.write_project_markdown(&project_id, output.as_deref())?;
                if cli.json {
                    print_json(&generated)?;
                } else {
                    println!("wrote {}", generated.path);
                }
            } else {
                print!("{}", store.render_project_markdown(&project_id));
            }
        }
        Command::RenderHtml {
            project_id,
            write,
            output,
        } => {
            if write || output.is_some() {
                let generated = store.write_project_html(&project_id, output.as_deref())?;
                if cli.json {
                    print_json(&generated)?;
                } else {
                    println!("wrote {}", generated.path);
                }
            } else {
                print!("{}", store.render_project_html(&project_id));
            }
        }
        Command::AppendEvaluation {
            project_id,
            id,
            relation,
            status,
            method,
            notes,
            n_min,
            n_max,
        } => {
            let written = store.append_evaluation(
                &project_id,
                EvaluationDraft {
                    id,
                    relation_id: relation,
                    status,
                    method,
                    notes,
                    checked_range: checked_range(n_min, n_max)?,
                    first_failure: None,
                    failure_reason: None,
                    timeout_seconds: None,
                    extra: BTreeMap::new(),
                },
            )?;
            print_written_evaluation(cli.json, &written)?;
        }
        Command::AppendGoal {
            project_id,
            id,
            statement,
            label,
            status,
            objects,
            motivation,
            current_best_route,
            depends_on,
        } => {
            let written = store.append_goal(
                &project_id,
                GoalDraft {
                    id,
                    label,
                    statement,
                    status,
                    objects,
                    motivation,
                    current_best_route,
                    depends_on,
                },
            )?;
            print_written_record(cli.json, &written)?;
        }
        Command::AppendConjecture {
            project_id,
            id,
            statement,
            label,
            status,
            relation,
            left,
            right,
            index_condition,
            depends_on,
        } => {
            let written = store.append_conjecture(
                &project_id,
                ConjectureDraft {
                    id,
                    label,
                    statement,
                    status,
                    relation,
                    left,
                    right,
                    index_condition,
                    depends_on,
                },
            )?;
            print_written_record(cli.json, &written)?;
        }
        Command::AppendImplication {
            project_id,
            id,
            from,
            to,
            label,
            status,
            explanation,
            proof_tags,
        } => {
            let written = store.append_implication(
                &project_id,
                ImplicationDraft {
                    id,
                    label,
                    status,
                    from,
                    to,
                    explanation,
                    proof_tags,
                },
            )?;
            print_written_record(cli.json, &written)?;
        }
        Command::AppendComputedRefinement {
            project_id,
            id,
            producer,
            output_kind,
            indices,
            label,
            status,
            description,
            command,
            source_files,
            depends_on,
        } => {
            let written = store.append_computed_refinement(
                &project_id,
                ComputedRefinementDraft {
                    id,
                    label,
                    status,
                    producer,
                    output_kind,
                    indices,
                    description,
                    command,
                    source_files,
                    depends_on,
                },
            )?;
            print_written_record(cli.json, &written)?;
        }
        Command::AppendCounterexample {
            project_id,
            id,
            relation,
            method,
            notes,
            n,
            first_failure_json,
            failure_reason,
        } => {
            let first_failure = first_failure_value(n, first_failure_json)?;
            let written = store.append_evaluation(
                &project_id,
                EvaluationDraft {
                    id,
                    relation_id: relation,
                    status: "counterexample_found".to_string(),
                    method,
                    notes,
                    checked_range: None,
                    first_failure: Some(first_failure),
                    failure_reason,
                    timeout_seconds: None,
                    extra: BTreeMap::new(),
                },
            )?;
            print_written_evaluation(cli.json, &written)?;
        }
        Command::AppendTimeout {
            project_id,
            id,
            relation,
            seconds,
            method,
            notes,
            n_min,
            n_max,
        } => {
            let written = store.append_evaluation(
                &project_id,
                EvaluationDraft {
                    id,
                    relation_id: relation,
                    status: "timeout".to_string(),
                    method,
                    notes,
                    checked_range: checked_range(n_min, n_max)?,
                    first_failure: None,
                    failure_reason: None,
                    timeout_seconds: Some(seconds),
                    extra: BTreeMap::new(),
                },
            )?;
            print_written_evaluation(cli.json, &written)?;
        }
    }
    Ok(())
}

fn checked_range(n_min: Option<i64>, n_max: Option<i64>) -> Result<Option<CheckedRange>> {
    match (n_min, n_max) {
        (Some(n_min), Some(n_max)) if n_min <= n_max => Ok(Some(CheckedRange { n_min, n_max })),
        (Some(_), Some(_)) => anyhow::bail!("expected --n-min <= --n-max"),
        (None, None) => Ok(None),
        _ => anyhow::bail!("expected both --n-min and --n-max, or neither"),
    }
}

fn first_failure_value(n: Option<i64>, first_failure_json: Option<String>) -> Result<Value> {
    match (n, first_failure_json) {
        (Some(n), None) => Ok(serde_json::json!({ "n": n })),
        (None, Some(text)) => serde_json::from_str(&text)
            .with_context(|| "failed to parse --first-failure-json as JSON"),
        (Some(_), Some(_)) => anyhow::bail!("use either --n or --first-failure-json, not both"),
        (None, None) => anyhow::bail!("counterexample records require --n or --first-failure-json"),
    }
}

fn cli_indexed_family(family_id: &str, offset: isize) -> String {
    match offset.cmp(&0) {
        std::cmp::Ordering::Greater => format!("{family_id}[n+{offset}]"),
        std::cmp::Ordering::Less => format!("{family_id}[n-{}]", offset.unsigned_abs()),
        std::cmp::Ordering::Equal => format!("{family_id}[n]"),
    }
}

fn print_written_evaluation(json: bool, written: &polynomial_lab::WrittenEvaluation) -> Result<()> {
    if json {
        print_json(written)
    } else {
        println!("wrote {}", written.path);
        Ok(())
    }
}

fn print_written_record(json: bool, written: &WrittenRecord) -> Result<()> {
    if json {
        print_json(written)
    } else {
        println!("wrote {}", written.path);
        Ok(())
    }
}

fn print_json(value: &impl Serialize) -> Result<()> {
    println!("{}", serde_json::to_string_pretty(value)?);
    Ok(())
}
