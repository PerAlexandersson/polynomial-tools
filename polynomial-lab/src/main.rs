use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use polynomial_lab::{
    default_lab_root, format_project_overviews, format_records, format_trace,
    format_validation_report, EvaluationFilter, LabStore,
};
use serde::Serialize;
use std::path::PathBuf;

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
    Validate,
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
    ListProofRules,
    ListSearchRecipes,
    TraceGoal {
        project_id: String,
        goal_id: String,
    },
    RenderMarkdown {
        project_id: String,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let root = cli.root.unwrap_or_else(default_lab_root);
    let store = LabStore::load(&root)
        .with_context(|| format!("failed to load polynomial lab at {}", root.display()))?;

    match cli.command {
        Command::Validate => {
            let report = store.validate();
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
        Command::RenderMarkdown { project_id } => {
            print!("{}", store.render_project_markdown(&project_id));
        }
    }
    Ok(())
}

fn print_json(value: &impl Serialize) -> Result<()> {
    println!("{}", serde_json::to_string_pretty(value)?);
    Ok(())
}
