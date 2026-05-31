use polynomial_lab::default_lab_root;
use polynomial_lab_mcp::PolynomialLabServer;
use rmcp::{transport::stdio, ServiceExt};
use std::path::PathBuf;

const BIN_NAME: &str = "poly-lab-mcp";
const CONTACT: &str = "Per Alexandersson <per.w.alexandersson@gmail.com>";

#[derive(Debug, Eq, PartialEq)]
enum StartupMode {
    Serve { root: PathBuf },
    Help,
    Version,
}

fn parse_args(args: &[String]) -> Result<StartupMode, String> {
    match args {
        [_program] => Ok(StartupMode::Serve {
            root: default_lab_root(),
        }),
        [_program, flag] if flag == "--help" || flag == "-h" => Ok(StartupMode::Help),
        [_program, flag] if flag == "--version" || flag == "-V" => Ok(StartupMode::Version),
        [_program, flag, root] if flag == "--root" => Ok(StartupMode::Serve {
            root: PathBuf::from(root),
        }),
        [_program, flag] if flag == "--root" => Err("--root requires a path".to_string()),
        [_program, flag] => Err(format!("unknown option: {flag}")),
        [_program, flag, ..] => Err(format!("unexpected extra arguments after: {flag}")),
        [] => Ok(StartupMode::Serve {
            root: default_lab_root(),
        }),
    }
}

fn version_text() -> String {
    format!("{} {}", BIN_NAME, env!("CARGO_PKG_VERSION"))
}

fn help_text() -> String {
    format!(
        "\
{name} {version}
{description}

Usage:
  {name}
  {name} --root PATH
  {name} --help
  {name} --version

Options:
      --root PATH    Lab root. Defaults to POLY_LAB_ROOT or /workspace/projects/polynomial-interlacing-lab
  -h, --help         Print this help text and exit
  -V, --version      Print version information and exit

MCP usage:
  Run in serve mode from an MCP client. The server uses stdio transport and
  exposes project-index tools only.

Contact:
  {contact}
",
        name = BIN_NAME,
        version = env!("CARGO_PKG_VERSION"),
        description = env!("CARGO_PKG_DESCRIPTION"),
        contact = CONTACT,
    )
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    match parse_args(&args) {
        Ok(StartupMode::Help) => {
            print!("{}", help_text());
            return Ok(());
        }
        Ok(StartupMode::Version) => {
            println!("{}", version_text());
            return Ok(());
        }
        Ok(StartupMode::Serve { root }) => {
            let service = PolynomialLabServer::with_root(root).serve(stdio()).await?;
            service.waiting().await?;
        }
        Err(message) => {
            eprintln!("error: {message}");
            eprintln!("Try '{BIN_NAME} --help' for usage.");
            std::process::exit(2);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| value.to_string()).collect()
    }

    #[test]
    fn parses_cli_flags() {
        assert!(matches!(
            parse_args(&args(&["poly-lab-mcp"])),
            Ok(StartupMode::Serve { .. })
        ));
        assert_eq!(
            parse_args(&args(&["poly-lab-mcp", "--help"])),
            Ok(StartupMode::Help)
        );
        assert_eq!(
            parse_args(&args(&["poly-lab-mcp", "-h"])),
            Ok(StartupMode::Help)
        );
        assert_eq!(
            parse_args(&args(&["poly-lab-mcp", "--version"])),
            Ok(StartupMode::Version)
        );
        assert_eq!(
            parse_args(&args(&["poly-lab-mcp", "--root", "/tmp/lab"])),
            Ok(StartupMode::Serve {
                root: PathBuf::from("/tmp/lab")
            })
        );
        assert!(parse_args(&args(&["poly-lab-mcp", "--bad"])).is_err());
        assert!(parse_args(&args(&["poly-lab-mcp", "--root"])).is_err());
    }

    #[test]
    fn help_mentions_stdio_and_root() {
        let help = help_text();
        assert!(help.contains("stdio transport"));
        assert!(help.contains("--root PATH"));
        assert!(help.contains(CONTACT));
    }
}
