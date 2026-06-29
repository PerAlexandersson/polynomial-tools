use std::fs;
use std::io::Write;
use std::process::{Command, Output, Stdio};

fn run_polytool(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn polytool");

    child
        .stdin
        .as_mut()
        .expect("polytool stdin")
        .write_all(input.as_bytes())
        .expect("write polytool stdin");

    child.wait_with_output().expect("wait for polytool")
}

#[test]
fn real_rooted_accepts_bigint_coefficients() {
    let output = run_polytool(&["real-rooted"], "-1000000000000000000000000000000,1\n");

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    assert!(stdout.contains("-1000000000000000000000000000000 + t: real-rooted"));
}

#[test]
fn gamma_expansion_json_accepts_bigint_coefficients() {
    let output = run_polytool(
        &["gamma-expansion", "--json"],
        "1,100000000000000000000,1\n",
    );

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    assert!(stdout.contains("\"coefficients\":[\"1\",\"100000000000000000000\",\"1\"]"));
    assert!(stdout.contains("\"gamma\":[\"1\",\"99999999999999999998\"]"));
    assert!(stdout.contains("\"expansion\":\"(1+t)^2 + 99999999999999999998 t\""));
}

#[test]
fn sequence_json_generates_bigint_coefficients() {
    let output = run_polytool(&["sequence", "chebyshev-t", "64", "--json"], "");

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    assert!(stdout.contains("\"degree\":64"));
    assert!(stdout.contains("\"9223372036854775808\""));
}

#[test]
fn resultant_and_discriminant_accept_bigint_coefficients() {
    let resultant = run_polytool(
        &["resultant"],
        "-100000000000000000000,1\n-100000000000000000001,1\n",
    );
    assert!(resultant.status.success());
    let stdout = String::from_utf8(resultant.stdout).expect("stdout is utf8");
    assert!(stdout.contains("= -1"));

    let discriminant = run_polytool(&["discriminant"], "1,0,1000000000000000000000000000000\n");
    assert!(discriminant.status.success());
    let stdout = String::from_utf8(discriminant.stdout).expect("stdout is utf8");
    assert!(stdout.contains("-4000000000000000000000000000000"));
}

#[test]
fn bench_interlacing_reports_tsv() {
    let output = run_polytool(
        &[
            "bench",
            "interlacing",
            "--sequence",
            "eulerian",
            "--max-n",
            "4",
            "--repeat",
            "1",
        ],
        "",
    );

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    assert!(
        stdout.starts_with("sequence\tleft_index\tright_index\tdegree\trepeat\tavg_us\tresult\n")
    );
    assert!(stdout.contains("eulerian\t"));
}

#[test]
fn bench_recurrence_fixture_reports_tsv() {
    let report_path = std::env::temp_dir().join(format!(
        "polytool-recurrence-bench-report-{}.md",
        std::process::id()
    ));
    let report_arg = report_path.to_string_lossy().into_owned();
    let output = run_polytool(
        &[
            "bench",
            "recurrence-fixtures",
            "--only",
            "01_scalar_geometric",
            "--repeat",
            "1",
            "--summary",
            "--report",
            &report_arg,
        ],
        "",
    );

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).expect("stdout is utf8");
    assert!(stdout.starts_with(
        "slug\trun\tfound\telapsed_ms\tcandidates\tunknowns\tweighted\tfit_rows\tverify_rows\trecurrence\n"
    ));
    assert!(stdout.contains("01_scalar_geometric\t1\ttrue\t"));
    assert!(stdout.contains("P(n) = 2 P(n-1)"));
    assert!(stdout.contains("# fixture_summary\n"));
    assert!(stdout.contains("# category_summary\n"));

    let report = fs::read_to_string(&report_path).expect("read benchmark report");
    assert!(report.contains("# Recurrence Fixture Benchmark Report"));
    assert!(report.contains("| synthetic | 1 | 1 | 1 |"));
    assert!(report.contains("| 01_scalar_geometric |"));
    let _ = fs::remove_file(report_path);
}
