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
