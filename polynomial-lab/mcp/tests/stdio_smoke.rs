use serde_json::{json, Value};
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

fn fixture_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../tests/fixtures/minimal_lab")
}

fn writable_fixture_root() -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time")
        .as_nanos();
    let root = std::env::temp_dir().join(format!("polynomial-lab-mcp-fixture-{nonce}"));
    copy_directory(&fixture_root(), &root).expect("copy fixture");
    root
}

fn copy_directory(source: &Path, destination: &Path) -> std::io::Result<()> {
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

fn start_server(root: &Path) -> (Child, ChildStdin, BufReader<std::process::ChildStdout>) {
    let mut child = Command::new(env!("CARGO_BIN_EXE_poly-lab-mcp"))
        .env("POLY_LAB_ROOT", root)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("spawn poly-lab-mcp");
    let stdin = child.stdin.take().expect("child stdin");
    let stdout = BufReader::new(child.stdout.take().expect("child stdout"));
    (child, stdin, stdout)
}

fn send(stdin: &mut ChildStdin, value: Value) {
    writeln!(stdin, "{}", serde_json::to_string(&value).unwrap()).expect("write request");
    stdin.flush().expect("flush request");
}

fn read_response(stdout: &mut BufReader<std::process::ChildStdout>, id: i64) -> Value {
    let mut line = String::new();
    loop {
        line.clear();
        let n = stdout.read_line(&mut line).expect("read response");
        assert!(n > 0, "server closed stdout before response {id}");
        let value: Value = serde_json::from_str(line.trim_end()).expect("json response");
        if value.get("id").and_then(Value::as_i64) == Some(id) {
            return value;
        }
    }
}

#[test]
fn lists_tools_and_traces_derangement_goal() {
    let root = writable_fixture_root();
    let (mut child, mut stdin, mut stdout) = start_server(&root);

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 1,
            "method": "initialize",
            "params": {
                "protocolVersion": "2025-03-26",
                "capabilities": {},
                "clientInfo": { "name": "polynomial-lab-mcp-test", "version": "0.1.0" }
            }
        }),
    );
    let init = read_response(&mut stdout, 1);
    assert!(init.get("result").is_some(), "initialize failed: {init}");

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "method": "notifications/initialized"
        }),
    );

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 2,
            "method": "tools/list",
            "params": {}
        }),
    );
    let tools = read_response(&mut stdout, 2);
    let tool_names: Vec<&str> = tools["result"]["tools"]
        .as_array()
        .expect("tool list")
        .iter()
        .filter_map(|tool| tool["name"].as_str())
        .collect();
    assert!(tool_names.contains(&"list_projects"));
    assert!(tool_names.contains(&"compute_family"));
    assert!(tool_names.contains(&"check_family_real_rooted"));
    assert!(tool_names.contains(&"check_family_interlacing"));
    assert!(tool_names.contains(&"create_project"));
    assert!(tool_names.contains(&"append_goal"));
    assert!(tool_names.contains(&"append_conjecture"));
    assert!(tool_names.contains(&"append_implication"));
    assert!(tool_names.contains(&"append_computed_refinement"));
    assert!(tool_names.contains(&"list_experiments"));
    assert!(tool_names.contains(&"trace_goal_support"));
    assert!(tool_names.contains(&"render_project_markdown"));
    assert!(tool_names.contains(&"append_timeout"));
    assert!(tool_names.contains(&"write_project_html"));

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 3,
            "method": "tools/call",
            "params": {
                "name": "trace_goal_support",
                "arguments": {
                    "project_id": "demo_project",
                    "goal_id": "demo_real_rootedness_goal"
                }
            }
        }),
    );
    let result = read_response(&mut stdout, 3);
    let structured = &result["result"]["structuredContent"];
    assert_eq!(structured["project_id"], "demo_project");
    assert_eq!(structured["goal_id"], "demo_real_rootedness_goal");
    assert_eq!(
        structured["incoming_implications"][0]["implication"]["id"],
        "demo_interlacing_implies_real_rootedness"
    );

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 4,
            "method": "tools/call",
            "params": {
                "name": "compute_family",
                "arguments": {
                    "family_id": "derangement_descent_polynomial",
                    "n": 4
                }
            }
        }),
    );
    let family = read_response(&mut stdout, 4);
    assert_eq!(
        family["result"]["structuredContent"]["coefficients"],
        json!(["0", "4", "4", "1"])
    );

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 5,
            "method": "tools/call",
            "params": {
                "name": "check_family_real_rooted",
                "arguments": {
                    "family_id": "derangement_descent_polynomial",
                    "n_min": 2,
                    "n_max": 6,
                    "append": true,
                    "project_id": "demo_project",
                    "relation_id": "demo_derangement_real_rooted"
                }
            }
        }),
    );
    let check = read_response(&mut stdout, 5);
    assert_eq!(
        check["result"]["structuredContent"]["report"]["all_real_rooted"],
        true
    );
    assert_eq!(
        check["result"]["structuredContent"]["written_evidence"]["path"],
        "projects/demo_project/evidence/demo_derangement_real_rooted_n_2_6.json"
    );

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 6,
            "method": "tools/call",
            "params": {
                "name": "create_project",
                "arguments": {
                    "project_id": "mcp_created_project",
                    "label": "MCP-created project",
                    "description": "Created by the stdio smoke test."
                }
            }
        }),
    );
    let created = read_response(&mut stdout, 6);
    assert_eq!(
        created["result"]["structuredContent"]["path"],
        "projects/mcp_created_project/project.toml"
    );

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 7,
            "method": "tools/call",
            "params": {
                "name": "append_goal",
                "arguments": {
                    "project_id": "mcp_created_project",
                    "id": "mcp_goal",
                    "statement": "F_n(t) is real-rooted.",
                    "current_best_route": "Find an interlacing refinement."
                }
            }
        }),
    );
    let goal = read_response(&mut stdout, 7);
    assert_eq!(
        goal["result"]["structuredContent"]["path"],
        "projects/mcp_created_project/goals/mcp_goal.toml"
    );

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 8,
            "method": "tools/call",
            "params": {
                "name": "append_conjecture",
                "arguments": {
                    "project_id": "mcp_created_project",
                    "id": "mcp_relation",
                    "statement": "F_n(t) interlaces G_n(t).",
                    "relation": "weak_interlaces",
                    "left": "mcp_family",
                    "right": "mcp_envelope"
                }
            }
        }),
    );
    let relation = read_response(&mut stdout, 8);
    assert_eq!(
        relation["result"]["structuredContent"]["path"],
        "projects/mcp_created_project/conjectures/mcp_relation.toml"
    );

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 9,
            "method": "tools/call",
            "params": {
                "name": "append_implication",
                "arguments": {
                    "project_id": "mcp_created_project",
                    "id": "mcp_relation_implies_goal",
                    "from": ["mcp_relation"],
                    "to": "mcp_goal",
                    "explanation": "Interlacing implies real-rootedness."
                }
            }
        }),
    );
    let implication = read_response(&mut stdout, 9);
    assert_eq!(
        implication["result"]["structuredContent"]["path"],
        "projects/mcp_created_project/implications/mcp_relation_implies_goal.toml"
    );

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 10,
            "method": "tools/call",
            "params": {
                "name": "append_computed_refinement",
                "arguments": {
                    "project_id": "mcp_created_project",
                    "id": "mcp_refinement",
                    "producer": "cargo run -p experiments --bin mcp_refinement",
                    "output_kind": "polynomial_vector_by_n",
                    "indices": ["n", "j"],
                    "description": "Fixture computed refinement.",
                    "depends_on": ["mcp_relation"]
                }
            }
        }),
    );
    let refinement = read_response(&mut stdout, 10);
    assert_eq!(
        refinement["result"]["structuredContent"]["path"],
        "projects/mcp_created_project/experiments/mcp_refinement.toml"
    );

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 11,
            "method": "tools/call",
            "params": {
                "name": "list_experiments",
                "arguments": {
                    "project_id": "mcp_created_project"
                }
            }
        }),
    );
    let experiments = read_response(&mut stdout, 11);
    assert_eq!(
        experiments["result"]["structuredContent"]["records"][0]["id"],
        "mcp_refinement"
    );

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 12,
            "method": "tools/call",
            "params": {
                "name": "append_timeout",
                "arguments": {
                    "project_id": "demo_project",
                    "id": "demo_mcp_timeout",
                    "relation_id": "demo_interlaces_envelope",
                    "seconds": 30,
                    "method": "stdio_smoke",
                    "checked_range": { "n_min": 5, "n_max": 7 }
                }
            }
        }),
    );
    let append = read_response(&mut stdout, 12);
    assert_eq!(
        append["result"]["structuredContent"]["path"],
        "projects/demo_project/evidence/demo_mcp_timeout.json"
    );
    assert!(root
        .join("projects/demo_project/evidence/demo_mcp_timeout.json")
        .exists());

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 13,
            "method": "tools/call",
            "params": {
                "name": "write_project_html",
                "arguments": {
                    "project_id": "demo_project"
                }
            }
        }),
    );
    let write = read_response(&mut stdout, 13);
    assert_eq!(
        write["result"]["structuredContent"]["path"],
        "projects/demo_project/generated/project-summary.html"
    );
    assert!(root
        .join("projects/demo_project/generated/project-summary.html")
        .exists());

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 14,
            "method": "tools/call",
            "params": {
                "name": "check_family_interlacing",
                "arguments": {
                    "left_family_id": "normalized_derangement_descent_polynomial",
                    "right_family_id": "reciprocal_eulerian_derivative_polynomial",
                    "n_min": 5,
                    "n_max": 6,
                    "mode": "weak",
                    "append": true,
                    "project_id": "demo_project",
                    "relation_id": "demo_interlaces_envelope"
                }
            }
        }),
    );
    let interlacing = read_response(&mut stdout, 14);
    assert_eq!(
        interlacing["result"]["structuredContent"]["report"]["all_interlacing"],
        true
    );
    assert_eq!(
        interlacing["result"]["structuredContent"]["written_evidence"]["path"],
        "projects/demo_project/evidence/demo_interlaces_envelope_weak_n_5_6.json"
    );
    assert!(root
        .join("projects/demo_project/evidence/demo_interlaces_envelope_weak_n_5_6.json")
        .exists());

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 15,
            "method": "tools/call",
            "params": {
                "name": "check_family_interlacing",
                "arguments": {
                    "left_family_id": "chebyshev_u_polynomial",
                    "right_family_id": "chebyshev_t_polynomial",
                    "left_offset": -1,
                    "n_min": 2,
                    "n_max": 3,
                    "mode": "strict"
                }
            }
        }),
    );
    let offset_interlacing = read_response(&mut stdout, 15);
    assert_eq!(
        offset_interlacing["result"]["structuredContent"]["report"]["items"][0]["left_n"],
        1
    );
    assert_eq!(
        offset_interlacing["result"]["structuredContent"]["report"]["all_interlacing"],
        true
    );

    drop(stdin);
    let _ = child.kill();
    let _ = child.wait();
    let _ = fs::remove_dir_all(root);
}
