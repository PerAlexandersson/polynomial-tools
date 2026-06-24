use serde_json::{json, Value};
use std::io::{BufRead, BufReader, Write};
use std::process::{Child, ChildStdin, Command, Stdio};

fn start_server() -> (Child, ChildStdin, BufReader<std::process::ChildStdout>) {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool-mcp"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("spawn polytool-mcp");
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
fn lists_tools_and_calls_properties() {
    let (mut child, mut stdin, mut stdout) = start_server();

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 1,
            "method": "initialize",
            "params": {
                "protocolVersion": "2025-03-26",
                "capabilities": {},
                "clientInfo": { "name": "polynomial-tools-mcp-test", "version": "0.1.0" }
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
    let tools_array = tools["result"]["tools"].as_array().expect("tool list");
    let tool_names: Vec<&str> = tools_array
        .iter()
        .filter_map(|tool| tool["name"].as_str())
        .collect();
    assert!(tool_names.contains(&"polynomial_properties"));
    assert!(tool_names.contains(&"check_polynomial_family"));
    assert!(tool_names.contains(&"find_recurrence"));
    let recurrence_schema = &tools_array
        .iter()
        .find(|tool| tool["name"] == "find_recurrence")
        .expect("find_recurrence tool")["inputSchema"];
    assert!(recurrence_schema["oneOf"].is_array());
    assert!(recurrence_schema["properties"]["coefficients"].is_object());

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 3,
            "method": "tools/call",
            "params": {
                "name": "polynomial_properties",
                "arguments": {
                    "text": "1, 11, 11, 1"
                }
            }
        }),
    );
    let result = read_response(&mut stdout, 3);
    let structured = &result["result"]["structuredContent"];
    assert_eq!(structured["items"][0]["real_rooted"], true);
    assert_eq!(structured["items"][0]["gamma_coefficients"], json!([1, 8]));

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 4,
            "method": "tools/call",
            "params": {
                "name": "check_polynomial_family",
                "arguments": {
                    "sequence": "eulerian",
                    "max_n": 4,
                    "options": {
                        "require_real_rooted": true,
                        "require_gamma_positive": true
                    }
                }
            }
        }),
    );
    let family = read_response(&mut stdout, 4);
    let structured = &family["result"]["structuredContent"];
    assert_eq!(structured["all_required_checks_passed"], true);
    assert!(structured["markdown"]
        .as_str()
        .expect("markdown")
        .contains("Polynomial family check"));

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 5,
            "method": "tools/call",
            "params": {
                "name": "find_recurrence",
                "arguments": {
                    "coefficients": [[1], [1], [2], [3], [5], [8]]
                }
            }
        }),
    );
    let recurrence = read_response(&mut stdout, 5);
    let structured = &recurrence["result"]["structuredContent"];
    assert_eq!(structured["found"], true);
    assert_eq!(structured["recurrence"], "P(n) = P(n-1) + P(n-2)");

    drop(stdin);
    let _ = child.kill();
    let _ = child.wait();
}
