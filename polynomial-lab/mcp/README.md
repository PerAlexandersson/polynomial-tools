# polynomial-lab-mcp

MCP server exposing the `polynomial-lab` project index through stdio tools.

Run without arguments from an MCP client:

```bash
POLY_LAB_ROOT=/workspace/projects/polynomial-interlacing-lab poly-lab-mcp
```

The server writes protocol messages only to stdout. Diagnostics must go to
stderr.

Example MCP registration:

```json
{
  "mcpServers": {
    "polynomial-lab": {
      "command": "poly-lab-mcp",
      "args": [
        "--root",
        "/workspace/projects/polynomial-interlacing-lab"
      ]
    }
  }
}
```

Exposed tools:

- `validate_lab`
- `list_projects`
- `get_project`
- `list_goals`
- `list_evaluations`
- `trace_goal_support`
- `create_project`
- `append_goal`
- `append_conjecture`
- `append_implication`
- `list_proof_rules`
- `list_search_recipes`
- `list_families`
- `compute_family`
- `check_family_real_rooted`
- `check_family_interlacing`
- `render_project_markdown`
- `render_project_html`
- `write_project_markdown`
- `write_project_html`
- `append_evaluation`
- `append_counterexample`
- `append_timeout`

The record-writing tools create new TOML or JSON files with create-new
semantics.  They never overwrite existing records.  Evidence append tools write
under `projects/<project-id>/evidence/<record-id>.json`.

`check_family_interlacing` supports optional `left_offset` and `right_offset`
integer fields, so one can check relations such as `F_{n-1} \ll G_n`.
