# polynomial-lab-mcp

MCP server exposing the `polynomial-lab` project index through stdio tools.

Run without arguments from an MCP client:

```bash
POLY_LAB_ROOT=/workspace/projects/polynomial-interlacing-lab poly-lab-mcp
```

The server writes protocol messages only to stdout. Diagnostics must go to
stderr.
