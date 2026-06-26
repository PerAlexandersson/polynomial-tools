# polynomial-tools MCP server

`polynomial-tools-mcp` is a local Model Context Protocol server for
`polynomial-tools`. It exposes exact univariate polynomial computations to MCP
clients over stdio.

The server is intentionally small:

- stdio transport only
- tools only; no resources, prompts, sampling, HTTP server, or filesystem access
- structured JSON results for agent-friendly use
- exact arithmetic inherited from `polynomial-tools`

## Build

From the Rust workspace root:

```sh
cargo build --release -p polynomial-tools-mcp
```

The binary is:

```text
target/release/polytool-mcp
```

Human-facing command-line help is available without starting the MCP server:

```sh
polytool-mcp --help
polytool-mcp --version
```

Running `polytool-mcp` with no arguments starts the stdio MCP server. In that
mode, stdout is reserved for MCP protocol messages.

## Install

For a local user install, run the installer from the repository root:

```sh
./polynomial-tools/mcp/install.sh
```

If `polynomial-tools` later lives in its own Git repository, the same script can
be run as:

```sh
./mcp/install.sh
```

By default the script:

- builds `polytool-mcp` in release mode
- uses `target` under the `polynomial-tools` directory as the Cargo target
  directory
- copies the binary to `$HOME/.local/bin/polytool-mcp`
- prints an MCP client configuration snippet with the absolute binary path

Useful options:

```sh
# Install somewhere else
./mcp/install.sh --prefix "$HOME/.local"

# Put the binary directly in a specific directory
./mcp/install.sh --bin-dir "$HOME/bin"

# Use a specific Cargo binary
./mcp/install.sh --cargo /home/user/.rustup/toolchains/stable-x86_64-unknown-linux-gnu/bin/cargo

# Reuse an already-built binary from BUILD_DIR/release/polytool-mcp
./mcp/install.sh --build-dir /tmp/polytool-mcp-target --skip-build
```

This machine currently needs the rustup toolchain binaries rather than the snap
wrappers, so a reliable local install command is:

```sh
env \
  RUSTC=/home/peal0658/.rustup/toolchains/stable-x86_64-unknown-linux-gnu/bin/rustc \
  RUSTDOC=/home/peal0658/.rustup/toolchains/stable-x86_64-unknown-linux-gnu/bin/rustdoc \
  ./polynomial-tools/mcp/install.sh \
  --cargo /home/peal0658/.rustup/toolchains/stable-x86_64-unknown-linux-gnu/bin/cargo
```

On another computer, install Rust with `rustup`, clone the repository, and run
the same installer. Cargo will fetch the Rust dependencies from crates.io.

If your shell resolves `cargo`, `rustc`, or `rustdoc` through a wrapper that does
not work in the current environment, use the rustup toolchain binaries directly.
On this computer, that build command is:

```sh
env \
  RUSTC=/home/peal0658/.rustup/toolchains/stable-x86_64-unknown-linux-gnu/bin/rustc \
  RUSTDOC=/home/peal0658/.rustup/toolchains/stable-x86_64-unknown-linux-gnu/bin/rustdoc \
  /home/peal0658/.rustup/toolchains/stable-x86_64-unknown-linux-gnu/bin/cargo \
  build --release -p polynomial-tools-mcp
```

## MCP client configuration

Use an absolute path to the compiled binary. A typical local configuration is:

```json
{
  "mcpServers": {
    "polynomial-tools": {
      "command": "/absolute/path/to/rust/target/release/polytool-mcp"
    }
  }
}
```

The server writes MCP protocol messages to stdout. Do not add debugging prints to
stdout; diagnostics should go to stderr.

For troubleshooting:

```sh
polytool-mcp --help
polytool-mcp --version
```

## Polynomial input format

All coefficient vectors are in ascending degree order:

```text
[a0, a1, a2] = a0 + a1 t + a2 t^2
```

Most tools accept one of these forms:

```json
{ "coefficients": [1, 11, 11, 1] }
```

```json
{ "expression": "1 + 11t + 11t^2 + t^3" }
```

`PolynomialInput` requires exactly one of `coefficients` or `expression`.

Batch tools accept exactly one of:

```json
{
  "polynomials": [
    { "coefficients": [1, 11, 11, 1] },
    { "expression": "1 + 26t + 66t^2 + 26t^3 + t^4" }
  ]
}
```

```json
{
  "text": "1, 11, 11, 1\n1 + 26t + 66t^2 + 26t^3 + t^4"
}
```

Text batches use the same parser as the library:

- comma-separated coefficients: `1, 2, 3`
- whitespace-separated coefficients: `1 2 3`
- bracketed coefficient lists: `[1, 2, 3]`
- expanded polynomial notation: `1 + 2t + 3t^2`

Blank lines and lines beginning with `#` are ignored in text batches.

## Output format

Every tool returns structured MCP output. `rmcp` also includes the same compact
JSON in a text content block for compatibility with clients that do not read
`structuredContent`.

Integer values that may exceed machine size are serialized as strings. Rational
values are serialized as strings such as `"3/2"` or `"1"`.

For batch tools, invalid polynomial lines become per-item errors. Invalid
top-level request shapes, such as providing both `text` and `polynomials`, are
reported as MCP invalid-params errors.

## Tools

### `parse_polynomials`

Parse and normalize one or more polynomials.

Example arguments:

```json
{
  "text": "1, 11, 11, 1\n1 + 2t + t^2"
}
```

Returns `items`, each with `index`, `ok`, and either normalized polynomial data
or an `error`.

### `check_polynomial_family`

Run the common research checks in one exact batch.  The tool accepts either
explicit polynomials, text input, or one of the built-in standard sequences.
It reports polynomial properties, consecutive interlacing data, optional finite
`Lace(A)` total-nonnegativity checks, and optional recurrence search.  It also
returns a Markdown summary suitable for a project note or `HANDOFF.md`.

Example with a generated sequence:

```json
{
  "sequence": "eulerian",
  "max_n": 6,
  "options": {
    "require_real_rooted": true,
    "require_gamma_positive": true,
    "check_consecutive_interlacing": true
  }
}
```

Example with the Athanasiadis--Wagner pairwise-but-not-fully-interlacing
triple:

```json
{
  "polynomials": [
    { "coefficients": [2, 1] },
    { "coefficients": [8, 6, 1] },
    { "coefficients": [3, 4, 1] }
  ],
  "lace": {
    "block_rows": 1,
    "block_cols": 3,
    "max_minor_size": 3,
    "include_matrix": true
  }
}
```

The `all_required_checks_passed` field uses conservative defaults:
real-rootedness is required, while gamma-positivity, log-concavity,
palindromicity, and consecutive interlacing are reported but only required
when their `require_*` flags are set.  A requested `lace` check is required to
be TNN, and a requested recurrence search is required to find a recurrence.

### `polynomial_properties`

Compute standard properties for each input polynomial:

- real-rooted
- simple roots
- palindromic, ignoring initial zeros
- gamma-positive, ignoring initial zeros
- gamma coefficients
- log-concave
- ultra-log-concave

Example arguments:

```json
{
  "polynomials": [
    { "coefficients": [1, 11, 11, 1] }
  ]
}
```

### `check_interlacing_pair`

Check strict and weak interlacing for two polynomials.

Example arguments:

```json
{
  "p": { "coefficients": [-2, 1] },
  "q": { "coefficients": [3, -4, 1] }
}
```

Returns normalized `p` and `q`, plus `strict`, `weak`, and a status string.
`strict` and `weak` can be `null` when degrees are incompatible or a required
real-rootedness condition fails.

### `check_interlacing_sequence`

Apply `check_interlacing_pair` to consecutive pairs in a batch.

Example arguments:

```json
{
  "text": "1, 4, 1\n1, 11, 11, 1\n1, 26, 66, 26, 1"
}
```

Returns parsed `items` and pairwise `pairs`.

### `check_interlacing_profile`

For each row, check previous rows backward and stop at the first failure.  This
is useful for quickly profiling how long the local interlacing chain persists.

Example arguments:

```json
{
  "text": "1\n1, 1\n1, 2, 1\n1, 0, 1"
}
```

Returns parsed `items` and a `profile` array.  Each profile row includes
`previous_count`, `checked_previous_count`, `interlacing_previous_count`,
strict/weak counts, the successful previous indices, and the checked pair
reports.  The checked pair reports are ordered backward from the immediately
previous row and include the first failure when there is one.

### `real_roots`

Return rational midpoint representatives from the existing Sturm isolation
helper, or `real_rooted: false`.

Example arguments:

```json
{
  "polynomials": [
    { "expression": "t^2 - 1" },
    { "expression": "t^2 + 1" }
  ]
}
```

### `find_recurrence`

Search adaptively for a recurrence in a polynomial sequence.

Provide exactly one input source:

- `coefficients`: dense coefficient vectors, one polynomial per row;
- `expressions`: polynomial expressions such as `"1 + 4*x + x^2"`;
- `polynomials`: explicit `PolynomialInput` objects;
- `text`: newline-separated polynomial data.

Compact coefficient example:

```json
{
  "coefficients": [[1], [1], [2], [3], [5], [8]]
}
```

Explicit polynomial example:

```json
{
  "polynomials": [
    { "coefficients": [1] },
    { "coefficients": [1] },
    { "coefficients": [2] },
    { "coefficients": [3] },
    { "coefficients": [5] },
    { "coefficients": [8] }
  ]
}
```

Important option fields:

- `skip_prefix`
- `min_rec_len`, `max_rec_len`
- `min_var_deg`, `max_var_deg`
- `min_idx_deg`, `max_idx_deg`
- `min_diff_deg`, `max_diff_deg`
- `try_inhomogeneous`
- `min_inhomo_var_deg`, `max_inhomo_var_deg`
- `min_inhomo_idx_deg`, `max_inhomo_idx_deg`
- `try_denominator`
- `try_alternating_sign`
- `max_denom_var_deg`, `max_denom_idx_deg`
- `min_margin`
- `fit_extra_rows`
- `no_verify`

`verbose` is intentionally not exposed through MCP, so the server never writes
search traces into the stdio protocol stream.

When a recurrence is found, the response includes:

- `recurrence`: compact plaintext recurrence
- `latex`
- `mathematica`
- `sage`
- `python`: exact standalone Python code using `fractions.Fraction`
- `recurrence_json`: a saved recurrence record with initial values, suitable
  for `generate_recurrence_rows` or `polytool recurrence-generate`
- search metadata such as `unknowns`, `equations`, `fit_polynomials`,
  `verification_polynomials`, and `candidates_tried`

### `generate_recurrence_rows`

Generate or extend coefficient rows from `recurrence_json` returned by
`find_recurrence`.

Example arguments:

```json
{
  "recurrence_json": "{ ... JSON string returned by find_recurrence ... }",
  "rows": 20
}
```

Use exactly one of:

- `rows`: total number of rows to return, including initial rows;
- `additional`: append this many rows after the initial rows.

The response contains `first_index`, `row_count`, and `polynomials`, where
coefficients are exact rational strings in ascending powers of `t`.

### `resultant`

Compute the exact resultant of two polynomials.

Example arguments:

```json
{
  "p": { "coefficients": [1, 0, 1] },
  "q": { "coefficients": [-1, 1] }
}
```

The resultant is returned as a string.

### `discriminant`

Compute exact discriminants for one or more polynomials.

Example arguments:

```json
{
  "text": "1, 0, -3, 1\n-1, 0, 1"
}
```

Discriminants are returned as strings.

### `ehrhart_hstar`

Convert between h*-vectors and Ehrhart polynomial coefficient vectors.

h*-vector to Ehrhart polynomial:

```json
{
  "mode": "hstar_to_ehrhart",
  "hstar": [1, 1, 0]
}
```

Ehrhart polynomial to h*-vector using rational strings:

```json
{
  "mode": "ehrhart_to_hstar",
  "ehrhart_coefficients": ["1", "2", "1"]
}
```

Ehrhart polynomial to h*-vector using a common denominator:

```json
{
  "mode": "ehrhart_to_hstar",
  "numerator_coefficients": [2, 3, 1],
  "denominator": 2
}
```

### `analyze_decomposition`

Analyze symmetric decomposition, R-transform data, interlacing checks, and
magic-basis coordinates.

Example arguments:

```json
{
  "polynomials": [
    { "coefficients": [1, 3, 1] }
  ]
}
```

### `generate_sequence`

Generate standard polynomial sequences.

Supported `sequence` values:

- `eulerian`
- `narayana`
- `type_b_eulerian`
- `chebyshev_t`
- `chebyshev_u`
- `hermite`

Example arguments:

```json
{
  "sequence": "eulerian",
  "max_n": 6
}
```

## Development

Run the MCP server tests:

```sh
cargo test -p polynomial-tools-mcp
```

Run the polynomial-tools baseline checks:

```sh
cargo test -p polynomial-tools
cargo test -p polynomial-tools-web
```

Format the MCP package:

```sh
cargo fmt -p polynomial-tools-mcp
```

## Repository notes

This package currently lives as a workspace member under the broader Rust
workspace. If `polynomial-tools` is split into its own Git repository, keep
these pieces together:

- the core crate in `polynomial-tools/`
- the MCP package in `polynomial-tools/mcp/`
- the web package in `polynomial-tools/web/`, if the browser UI should remain
- a workspace manifest that includes `mcp` and `web`
- `Cargo.lock`, if release binaries or reproducible local MCP installs matter
- CI jobs for `cargo test -p polynomial-tools`, `cargo test -p polynomial-tools-mcp`,
  and `cargo test -p polynomial-tools-web`

Before publishing, also fill in repository metadata in the package manifests.

## Support

Report issues in the Git repository, or contact Per Alexandersson
(@PerAlexandersson, <per.w.alexandersson@gmail.com>).
