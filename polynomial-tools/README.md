# polynomial-tools

Dense univariate polynomial toolkit for combinatorial research. Provides
real-rootedness testing by primitive integer root counting, Bézout-matrix
interlacing checks, log-concavity, gamma-positivity, resultants,
Ehrhart/h\*-vector conversion, recurrence search for polynomial sequences,
finite Athanasiadis--Wagner interlacing matrices, and standard sequences —
all with exact arithmetic.

## Installation

From the current workspace root:

```sh
cargo build --release -p polynomial-tools
```

The CLI binary is at `target/release/polytool`.

In a future standalone `polynomial-tools` repository, the equivalent command is:

```sh
cargo build --release
```

## Using as a crate

The library crate is usable today from another Rust project by depending on
this directory:

```toml
[dependencies]
polynomial-tools = { path = "/home/paxinum/Dropbox/AI-projects/rust/polynomial-tools" }
```

From a public Git repository containing this Rust workspace, Cargo can also
depend on the package by name:

```toml
[dependencies]
polynomial-tools = { git = "https://github.com/USER/REPO.git", package = "polynomial-tools" }
```

If this directory is pushed as its own standalone GitHub repository, the
dependency can instead be:

```toml
[dependencies]
polynomial-tools = { git = "https://github.com/USER/polynomial-tools.git" }
```

After publishing to crates.io, the dependency would be:

```toml
[dependencies]
polynomial-tools = "0.1"
```

The public API is re-exported from `polynomial_tools`, so a downstream crate
can write:

```rust
use polynomial_tools::{check_weak_interlacing, is_real_rooted};

assert!(is_real_rooted(&[1, 11, 11, 1]));
assert_eq!(check_weak_interlacing(&[2, -3, 1], &[-1, 1]), Some(true));
```

Coefficient vectors are always in ascending degree order. Most convenience
functions accept `&[i64]`; use the `*_bigint_coeffs` variants when coefficients
may exceed `i64`. Interlacing functions return `Option<bool>`: `None` means the
directed degree relation is not valid for that test, not that interlacing failed.

For Hermite--Biehler or Euclidean-chain experiments, the library also exposes a
signed Sturm continued-fraction certificate checker over exact rationals:

```rust
use num_bigint::BigInt;
use num_rational::BigRational;
use polynomial_tools::check_sturm_continued_fraction_certificate;

let q = |n: i64| BigRational::from_integer(BigInt::from(n));
let p0 = vec![q(1), q(3), q(1)]; // x^2 + 3x + 1
let p1 = vec![q(2), q(1)];       // x + 2
let certificate =
    check_sturm_continued_fraction_certificate(&p0, &p1).expect("Sturm-CF certificate");
assert_eq!(certificate.steps, 2);
```

To build the MCP server from the workspace root:

```sh
cargo build --release -p polynomial-tools-mcp
```

The MCP binary is at `target/release/polytool-mcp`.

To install the MCP server into a user bin directory and print a client
configuration snippet:

```sh
./polynomial-tools/mcp/install.sh
```

In a future standalone `polynomial-tools` repository, run:

```sh
./mcp/install.sh
```

## CLI usage

### Input format

Polynomials are given as comma-separated integer coefficients in
**ascending degree order**: `a_0, a_1, ..., a_d` represents
`a_0 + a_1 t + ... + a_d t^d`.

A text file `polys.txt` might look like:

```
1, 11, 11, 1
1, 26, 66, 26, 1
1, 57, 302, 302, 57, 1
```

### Check real-rootedness

```sh
polytool real-rooted < polys.txt
```

Output (one line per polynomial):

```
1 + 11t + 11t^2 + t^3: real-rooted
1 + 26t + 66t^2 + 26t^3 + t^4: real-rooted
1 + 57t + 302t^2 + 302t^3 + 57t^4 + t^5: real-rooted
```

### Check interlacing

Test whether consecutive polynomials interlace:

```sh
polytool interlacing < polys.txt
```

### Check unimodality, log-concavity, palindromicity, gamma-positivity

```sh
polytool properties < polys.txt
```

Output includes all properties for each polynomial:

```
1 + 11t + 11t^2 + t^3: real-rooted, palindromic, gamma-positive [1, 8], unimodal, log-concave
```

For a palindromic polynomial, print the full gamma expansion:

```sh
polytool gamma-expansion < polys.txt
```

Example output:

```text
1 + 11t + 11t^2 + t^3: gamma [1, 8]; expansion: (1+t)^3 + 8 t (1+t)
```

### Search for a recurrence

Given a sequence of polynomials (one per line), search for a polynomial
recurrence `f(n,t) P_n(t) = Σ c_{r,d}(n,t) D^d P_{n-r}(t)`. With
`--alternating-sign`, the search also allows terms
`(-1)^n c_{r,d}(n,t) D^d P_{n-r}(t)`:

```sh
polytool recurrence < polys.txt
```

Options:

```
--skip-prefix <k>    Ignore the first k input polynomials before searching
--full-depth         Require all offsets 1..depth to appear in the recurrence
--min-rec-len <k>    Minimum recurrence depth to try (default: 1)
--max-rec-len <k>    Maximum recurrence depth (default: 5)
--min-var-deg <d>    Minimum degree in t for coefficients (default: 0)
--max-var-deg <d>    Maximum degree in t for coefficients (default: 3)
--min-idx-deg <d>    Minimum degree in n for coefficients (default: 0)
--max-idx-deg <d>    Maximum degree in n for coefficients (default: 3)
--min-diff-deg <d>   Minimum derivative order (default: 0)
--max-diff-deg <d>   Maximum derivative order (default: 2)
--inhomogeneous      Also try inhomogeneous recurrences
--min-inhomo-var-deg Minimum degree in t for the inhomogeneous term
--max-inhomo-var-deg Maximum degree in t for the inhomogeneous term
--min-inhomo-idx-deg Minimum degree in n for the inhomogeneous term
--max-inhomo-idx-deg Maximum degree in n for the inhomogeneous term
--denominator        Allow a nontrivial LHS factor f(n,t)
--alternating-sign   Also allow right-hand-side terms multiplied by (-1)^n
--max-denom-var-deg  Max degree in t for f(n,t) (default: 2, implies --denominator)
--max-denom-idx-deg  Max degree in n for f(n,t) (default: 2, implies --denominator)
--verbose            Print each candidate tried
```

### Scout BKW equal-modulus loci

For a polynomial recurrence with characteristic symbol

```text
F(x,z) = a_0(x) + a_1(x) z + ... + a_r(x) z^r,
```

`bkw-scout` numerically scans a complex rectangle for points where the two
dominant characteristic roots have nearly equal modulus.  This is a scout for
Beraha--Kahane--Weiss accumulation obstructions; it does not certify dominance,
amplitudes, or eventual non-real-rootedness.

Give the symbol as the z-coefficient polynomials in ascending z-degree:

```sh
polytool bkw-scout \
  --symbol '1; -x; 1' \
  --box -3 3 -2 2 \
  --grid 61 \
  --top 10
```

The same symbol can be supplied on stdin, one coefficient polynomial per line.
Use `--format json` for machine-readable output, and `--mathematica` to print a
`Reduce` skeleton for exact equal-modulus follow-up.

### Compute resultant/discriminant

```sh
# Resultant of two polynomials (given as two lines)
echo "1, 0, 1
1, -1" | polytool resultant

# Discriminant of a single polynomial
echo "1, 0, -3, 1" | polytool discriminant
```

### Ehrhart ↔ h\*-vector conversion

```sh
# h*-vector → Ehrhart polynomial
echo "1, 8, 35, 32, 9" | polytool hstar-to-ehrhart

# Ehrhart polynomial → h*-vector (coefficients as rationals: num/den)
echo "1/1, 2/1, 1/1" | polytool ehrhart-to-hstar
```

## Library usage

Add to your `Cargo.toml`:

```toml
[dependencies]
polynomial-tools = { path = "../polynomial-tools" }
```

### Real-rootedness and interlacing

```rust
use polynomial_tools::*;

// Coefficients in ascending degree order: coeffs[i] = coeff of t^i
let eulerian_4 = [1, 11, 11, 1];

// Default exact PRS/root-counting method
assert!(is_real_rooted(&eulerian_4));

// Sturm chain method (slower, gives root locations)
assert!(is_real_rooted_sturm(&eulerian_4));

// Strict interlacing: deg(p) = deg(q) + 1
let p = [-15, 23, -9, 1];  // (t-1)(t-3)(t-5)
let q = [8, -6, 1];         // (t-2)(t-4)
assert_eq!(check_interlacing(&p, &q), Some(true));

// Weak interlacing (allows shared roots)
let f = [2, -3, 1];   // (t-1)(t-2)
let g = [-1, 1];       // (t-1)
assert_eq!(check_weak_interlacing(&f, &g), Some(true));

// Same-degree interlacing: roots alternate on the real line
// (t-1)(t-3) roots {1,3} and (t-2)(t-4) roots {2,4}: 1 < 2 < 3 < 4
assert_eq!(check_weak_interlacing(&[3, -4, 1], &[8, -6, 1]), Some(true));

// Nested roots do NOT interlace: (t-1)(t-4) vs (t-2)(t-3)
assert_eq!(check_weak_interlacing(&[4, -5, 1], &[6, -5, 1]), Some(false));
```

For large exact coefficients, use the `BigInt` APIs:

```rust
use num_bigint::BigInt;
use polynomial_tools::check_interlacing_bigint_coeffs;

let center = BigInt::from(10).pow(20);
let f = vec![-&center, BigInt::from(1)]; // t - center
let g = vec![
    (&center - 1u32) * (&center + 1u32),
    -BigInt::from(2) * &center,
    BigInt::from(1),
];

assert_eq!(check_interlacing_bigint_coeffs(&f, &g), Some(true));
```

### Athanasiadis--Wagner interlacing matrices

```rust
use polynomial_tools::*;

// A column vector of polynomials P, Q, R in ascending coefficient order.
let polys = vec![
    vec![2, 1],       // 2 + t
    vec![8, 6, 1],    // (2+t)(4+t)
    vec![3, 4, 1],    // (1+t)(3+t)
];

// One finite truncation of the infinite Lace(A) matrix.
let lace = lace_matrix_sequence_i64(&polys, 1, 3).unwrap();
assert_eq!(lace, vec![vec![2, 1, 0], vec![8, 6, 1], vec![3, 4, 1]]);

// This is Athanasiadis--Wagner's pairwise-but-not-fully-interlacing example:
// the 3 x 3 determinant is negative.
assert!(check_lace_sequence_total_nonnegative_i64(&polys, 1, 3, 3).is_err());
```

For a `p x q` polynomial matrix `A`, `lace_matrix_i64(A, r, c)` returns the
`p r` by `q c` truncation whose `(p rb+i, q cb+j)` entry is the coefficient
of `x^(cb-rb)` in `A_ij(x)`.  This follows Definition 3.5 of
Athanasiadis--Wagner, *Veronese sections and interlacing matrices of formal
power series*.  A finite TNN check is computational evidence for the infinite
matrix, not a proof of full interlacing unless a separate finite criterion
applies.

### Polynomial arithmetic

```rust
use polynomial_tools::Polynomial;

let p = Polynomial::<i64>::from_i64_coeffs(&[1, 1]);  // 1 + t
let q = p.clone() * p.clone();                         // 1 + 2t + t^2
assert_eq!(q.evaluate(3), 16);
assert!(q.is_palindromic());

let deriv = q.derivative();  // 2 + 2t
let gcd = Polynomial::<i64>::gcd(&q, &deriv);
```

### Recurrence search

```rust
use polynomial_tools::recurrence::*;

// Eulerian polynomials A_1, A_2, ..., A_10
let polys: Vec<Vec<i64>> = vec![
    vec![1],
    vec![1, 1],
    vec![1, 4, 1],
    vec![1, 11, 11, 1],
    vec![1, 26, 66, 26, 1],
    vec![1, 57, 302, 302, 57, 1],
    vec![1, 120, 1191, 2416, 1191, 120, 1],
    vec![1, 247, 4293, 15619, 15619, 4293, 247, 1],
    // ...
];

// Adaptive search: tries small parameter spaces first
let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default());
if let Some(res) = result {
    println!("{}", res.recurrence);
}

// Or search with specific parameters
let opts = RecurrenceOptions {
    rec_len: 2,
    var_deg: 1,
    idx_deg: 1,
    diff_deg: 1,
    ..Default::default()
}.with_alternating_sign(true);
if let Some(rec) = find_polynomial_recurrence(&polys, &opts) {
    println!("{}", rec);
}
```

### Other properties

```rust
use polynomial_tools::*;

let coeffs = [1, 11, 11, 1];

assert!(is_palindromic(&coeffs));
assert!(is_unimodal(&coeffs));
assert!(is_log_concave(&coeffs));
assert!(is_ultra_log_concave(&coeffs));
assert!(is_gamma_positive(&coeffs));

// Gamma coefficients: p(t) = Σ γ_i t^i (1+t)^{d-2i}
let gamma = gamma_coefficients(&coeffs).unwrap();
assert_eq!(gamma, vec![1, 8]);

// Resultant and discriminant
let disc = discriminant(&[1, 0, -3, 1]);
```

## MCP server

The repository also contains a local Model Context Protocol server in
`mcp/`. It exposes the exact polynomial routines to MCP clients over stdio.

Build it from the workspace root:

```sh
cargo build --release -p polynomial-tools-mcp
```

Install it locally:

```sh
./polynomial-tools/mcp/install.sh
```

Check the installed binary:

```sh
polytool-mcp --help
polytool-mcp --version
```

Example MCP client configuration:

```json
{
  "mcpServers": {
    "polynomial-tools": {
      "command": "/absolute/path/to/rust/target/release/polytool-mcp"
    }
  }
}
```

The MCP server is tools-only: no resources, prompts, HTTP server, sampling, or
filesystem access. It returns structured JSON and also includes the same compact
JSON as text content for client compatibility.

Available MCP tools:

- `parse_polynomials`
- `check_polynomial_family`
- `polynomial_properties`
- `check_interlacing_pair`
- `check_interlacing_sequence`
- `real_roots`
- `find_recurrence`
- `resultant`
- `discriminant`
- `ehrhart_hstar`
- `analyze_decomposition`
- `generate_sequence`

See [`mcp/README.md`](mcp/README.md) for request schemas, examples, and
development notes.

## Development

Run the core tests:

```sh
cargo test -p polynomial-tools
```

Run the MCP server tests:

```sh
cargo test -p polynomial-tools-mcp
```

Run the web wrapper tests:

```sh
cargo test -p polynomial-tools-web
```

If this directory is split out into its own Git repository, keep the core crate,
`mcp/`, and `web/` under one workspace manifest, and keep `Cargo.lock` if
reproducible binary builds are important.

## Support

Report issues in the Git repository, or contact Per Alexandersson
(@PerAlexandersson, <per.w.alexandersson@gmail.com>).

## Algorithm notes

### Real-rootedness algorithms

The default real-rootedness check uses the primitive integer PRS path from
`root_count`.  It replaces the polynomial by its square-free part and counts
distinct real roots exactly over `BigInt`; the polynomial is real-rooted
precisely when this count equals the square-free degree.  For one-signed
coefficient polynomials, the backend first tries cheap coefficient filters and
then counts positive roots of `f(-t)`, since such polynomials can only have
non-positive real roots.

Explicit matrix comparison paths remain available as
`is_real_rooted_bezout`, `is_real_rooted_bezout_squarefree`, and
`is_real_rooted_hermite`.

All interlacing checks use the **Bézout matrix** (Fisk, Cor. 9.145). For
polynomials f (degree d) and g (degree d−1), the Bézout matrix B(f,g) is the
d×d symmetric matrix with (i,j) entry equal to the coefficient of x^i y^j in
`(f(x)g(y) - f(y)g(x)) / (x-y)`.

**Theorem:** f and g are both real-rooted and g strictly interlaces f
if and only if B(f,g) is positive definite.

This reduces interlacing to a single exact matrix definiteness check, avoiding
root isolation entirely. It is 100–400× faster than Sturm chains at degree 15+
for interlacing checks.  The definiteness implementation uses fraction-free
BigInt Bareiss elimination below dimension `30`, and switches to a
CRT-over-prime-fields path for larger matrices.  The modular path reconstructs
the leading principal minors exactly using Hadamard bounds.  For speed
comparisons, the explicit entry points are `is_positive_definite_bareiss`,
`is_positive_definite_modular`, and `modular_leading_principal_minors_bigint`.

Run the local comparison example with:

```bash
cargo run --release -p polynomial-tools --example bench_bareiss_vs_modular_linalg
```

For **same-degree** polynomials, `check_weak_interlacing` reduces to the
deg+1 case by extending one polynomial with a root far to the right
(using the Cauchy bound for root radius). When all coefficients are
positive (all roots negative), multiplying by t suffices, skipping the
Cauchy bound computation.

See `doc/bezout-interlacing.md` for the Mathematica reference implementation.

### Interlacing matrices

Athanasiadis--Wagner's interlacing matrix construction packages a polynomial
sequence or polynomial matrix into an infinite block Toeplitz matrix `Lace(A)`.
For a column vector, total nonnegativity of this infinite matrix is called
full interlacing and implies pairwise interlacing of the entries.  The converse
fails, so this gives a stronger finite-experiment target than checking
ordinary pairwise interlacing.

The crate implements finite truncations through `lace_matrix_i64`,
`lace_matrix_bigint`, and the sequence wrappers
`lace_matrix_sequence_i64` / `lace_matrix_sequence_bigint`.  Exact finite TNN
checks are available with `check_lace_total_nonnegative_i64` and the faster
Neville-elimination wrappers.
