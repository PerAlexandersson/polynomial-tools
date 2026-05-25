# polynomial-tools refactoring notes

These notes record cleanup opportunities after the May 25, 2026 pass.  The
crate is currently test-clean and clippy-clean for `polynomial-tools`.

## Current structure

- `real_rootedness.rs` owns the public real-rootedness, interlacing,
  gamma-positivity, resultant, discriminant, and Ehrhart helpers.  It is the
  main file that could use future splitting.
- `root_count.rs` owns the primitive integer PRS/Sturm implementation.  This is
  the preferred exact path for one-signed combinatorial polynomials.
- `linalg.rs` owns exact determinant, definiteness, and TNN checks.  The
  index-heavy loops are intentional matrix elimination code; the module has a
  local clippy allowance for `needless_range_loop`.
- `basis.rs` and `decomposition.rs` are small, coherent modules and should stay
  separate from `real_rootedness.rs`.

## Near-term refactors

- Split `real_rootedness.rs` into smaller internal modules:
  `formatting`, `bezout`, `interlacing`, `qpoly`, `ehrhart`, and
  `resultants`.  Re-export the same public API from `real_rootedness`.
- Move the private rational-polynomial helpers in `real_rootedness.rs`
  (`poly_gcd_q`, `poly_rem_q`, `poly_exact_div_q`, `q_poly_to_i64`) into a
  small internal module.  They are not conceptually tied to interlacing.
- Consider using `root_count::is_real_rooted_fast_bigint_coeffs` as the general
  mixed-sign fallback for default real-rootedness once it has benchmark coverage
  against the Bézout path on the examples that matter for current projects.
- Add focused benchmarks for positive-coefficient families of degree 20-80:
  Eulerian, Narayana, h*-vectors, and known counterexamples.  The existing
  examples are useful, but they do not yet define a regression benchmark suite.

## Things not to change casually

- Do not replace the exact PRS/Bézout checks with floating-point root finding.
- Keep coefficient vectors in ascending degree order throughout the public API.
- Keep the zero polynomial conventions documented in tests; several research
  callers rely on them.
