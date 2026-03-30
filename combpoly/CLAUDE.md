# Project: combpoly — Combinatorial generating polynomials

## Purpose

Library and CLI for generating polynomials from combinatorial objects
(permutations, words, parking functions) by counting with a statistic.

Core operation: given a set of objects and a statistic,
`build_generating_polynomial` returns the coefficient vector where
`coeffs[k] = #{w : stat(w) = k}`.

## CLI tool: `combpoly`

```bash
cargo run --release -- poly --perms 7 --avoiding 312 --stat des --real-rooted
cargo run --release -- scan --size 7 --avoiding 312 --ideal bruhat --stat exc
cargo run --release -- list --perms 5 --avoiding 312
cargo run --release -- recurrence --perms 3:12 --avoiding 312 --stat des --auto
```

## Library modules (`src/`)

| Module | Contents |
|--------|----------|
| `permutation.rs` | Generation, pattern avoidance, backtrack_image, filtered_permutations |
| `statistics.rs` | 18 statistics: des, exc, peak, valley, inv, maj, lrmax, fix, cyc, ... |
| `polynomial_builder.rs` | `build_generating_polynomial(objects, stat)` |
| `order.rs` | Bruhat/weak order ideals |
| `word.rs` | Words on Ferrers boards |
| `parking.rs` | Parking functions, run-sorted variants |
| `catalan.rs` | Catalan/Dyck path utilities |
| `cayley.rs` | Cayley permutations |

## Dependencies

- `polynomial-tools` for real-rootedness, interlacing, format_poly, recurrence search
- Polynomial analysis is NOT in this crate — use `polynomial_tools` directly

## Exploration binaries

Research binaries are in the `experiments/` crate (workspace sibling), not here.

## Paper: Backtrack permutations

Associated research paper in `paper/Backtrack-permutations.tex`.
See the paper directory for current status.

## OEIS connections

| Triangle | OEIS | Context |
|----------|------|---------|
| Narayana (des on Av_132 etc.) | A001263 | h-vector of associahedron |
| Av_321 + des | A091156 | Dyck paths by long ascents |
| Av_132 + des | A048994 | Unsigned Stirling 1st kind |
| Peak on Av_132 (Class A) | A091894 | Dyck paths by ddu's |
| Peak on Av_312 (Class B) | A236406 | 123-avoiding by peaks |
| Exc on Av_312 | NOT IN OEIS | New, fails real-rootedness |
| Exc on Av_231 | NOT IN OEIS | New, real-rooted, no recurrence |

## Style notes

- Paper uses amsart, explicit `\ref` (no cleveref).
- LaTeX macros: `\backtrack`, `\symS`, `\Av`, `\des`, `\exc`, `\peak`, `\ltr`, `\oeis{Annnnnn}`.
