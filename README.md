# Rust combinatorics workspace

This Cargo workspace collects combinatorics and polynomial code used for
experiments, checked examples, and `symmetricfunctions.com` support.

## Library Map

- `combinatoric-core`: foundational partitions, compositions, permutations,
  graphs, and posets.
- `sym-poly/core`: shared symmetric-function data structures, including
  tableaux, skew tableaux, semi-standard augmented fillings, crystals, and
  exact arithmetic helpers.
- `sym-poly/sym`, `sym-poly/qsym`, and `sym-poly/multipoly`: symmetric,
  quasisymmetric, and multivariate polynomial code, including Schur,
  quasisymmetric Schur, key, atom, slide/glide, lock, Kohnert/Lascoux,
  Schubert, and nonsymmetric Macdonald filling-formula routines.
- `polytool`: exact univariate polytool for real-rootedness,
  interlacing, gamma-positivity, recurrence searches, and Ehrhart/h*-vectors.
- `combpoly`: command-line exploration of combinatorial polynomials from
  permutations, words, and parking functions.
- `kostka`: Kostka coefficients, Gelfand--Tsetlin Ehrhart polynomials, and
  h*-vectors.

For website examples that need tableau, SSAF, or key-polynomial verification,
start in `sym-poly/README.md` and add a small checked example under the
relevant crate's `examples/` directory.

## Polytool repository layout

The `master` branch is the canonical history for this entire workspace.
Polytool development belongs under `polytool/` on `master`; do not commit
directly to `main`.

The `main` branch is a generated standalone view of `polytool/`, maintained
with `git subtree split`. After committing and pushing polytool work on
`master`, publish the standalone branch with:

```sh
./scripts/sync-polytool-main.sh
```

The script requires a clean `polytool/` tree and refuses to overwrite a
diverged `main`. This keeps the standalone crate available without maintaining
a second implementation or history by hand.
