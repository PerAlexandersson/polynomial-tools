# sym-poly

`sym-poly` is the multi-crate symmetric-function workspace used for checked
examples on `symmetricfunctions.com`.

## Tableau Breadcrumbs

- `core/src/tableau.rs` contains the shared `Tableau` and `SkewTableau` types,
  standard and semistandard generators, reading words, descents, RSK, promotion,
  evacuation, crystals, charge, and key-tableau helpers.
- `core/src/ssaf.rs` contains semi-standard augmented fillings, including atom
  fillings, key fillings, permuted basements, and Mason's SSYT-to-SSAF map.
- `qsym/src/schur_qsym.rs` contains composition-tableau and immaculate-tableau
  generators for quasisymmetric Schur and dual immaculate functions.
- `qsym/README.md` tracks which Mason/Assaf/Searles-style QSym Schur variants
  are implemented, missing, or need definition-level verification.
- `multipoly/src/key_polynomial.rs` computes key polynomials and tests them
  against the SSAF weight enumerator.

When adding website examples, prefer a small checked Rust example under the
relevant crate's `examples/` directory and cite it from TeX with a short
`% Related Rust:` comment near the example.
