# sym-poly-qsym roadmap

This crate contains checked quasisymmetric-function code used by
`symmetricfunctions.com` examples.

## Implemented and Usable

- `QSymFunction<C>` with monomial, fundamental, quasisymmetric Schur,
  dual immaculate, and two quasisymmetric power-sum bases.
- `p_partition.rs`: Stanley `(P,w)`-partition generating functions and
  labeled linear-extension data.
- `chromatic_qsym.rs`: chromatic quasisymmetric functions and asc-weighted
  variants.
- `power_sum.rs`: `Psi` and `Phi` quasisymmetric power-sum bases.
- `schur_qsym.rs`: composition-tableau helpers for quasisymmetric Schur
  functions, standard immaculate tableaux for dual immaculate functions,
  Young quasisymmetric Schur functions via the rho/reverse-index relation,
  row-strict quasisymmetric Schur functions via the omega relation,
  row-strict Young quasisymmetric Schur functions via the omega relation,
  and row-strict dual immaculate functions in monomial and fundamental bases.
- Quasisymmetric Schur tests include the Tewari--van Willigenburg
  `qSchur_(2,1,3)` fundamental expansion example from
  `tex-source/qsymSchur.tex`.
- Young quasisymmetric Schur tests include the
  `yqSchur_(3,1,2) = rho(qSchur_(2,1,3))` source-relation example.
- `examples/qsym_schur_degree4.rs`: degree-four quasisymmetric Schur
  fundamental table for `tex-source/qsymSchur.tex`.
- `examples/young_qsym_schur_degree4.rs`: degree-four Young
  quasisymmetric Schur fundamental table for `tex-source/qsymSchur.tex`.
- `examples/row_strict_qsym_schur_degree4.rs`: degree-four row-strict
  quasisymmetric Schur fundamental table for `tex-source/qsymSchur.tex`.
- `examples/row_strict_young_qsym_schur_degree4.rs`: degree-four row-strict
  Young quasisymmetric Schur fundamental table for
  `tex-source/qsymSchur.tex`.
- `examples/p_partition_site_example.rs`: checked example for
  `tex-source/pPartitions.tex`.
- `examples/dual_immaculate_site_example.rs`: checked example for
  `tex-source/qsymSchur.tex`.
- `examples/dual_immaculate_degree4.rs`: degree-four dual immaculate `qmonom`
  table for `tex-source/qsymSchur.tex`.
- `examples/row_strict_dual_immaculate_degree4.rs`: degree-four
  row-strict dual immaculate `qmonom` table for `tex-source/qsymSchur.tex`.
- Dual immaculate monomial tests compare packed semistandard tableaux with the
  standard-tableau fundamental expansion through degree 4.
- Row-strict dual immaculate tests include the paper example
  `RS*_(1,2) = M_(2,1) + M_(1,1,1)` and the identity
  `RS*_alpha = psi(S*_alpha)`.

## Website Families Still Missing or Needing Verification

The site page `tex-source/qsymSchur.tex` also discusses several Schur-like
families that should get explicit Rust support before new computed examples are
added.

- Quasisymmetric Schur Q-functions (`qSchurQ`).
- Peak Young quasisymmetric Schur functions (`peakYqSchur`).
- Extended Schur functions (`extSchur`) and their lock-polynomial stable-limit
  interpretation.

For each new family, prefer:

- one generator with a named tableau/filling type when a tableau model exists;
- a small example under `examples/` whose output can be copied into TeX;
- tests matching examples from the source paper, and Sage comparisons when a
  Sage implementation is available.
