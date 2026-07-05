# sym-poly-qsym roadmap

This crate contains checked quasisymmetric-function code used by
`symmetricfunctions.com` examples.

## Implemented and Usable

- `QSymFunction<C>` with monomial, fundamental, quasisymmetric Schur,
  dual immaculate, extended Schur, row-strict extended Schur, flipped
  extended Schur, backward extended Schur, and two quasisymmetric power-sum
  bases.
- `p_partition.rs`: Stanley `(P,w)`-partition generating functions and
  labeled linear-extension data.
- `chromatic_qsym.rs`: chromatic quasisymmetric functions and asc-weighted
  variants.
- `peak.rs`: Stembridge peak quasisymmetric functions in the fundamental
  basis, plus quasisymmetric Schur Q-functions via standard peak composition
  tableaux and peak Young quasisymmetric Schur functions via standard peak
  Young composition tableaux.
- `power_sum.rs`: `Psi` and `Phi` quasisymmetric power-sum bases, plus the
  `p_alpha` and reverse `p^r_alpha` P-partition combinatorial power-sum
  bases of Aliniaeifard--Wang--van Willigenburg.
- `schur_qsym.rs`: composition-tableau helpers for quasisymmetric Schur
  functions, standard immaculate tableaux for dual immaculate functions,
  Young quasisymmetric Schur functions via the rho/reverse-index relation,
  row-strict quasisymmetric Schur functions via the omega relation,
  row-strict Young quasisymmetric Schur functions via the omega relation,
  and row-strict dual immaculate functions in monomial and fundamental bases.
- `kohnert_qsym.rs`: stable Kohnert limits for finite diagrams, stable key
  limits from left-justified key diagrams, extended Schur functions from
  right-justified lock diagrams, and row-strict, flipped, and backward
  extended Schur functions via the source-paper involution relations.
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
- `examples/peak_quasisymmetric_degree4.rs`: degree-four peak
  quasisymmetric fundamental table for `tex-source/peakQuasisymmetric.tex`.
- `examples/peak_tableau_site_examples.rs`: checked `(4,3,1)` SPCT/SPYCT
  `ytableau` example for `tex-source/qsymSchur.tex`.
- `examples/qsym_schur_q_degree4.rs`: degree-four quasisymmetric Schur
  Q-function peak-basis table for `tex-source/qsymSchur.tex`.
- `examples/peak_young_qsym_schur_degree4.rs`: degree-four peak Young
  quasisymmetric Schur peak-basis table for `tex-source/qsymSchur.tex`.
- `examples/dual_immaculate_site_example.rs`: checked example for
  `tex-source/qsymSchur.tex`.
- `examples/dual_immaculate_degree4.rs`: degree-four dual immaculate `qmonom`
  table for `tex-source/qsymSchur.tex`.
- `examples/row_strict_dual_immaculate_degree4.rs`: degree-four
  row-strict dual immaculate `qmonom` table for `tex-source/qsymSchur.tex`.
- `examples/extended_schur_degree4.rs`: degree-four extended Schur `qmonom`
  and fundamental tables for future `tex-source/qsymSchur.tex` examples.
- `examples/combinatorial_power_sum_degree4.rs`: degree-four `p_alpha` and
  reverse `p^r_alpha` monomial tables for future
  `tex-source/pPartitions.tex` or `tex-source/qsymSchur.tex` examples.
- `examples/kohnert_key_stable_degree4.rs`: degree-four stable key Kohnert
  monomial and fundamental tables for future `tex-source/kohnert.tex`
  examples.
- Dual immaculate monomial tests compare packed semistandard tableaux with the
  standard-tableau fundamental expansion through degree 4.
- Row-strict dual immaculate tests include the paper example
  `RS*_(1,2) = M_(2,1) + M_(1,1,1)` and the identity
  `RS*_alpha = psi(S*_alpha)`.
- Extended Schur tests check Assaf--Searles paper examples and properties:
  the published `E_(2,1,2)` fundamental expansion, the adjacent-swap
  `F`-positive difference example, reverse-hook single-fundamental examples,
  stable Kohnert invariance under added empty rows, monomial and fundamental
  positivity through degree 4, agreement with Schur functions for partition
  shapes through degree 4, Daugherty's involution/partition specializations
  for the four extended families, and round-trip basis conversion through
  both fundamental and monomial bases in degree 4.
- P-partition combinatorial power-sum tests check the Aliniaeifard--Wang--van
  Willigenburg degree-four monomial table, the reverse-basis monomial formula,
  the `p_(2,1,1)` partition-refinement example, the product example
  `p_(1,2)p_(1) = (1/2)p_(1,2,1)+p_(1,1,2)`, the involution examples for
  `p_(1,1,2)`, and round-trip basis conversion for `p_alpha` and
  `p^r_alpha` in degree 4.
- Stable key Kohnert tests check the left-justified key-diagram convention,
  agreement with Schur functions for partition shapes through degree 4, and
  the Demazure-character stable-limit behavior on a non-dominant weak
  composition with an internal zero row.

## Website Families Still Missing or Needing Verification

The site discusses several QSym-related families that should get explicit Rust
support before new computed examples are added.

- Quasisymmetric Macdonald `E` and `J` polynomials.  These need a
  source-checked HHL/nonattacking-filling implementation, preferably sharing
  conventions with the existing nonsymmetric Macdonald code in
  `sym-poly/multipoly`.
- Quasisymmetric hook Schur functions.  These need the Mason--Niese super
  Gessel/fundamental model and should probably live in a superspace-aware
  module rather than overloading ordinary `QSymFunction`.  The ordinary
  hook-tableau generator is now available in `sym-poly/sym/src/hook_schur.rs`
  and can be reused as the supersymmetric half of this implementation.
- General Kohnert quasisymmetric examples beyond key and lock diagrams.  The
  stable Kohnert engine is present; more named diagram constructors and
  checked examples can now be added incrementally.

For each new family, prefer:

- one generator with a named tableau/filling type when a tableau model exists;
- a small example under `examples/` whose output can be copied into TeX;
- tests matching examples from the source paper, and Sage comparisons when a
  Sage implementation is available.
