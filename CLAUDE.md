# Workspace: Rust combinatorics projects

## Layout

```
rust/
  Cargo.toml                  ← workspace root (shared Cargo.lock)
  sym-poly/core/              ← Ring, Partition, Composition, BasisIndex, FormalSum
  sym-poly/sym/               ← Symmetric functions (6 classical bases)
  sym-poly/qsym/              ← Quasisymmetric functions (M, F, Ψ, Φ bases)
  sym-poly/multipoly/         ← Multivariate polynomials, divided differences, key polys
  polynomial-tools/           ← Polynomial toolkit (will be public)
  combinatoric-core/          ← Graphs, posets, chromatic; re-exports sym-poly
  combpoly/                   ← Permutation statistics + combpoly CLI
  kostka/                     ← Public crate: Kostka/LR coefficients, GT polytopes
  experiments/                ← 44 research exploration binaries
```

### sym-poly/ (multi-crate)
Modular workspace for symmetric polynomial algebras. Four crates, all generic
over `Ring` coefficients (i64, BigInt, Ratio<BigInt>, Ratio<i64>).

**sym-poly/core/** — Shared foundation:
- `ring` — Ring trait and implementations
- `partition`, `composition` — Indexing types with full combinatorial utilities
- `index` — `BasisIndex` trait (generic over Partition/Composition)
- `formal_sum` — `FormalSum<I, C>` generic linear combination
- `matrix` — Integer matrix utilities (multiply, transpose, invert)
- `transition_cache` — `TransitionCache<B>` per-algebra cached transition matrices

**sym-poly/sym/** — Symmetric functions:
- `SymmetricFunction<C>` with 6 bases (m, e, h, p, s, f), all 30 conversions
- Kostka coefficients, S_n characters, omega involution, multiplication
- Hall inner product, skew Schur (Jacobi-Trudi), plethysm, specializations

**sym-poly/qsym/** — Quasisymmetric functions:
- `QSymFunction<C>` with 4 bases: M (monomial), F (fundamental), Ψ (type 1), Φ (type 2)
- Power sum bases Ψ, Φ per Ballantine--Daugherty--Hicks--Mason--Niese (JCTA 2020)
- Normalized Ψ̃ = Ψ/z, Φ̃ = Φ/z utilities
- Omega involution (ω² = id, ω(Ψ_α) = (-1)^{n-ℓ(α)} Ψ_{α^r})
- Quasi-shuffle multiplication, Sym↔QSym maps
- (P,w)-partitions (Ψ̃-positive for naturally labeled posets)
- Chromatic quasisymmetric functions

**sym-poly/multipoly/** — Multivariate polynomials:
- `MultiPoly<C>` sparse multivariate polynomial type
- Divided difference operators ∂_i and π_i (Demazure)
- Key polynomials κ_α via operator recursion

### polynomial-tools/
Univariate polynomial toolkit for combinatorial research. Uses Bézout matrices
(not Sturm chains) as the default for real-rootedness and interlacing (100-400x faster).

**Modules:**
- `polynomial` — `Polynomial<C>` with `CoeffRing`/`FieldRing` traits, arithmetic,
  derivative, evaluate, shift, reverse, dilate, GCD, division, Lagrange interpolation
- `real_rootedness` — Bézout/Sturm real-rootedness, strict/weak interlacing,
  log-concavity, ultra-log-concavity, palindromic, gamma-positivity,
  resultant, discriminant, Ehrhart h*-vector conversion, format_poly
- `sturm` — Sturm chains (internal, used as fallback)
- `recurrence` — Adaptive recurrence search for polynomial sequences
- `sequences` — Eulerian, Narayana, type B Eulerian, Chebyshev T/U, Hermite

### combinatoric-core/
Facade crate: re-exports sym-poly-core and sym-poly-sym for backward compatibility.
Retains local modules for combinatorial structures not in the function algebra hierarchy.

**Re-exported** (from sym-poly): `Ring`, `Partition`, `Composition`, `Basis`,
`SymmetricFunction`, `kostka`, `transition` — all existing `use combinatoric_core::*`
paths continue to work.

**Local modules:**
- `graph` — Graph type with generators (complete, path, cycle, bipartite,
  multipartite, Ferrers board, unit interval, Petersen), operations (complement,
  induced subgraph, delete/contract edge), predicates (connected, bipartite,
  claw-free), matchings, independence number, graph6 parser
- `poset` — Poset type (Hasse diagram). Constructors: chain, antichain, fence,
  k-alternating, from Young/skew diagram, dual. Algorithms: linear extensions,
  order-preserving maps (backtracking + frontier DP), order polytope Ehrhart
  polynomial and h*-vector (BigRational, no overflow), P-Eulerian polynomial,
  natural relabeling. Frontier DP gives 100-18000x speedup over backtracking.
- `chromatic` — Chromatic symmetric function (Sym version, via deletion-contraction)
- `key_polynomial` — Key polynomials κ_{λ,σ} via Kogan faces of GT-polytopes,
  Ehrhart polynomials for key-Kostka coefficients

**Numeric policy / upgrade path:**
- Prefer one exact implementation plus one ergonomic wrapper, not two separate algorithms.
- For integer-valued combinatorial functions whose coefficients may grow, the target pattern is:
  `foo_bigint(...) -> Vec<BigInt>` as the canonical implementation, and
  `foo(...) -> Vec<i64>` as a checked convenience wrapper calling `to_i64().expect(...)`.
- For genuinely rational algorithms (for example Ehrhart interpolation), keep the exact
  `BigRational` version as the primary API.
- Avoid making every public API fully generic over coefficient type unless the generic
  version clearly improves reuse; generic helper layers are preferred over genericizing
  every algorithm.

### combpoly/
Permutation statistics and generating polynomials. Provides the `combpoly` CLI.

**Library modules:**
- `permutation` — Generation, pattern avoidance, backtrack image, filtered generation
- `statistics` — 18 statistics: des, exc, peak, valley, inv, maj, lrmax, fix, cyc, ...
- `polynomial_builder` — `build_generating_polynomial(objects, stat)`
- `order` — Bruhat/weak order ideals
- `word` — Words on Ferrers boards
- `parking` — Parking functions
- `catalan` — Catalan/Dyck path utilities
- `cayley` — Cayley permutations

### kostka/
Public crate. Kostka coefficients (DP), LR coefficients (GT polytope),
Ehrhart polynomials and h*-vectors of Gelfand-Tsetlin polytopes. Has CLI.

### experiments/
44 research exploration binaries (peak_*, cayley_*, backtrack_*, pf_*, stembridge_*, etc.).
Depends on combpoly and polynomial-tools. Flat layout with prefix-based grouping.

## Naming conventions

- Use verbose, Mathematica-style names: `conjugate_partition`, `kostka_coefficient`,
  `to_monomial_basis`, `schur_symmetric`, `build_generating_polynomial`.
- CLI should be AI-friendly: use `clap` with detailed `--help`, support `--format json` output.

## Mathematica reference packages

Located in `~/AI-projects/combinatoric-tools/mathematica-packages/` (shared, read-only):
- `SymmetricFunctions.m` (~2000 lines) — full symmetric function package
- `CombinatoricTools.m` (~1600 lines) — partitions, compositions, tableaux, charge
- `NewTableaux.m` — SSYT/SYT generation, cylindric tableaux

Also: `bezout-interlacing.md` in workspace root has Mathematica code for the Bézout
interlacing algorithm, ready to add to the symmetric functions package.

## What's next

### sym-poly/sym — Tier 2-4
- [ ] SSYT generation
- [ ] Ribbon Schur functions
- [ ] Cylindric Schur functions
- [ ] Hall-Littlewood P, T (needs `Polynomial<BigInt>` coefficients)
- [ ] Jack P, J (needs rational parameter)
- [ ] Schur Q, P (shifted symmetric functions)
- [ ] Macdonald P, H, J
- [ ] LLT polynomials, Delta/Nabla operators

### sym-poly/qsym
- [ ] Young quasisymmetric Schur (YQS) basis
- [ ] General (P,w)-partitions (non-natural labelings)
- [ ] q-chromatic quasisymmetric function (Shareshian-Wachs, polynomial coefficients in q)
- [ ] NCSym and Hopf algebra duality with QSym

### sym-poly/multipoly
- [ ] Macdonald E-polynomials (needs Ring impl for Polynomial<Ratio<BigInt>>)
- [ ] Permuted basement Macdonalds
- [ ] Schubert polynomials via divided differences

## References
- Stanley, EC2 for symmetric function foundations
- Fisk, "Polynomials, roots, and interlacing" for Bézout matrix theory
- Postnikov (2005) for cylindric Schur functions
- Ballantine--Daugherty--Hicks--Mason--Niese, JCTA 2020, for QSym power sums Ψ and Φ
