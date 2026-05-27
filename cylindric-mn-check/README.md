# Cylindric MN Inclusion-Exclusion Checker

Small brute-force checker for the appendix proof in
`Cylindric-Schur-functions.tex`.

The checker models a cylindric diagram as boxes on a one-period
cylinder, forms the associated poset, and computes

```tex
\sum_{S\subseteq E}(-1)^{|S|}m_S(P)
```

by exhaustive enumeration of subsets of strict edges.

For a fixed `S`, the forced equality classes are computed as strongly
connected components after adding reverse edges for the strict edges in
`S`.  The value `m_S(P)` is the size of the unique source component in
the condensation DAG, or zero if there is not a unique source.

Run:

```bash
cargo run -- --max-ordinary-moves 10 --max-loop-len 9 --bad-grid 4
```

Current checked output:

```text
ordinary ribbons checked: 2047, failures: 0
loop ribbons checked: 1004, failures: 0
anchored residual ordinary ribbons checked: 2047, failures: 0
anchored residual loop ribbons checked: 1004, failures: 0
anchored residual non-ribbon grid subsets checked: 372, failures: 0
ordinary non-ribbon grid subsets checked: 61780, failures: 0
```

The anchored residual checks model the point in the appendix after the
loop layers have been collapsed.  They add one weighted source node for
the collapsed stack and strict connector edges from that node to the
minimal boxes of the residual shape `F`.  The checker verifies that the
residual contribution is independent of the source weight and has the
sign predicted by the open/closed edge count.

Use `--residual-report` to print the individual anchored residual
cases instead of only the summary.

The `--picture-report` flag prints diagnostics for a direct encoding
of the local shapes in `toggleProofPics.tex`.  This is only a notation
sanity check for the proof pictures, not a full stacked-ribbon test.

The `--stacked-report` flag prints diagnostic data for naive shifted
loop bands.  These are not currently claimed to be faithful models of
the paper's stacked-ribbon peeling.

The `--extended-report` flag runs an exploratory scan of longer
path-generated cylindric shapes.  This is deliberately not part of the
pass/fail result, because the generator produces many path-like shapes
that are connected and have no `2 x 2` block but are not the stacked
ribbons defined in the draft.
