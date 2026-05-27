# flagged-lorentzian

Reusable Rust tools for the flagged skew Schur normalized Lorentzian project.

The crate is deliberately separate from `experiments/`: code here should be
stable enough to reuse across sessions, with typed APIs and regression tests for
counterexamples and promising refinements.

Current focus:

- enumerate row-flagged skew semistandard tableaux;
- group tableaux by content and descent-type statistics;
- check the \(2\times2\) descent-fiber inequalities
  \[
  \mathcal T_\beta(i,i)\times \mathcal T_\beta(j,j)
  \hookrightarrow
  \mathcal T_\beta(i,j)\times \mathcal T_\beta(i,j)
  \]
  fiberwise with respect to unordered pair statistics.
- work with skew GT patterns and pair-sum GT arrays for the active-row /
  pair-envelope retile model.

The first durable statistic implemented here is componentwise descent data:
connected components of a skew shape are read separately, avoiding artificial
descents created by concatenating disconnected components in one global row
word.

## CLI

Scan a single ordinary skew shape:

```text
cargo run -q -p flagged-lorentzian --bin descent_fiber_scan -- \
  --lambda 4,3,1 --mu 3,1 --alphabet 5 --stat global
```

The same shape with componentwise descents:

```text
cargo run -q -p flagged-lorentzian --bin descent_fiber_scan -- \
  --lambda 4,3,1 --mu 3,1 --alphabet 5 --stat componentwise
```

Scan an ordinary family:

```text
cargo run -q -p flagged-lorentzian --bin descent_fiber_scan -- \
  --max-skew-size 7 --max-outer-extra 8 --alphabet 6 --stat componentwise
```

Recent stress-test commands:

```text
cargo run -q -p flagged-lorentzian --bin descent_fiber_scan -- \
  --max-skew-size 8 --max-outer-extra 5 --alphabet 6 \
  --row-flags 3,4,5,6 --stat componentwise --tableau-limit 500000
```

```text
cargo run -q -p flagged-lorentzian --bin descent_fiber_scan -- \
  --max-skew-size 8 --max-outer-extra 5 --alphabet 7 \
  --row-flags 4,5,6,7 --stat componentwise --tableau-limit 500000
```

Both passed in the current workspace run.  These are evidence-gathering scans,
not proofs; keep the exact command and caps with any reported result.

## GT Retile Scans

The current leading \(2\times2\) model uses the binary
`gt_exchange_matching_scan`.  It compares negative pairs with mixed pairs after
compressing each pair to

```text
(active GT row sum, pair envelope, optional descent data, pair-sum GT array).
```

The option `--matching-mode flow` checks existence by max-flow.  The option
`--matching-mode greedy` uses a deterministic least-reachable target rule:
negative keys and positive keys are sorted, and each negative key is assigned
to the first still-available mixed key with the same active row and an allowed envelope
reachable within the allowed elementary exchange depth.

Use `--envelope-mode exact`, `--envelope-mode nonincrease`, or
`--active-row-only` to control the flag-envelope condition.  The current best
flag-compatible model is `nonincrease`, i.e. the target may drop to a lower
sharp-envelope layer but may not move upward.

Use `--invariant-level R` to preserve a different pair-sum GT row, or
`--upper-invariant` to preserve the row above the active adjacent labels.  This
is useful for lower flagged labels, where preserving the lower active row can
fail.

Use `--descent-mode componentwise` or `--descent-mode global` to include the
unordered pair of descent data in the exchange key.  The `active-*` variants
are available as negative tests; active-subword descents are known to fail
already for the straight shape \((2,1)\).

Representative command:

```text
cargo run --release -q -p flagged-lorentzian --bin gt_exchange_matching_scan -- \
  --max-skew-size 6 --max-outer-extra 7 --alphabet 5 --lower-label 4 \
  --max-exchange-depth 2 --matching-mode greedy
```

In the current run this passed through skew size \(6\): \(4650\) shapes,
\(114885\) fibers, and \(4513111\) negative pairs.  The same greedy depth-\(2\)
rule also passed a size-\(7\) ordinary family with outer bound
`skew_size + 5`, and greedy depth \(1\) passed the flagged shape
\((7,6,5)/(3,2)\) with flags \((4,5,5)\).

The same ordinary size-\(6\) and size-\(7\) exact-envelope scans pass with
`--descent-mode componentwise`; connected ordinary shapes pass through size
\(7\) with `--descent-mode global`.  Exact envelope first fails in the tested
range at \((4,4,3,2)/(2,2,1)\), \(\beta=(2,3,1,0,0)\).  With
`--envelope-mode nonincrease`, greedy depth \(1\) passes connected ordinary
shapes through skew size \(8\), outer bound `skew_size + 7`, with global
descents.  It also passes the ordinary skew-size-\(8\) family with outer bound
`skew_size + 5` and componentwise descents.

The diagnostic binary `gt_envelope_drop_dump` prints the first exact-envelope
failure and its one-step lower-envelope repair.

The diagnostic binary `gt_lower_label_failure_dump` prints the first flagged
lower-label failure of the \(Z_i\)-preserving rule.  It is a reminder that the
current strong evidence is for the top-label/local-top mechanism; lower flagged
labels probably need an interval version of the GT retile.
