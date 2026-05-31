# polynomial-lab

`polynomial-lab` is a structured index over project data for polynomial
real-rootedness and interlacing experiments.

The default lab root is:

```text
/workspace/projects/polynomial-interlacing-lab
```

Set `POLY_LAB_ROOT` or pass `--root` to point at another lab directory.

Examples:

```bash
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- validate
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- validate --strict
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- list-projects
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- list-families
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- \
  compute-family derangement_descent_polynomial --n 8
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- \
  check-family-real-rooted derangement_descent_polynomial \
  --n-min 2 --n-max 30 \
  --project derangement_descents \
  --relation derangement_descent_real_rootedness_checked \
  --append
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- \
  trace-goal derangement_descents derangement_descent_real_rootedness
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- \
  render-markdown derangement_descents
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- \
  render-markdown derangement_descents --write
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- \
  render-html derangement_descents --write
```

Append-only evidence records:

```bash
timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- \
  append-evaluation derangement_descents my_relation_n_1_20 \
  --relation my_relation \
  --status holds_for_checked_domain \
  --method exact_bezout_check \
  --n-min 1 --n-max 20

timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- \
  append-counterexample derangement_descents my_false_strengthening_n_4 \
  --relation my_false_strengthening \
  --method exact_bezout_check \
  --n 4 \
  --failure-reason "The proposed strengthening fails at n=4."

timeout 60s nice -n 10 cargo run -p polynomial-lab --bin poly-lab -- \
  append-timeout derangement_descents my_relation_timeout \
  --relation my_relation \
  --seconds 60 \
  --method exact_bezout_check \
  --n-min 30 --n-max 40
```

Tests use `tests/fixtures/minimal_lab` as a portable complete example with a
goal, definitions, a conjecture, an implication, a success record, a failure
record, a proof rule, and a search recipe.  There is also a small integration
test that checks the real derangement lab when `/workspace/projects/...` is
available.

## Family Registry

The family registry deliberately lives in `polynomial-lab`, not
`polynomial-tools`.  Stable standard families are backed by `polynomial-tools`;
project-specific families remain in the lab namespace until they become broadly
useful.

Currently registered families:

- `eulerian_polynomial`
- `narayana_polynomial`
- `type_b_eulerian_polynomial`
- `chebyshev_t_polynomial`
- `chebyshev_u_polynomial`
- `hermite_polynomial`
- `derangement_descent_polynomial`
- `normalized_derangement_descent_polynomial`
- `reciprocal_eulerian_derivative_polynomial`
