# Recurrence Benchmark Fixtures

Generated exact coefficient rows for adaptive recurrence-search benchmarks.
Each `rows/*.txt` file contains one dense coefficient list per line, with
coefficients in ascending powers of `t`. These files contain no headers, so
they can be piped directly into `polytool recurrence`.
Metadata lives separately in `manifest.tsv` and in the table below.
The matching `json/*.json` files are recurrence JSON records emitted by
`polytool recurrence --json`; they include minimal initial conditions and
can regenerate or extend the raw row files with `recurrence-generate`.
Fixtures `01`--`23` are synthetic stress tests. Fixtures `24` onward are
natural OEIS-derived cases copied from the curated
`projects/real-rooted-oeis/sequences` queue.

Regenerate from the Rust workspace root with:

```sh
cargo run -p polynomial-tools --example generate_recurrence_benchmarks
bash polynomial-tools/fixtures/recurrence-benchmarks/regenerate-json.sh
```

Run the built-in timing suite:

```sh
polytool bench recurrence-fixtures --repeat 3
polytool bench recurrence-fixtures --only 23_sparse --repeat 5
```

Example one-off timing command:

```sh
time cargo run -q -p polynomial-tools --bin polytool -- recurrence \
  --max-rec-len 1 --max-var-deg 1 --max-idx-deg 0 --max-diff-deg 0 \
  < polynomial-tools/fixtures/recurrence-benchmarks/rows/03_binomial_powers.txt
```

Example regeneration command:

```sh
polytool recurrence-generate \
  --recurrence polynomial-tools/fixtures/recurrence-benchmarks/json/03_binomial_powers.json \
  --rows 50
```

| slug | features | suggested args | recurrence |
|---|---|---|---|
| `01_scalar_geometric` | constant coefficient, first order | `--max-rec-len 1 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0` | `P_n = 2 P_{n-1}` |
| `02_scalar_fibonacci` | constant coefficient, second order | `--max-rec-len 2 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0` | `P_n = P_{n-1} + P_{n-2}` |
| `03_binomial_powers` | t-dependent coefficient | `--max-rec-len 1 --max-var-deg 1 --max-idx-deg 0 --max-diff-deg 0` | `P_n = (1+t) P_{n-1}` |
| `04_chebyshev_t` | t-dependent coefficient, second order | `--max-rec-len 2 --max-var-deg 1 --max-idx-deg 0 --max-diff-deg 0` | `P_n = 2t P_{n-1} - P_{n-2}` |
| `05_factorial_index` | n-dependent coefficient | `--max-rec-len 1 --max-var-deg 0 --max-idx-deg 1 --max-diff-deg 0` | `P_n = n P_{n-1}` |
| `06_affine_product` | n- and t-dependent coefficient | `--max-rec-len 1 --max-var-deg 1 --max-idx-deg 1 --max-diff-deg 0` | `P_n = (n+t) P_{n-1}` |
| `07_hermite_indexed_second_order` | n-dependent coefficient, second order | `--max-rec-len 2 --max-var-deg 1 --max-idx-deg 1 --max-diff-deg 0` | `P_n = 2t P_{n-1} - 2(n-2) P_{n-2}` |
| `08_inhomogeneous_linear` | inhomogeneous, degree one in n and t | `--inhomogeneous --max-rec-len 1 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0 --max-inhomo-var-deg 1 --max-inhomo-idx-deg 1` | `P_n = P_{n-1} + n + t` |
| `09_inhomogeneous_quadratic` | inhomogeneous, degree two in n and t | `--inhomogeneous --max-rec-len 1 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0 --max-inhomo-var-deg 2 --max-inhomo-idx-deg 2` | `P_n = P_{n-1} + n^2 + nt + t^2` |
| `10_alternating_scalar` | alternating sign | `--alternating-sign --max-rec-len 1 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0` | `P_n = (-1)^n P_{n-1}` |
| `11_alternating_fibonacci` | alternating sign, second order | `--alternating-sign --max-rec-len 2 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0` | `P_n = P_{n-1} + (-1)^n P_{n-2}` |
| `12_eulerian_derivative` | first derivative, n- and t-dependent coefficients | `--max-rec-len 1 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 1` | `P_n = (1+(n-2)t)P_{n-1} + (t-t^2)P'_{n-1}` |
| `13_derivative_appell` | first derivative, n- and t-dependent coefficient | `--max-rec-len 1 --max-var-deg 1 --max-idx-deg 1 --max-diff-deg 1` | `P_n = (n+t)P_{n-1} + P'_{n-1}` |
| `14_second_derivative` | second derivative | `--max-rec-len 1 --max-var-deg 2 --max-idx-deg 0 --max-diff-deg 2` | `P_n = (1+t)P_{n-1} + t^2 P''_{n-1}` |
| `15_mixed_derivative_second_order` | first derivative, n- and t-dependent coefficient, second order | `--max-rec-len 2 --max-var-deg 1 --max-idx-deg 1 --max-diff-deg 1` | `P_n = (1+nt)P_{n-1} + tP'_{n-1} + P_{n-2}` |
| `16_denominator_linear_index` | LHS denominator in n | `--denominator --max-rec-len 1 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0 --max-denom-idx-deg 1` | `(n+1)P_n = P_{n-1}` |
| `17_denominator_quadratic_index` | LHS denominator, quadratic in n | `--denominator --max-rec-len 1 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0 --max-denom-idx-deg 2` | `(1+n+n^2)P_n = P_{n-1}` |
| `18_denominator_with_t_rhs` | LHS denominator, t-dependent RHS | `--denominator --max-rec-len 1 --max-var-deg 1 --max-idx-deg 0 --max-diff-deg 0 --max-denom-idx-deg 1` | `(n+1)P_n = (1+t)P_{n-1}` |
| `19_denominator_derivative` | LHS denominator, first derivative | `--denominator --max-rec-len 1 --max-var-deg 1 --max-idx-deg 0 --max-diff-deg 1 --max-denom-idx-deg 1` | `(n+1)P_n = P_{n-1} + tP'_{n-1}` |
| `20_complex_mixed_alternating_derivative` | alternating sign, first derivative, n- and t-dependent coefficients, second order | `--alternating-sign --max-rec-len 2 --max-var-deg 1 --max-idx-deg 1 --max-diff-deg 1` | `P_n = (n+t)P_{n-1} + (1-t)P'_{n-1} + (-1)^n tP_{n-2}` |
| `21_fifth_order_mixed_coefficients` | fifth-order recurrence, n^2- and t^5-dependent coefficients | `--max-rec-len 5 --max-var-deg 5 --max-idx-deg 2 --max-diff-deg 0` | `P_n = tP_{n-1} + (t^2-t-n)P_{n-2} + (t^2+tn+1)P_{n-3} + (t^4+t^3n^2+1)P_{n-4} - t^5P_{n-5}` |
| `22_fifth_derivative_generic` | second-order recurrence, derivatives through order five, n- and t-dependent coefficients | `--max-rec-len 2 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 5` | `P_n = (2+t-nt^2)P_{n-1} + (-1+3t)P'_{n-1} + (1+t^2)P''_{n-2} + (4-2n+t)P'''_{n-2} + (-3+t+nt^2)P^{(4)}_{n-2} + (5+t-2nt)P^{(5)}_{n-2}` |
| `23_sparse_second_derivative_lag` | second derivative, third-order recurrence, sparse t^2 coefficients | `--max-rec-len 3 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 2 --fit-extra-rows 2` | `P_n = tP_{n-1} + t^2P''_{n-1} - t^2P_{n-3} + nt^2P''_{n-3}` |
| `24_oeis_a114655_denominator_derivative` | natural OEIS recurrence, denominator, first derivative | `--min-rec-len 1 --max-rec-len 1 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1` | `(n+1) P(n) = (3nt+2n-3t+2) P(n-1) + 2t(2-t) P'(n-1)` |
| `25_oeis_a110319_skip_denominator_derivative` | natural OEIS recurrence, skip prefix, denominator, two-lag derivative | `--skip-prefix 1 --min-rec-len 2 --max-rec-len 2 --max-var-deg 2 --max-idx-deg 2 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 2` | `(n+2) P(n) = (3n-2)t P(n-1) - 2t^2 P'(n-1) + (14-3n)t P(n-2) + 8t^2 P'(n-2)` |
| `26_oeis_a156289_second_derivative` | natural OEIS recurrence, second derivative | `--min-rec-len 1 --max-rec-len 1 --max-var-deg 2 --max-idx-deg 0 --max-diff-deg 2` | `P(n) = (1 + 3t) P(n-1) + (3t + 2t^2) P'(n-1) + t^2 P^(2)(n-1)` |
| `27_oeis_a219836_two_lag_derivative` | natural OEIS recurrence, two-lag derivative | `--min-rec-len 2 --max-rec-len 2 --max-var-deg 3 --max-idx-deg 1 --max-diff-deg 1` | `P(n) = (1 - 3t + nt) P(n-1) + (t - t^2) P'(n-1) + (t - 3t^2 + nt^2) P(n-2) + (t^2 - t^3) P'(n-2)` |
| `28_oeis_a089627_denominator_second_order` | natural OEIS recurrence, denominator, second order | `--min-rec-len 2 --max-rec-len 2 --max-var-deg 1 --max-idx-deg 1 --max-diff-deg 0 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1` | `(1 - n) P(n) = (3 - 2n) P(n-1) + (-2 + n + 8t - 4nt) P(n-2)` |
| `29_oeis_a101920_three_lag_denominator_derivative` | natural OEIS recurrence, denominator, three-lag derivative | `--min-rec-len 3 --max-rec-len 3 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1` | `(1 + n) P(n) = (4n - 3t + 3nt) P(n-1) + (2t - 2t^2) P'(n-1) + (8 - 4n + 3t - 3nt) P(n-2) + (2t - 2t^2) P'(n-2) + (-3 + n) P(n-3)` |
| `30_oeis_a177970_order8_closed_form` | natural OEIS recurrence, closed-form-derived, order 8 | `--min-rec-len 8 --max-rec-len 8 --min-var-deg 6 --max-var-deg 6 --max-idx-deg 0 --max-diff-deg 0` | `P(n) = (6 + 6t) P(n-1) + (-15t^2 - 24t - 15) P(n-2) + ... + (-t^6 + 4t^5 - 6t^4 + 4t^3 - t^2) P(n-8)` |
| `31_oeis_a321331_signed_derivative_lag` | natural OEIS recurrence, signed terms, derivative lag | `--min-rec-len 2 --max-rec-len 2 --max-var-deg 2 --max-idx-deg 0 --max-diff-deg 1` | `P(n) = (1 + 2t) P(n-1) + t P'(n-1) - t^2 P(n-2) - t^2 P'(n-2)` |
| `32_oeis_a390433_quadratic_index_second_derivative` | natural OEIS recurrence, quadratic index, second derivative | `--min-rec-len 1 --max-rec-len 1 --max-var-deg 2 --max-idx-deg 2 --max-diff-deg 2` | `P(n) = (4 - 4n + n^2 + t) P(n-1) + (-3t + 2nt) P'(n-1) + t^2 P^(2)(n-1)` |
