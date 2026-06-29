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
`projects/real-rooted-oeis/sequences` queue and related recurrence workbench
evidence. Fixtures `46` onward are deliberately harder high-weight OEIS
recurrences.

Regenerate from the Rust workspace root with:

```sh
cargo run -p polynomial-tools --example generate_recurrence_benchmarks
bash polynomial-tools/fixtures/recurrence-benchmarks/regenerate-json.sh
```

Run the built-in timing suite:

```sh
polytool bench recurrence-fixtures --repeat 3
polytool bench recurrence-fixtures --only 23_sparse --repeat 5
polytool bench recurrence-fixtures --only oeis --repeat 3 \
  --summary --report bench-results/recurrence-fixtures/oeis.md
polytool bench recurrence-fixtures --only oeis --repeat 3 --format json \
  > bench-results/recurrence-fixtures/oeis.json
polytool bench compare old.json new.json --top 10
```

The default output is per-run TSV. The optional `--summary` flag appends
fixture-level and category-level TSV summaries. The optional `--report <path>`
writes the same benchmark run as a Markdown report suitable for checked-in
benchmark notes or project handoffs.
The JSON output includes per-run records, fixture/category summaries, and
adaptive search diagnostics such as generated candidates, exact solve attempts,
modular-prefilter rejections, and held-out verification failures. The compare
subcommand reads two JSON benchmark outputs and reports speedups, slowdowns,
new fixtures, removed fixtures, and worst regressions.

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
| `33_oeis_a008292_eulerian_variant_derivative` | natural OEIS recurrence, Eulerian-like first derivative | `--min-rec-len 1 --max-rec-len 1 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 1` | `P(n) = (1 - t + nt) P(n-1) + (t - t^2) P'(n-1)` |
| `34_oeis_a019538_derivative_appell` | natural OEIS recurrence, first derivative | `--min-rec-len 1 --max-rec-len 1 --max-var-deg 2 --max-idx-deg 0 --max-diff-deg 1` | `P(n) = (1 + 2t) P(n-1) + (t + t^2) P'(n-1)` |
| `35_oeis_a035607_two_lag_constant_t` | natural OEIS recurrence, two-lag t-dependent coefficients | `--min-rec-len 2 --max-rec-len 2 --max-var-deg 1 --max-idx-deg 0 --max-diff-deg 0` | `P(n) = (1 + t) P(n-1) + t P(n-2)` |
| `36_oeis_a049218_quadratic_index_lag` | natural OEIS recurrence, quadratic index, second order | `--min-rec-len 2 --max-rec-len 2 --max-var-deg 1 --max-idx-deg 2 --max-diff-deg 0` | `P(n) = t P(n-1) + (-n - n^2) P(n-2)` |
| `37_oeis_a088729_second_derivative_sturm` | natural OEIS recurrence, second derivative | `--min-rec-len 1 --max-rec-len 1 --max-var-deg 1 --max-idx-deg 0 --max-diff-deg 2` | `P(n) = (3 + t) P(n-1) + (4 + 3t) P'(n-1) + 2t P^(2)(n-1)` |
| `38_oeis_a104684_denominator_second_order` | natural OEIS recurrence, denominator, second order | `--min-rec-len 2 --max-rec-len 2 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 0 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1` | `(1 - n) P(n) = (6 - 4n + 3t - 2nt) P(n-1) + (-2t^2 + nt^2) P(n-2)` |
| `39_oeis_a130749_denominator_two_lag_derivative` | natural OEIS recurrence, denominator, two-lag derivative | `--min-rec-len 2 --max-rec-len 2 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1` | `(1 + 1/2n) P(n) = (1 + 3/2n + t + 3/2nt) P(n-1) + (t - t^2) P'(n-1) + (-n - 3/2nt) P(n-2)` |
| `40_oeis_a257142_denominator_three_lag_derivative` | natural OEIS recurrence, denominator, three-lag derivative | `--min-rec-len 3 --max-rec-len 3 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1` | `(1 + 1/3n) P(n) = (7/3 + 4/3n + t + nt) P(n-1) + ... + (-2/3 + 1/3n) P(n-3)` |
| `41_oeis_a348576_two_lag_derivative_dense` | natural OEIS recurrence, dense two-lag derivative | `--min-rec-len 2 --max-rec-len 2 --max-var-deg 3 --max-idx-deg 1 --max-diff-deg 1` | `P(n) = (1 + 2t + nt) P(n-1) + (t + t^2) P'(n-1) + (-t^2 - nt^2) P(n-2) + (-t^2 - nt^2 - t^3 - nt^3) P'(n-2)` |
| `42_oeis_a375853_denominator_derivative` | natural OEIS recurrence, denominator, first derivative | `--min-rec-len 1 --max-rec-len 1 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1` | `(1 - n) P(n) = (-2 - n + 2t - 3nt) P(n-1) + (-2t + 2t^2) P'(n-1)` |
| `43_oeis_a395454_two_lag_index_derivative` | natural OEIS recurrence, two-lag index derivative | `--min-rec-len 2 --max-rec-len 2 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 1` | `P(n) = (-2 + n - 2t + nt) P(n-1) + (t - t^2) P'(n-1) + (-4t + 2nt) P(n-2)` |
| `44_oeis_a109954_third_order_constant` | natural OEIS recurrence, third order t-dependent coefficients | `--min-rec-len 3 --max-rec-len 3 --max-var-deg 1 --max-idx-deg 0 --max-diff-deg 0` | `P(n) = (-3 + t) P(n-1) + (-3 + t) P(n-2) - P(n-3)` |
| `45_oeis_a176231_second_derivative_signed` | natural OEIS recurrence, signed second derivative | `--min-rec-len 1 --max-rec-len 1 --max-var-deg 1 --max-idx-deg 0 --max-diff-deg 2` | `P(n) = (-1 + t) P(n-1) + (2 - 4t) P'(n-1) + 4t P^(2)(n-1)` |
| `46_oeis_a173882_high_weight_derivative` | natural OEIS recurrence, high weighted unknown count, sixth-order derivative denominator | `--min-rec-len 6 --max-rec-len 6 --max-var-deg 3 --max-idx-deg 2 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 2 --fit-extra-rows 20` | `(1 + 2/3n - 1/3n^2) P(n) = (-1/9 + 14/3n - 11/9n^2 - 343/27t + 436/27nt - 131/27n^2t) P(n-1) + (17/9t - 17/9nt - 289/27t^2 + 67/9nt^2) P'(n-1) + (82/9 - 28/3n + 14/9n^2 + 2642/27t - 677/9nt + 131/9n^2t - 55/3t^2 + 569/27nt^2 - 146/27n^2t^2) P(n-2) + (-116/9t + 56/9nt + 1415/27t^2 - 547/27nt^2 - 101/27t^3 + 23/9nt^3) P'(n-2) + (-14 + 20/3n - 2/3n^2 - 4840/27t + 2738/27nt - 400/27n^2t + 1246/9t^2 - 878/9nt^2 + 148/9n^2t^2 - 16/27t^3 + 4/9nt^3 - 2/27n^2t^3) P(n-3) + (70/3t - 22/3nt - 1720/27t^2 + 430/27nt^2 + 530/27t^3 - 212/27nt^3) P'(n-3) + (55/9 - 2/3n - 1/9n^2 + 2578/27t - 1216/27nt + 152/27n^2t - 794/3t^2 + 1216/9nt^2 - 152/9n^2t^2 + 10/3t^3 - 16/9nt^3 + 2/9n^2t^3) P(n-4) + (-140/9t + 32/9nt + 316/27t^2 - 2/3nt^2 - 88/3t^3 + 74/9nt^3) P'(n-4) + (-1/9 - 2/3n + 1/9n^2 - 317/27t + 6nt - 7/9n^2t + 1406/9t^2 - 1670/27nt^2 + 164/27n^2t^2 - 16/3t^3 + 20/9nt^3 - 2/9n^2t^3) P(n-5) + (29/9t - 5/9nt + 313/27t^2 - 71/27nt^2 + 398/27t^3 - 28/9nt^3) P'(n-5) + (280/27t - 89/27nt + 7/27n^2t - 35/3t^2 + 29/9nt^2 - 2/9n^2t^2 + 70/27t^3 - 8/9nt^3 + 2/27n^2t^3) P(n-6) + (-35/27t^2 + 5/27nt^2 - 35/27t^3 + 5/27nt^3) P'(n-6)` |
| `47_oeis_a174148_high_weight_second_derivative` | natural OEIS recurrence, high weighted unknown count, second derivatives, denominator | `--min-rec-len 2 --max-rec-len 2 --max-var-deg 3 --max-idx-deg 3 --max-diff-deg 2 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 2 --fit-extra-rows 20` | `(1 + n + 1/4n^2) P(n) = (1 - 15/4n - 4n^2 + t + 21/4nt + 23/4n^2t + 3/2n^3t) P(n-1) + (69/2t + 63/2nt - 63/4n^2t - 3/2t^2 - 9nt^2 - 15/4n^2t^2) P'(n-1) + (45/2t^2 + 117/4nt^2 + 9/2t^3 + 9/4nt^3) P^(2)(n-1) + (19/4n + 17/4n^2 + 143/4nt - 71n^2t + 21/4n^3t - 1/2nt^2 + 11/4n^2t^2 + 3/4n^3t^2) P(n-2) + (117/4nt + 63/4n^2t + 63nt^2 - 39n^2t^2 - 9/4nt^3 - 3/4n^2t^3) P'(n-2) + (18nt^2 + 18nt^3) P^(2)(n-2)` |
| `48_oeis_a176200_high_weight_eulerian_transform` | natural OEIS recurrence, high weighted unknown count, shifted Eulerian transform, derivative denominator | `--min-rec-len 3 --max-rec-len 3 --max-var-deg 3 --max-idx-deg 3 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 2 --fit-extra-rows 20` | `(1 - 3/2n + 1/2n^2) P(n) = (1 - 5/2n + n^2 - t + 3/2nt - 3/2n^2t + 1/2n^3t) P(n-1) + (t - 3/2nt + 1/2n^2t - t^2 + 3/2nt^2 - 1/2n^2t^2) P'(n-1) + (n - 1/2n^2 - 1/2nt + n^2t - 1/2n^3t - 2nt^2 + 2n^2t^2 - 1/2n^3t^2) P(n-2) + (nt - 1/2n^2t - nt^3 + 1/2n^2t^3) P'(n-2) + (-1/2nt + 1/2n^2t + 3/2nt^2 - 2n^2t^2 + 1/2n^3t^2) P(n-3) + (-1/2nt^2 + 1/2n^2t^2 + 1/2nt^3 - 1/2n^2t^3) P'(n-3)` |
| `49_oeis_a176204_high_weight_shifted_eulerian_transform` | natural OEIS recurrence, high weighted unknown count, shifted Eulerian transform, derivative denominator | `--max-rec-len 8 --max-var-deg 3 --max-idx-deg 3 --max-diff-deg 2 --denominator --max-denom-var-deg 1 --max-denom-idx-deg 2 --inhomogeneous --max-inhomo-var-deg 3 --max-inhomo-idx-deg 3 --fit-extra-rows 20` | `(1 - 3/2n + 1/2n^2) P(n) = (1 - 5/2n + n^2 - t + 3/2nt - 3/2n^2t + 1/2n^3t) P(n-1) + (t - 3/2nt + 1/2n^2t - t^2 + 3/2nt^2 - 1/2n^2t^2) P'(n-1) + (n - 1/2n^2 - 1/2nt + n^2t - 1/2n^3t - 2nt^2 + 2n^2t^2 - 1/2n^3t^2) P(n-2) + (nt - 1/2n^2t - nt^3 + 1/2n^2t^3) P'(n-2) + (-1/2nt + 1/2n^2t + 3/2nt^2 - 2n^2t^2 + 1/2n^3t^2) P(n-3) + (-1/2nt^2 + 1/2n^2t^2 + 1/2nt^3 - 1/2n^2t^3) P'(n-3)` |
| `50_oeis_a155495_high_weight_binomial_product` | natural OEIS recurrence, high weighted unknown count, binomial product, derivative denominator | `--min-rec-len 3 --max-rec-len 3 --max-var-deg 3 --max-idx-deg 2 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 2 --fit-extra-rows 20` | `(1 + 3n + 2n^2) P(n) = (1 - 11n + 30n^2 + t + 37nt - 24n^2t) P(n-1) + (48t - 54nt - 48t^2 + 54nt^2) P'(n-1) + (-12 + 34n - 36n^2 + 24t - 124nt + 144n^2t + 36t^2 - 14nt^2 - 36n^2t^2) P(n-2) + (-48t + 48t^3) P'(n-2) + (12 - 20n + 8n^2 + 36t - 60nt + 24n^2t + 36t^2 - 60nt^2 + 24n^2t^2 + 12t^3 - 20nt^3 + 8n^2t^3) P(n-3)` |
| `51_oeis_a168287_high_weight_pascal_transform` | natural OEIS recurrence, high weighted unknown count, Pascal-transform family, derivative denominator | `--min-rec-len 3 --max-rec-len 3 --max-var-deg 3 --max-idx-deg 2 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1 --fit-extra-rows 20` | `(1 - 1/3n) P(n) = (5/3 - 2/3n - 1/3t + nt - 1/3n^2t) P(n-1) + (t - 1/3nt - t^2 + 1/3nt^2) P'(n-1) + (-2/3 + 1/3n + 8/3t - 8/3nt + 2/3n^2t + 4/3t^2 - 4/3nt^2 + 1/3n^2t^2) P(n-2) + (-2/3t + 1/3nt + 2/3t^3 - 1/3nt^3) P'(n-2) + (-2t + 5/3nt - 1/3n^2t - 2t^2 + 5/3nt^2 - 1/3n^2t^2) P(n-3)` |
| `52_oeis_a168288_high_weight_pascal_transform` | natural OEIS recurrence, high weighted unknown count, Pascal-transform family, derivative denominator | `--min-rec-len 3 --max-rec-len 3 --max-var-deg 3 --max-idx-deg 2 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1 --fit-extra-rows 20` | `(1 - 1/3n) P(n) = (5/3 - 2/3n - 1/3t + nt - 1/3n^2t) P(n-1) + (t - 1/3nt - t^2 + 1/3nt^2) P'(n-1) + (-2/3 + 1/3n + 8/3t - 8/3nt + 2/3n^2t + 4/3t^2 - 4/3nt^2 + 1/3n^2t^2) P(n-2) + (-2/3t + 1/3nt + 2/3t^3 - 1/3nt^3) P'(n-2) + (-2t + 5/3nt - 1/3n^2t - 2t^2 + 5/3nt^2 - 1/3n^2t^2) P(n-3)` |
| `53_oeis_a168289_high_weight_pascal_transform` | natural OEIS recurrence, high weighted unknown count, Pascal-transform family, derivative denominator | `--min-rec-len 3 --max-rec-len 3 --max-var-deg 3 --max-idx-deg 2 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1 --fit-extra-rows 20` | `(1 - 1/3n) P(n) = (5/3 - 2/3n - 1/3t + nt - 1/3n^2t) P(n-1) + (t - 1/3nt - t^2 + 1/3nt^2) P'(n-1) + (-2/3 + 1/3n + 8/3t - 8/3nt + 2/3n^2t + 4/3t^2 - 4/3nt^2 + 1/3n^2t^2) P(n-2) + (-2/3t + 1/3nt + 2/3t^3 - 1/3nt^3) P'(n-2) + (-2t + 5/3nt - 1/3n^2t - 2t^2 + 5/3nt^2 - 1/3n^2t^2) P(n-3)` |
| `54_oeis_a168290_high_weight_pascal_transform` | natural OEIS recurrence, high weighted unknown count, Pascal-transform family, derivative denominator | `--min-rec-len 3 --max-rec-len 3 --max-var-deg 3 --max-idx-deg 2 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1 --fit-extra-rows 20` | `(1 - 1/3n) P(n) = (5/3 - 2/3n - 1/3t + nt - 1/3n^2t) P(n-1) + (t - 1/3nt - t^2 + 1/3nt^2) P'(n-1) + (-2/3 + 1/3n + 8/3t - 8/3nt + 2/3n^2t + 4/3t^2 - 4/3nt^2 + 1/3n^2t^2) P(n-2) + (-2/3t + 1/3nt + 2/3t^3 - 1/3nt^3) P'(n-2) + (-2t + 5/3nt - 1/3n^2t - 2t^2 + 5/3nt^2 - 1/3n^2t^2) P(n-3)` |
| `55_oeis_a319251_narayana_correction_derivative_denominator` | natural OEIS recurrence, Narayana correction, three-lag derivative denominator | `--min-rec-len 3 --max-rec-len 3 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1 --fit-extra-rows 5` | `(1 + n) P(n) = (19/2 - 5/2n - 5/2t + 5/2nt) P(n-1) + (11t - t^2) P'(n-1) + (-33/2 + 11/2n + 19t - 7nt + 3/2t^2 - 1/2nt^2) P(n-2) + (-6t + 2t^2) P'(n-2) + (8 - 2n - 8t + 2nt) P(n-3)` |
| `56_oeis_a162303_three_lag_derivatives_denominator` | natural OEIS recurrence, three lags with derivatives, denominator | `--min-rec-len 3 --max-rec-len 3 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1 --fit-extra-rows 5` | `(1 + n) P(n) = (-1 - 3t + 3nt) P(n-1) + (2t - 2t^2) P'(n-1) + (3 + 2n - 9t + 6nt) P(n-2) + (6t - 6t^2) P'(n-2) + (2 + n - 6t + 3nt) P(n-3) + (4t - 4t^2) P'(n-3)` |
| `57_oeis_a119307_second_derivative_denominator` | natural OEIS recurrence, second derivatives, denominator | `--min-rec-len 2 --max-rec-len 2 --max-var-deg 2 --max-idx-deg 2 --max-diff-deg 2 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 2 --fit-extra-rows 8` | `(1 + 3/2n + 1/2n^2) P(n) = (1 + 5n + 5/2n^2 + t + 3nt + 3/2n^2t) P(n-1) + (-3/2 - 7t - 2nt - 3/2t^2 - nt^2) P'(n-1) + (-3/2t - 3/2t^2) P^(2)(n-1) + (-3/2 - 2n - 1/2n^2 + 7/2t + 5nt + 3/2n^2t) P(n-2) + (-t - nt + t^2 + nt^2) P'(n-2)` |
| `58_oeis_a108558_type_d_no_derivative_denominator` | natural OEIS recurrence, type-D h-vector, no derivatives, denominator | `--max-rec-len 4 --max-var-deg 3 --max-idx-deg 2 --max-diff-deg 0 --denominator --fit-extra-rows 16` | `(1 - 9/20n + 1/20n^2) P(n) = (11/4 - 13/10n + 3/20n^2 + 11/4t - 13/10nt + 3/20n^2t) P(n-1) + (-5/2 + 5/4n - 3/20n^2 - t + 7/10nt - 1/10n^2t - 5/2t^2 + 5/4nt^2 - 3/20n^2t^2) P(n-2) + (3/4 - 2/5n + 1/20n^2 - 3/4t + 2/5nt - 1/20n^2t - 3/4t^2 + 2/5nt^2 - 1/20n^2t^2 + 3/4t^3 - 2/5nt^3 + 1/20n^2t^3) P(n-3)` |
| `59_oeis_a142706_derivative_coefficients_eulerian` | natural OEIS recurrence, derivative coefficients of Eulerian polynomials, denominator | `--min-rec-len 2 --max-rec-len 2 --max-var-deg 3 --max-idx-deg 3 --max-diff-deg 1 --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1 --fit-extra-rows 5` | `(1 + n) P(n) = (4 + 3n + 2t + 4nt + 2n^2t) P(n-1) + (t + nt - t^2 - nt^2) P'(n-1) + (-4 - 2n - 4nt - 2n^2t - 2nt^2 - 3n^2t^2 - n^3t^2) P(n-2) + (-2t - nt - 2nt^2 - n^2t^2 + 2t^3 + 3nt^3 + n^2t^3) P'(n-2)` |
