# Non-Nesting Rooks Worklog

This note tracks what currently works, what fails, and what looks like the
best next proof targets.

## Main objects

- `R_mu(t)`: non-nesting rook polynomial on a Ferrers board `mu`.
- `F_c^nu = R_{nu_<c>}`: column-strip polynomials.
- `Delta_c^nu = (F_c^nu - F_{c+1}^nu)/t`.
- For `mu = (mu', m)`:
  - `G(c,s) = F_c^{mu'} + t * sum_{j=c+1}^s F_j^{mu'}`.
  - `R_{mu_<c>} = G(c,m)`.

## What works

### 1. Strict Ferrers-poset reformulation

- A non-nesting rook placement is exactly a chain in the strict Ferrers-cell
  poset
  `P_mu = {(i,j)}`
  with order `(i,j) < (i',j')` iff `i < i'` and `j > j'`.
- So `R_mu(t)` is the chain polynomial of `P_mu`.

### 2. Fixed-row LGV determinant

- For fixed rows `I = {i_1 < ... < i_k}`, the number of placements using
  exactly those rows is
  `d_I(mu) = det[ binom(mu_{i_{k+1-b}} - b + 1, a - b + 1) ]_(a,b=1)^k`.
- Hence
  `r_k(mu) = sum_I d_I(mu)`.
- Verified directly by a brute-force Python check for all partitions of size
  `<= 10`.

### 3. Delta-grid package

The clean recurrence package is:

- `Delta_c^mu = Delta_c^{mu'} + F_{c+1}^{mu'}`.
- `G(c,s) = G(c+1,s) + t Delta_c^mu`.
- `G(c,s+1) = G(c,s) + t F_{s+1}^{mu'}`.

Empirically this is very strong.

- `nn_rook_grid.rs` passed for all partitions with `|mu| <= 20`.
- In particular:
  - `Delta_c^mu` is real-rooted and coefficientwise nonnegative.
  - `Delta_c^mu << G(c+1,s)`.
  - `Delta_c^mu << G(c,s)`.
  - `G(c+1,s) << G(c,s)`.
  - Both diagonal relations also hold.

### 4. q/deformation evidence

- `nn_rook_qdeform2.rs`:
  - for `q in {0, 1/4, 1/2, 3/4, 1}`, the `q^{nest}` specialization stayed
    real-rooted in all tested cases;
  - `G(c+1,s) << G(c,s)` held for all tested `q`.
- `nn_rook_deform.rs`:
  - along the path
    `A(alpha) = G(c+1,s) + alpha * t * R`,
    `B(alpha) = G(c,s) + alpha * t * R`,
    interlacing survived at all sampled points (`33633/33633` checks).

### 5. 2x2 square structure

With

- `A = G(c+1,s)`,
- `B = G(c,s) = A + t Delta`,
- `R = F_{s+1}^{mu'}`,
- `A' = A + tR = G(c+1,s+1)`,
- `B' = B + tR = G(c,s+1)`,

the following always held in `nn_rook_square.rs` for all squares with
`|mu| <= 14`:

- rows:
  - `A << B`,
  - `A' << B'`;
- diagonal:
  - `A << B'`;
- anti-diagonal:
  - `A' << B`;
- input relations:
  - `Delta << A`,
  - `Delta << B`,
  - `Delta + R << A`,
  - `Delta - R` has nonnegative coefficients,
  - `Delta - R << A'`,
  - `Delta - R << B`.

This gives two experimentally clean side-lemmas:

- diagonal:
  `B' = A + t(Delta + R)` and `Delta + R << A`;
- anti-diagonal:
  `B = A' + t(Delta - R)` and `Delta - R << A'`.

Equivalently,

- `(G(c,s+1) - G(c+1,s))/t = Delta_c^mu + F_{s+1}^{mu'}`,
- `(G(c,s) - G(c+1,s+1))/t = Delta_c^mu - F_{s+1}^{mu'}`.

## What fails

### 1. Easy Branden-family proof

- The naive pairwise-interlacing / common-interlacing route does not work.
- Degree gaps break the termwise Wagner argument.

### 2. Easy graph-class proof

- The conflict / incomparability graph is not claw-free in general.
- It is not interval/chordal either.

### 3. Easy vertical square explanation

- `R << A` and `R << B` fail often.
- Accordingly, the vertical edges
  - `A << A'`,
  - `B << B'`
  fail often too.

So the square is not controlled by a naive repeated application of the standard
lemma `f << g => g << g + t f`.

## Best next proof targets

### Route A: Global LGV/TNN packaging

Turn the fixed-row determinants into one global total-nonnegative object.
This would give a structural proof from the poset/chain-polynomial side.

### Route B: Coupled square lemmas

Try to prove the two diagonal quotient statements:

- `Delta_c^mu + F_{s+1}^{mu'} << G(c+1,s)`,
- `Delta_c^mu - F_{s+1}^{mu'} << G(c+1,s+1)`.

If both are true, the diagonal and anti-diagonal edges become immediate from
the standard additive step.

### Route C: q-family deformation

Try to show that the interlacing property persists for the `q^{nest}` family
and can be specialized from `q = 1` (standard rook) to `q = 0`
(non-nesting rook).

## Related Touchard / big-block thread

- The staircase standard rook polynomial is the reversed Touchard polynomial:
  `R^{std}_{delta_{n-1}}(t) = t^n B_n(1/t)`.
- This motivated adding a related section to the draft on:
  - Touchard/Bell polynomials,
  - the multivariate Bell polynomial,
  - big-block polynomials `P_{n,j}(t)`,
  - and the derivative-style recurrence for `Q_{n,j}(t) = t P_{n,j}(t)`.

What is now proved/documented in the draft:

- `P_{n+1,j}(t) = sum_{k=0}^{j-2} binom(n,k) P_{n-k,j}(t)
  + t sum_{k=j-1}^{n} binom(n,k) P_{n-k,j}(t)`.
- If `Q_{n,j}(t) = t P_{n,j}(t)`, then
  `Q_{n+1,j}(t) = t Q'_{n,j}(t)
  + sum_{k=1}^{j-2} binom(n,k) Q_{n-k,j}(t)
  + t binom(n,j-1) Q_{n-j+1,j}(t)`.
- The clean proof of the `Q`-recurrence uses the exponential generating
  function
  `sum_n P_{n,j}(t) z^n/n! = exp(sum_{r=1}^{j-1} z^r/r! + t sum_{r=j}^\infty z^r/r!)`.

Checked separately by brute-force set partitions:

- both `P`- and `Q`-recurrences hold for `2 <= j <= 5` and `n <= 7`.

Open direction:

- find a natural rook statistic on staircase placements whose distribution is
  `P_{n,j}` or `Q_{n,j}`.

## Useful commands

```bash
cargo check --offline -p experiments --bin nn_rook_grid
CARGO_TARGET_DIR=/tmp/rust-target cargo run --offline --quiet -p experiments --bin nn_rook_grid -- 20

cargo check --offline -p experiments --bin nn_rook_square
CARGO_TARGET_DIR=/tmp/rust-target cargo run --offline --quiet -p experiments --bin nn_rook_square

CARGO_TARGET_DIR=/tmp/rust-target cargo run --offline --quiet -p experiments --bin nn_rook_qdeform2
CARGO_TARGET_DIR=/tmp/rust-target cargo run --offline --quiet -p experiments --bin nn_rook_deform

latexmk -pdf -interaction=nonstopmode -halt-on-error documents/non-nesting-rooks.tex
```
