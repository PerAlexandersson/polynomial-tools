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

What now also works:

- There is an explicit staircase-rook bijection behind the Touchard identity.
  After reversing columns, the staircase board is the triangular board
  `T_n = {(i,j) : 1 <= i < j <= n}`.
- A standard rook placement on `T_n` is a set of directed edges `i -> j` with
  no repeated source and no repeated target, hence a disjoint union of directed
  paths.
- The vertex sets of those path components form a set partition of `[n]`, and
  the inverse map sends a block `{b_1 < ... < b_m}` to the rook chain
  `(b_1,b_2), (b_2,b_3), ..., (b_{m-1},b_m)`.
- Under this bijection, `bb_j` is exactly the number of path components with at
  least `j` vertices, equivalently the number of maximal rook chains with at
  least `j-1` rooks.
- So `P_{n,j}(t)` is already a natural staircase standard-rook generating
  polynomial; the old "find a rook statistic" problem is solved.

Verified directly by the dedicated checker `touchard_staircase.rs`:

- bijection counts agree for `n <= 8`,
- the map is onto and the inverse is correct on all tested placements,
- the induced big-block distributions match set partitions for `2 <= j <= 5`
  and `n <= 8`.

Open direction:

- turn the big-block recurrence / rook-chain model into an actual
  real-rootedness or interlacing proof.

## General Ferrers path-component extension

Tried extending the staircase "union of directed paths" model to an arbitrary
Ferrers board `mu = (mu_1 >= ... >= mu_ell)` as follows:

- set `n(mu) = max_i (i + mu_i)`,
- reverse columns inside the width-`n(mu)-1` rectangle,
- interpret a rook in row `i` and original column `c` as an edge
  `i -> n(mu)-c+1`.

Then every edge satisfies `i < n(mu)-c+1`, so every standard rook placement
becomes a directed graph on `[n(mu)]` with indegree/outdegree at most `1`,
hence a disjoint union of directed paths.

Define `P_{mu,j}(t)` to be the generating polynomial for the statistic

- `rch_j(rho)` = number of path components with at least `j` vertices.

What the experiments say:

- For `j=2`, all scanned `P_{mu,2}(t)` are real-rooted for every Ferrers shape
  with `|mu| <= 20` (`2713/2713`).
- For `j=3` and `j=4`, all scanned `P_{mu,j}(t)` are real-rooted for every
  Ferrers shape with `|mu| <= 18` (`1596/1596` in both cases).
- Row-deletion interlacing fails early, so the real-rootedness does not seem to
  come from a naive Sturm sequence under `mu -> most(mu)`.
  Small failures:
  - `j=2`: `mu=[2,2,1]`, where `P_[2,2]=1+4t+2t^2` and
    `P_[2,2,1]=1+7t+2t^2` do not interlace.
  - `j=3`: `mu=[2,1,1]`, where `4+t` and `6+t` do not interlace in the
    directed same-degree sense.
- The statistic is often genuinely different from the ordinary rook polynomial:
  only `66/2713` shapes with `|mu| <= 20` had `P_{mu,2}` equal to the standard
  rook polynomial.
- But it also becomes trivial for wide/rectangular-ish shapes. Example:
  `mu=[4,4,4]` has `P_{mu,3}(t)=73`, so no placement supports a path component
  with at least `3` vertices.

So the extension looks much better as a real-rootedness phenomenon than as an
interlacing/Sturm phenomenon.

## Ordered set partitions with big blocks

Also tested the ordered-set-partition analogue:

- `O_{n,j}(t) = sum t^{bb_j(pi)}` over ordered set partitions of `[n]`.

Useful generating function:

- if `A_j(z,t) = sum_{r=1}^{j-1} z^r/r! + t sum_{r=j}^\infty z^r/r!`, then
  `sum_n O_{n,j}(t) z^n/n! = 1 / (1 - A_j(z,t))`.

Experimentally this looks even stronger than the unordered case.
Using the refined-by-(blocks, big blocks) dynamic program in
`ordered_big_blocks.rs`, the exact Bézout checks give:

- for each `j = 2,3,4,5,6`, all rows `O_{n,j}(t)` with `n <= 18` are real-rooted,
- for each `j = 2,3,4,5,6`, the consecutive rows `O_{n,j} << O_{n+1,j}` hold
  for every tested `n <= 17`.

Sample rows:

- `j=2`:
  - `n=4`: `24 + 45t + 6t^2`
  - `n=5`: `120 + 311t + 110t^2`
  - `n=6`: `720 + 2383t + 1490t^2 + 90t^3`
- `j=3`:
  - `n=4`: `66 + 9t`
  - `n=5`: `450 + 91t`
  - `n=6`: `3690 + 973t + 20t^2`

At the moment this is only computational evidence, but it strongly suggests
that ordered set partitions may be the cleaner family from the real-rootedness
point of view.

### Ordered refinements by block count and by position of `n`

I then pushed the ordered family one level further in two different ways.

1. Exact number of blocks:

- Define `O_{n,j}^{<m>}(t)` to count ordered set partitions of `[n]` with
  exactly `m` blocks, weighted by `bb_j`.
- Then the block containing `n+1` gives the clean recurrence
  `O_{n+1,j}^{<m>}(t)
   = m (sum_{k=0}^{j-2} binom(n,k) O_{n-k,j}^{<m-1>}(t)
      + t sum_{k=j-1}^n binom(n,k) O_{n-k,j}^{<m-1>}(t))`.
- In bivariate form, for
  `G_{n,j}(u,t) = sum_m O_{n,j}^{<m>}(t) u^m`,
  this becomes
  `G_{n+1,j}(u,t)
   = u(1+u d/du)(sum_{k=0}^{j-2} binom(n,k) G_{n-k,j}(u,t)
      + t sum_{k=j-1}^n binom(n,k) G_{n-k,j}(u,t))`.

2. Position of the block containing `n`:

- Let `H_{n,j,p}(t)` count ordered set partitions where the block containing
  `n` is in position `p`, and put
  `H_{n,j}(x,t) = sum_{p>=1} H_{n,j,p}(t) x^{p-1}`.
- With
  `A_j(z,t) = sum_{r=1}^{j-1} z^r/r! + t sum_{r=j}^\infty z^r/r!`,
  the marked-block decomposition gives
  `sum_{n>=1} H_{n,j}(x,t) z^{n-1}/(n-1)!
   = A'_j(z,t) / ((1-A_j(z,t))(1-xA_j(z,t)))`.

The good news is that these really are cleaner than the raw `O_{n,j}` family.
The bad news is that they do not seem to become stable in any naive bivariate
sense.

Using `ordered_big_blocks_position.rs`, I checked:

- for each `j = 2,3,4,5,6`, all coordinate specializations
  `G_{n,j}(c,t)` and `H_{n,j}(c,t)` with `c = 1,2,3` are real-rooted for
  every tested `n <= 14`;
- but same-phase line tests fail very early.

Small explicit failures:

- `H_{3,2}(s,s) = 2 + 6s + 5s^2`, which has negative discriminant;
- `G_{4,2}(s,s) = s^2 + 8s^3 + 66s^4 = s^2(1 + 8s + 66s^2)`, whose quadratic
  factor is not real-rooted.

So the ordered refinements look promising as recurrence packages, but not as
direct stable-polynomial lifts.

### Matrix form and interlacing tests for the ordered refinements

The exact-block refinement does admit a clean matrix/operator form.
If
`v_n = (O_{n,j}^{<0>}(t), O_{n,j}^{<1>}(t), O_{n,j}^{<2>}(t), ...)^\top`,
then

- `v_{n+1} = J ( sum_{k=0}^{j-2} binom(n,k) v_{n-k}
             + t sum_{k=j-1}^n binom(n,k) v_{n-k} )`,

where `J` is the weighted shift matrix with entries `J_{m,m-1} = m`.
Equivalently, on the ordinary generating variable `u`, `J` is the operator
`u(1+u d/du)`.

So there is definitely a matrix picture, but it is not the same as the usual
Brändén staircase-matrix setup: one step in `n` uses a whole binomially weighted
sum of earlier vectors, not just one previous vector.

I also checked whether the natural coefficient vectors are interlacing families.

- For the exact-block vector in the natural order `m=1,2,...`, adjacent
  interlacing fails often.
- If one reverses the order to `m=n,n-1,...,1`, then every *eligible adjacent*
  pair interlaces in the tested range:
  - `j=2`: `55/55`
  - `j=3`: `73/73`
  - `j=4`: `81/81`
  - `j=5`: `86/86`
  - `j=6`: `88/88`
- But this still does not make the full vector a Brändén interlacing sequence,
  because the degrees are not globally compatible.
  Example: for `j=2, n=5`, the reversed row contains
  `120, 240t, 60t+90t^2, 10t+20t^2, t`;
  the constant polynomial `120` and the quadratic `60t+90t^2` differ in degree
  by `2`, so the full row cannot be pairwise interlacing in the usual sense.

For the position-of-`n` refinement the situation is weaker.

- In the natural order `p=1,2,...`, adjacent interlacing fails badly.
- In the reversed order `p=n,n-1,...,1`, things improve but still do not become
  clean:
  - `j=2`: `21/91`
  - `j=3`: `48/91`
  - `j=4`: `62/91`
  - `j=5`: `72/91`
  - `j=6`: `82/91`

So the current picture is:

- exact-block refinement: clean matrix/operator, strong reversed-adjacent
  interlacing, but not a full interlacing family;
- position refinement: useful formulas, but not a clean interlacing family;
- therefore no direct application yet of Brändén's matrix-preservation theorem.

### Cumulative tails of the position vector

The next natural transform was to replace the raw position coefficients by
cumulative tails

- `T_{n,j,p}(t) = sum_{q >= p} H_{n,j,q}(t)`.

This has a clean generating-function description: since
`sum_n H_{n,j,p}(t) z^{n-1}/(n-1)! = A'_j(z,t) A_j(z,t)^{p-1} / (1-A_j(z,t))`,
the tail sums satisfy

- `sum_n T_{n,j,p}(t) z^{n-1}/(n-1)! = A'_j(z,t) A_j(z,t)^{p-1} / (1-A_j(z,t))^2`.

Empirically this is the first really promising position-based family.

For the reversed tail vector
`(T_{n,j,n}, T_{n,j,n-1}, ..., T_{n,j,1})`,
adjacent interlacing checks gave:

- `n <= 14`:
  - `j=2`: `83/91`
  - `j=3`: `88/91`
  - `j=4`: `91/91`
  - `j=5`: `91/91`
  - `j=6`: `91/91`
- `n <= 18`:
  - `j=2`: `128/153`
  - `j=3`: `138/153`
  - `j=4`: `145/153`
  - `j=5`: `153/153`
  - `j=6`: `153/153`
  - `j=7`: `153/153`
  - `j=8`: `153/153`

So the large-`j` behavior is strikingly better here.

Smallest failures:

- `j=4`: first failure at `n=15`, comparing
  `T_{15,4,7}` and `T_{15,4,6}`;
- `j=3`: first failure at `n=13`, comparing
  `T_{13,3,8}` and `T_{13,3,7}`;
- `j=2`: failures already appear earlier and more often.

The head sums `sum_{q <= p} H_{n,j,q}(t)` were much worse and do not look
useful.

At the moment, the cumulative tails are the best lead on the position side.

One important correction: the strong positive statement above is only about
adjacent pairs in the reversed tail vector.
If one asks for a full interlacing sequence in the usual sense, i.e. for
`T_{n,j,p} << T_{n,j,q}` whenever `p > q`, that is still false in the tested
range.

For `n <= 18`, the pairwise results on eligible pairs were:

- `j=2`: `368/393` passes, with `576` ineligible pairs;
- `j=3`: `385/449` passes, with `520` ineligible pairs;
- `j=4`: `504/567` passes, with `402` ineligible pairs;
- `j=5`: `638/681` passes, with `288` ineligible pairs;
- `j=6`: `744/779` passes, with `190` ineligible pairs;
- `j=7`: `827/849` passes, with `120` ineligible pairs;
- `j=8`: `900/906` passes, with `63` ineligible pairs.

So the reversed tail vector is much closer to an interlacing sequence than the
raw position vector, but we still do not have full pairwise interlacing.

### Clean unrefined recurrence for the ordered family

There is a nice closed formula for the unrefined ordered polynomials after all,
but it is not a derivative-in-`t` recurrence like the unordered one.

If

- `O_j(z,t) = sum_{n>=0} O_{n,j}(t) z^n/n! = 1/(1-A_j(z,t))`,
- `A_j(z,t) = sum_{r=1}^{j-1} z^r/r! + t sum_{r>=j} z^r/r!`,

then

- `A'_j(z,t) = A_j(z,t) + 1 + (t-1) z^{j-1}/(j-1)!`,

so

- `O'_j(z,t) = (2 + (t-1) z^{j-1}/(j-1)!) O_j(z,t)^2 - O_j(z,t)`.

Thus the natural companion family is not `t O'_{n,j}(t)` but the square
convolution

- `C_{n,j}(t) = sum_{r=0}^n binom(n,r) O_{r,j}(t) O_{n-r,j}(t)`.

Coefficient extraction gives

- `O_{n+1,j}(t)
   = 2 C_{n,j}(t) - O_{n,j}(t)
     + (t-1) binom(n,j-1) C_{n-j+1,j}(t)`.

For `j=3` this specializes to

- `O_{n+1,3}(t)
   = 2 C_{n,3}(t) - O_{n,3}(t)
     + (t-1) binom(n,2) C_{n-2,3}(t)`.

So the ordered analogue of the unordered derivative recurrence is really a
Riccati / quadratic-convolution recurrence.

There is also a clean fixed-step recurrence for the companion itself.
Writing `d=j-1`, one gets by rearranging:

- `C_{n,j}(t)
   = (O_{n+1,j}(t)+O_{n,j}(t))/2
     - ((t-1)/2) binom(n,d) C_{n-d,j}(t)` for `n >= d`,
- and for `0 <= n < d`,
  `C_{n,j}(t) = (O_{n+1,j}(t)+O_{n,j}(t))/2`.

So `C_{n,j}` splits into `d` residue classes modulo `d`, each satisfying a
first-order recurrence with varying coefficients.
For `j=3`, this means the even and odd subsequences satisfy independent
first-order recurrences.

This is finite-length for the companion, but it still does not obviously turn
into a fixed finite-memory recurrence for the unrefined ordered sequence
`O_{n,j}(t)` alone.

For `j=3`, the cleanest matrix-like object is actually an infinite hierarchy.
If

- `D_n^{(m)}(t) = [z^n/n!] O_3(z,t)^m`,

so `D_n^{(1)} = O_{n,3}` and `D_n^{(2)} = C_{n,3}`, then

- `D_{n+1}^{(m)}(t)
   = -m D_n^{(m)}(t)
     + 2m D_n^{(m+1)}(t)
     + m(t-1) binom(n,2) D_{n-2}^{(m+1)}(t)`.

So there is a very clean linear recursion in the power index `m`, but it is an
infinite upper-shift hierarchy, not a closed finite-dimensional package on just
`O`, `C`, and a few of their predecessors.

## Ordered Ferrers rook-path model

Circled back to the rook side and wrote down the uniform ordered model that
packages together:

- a Ferrers board `mu`,
- a nesting weight `q`,
- the path-component big-block statistic `rch_j`,
- and an order on the path components.

For a standard rook placement `rho` on `mu`, build the reversed-column digraph
`D_mu(rho)` on `[N(mu)]` as in the Ferrers path model. Let

- `nest(rho)` be the usual number of nesting pairs,
- `kappa_mu(rho)` be the number of path components of `D_mu(rho)`,
- `rch_j^mu(rho)` be the number of path components with at least `j` vertices.

Since `D_mu(rho)` is a forest of directed paths on `N(mu)` vertices with
`|rho|` edges, one has

- `kappa_mu(rho) = N(mu) - |rho|`.

The ordered generating function is therefore

- `G_{mu,j}(u,q,t)
   = sum_rho kappa_mu(rho)! u^{kappa_mu(rho)}
       q^{nest(rho)} t^{rch_j^mu(rho)}`,

equivalently, it counts pairs `(rho, sigma)` where `sigma` is a total order on
the path components of `D_mu(rho)`.

This is exactly the right staircase extension:

- for `mu = delta_{n-1}` and `q=1`, it agrees with the ordered set-partition
  family `G_{n,j}(u,t)`,
- at `u=1`, this gives `O_{n,j}(t)`,
- at `q=0`, it gives an ordered non-nesting rook-path model on arbitrary
  Ferrers boards.

Added checker:

- `experiments/src/bin/ordered_ferrers_rook_paths.rs`

What it verifies:

- staircase sanity check: for `n <= 8` and `j = 2,3,4,5`, the `q=1`
  staircase specialization agrees exactly with the ordered set-partition
  big-block polynomials,
- endpoint real-rootedness:
  - `q=0`: all shapes `|mu| <= 16` passed for `j = 2,3,4` (`914/914` each),
  - `q=1`: all shapes `|mu| <= 16` passed for `j = 2,3,4` (`914/914` each).

Sample staircase rows:

- `j=2`:
  - `q=1`: `1`, `2+t`, `6+7t`, `24+45t+6t^2`, `120+311t+110t^2`, ...
  - `q=0`: `1`, `2+t`, `6+7t`, `24+45t+4t^2`, `120+311t+70t^2`, ...
- `j=3`:
  - `q=1`: `1`, `3`, `12+t`, `66+9t`, `450+91t`, ...
  - `q=0`: `1`, `3`, `12+t`, `64+9t`, `420+81t`, ...

This looks like a genuinely useful bridge object: it contains the ordered
set-partition family at one corner and an ordered non-nesting rook family at
another, while staying entirely inside the Ferrers-rook language.

## Useful commands

```bash
cargo check --offline -p experiments --bin nn_rook_grid
CARGO_TARGET_DIR=/tmp/rust-target cargo run --offline --quiet -p experiments --bin nn_rook_grid -- 20

cargo check --offline -p experiments --bin nn_rook_square
CARGO_TARGET_DIR=/tmp/rust-target cargo run --offline --quiet -p experiments --bin nn_rook_square

CARGO_TARGET_DIR=/tmp/rust-target cargo run --offline --quiet -p experiments --bin nn_rook_qdeform2
CARGO_TARGET_DIR=/tmp/rust-target cargo run --offline --quiet -p experiments --bin nn_rook_deform

cargo check --offline -p experiments --bin touchard_staircase
CARGO_TARGET_DIR=/tmp/rust-target cargo run --offline --quiet -p experiments --bin touchard_staircase

cargo check --offline -p experiments --bin ferrers_rook_paths
CARGO_TARGET_DIR=/tmp/rust-target cargo run --offline --quiet -p experiments --bin ferrers_rook_paths -- 18 4
CARGO_TARGET_DIR=/tmp/rust-target cargo run --offline --quiet -p experiments --bin ferrers_rook_paths -- 20 2

cargo check --offline -p experiments --bin ordered_big_blocks
CARGO_TARGET_DIR=/tmp/rust-target cargo run --offline --quiet -p experiments --bin ordered_big_blocks -- 24 6

cargo check --offline -p experiments --bin ordered_big_blocks_position
CARGO_TARGET_DIR=/tmp/rust-target cargo run --offline --quiet -p experiments --bin ordered_big_blocks_position -- 14 6

cargo check --offline -p experiments --bin ordered_ferrers_rook_paths
CARGO_TARGET_DIR=/tmp/rust-target cargo run --offline --quiet -p experiments --bin ordered_ferrers_rook_paths -- 16 4

latexmk -pdf -interaction=nonstopmode -halt-on-error documents/non-nesting-rooks.tex
```
