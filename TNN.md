# Totally Nonnegative Matrices, Planar Networks, and Reconstruction

This note is for the implementation problem:

- input: a finite prefix of a monic polynomial sequence
  `P_0(t), P_1(t), ..., P_N(t)` with `deg(P_n) = n`,
- form the coefficient matrix
  `R = (r_{n,k})_{0 <= k <= n <= N}` via
  `P_n(t) = sum_{k=0}^n r_{n,k} t^k`,
- if `R` is lower-unitriangular and totally nonnegative (TNN), reconstruct a
  planar network with nonnegative edge weights whose path matrix is `R`.

The short version is:

- LGV gives `planar network => TNN`.
- For lower-unitriangular matrices there is also a converse: every TNN matrix
  comes from a canonical triangular planar network.
- The canonical network can be reconstructed recursively by a Whitney-type
  reduction on the first column.
- This is the right theorem to implement in Rust.

There is one important caveat up front: the statement

> "if the coefficient matrix of the polynomials is TNN, then the row
> polynomials themselves are real-rooted"

is false in general.

For example,

```text
R = [ [1,0,0],
      [1,1,0],
      [1,1,1] ]
```

is lower-unitriangular and TNN, but its third row polynomial is

```text
1 + t + t^2
```

which is not real-rooted.

So there are really two different mechanisms in the literature:

1. `lower-unitriangular TNN => canonical network => resolvability => chain polynomials are real-rooted`
2. `Toeplitz matrix of a coefficient sequence is TNN => that one polynomial is real-rooted`

The first mechanism is the one that gives the network reconstruction from a
lower-unitriangular matrix. The second is the Pólya frequency / Toeplitz route.
Some recent papers build a single "unified" network that realizes both the
triangle and the Toeplitz matrices of its rows, but that is extra structure,
not automatic for an arbitrary TNN coefficient matrix.

There is also a third, especially useful, certificate due to Brenti:

- if the row polynomials come from a planar weakly `y`-invariant digraph in the
  sense of Brenti, then every row is automatically a Pólya frequency sequence,
  so every row polynomial is real-rooted with nonpositive zeros.

This is stronger than merely reconstructing a canonical TNN network from a
matrix.

## 1. The matrix attached to a monic polynomial sequence

Write

```math
P_n(t) = \sum_{k=0}^n r_{n,k} t^k,
\qquad r_{n,n} = 1.
```

If we index columns by increasing degree, then the coefficient matrix

```math
R = (r_{n,k})_{n,k >= 0}
```

is lower-unitriangular:

- `r_{n,k} = 0` for `k > n`,
- `r_{n,n} = 1`.

So the natural input class for a monic polynomial sequence is exactly the class
treated by Branden--Saud Maia Leite.

For a finite computation we only need the principal truncation

```math
R^{(N)} = (r_{n,k})_{0 <= k <= n <= N}.
```

If the user gives `P_1, P_2, ...`, the right convention is to prepend
`P_0(t) = 1`.

## 2. LGV gives path matrix => TNN

Let `Gamma_N` be the triangular directed graph with vertices

```math
(i,j), \qquad 0 <= j <= i <= N,
```

with edges

- horizontal: `(i,j) -> (i,j+1)`,
- vertical: `(i+1,j) -> (i,j)`.

Attach weight `1` to each horizontal edge, and attach a nonnegative weight
`lambda_{i,j}` to each vertical edge.

Let `r_{n,k}(Gamma_N, lambda)` be the total weight of all directed paths from
`(n,0)` to `(k,k)`. Then the path matrix

```math
R(Gamma_N, lambda) = (r_{n,k}(Gamma_N, lambda))
```

is TNN by the Lindstrom--Gessel--Viennot lemma.

This is the easy direction:

```text
planar network with nonnegative weights => TNN matrix
```

and it is the reason planar networks are such a convenient certificate.

## 3. The lower-unitriangular converse

For the implementation problem, the crucial result is the converse in the
lower-unitriangular setting. In Branden--Saud Maia Leite this is Theorem 2.1.

Let `R = (r_{n,k})_{0 <= k <= n <= N}` be lower-unitriangular. Then:

- `R` is TNN if and only if there exists a nonnegative triangular array
  `lambda = (lambda_{n,k})_{0 <= k <= n < N}` such that
  `R = R(Gamma_N, lambda)`.
- We may require the zero-propagation condition

```math
lambda_{n,k} = 0 \Longrightarrow lambda_{n+1,k} = 0.
```

- Under that condition, `lambda` is unique.

This gives a canonical network attached to `R`, not just some network.

That uniqueness is exactly what makes the theorem implementable: we are not
solving an underdetermined inverse problem, we are reconstructing the canonical
`lambda` array.

## 4. Whitney reduction and the recursive reconstruction

The reconstruction is recursive on the size of the matrix.

Let

```math
R = (r_{n,k})_{0 <= k <= n <= N}.
```

Assume `R` is lower-unitriangular and TNN.

Let `m+1` be the smallest index such that `r_{m+1,0} = 0`. If there is no such
index, set `m = N`.

Define the first-column weights

```math
mu_n =
\begin{cases}
r_{n+1,0} / r_{n,0}, & n <= m, \\
0, & n > m.
\end{cases}
```

Then define the reduced matrix `R_tilde = (r_tilde_{n,k})_{0 <= k <= n <= N-1}`
by

```math
r_tilde_{n,k} = r_{n+1,k+1} - mu_n r_{n,k+1}.
```

Whitney reduction says:

- `R` is TNN if and only if all entries below the first zero in the first column
  also vanish, and `R_tilde` is TNN.
- If `gamma` is the canonical network array for `R_tilde`, then the canonical
  array for `R` is obtained by adding the first column `mu` to the left:

```text
lambda_{n,0} = mu_n
lambda_{n,k} = gamma_{n-1,k-1}   for 1 <= k <= n
```

In words:

- first recover the whole first column of vertical weights from the first column
  of the matrix,
- peel that column off,
- recurse on the smaller lower-unitriangular matrix obtained by shifting rows
  and columns by one.

This is the algorithm we want in code.

### Why this formula makes combinatorial sense

A path from `(n,0)` to `(k,k)` can do one of two things:

- stay on the first column for a while, accumulating the first-column vertical
  weights, then move right once into the shifted subnetwork,
- or, in the case `k = 0`, stay entirely in the first column.

So the first column separates cleanly from the rest of the network. The reduced
matrix exactly records what remains after removing the contribution from that
first column.

## 5. Equivalent operator picture: resolvability

Branden--Saud Maia Leite package the same structure in a second way.

A lower-unitriangular matrix `R` is called *resolvable* if there exist:

- nonnegative numbers `lambda_{n,k}`,
- monic polynomials `R_{n,k}(t)` with `t^k | R_{n,k}(t)`,

such that

```math
R_{n,0}(t) = R_n(t) = \sum_{k=0}^n r_{n,k} t^k,
\qquad
R_{n,n}(t) = t^n,
```

and

```math
R_{n+1,k}(t) = R_{n+1,k+1}(t) + lambda_{n,k} R_{n,k}(t).
```

They prove the equivalence

```text
R is resolvable <=> R is lower-unitriangular TNN
```

This is Theorem 2.6 in their paper.

and also the diagonal-operator factorization

```math
R_n(t) = (t + alpha_1)(t + alpha_2) \cdots (t + alpha_n) 1,
```

where each `alpha_i` is a diagonal operator:

```math
alpha_i(t^k) = alpha_{i,k} t^k, \qquad alpha_{i,k} >= 0.
```

The two parametrizations are equivalent:

```math
lambda_{n,k} = alpha_{n+1-k,k}.
```

This operator viewpoint is very useful conceptually, but for implementation the
Whitney reduction above is the cleaner reconstruction algorithm.

## 6. What real-rootedness result is actually true

The raw row polynomials

```math
R_n(t) = \sum_{k=0}^n r_{n,k} t^k
```

need not be real-rooted for an arbitrary lower-unitriangular TNN matrix.

What Branden--Saud Maia Leite prove is that the associated *chain polynomials*

```math
p_0(t) = 1,
\qquad
p_n(t) = t \sum_{k=0}^{n-1} r_{n,k} p_k(t)
```

are real-rooted with zeros in `[-1,0]`, and in fact form an interlacing
sequence.

That is their Theorem 3.7.

So if we implement `matrix -> canonical network`, we should keep the
real-rootedness claims phrased in terms of:

- chain polynomials, or
- Toeplitz/PF matrices of individual rows,

not in terms of arbitrary raw row polynomials.

## 7. Where the Toeplitz/PF route enters

For a single polynomial

```math
f(t) = \sum_{k=0}^d a_k t^k
```

with `a_k >= 0`, real-rootedness is equivalent to total nonnegativity of the
Toeplitz matrix `(a_{i-j})`. This is the Aissen--Schoenberg--Whitney theorem on
Pólya frequency sequences; in Branden--Saud Maia Leite this appears as their
Theorem 4.1.

This is different from saying that the coefficient *triangle*
`(r_{n,k})` is TNN.

The recent Chen--Fu--Ruan framework is interesting exactly because it builds a
single planar network that can realize:

- the lower-triangular matrix itself,
- its reversal,
- and Toeplitz matrices attached to the rows,

by choosing different source/sink sets in one common graph.

That is why their method can prove real-rootedness of actual row generating
polynomials for certain structured combinatorial triangles. But this is extra
structure; it is not a generic consequence of lower-unitriangular TNN alone.

For this project, that suggests two future implementation paths:

1. generic path: reconstruct the canonical `Gamma_N` network from an arbitrary
   lower-unitriangular TNN matrix;
2. structured-family path: if a matrix comes from a recurrence / production
   matrix of the Chen--Fu--Ruan type, build the unified network directly instead
   of reconstructing it from each finite prefix.

## 8. Brenti's planar weakly `y`-invariant digraph theorem

Brenti uses the term "totally positive" in the older sense:

```text
every minor has nonnegative determinant
```

so in our terminology this is total nonnegativity.

His corollary says that if `D` is:

- a locally finite nonnegative digraph on `N x N`,
- planar,
- weakly `y`-invariant, meaning the outgoing edges from `(m,k)` do not depend on
  `m`,

and

```math
M_{n,k} = P_D((0,0),(n,k)),
```

then:

- `M` is TNN,
- every row `(M_{n,k})_{k >= 0}` is a PF sequence,
- hence every row polynomial `sum_k M_{n,k} t^k` has only real nonpositive
  zeros.

This is exactly the kind of certificate one wants when the end goal is not only
to prove TNN of a matrix, but to prove real-rootedness of the actual row
polynomials.

That motivates the following definition.

### Brenti sequence

Call a polynomial sequence `P_1(t), P_2(t), ...` a *Brenti sequence* if there
exists a planar weakly `y`-invariant nonnegative digraph `D` on `N x N` and row
indices `r_1, r_2, ...` such that

```math
P_i(t) = \sum_{k >= 0} M_{r_i,k} t^k,
\qquad
M_{n,k} = P_D((0,0),(n,k)).
```

Exhibiting such a digraph is then a proof of real-rootedness for the sequence.

### What we can implement cleanly

The full class of Brenti digraphs is broad. For code, the most practical first
step is a restricted but theorem-compatible model:

- vertices `(x,y)` with `x,y in N`,
- edges always go from column `x` to column `x+1`,
- the outgoing transitions depend only on the height `y`,
- planarity is guaranteed by an explicit monotonicity condition on target
  heights.

This "strip digraph" model still covers many natural lattice-path and transfer
matrix examples, including Pascal-type networks, and is straightforward to
compute with exactly.

## 9. Rust implementation plan

The repo already has the exact arithmetic and TNN certification pieces we need:

- `polytool/src/linalg.rs` has exact rational linear algebra and
  `check_tnn_neville_bigint`.
- `num-rational` is already in the workspace, so returning `BigRational`
  weights is natural.

### Suggested data type

```rust
pub struct CanonicalPlanarNetwork<T> {
    /// lambda[n][k] for 0 <= k <= n < num_rows - 1
    pub lambda: Vec<Vec<T>>,
}
```

For exact reconstruction from integer input, use
`T = num_rational::Ratio<num_bigint::BigInt>`.

### Suggested API

```rust
pub fn reconstruct_canonical_tnn_network(
    rows: &[Vec<BigInt>],
) -> Result<CanonicalPlanarNetwork<BigRational>, NetworkError>;
```

and then a convenience wrapper

```rust
pub fn reconstruct_from_monic_polynomials(
    polys: &[Polynomial<BigInt>],
) -> Result<CanonicalPlanarNetwork<BigRational>, NetworkError>;
```

where the wrapper:

- prepends `P_0 = 1` if the caller supplied `P_1, P_2, ...`,
- builds the lower-unitriangular coefficient matrix in increasing powers,
- calls the matrix routine.

### Preconditions to validate

Before reconstructing:

1. matrix is triangular: `r_{n,k} = 0` for `k > n`
2. matrix is unitriangular: `r_{n,n} = 1`
3. entries are exact and nonnegative
4. TNN holds, using `check_tnn_neville_bigint`

If any of these fail, return a descriptive error.

### Reconstruction pseudocode

```text
recover_lambda(R):
    # R has rows 0..N, row n has entries 0..n
    if N == 0:
        return []

    find z = smallest i in {1,...,N} with R[i][0] == 0
    if no such i exists:
        z = N + 1

    verify R[i][0] == 0 for all i >= z

    for n = 0..N-1:
        if n + 1 <= z:
            mu[n] = R[n+1][0] / R[n][0]
        else:
            mu[n] = 0

    build R_tilde of size N:
        R_tilde[n][k] = R[n+1][k+1] - mu[n] * R[n][k+1]
        for 0 <= k <= n <= N-1

    recursively compute gamma = recover_lambda(R_tilde)

    build lambda:
        lambda[n][0] = mu[n]
        lambda[n][k] = gamma[n-1][k-1] for 1 <= k <= n

    return lambda
```

### Remarks for the implementation

- Use exact rationals internally even if the input matrix is integral. The
  canonical `lambda` weights need not remain integers.
- For a finite prefix `R^{(N)}`, this reconstructs the canonical network for
  that truncation. For an infinite family, stabilization across `N` is a
  separate compatibility question.
- Neville elimination is excellent for checking TNN, but the reconstruction
  itself should follow Whitney reduction, not the Neville row-subtraction path.

## 10. Sanity-check example: Pascal

For the Pascal matrix

```math
r_{n,k} = \binom{n}{k},
```

the reconstruction gives

```math
lambda_{n,k} = 1
```

for all `0 <= k <= n`.

So the canonical triangular network for Pascal has weight `1` on every vertical
edge and weight `1` on every horizontal edge. This is a good first regression
test for the Rust implementation.

## 11. Recommended implementation order

1. Add a new module in `polytool` for canonical network reconstruction.
2. Start with matrix input over `BigInt`, convert to `BigRational` internally.
3. Reuse `check_tnn_neville_bigint` as the validator.
4. Implement the recursive Whitney reduction exactly as above.
5. Add a path-matrix evaluator for `Gamma_N(lambda)` so the tests can verify
   `reconstruct -> evaluate -> original matrix`.
6. Add tests:
   - identity matrix
   - Pascal matrix
   - a matrix with a zero tail in the first column
   - a non-TNN input that should be rejected
7. Only after that add a polynomial-sequence wrapper.

At that point add a second module for Brenti certificates:

1. represent a weakly `y`-invariant strip digraph,
2. compute rows `M_{n,k}`,
3. build row polynomials,
4. verify that a given polynomial sequence appears among those rows,
5. use Brenti's theorem as the real-rootedness certificate.

## References

- Bernt Lindstrom, "On the Vector Representations of Induced Matroids"
  (1973). DOI: https://doi.org/10.1112/blms/5.1.85
- Ira Gessel and Gerard Viennot, "Binomial determinants, paths, and hook
  length formulae" (1985). DOI: https://doi.org/10.1016/0001-8708(85)90121-5
- Petter Branden and Leonardo Saud Maia Leite, "Totally nonnegative matrices,
  chain enumeration and zeros of polynomials" (arXiv:2412.06595).
  https://arxiv.org/abs/2412.06595
- Francesco Brenti, "Combinatorics and total positivity" (1995).
  DOI: https://doi.org/10.1016/0097-3165(95)90000-4
- Samuel Karlin, *Total Positivity, Vol. I* (1968).
- Shaun M. Fallat and Charles R. Johnson, *Totally Nonnegative Matrices*
  (2011). DOI: https://doi.org/10.1515/9781400839018
- Xi Chen, Lang Fu, and Jiajie Ruan, "A unified planar network approach to
  total positivity of combinatorial matrices and real-rootedness of
  polynomials" (arXiv:2512.08369). https://arxiv.org/abs/2512.08369
