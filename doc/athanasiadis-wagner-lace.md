# Athanasiadis--Wagner interlacing matrices

Reference: C.~A. Athanasiadis and V.~Wagner, *Veronese sections and
interlacing matrices of formal power series*.

For a `p x q` matrix `A` of formal power series

```text
A_ij(x) = sum_n a_ij(n) x^n,
```

Definition 3.5 constructs the infinite matrix `Lace(A)` indexed by integers.
If

```text
u = p u' + i,    v = q v' + j,
0 <= i < p,     0 <= j < q,
```

then

```text
Lace(A)_{u,v} = a_ij(v' - u').
```

For a column vector of series, Athanasiadis--Wagner call the vector fully
interlacing when this infinite `Lace(A)` matrix is totally nonnegative.  This
implies ordinary pairwise interlacing of the entries, but the converse fails.

## Pairwise but not fully interlacing

Their Example 3.4 takes

```text
P(x) = t + x,
Q(x) = (b+x)(d+x),
R(x) = (a+x)(c+x).
```

The triple is pairwise interlacing when `a <= b <= t <= c <= d`.  However,
the finite `3 x 3` truncation

```text
[ t   1     0 ]
[ bd  b+d   1 ]
[ ac  a+c   1 ]
```

has determinant

```text
t(b+d-a-c) - (bd-ac),
```

which can be negative.  For example, `a=1`, `b=2`, `t=2`, `c=3`, `d=4`
gives determinant `-1`.

The Rust module `interlacing_matrix` implements finite truncations of
`Lace(A)` and exact finite total-nonnegativity checks.  These truncations are
useful for experiments; proving full interlacing of the infinite matrix still
requires a theorem or finite criterion that controls all minors.
