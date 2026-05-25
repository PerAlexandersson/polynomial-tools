# Bézout matrix for interlacing and real-rootedness

## Theory

Given polynomials f (degree d) and g (degree d-1), the **Bézout matrix** B(f,g) is the d x d symmetric matrix with entries

```
B[i,j] = coefficient of x^i y^j in (f(x)g(y) - f(y)g(x)) / (x - y)
```

There are two standard oriented versions of the criterion.

- If `deg(f) = deg(g) + 1`, then `g` interlaces `f` iff `B(f,g)` is
  positive semidefinite; in the coprime/strict case this is positive definite.
  See Kummer--Naldi--Plaumann, Theorem 2.13, citing Krein--Naimark, Section 2.2.
- If `deg(f) = deg(g)`, the same-degree alternation criterion uses the
  same Bezoutian with the orientation fixed by the argument order. Fisk,
  Section 9.21, Corollary 9.145 gives the positive definite/no-common-root
  form.

For real-rootedness: f is real-rooted iff B(f, f') is positive semi-definite (semi-definite because repeated roots make B singular).

When f and g share roots, B(f,g) is singular. Divide out gcd(f,g) first,
verify that the gcd is real-rooted, then check strict interlacing of the
reduced polynomials.  The gcd check matters: positive semidefiniteness of the
Bezoutian alone does not certify that a common factor has only real roots.

## Mathematica code

```mathematica
(* --- Bézout matrix -------------------------------------------------- *)

BezoutMatrix[f_, g_, t_] := Module[
  {df = Exponent[f, t], dg = Exponent[g, t], h, b},
  (* Require deg f = deg g + 1; swap externally if needed *)
  h = Cancel[(f /. t -> x) (g /. t -> y) - (f /. t -> y) (g /. t -> x)] / (x - y) // Expand;
  b = Table[
    Coefficient[Coefficient[h, x, i], y, j],
    {i, 0, df - 1}, {j, 0, df - 1}
  ];
  b
]

(* --- Strict interlacing via Bézout matrix --------------------------- *)

CheckInterlacingBezout[f_, g_, t_] := Module[
  {df = Exponent[f, t], dg = Exponent[g, t], ff, gg, b},
  (* Auto-swap so ff has the higher degree *)
  If[df == dg + 1,
    {ff, gg} = {f, g},
    If[dg == df + 1,
      {ff, gg} = {g, f},
      Return[Indeterminate] (* degree diff is not 1 *)
    ]
  ];
  (* Align leading coefficients to be positive *)
  If[Coefficient[ff, t, Exponent[ff, t]] < 0, ff = -ff];
  If[Coefficient[gg, t, Exponent[gg, t]] < 0, gg = -gg];
  b = BezoutMatrix[ff, gg, t];
  PositiveDefiniteMatrixQ[b]
]

(* --- Weak interlacing (shared roots allowed) ------------------------ *)

CheckWeakInterlacingBezout[f_, g_, t_] := Module[
  {d, fred, gred},
  d = PolynomialGCD[f, g, t];
  (* Production code should also verify that d is real-rooted.  The Rust
     implementation does this before checking the reduced pair. *)
  fred = Cancel[f / d];
  gred = Cancel[g / d];
  (* If both reduce to constants, all roots shared: trivially interlacing *)
  If[Exponent[fred, t] == 0 && Exponent[gred, t] == 0,
    Return[True]
  ];
  CheckInterlacingBezout[fred, gred, t]
]

(* --- Real-rootedness via Bézout matrix ------------------------------ *)

RealRootedBezoutQ[f_, t_] := Module[
  {d = Exponent[f, t], ff, fp, b},
  If[d <= 1, Return[True]];
  ff = f;
  fp = D[f, t];
  (* Align leading coefficients *)
  If[Coefficient[ff, t, d] Coefficient[fp, t, d - 1] < 0, fp = -fp];
  b = BezoutMatrix[ff, fp, t];
  PositiveSemidefiniteMatrixQ[b]
]
```

## Examples

```mathematica
(* Eulerian polynomial A_4(t) = 1 + 11t + 11t^2 + t^3 *)
f = 1 + 11 t + 11 t^2 + t^3;
RealRootedBezoutQ[f, t]
(* True *)

(* Interlacing: (t-1)(t-3) and (t-2) *)
f = (t - 1)(t - 3) // Expand;
g = (t - 2);
CheckInterlacingBezout[f, g, t]
(* True *)

(* Weak interlacing with shared root *)
f = (t - 1)^2 (t - 3) // Expand;
g = (t - 1)(t - 2) // Expand;
CheckWeakInterlacingBezout[f, g, t]
(* True *)

(* Inspect the Bézout matrix directly *)
BezoutMatrix[(t - 1)(t - 3)(t - 5) // Expand, (t - 2)(t - 4) // Expand, t] // MatrixForm
```

## Notes

- Mathematica's `PositiveDefiniteMatrixQ` / `PositiveSemidefiniteMatrixQ` work with exact rationals, so this is fully rigorous.
- The Bézout matrix approach is much faster than root isolation for checking interlacing — it avoids computing roots entirely.
- For numerical (approximate) coefficients, use `PositiveDefiniteMatrixQ` with a tolerance, or convert to exact rationals first via `Rationalize[..., 0]`.
- References:
  - S. Fisk, *Polynomials, roots, and interlacing*, arXiv:math/0612833,
    Section 9.21, Corollary 9.145.
  - M. Kummer, S. Naldi, and D. Plaumann, *Spectrahedral representations of
    plane hyperbolic curves*, arXiv:1807.10901, Theorem 2.13.
  - M. G. Krein and M. A. Naimark, *The method of symmetric and Hermitian
    forms in the theory of the separation of the roots of algebraic
    equations*, Linear and Multilinear Algebra 10 (1981), Section 2.2.
  - MathOverflow discussion of the common-factor caveat:
    https://mathoverflow.net/questions/403708/bezout-matrices-and-interlacing-roots
