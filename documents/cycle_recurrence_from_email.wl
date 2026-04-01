(* Mathematica script encoding the cycle recurrence from the email.
   It uses the formula verbatim, with the Lucas-type sequence
   1, 3, 4, 7, 11, 18, ... as stated in the message. *)

ClearAll[t, c, coeffList, generateData, forcingTerm, emailLucas];

(* Lucas-type sequence from the email, not Mathematica's built-in LucasL. *)
emailLucas[1] = 1;
emailLucas[2] = 3;
emailLucas[n_Integer?Positive] /; n >= 3 :=
  emailLucas[n] = emailLucas[n - 1] + emailLucas[n - 2];

(* Initial values. *)
c[3] = 1 + 3 t + 2 t^2 + t^3;
c[4] = 1 + 4 t + 7 t^2 + 3 t^3;
c[5] = 1 + 5 t + 14 t^2 + 10 t^3 + t^4;

(* Inhomogeneous term from the email. *)
forcingTerm[n_Integer] /; EvenQ[n] && n >= 6 :=
  t^3 + (2^(n/2 - 2) + 1) t^(n/2) (1 - t);

forcingTerm[n_Integer] /; OddQ[n] && n >= 7 :=
  t^3 + (emailLucas[n - 4] - 1 + 2 t) t^((n - 1)/2) (1 - t);

(* Recurrence:
   c_n - c_(n-1) - 3 c_(n-2) + 2 c_(n-3) = forcingTerm[n]. *)
c[n_Integer] /; n >= 6 :=
  c[n] = Expand[c[n - 1] + 3 c[n - 2] - 2 c[n - 3] + forcingTerm[n]];

coeffList[poly_] := CoefficientList[Expand[poly], t];
coeffList[n_Integer] := coeffList[c[n]];

generateData[nMax_Integer?Positive, nMin_Integer:3] :=
  Table[n -> coeffList[n], {n, nMin, nMax}];

(* Default output when run as a script. *)
Print["Cycle recurrence from email"];
Print["Lucas-type sequence used: ", Table[emailLucas[k], {k, 1, 10}]];
Print["Coefficient data for c_n, n = 3..15:"];
Do[
  Print["c_", n, " = ", coeffList[n]],
  {n, 3, 15}
];

