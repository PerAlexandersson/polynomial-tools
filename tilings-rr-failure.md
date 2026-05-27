# Real-rootedness failures for the `bb` / `ff` family

We consider

\[
bb(x,t) = x + x^s \frac{(1 + t x^r ((x^s - 1)/(x - 1)))^d - 1}{(x^s - 1)/(x - 1)},
\qquad
ff(x,t) = \frac{1}{1 - bb(x,t)}.
\]

Write

\[
P_{d,r,s,n}(t) = [x^n] ff(x,t).
\]

For fixed `(d,r,s)`, the sequence `P_{d,r,s,n}(t)` in `n` satisfies a linear recurrence, and the Rust experiment in
[experiments/src/bin/bb_rr_table.rs](/home/paxinum/Dropbox/AI-projects/rust/experiments/src/bin/bb_rr_table.rs)
searches for the first `n` where real-rootedness fails.
The completed table below was also cross-checked with Wolfram `CountRoots`.

## Search box

- `d = 1, ..., 10`
- `r = 2, ..., 10`
- `s = 1, ..., 200`
- `n = 1, ..., 200`

In the table below, an entry `a->b` means:

- `a` is the smallest `s` for which a failure was found in the search box, and
- `b` is then the smallest `n <= 200` such that `P_{d,r,a,b}(t)` is not real-rooted.

Also:

- `NF` means no failure was found for that `(d,r)` in the full search box.

## Completed table

| d\\r | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | NF | NF | NF | NF | NF | NF | NF | NF | NF |
| 2 | 4->98 | 5->190 | 8->164 | 11->192 | 15->195 | 22->189 | 29->198 | NF | NF |
| 3 | 5->111 | 7->161 | 10->175 | 13->196 | 19->186 | 25->189 | 34->188 | 57->185 | NF |
| 4 | 3->88 | 9->164 | 12->190 | 18->182 | 24->185 | 30->200 | 45->183 | 54->193 | NF |
| 5 | 3->37 | 11->174 | 15->199 | 22->160 | 27->196 | 34->158 | 39->182 | NF | NF |
| 6 | 3->28 | 4->66 | 7->147 | 11->195 | 15->200 | 19->185 | 31->166 | 35->188 | NF |
| 7 | 3->17 | 4->47 | 5->132 | 8->125 | 10->196 | 13->175 | 16->200 | 20->171 | 22->190 |
| 8 | 3->17 | 4->28 | 5->74 | 6->184 | 9->160 | 11->185 | 14->153 | 16->166 | 18->186 |
| 9 | 3->17 | 4->28 | 5->78 | 6->102 | 8->172 | 10->144 | 12->161 | 14->183 | 16->197 |
| 10 | 3->17 | 4->28 | 5->45 | 6->107 | 8->127 | 9->160 | 11->172 | 13->189 | 16->190 |

## Immediate observations

1. The entire `d=1` row looks safe in this search box.

2. The failure region is not monotone in `r`.
   For example, for `d=5` there are failures through `r=8`, but `r=9,10` are `NF` in the current box.

3. For moderate `d` (say `d=3,4,5,6`) the first failures at larger `r` often occur very near the cutoff `n=200`.

4. For larger `d`, failures can occur much earlier.
   For example:
   - `(d,r) = (7,2)` has first witness `3->17`,
   - `(d,r) = (10,4)` has first witness `5->45`,
   - `(d,r) = (9,5)` has first witness `6->102`.

5. The minimal `s` tends to drop again for larger `d`.
   Compare the rows `d=3,4,5` with `d=8,9,10`.

## A few especially striking cells

- `(7,2)`, `(8,2)`, `(9,2)`, and `(10,2)` all give `3->17`.
- `(5,2)` gives `3->37`.
- `(8,4)` gives `5->74`.
- `(9,4)` gives `5->78`.
- `(10,4)` gives `5->45`.

These suggest that for larger `d`, the first non-real-rooted examples may occur at quite small `s` and `n`.
