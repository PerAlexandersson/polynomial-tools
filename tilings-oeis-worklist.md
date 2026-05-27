# OEIS worklist for partial tilings

This file collects the sequences from `main.tex` that should either be
submitted to OEIS as new entries or used to update existing entries with the
tiling interpretation from the paper.

Status notes:

- Terms below were generated from the paper formulas and the Rust experiment
  `experiments/src/bin/tiling_oeis.rs`.
- Entries marked `New` in the paper should still get one final exact OEIS
  duplicate search before submission.
- For Ferrers shapes of the form `(h,...)`, the first column height is ignored
  by the anchor-word inequalities once it is large enough; the listed terms use
  the smallest representative such as `(4,2,2,1)` for `(h,2,2,1)`.

## Submit as new entries

### Hook / L-polyomino anchor counts

General form:

```text
H_d(n) = [x^n] 1/(1 - x(1+x)^d), n >= 0.
```

These are the counts `|anchor(d,n)|` in Table `tab:LCountTable`.

#### `H_6(n)`

Action: submit new.

Formula:

```text
OGF = 1/(1 - x(1+x)^6)
```

Terms, `n = 0,1,2,...`:

```text
1, 1, 7, 28, 105, 413, 1624, 6349, 24851, 97315, 380989, 1491567, 5839638, 22862658, 89508951, 350434385, 1371977475, 5371396171, 21029425081, 82331801783, 322335276473, 1261967164192, 4940697593195, 19343207491818, 75730131016730
```

#### `H_7(n)`

Action: submit new.

Formula:

```text
OGF = 1/(1 - x(1+x)^7)
```

Terms, `n = 0,1,2,...`:

```text
1, 1, 8, 36, 148, 638, 2766, 11908, 51284, 221049, 952613, 4104980, 17689720, 76231176, 328505052, 1415636084, 6100444144, 26288828104, 113287234449, 488192085677, 2103780839028, 9065886015180, 39067895161612, 168356455227146, 725503534206186
```

### Ferrers tiling counts `P_n(1)`

#### `P_{(3,1,1),3,n}(1)`

Action: submit new. This is the sequence displayed near the example computing
`B_{311,3}(x,t)`.

Formula:

```text
OGF = 1/(1 - x - 3*x^3 - 3*x^4 - 4*x^5 - 2*x^6 - x^7)
```

Terms, `n = 0,1,2,...`:

```text
1, 1, 1, 4, 10, 20, 41, 90, 199, 431, 928, 2009, 4361, 9455, 20478, 44361, 96132, 208321, 451389, 978051, 2119265, 4592123, 9950352, 21560630, 46718097, 101230073
```

#### `P_{(3,2,1),3,n}(1)`

Action: submit new.

Formula:

```text
OGF = 1/(1 - x - 3*x^3 - x^4 - 3*x^5 - x^7)
```

Terms, `n = 0,1,2,...`:

```text
1, 1, 1, 4, 8, 15, 31, 63, 129, 262, 531, 1082, 2201, 4474, 9100, 18507, 37638, 76546, 155671, 316593, 643864, 1309437, 2663032, 5415868, 11014368, 22400164
```

#### `P_{(h,2,2,1),4,n}(1)`

Action: submit new. Terms use representative `mu=(4,2,2,1)`.

Formula:

```text
OGF = 1/(1 - x - 4*x^4 - 3*x^5 - 3*x^6 - 6*x^7 - 2*x^8 - 2*x^9 - 4*x^10 - x^13)
```

Terms, `n = 0,1,2,...`:

```text
1, 1, 1, 1, 5, 12, 22, 38, 72, 148, 301, 593, 1149, 2242, 4423, 8743, 17218, 33815, 66436, 130699, 257247, 506154, 995529, 1957987, 3851477, 7576734
```

#### `P_{(3,3,2),3,n}(1)`

Action: submit new.

Formula:

```text
OGF = 1/(1 - x - 3*x^3 - x^5)
```

Terms, `n = 0,1,2,...`:

```text
1, 1, 1, 4, 7, 11, 24, 46, 83, 162, 311, 584, 1116, 2132, 4046, 7705, 14685, 27939, 53186, 101287, 192809, 367052, 698852, 1330465, 2532908, 4822273
```

### Dense tiling counts `Q_n(1)`

#### `Q_{(3,1,1),3,n}(1)`

Action: submit new.

Formula:

```text
OGF = 1/(1 - 3*x^3 - 3*x^4 - 4*x^5 - 2*x^6 - x^7)
```

Terms, `n = 0,1,2,...`:

```text
1, 0, 0, 3, 3, 4, 11, 19, 33, 63, 115, 211, 390, 715, 1315, 2422, 4452, 8187, 15062, 27702, 50950, 93714, 172366, 317030, 583111, 1072506
```

#### `Q_{(h,2,2,1),4,n}(1)`

Action: submit new. Terms use representative `mu=(4,2,2,1)`.

Formula:

```text
OGF = 1/(1 - 4*x^4 - 3*x^5 - 3*x^6 - 6*x^7 - 2*x^8 - 2*x^9 - 4*x^10 - x^13)
```

Terms, `n = 0,1,2,...`:

```text
1, 0, 0, 0, 4, 3, 3, 6, 18, 26, 37, 66, 125, 209, 344, 591, 1025, 1747, 2975, 5086, 8695, 14850, 25394, 43403, 74119, 126626
```

### Cylindric L-shape tiling counts

General form for `mu=(h,1)` and `d=h`:

```text
C_d(x) = 1 + x*(d*x*(1+x)^(d-1) + (1+x)^d)/(1 - x*(1+x)^d).
```

The paper table lists `c_d(n)` for `n = 1,2,3,...`, i.e. the constant term
`c_d(0)=1` is omitted from the displayed sequence.

#### `c_3(n)`

Action: submit new.

Terms, `n = 1,2,3,...`:

```text
1, 7, 19, 47, 126, 331, 869, 2287, 6013, 15812, 41581, 109343, 287535, 756119, 1988334, 5228639, 13749533, 36156571, 95079421, 250026372, 657483881, 1728957831, 4546568011, 11955919519, 31439980926
```

#### `c_4(n)`

Action: submit new.

Terms, `n = 1,2,3,...`:

```text
1, 9, 31, 89, 276, 855, 2626, 8089, 24916, 76724, 236281, 727655, 2240876, 6900994, 21252276, 65448409, 201554636, 620706780, 1911525876, 5886726724, 18128737861, 55829181769, 171931303841, 529478890855, 1630580875001
```

#### `c_5(n)`

Action: submit new.

Terms, `n = 1,2,3,...`:

```text
1, 11, 46, 151, 526, 1862, 6518, 22839, 80110, 280886, 984842, 3453214, 12108097, 42454836, 148860421, 521952807, 1830135154, 6417045458, 22500235785, 78893099626, 276624709127, 969935648264, 3400907910453, 11924682463686, 41811791323776
```

#### `c_6(n)`

Action: submit new.

Terms, `n = 1,2,3,...`:

```text
1, 13, 64, 237, 911, 3604, 14113, 55181, 216091, 846103, 3312387, 12968220, 50771904, 198775653, 778221389, 3046796781, 11928441715, 46700756563, 182837020787, 715820868167, 2802493231336, 10971974511203, 42956116132980, 168176467326476, 658423682322911
```

#### `c_7(n)`

Action: submit new.

Terms, `n = 1,2,3,...`:

```text
1, 15, 85, 351, 1471, 6399, 27637, 118911, 512410, 2208595, 9517289, 41012271, 176736665, 761617795, 3282058705, 14143474111, 60948907145, 262648952136, 1131841009133, 4877476515311, 21018656212093, 90576327070959, 390323289449083, 1682031886574271, 7248430579431021
```

#### `c_8(n)`

Action: submit new.

Terms, `n = 1,2,3,...`:

```text
1, 17, 109, 497, 2251, 10637, 50107, 234737, 1100872, 5166027, 24237082, 113705261, 533453428, 2502727979, 11741634154, 55086294257, 258439478200, 1212478657352, 5688389395261, 26687294254827, 125204452428682, 587401431332282, 2755816064865757, 12929015461655213, 60656965802389576
```

## Update existing OEIS entries with tiling interpretations

### Hook / L-polyomino family

These all belong to the same family

```text
H_d(n) = [x^n] 1/(1 - x(1+x)^d), n >= 0.
```

- `A000045`: `H_1(n)`, Fibonacci numbers. Add interpretation as anchor words
  `|anchor(1,n)|`, equivalently tilings by the `d=1` L-family specialization.
- `A002478`: `H_2(n)`. Add interpretation as tilings of the `3 x n` board by
  L-triominoes and unit squares, and also as `P_{(h,1),2,n}(1)`.
- `A099234`: `H_3(n)`, also `P_{(h,1),3,n}(1)`.
- `A099235`: `H_4(n)`, also `P_{(h,1),4,n}(1)`.
- `A360090`: `H_5(n)`, also `P_{(h,1),5,n}(1)`.

### Other total tiling counts

- `A097076`: update with `P_{(2,2),3,n}(1)` and
  `OGF = 1/(1 - x - 3*x^2 - x^3)`.
- `A193147`: update with `P_{(3,2,1),2,n}(1)` and
  `OGF = 1/(1 - x - 2*x^3 - x^5)`.

### Dense tiling counts

- `A008346`: update with dense tilings `Q_{(h,1),2,n}(1)` and
  `OGF = 1/(1 - 2*x^2 - x^3)`.

### Row-refined and row-avoiding counts

- `A000045`: also appears as `S_n(0)` for `(h,1), d=2`, with
  `OGF = 1/(1 - x - x^2)`.
- `A001045`: appears as `S_n(0)` for `(k,k), d=3`, and also as
  `S_{(k,k),2,n}(1)`, with `OGF = 1/(1 - x - 2*x^2)`.
- `A006130`: update with `S_{(k,k),3,n}(1)` and
  `OGF = 1/(1 - x - 3*x^2)`.
- `A006131`: update with `S_{(k,k),4,n}(1)` and
  `OGF = 1/(1 - x - 4*x^2)`.
- `A000930`: update with `S_n(0)` for `(h,1,1), d=2`, with
  `OGF = 1/(1 - x - x^3)`.
- `A023610`: update with the coefficient of `s^1` in
  `S_{(h,1),2,n}(s,1)`. Terms, `n=0,1,2,...`:

```text
0, 0, 1, 3, 7, 15, 30, 58, 109, 201, 365, 655, 1164, 2052, 3593, 6255, 10835, 18687, 32106, 54974, 93845, 159765, 271321, 459743, 777432, 1312200
```

- `A001629`: update with the coefficient of `s^1` in
  `S_{(k,k),2,n}(s,1)`, equivalently the self-convolution of Fibonacci.
  Terms, `n=0,1,2,...`:

```text
0, 0, 1, 2, 5, 10, 20, 38, 71, 130, 235, 420, 744, 1308, 2285, 3970, 6865, 11822, 20284, 34690, 59155, 100610, 170711, 289032, 488400
```

- `A073371`: update with the coefficient of `s^1` in
  `S_{(k,k),3,n}(s,1,1)`, equivalently the self-convolution of the
  Jacobsthal sequence. Terms, `n=0,1,2,...`:

```text
0, 0, 1, 2, 7, 16, 41, 94, 219, 492, 1101, 2426, 5311, 11528, 24881, 53398, 114083, 242724, 514581, 1087410, 2291335, 4815680, 10097401, 21126862, 44117867
```

### Cylindric tiling counts

These are from the cylindric section, for `mu=(h,1)` and `d=h`.

- `A000032`: update with cylindric counts for `d=1`, listed for `n>=1`:

```text
1, 3, 4, 7, 11, 18, 29, 47, 76, 123, 199, 322, 521, 843, 1364, 2207, 3571, 5778, 9349, 15127, 24476, 39603, 64079, 103682, 167761
```

- `A286910`: update with cylindric counts for `d=2`, listed for `n>=1`:

```text
1, 5, 10, 21, 46, 98, 211, 453, 973, 2090, 4489, 9642, 20710, 44483, 95545, 205221, 440794, 946781, 2033590, 4367946, 9381907, 20151389, 43283149, 92967834, 199685521
```

## Useful OEIS wording

Possible title pattern for new entries:

```text
Number of anchor-word tilings of width n for the Ferrers tile mu=... in strip parameter d=...
```

Possible comment:

```text
Equivalently, this is P_{mu,d,n}(1), where P_{mu,d,n}(t) is the tiling
polynomial counting anchor words by number of big tiles.
```

For dense sequences:

```text
Dense version: Q_n(1) = [x^n] 1/(1 - (B_{mu,d}(x,1)-x)), where B_{mu,d}
is the fault-free block generating function.
```

For cylindric sequences:

```text
Counts cylindric anchor words for mu=(h,1) and d=h, where the left and right
boundaries of the board are identified.
```
