//! Partially ordered sets (posets) on vertex set {0, 1, ..., n-1}.
//!
//! A [`Poset`] is stored as a Hasse diagram (cover relations). It supports
//! enumeration of linear extensions, order-preserving maps (backtracking and
//! frontier DP), order polytope Ehrhart polynomials and h*-vectors, and
//! P-Eulerian polynomials.
//!
//! # Order-preserving map counting
//!
//! Two algorithms are provided:
//! - **Backtracking** (`count_weak_order_preserving`): O(k^n), simple, good for n ≤ 8
//! - **Frontier DP** (`count_weak_order_preserving_dp`): O(n·k^w) where w is the max
//!   frontier width (number of live vertices with unprocessed children). Gives
//!   100–18000x speedup for chains, antichains, fences, and other narrow posets.
//!
//! Both auto-relabel to natural labeling if needed.
//!
//! # Ehrhart theory
//!
//! `order_polytope_ehrhart()` computes the Ehrhart polynomial of the order polytope
//! O(P) using the frontier DP and Stanley reciprocity (strict maps give negative
//! evaluation points for free). Uses `BigRational` for exact arithmetic — no overflow.
//!
//! `order_polytope_hstar()` returns the h*-vector. By Stanley's theorem, for
//! naturally labeled posets this equals the P-Eulerian polynomial.
//!
//! # Examples
//!
//! ```
//! use combinatoric_core::poset::Poset;
//!
//! // Diamond poset: 0 < 1, 0 < 2, 1 < 3, 2 < 3
//! let diamond = Poset::new(4, &[(0, 1), (0, 2), (1, 3), (2, 3)]);
//! assert_eq!(diamond.num_linear_extensions(), 2);
//! assert_eq!(diamond.p_eulerian_polynomial(), vec![1, 1]);
//!
//! // Order-preserving maps: Ω(diamond, 5) = 105
//! assert_eq!(diamond.count_weak_order_preserving_dp(5), 105);
//!
//! // h*-vector = P-Eulerian polynomial for naturally labeled posets
//! assert_eq!(diamond.order_polytope_hstar(), vec![1, 1]);
//! ```

use crate::partition::Partition;

// ---------------------------------------------------------------------------
// Poset type
// ---------------------------------------------------------------------------

/// A finite poset on vertices `0..n`, stored as a Hasse diagram.
///
/// Cover relations are directed edges `(a, b)` meaning a < b (a is covered by b).
/// The poset is assumed to have no cycles (it is a DAG).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Poset {
    n: usize,
    /// Cover relations (Hasse diagram): (a, b) means a < b.
    covers: Vec<(usize, usize)>,
    /// Adjacency list: children[a] = {b : a < b is a cover}.
    children: Vec<Vec<usize>>,
    /// Reverse adjacency: parents[b] = {a : a < b is a cover}.
    parents: Vec<Vec<usize>>,
}

impl Poset {
    // -- Constructors -------------------------------------------------------

    /// Create a poset from cover relations.
    ///
    /// Each pair `(a, b)` means a < b. The input should be a Hasse diagram
    /// (no transitive edges), though redundant edges are harmless.
    pub fn new(n: usize, covers: &[(usize, usize)]) -> Self {
        let mut children = vec![Vec::new(); n];
        let mut parents = vec![Vec::new(); n];
        let mut deduped = Vec::new();
        let mut seen = std::collections::BTreeSet::new();

        for &(a, b) in covers {
            if a >= n || b >= n || a == b {
                continue;
            }
            if seen.insert((a, b)) {
                children[a].push(b);
                parents[b].push(a);
                deduped.push((a, b));
            }
        }

        Poset {
            n,
            covers: deduped,
            children,
            parents,
        }
    }

    /// Chain poset: 0 < 1 < 2 < ... < (n-1).
    pub fn chain(n: usize) -> Self {
        let covers: Vec<_> = (0..n.saturating_sub(1)).map(|i| (i, i + 1)).collect();
        Poset::new(n, &covers)
    }

    /// Antichain poset on n elements (no relations).
    pub fn antichain(n: usize) -> Self {
        Poset::new(n, &[])
    }

    /// Fence (zigzag) poset: 0 < 1 > 2 < 3 > 4 < ...
    ///
    /// Covers: (0,1), (2,1), (2,3), (4,3), (4,5), ...
    pub fn fence(n: usize) -> Self {
        let mut covers = Vec::new();
        for i in 0..n.saturating_sub(1) {
            if i % 2 == 0 {
                covers.push((i, i + 1));
            } else {
                covers.push((i + 1, i));
            }
        }
        Poset::new(n, &covers)
    }

    /// Poset from a (skew) Young diagram λ/μ.
    ///
    /// Boxes are labeled 0..k in reading order (left to right, top to bottom).
    /// Box A < Box B if A is directly left of B or directly above B.
    pub fn from_skew_shape(lambda: &Partition, mu: &Partition) -> Self {
        let boxes = lambda.skew_diagram_boxes(mu);
        let n = boxes.len();
        if n == 0 {
            return Poset::new(0, &[]);
        }

        // Map (row, col) -> index
        let mut box_idx = std::collections::HashMap::new();
        for (i, &(r, c)) in boxes.iter().enumerate() {
            box_idx.insert((r, c), i);
        }

        let mut covers = Vec::new();
        for (i, &(r, c)) in boxes.iter().enumerate() {
            // Right neighbor
            if let Some(&j) = box_idx.get(&(r, c + 1)) {
                covers.push((i, j));
            }
            // Below neighbor
            if let Some(&j) = box_idx.get(&(r + 1, c)) {
                covers.push((i, j));
            }
        }

        Poset::new(n, &covers)
    }

    /// Poset from a straight shape λ (special case of skew with empty inner).
    pub fn from_shape(lambda: &Partition) -> Self {
        Self::from_skew_shape(lambda, &Partition::empty())
    }

    // -- Accessors ----------------------------------------------------------

    /// Number of elements.
    pub fn num_elements(&self) -> usize {
        self.n
    }

    /// Cover relations.
    pub fn covers(&self) -> &[(usize, usize)] {
        &self.covers
    }

    /// Minimal elements (no parents in Hasse diagram).
    pub fn minimal_elements(&self) -> Vec<usize> {
        (0..self.n).filter(|&v| self.parents[v].is_empty()).collect()
    }

    /// Maximal elements (no children in Hasse diagram).
    pub fn maximal_elements(&self) -> Vec<usize> {
        (0..self.n).filter(|&v| self.children[v].is_empty()).collect()
    }

    /// In-degree of element v in the Hasse diagram.
    pub fn in_degree(&self, v: usize) -> usize {
        self.parents[v].len()
    }

    // -- Linear extensions --------------------------------------------------

    /// All linear extensions of the poset.
    ///
    /// A linear extension is a permutation `σ` of `0..n` such that
    /// a < b in the poset implies σ^{-1}(a) < σ^{-1}(b) (i.e., a appears
    /// before b in the sequence).
    ///
    /// Uses Kahn's algorithm variant: at each step, choose any element
    /// whose in-degree (among remaining elements) is 0.
    pub fn linear_extensions(&self) -> Vec<Vec<usize>> {
        let mut result = Vec::new();
        let mut in_deg: Vec<usize> = (0..self.n).map(|v| self.parents[v].len()).collect();
        let mut current = Vec::with_capacity(self.n);
        let mut used = vec![false; self.n];

        self.linext_rec(&mut in_deg, &mut current, &mut used, &mut result);
        result
    }

    fn linext_rec(
        &self,
        in_deg: &mut Vec<usize>,
        current: &mut Vec<usize>,
        used: &mut Vec<bool>,
        result: &mut Vec<Vec<usize>>,
    ) {
        if current.len() == self.n {
            result.push(current.clone());
            return;
        }

        for v in 0..self.n {
            if !used[v] && in_deg[v] == 0 {
                used[v] = true;
                current.push(v);
                // Decrease in-degree of children
                for &c in &self.children[v] {
                    in_deg[c] -= 1;
                }

                self.linext_rec(in_deg, current, used, result);

                // Restore
                current.pop();
                used[v] = false;
                for &c in &self.children[v] {
                    in_deg[c] += 1;
                }
            }
        }
    }

    /// Number of linear extensions.
    pub fn num_linear_extensions(&self) -> usize {
        self.linear_extensions().len()
    }

    // -- Order-preserving maps ----------------------------------------------

    /// Count weakly order-preserving maps from the poset to {1, ..., k}.
    ///
    /// A map f: P → {1,...,k} is weakly order-preserving if a ≤ b in P implies f(a) ≤ f(b).
    /// The order polynomial Ω(P, k) counts these maps.
    ///
    /// Automatically relabels to natural labeling if needed. For better
    /// performance on larger posets, use [`count_weak_order_preserving_dp`].
    pub fn count_weak_order_preserving(&self, k: usize) -> usize {
        let p = self.natural_relabeling();
        if p.n == 0 {
            return 1;
        }
        let mut count = 0;
        let mut coloring = vec![0usize; p.n];
        p.order_map_rec(0, k, true, &mut coloring, &mut count);
        count
    }

    /// Count strictly order-preserving maps from the poset to {1, ..., k}.
    ///
    /// A map f: P → {1,...,k} is strictly order-preserving if a < b in P implies f(a) < f(b).
    /// The strict order polynomial Ω̄(P, k) counts these maps.
    ///
    /// Automatically relabels to natural labeling if needed. For better
    /// performance on larger posets, use [`count_strict_order_preserving_dp`].
    pub fn count_strict_order_preserving(&self, k: usize) -> usize {
        let p = self.natural_relabeling();
        if p.n == 0 {
            return 1;
        }
        let mut count = 0;
        let mut coloring = vec![0usize; p.n];
        p.order_map_rec(0, k, false, &mut coloring, &mut count);
        count
    }

    /// Backtracking helper. Requires natural labeling (all parents have index < v).
    fn order_map_rec(
        &self,
        v: usize,
        k: usize,
        weak: bool,
        coloring: &mut Vec<usize>,
        count: &mut usize,
    ) {
        if v == self.n {
            *count += 1;
            return;
        }

        // Determine lower bound from parents (all have index < v by natural labeling)
        let mut lo = 1;
        for &p in &self.parents[v] {
            let parent_val = coloring[p];
            let bound = if weak { parent_val } else { parent_val + 1 };
            if bound > lo {
                lo = bound;
            }
        }

        for c in lo..=k {
            coloring[v] = c;
            self.order_map_rec(v + 1, k, weak, coloring, count);
        }
    }

    // -- Frontier DP for order-preserving maps --------------------------------

    /// Count order-preserving maps using frontier DP (requires natural labeling).
    ///
    /// Processes vertices 0..n in order, tracking only "live" vertices
    /// (those with unprocessed children whose values constrain future vertices).
    /// States are merged when they agree on live vertex values, giving
    /// exponential speedup over backtracking for narrow posets (chains,
    /// trees, Young diagrams).
    ///
    /// Complexity: O(n · k^w) where w = max frontier width ≤ poset width.
    fn count_order_preserving_dp(&self, k: usize, weak: bool) -> usize {
        use std::collections::HashMap;

        let n = self.n;
        if n == 0 {
            return 1;
        }
        if k == 0 {
            return 0;
        }

        debug_assert!(self.is_naturally_labeled());

        let delta: usize = if weak { 0 } else { 1 };

        // max_child[v] = largest child of v in the Hasse diagram
        let mut max_child: Vec<Option<usize>> = vec![None; n];
        for &(a, b) in &self.covers {
            match max_child[a] {
                None => max_child[a] = Some(b),
                Some(cur) if b > cur => max_child[a] = Some(b),
                _ => {}
            }
        }

        // Live set: vertices u < v whose assigned value still constrains
        // some future vertex (i.e., max_child[u] >= v). Maintained incrementally.
        let mut live: Vec<usize> = Vec::new();
        let mut dp: HashMap<Vec<usize>, usize> = HashMap::new();
        dp.insert(Vec::new(), 1);

        for v in 0..n {
            // Positions of v's parents in the live vector
            let parent_positions: Vec<usize> = self.parents[v]
                .iter()
                .map(|&p| live.iter().position(|&u| u == p).unwrap())
                .collect();

            // Which live slots die after processing v (their max child is v)
            let dying_positions: Vec<usize> = live
                .iter()
                .enumerate()
                .filter(|&(_, &u)| max_child[u] == Some(v))
                .map(|(i, _)| i)
                .collect();

            // Positions to keep in the state vector
            let keep_positions: Vec<usize> = (0..live.len())
                .filter(|i| !dying_positions.contains(i))
                .collect();

            // Does v enter the next live set? (only if v has children)
            let v_enters = max_child[v].is_some();

            let mut new_dp: HashMap<Vec<usize>, usize> = HashMap::new();

            for (state, &count) in &dp {
                // Lower bound for f(v) from parent constraints
                let mut lo = 1usize;
                for &pos in &parent_positions {
                    let bound = state[pos] + delta;
                    if bound > lo {
                        lo = bound;
                    }
                }

                if lo > k {
                    continue;
                }

                if !v_enters {
                    // v has no children: all values c = lo..=k yield the same
                    // successor state, so multiply instead of iterating.
                    let num_choices = k - lo + 1;
                    let new_state: Vec<usize> =
                        keep_positions.iter().map(|&i| state[i]).collect();
                    *new_dp.entry(new_state).or_insert(0) += count * num_choices;
                } else {
                    // v enters the live set; each value creates a distinct state.
                    for c in lo..=k {
                        let mut new_state: Vec<usize> =
                            keep_positions.iter().map(|&i| state[i]).collect();
                        new_state.push(c); // v is largest, appended at end
                        *new_dp.entry(new_state).or_insert(0) += count;
                    }
                }
            }

            // Update the live set
            live.retain(|u| max_child[*u] != Some(v));
            if v_enters {
                live.push(v);
            }

            dp = new_dp;
        }

        dp.values().sum()
    }

    /// Count weakly order-preserving maps P → {1,...,k} using frontier DP.
    ///
    /// Automatically relabels the poset to natural labeling if needed.
    pub fn count_weak_order_preserving_dp(&self, k: usize) -> usize {
        let p = self.natural_relabeling();
        p.count_order_preserving_dp(k, true)
    }

    /// Count strictly order-preserving maps P → {1,...,k} using frontier DP.
    ///
    /// Automatically relabels the poset to natural labeling if needed.
    pub fn count_strict_order_preserving_dp(&self, k: usize) -> usize {
        let p = self.natural_relabeling();
        p.count_order_preserving_dp(k, false)
    }

    // -- Order polynomial ---------------------------------------------------

    /// Order polynomial as a vector of values [Ω(P,0), Ω(P,1), ..., Ω(P,n)].
    ///
    /// Since Ω(P,k) is a polynomial of degree n in k, evaluating at n+1 points
    /// determines it completely. Returns the values for Lagrange interpolation.
    pub fn order_polynomial_values(&self, max_k: usize) -> Vec<usize> {
        (0..=max_k).map(|k| self.count_weak_order_preserving(k)).collect()
    }

    // -- P-Eulerian polynomial ----------------------------------------------

    /// The P-Eulerian polynomial (descent polynomial over linear extensions).
    ///
    /// For a naturally labeled poset P, the P-Eulerian polynomial is
    ///
    /// ```text
    /// W_P(t) = Σ_{σ ∈ L(P)} t^{des(σ)}
    /// ```
    ///
    /// where L(P) is the set of linear extensions and des(σ) counts descents
    /// (positions i where σ(i) > σ(i+1)).
    pub fn p_eulerian_polynomial(&self) -> Vec<i64> {
        if self.n == 0 {
            return vec![1];
        }

        let max_des = self.n.saturating_sub(1);
        let mut coeffs = vec![0i64; max_des + 1];
        let mut in_deg: Vec<usize> = (0..self.n).map(|v| self.parents[v].len()).collect();
        let mut current = Vec::with_capacity(self.n);
        let mut used = vec![false; self.n];

        self.p_eulerian_rec(&mut in_deg, &mut current, &mut used, &mut coeffs);

        // Trim trailing zeros
        while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
            coeffs.pop();
        }
        coeffs
    }

    /// Streaming backtracking: count descents per extension without storing them.
    fn p_eulerian_rec(
        &self,
        in_deg: &mut Vec<usize>,
        current: &mut Vec<usize>,
        used: &mut Vec<bool>,
        coeffs: &mut Vec<i64>,
    ) {
        if current.len() == self.n {
            let des = (0..self.n.saturating_sub(1))
                .filter(|&i| current[i] > current[i + 1])
                .count();
            coeffs[des] += 1;
            return;
        }

        for v in 0..self.n {
            if !used[v] && in_deg[v] == 0 {
                used[v] = true;
                current.push(v);
                for &c in &self.children[v] {
                    in_deg[c] -= 1;
                }

                self.p_eulerian_rec(in_deg, current, used, coeffs);

                current.pop();
                used[v] = false;
                for &c in &self.children[v] {
                    in_deg[c] += 1;
                }
            }
        }
    }

    // -- Order polytope Ehrhart theory ----------------------------------------

    /// Compute the Ehrhart polynomial of the order polytope O(P).
    ///
    /// The order polytope O(P) ⊂ [0,1]^n has the property that its
    /// lattice points at dilation t are order-preserving maps P → {0,...,t},
    /// so Ehr(O(P), t) = Ω(P, t+1) (the order polynomial shifted by 1).
    ///
    /// Uses Stanley reciprocity: Ω̄(P, k) = (-1)^n Ω(P, -k), where Ω̄
    /// counts strict order-preserving maps. This provides evaluation
    /// points at negative dilations for free, halving the computation.
    ///
    /// Uses frontier DP internally for fast lattice point counting.
    ///
    /// Returns rational coefficients as `Ratio<BigInt>` in ascending degree order.
    pub fn order_polytope_ehrhart(&self) -> Vec<BigRat> {
        let d = self.n; // dimension of O(P)
        if d == 0 {
            return vec![BigRat::one()];
        }

        // Relabel once for the DP
        let p = self.natural_relabeling();

        let num_positive = (d + 2) / 2; // ceil((d+1)/2)
        let num_negative = (d + 1) / 2; // floor((d+1)/2)
        let sign = if d % 2 == 0 { 1i64 } else { -1i64 };

        let mut points: Vec<i64> = Vec::with_capacity(d + 1);
        let mut values: Vec<i64> = Vec::with_capacity(d + 1);

        // Positive evaluations: Ehr(t) = Ω(P, t+1)
        for t in 0..num_positive {
            let val = p.count_order_preserving_dp(t + 1, true);
            points.push(t as i64);
            values.push(val as i64);
        }

        // Negative evaluations via Stanley reciprocity:
        // Ehr(-k) = (-1)^n * Ω̄(P, k-1) for k ≥ 1.
        for k in 1..=num_negative {
            let strict_val = p.count_order_preserving_dp(k - 1, false);
            let ehr_neg = sign * strict_val as i64;
            points.push(-(k as i64));
            values.push(ehr_neg);
        }

        lagrange_interpolation_big(&points, &values)
    }

    /// Compute the h*-vector of the order polytope O(P).
    ///
    /// The h*-vector is obtained from the Ehrhart polynomial via the
    /// standard conversion. By Stanley's theorem, for naturally labeled
    /// posets, h*_i equals the number of linear extensions with exactly
    /// i descents (i.e., the P-Eulerian polynomial).
    pub fn order_polytope_hstar(&self) -> Vec<i64> {
        let ehrhart = self.order_polytope_ehrhart();
        let d = self.n;
        if d == 0 {
            return vec![1];
        }

        // h*_i = sum_{k=0}^{i} (-1)^k * C(d+1, k) * Ehr(i-k)
        let mut hstar = Vec::with_capacity(d + 1);
        for i in 0..=d {
            let mut val = BigRat::zero();
            for k in 0..=i {
                let binom = BigRat::from(binomial_big(d + 1, k));
                let sign = if k % 2 == 0 {
                    BigRat::one()
                } else {
                    -BigRat::one()
                };
                let ehr = eval_big_poly(&ehrhart, (i - k) as i64);
                val += sign * binom * ehr;
            }
            assert!(val.is_integer(), "h*-vector entry not integer at i={}", i);
            hstar.push(val.to_integer().to_i64().expect("h* entry too large for i64"));
        }
        // Trim trailing zeros
        while hstar.len() > 1 && *hstar.last().unwrap() == 0 {
            hstar.pop();
        }
        hstar
    }

    // -- k-alternating poset --------------------------------------------------

    /// Construct the poset whose linear extensions are exactly the
    /// k-alternating permutations of [n].
    ///
    /// A k-alternating permutation has descents exactly at positions
    /// divisible by k. The poset enforces:
    /// - At ascent positions j (k ∤ j): element at position j < element at position j+1
    /// - At descent positions j (k | j): element at position j > element at position j+1
    ///
    /// This is NOT a poset on [n] directly; it's the "zig-zag" poset
    /// on positions {1,...,n} where covers encode the required comparisons.
    /// Linear extensions of this poset correspond to k-alternating permutations.
    pub fn k_alternating(n: usize, k: usize) -> Self {
        // The poset is on positions {0, ..., n-1} (0-indexed).
        // For each adjacent pair (j, j+1):
        //   If position (j+1) is NOT a descent (ascent): j < j+1 in poset
        //   If position (j+1) IS a descent (k | (j+1)): j+1 < j in poset
        let mut covers = Vec::new();
        for j in 0..(n - 1) {
            let pos = j + 1; // 1-indexed position
            if pos % k == 0 {
                // descent at position pos: value at pos > value at pos+1
                // so position j+1 < position j in the poset
                covers.push((j + 1, j));
            } else {
                // ascent: position j < position j+1
                covers.push((j, j + 1));
            }
        }
        Poset::new(n, &covers)
    }

    // -- Dual / opposite poset ----------------------------------------------

    /// The dual (opposite) poset: reverse all relations.
    pub fn dual(&self) -> Self {
        let reversed: Vec<_> = self.covers.iter().map(|&(a, b)| (b, a)).collect();
        Poset::new(self.n, &reversed)
    }

    // -- Natural labeling ---------------------------------------------------

    /// Check if the poset is naturally labeled: a < b implies a < b as integers.
    pub fn is_naturally_labeled(&self) -> bool {
        self.covers.iter().all(|&(a, b)| a < b)
    }

    /// Return an isomorphic poset with a natural labeling (every cover a < b
    /// has a < b as integers). Uses topological sort.
    ///
    /// If the poset is already naturally labeled, returns a clone.
    pub fn natural_relabeling(&self) -> Self {
        if self.is_naturally_labeled() {
            return self.clone();
        }
        let topo = self.topological_order();
        let mut new_label = vec![0usize; self.n];
        for (pos, &v) in topo.iter().enumerate() {
            new_label[v] = pos;
        }
        let new_covers: Vec<_> = self
            .covers
            .iter()
            .map(|&(a, b)| (new_label[a], new_label[b]))
            .collect();
        Poset::new(self.n, &new_covers)
    }

    /// Reduce to Hasse diagram by removing transitive edges.
    ///
    /// If the input has edges a→b and a→c→b, the edge a→b is transitive
    /// and is removed.
    pub fn to_hasse_diagram(&self) -> Self {
        // Compute transitive closure via reachability
        let mut reach = vec![vec![false; self.n]; self.n];
        // Topological order
        let topo = self.topological_order();

        // Process in reverse topological order
        for &v in topo.iter().rev() {
            for &c in &self.children[v] {
                reach[v][c] = true;
                for u in 0..self.n {
                    if reach[c][u] {
                        reach[v][u] = true;
                    }
                }
            }
        }

        // An edge (a, b) is a cover iff there is no c with a < c < b
        let hasse_covers: Vec<_> = self
            .covers
            .iter()
            .copied()
            .filter(|&(a, b)| {
                !self.children[a]
                    .iter()
                    .any(|&c| c != b && reach[c][b])
            })
            .collect();

        Poset::new(self.n, &hasse_covers)
    }

    fn topological_order(&self) -> Vec<usize> {
        let mut in_deg: Vec<usize> = (0..self.n).map(|v| self.parents[v].len()).collect();
        let mut queue: Vec<usize> = (0..self.n).filter(|&v| in_deg[v] == 0).collect();
        let mut order = Vec::with_capacity(self.n);

        while let Some(v) = queue.pop() {
            order.push(v);
            for &c in &self.children[v] {
                in_deg[c] -= 1;
                if in_deg[c] == 0 {
                    queue.push(c);
                }
            }
        }
        order
    }

    // -- Comparability graph ------------------------------------------------

    /// The comparability graph: vertices are poset elements, edges connect
    /// comparable pairs (a < b or b < a in the transitive closure).
    pub fn comparability_graph(&self) -> crate::graph::Graph {
        let mut reach = vec![vec![false; self.n]; self.n];
        let topo = self.topological_order();

        for &v in topo.iter().rev() {
            for &c in &self.children[v] {
                reach[v][c] = true;
                for u in 0..self.n {
                    if reach[c][u] {
                        reach[v][u] = true;
                    }
                }
            }
        }

        let mut edges = Vec::new();
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                if reach[i][j] || reach[j][i] {
                    edges.push((i, j));
                }
            }
        }
        crate::graph::Graph::new(self.n, &edges)
    }

    /// The incomparability graph: edges connect incomparable pairs.
    /// This is the complement of the comparability graph.
    pub fn incomparability_graph(&self) -> crate::graph::Graph {
        self.comparability_graph().complement()
    }
}

impl std::fmt::Display for Poset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Poset(n={}, covers={})",
            self.n,
            self.covers.len()
        )
    }
}

// ---------------------------------------------------------------------------
// Helper functions for Ehrhart computation (BigRational — no overflow)
// ---------------------------------------------------------------------------

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{Zero, One, ToPrimitive};

type BigRat = Ratio<BigInt>;

/// Evaluate a polynomial with BigRational coefficients at an integer point.
fn eval_big_poly(coeffs: &[BigRat], x: i64) -> BigRat {
    let x = BigRat::from(BigInt::from(x));
    let mut result = BigRat::zero();
    let mut x_pow = BigRat::one();
    for c in coeffs {
        result += c * &x_pow;
        x_pow *= &x;
    }
    result
}

/// Lagrange interpolation with integer points and values.
/// Returns exact rational coefficients using BigRational (no overflow).
fn lagrange_interpolation_big(points: &[i64], values: &[i64]) -> Vec<BigRat> {
    let n = points.len();
    let pts: Vec<BigInt> = points.iter().map(|&x| BigInt::from(x)).collect();
    let mut dd: Vec<BigRat> = values
        .iter()
        .map(|&v| BigRat::from(BigInt::from(v)))
        .collect();

    // Newton's divided differences
    for j in 1..n {
        for i in (j..n).rev() {
            let diff = &dd[i] - &dd[i - 1];
            let x_diff = BigRat::from(&pts[i] - &pts[i - j]);
            dd[i] = diff / x_diff;
        }
    }

    // Convert from Newton form to standard polynomial
    // p(x) = dd[0] + dd[1]*(x-x0) + dd[2]*(x-x0)*(x-x1) + ...
    let mut coeffs: Vec<BigRat> = vec![BigRat::zero(); n];
    coeffs[0] = dd[n - 1].clone();
    for i in (0..n - 1).rev() {
        let neg_pt = BigRat::from(-&pts[i]);
        for j in (0..n - 1).rev() {
            let c = coeffs[j].clone();
            coeffs[j + 1] = &coeffs[j + 1] + &c;
            coeffs[j] = c * &neg_pt;
        }
        coeffs[0] = &coeffs[0] + &dd[i];
    }
    coeffs
}

fn binomial_big(n: usize, k: usize) -> BigInt {
    if k > n {
        return BigInt::zero();
    }
    let k = k.min(n - k);
    let mut result = BigInt::one();
    for i in 0..k {
        result = result * BigInt::from(n - i) / BigInt::from(i + 1);
    }
    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- Linear extension tests (Mathematica-verified) --

    #[test]
    fn test_chain_linext() {
        // Chain 0 < 1 < 2: one linear extension [0, 1, 2]
        let p = Poset::chain(3);
        assert_eq!(p.num_linear_extensions(), 1);
    }

    #[test]
    fn test_antichain_linext() {
        // Antichain on 3: all 6 permutations
        let p = Poset::antichain(3);
        assert_eq!(p.num_linear_extensions(), 6);
    }

    #[test]
    fn test_diamond_linext() {
        // Diamond: 0<1, 0<2, 1<3, 2<3 → 2 linear extensions
        let p = Poset::new(4, &[(0, 1), (0, 2), (1, 3), (2, 3)]);
        assert_eq!(p.num_linear_extensions(), 2);
    }

    #[test]
    fn test_v_poset_linext() {
        // V: 0<1, 0<2 → 2 linear extensions (but Mathematica says 2 with 1-indexed)
        // Actually: [0,1,2] and [0,2,1]
        let p = Poset::new(3, &[(0, 1), (0, 2)]);
        let exts = p.linear_extensions();
        assert_eq!(exts.len(), 2);
    }

    #[test]
    fn test_n_poset_linext() {
        // N-poset: 0<1, 0<2, 2<3 → 3 linear extensions
        let p = Poset::new(4, &[(0, 1), (0, 2), (2, 3)]);
        assert_eq!(p.num_linear_extensions(), 3);
    }

    #[test]
    fn test_fence_linext() {
        // Fence 4: 0<1, 2<1, 2<3 → 5 linear extensions
        let p = Poset::new(4, &[(0, 1), (2, 1), (2, 3)]);
        assert_eq!(p.num_linear_extensions(), 5);
    }

    // -- Order-preserving map tests (Mathematica-verified) --

    #[test]
    fn test_chain_order_poly() {
        // Chain 0<1<2, k=3: C(3+2,3) = 10 weakly order-preserving maps
        let p = Poset::chain(3);
        assert_eq!(p.count_weak_order_preserving(3), 10);
        assert_eq!(p.count_weak_order_preserving(4), 20);
    }

    #[test]
    fn test_antichain_order_poly() {
        // Antichain on 3, k=2: 2^3 = 8
        let p = Poset::antichain(3);
        assert_eq!(p.count_weak_order_preserving(2), 8);
    }

    #[test]
    fn test_diamond_order_poly() {
        // Diamond: Ω(P,k) values from Mathematica: 1, 6, 20, 50, 105
        let p = Poset::new(4, &[(0, 1), (0, 2), (1, 3), (2, 3)]);
        assert_eq!(p.count_weak_order_preserving(1), 1);
        assert_eq!(p.count_weak_order_preserving(2), 6);
        assert_eq!(p.count_weak_order_preserving(3), 20);
        assert_eq!(p.count_weak_order_preserving(4), 50);
        assert_eq!(p.count_weak_order_preserving(5), 105);
    }

    #[test]
    fn test_strict_order_poly() {
        // Chain 3 strict k=3: C(3,3) = 1
        let p = Poset::chain(3);
        assert_eq!(p.count_strict_order_preserving(3), 1);
        assert_eq!(p.count_strict_order_preserving(4), 4);

        // Diamond strict k=3: 1
        let d = Poset::new(4, &[(0, 1), (0, 2), (1, 3), (2, 3)]);
        assert_eq!(d.count_strict_order_preserving(3), 1);
    }

    // -- P-Eulerian polynomial tests (Mathematica-verified) --

    #[test]
    fn test_chain_p_eulerian() {
        // Chain: one extension with 0 descents → [1]
        assert_eq!(Poset::chain(3).p_eulerian_polynomial(), vec![1]);
    }

    #[test]
    fn test_antichain_p_eulerian() {
        // Antichain 3 = all 6 permutations → Eulerian [1, 4, 1]
        assert_eq!(Poset::antichain(3).p_eulerian_polynomial(), vec![1, 4, 1]);
    }

    #[test]
    fn test_diamond_p_eulerian() {
        // Diamond → [1, 1]
        let p = Poset::new(4, &[(0, 1), (0, 2), (1, 3), (2, 3)]);
        assert_eq!(p.p_eulerian_polynomial(), vec![1, 1]);
    }

    #[test]
    fn test_v_poset_p_eulerian() {
        // V-poset → [1, 1]
        let p = Poset::new(3, &[(0, 1), (0, 2)]);
        assert_eq!(p.p_eulerian_polynomial(), vec![1, 1]);
    }

    #[test]
    fn test_fence_p_eulerian() {
        // Fence 4 → [0, 4, 1] (Mathematica)
        let p = Poset::new(4, &[(0, 1), (2, 1), (2, 3)]);
        assert_eq!(p.p_eulerian_polynomial(), vec![0, 4, 1]);
    }

    #[test]
    fn test_n_poset_p_eulerian() {
        // N-poset → [1, 2]
        let p = Poset::new(4, &[(0, 1), (0, 2), (2, 3)]);
        assert_eq!(p.p_eulerian_polynomial(), vec![1, 2]);
    }

    // -- Constructor tests --

    #[test]
    fn test_fence_constructor() {
        let f = Poset::fence(4);
        // 0<1, 2<1, 2<3
        assert_eq!(f.num_elements(), 4);
        assert_eq!(f.covers().len(), 3);
        assert_eq!(f.num_linear_extensions(), 5);
    }

    #[test]
    fn test_from_shape() {
        // λ = (2,1): boxes (0,0), (0,1), (1,0)
        // Covers: (0,0)<(0,1) and (0,0)<(1,0)
        let p = Poset::from_shape(&Partition::new(vec![2, 1]));
        assert_eq!(p.num_elements(), 3);
        // Linear extensions: [0,1,2] and [0,2,1]
        assert_eq!(p.num_linear_extensions(), 2);
    }

    // -- Structural tests --

    #[test]
    fn test_minimal_maximal() {
        let p = Poset::new(4, &[(0, 1), (0, 2), (1, 3), (2, 3)]);
        assert_eq!(p.minimal_elements(), vec![0]);
        assert_eq!(p.maximal_elements(), vec![3]);

        let a = Poset::antichain(3);
        assert_eq!(a.minimal_elements(), vec![0, 1, 2]);
        assert_eq!(a.maximal_elements(), vec![0, 1, 2]);
    }

    #[test]
    fn test_dual() {
        let p = Poset::chain(3);
        let d = p.dual();
        // Dual of 0<1<2 is 2<1<0
        assert_eq!(d.num_linear_extensions(), 1);
        let exts = d.linear_extensions();
        assert_eq!(exts[0], vec![2, 1, 0]);
    }

    #[test]
    fn test_naturally_labeled() {
        assert!(Poset::chain(3).is_naturally_labeled());
        // Fence has covers (2,1) which violates natural labeling
        assert!(!Poset::fence(4).is_naturally_labeled());
    }

    #[test]
    fn test_comparability_graph() {
        // Chain 0<1<2: all pairs comparable → complete graph K3
        let g = Poset::chain(3).comparability_graph();
        assert_eq!(g.num_edges(), 3);

        // Antichain: no edges
        let g2 = Poset::antichain(3).comparability_graph();
        assert_eq!(g2.num_edges(), 0);
    }

    #[test]
    fn test_hasse_reduction() {
        // Add a transitive edge 0→2 to chain 0<1<2
        let p = Poset::new(3, &[(0, 1), (1, 2), (0, 2)]);
        let h = p.to_hasse_diagram();
        assert_eq!(h.covers().len(), 2); // 0→2 removed
    }

    #[test]
    fn test_empty_poset() {
        let p = Poset::new(0, &[]);
        assert_eq!(p.num_linear_extensions(), 1); // empty permutation
        assert_eq!(p.count_weak_order_preserving(5), 1);
    }

    // -- Natural relabeling tests --

    #[test]
    fn test_natural_relabeling_fence() {
        let f = Poset::fence(4); // covers (0,1),(2,1),(2,3) — not natural
        assert!(!f.is_naturally_labeled());
        let g = f.natural_relabeling();
        assert!(g.is_naturally_labeled());
        // Isomorphic: same number of linear extensions
        assert_eq!(g.num_linear_extensions(), f.num_linear_extensions());
    }

    #[test]
    fn test_natural_relabeling_already_natural() {
        let p = Poset::chain(4);
        assert!(p.is_naturally_labeled());
        let q = p.natural_relabeling();
        assert_eq!(p, q);
    }

    #[test]
    fn test_natural_relabeling_dual() {
        // Dual of chain 0<1<2<3 has covers (3,2),(2,1),(1,0) — all reversed
        let d = Poset::chain(4).dual();
        assert!(!d.is_naturally_labeled());
        let g = d.natural_relabeling();
        assert!(g.is_naturally_labeled());
        assert_eq!(g.num_linear_extensions(), 1);
    }

    // -- DP vs backtracking agreement tests --

    #[test]
    fn test_dp_vs_backtrack_chain() {
        let p = Poset::chain(5);
        for k in 0..=8 {
            assert_eq!(
                p.count_weak_order_preserving_dp(k),
                p.count_weak_order_preserving(k),
                "weak chain(5) k={}",
                k
            );
            assert_eq!(
                p.count_strict_order_preserving_dp(k),
                p.count_strict_order_preserving(k),
                "strict chain(5) k={}",
                k
            );
        }
    }

    #[test]
    fn test_dp_vs_backtrack_antichain() {
        let p = Poset::antichain(4);
        for k in 0..=6 {
            assert_eq!(
                p.count_weak_order_preserving_dp(k),
                p.count_weak_order_preserving(k),
                "weak antichain(4) k={}",
                k
            );
            assert_eq!(
                p.count_strict_order_preserving_dp(k),
                p.count_strict_order_preserving(k),
                "strict antichain(4) k={}",
                k
            );
        }
    }

    #[test]
    fn test_dp_vs_backtrack_diamond() {
        let p = Poset::new(4, &[(0, 1), (0, 2), (1, 3), (2, 3)]);
        for k in 0..=8 {
            assert_eq!(
                p.count_weak_order_preserving_dp(k),
                p.count_weak_order_preserving(k),
                "weak diamond k={}",
                k
            );
            assert_eq!(
                p.count_strict_order_preserving_dp(k),
                p.count_strict_order_preserving(k),
                "strict diamond k={}",
                k
            );
        }
    }

    #[test]
    fn test_dp_vs_backtrack_fence() {
        // Fence is NOT naturally labeled — the old backtracking silently
        // gives wrong results (skips parents with index > v). Compare the
        // DP (which auto-relabels) against the relabeled backtracking.
        let p = Poset::fence(5);
        assert!(!p.is_naturally_labeled());
        let q = p.natural_relabeling();
        assert!(q.is_naturally_labeled());
        for k in 0..=7 {
            assert_eq!(
                p.count_weak_order_preserving_dp(k),
                q.count_weak_order_preserving(k),
                "weak fence(5) k={}",
                k
            );
            assert_eq!(
                p.count_strict_order_preserving_dp(k),
                q.count_strict_order_preserving(k),
                "strict fence(5) k={}",
                k
            );
        }
    }

    #[test]
    fn test_dp_vs_backtrack_young_diagram() {
        // Shape (3,2,1) = 6 boxes
        let p = Poset::from_shape(&Partition::new(vec![3, 2, 1]));
        for k in 0..=8 {
            assert_eq!(
                p.count_weak_order_preserving_dp(k),
                p.count_weak_order_preserving(k),
                "weak (3,2,1) k={}",
                k
            );
            assert_eq!(
                p.count_strict_order_preserving_dp(k),
                p.count_strict_order_preserving(k),
                "strict (3,2,1) k={}",
                k
            );
        }
    }

    // -- Ehrhart polynomial and h*-vector tests (Mathematica-verified) --
    //
    // Mathematica verification:
    //   Needs["Combinatorica`"];
    //   (* Chain(3): O(P) = {0 ≤ x1 ≤ x2 ≤ x3 ≤ 1} = standard simplex *)
    //   (* Ehr(t) = C(t+3,3) = (t+1)(t+2)(t+3)/6 *)
    //   (* Antichain(3): O(P) = [0,1]^3, Ehr(t) = (t+1)^3 *)

    fn rat(n: i64, d: i64) -> BigRat {
        BigRat::new(BigInt::from(n), BigInt::from(d))
    }

    #[test]
    fn test_ehrhart_chain3() {
        // Chain 0<1<2: O(P) is the 3-simplex {0 ≤ x1 ≤ x2 ≤ x3 ≤ 1}
        // Ehr(t) = C(t+3,3) = (t+1)(t+2)(t+3)/6 = t^3/6 + t^2 + 11t/6 + 1
        let p = Poset::chain(3);
        let ehr = p.order_polytope_ehrhart();
        assert_eq!(ehr.len(), 4); // degree 3
        assert_eq!(ehr[0], rat(1, 1)); // constant = 1
        assert_eq!(ehr[3], rat(1, 6)); // leading = 1/6
        // Verify at t=1: C(4,3)=4
        assert_eq!(eval_big_poly(&ehr, 1), rat(4, 1));
        // Verify at t=2: C(5,3)=10
        assert_eq!(eval_big_poly(&ehr, 2), rat(10, 1));
    }

    #[test]
    fn test_ehrhart_antichain3() {
        // Antichain(3): O(P) = [0,1]^3, Ehr(t) = (t+1)^3
        let p = Poset::antichain(3);
        let ehr = p.order_polytope_ehrhart();
        // (t+1)^3 = t^3 + 3t^2 + 3t + 1
        assert_eq!(ehr.len(), 4);
        assert_eq!(ehr[0], rat(1, 1));
        assert_eq!(ehr[1], rat(3, 1));
        assert_eq!(ehr[2], rat(3, 1));
        assert_eq!(ehr[3], rat(1, 1));
    }

    #[test]
    fn test_hstar_chain3() {
        // Chain(3): h* = [1] (simplex has h*=(1,0,0,...))
        // Actually h* = [1, 0, 0] but we trim trailing zeros → [1]
        let p = Poset::chain(3);
        assert_eq!(p.order_polytope_hstar(), vec![1]);
    }

    #[test]
    fn test_hstar_antichain3() {
        // Antichain(3): h* = Eulerian numbers A(3,k) = [1, 4, 1]
        let p = Poset::antichain(3);
        assert_eq!(p.order_polytope_hstar(), vec![1, 4, 1]);
    }

    #[test]
    fn test_hstar_diamond() {
        // Diamond: h* = P-Eulerian = [1, 1]
        let p = Poset::new(4, &[(0, 1), (0, 2), (1, 3), (2, 3)]);
        assert_eq!(p.order_polytope_hstar(), vec![1, 1]);
    }

    #[test]
    fn test_hstar_equals_p_eulerian_naturally_labeled() {
        // Stanley's theorem: for naturally labeled posets, h* = P-Eulerian poly
        let posets = vec![
            Poset::chain(4),
            Poset::antichain(4),
            Poset::new(4, &[(0, 1), (0, 2), (1, 3), (2, 3)]), // diamond
            Poset::new(4, &[(0, 1), (0, 2), (2, 3)]),          // N-poset
            Poset::from_shape(&Partition::new(vec![3, 2])),
            Poset::from_shape(&Partition::new(vec![2, 2, 1])),
        ];
        for p in posets {
            assert!(p.is_naturally_labeled());
            assert_eq!(
                p.order_polytope_hstar(),
                p.p_eulerian_polynomial(),
                "h* != P-Eulerian for {:?}",
                p
            );
        }
    }

    #[test]
    fn test_ehrhart_shape_22() {
        // λ = (2,2): this is the diamond poset.
        // Ω(P,k) for diamond at k=1..5: 1, 6, 20, 50, 105
        // Ehr(t) = Ω(P, t+1). Verify evaluations.
        let p = Poset::from_shape(&Partition::new(vec![2, 2]));
        let ehr = p.order_polytope_ehrhart();
        assert_eq!(eval_big_poly(&ehr, 0), rat(1, 1));
        assert_eq!(eval_big_poly(&ehr, 1), rat(6, 1));
        assert_eq!(eval_big_poly(&ehr, 2), rat(20, 1));
        assert_eq!(eval_big_poly(&ehr, 3), rat(50, 1));
    }

    #[test]
    fn test_ehrhart_fence10() {
        // fence(10): n=10, previously overflowed with i64 arithmetic.
        // Verify h* = P-Eulerian of the naturally relabeled poset.
        let p = Poset::fence(10);
        let hstar = p.order_polytope_hstar();
        let nat = p.natural_relabeling();
        let pe = nat.p_eulerian_polynomial();
        assert_eq!(hstar, pe);
    }
}
