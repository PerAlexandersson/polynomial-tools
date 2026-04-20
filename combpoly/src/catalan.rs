/// Enumerate Dyck paths (as area sequences) and compute peak-nesting statistics.
///
/// An area sequence is a list (a_1, ..., a_n) with a_i < i and a_{i+1} <= a_i + 1.
/// These are in bijection with Dyck paths of size n.
/// Visit all area sequences of length `n` without storing them all in memory.
pub fn for_each_area_sequence<F>(n: usize, mut visitor: F)
where
    F: FnMut(&[u8]),
{
    let mut current = Vec::with_capacity(n);
    generate_area_sequences_with(n, &mut current, &mut visitor);
}

fn generate_area_sequences_with<F>(n: usize, current: &mut Vec<u8>, visitor: &mut F)
where
    F: FnMut(&[u8]),
{
    let i = current.len(); // 0-indexed position, so a_{i+1} in 1-indexed
    if i == n {
        visitor(current);
        return;
    }

    // a_{i+1} < i+1, and if i > 0 then a_{i+1} <= a_i + 1
    let max_val = if i == 0 {
        0
    } else {
        current[i - 1] as usize + 1
    };
    let max_val = max_val.min(i); // a_{i+1} < i+1 means a_{i+1} <= i

    for v in 0..=max_val {
        current.push(v as u8);
        generate_area_sequences_with(n, current, visitor);
        current.pop();
    }
}

/// Generate all area sequences of length `n`.
pub fn all_area_sequences(n: usize) -> Vec<Vec<u8>> {
    let mut result = Vec::new();
    for_each_area_sequence(n, |a| result.push(a.to_vec()));
    result
}

/// Convert an area sequence to its Dyck word in the alphabet `{N,E}`.
pub fn area_sequence_to_dyck_word(a: &[u8]) -> String {
    let n = a.len();
    let mut north_counts = vec![0usize; n];

    for (i, &v) in a.iter().enumerate() {
        let east_before = i - v as usize;
        north_counts[east_before] += 1;
    }

    let mut word = String::with_capacity(2 * n);
    for &count in &north_counts {
        for _ in 0..count {
            word.push('N');
        }
        word.push('E');
    }

    word
}

/// Compute the peak-cliques of an area sequence.
///
/// Returns a list of (start, end) pairs where each peak-clique is {start, start+1, ..., end}
/// using 1-based indexing.
pub fn peak_cliques(a: &[u8]) -> Vec<(usize, usize)> {
    let n = a.len();
    let mut cliques = Vec::new();

    for j in 0..n {
        // 1-indexed: position j+1, value a[j]
        // Peak at j+1 iff j+1 == n or a[j+1] <= a[j]
        let is_peak = if j + 1 == n { true } else { a[j + 1] <= a[j] };

        if is_peak {
            // Peak-clique: {(j+1) - a[j], ..., j+1} in 1-based
            let start = (j + 1) - a[j] as usize; // 1-based
            let end = j + 1; // 1-based
            cliques.push((start, end));
        }
    }

    cliques
}

/// Return a triple of peak-intervals witnessing the `Q6 \cong 332/1` obstruction.
///
/// In the interval convention used for Dyck-path peaks, this means finding
/// three peak-intervals that standardize to `12`, `234`, and `45`.
pub fn q6_peak_witness(a: &[u8]) -> Option<[(usize, usize); 3]> {
    let peaks = peak_cliques(a);

    for i in 0..peaks.len() {
        let (a1, b1) = peaks[i];
        for j in (i + 1)..peaks.len() {
            let (a2, b2) = peaks[j];

            if !(a1 < a2 && a2 <= b1) {
                continue;
            }

            for &(a3, b3) in peaks.iter().skip(j + 1) {
                if a3 <= b1 + 1 {
                    continue;
                }
                if a3 <= b2 && b2 < b3 {
                    return Some([(a1, b1), (a2, b2), (a3, b3)]);
                }
            }
        }
    }

    None
}

/// Check whether the Dyck path contains the `Q6 \cong 332/1` peak-subpath.
pub fn has_q6_peak_subpath(a: &[u8]) -> bool {
    q6_peak_witness(a).is_some()
}

/// Count `Q6 \cong 332/1` peak-subpaths in the Dyck path.
///
/// Each instance is a triple of peaks whose intervals standardize to `12`, `234`, and `45`.
pub fn q6_peak_subpath_count(a: &[u8]) -> usize {
    let peaks = peak_cliques(a);
    let mut count = 0usize;

    for i in 0..peaks.len() {
        let (a1, b1) = peaks[i];
        for j in (i + 1)..peaks.len() {
            let (a2, b2) = peaks[j];

            if !(a1 < a2 && a2 <= b1) {
                continue;
            }

            for &(a3, b3) in peaks.iter().skip(j + 1) {
                if a3 <= b1 + 1 {
                    continue;
                }
                if a3 <= b2 && b2 < b3 {
                    count += 1;
                }
            }
        }
    }

    count
}

/// Number of peaks of the Dyck path, equivalently the number of peak-intervals.
pub fn num_peaks(a: &[u8]) -> usize {
    peak_cliques(a).len()
}

/// Count Dyck paths of semilength `n` that contain the `Q6` obstruction.
///
/// Returns `(containing, total)`.
pub fn q6_counts(n: usize) -> (i64, i64) {
    let mut total = 0i64;
    let mut containing = 0i64;

    for_each_area_sequence(n, |a| {
        total += 1;
        if has_q6_peak_subpath(a) {
            containing += 1;
        }
    });

    (containing, total)
}

/// Peak-generating polynomial for `Q6`-avoiding Dyck paths of semilength `n`.
///
/// Returns coefficients `[c_0, c_1, ..., c_n]` where
/// `c_k = #{Dyck paths of semilength n avoiding Q6 with exactly k peaks}`.
pub fn q6_avoiding_peak_poly(n: usize) -> Vec<i64> {
    let mut coeffs = vec![0i64; 1];

    for_each_area_sequence(n, |a| {
        if !has_q6_peak_subpath(a) {
            let peaks = num_peaks(a);
            while coeffs.len() <= peaks {
                coeffs.push(0);
            }
            coeffs[peaks] += 1;
        }
    });

    while coeffs.len() <= n {
        coeffs.push(0);
    }
    coeffs
}

/// Generating polynomial for Dyck paths of semilength `n` by number of `Q6` peak-subpaths.
///
/// Returns coefficients `[c_0, c_1, ...]` where
/// `c_k = #{Dyck paths of semilength n with exactly k Q6 peak-subpaths}`.
pub fn q6_peak_count_poly(n: usize) -> Vec<i64> {
    let mut coeffs = vec![0i64; 1];

    for_each_area_sequence(n, |a| {
        let q6_count = q6_peak_subpath_count(a);
        while coeffs.len() <= q6_count {
            coeffs.push(0);
        }
        coeffs[q6_count] += 1;
    });

    coeffs
}

/// Dyck paths of semilength `n` maximizing the number of `Q6` peak-subpaths.
///
/// Returns `(max_count, area_sequences)`.
pub fn q6_peak_maximizers(n: usize) -> (usize, Vec<Vec<u8>>) {
    let mut max_count = 0usize;
    let mut maximizers = Vec::new();

    for_each_area_sequence(n, |a| {
        let q6_count = q6_peak_subpath_count(a);
        if q6_count > max_count {
            max_count = q6_count;
            maximizers.clear();
            maximizers.push(a.to_vec());
        } else if q6_count == max_count {
            maximizers.push(a.to_vec());
        }
    });

    (max_count, maximizers)
}

/// Compute the peak-nesting of an area sequence (maximum over all vertices).
pub fn peak_nesting(a: &[u8]) -> usize {
    let n = a.len();
    if n == 0 {
        return 0;
    }

    let cliques = peak_cliques(a);

    // Count how many peak-cliques contain each vertex
    let mut nesting = vec![0usize; n + 1]; // 1-indexed
    for &(start, end) in &cliques {
        for nesting_v in &mut nesting[start..=end] {
            *nesting_v += 1;
        }
    }

    *nesting.iter().max().unwrap()
}

/// Check whether the unit interval graph of an area sequence is connected.
///
/// Connected iff a_j >= 1 for all j >= 2 (i.e., `a[j] >= 1` for j >= 1 in 0-indexed).
pub fn is_connected(a: &[u8]) -> bool {
    // Empty graph is not connected by convention.
    // a[0] = 0 always. Check a[1], a[2], ..., a[n-1] >= 1.
    !a.is_empty() && (a.len() == 1 || a[1..].iter().all(|&v| v >= 1))
}

/// Compute the peak-nesting polynomial pnest_n(t) over all area sequences of length n.
///
/// Returns coefficients [c_0, c_1, ..., c_d] where the polynomial is c_0 + c_1*t + ... + c_d*t^d.
pub fn peak_nesting_poly(n: usize) -> Vec<i64> {
    let (all, _conn) = peak_nesting_polys(n);
    all
}

/// Compute the connected peak-nesting polynomial pnest_n^c(t).
pub fn peak_nesting_poly_connected(n: usize) -> Vec<i64> {
    let (_all, conn) = peak_nesting_polys(n);
    conn
}

/// Compute both pnest_n(t) and pnest_n^c(t) in a single streaming pass.
///
/// Avoids storing all area sequences in memory.
pub fn peak_nesting_polys(n: usize) -> (Vec<i64>, Vec<i64>) {
    let mut all_coeffs: Vec<i64> = Vec::new();
    let mut conn_coeffs: Vec<i64> = Vec::new();
    for_each_area_sequence(n, |a| {
        let pn = peak_nesting(a);
        while all_coeffs.len() <= pn {
            all_coeffs.push(0);
        }
        all_coeffs[pn] += 1;

        if is_connected(a) {
            while conn_coeffs.len() <= pn {
                conn_coeffs.push(0);
            }
            conn_coeffs[pn] += 1;
        }
    });
    if all_coeffs.is_empty() {
        all_coeffs.push(0);
    }
    if conn_coeffs.is_empty() {
        conn_coeffs.push(0);
    }
    (all_coeffs, conn_coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_catalan_counts() {
        // Catalan numbers: 1, 1, 2, 5, 14, 42, 132, 429
        assert_eq!(all_area_sequences(0).len(), 1);
        assert_eq!(all_area_sequences(1).len(), 1);
        assert_eq!(all_area_sequences(2).len(), 2);
        assert_eq!(all_area_sequences(3).len(), 5);
        assert_eq!(all_area_sequences(4).len(), 14);
        assert_eq!(all_area_sequences(5).len(), 42);
        assert_eq!(all_area_sequences(6).len(), 132);
        assert_eq!(all_area_sequences(7).len(), 429);
    }

    #[test]
    fn test_peak_cliques_example() {
        // From the paper: area sequence (0,1,1,2,3,3,2,0)
        // Peak-cliques: {1,2}, {2,3,4,5}, {3,4,5,6}, {5,6,7}, {8}
        let a = vec![0, 1, 1, 2, 3, 3, 2, 0];
        let cliques = peak_cliques(&a);
        assert_eq!(cliques, vec![(1, 2), (2, 5), (3, 6), (5, 7), (8, 8)]);
    }

    #[test]
    fn test_area_sequence_to_dyck_word() {
        let a = vec![0, 1, 1, 2, 1];
        assert_eq!(area_sequence_to_dyck_word(&a), "NNENNEENEE");
    }

    #[test]
    fn test_q6_minimal_witness() {
        let a = vec![0, 1, 1, 2, 1];
        assert_eq!(peak_cliques(&a), vec![(1, 2), (2, 4), (4, 5)]);
        assert_eq!(q6_peak_witness(&a), Some([(1, 2), (2, 4), (4, 5)]));
        assert!(has_q6_peak_subpath(&a));
        assert_eq!(q6_peak_subpath_count(&a), 1);
    }

    #[test]
    fn test_q6_absent_below_size_five() {
        for n in 0..5 {
            let (containing, total) = q6_counts(n);
            assert_eq!(containing, 0, "unexpected Q6 witness at n={n}");
            assert!(total >= 1);
        }
    }

    #[test]
    fn test_q6_avoiding_peak_poly_below_size_five() {
        assert_eq!(q6_avoiding_peak_poly(1), vec![0, 1]);
        assert_eq!(q6_avoiding_peak_poly(2), vec![0, 1, 1]);
        assert_eq!(q6_avoiding_peak_poly(3), vec![0, 1, 3, 1]);
        assert_eq!(q6_avoiding_peak_poly(4), vec![0, 1, 6, 6, 1]);
    }

    #[test]
    fn test_q6_peak_count_poly_small() {
        assert_eq!(q6_peak_count_poly(1), vec![1]);
        assert_eq!(q6_peak_count_poly(2), vec![2]);
        assert_eq!(q6_peak_count_poly(3), vec![5]);
        assert_eq!(q6_peak_count_poly(4), vec![14]);
        assert_eq!(q6_peak_count_poly(5), vec![41, 1]);
    }

    #[test]
    fn test_peak_nesting_example() {
        // pnest((0,1,1,2,3,3,2,0)) = 3
        let a = vec![0, 1, 1, 2, 3, 3, 2, 0];
        assert_eq!(peak_nesting(&a), 3);
    }

    #[test]
    fn test_connected_example() {
        assert!(is_connected(&[0, 1, 1, 2, 3, 3, 2])); // all a[j]>=1 for j>=1
        assert!(!is_connected(&[0, 1, 1, 2, 3, 3, 2, 0])); // a[7]=0
        assert!(is_connected(&[0])); // single vertex
        assert!(!is_connected(&[])); // empty: not connected
    }

    #[test]
    fn test_n4_all_pnest() {
        // From the paper's detailed example:
        // 8 area sequences with pnest=1, 6 with pnest=2
        // pnest_4(t) = 8t + 6t^2
        let poly = peak_nesting_poly(4);
        assert_eq!(poly, vec![0, 8, 6]);
    }

    #[test]
    fn test_n4_connected_pnest() {
        // pnest_4^c(t) = t + 4t^2
        let poly = peak_nesting_poly_connected(4);
        assert_eq!(poly, vec![0, 1, 4]);
    }

    #[test]
    fn test_pnest_poly_table() {
        // Verify against the paper's Table for n=0..10
        let expected_all: Vec<Vec<i64>> = vec![
            vec![1],                           // n=0
            vec![0, 1],                        // n=1
            vec![0, 2],                        // n=2
            vec![0, 4, 1],                     // n=3
            vec![0, 8, 6],                     // n=4
            vec![0, 16, 25, 1],                // n=5
            vec![0, 32, 90, 10],               // n=6
            vec![0, 64, 301, 63, 1],           // n=7
            vec![0, 128, 966, 322, 14],        // n=8
            vec![0, 256, 3025, 1463, 117, 1],  // n=9
            vec![0, 512, 9330, 6174, 762, 18], // n=10
        ];

        let expected_conn: Vec<Vec<i64>> = vec![
            vec![0],                         // n=0
            vec![0, 1],                      // n=1
            vec![0, 1],                      // n=2
            vec![0, 1, 1],                   // n=3
            vec![0, 1, 4],                   // n=4
            vec![0, 1, 12, 1],               // n=5
            vec![0, 1, 33, 8],               // n=6
            vec![0, 1, 88, 42, 1],           // n=7
            vec![0, 1, 232, 184, 12],        // n=8
            vec![0, 1, 609, 731, 88, 1],     // n=9
            vec![0, 1, 1596, 2737, 512, 16], // n=10
        ];

        for n in 0..=10 {
            let poly = peak_nesting_poly(n);
            assert_eq!(poly, expected_all[n], "pnest_{n}(t) mismatch");

            let poly_c = peak_nesting_poly_connected(n);
            assert_eq!(poly_c, expected_conn[n], "pnest_{n}^c(t) mismatch");
        }
    }

    #[test]
    fn test_n4_detailed() {
        // Verify each area sequence of length 4 and its peak-nesting
        let expected: Vec<(Vec<u8>, usize)> = vec![
            (vec![0, 0, 0, 0], 1),
            (vec![0, 1, 0, 0], 1),
            (vec![0, 0, 1, 0], 1),
            (vec![0, 1, 1, 0], 2),
            (vec![0, 1, 2, 0], 1),
            (vec![0, 0, 0, 1], 1),
            (vec![0, 1, 0, 1], 1),
            (vec![0, 0, 1, 1], 2),
            (vec![0, 0, 1, 2], 1),
            (vec![0, 1, 1, 1], 2),
            (vec![0, 1, 2, 1], 2),
            (vec![0, 1, 1, 2], 2),
            (vec![0, 1, 2, 2], 2),
            (vec![0, 1, 2, 3], 1),
        ];

        for (a, exp_pn) in &expected {
            let pn = peak_nesting(a);
            assert_eq!(
                pn, *exp_pn,
                "peak_nesting({:?}) = {} (expected {})",
                a, pn, exp_pn
            );
        }
    }
}
