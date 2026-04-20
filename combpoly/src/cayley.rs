//! Cayley permutations: surjections onto initial segments.
//!
//! A Cayley permutation of size n is a sequence (a_1, ..., a_n) with positive
//! integer entries whose image is {1, ..., m} for some m, i.e., every value
//! from 1 to the maximum appears. The count is the Fubini number (A000670).

/// Generate all Cayley permutations of size n.
pub fn all_cayley_permutations(n: u8) -> Vec<Vec<u8>> {
    let m = n as usize;
    let mut result = Vec::new();
    if m == 0 {
        result.push(Vec::new());
        return result;
    }
    let mut current = Vec::with_capacity(m);
    let mut freq = vec![0usize; m]; // freq[v-1] = how many times value v appears
    gen_cayley(m, 0, &mut current, &mut freq, &mut result);
    result
}

fn gen_cayley(
    n: usize,
    max_used: u8, // current maximum value used
    current: &mut Vec<u8>,
    freq: &mut Vec<usize>,
    result: &mut Vec<Vec<u8>>,
) {
    let placed = current.len();

    if placed == n {
        result.push(current.clone());
        return;
    }

    let remaining = n - placed - 1;

    // We can use any value from 1..=max_used (existing), or introduce new
    // values up to max_used + remaining + 1 (the remaining positions can
    // cover at most `remaining` new gaps).
    let upper = (max_used as usize + remaining + 1).min(n) as u8;
    for v in 1..=upper {
        let vi = v as usize - 1;
        let new_max = if v > max_used { v } else { max_used };
        freq[vi] += 1;

        // Count values in 1..=new_max that still have freq 0
        let uncovered = freq[..new_max as usize].iter().filter(|&&f| f == 0).count();

        if uncovered <= remaining {
            current.push(v);
            gen_cayley(n, new_max, current, freq, result);
            current.pop();
        }

        freq[vi] -= 1;
    }
}

/// Iterate over all fixed-point-free Cayley permutations of size n,
/// calling `f` on each one without storing them all in memory.
/// A Cayley permutation w is fixed-point-free if w_i != i for all i.
pub fn for_each_fpf_cayley(n: u8, mut f: impl FnMut(&[u8])) {
    let m = n as usize;
    if m == 0 {
        f(&[]); // empty permutation is vacuously FPF
        return;
    }
    let mut current = Vec::with_capacity(m);
    let mut freq = vec![0usize; m];
    gen_fpf_cayley(m, 0, &mut current, &mut freq, &mut f);
}

fn gen_fpf_cayley(
    n: usize,
    max_used: u8,
    current: &mut Vec<u8>,
    freq: &mut Vec<usize>,
    f: &mut impl FnMut(&[u8]),
) {
    let placed = current.len();

    if placed == n {
        f(current);
        return;
    }

    let remaining = n - placed - 1;
    let upper = (max_used as usize + remaining + 1).min(n) as u8;
    let pos_1based = (placed + 1) as u8; // 1-indexed position

    for v in 1..=upper {
        if v == pos_1based {
            continue; // skip fixed points
        }
        let vi = v as usize - 1;
        let new_max = if v > max_used { v } else { max_used };
        freq[vi] += 1;

        let uncovered = freq[..new_max as usize].iter().filter(|&&f| f == 0).count();

        if uncovered <= remaining {
            current.push(v);
            gen_fpf_cayley(n, new_max, current, freq, f);
            current.pop();
        }

        freq[vi] -= 1;
    }
}

/// Check if a sequence is a Cayley permutation.
pub fn is_cayley_permutation(a: &[u8]) -> bool {
    if a.is_empty() {
        return true;
    }
    let max_val = *a.iter().max().unwrap() as usize;
    let mut seen = vec![false; max_val];
    for &v in a {
        if v == 0 {
            return false;
        }
        seen[v as usize - 1] = true;
    }
    seen.iter().all(|&s| s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cayley_counts() {
        // Fubini numbers: 1, 1, 3, 13, 75, 541
        assert_eq!(all_cayley_permutations(0).len(), 1);
        assert_eq!(all_cayley_permutations(1).len(), 1);
        assert_eq!(all_cayley_permutations(2).len(), 3);
        assert_eq!(all_cayley_permutations(3).len(), 13);
        assert_eq!(all_cayley_permutations(4).len(), 75);
        assert_eq!(all_cayley_permutations(5).len(), 541);
    }

    #[test]
    fn test_is_cayley_permutation() {
        assert!(is_cayley_permutation(&[1, 1, 2]));
        assert!(is_cayley_permutation(&[2, 1, 3]));
        assert!(is_cayley_permutation(&[1, 1, 1]));
        assert!(!is_cayley_permutation(&[1, 3, 3])); // missing 2
        assert!(!is_cayley_permutation(&[2, 2, 2])); // missing 1
        assert!(is_cayley_permutation(&[1]));
        assert!(is_cayley_permutation(&[]));
    }

    #[test]
    fn test_all_are_valid() {
        for cp in all_cayley_permutations(5) {
            assert!(is_cayley_permutation(&cp), "Invalid Cayley perm: {:?}", cp);
        }
    }
}
