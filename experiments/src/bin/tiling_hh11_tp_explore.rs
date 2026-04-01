/// Investigate total positivity of M_{n,j} = sum_r C(hr, j) * C(j+1, n-2j-2r)
/// for the (h,h,1,1) tiling family.
///
/// 1. Check if the transfer matrix for the reduced state machine is TNN.
/// 2. Check individual term contributions in 2x2 minors.
/// 3. Print M_{n,j} for h=2 and look for closed-form patterns.

fn binom(n: i128, k: i128) -> i128 {
    if k < 0 || k > n || n < 0 {
        return 0;
    }
    let k = k.min(n - k);
    if k == 0 {
        return 1;
    }
    let mut result: i128 = 1;
    for i in 0..k {
        result = result * (n - i) / (i + 1);
    }
    result
}

/// M_{n,j} = sum_r C(h*r, j) * C(j+1, n - 2*j - 2*r)
fn m_entry(h: i128, n: i128, j: i128) -> i128 {
    let mut total: i128 = 0;
    // r ranges: we need h*r >= 0 and n - 2*j - 2*r >= 0
    // so r >= 0 and r <= (n - 2*j) / 2
    let r_max = if n >= 2 * j { (n - 2 * j) / 2 } else { -1 };
    for r in 0..=r_max {
        total += binom(h * r, j) * binom(j + 1, n - 2 * j - 2 * r);
    }
    total
}

/// Individual term in M_{n,j} for given r
fn m_term(h: i128, n: i128, j: i128, r: i128) -> i128 {
    binom(h * r, j) * binom(j + 1, n - 2 * j - 2 * r)
}

// ============================================================
// Part 1: Transfer matrix TNN check
// ============================================================

/// States for (h,h,1,1) with d=h:
/// 0 = empty
/// For v=1..h, d=1,2,3: state index = 1 + (v-1)*3 + (d-1)
/// Total states: 1 + 3*h

fn state_index(v: usize, d: usize) -> usize {
    // v is 1-indexed, d is 1-indexed
    1 + (v - 1) * 3 + (d - 1)
}

fn state_name(idx: usize, h: usize) -> String {
    if idx == 0 {
        return "empty".to_string();
    }
    let idx0 = idx - 1;
    let v = idx0 / 3 + 1;
    let d = idx0 % 3 + 1;
    format!("({},{})", v, d)
}

/// Build transfer matrix for given h, with weight t for tile starts.
/// Returns matrix as Vec<Vec<f64>> in a given state ordering.
fn build_transfer_matrix(h: usize, t: f64, ordering: &[usize]) -> Vec<Vec<f64>> {
    let n = 1 + 3 * h;
    // Build in natural ordering first
    let mut mat = vec![vec![0.0f64; n]; n];

    // From empty (state 0):
    //   to (v,1) with weight t for any v=1..h
    //   to empty with weight 1
    mat[0][0] = 1.0; // empty -> empty
    for v in 1..=h {
        mat[0][state_index(v, 1)] = t; // empty -> (v,1)
    }

    for v in 1..=h {
        // From (v,1): to (v,2) with weight 1
        mat[state_index(v, 1)][state_index(v, 2)] = 1.0;

        // From (v,2): to (v',1) with weight t for v'>v, or to (v,3) with weight 1
        mat[state_index(v, 2)][state_index(v, 3)] = 1.0;
        for vp in (v + 1)..=h {
            mat[state_index(v, 2)][state_index(vp, 1)] = t;
        }

        // From (v,3): to (v',1) with weight t for v'>v, or to empty with weight 1
        mat[state_index(v, 3)][0] = 1.0;
        for vp in (v + 1)..=h {
            mat[state_index(v, 3)][state_index(vp, 1)] = t;
        }
    }

    // Now reorder according to the given ordering
    let m = ordering.len();
    let mut reordered = vec![vec![0.0f64; m]; m];
    for i in 0..m {
        for j in 0..m {
            reordered[i][j] = mat[ordering[i]][ordering[j]];
        }
    }
    reordered
}

/// Check all minors up to size k of a matrix for non-negativity.
/// Returns (all_nonneg, first_negative_minor_info)
fn check_tnn(mat: &[Vec<f64>], max_minor_size: usize) -> (bool, String) {
    let n = mat.len();
    let mut all_nonneg = true;
    let mut info = String::new();

    for size in 1..=max_minor_size.min(n) {
        let row_subsets = subsets(n, size);
        let col_subsets = subsets(n, size);
        for rows in &row_subsets {
            for cols in &col_subsets {
                let det = minor_det(mat, rows, cols);
                if det < -1e-10 {
                    if all_nonneg {
                        info = format!(
                            "Negative {}x{} minor: rows={:?}, cols={:?}, det={:.6}",
                            size, size, rows, cols, det
                        );
                    }
                    all_nonneg = false;
                }
            }
        }
    }
    (all_nonneg, info)
}

fn subsets(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut current = Vec::new();
    subsets_helper(n, k, 0, &mut current, &mut result);
    result
}

fn subsets_helper(
    n: usize,
    k: usize,
    start: usize,
    current: &mut Vec<usize>,
    result: &mut Vec<Vec<usize>>,
) {
    if current.len() == k {
        result.push(current.clone());
        return;
    }
    for i in start..n {
        current.push(i);
        subsets_helper(n, k, i + 1, current, result);
        current.pop();
    }
}

fn minor_det(mat: &[Vec<f64>], rows: &[usize], cols: &[usize]) -> f64 {
    let n = rows.len();
    if n == 1 {
        return mat[rows[0]][cols[0]];
    }
    if n == 2 {
        return mat[rows[0]][cols[0]] * mat[rows[1]][cols[1]]
            - mat[rows[0]][cols[1]] * mat[rows[1]][cols[0]];
    }
    // General case: Leibniz formula via permutations
    let perms = permutations(n);
    let mut det = 0.0;
    for (perm, sign) in &perms {
        let mut prod = *sign as f64;
        for i in 0..n {
            prod *= mat[rows[i]][cols[perm[i]]];
        }
        det += prod;
    }
    det
}

fn permutations(n: usize) -> Vec<(Vec<usize>, i32)> {
    let mut result = Vec::new();
    let mut perm: Vec<usize> = (0..n).collect();
    loop {
        let sign = perm_sign(&perm);
        result.push((perm.clone(), sign));
        if !next_permutation(&mut perm) {
            break;
        }
    }
    result
}

fn perm_sign(perm: &[usize]) -> i32 {
    let n = perm.len();
    let mut inversions = 0;
    for i in 0..n {
        for j in (i + 1)..n {
            if perm[i] > perm[j] {
                inversions += 1;
            }
        }
    }
    if inversions % 2 == 0 {
        1
    } else {
        -1
    }
}

fn next_permutation(arr: &mut Vec<usize>) -> bool {
    let n = arr.len();
    if n <= 1 {
        return false;
    }
    let mut i = n - 1;
    while i > 0 && arr[i - 1] >= arr[i] {
        i -= 1;
    }
    if i == 0 {
        return false;
    }
    let mut j = n - 1;
    while arr[j] <= arr[i - 1] {
        j -= 1;
    }
    arr.swap(i - 1, j);
    arr[i..].reverse();
    true
}

// ============================================================
// Part 2: Individual term analysis in 2x2 minors
// ============================================================

fn check_minor_terms(h: i128, n1: i128, n2: i128, j1: i128, j2: i128) {
    println!(
        "\n--- Individual term analysis for h={}, (n1,n2,j1,j2)=({},{},{},{}) ---",
        h, n1, n2, j1, j2
    );

    // M_{n1,j1}*M_{n2,j2} - M_{n1,j2}*M_{n2,j1}
    // = sum_{r,s} [T(n1,j1,r)*T(n2,j2,s) - T(n1,j2,r)*T(n2,j2,s)]
    // Wait, we need to be more careful. The minor is:
    // sum_r T(n1,j1,r) * sum_s T(n2,j2,s) - sum_r T(n1,j2,r) * sum_s T(n2,j1,s)
    // = sum_{r,s} [T(n1,j1,r)*T(n2,j2,s) - T(n1,j2,r)*T(n2,j1,s)]

    let r1_max = if n1 >= 2 * j1 {
        ((n1 - 2 * j1) / 2) as usize
    } else {
        0
    };
    let r2_max = if n1 >= 2 * j2 {
        ((n1 - 2 * j2) / 2) as usize
    } else {
        0
    };
    let s1_max = if n2 >= 2 * j2 {
        ((n2 - 2 * j2) / 2) as usize
    } else {
        0
    };
    let s2_max = if n2 >= 2 * j1 {
        ((n2 - 2 * j1) / 2) as usize
    } else {
        0
    };

    let rmax = r1_max.max(r2_max);
    let smax = s1_max.max(s2_max);

    let mut total_pos: i128 = 0;
    let mut total_neg: i128 = 0;
    let mut neg_count = 0;
    let mut pos_count = 0;
    let mut zero_count = 0;

    println!("  (r, s) -> T(n1,j1,r)*T(n2,j2,s) - T(n1,j2,r)*T(n2,j1,s)");
    for r in 0..=(rmax as i128) {
        for s in 0..=(smax as i128) {
            let a = m_term(h, n1, j1, r) * m_term(h, n2, j2, s);
            let b = m_term(h, n1, j2, r) * m_term(h, n2, j1, s);
            let diff = a - b;
            if diff != 0 {
                println!("  r={}, s={}: {} - {} = {}", r, s, a, b, diff);
            }
            if diff > 0 {
                total_pos += diff;
                pos_count += 1;
            } else if diff < 0 {
                total_neg += diff;
                neg_count += 1;
            } else {
                zero_count += 1;
            }
        }
    }

    let minor_val =
        m_entry(h, n1, j1) * m_entry(h, n2, j2) - m_entry(h, n1, j2) * m_entry(h, n2, j1);
    println!(
        "  Summary: pos_terms={}, neg_terms={}, zero_terms={}",
        pos_count, neg_count, zero_count
    );
    println!(
        "  total_pos={}, total_neg={}, net={}",
        total_pos,
        total_neg,
        total_pos + total_neg
    );
    println!(
        "  M[{},{}]*M[{},{}] - M[{},{}]*M[{},{}] = {}*{} - {}*{} = {}",
        n1,
        j1,
        n2,
        j2,
        n1,
        j2,
        n2,
        j1,
        m_entry(h, n1, j1),
        m_entry(h, n2, j2),
        m_entry(h, n1, j2),
        m_entry(h, n2, j1),
        minor_val
    );

    if neg_count > 0 {
        println!("  ==> CANCELLATIONS present! Not term-by-term non-negative.");
    } else {
        println!("  ==> ALL terms non-negative! No cancellation needed.");
    }
}

// ============================================================
// Part 3: Print M_{n,j} table and look for patterns
// ============================================================

fn print_m_table(h: i128, n_max: i128, j_max: i128) {
    println!("\n=== M_{{n,j}} table for h={} ===", h);
    print!("{:>4}", "n\\j");
    for j in 0..=j_max {
        print!("{:>10}", j);
    }
    println!();
    for n in 0..=n_max {
        print!("{:>4}", n);
        for j in 0..=j_max {
            print!("{:>10}", m_entry(h, n, j));
        }
        println!();
    }
}

fn main() {
    println!("=============================================================");
    println!("  Investigating total positivity of M_{{n,j}} for (h,h,1,1)");
    println!("  M_{{n,j}} = sum_r C(h*r, j) * C(j+1, n-2j-2r)");
    println!("=============================================================");

    // ============================================================
    // Part 1: Transfer matrix TNN check
    // ============================================================
    println!("\n########################################");
    println!("# PART 1: Transfer Matrix TNN Check");
    println!("########################################");

    for h in 2..=3 {
        let n = 1 + 3 * h;
        println!("\n--- h={}, {} states ---", h, n);

        // Print state names
        for i in 0..n {
            println!("  State {}: {}", i, state_name(i, h));
        }

        // Natural ordering
        let natural: Vec<usize> = (0..n).collect();
        let mat = build_transfer_matrix(h, 1.0, &natural);

        println!("\n  Transfer matrix (natural ordering, t=1):");
        for i in 0..n {
            print!("    ");
            for j in 0..n {
                print!("{:5.1}", mat[i][j]);
            }
            println!("  <- {}", state_name(natural[i], h));
        }

        let (tnn, info) = check_tnn(&mat, n.min(5));
        println!(
            "  TNN (natural ordering): {} {}",
            tnn,
            if tnn { "" } else { &info }
        );

        // Try ordering: empty last
        let mut ordering_empty_last: Vec<usize> = (1..n).collect();
        ordering_empty_last.push(0);
        let mat2 = build_transfer_matrix(h, 1.0, &ordering_empty_last);
        let (tnn2, info2) = check_tnn(&mat2, n.min(5));
        println!(
            "  TNN (empty last): {} {}",
            tnn2,
            if tnn2 { "" } else { &info2 }
        );

        // Try ordering: group by d: all (v,1), then all (v,2), then all (v,3), then empty
        let mut ordering_by_d = Vec::new();
        for d in 1..=3 {
            for v in 1..=h {
                ordering_by_d.push(state_index(v, d));
            }
        }
        ordering_by_d.push(0);
        let mat3 = build_transfer_matrix(h, 1.0, &ordering_by_d);
        let (tnn3, info3) = check_tnn(&mat3, n.min(5));
        print!("  TNN (group by d, then empty): {} ", tnn3);
        if !tnn3 {
            println!("{}", info3);
        } else {
            println!();
        }

        // Try ordering: empty, then (1,1),(1,2),(1,3),(2,1),(2,2),(2,3),...
        // This is the natural ordering, already done.

        // Try ordering: empty first, then reverse v order
        let mut ordering_rev_v = vec![0usize];
        for v in (1..=h).rev() {
            for d in 1..=3 {
                ordering_rev_v.push(state_index(v, d));
            }
        }
        let mat4 = build_transfer_matrix(h, 1.0, &ordering_rev_v);
        let (tnn4, info4) = check_tnn(&mat4, n.min(5));
        print!("  TNN (empty, then rev v): {} ", tnn4);
        if !tnn4 {
            println!("{}", info4);
        } else {
            println!();
        }

        // Try: empty, (h,3),(h,2),(h,1),...,(1,3),(1,2),(1,1)
        let mut ordering_rev_all = vec![0usize];
        for v in (1..=h).rev() {
            for d in (1..=3).rev() {
                ordering_rev_all.push(state_index(v, d));
            }
        }
        let mat5 = build_transfer_matrix(h, 1.0, &ordering_rev_all);
        let (tnn5, info5) = check_tnn(&mat5, n.min(5));
        print!("  TNN (empty, rev v rev d): {} ", tnn5);
        if !tnn5 {
            println!("{}", info5);
        } else {
            println!();
        }

        // Try: (1,3),(1,2),(1,1),(2,3),(2,2),(2,1),..., empty
        let mut ordering_rev_d = Vec::new();
        for v in 1..=h {
            for d in (1..=3).rev() {
                ordering_rev_d.push(state_index(v, d));
            }
        }
        ordering_rev_d.push(0);
        let mat6 = build_transfer_matrix(h, 1.0, &ordering_rev_d);
        let (tnn6, info6) = check_tnn(&mat6, n.min(5));
        print!("  TNN (rev d within v, empty last): {} ", tnn6);
        if !tnn6 {
            println!("{}", info6);
        } else {
            println!();
        }

        // For h=2, try ALL permutations of 7 states (7! = 5040, feasible)
        if h == 2 {
            println!(
                "\n  Brute-force search over all {}! = {} orderings for h={}...",
                n,
                (1..=n).product::<usize>(),
                h
            );
            let mut best_min_minor = f64::NEG_INFINITY;
            let mut best_ordering: Vec<usize> = Vec::new();
            let mut tnn_count = 0;

            let mut perm: Vec<usize> = (0..n).collect();
            let total_perms: usize = (1..=n).product();
            let mut count = 0;
            loop {
                let mat_p = build_transfer_matrix(h, 1.0, &perm);
                let (is_tnn, _) = check_tnn(&mat_p, n);
                if is_tnn {
                    tnn_count += 1;
                    println!(
                        "    TNN ordering found: {:?} = [{}]",
                        perm,
                        perm.iter()
                            .map(|&i| state_name(i, h))
                            .collect::<Vec<_>>()
                            .join(", ")
                    );
                }

                // Track the ordering with best (least negative) min minor
                let min_minor = all_minors_min(&mat_p, n);
                if min_minor > best_min_minor {
                    best_min_minor = min_minor;
                    best_ordering = perm.clone();
                }

                count += 1;
                if !next_permutation(&mut perm) {
                    break;
                }
            }
            println!("  Checked {} orderings. TNN count: {}", count, tnn_count);
            if tnn_count == 0 {
                println!(
                    "  Best (closest to TNN) ordering: {:?} = [{}], min minor = {:.6}",
                    best_ordering,
                    best_ordering
                        .iter()
                        .map(|&i| state_name(i, h))
                        .collect::<Vec<_>>()
                        .join(", "),
                    best_min_minor
                );
            }
        }

        // Also try with symbolic t: check TNN for several t values
        println!("\n  Checking TNN for various t values (natural ordering):");
        for &t_val in &[0.5, 1.0, 2.0, 0.1] {
            let mat_t = build_transfer_matrix(h, t_val, &natural);
            let (tnn_t, info_t) = check_tnn(&mat_t, n.min(5));
            println!(
                "    t={}: TNN={} {}",
                t_val,
                tnn_t,
                if tnn_t { "" } else { &info_t }
            );
        }
    }

    // ============================================================
    // Part 2: Individual term analysis
    // ============================================================
    println!("\n########################################");
    println!("# PART 2: Individual Term Analysis in 2x2 Minors");
    println!("########################################");

    let h: i128 = 2;
    // Check several 2x2 minors
    check_minor_terms(h, 8, 12, 1, 2);
    check_minor_terms(h, 6, 10, 0, 1);
    check_minor_terms(h, 10, 14, 2, 3);
    check_minor_terms(h, 4, 8, 0, 2);
    check_minor_terms(h, 6, 12, 1, 3);

    // Check adjacent row/column minors (these are key for TP)
    println!("\n--- Systematic check of ADJACENT 2x2 minors for h=2 ---");
    let mut any_negative = false;
    for n1 in 0..=16i128 {
        for n2 in (n1 + 1)..=18 {
            for j1 in 0..=5i128 {
                for j2 in (j1 + 1)..=6 {
                    let minor = m_entry(h, n1, j1) * m_entry(h, n2, j2)
                        - m_entry(h, n1, j2) * m_entry(h, n2, j1);
                    if minor < 0 {
                        println!(
                            "  NEGATIVE minor! n1={}, n2={}, j1={}, j2={}: {}",
                            n1, n2, j1, j2, minor
                        );
                        any_negative = true;
                    }
                }
            }
        }
    }
    if !any_negative {
        println!("  All 2x2 minors non-negative for n in 0..18, j in 0..6 (h=2).");
    }

    // Also check for h=3
    println!("\n--- Systematic 2x2 minor check for h=3 ---");
    let h3: i128 = 3;
    any_negative = false;
    for n1 in 0..=20i128 {
        for n2 in (n1 + 1)..=22 {
            for j1 in 0..=5i128 {
                for j2 in (j1 + 1)..=6 {
                    let minor = m_entry(h3, n1, j1) * m_entry(h3, n2, j2)
                        - m_entry(h3, n1, j2) * m_entry(h3, n2, j1);
                    if minor < 0 {
                        println!(
                            "  NEGATIVE minor! n1={}, n2={}, j1={}, j2={}: {}",
                            n1, n2, j1, j2, minor
                        );
                        any_negative = true;
                    }
                }
            }
        }
    }
    if !any_negative {
        println!("  All 2x2 minors non-negative for n in 0..22, j in 0..6 (h=3).");
    }

    // ============================================================
    // Part 3: Print table and look for patterns
    // ============================================================
    println!("\n########################################");
    println!("# PART 3: M_{{n,j}} Table and Pattern Search");
    println!("########################################");

    print_m_table(2, 20, 6);
    print_m_table(3, 20, 6);

    // Look for closed forms for h=2
    println!("\n--- Pattern search for h=2 ---");
    println!("Checking if M_{{n,j}} for h=2 matches simple binomial expressions:");
    let h2: i128 = 2;
    println!("\n  j=0: M_{{n,0}} = sum_r C(0,0)*C(1, n-2r) = sum_r C(1, n-2r)");
    println!("  Since C(1,k) = k for k=0,1 and 0 otherwise:");
    for n in 0..=10 {
        let val = m_entry(h2, n, 0);
        // C(1, n-2r) is nonzero only when n-2r in {0,1}
        // r = n/2 gives C(1,0)=1; r=(n-1)/2 gives C(1,1)=1 if n odd
        println!("    M_{{{},0}} = {}", n, val);
    }

    println!("\n  j=1: M_{{n,1}} = sum_r C(2r, 1)*C(2, n-2-2r) = sum_r 2r * C(2, n-2-2r)");
    for n in 0..=14 {
        let val = m_entry(h2, n, 1);
        // Try to match: maybe floor((n-1)^2/8) or similar?
        let candidate1 = if n >= 2 { (n - 1) * (n - 2) / 4 } else { 0 };
        // Actually let's check: C(2, n-2-2r) is nonzero for n-2-2r in {0,1,2}
        // i.e., r = (n-2)/2, (n-3)/2, (n-4)/2
        println!(
            "    M_{{{},1}} = {}, guess (n-1)(n-2)/4 = {}",
            n, val, candidate1
        );
    }

    println!("\n  Trying to express M_{{n,j}} as single binomial for h=2:");
    for j in 0..=4i128 {
        println!("  j={}:", j);
        for n in 0..=16i128 {
            let val = m_entry(h2, n, j);
            if val > 0 {
                // Check C(n-j, j), C(n-j-1, j), C(n-j, 2j), etc.
                let guess_a = binom(n - j, j);
                let guess_b = binom(n - j - 1, j);
                let guess_c = binom(n - j, 2 * j);
                let guess_d = binom(n, j + 1) - binom(n - 2, j + 1);
                // More: C(n-j, j+1)
                let guess_e = binom(n - j, j + 1);
                let guess_f = binom(n - j + 1, j + 1);

                let mut matches = Vec::new();
                if val == guess_a {
                    matches.push(format!("C({},{})", n - j, j));
                }
                if val == guess_b {
                    matches.push(format!("C({},{})", n - j - 1, j));
                }
                if val == guess_c {
                    matches.push(format!("C({},{})", n - j, 2 * j));
                }
                if val == guess_d {
                    matches.push(format!("C({},{}) - C({},{})", n, j + 1, n - 2, j + 1));
                }
                if val == guess_e {
                    matches.push(format!("C({},{})", n - j, j + 1));
                }
                if val == guess_f {
                    matches.push(format!("C({},{})", n - j + 1, j + 1));
                }

                if !matches.is_empty() {
                    println!("    n={}: M={}, matches: {}", n, val, matches.join(" or "));
                } else {
                    println!("    n={}: M={}, no simple match found", n, val);
                }
            }
        }
    }

    // Check 3x3 and 4x4 minors for h=2
    println!("\n--- Checking 3x3 minors for h=2 ---");
    any_negative = false;
    let mut checked_3 = 0u64;
    'outer3: for n1 in 0..=14i128 {
        for n2 in (n1 + 1)..=16 {
            for n3 in (n2 + 1)..=18 {
                for j1 in 0..=4i128 {
                    for j2 in (j1 + 1)..=5 {
                        for j3 in (j2 + 1)..=6 {
                            // 3x3 determinant
                            let det = det3x3(h2, n1, n2, n3, j1, j2, j3);
                            checked_3 += 1;
                            if det < 0 {
                                println!(
                                    "  NEGATIVE 3x3 minor! n=({},{},{}), j=({},{},{}): {}",
                                    n1, n2, n3, j1, j2, j3, det
                                );
                                any_negative = true;
                                if checked_3 > 100000 {
                                    break 'outer3;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    println!(
        "  Checked {} 3x3 minors for h=2. Any negative: {}",
        checked_3, any_negative
    );
}

fn det3x3(h: i128, n1: i128, n2: i128, n3: i128, j1: i128, j2: i128, j3: i128) -> i128 {
    let a = m_entry(h, n1, j1);
    let b = m_entry(h, n1, j2);
    let c = m_entry(h, n1, j3);
    let d = m_entry(h, n2, j1);
    let e = m_entry(h, n2, j2);
    let f = m_entry(h, n2, j3);
    let g = m_entry(h, n3, j1);
    let hh = m_entry(h, n3, j2);
    let i = m_entry(h, n3, j3);
    a * (e * i - f * hh) - b * (d * i - f * g) + c * (d * hh - e * g)
}

fn all_minors_min(mat: &[Vec<f64>], n: usize) -> f64 {
    let mut min_val = f64::INFINITY;
    for size in 1..=n {
        let row_subsets = subsets(n, size);
        let col_subsets = subsets(n, size);
        for rows in &row_subsets {
            for cols in &col_subsets {
                let det = minor_det(mat, rows, cols);
                if det < min_val {
                    min_val = det;
                }
            }
        }
    }
    min_val
}
