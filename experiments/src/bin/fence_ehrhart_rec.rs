use combinatoric_core::poset::Poset;

fn binom(n: i64, k: i64) -> i64 {
    if k < 0 || k > n {
        return 0;
    }
    let mut r = 1i64;
    for i in 0..k {
        r = r * (n - i) / (i + 1);
    }
    r
}

fn a109466(n: i64, k: i64) -> i64 {
    let sign = if (n - k) % 2 == 0 { 1 } else { -1 };
    sign * binom(k, n - k)
}

fn main() {
    // Claim: the GF for the char polys is
    //   sum_{s>=0} p_s(x) y^s = (1+y) / (1 + (2-x)y + y^2)
    //
    // Recurrence: p_s(x) + (2-x)*p_{s-1}(x) + p_{s-2}(x) = 0 for s >= 2,
    // with p_0 = 1, p_1 = x - 1.

    println!("=== Verifying GF = (1+y) / (1 + (2-x)y + y^2) ===\n");

    let max_s = 12;
    let max_deg = max_s + 1;

    // Compute p_s via recurrence: p_s = -(2-x)*p_{s-1} - p_{s-2}
    // = (x-2)*p_{s-1} - p_{s-2}
    let mut q: Vec<Vec<i64>> = Vec::new();
    q.push(vec![1]); // p_0 = 1
    q.push(vec![-1, 1]); // p_1 = -1 + x

    for s in 2..=max_s {
        let mut ps = vec![0i64; max_deg];
        // (x-2)*p_{s-1}: -2*p_{s-1}[j] + p_{s-1}[j-1]
        for j in 0..max_deg {
            if j < q[s - 1].len() {
                ps[j] += -2 * q[s - 1][j];
            }
            if j > 0 && j - 1 < q[s - 1].len() {
                ps[j] += q[s - 1][j - 1];
            }
            // -p_{s-2}
            if j < q[s - 2].len() {
                ps[j] -= q[s - 2][j];
            }
        }
        while ps.last() == Some(&0) {
            ps.pop();
        }
        q.push(ps);
    }

    // Compare with A109466
    for s in 0..=max_s as i64 {
        let mut a_coeffs = Vec::new();
        for j in 0..=s {
            a_coeffs.push(a109466(2 * s, s + j));
        }

        let match_ok = a_coeffs == q[s as usize];
        println!("s={:2}: A109466={:?}", s, a_coeffs);
        println!("       GF_rec ={:?}  match={}", q[s as usize], match_ok);
    }

    // Cross-check with actual fence data
    println!("\n=== Cross-check with fence Ehrhart evaluations ===\n");
    let max_n = 25;
    let mut all_evals: Vec<Vec<i64>> = Vec::new();
    for n in 1..=max_n {
        let p = Poset::fence(n);
        let mut row = Vec::new();
        for s in 0..=max_s as usize {
            row.push(p.count_weak_order_preserving_dp(s) as i64);
        }
        all_evals.push(row);
    }

    for s in 2..=max_s {
        let odd_col: Vec<i64> = all_evals
            .iter()
            .enumerate()
            .filter(|(i, _)| i % 2 == 0)
            .map(|(_, r)| r[s])
            .collect();

        let char_poly = &q[s];
        let order = char_poly.len() - 1;
        let mut ok = true;
        for m in order..odd_col.len() {
            let mut sum = 0i64;
            for (j, &c) in char_poly.iter().enumerate() {
                sum += c * odd_col[m - j];
            }
            if sum != 0 {
                ok = false;
                break;
            }
        }
        println!("s={:2}: odd-n recurrence verified={}", s, ok);
    }

    println!("\n=== Summary ===");
    println!("The characteristic polynomials p_s(x) for the odd/even");
    println!("fence subsequence recurrences satisfy:");
    println!("  sum_{{s>=0}} p_s(x) y^s = (1+y) / (1 + (2-x)y + y^2)");
    println!("Equivalently: p_s + (2-x) p_{{s-1}} + p_{{s-2}} = 0 for s >= 2.");
}
