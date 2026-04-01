// Test: does Σ_{i<j} U_i ≪ Σ_i W_i hold for all j and all boards?
// This is the KEY condition for forward AW preservation.
// Also test: does Σ_{i<j} U_i ≪ Σ_{i≥j} W_i + t·Σ_{i<j} D_i hold? (= R decomposition)
use combpoly::order::bruhat_lower_ideal;
use experiments::peak_utils::{
    board_to_perm, gen_boards, is_312_avoiding, pa, peak_count, pmt, pt, pz,
};
use polynomial_tools::real_rootedness::check_weak_interlacing;
fn interlaces(f: &[i64], g: &[i64]) -> bool {
    let f = pt(f);
    let g = pt(g);
    if pz(&f) {
        return true;
    }
    if pz(&g) {
        return false;
    }
    check_weak_interlacing(&f, &g).unwrap_or(false)
}
fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);
    let mut sum_u_sum_w = [0u64; 2]; // Σ_{i<j} U_i ≪ Σ_i W_i
    let mut sum_u_r = [0u64; 2]; // Σ_{i<j} U_i ≪ R = Σ_{i≥j} W_i + t·Σ_{i<j} D_i
    let mut each_u_sum_w = [0u64; 2]; // each U_i ≪ Σ W_i (for i < j)
    let mut sum_u_each_w = [0u64; 2]; // Σ U_i ≪ each W_{i'} (for i' ≥ 1)
                                      // Also: does Σ_{i<j} U_i ≪ tΣ_{i<j} D_i hold? (sub-condition)
    let mut sum_u_t_sum_d = [0u64; 2];
    for n in 1..=max_n {
        for board in &gen_boards(n) {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            let m = board[0] as usize;
            if n <= 1 {
                continue;
            }
            let ideal = bruhat_lower_ideal(&perm);
            let mut d_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            let mut u_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            for pi in &ideal {
                let j = pi[0] as usize;
                if j > m {
                    continue;
                }
                let l = *pi.last().unwrap() as usize;
                if l > m {
                    continue;
                }
                let pk = peak_count(pi);
                let poly = if pi.len() >= 2 && pi[0] > pi[1] {
                    &mut d_jl[j][l]
                } else {
                    &mut u_jl[j][l]
                };
                while poly.len() <= pk {
                    poly.push(0);
                }
                poly[pk] += 1;
            }
            let mut w_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            for j in 1..=m {
                for l in 1..=m {
                    w_jl[j][l] = pa(&pmt(&d_jl[j][l]), &u_jl[j][l]);
                }
            }
            for l in 1..=m {
                // Compute Σ W_i for all i
                let mut sum_w = vec![0i64];
                for i in 1..=m {
                    sum_w = pa(&sum_w, &w_jl[i][l]);
                }
                for j in 2..=m {
                    // j is the cutoff: test Σ_{i<j} U_i
                    let mut sum_u = vec![0i64];
                    for i in 1..j {
                        sum_u = pa(&sum_u, &u_jl[i][l]);
                    }
                    if pz(&sum_u) {
                        continue;
                    }
                    // Test Σ U ≪ Σ W
                    if !pz(&sum_w) {
                        sum_u_sum_w[0] += 1;
                        if !interlaces(&sum_u, &sum_w) {
                            sum_u_sum_w[1] += 1;
                        }
                    }
                    // Test Σ U ≪ R = Σ_{i≥j} W_i + t·Σ_{i<j} D_i
                    let mut r = vec![0i64];
                    for i in j..=m {
                        r = pa(&r, &w_jl[i][l]);
                    }
                    let mut td = vec![0i64];
                    for i in 1..j {
                        td = pa(&td, &pmt(&d_jl[i][l]));
                    }
                    r = pa(&r, &td);
                    if !pz(&r) {
                        sum_u_r[0] += 1;
                        if !interlaces(&sum_u, &r) {
                            sum_u_r[1] += 1;
                        }
                    }
                    // Test Σ U ≪ tΣ D (sub-condition)
                    if !pz(&td) {
                        sum_u_t_sum_d[0] += 1;
                        if !interlaces(&sum_u, &td) {
                            sum_u_t_sum_d[1] += 1;
                        }
                    }
                    // Test each U_i ≪ Σ W
                    for i in 1..j {
                        if !pz(&u_jl[i][l]) && !pz(&sum_w) {
                            each_u_sum_w[0] += 1;
                            if !interlaces(&u_jl[i][l], &sum_w) {
                                each_u_sum_w[1] += 1;
                            }
                        }
                    }
                }
                // Test Σ_{i<j} U_i ≪ each W_{i'}
                for j in 2..=m {
                    let mut sum_u = vec![0i64];
                    for i in 1..j {
                        sum_u = pa(&sum_u, &u_jl[i][l]);
                    }
                    if pz(&sum_u) {
                        continue;
                    }
                    for ip in 1..=m {
                        if !pz(&w_jl[ip][l]) {
                            sum_u_each_w[0] += 1;
                            if !interlaces(&sum_u, &w_jl[ip][l]) {
                                sum_u_each_w[1] += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    println!(
        "Σ U_{{i<j}} ≪ Σ W_i (full sum):     {}/{} pass ({} FAIL)",
        sum_u_sum_w[0] - sum_u_sum_w[1],
        sum_u_sum_w[0],
        sum_u_sum_w[1]
    );
    println!(
        "Σ U_{{i<j}} ≪ R (= Σ W_{{i≥j}} + tΣ D_{{i<j}}): {}/{} pass ({} FAIL)",
        sum_u_r[0] - sum_u_r[1],
        sum_u_r[0],
        sum_u_r[1]
    );
    println!(
        "Σ U_{{i<j}} ≪ tΣ D_{{i<j}}:           {}/{} pass ({} FAIL)",
        sum_u_t_sum_d[0] - sum_u_t_sum_d[1],
        sum_u_t_sum_d[0],
        sum_u_t_sum_d[1]
    );
    println!(
        "each U_i ≪ Σ W (individual):    {}/{} pass ({} FAIL)",
        each_u_sum_w[0] - each_u_sum_w[1],
        each_u_sum_w[0],
        each_u_sum_w[1]
    );
    println!(
        "Σ U_{{i<j}} ≪ each W_{{i'}}:           {}/{} pass ({} FAIL)",
        sum_u_each_w[0] - sum_u_each_w[1],
        sum_u_each_w[0],
        sum_u_each_w[1]
    );
}
