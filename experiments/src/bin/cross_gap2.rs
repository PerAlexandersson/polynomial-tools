//! Test: U_{1,l'} ≪ W_{i,l'-1} (cross-column U1≪W) and
//! D_{k'}^+(col l') ≪ W_{i,l'-1} (cross-column D+≪W).
//! If BOTH hold → cross-group DU by right-cone + transitivity.
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
    // Strategy: D_k^+(col1) ≪ U_1(col1) ≪ W_i(col2) → D_k^+ ≪ W_i → D_k^+ ≪ U_j^+(col2)
    let mut u1_w_cross = [0u64; 2]; // U_{1,l'} ≪ W_{i,l'-1} for all i
    let mut dk_w_cross = [0u64; 2]; // D_k^+(col1) ≪ W_{i,col2} for all cross-group pairs
    let mut dk_u1_same = [0u64; 2]; // D_k^+(col1) ≪ U_{1,col1} (= A_{1,col1})
                                    // Also: A_{1,l'} ≪ A_{i,l'-1} for all i (cross-col A1 vs A)
    let mut a1_a_cross = [0u64; 2];
    // A_{1,l'} ≪ W_{i,l'-1} for all i
    let mut a1_w_cross = [0u64; 2];
    // Also test the OTHER cross direction: D_j^+(col2) ≪ U_1(col2) ≪ W_i(col1)?
    // For the other cross-group: k' < l', j' > l'. D uses col l'-1, U uses col l'.
    let mut u1_w_cross2 = [0u64; 2]; // U_{1,l'-1} ≪ W_{i,l'} for all i
    for n in 1..=max_n {
        for board in &gen_boards(n) {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) {
                continue;
            }
            let m = board[0] as usize;
            if n <= 2 {
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
            let mut a_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            let mut w_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            for j in 1..=m {
                for l in 1..=m {
                    a_jl[j][l] = pa(&d_jl[j][l], &u_jl[j][l]);
                    w_jl[j][l] = pa(&pmt(&d_jl[j][l]), &u_jl[j][l]);
                }
            }
            // Test cross-column conditions between adjacent columns
            for l in 2..=m {
                let col1 = l; // "upper" column
                let col2 = l - 1; // "lower" column
                                  // U_{1,col1} ≪ W_{i,col2} for all i
                for i in 1..=m {
                    if !pz(&u_jl[1][col1]) && !pz(&w_jl[i][col2]) {
                        u1_w_cross[0] += 1;
                        if !interlaces(&u_jl[1][col1], &w_jl[i][col2]) {
                            u1_w_cross[1] += 1;
                        }
                    }
                    // A_{1,col1} ≪ A_{i,col2}
                    if !pz(&a_jl[1][col1]) && !pz(&a_jl[i][col2]) {
                        a1_a_cross[0] += 1;
                        if !interlaces(&a_jl[1][col1], &a_jl[i][col2]) {
                            a1_a_cross[1] += 1;
                        }
                    }
                    // A_{1,col1} ≪ W_{i,col2}
                    if !pz(&a_jl[1][col1]) && !pz(&w_jl[i][col2]) {
                        a1_w_cross[0] += 1;
                        if !interlaces(&a_jl[1][col1], &w_jl[i][col2]) {
                            a1_w_cross[1] += 1;
                        }
                    }
                    // U_{1,col2} ≪ W_{i,col1}
                    if !pz(&u_jl[1][col2]) && !pz(&w_jl[i][col1]) {
                        u1_w_cross2[0] += 1;
                        if !interlaces(&u_jl[1][col2], &w_jl[i][col1]) {
                            u1_w_cross2[1] += 1;
                        }
                    }
                }
                // D_k^+(col1) ≪ W_{i,col2} for cross-group pairs
                let mp = m + 1;
                for kp in (l + 1)..=mp {
                    // D_k^+(col1) = Σ_{j<kp} A_{j,col1}
                    let mut dk = vec![0i64];
                    for j in 1..kp.min(m + 1) {
                        dk = pa(&dk, &a_jl[j][col1]);
                    }
                    for i in 1..=m {
                        if !pz(&dk) && !pz(&w_jl[i][col2]) {
                            dk_w_cross[0] += 1;
                            if !interlaces(&dk, &w_jl[i][col2]) {
                                dk_w_cross[1] += 1;
                            }
                        }
                    }
                    // D_k^+(col1) ≪ U_{1,col1}
                    if !pz(&dk) && !pz(&u_jl[1][col1]) {
                        dk_u1_same[0] += 1;
                        if !interlaces(&dk, &u_jl[1][col1]) {
                            dk_u1_same[1] += 1;
                        }
                    }
                }
            }
        }
    }
    println!("=== Cross-group strategy: transitivity via U_1 ===");
    let show = |name: &str, c: [u64; 2]| {
        if c[0] == 0 {
            println!("  {}: (no data)", name);
        } else if c[1] == 0 {
            println!("  {}: {}/{} ALL PASS <<<", name, c[0], c[0]);
        } else {
            println!("  {}: {}/{} pass ({} FAIL)", name, c[0] - c[1], c[0], c[1]);
        }
    };
    show("D_k^+(col1) ≪ U_{1,col1}", dk_u1_same);
    show("U_{1,l'} ≪ W_{i,l'-1}", u1_w_cross);
    show("A_{1,l'} ≪ W_{i,l'-1}", a1_w_cross);
    show("A_{1,l'} ≪ A_{i,l'-1}", a1_a_cross);
    show("D_k^+(col1) ≪ W_{i,col2}", dk_w_cross);
    show("U_{1,l'-1} ≪ W_{i,l'} (other dir)", u1_w_cross2);
}
