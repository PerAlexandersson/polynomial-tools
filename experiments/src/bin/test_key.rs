// The KEY identity: x*D_k^+ + y*U_l^+ = Σ b_j A_j + y(t-1)·Σ_{j≥l} D_j
// where b_j = x*[j<k] + y*[j≥l].
// So DU at λ^+ iff for all k,l: Σ_{j≥l} D_j ≪ Σ b_j A_j and boundary holds.
// Test: does Σ_{j≥l} D_j ≪ Σ b_j A_j hold for ALL k, l, x=1, y=1?
// (This is the shift lemma condition for DU.)
use combpoly::order::bruhat_lower_ideal;
use experiments::peak_utils::{board_to_perm, gen_boards, is_312_avoiding, pa, peak_count, pt, pz};
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
    // Key identity: D_k^+ + U_l^+ = Σ b_j A_j + (t-1)·Σ_{j≥l} D_j
    // where b_j = [j<k] + [j≥l].
    // Shift lemma: D_k^+ ≪ D_k^+ + U_l^+ (i.e., DU) iff
    //   (1) Σ_{j≥l} D_j ≪ D_k^+ + U_l^+ and (2) boundary holds.
    // But we can also write: D_k^+ + U_l^+ = P + (t-1)h where
    //   P = D_k^+ + U_l^+ - (t-1)h = Σ b_j A_j and h = Σ_{j≥l} D_j.
    // Actually: for the DU condition D_k^+ ≪ U_l^+, we consider
    //   U_l^+ = D_k^+ + [(t-1) partial D terms + rest].
    // Let me instead directly verify: h = Σ_{j≥l} D_j ≪ Σ b_j A_j
    // where b_j = [j<k] + [j≥l].
    let mut h_ba = [0u64; 2]; // h ≪ Σ b_j A_j for all k, l
    let mut h_ba_overlap = [0u64; 2]; // same but only k > l (overlapping case)
                                      // Also test D_j ≪ Σ b_j A_j for each individual D_j (j ≥ l)
    let mut each_d_ba = [0u64; 2];
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
            let mut a_jl: Vec<Vec<Vec<i64>>> = vec![vec![vec![0i64]; m + 1]; m + 1];
            for j in 1..=m {
                for l in 1..=m {
                    a_jl[j][l] = pa(&d_jl[j][l], &u_jl[j][l]);
                }
            }
            for col in 1..=m {
                for k in 1..=m + 1 {
                    // k ranges 1..=m+1 for the output
                    for l in 1..=m + 1 {
                        if k == l {
                            continue;
                        }
                        // b_j = [j<k] + [j≥l] for j in 1..=m
                        let mut ba = vec![0i64]; // Σ b_j A_j
                        for j in 1..=m {
                            let b = (if j < k { 1 } else { 0 }) + (if j >= l { 1 } else { 0 });
                            if b > 0 {
                                let scaled: Vec<i64> =
                                    a_jl[j][col].iter().map(|&c| c * b as i64).collect();
                                ba = pa(&ba, &scaled);
                            }
                        }
                        // h = Σ_{j≥l} D_j
                        let mut h = vec![0i64];
                        for j in l.max(1)..=m {
                            h = pa(&h, &d_jl[j][col]);
                        }
                        if pz(&h) || pz(&ba) {
                            continue;
                        }
                        h_ba[0] += 1;
                        if !interlaces(&h, &ba) {
                            h_ba[1] += 1;
                        }
                        if k > l {
                            h_ba_overlap[0] += 1;
                            if !interlaces(&h, &ba) {
                                h_ba_overlap[1] += 1;
                            }
                        }
                        // Each D_j (j >= l) ≪ Σ b_j A_j
                        for j in l.max(1)..=m {
                            if !pz(&d_jl[j][col]) && !pz(&ba) {
                                each_d_ba[0] += 1;
                                if !interlaces(&d_jl[j][col], &ba) {
                                    each_d_ba[1] += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    println!(
        "h = Σ D_{{j≥l}} ≪ Σ b_j A_j (ALL k,l):    {}/{} pass ({} FAIL)",
        h_ba[0] - h_ba[1],
        h_ba[0],
        h_ba[1]
    );
    println!(
        "h ≪ Σ b_j A_j (overlap k>l only):   {}/{} pass ({} FAIL)",
        h_ba_overlap[0] - h_ba_overlap[1],
        h_ba_overlap[0],
        h_ba_overlap[1]
    );
    println!(
        "each D_j ≪ Σ b_j A_j:               {}/{} pass ({} FAIL)",
        each_d_ba[0] - each_d_ba[1],
        each_d_ba[0],
        each_d_ba[1]
    );
}
