// Quick test: does reversed UU hold? U_{j,l} << U_{j',l} for j > j'?
// And: does each U_i << U_1 (= A_1) for all i?
use combpoly::order::bruhat_lower_ideal;
use experiments::peak_utils::{
    board_to_perm, gen_boards, is_312_avoiding, peak_count, pmt, pt, pz,
};
use polynomial_tools::real_rootedness::check_weak_interlacing;
fn pdeg(p: &[i64]) -> Option<usize> {
    let v = pt(p);
    if pz(&v) {
        None
    } else {
        Some(v.len() - 1)
    }
}
fn interlaces(f: &[i64], g: &[i64]) -> bool {
    let f = pt(f);
    let g = pt(g);
    if pz(&f) {
        return true;
    }
    if pz(&g) {
        return false;
    }
    match check_weak_interlacing(&f, &g) {
        Some(true) => true,
        Some(false) => false,
        None => match (pdeg(&f), pdeg(&g)) {
            (Some(df), Some(dg)) if df == dg => {
                let tf = pmt(&f);
                check_weak_interlacing(&g, &tf).unwrap_or(false)
            }
            _ => false,
        },
    }
}

fn main() {
    let max_n: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);
    let mut rev_uu = [0u64; 2]; // U_j << U_j' for j > j'
    let mut fwd_uu = [0u64; 2]; // U_j << U_j' for j < j'
    let mut u_to_u1 = [0u64; 2]; // U_j << U_1 for j > 1
    let mut all_uu = [0u64; 2]; // U_j << U_j' for all j != j'
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
                if pi.len() >= 2 && pi[0] < pi[1] {
                    while u_jl[j][l].len() <= pk {
                        u_jl[j][l].push(0);
                    }
                    u_jl[j][l][pk] += 1;
                }
            }
            for l in 1..=m {
                for j in 1..=m {
                    for jp in 1..=m {
                        if j == jp {
                            continue;
                        }
                        if !pz(&u_jl[j][l]) && !pz(&u_jl[jp][l]) {
                            all_uu[0] += 1;
                            let ok = interlaces(&u_jl[j][l], &u_jl[jp][l]);
                            if !ok {
                                all_uu[1] += 1;
                            }
                            if j > jp {
                                rev_uu[0] += 1;
                                if !ok {
                                    rev_uu[1] += 1;
                                }
                            }
                            if j < jp {
                                fwd_uu[0] += 1;
                                if !ok {
                                    fwd_uu[1] += 1;
                                }
                            }
                            if jp == 1 && j > 1 {
                                u_to_u1[0] += 1;
                                if !ok {
                                    u_to_u1[1] += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    println!(
        "Forward UU (j < j'): {}/{} pass ({} FAIL)",
        fwd_uu[0] - fwd_uu[1],
        fwd_uu[0],
        fwd_uu[1]
    );
    println!(
        "Reversed UU (j > j'): {}/{} pass ({} FAIL)",
        rev_uu[0] - rev_uu[1],
        rev_uu[0],
        rev_uu[1]
    );
    println!(
        "U_j << U_1 (j > 1): {}/{} pass ({} FAIL)",
        u_to_u1[0] - u_to_u1[1],
        u_to_u1[0],
        u_to_u1[1]
    );
    println!(
        "All UU (j != j'): {}/{} pass ({} FAIL)",
        all_uu[0] - all_uu[1],
        all_uu[0],
        all_uu[1]
    );
}
