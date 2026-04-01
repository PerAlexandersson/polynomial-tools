use combpoly::order::bruhat_lower_ideal;
use experiments::peak_utils::{
    board_to_perm, gen_boards, is_312_avoiding, pa, peak_count, pmt, pt, pz,
};
use polynomial_tools::real_rootedness::{check_weak_interlacing, format_poly};
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
    // Find the 1 failure of D_k^+(col1) ≪ W_{i,col2}
    for n in 1..=8usize {
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
            let mp = m + 1;
            for l in 2..=m {
                let col1 = l;
                let col2 = l - 1;
                for kp in (l + 1)..=mp {
                    let mut dk = vec![0i64];
                    for j in 1..kp.min(m + 1) {
                        dk = pa(&dk, &a_jl[j][col1]);
                    }
                    for i in 1..=m {
                        if !pz(&dk) && !pz(&w_jl[i][col2]) && !interlaces(&dk, &w_jl[i][col2]) {
                            println!(
                                "FAIL: board={:?} col1={} col2={} k'={} i={}",
                                board, col1, col2, kp, i
                            );
                            println!("  D_k^+ = {}", format_poly(&pt(&dk)));
                            println!("  W_i   = {}", format_poly(&pt(&w_jl[i][col2])));
                            println!("  U_1(col1) = {}", format_poly(&pt(&u_jl[1][col1])));
                            // Check transitivity chain
                            println!("  D_k^+ ≪ U_1(col1)? {}", interlaces(&dk, &u_jl[1][col1]));
                            println!(
                                "  U_1(col1) ≪ W_i(col2)? {}",
                                interlaces(&u_jl[1][col1], &w_jl[i][col2])
                            );
                        }
                    }
                }
            }
        }
    }
    println!("Done checking.");
}
