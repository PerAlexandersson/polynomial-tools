use combpoly::parking::all_parking_functions;
use polynomial_tools::real_rootedness as polynomial;
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};

fn des_word(w: &[u8]) -> usize {
    let mut d = 0;
    for i in 0..w.len().saturating_sub(1) {
        if w[i] > w[i + 1] { d += 1; }
    }
    d
}

fn main() {
    let max_n: u8 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(9);
    let mut polys: Vec<Vec<i64>> = Vec::new();
    for n in 1u8..=max_n {
        let pfs = all_parking_functions(n);
        let max_des = (n as usize).saturating_sub(1);
        let mut coeffs = vec![0i64; max_des + 1];
        for pf in &pfs { let d = des_word(pf); if d < coeffs.len() { coeffs[d] += 1; } }
        while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 { coeffs.pop(); }
        polys.push(coeffs);
    }

    // Try various recurrence shapes
    let configs = [
        ("rec_len=1, diff=1", AdaptiveSearchOptions { max_rec_len: 1, max_var_deg: 3, max_idx_deg: 3, max_diff_deg: 1, verbose: false, ..Default::default() }),
        ("rec_len=1, diff=2", AdaptiveSearchOptions { max_rec_len: 1, max_var_deg: 3, max_idx_deg: 3, max_diff_deg: 2, verbose: false, ..Default::default() }),
        ("rec_len=2, diff=0", AdaptiveSearchOptions { max_rec_len: 2, max_var_deg: 2, max_idx_deg: 3, max_diff_deg: 0, verbose: false, ..Default::default() }),
        ("rec_len=2, diff=1", AdaptiveSearchOptions { max_rec_len: 2, max_var_deg: 2, max_idx_deg: 2, max_diff_deg: 1, verbose: false, ..Default::default() }),
        ("rec_len=2, diff=2", AdaptiveSearchOptions { max_rec_len: 2, max_var_deg: 2, max_idx_deg: 2, max_diff_deg: 2, verbose: false, ..Default::default() }),
    ];

    for (label, opts) in &configs {
        eprintln!("Trying {}...", label);
        match find_recurrence_adaptive(&polys, opts) {
            Some(r) => println!("{}: {}", label, r.recurrence),
            None => println!("{}: not found", label),
        }
    }
}
