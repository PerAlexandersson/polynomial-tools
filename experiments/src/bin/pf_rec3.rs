use combpoly::parking::is_parking_function;
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};

fn des_word(w: &[u8]) -> usize {
    (0..w.len().saturating_sub(1)).filter(|&i| w[i] > w[i + 1]).count()
}

fn pf_max_entry(n: u8, k: u8) -> Vec<Vec<u8>> {
    let mut result = Vec::new();
    let mut current = Vec::with_capacity(n as usize);
    gen(n, k, &mut current, &mut result);
    result
}

fn gen(n: u8, k: u8, cur: &mut Vec<u8>, res: &mut Vec<Vec<u8>>) {
    if cur.len() == n as usize {
        if is_parking_function(cur) { res.push(cur.clone()); }
        return;
    }
    for v in 1..=k { cur.push(v); gen(n, k, cur, res); cur.pop(); }
}

fn main() {
    for max_val in [3u8, 4] {
        let max_n = if max_val == 3 { 13u8 } else { 10 };
        eprintln!("Computing max_entry={} up to n={}...", max_val, max_n);
        let mut polys = Vec::new();
        for n in 1..=max_n {
            let pfs = pf_max_entry(n, max_val);
            let mut coeffs = vec![0i64; n as usize];
            for pf in &pfs { let d = des_word(pf); if d < coeffs.len() { coeffs[d] += 1; } }
            while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 { coeffs.pop(); }
            eprintln!("  n={}: |PF|={}", n, pfs.len());
            polys.push(coeffs);
        }

        println!("\n=== max_entry={} recurrence search ===", max_val);
        for &(label, rl, vd, id, dd) in &[
            ("r1d0", 1, 3, 3, 0), ("r1d1", 1, 3, 3, 1),
            ("r2d0", 2, 3, 3, 0), ("r2d1", 2, 3, 3, 1),
            ("r3d0", 3, 2, 2, 0), ("r3d1", 3, 2, 2, 1),
        ] {
            let opts = AdaptiveSearchOptions {
                max_rec_len: rl, max_var_deg: vd, max_idx_deg: id, max_diff_deg: dd,
                verbose: false, ..Default::default()
            };
            match find_recurrence_adaptive(&polys, &opts) {
                Some(r) => { println!("  {}: {}", label, r.recurrence); break; }
                None => println!("  {}: not found", label),
            }
        }
    }
}
