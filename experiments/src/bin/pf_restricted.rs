use combpoly::parking::is_parking_function;
use polynomial_tools::real_rootedness as polynomial;
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};

fn des_word(w: &[u8]) -> usize {
    let mut d = 0;
    for i in 0..w.len().saturating_sub(1) {
        if w[i] > w[i + 1] {
            d += 1;
        }
    }
    d
}

/// Generate all words on {1,...,k}^n that are parking functions
fn pf_max_entry(n: u8, k: u8) -> Vec<Vec<u8>> {
    let mut result = Vec::new();
    let mut current = Vec::with_capacity(n as usize);
    gen_words(n, k, &mut current, &mut result);
    result
}

fn gen_words(n: u8, k: u8, current: &mut Vec<u8>, result: &mut Vec<Vec<u8>>) {
    if current.len() == n as usize {
        if is_parking_function(current) {
            result.push(current.clone());
        }
        return;
    }
    for v in 1..=k {
        current.push(v);
        gen_words(n, k, current, result);
        current.pop();
    }
}

fn run_for_max_entry(k: u8, max_n: u8) {
    println!("\n=== PF with max entry ≤ {}, descent polynomial ===\n", k);
    let mut polys: Vec<Vec<i64>> = Vec::new();

    for n in 1..=max_n {
        let pfs = pf_max_entry(n, k);
        let max_des = (n as usize).saturating_sub(1);
        let mut coeffs = vec![0i64; max_des + 1];
        for pf in &pfs {
            let d = des_word(pf);
            if d < coeffs.len() {
                coeffs[d] += 1;
            }
        }
        while coeffs.len() > 1 && *coeffs.last().unwrap() == 0 {
            coeffs.pop();
        }
        let rr = polynomial::is_real_rooted(&coeffs);
        print!("n={}: [", n);
        for (i, &c) in coeffs.iter().enumerate() {
            if i > 0 {
                print!(", ");
            }
            print!("{}", c);
        }
        print!("]  sum={}", coeffs.iter().sum::<i64>());
        println!("  RR={}", if rr { "yes" } else { "NO" });
        polys.push(coeffs);
    }

    // Recurrence search
    println!("\nRecurrence search:");
    for (label, opts) in &[
        (
            "rec1,d0",
            AdaptiveSearchOptions {
                max_rec_len: 1,
                max_var_deg: 2,
                max_idx_deg: 2,
                max_diff_deg: 0,
                verbose: false,
                ..Default::default()
            },
        ),
        (
            "rec1,d1",
            AdaptiveSearchOptions {
                max_rec_len: 1,
                max_var_deg: 2,
                max_idx_deg: 2,
                max_diff_deg: 1,
                verbose: false,
                ..Default::default()
            },
        ),
        (
            "rec2,d0",
            AdaptiveSearchOptions {
                max_rec_len: 2,
                max_var_deg: 2,
                max_idx_deg: 2,
                max_diff_deg: 0,
                verbose: false,
                ..Default::default()
            },
        ),
        (
            "rec2,d1",
            AdaptiveSearchOptions {
                max_rec_len: 2,
                max_var_deg: 2,
                max_idx_deg: 2,
                max_diff_deg: 1,
                verbose: false,
                ..Default::default()
            },
        ),
    ] {
        match find_recurrence_adaptive(&polys, opts) {
            Some(r) => println!("  {}: {}", label, r.recurrence),
            None => println!("  {}: not found", label),
        }
    }
}

fn main() {
    let max_n_k2: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(14);
    let max_n_k3: u8 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(10);

    run_for_max_entry(2, max_n_k2);
    run_for_max_entry(3, max_n_k3);
    run_for_max_entry(4, 9);
}
