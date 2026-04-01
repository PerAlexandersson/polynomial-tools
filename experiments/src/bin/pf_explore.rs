use combpoly::parking::all_parking_functions;
use polynomial_tools::real_rootedness as polynomial;

fn des_word(w: &[u8]) -> usize {
    let mut d = 0;
    for i in 0..w.len().saturating_sub(1) {
        if w[i] > w[i + 1] {
            d += 1;
        }
    }
    d
}

fn find_roots_f64(coeffs: &[f64]) -> Option<Vec<f64>> {
    let n = coeffs.len() - 1;
    if n == 0 {
        return Some(vec![]);
    }
    if n == 1 {
        return Some(vec![-coeffs[0] / coeffs[1]]);
    }
    let mut roots = Vec::new();
    let mut p = coeffs.to_vec();
    for _ in 0..n {
        if p.len() <= 1 {
            break;
        }
        let r = newton_root(&p)?;
        roots.push(r);
        let deg = p.len() - 1;
        let mut q = vec![0.0; deg];
        q[deg - 1] = p[deg];
        for j in (0..deg - 1).rev() {
            q[j] = p[j + 1] + r * q[j + 1];
        }
        p = q;
    }
    Some(roots)
}

fn newton_root(coeffs: &[f64]) -> Option<f64> {
    for &start in &[
        -1000.0, -100.0, -10.0, -5.0, -2.0, -1.0, -0.5, -0.1, 0.0, 0.5,
    ] {
        let mut x = start;
        for _ in 0..500 {
            let (val, deriv) = eval_pd(coeffs, x);
            if val.abs() < 1e-12 {
                return Some(x);
            }
            if deriv.abs() < 1e-15 {
                break;
            }
            x -= val / deriv;
        }
        let (val, _) = eval_pd(coeffs, x);
        if val.abs() < 1e-8 {
            return Some(x);
        }
    }
    None
}

fn eval_pd(coeffs: &[f64], x: f64) -> (f64, f64) {
    let mut val = 0.0;
    let mut deriv = 0.0;
    for i in (0..coeffs.len()).rev() {
        deriv = deriv * x + val;
        val = val * x + coeffs[i];
    }
    (val, deriv)
}

fn check_interlaces(fc: &[i64], gc: &[i64]) -> bool {
    let ff: Vec<f64> = fc.iter().map(|&c| c as f64).collect();
    let gf: Vec<f64> = gc.iter().map(|&c| c as f64).collect();
    let (Some(mut fr), Some(mut gr)) = (find_roots_f64(&ff), find_roots_f64(&gf)) else {
        return false;
    };
    fr.sort_by(|a, b| a.partial_cmp(b).unwrap());
    gr.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let (df, dg) = (fr.len(), gr.len());
    if df == dg {
        fr.iter().zip(gr.iter()).all(|(fi, gi)| fi <= &(gi + 1e-8))
    } else if df + 1 == dg {
        (0..df).all(|i| gr[i] <= fr[i] + 1e-8 && fr[i] <= gr[i + 1] + 1e-8)
    } else {
        false
    }
}

fn main() {
    let max_n: u8 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);
    println!("=== Descent-generating polynomials for parking functions ===\n");
    let mut polys: Vec<Vec<i64>> = Vec::new();
    for n in 1u8..=max_n {
        let pfs = all_parking_functions(n);
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
        print!("n={}: ", n);
        for (i, &c) in coeffs.iter().enumerate() {
            if i > 0 {
                print!("\t");
            }
            print!("{}", c);
        }
        println!("  (RR={})", if rr { "yes" } else { "NO" });
        polys.push(coeffs);
    }
    println!("\n=== Interlacing: P_n ≼ P_{{n+1}} ===");
    for i in 0..polys.len() - 1 {
        println!(
            "P_{} ≼ P_{}: {}",
            i + 1,
            i + 2,
            check_interlaces(&polys[i], &polys[i + 1])
        );
    }
    println!("\n=== Reversed: P_{{n+1}} ≼ P_n ===");
    for i in 0..polys.len() - 1 {
        println!(
            "P_{} ≼ P_{}: {}",
            i + 2,
            i + 1,
            check_interlaces(&polys[i + 1], &polys[i])
        );
    }
    println!("\n=== Recurrence search ===");
    use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};
    let search = AdaptiveSearchOptions {
        max_rec_len: 3,
        max_var_deg: 3,
        max_idx_deg: 3,
        max_diff_deg: 2,
        verbose: false,
        ..Default::default()
    };
    match find_recurrence_adaptive(&polys, &search) {
        Some(r) => println!("RECURRENCE: {}", r.recurrence),
        None => println!("No simple recurrence found."),
    }
}
