use combinatoric_core::partition::Partition;
use combinatoric_core::key_polynomial::*;

fn main() {
    let cases: Vec<(&str, Vec<u32>, Vec<usize>)> = vec![
        ("id",       vec![2, 1], vec![1, 2, 3]),
        ("s1",       vec![2, 1], vec![2, 1, 3]),
        ("s2",       vec![2, 1], vec![1, 3, 2]),
        ("s1s2",     vec![2, 1], vec![3, 1, 2]),
        ("s2s1",     vec![2, 1], vec![2, 3, 1]),
        ("w0",       vec![2, 1], vec![3, 2, 1]),
        ("31,s1",    vec![3, 1], vec![2, 1, 3]),
        ("31,w0",    vec![3, 1], vec![3, 2, 1]),
        // n=4
        ("n4,paper", vec![2, 1], vec![2, 4, 3, 1]),
        ("n4,w0",    vec![2, 1], vec![4, 3, 2, 1]),
        ("n4,s1",    vec![2, 1], vec![2, 1, 3, 4]),
    ];

    println!("{:<12} {:>6} {:>6}  {:<30}  {:<30}  {:>5}",
             "case", "deg_p", "deg_r", "plain", "recip", "match");
    println!("{}", "-".repeat(100));

    for (name, lam, sigma) in &cases {
        let lambda = Partition::new(lam.clone());
        let n = sigma.len();
        let d = n * (n - 1) / 2;

        let ep_plain = key_ehrhart_polynomial_positive_only(&lambda, sigma, Some(d));
        let ep_recip = key_ehrhart_polynomial(&lambda, sigma, Some(d));

        // Compare coefficients
        let max_deg = ep_plain.degree.max(ep_recip.degree);
        let mut agree = true;
        for i in 0..=max_deg {
            let cp = if i < ep_plain.coeffs.len() { &ep_plain.coeffs[i] } else { &num_traits::Zero::zero() };
            let cr = if i < ep_recip.coeffs.len() { &ep_recip.coeffs[i] } else { &num_traits::Zero::zero() };
            if cp != cr {
                agree = false;
                break;
            }
        }

        // Also verify both predict the same values at k=0..d+2
        let mut values_agree = true;
        for k in 0..=(d as i64 + 2) {
            if ep_plain.eval(k) != ep_recip.eval(k) {
                values_agree = false;
                break;
            }
        }

        println!("{:<12} {:>6} {:>6}  {:<30}  {:<30}  {:>5}",
                 name,
                 ep_plain.degree, ep_recip.degree,
                 &ep_plain.display()[..ep_plain.display().len().min(30)],
                 &ep_recip.display()[..ep_recip.display().len().min(30)],
                 if agree && values_agree { "OK" } else { "FAIL" });
    }
}
