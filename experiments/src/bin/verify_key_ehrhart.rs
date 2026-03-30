/// Verify key-Kostka Ehrhart polynomials against the Osaka paper.
///
/// Paper Table 2 gives Ehrhart polys for λ=(a+b, a, 0), all σ ∈ S₃.
/// We test several (a,b) pairs and compare.
use combinatoric_core::partition::Partition;
use combinatoric_core::key_polynomial::*;

fn main() {
    println!("=== Verify Ehrhart polynomials against Osaka paper ===\n");

    // Paper formulas for λ=(a+b, a, 0), σ ∈ S₃:
    // [1,2,3]: P(k) = 1
    // [2,1,3]: P(k) = 1 + bk
    // [1,3,2]: P(k) = 1 + ak
    // [3,1,2]: P(k) = (bk+1)(2ak+bk+2)/2
    // [2,3,1]: P(k) = (ak+1)(2bk+ak+2)/2
    // [3,2,1]: P(k) = (ak+1)(bk+1)(2+ak+bk)/2

    let test_abs: Vec<(u32, u32)> = vec![
        (1, 1),  // λ=(2,1,0)
        (2, 1),  // λ=(3,2,0)
        (1, 2),  // λ=(3,1,0)
        (2, 2),  // λ=(4,2,0)
        (3, 1),  // λ=(4,3,0)
        (1, 3),  // λ=(4,1,0)
    ];

    let perms: Vec<(&str, Vec<usize>)> = vec![
        ("[1,2,3]", vec![1,2,3]),
        ("[2,1,3]", vec![2,1,3]),
        ("[1,3,2]", vec![1,3,2]),
        ("[3,1,2]", vec![3,1,2]),
        ("[2,3,1]", vec![2,3,1]),
        ("[3,2,1]", vec![3,2,1]),
    ];

    let mut all_ok = true;

    for &(a, b) in &test_abs {
        let lambda = Partition::new(vec![a + b, a]);
        println!("--- λ=({},{},0), a={}, b={} ---", a+b, a, a, b);

        for (name, sigma) in &perms {
            let ep = key_ehrhart_polynomial(&lambda, sigma, Some(3));

            // Compute paper's formula at k=0..5
            let paper_values: Vec<i64> = (0..=5)
                .map(|k| paper_formula(a, b, sigma, k))
                .collect();

            let rust_values: Vec<i64> = (0..=5)
                .map(|k| {
                    let v = ep.eval(k);
                    // Should be integer
                    assert!(v.is_integer(), "Non-integer at k={}: {}", k, v);
                    v.to_integer().try_into().unwrap()
                })
                .collect();

            let matches = paper_values == rust_values;
            if !matches { all_ok = false; }

            println!(
                "  {:<8} poly: {:<35}  {}",
                name,
                ep.display(),
                if matches { "OK" } else { "MISMATCH" }
            );
            if !matches {
                println!("    paper: {:?}", paper_values);
                println!("    rust:  {:?}", rust_values);
            }
        }
    }

    println!("\n{}", if all_ok { "ALL MATCH" } else { "SOME MISMATCHES" });
}

fn paper_formula(a: u32, b: u32, sigma: &[usize], k: i64) -> i64 {
    let ak = a as i64 * k;
    let bk = b as i64 * k;
    match sigma {
        [1,2,3] => 1,
        [2,1,3] => 1 + bk,
        [1,3,2] => 1 + ak,
        [3,1,2] => (bk + 1) * (2*ak + bk + 2) / 2,
        [2,3,1] => (ak + 1) * (2*bk + ak + 2) / 2,
        [3,2,1] => (ak + 1) * (bk + 1) * (2 + ak + bk) / 2,
        _ => panic!("unknown sigma"),
    }
}
