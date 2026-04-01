//! Explore lattice path interpretations for M_{n,j} = sum_r C(hr,j)*C(j+1, n-2j-2r)
//! for the (h,h,1,1) family, looking for an LGV-type proof of total positivity.
//!
//! Key idea: can M_{n,j} be expressed as the number of lattice paths from A_n to B_j
//! in some planar graph where non-crossing only occurs for the identity matching?

fn binom(n: i64, k: i64) -> i64 {
    if k < 0 || k > n || n < 0 {
        return 0;
    }
    let k = k.min(n - k) as usize;
    let n = n as usize;
    let mut result = 1i64;
    for i in 0..k {
        result = result * (n - i) as i64 / (i + 1) as i64;
    }
    result
}

fn m_entry(h: i64, n: i64, j: i64) -> i64 {
    let mut total = 0i64;
    for r in 0..=100 {
        let s = n - 2 * j - 2 * r;
        if s < 0 {
            break;
        }
        if s > j + 1 {
            continue;
        }
        total += binom(h * r, j) * binom(j + 1, s);
    }
    total
}

fn main() {
    // Print M_{n,j} for h=2 and look for patterns
    println!("=== h=2: M_{{n,j}} ===");
    println!(
        "{:>4} | {:>6} {:>6} {:>6} {:>6} {:>6} {:>6}",
        "n", "j=0", "j=1", "j=2", "j=3", "j=4", "j=5"
    );
    for n in 0..=20 {
        let row: Vec<String> = (0..=5)
            .map(|j| format!("{:>6}", m_entry(2, n, j)))
            .collect();
        println!("{:>4} | {}", n, row.join(" "));
    }

    // For the rigid case, M_{n,j} = C(d(n-j*ell), j) = C(d*n - j*(d*ell+1) + j, j)
    // which counts paths from A_n = (d*n, 0) to B_j = ((d*ell+1)*j, j)
    // using West(-1,0) and North(0,1) steps.
    //
    // For (h,h,1,1): M_{n,j} = sum_r C(hr,j) C(j+1, n-2j-2r)
    //
    // Let's see if we can interpret this as paths on a lattice with two types of steps.
    //
    // Model: paths from A_n to B_j using:
    //   Step type 1: "double-West" = (-2, 0) -- one unit of x^2 spacing
    //   Step type 2: "North" = (0, +1) -- one tile
    //   Step type 3: "optional-West" = (-1, 0) -- the (1+x) factor
    //
    // A path from (n, 0) to (0, j) would use:
    //   - j North steps
    //   - some number of double-West and optional-West steps to reach x=0
    //   Total: j + (number of double-Wests) + (number of optional-Wests) = ??
    //
    // The constraint is: n = 2*(double-Wests) + (optional-Wests)
    // And: optional-Wests <= j+1 (from C(j+1,...))
    //
    // With r = number of double-Wests after accounting for the 2j from North steps:
    //   Total x consumed = 2j (from North steps taking x^2 each) + 2r (double-Wests) + s (optional)
    //   So n = 2j + 2r + s, s = n - 2j - 2r
    //   C(j+1, s) ways to place the optional-Wests among j+1 slots
    //   C(hr, j) ways to place the North steps...
    //
    // Wait, C(hr, j) doesn't directly count "placements" in a path.
    // C(hr, j) = C(hr, j) counts choosing j items from hr.
    // If we think of r "blocks" each of height h, and we choose j of the hr total "slots",
    // that's like choosing which rows to place tiles in.
    //
    // So the model could be:
    // 1. The tiling has r "super-blocks", each providing h potential tile positions
    // 2. We choose j of the hr positions for tiles: C(hr, j)
    // 3. The j tiles create j+1 "gaps" (before first, between consecutive, after last)
    // 4. Each gap optionally adds 1 to the width: C(j+1, s) for s optional extensions
    // 5. Total width = base + 2r + s = n

    // But this doesn't directly give a planar lattice path model for LGV.
    //
    // Alternative: look at the "weighted" version.
    // Define sources A_n and sinks B_j on a lattice:

    println!("\n=== Trying to match with path model ===");
    // For h=2, rigid ell=3: A_n = (2n, 0), B_j = (7j, j) gives C(2n-7j+j, j) = C(2(n-3j), j)
    // But our M_{n,j} is NOT C(2(n-3j), j).
    // For n=8, j=2: C(2(8-6),2) = C(4,2) = 6, but M_{8,2} = 9. Different!

    // So the rigid lattice path model doesn't apply directly.
    // We need a different geometry.

    // Let me check: is there a single-binomial formula hidden here?
    // For h=2: M_{n,j} = sum_r C(2r,j)*C(j+1, n-2j-2r)
    //
    // Let's compute for h=2, j=2:
    println!("\nh=2, j=2 detailed:");
    for n in 6..=16 {
        let mut terms = Vec::new();
        for r in 0..=20 {
            let s = n - 2 * 2 - 2 * r;
            if s < 0 {
                break;
            }
            if s > 3 {
                continue;
            }
            let val = binom(2 * r, 2) * binom(3, s);
            if val > 0 {
                terms.push(format!(
                    "C({},2)*C(3,{})={}*{}={}",
                    2 * r,
                    s,
                    binom(2 * r, 2),
                    binom(3, s),
                    val
                ));
            }
        }
        println!(
            "  M_{{{},2}} = {} = {}",
            n,
            m_entry(2, n, 2),
            terms.join(" + ")
        );
    }

    // Check if the M matrix can be written as a product of two TP matrices
    // M = A * B where A_{n,r} and B_{r,j} are TP
    // A_{n,r} would need to depend only on n,r (not j) -- but C(j+1, n-2j-2r) depends on j!
    // So direct factorization fails.

    // However, there IS a factorization if we work with a "lifted" matrix.
    // Consider the generating function approach:
    // [t^j] F = (1+x)^{j+1} * x^{2j} * G_j(x^2)
    // where G_j(z) = sum_r C(hr,j) z^r
    //
    // The factor (1+x)^{j+1} is the key difference from the rigid case.
    // In the rigid case, this factor was (1+x)^0 = 1 (trivially TP).
    // Here it grows with j.

    // Total positivity of the PRODUCT of a TP matrix with a specific structured matrix
    // might follow from the Cauchy-Binet identity if both factors are TP.

    // Let's check: define
    //   N_{n,j} = C(j+1, n-2j) -- the "Pascal-like" part
    //   P_{n,j} = sum_r C(hr,j) * [x^{n-2r}] ... hmm, this doesn't factor cleanly.

    println!("\n=== Checking if (1+x) substitution preserves TP ===");
    println!("For rigid (h,h,1), d=h, ell=2:");
    println!("  M_{{n,j}}^rigid = C(h(n-2j), j) -- this is TP (proved in paper)");
    println!("\nFor (h,h,1,1), d=h:");
    println!("  M_{{n,j}} = sum_r C(hr,j) C(j+1, n-2j-2r)");
    println!("          = [x^{{n-2j}}] (1+x)^{{j+1}} sum_r C(hr,j) x^{{2r}}");
    println!("\nThe operation f(x) -> [x^n](1+x)^{{j+1}} f(x) applied to the rigid GF");
    println!("is a 'binomial transform' in n for each fixed j.");
    println!("If the rigid matrix is TP, does this transform preserve TP?");

    // The binomial transform of a TP matrix is NOT always TP.
    // But the specific structure here might help.

    // Let's verify: rigid (h,h,1) with d=h has M^rig_{n,j} = C(h*n - j*(2h+1) + j, j) = C(h(n-2j), j)
    // Wait: ell=2, d=h. Rigid formula: C(d(n-j*ell), j) = C(h(n-2j), j).
    // For h=2: C(2(n-2j), j) = C(2n-4j, j).

    println!("\n=== Comparison: rigid vs (h,h,1,1) for h=2 ===");
    println!(
        "{:>4} | {:>8} {:>8} {:>8} | {:>8} {:>8} {:>8}",
        "n", "M^rig_0", "M^rig_1", "M^rig_2", "M_0", "M_1", "M_2"
    );
    for n in 0..=14 {
        let rig: Vec<String> = (0..=2)
            .map(|j| format!("{:>8}", binom(2 * (n - 2 * j), j)))
            .collect();
        let our: Vec<String> = (0..=2)
            .map(|j| format!("{:>8}", m_entry(2, n, j)))
            .collect();
        println!("{:>4} | {} | {}", n, rig.join(" "), our.join(" "));
    }
}
