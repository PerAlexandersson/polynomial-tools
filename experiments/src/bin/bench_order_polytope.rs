/// Benchmark: frontier DP vs backtracking for order-preserving map counting.
use combinatoric_core::partition::Partition;
use combinatoric_core::poset::Poset;
use std::time::Instant;

fn bench(label: &str, poset: &Poset, k: usize) {
    let nat = poset.natural_relabeling();

    let t0 = Instant::now();
    let val_dp = poset.count_weak_order_preserving_dp(k);
    let dur_dp = t0.elapsed();

    // Only run backtracking if DP was fast (otherwise it would be hours)
    let (val_bt, dur_bt, speedup) = if dur_dp.as_millis() < 50 {
        let t1 = Instant::now();
        let v = nat.count_weak_order_preserving(k);
        let d = t1.elapsed();
        assert_eq!(v, val_dp, "{} k={}: mismatch", label, k);
        (v, format!("{:>10.2?}", d), format!("{:.1}x", d.as_secs_f64() / dur_dp.as_secs_f64().max(1e-9)))
    } else {
        (val_dp, "  (skipped)".into(), "   >>".into())
    };

    println!(
        "  {:<28} k={:<3} Ω={:<14} bt={:<14} dp={:>10.2?}  speedup={}",
        label, k, val_dp, dur_bt, dur_dp, speedup
    );
}

fn main() {
    println!("=== DP vs backtracking ===\n");

    println!("--- Chains (frontier width 1) ---");
    for n in [5, 10, 13, 15] {
        bench(&format!("chain({})", n), &Poset::chain(n), n);
    }
    // chain(20): DP only (backtracking would take hours)
    let t0 = Instant::now();
    let v = Poset::chain(20).count_weak_order_preserving_dp(20);
    println!("  {:<28} k={:<3} Ω={:<14} bt={:<14} dp={:>10.2?}  speedup=>>", "chain(20)", 20, v, "(hours)", t0.elapsed());

    println!("\n--- Antichains (frontier width 0) ---");
    for n in [5, 8, 10, 14] {
        bench(&format!("antichain({})", n), &Poset::antichain(n), 3);
    }

    println!("\n--- Young diagram posets ---");
    for s in &[vec![3,2,1], vec![4,3,2,1], vec![4,4,4], vec![3,3,3,3]] {
        let p = Poset::from_shape(&Partition::new(s.clone()));
        bench(&format!("shape {:?}", s), &p, p.num_elements() / 2 + 1);
    }

    println!("\n--- Fences (frontier width 1) ---");
    for n in [6, 10] {
        bench(&format!("fence({})", n), &Poset::fence(n), n);
    }
    let t0 = Instant::now();
    let v = Poset::fence(20).count_weak_order_preserving_dp(20);
    println!("  {:<28} k={:<3} Ω={:<14} bt={:<14} dp={:>10.2?}  speedup=>>", "fence(20)", 20, v, "(hours)", t0.elapsed());

    println!("\n--- Ehrhart polynomial (BigRational — no overflow) ---");
    for (label, p) in &[
        ("chain(8)", Poset::chain(8)),
        ("antichain(5)", Poset::antichain(5)),
        ("shape [4,3,2,1]", Poset::from_shape(&Partition::new(vec![4,3,2,1]))),
        ("fence(10)", Poset::fence(10)),
    ] {
        let t0 = Instant::now();
        let hstar = p.order_polytope_hstar();
        let dur = t0.elapsed();
        let nat = p.natural_relabeling();
        let pe = nat.p_eulerian_polynomial();
        assert_eq!(hstar, pe, "{}: h* != P-Eulerian", label);
        println!("  {:<28} n={:<3} time={:>10.2?}  h*={:?}", label, p.num_elements(), dur, hstar);
    }
}
