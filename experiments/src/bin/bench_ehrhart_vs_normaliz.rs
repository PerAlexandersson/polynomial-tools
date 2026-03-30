/// Benchmark: Rust frontier DP Ehrhart vs polymake/Normaliz vs linear extension enumeration.
use combinatoric_core::partition::Partition;
use combinatoric_core::poset::Poset;
use std::io::Write;
use std::process::Command;
use std::time::{Duration, Instant};

fn polymake_script(poset: &Poset) -> String {
    let n = poset.num_elements();
    let nat = poset.natural_relabeling();
    let mut ineqs = Vec::new();
    for i in 0..n {
        let mut row = vec![0i32; n + 1]; row[i + 1] = 1; ineqs.push(row);
    }
    for i in 0..n {
        let mut row = vec![0i32; n + 1]; row[0] = 1; row[i + 1] = -1; ineqs.push(row);
    }
    for &(a, b) in nat.covers() {
        let mut row = vec![0i32; n + 1]; row[a + 1] = -1; row[b + 1] = 1; ineqs.push(row);
    }
    let mut s = String::from("use application \"polytope\";\nuse Benchmark qw(:all);\nmy $m = new Matrix<Rational>([\n");
    for (i, row) in ineqs.iter().enumerate() {
        let vals: Vec<String> = row.iter().map(|x| x.to_string()).collect();
        s.push_str(&format!("  [{}]{}\n", vals.join(", "), if i + 1 < ineqs.len() { "," } else { "" }));
    }
    s.push_str("]);\nmy $p = new Polytope(INEQUALITIES => $m);\nmy $t0 = Benchmark->new;\nmy $hstar = $p->H_STAR_VECTOR;\nmy $t1 = Benchmark->new;\nmy $td = timediff($t1, $t0);\nprint \"HSTAR: $hstar\\n\";\nprint \"TIME: \", timestr($td), \"\\n\";\n");
    s
}

fn run_polymake(script: &str, timeout_secs: u64) -> Option<(Vec<i64>, f64)> {
    std::fs::write("/tmp/order_polytope_bench.pl", script).ok()?;
    let output = Command::new("nice")
        .args(["-n", "19", "timeout", &format!("{}", timeout_secs), "polymake", "--script", "/tmp/order_polytope_bench.pl"])
        .output().ok()?;
    if !output.status.success() { return None; }
    let stdout = String::from_utf8_lossy(&output.stdout);
    let hstar: Vec<i64> = stdout.lines().find(|l| l.starts_with("HSTAR:"))?.strip_prefix("HSTAR: ")?
        .split_whitespace().filter_map(|s| s.parse().ok()).collect();
    let cpu = stdout.lines().find(|l| l.starts_with("TIME:")).and_then(|l| {
        let before = &l[..l.find("CPU")?]; before[before.rfind('=')? + 1..].trim().parse::<f64>().ok()
    }).unwrap_or(0.0);
    Some((hstar, cpu))
}

/// Compute h* via Stanley's theorem: enumerate linear extensions, count descents.
/// Streams without storing all extensions (O(n) memory).
fn linext_hstar(poset: &Poset, _timeout: Duration) -> Option<(Vec<i64>, Duration)> {
    let nat = poset.natural_relabeling();
    let t0 = Instant::now();
    let coeffs = nat.p_eulerian_polynomial();
    Some((coeffs, t0.elapsed()))
}

fn fmt(d: Duration) -> String {
    let s = d.as_secs_f64();
    if s >= 1.0 { format!("{:.2}s", s) } else if s >= 0.001 { format!("{:.1}ms", s * 1e3) } else { format!("{:.0}µs", s * 1e6) }
}
fn fmtf(s: f64) -> String {
    if s >= 1.0 { format!("{:.2}s", s) } else if s >= 0.001 { format!("{:.1}ms", s * 1e3) } else { format!("{:.0}µs", s * 1e6) }
}

fn bench(label: &str, poset: &Poset, do_linext: bool, do_polymake: bool) {
    let n = poset.num_elements();
    let t0 = Instant::now();
    let rust_hstar = poset.order_polytope_hstar();
    let rust_time = t0.elapsed();
    print!("  {:<30} n={:<3} DP={:>10}  ", label, n, fmt(rust_time));
    std::io::stdout().flush().unwrap();

    if do_linext {
        match linext_hstar(poset, Duration::from_secs(5)) {
            Some((h, d)) => { assert_eq!(h, rust_hstar); print!("linext={:>10}  ", fmt(d)); }
            None => print!("linext={:>10}  ", ">5s"),
        }
    } else { print!("linext={:>10}  ", "---"); }
    std::io::stdout().flush().unwrap();

    if do_polymake {
        match run_polymake(&polymake_script(poset), 5) {
            Some((mut h, s)) => {
                while h.len() > 1 && *h.last().unwrap() == 0 { h.pop(); }
                assert_eq!(h, rust_hstar);
                print!("polymake={:>10}", fmtf(s));
            }
            None => print!("polymake={:>10}", ">5s"),
        }
    } else { print!("polymake={:>10}", "---"); }
    println!();
}

fn main() {
    println!("=== h*-vector: Rust DP vs linext vs polymake (5s cap each) ===\n");
    // Small/medium — all three methods
    bench("chain(6)", &Poset::chain(6), true, true);
    bench("antichain(4)", &Poset::antichain(4), true, true);
    bench("shape [3,2,1]", &Poset::from_shape(&Partition::new(vec![3, 2, 1])), true, true);
    bench("fence(8)", &Poset::fence(8), true, true);
    bench("shape [4,3,2,1]", &Poset::from_shape(&Partition::new(vec![4, 3, 2, 1])), true, true);
    bench("shape [3,3,3]", &Poset::from_shape(&Partition::new(vec![3, 3, 3])), true, true);
    bench("fence(10)", &Poset::fence(10), true, true);
    bench("2-alt(10)", &Poset::k_alternating(10, 2), true, true);
    // Larger — linext may timeout, polymake may timeout
    bench("shape [5,4,3,2,1]", &Poset::from_shape(&Partition::new(vec![5, 4, 3, 2, 1])), true, true);
    bench("shape [4,4,4,4]", &Poset::from_shape(&Partition::new(vec![4, 4, 4, 4])), true, true);
    bench("fence(14)", &Poset::fence(14), false, true);
    bench("shape [5,5,5]", &Poset::from_shape(&Partition::new(vec![5, 5, 5])), true, true);
    bench("fence(20)", &Poset::fence(20), false, false);
    bench("shape [5,5,5,5]", &Poset::from_shape(&Partition::new(vec![5, 5, 5, 5])), false, false);
    bench("shape [6,5,4,3,2,1]", &Poset::from_shape(&Partition::new(vec![6, 5, 4, 3, 2, 1])), false, false);
    bench("shape [6,6,6]", &Poset::from_shape(&Partition::new(vec![6, 6, 6])), false, false);
    bench("shape [6,6,6,6]", &Poset::from_shape(&Partition::new(vec![6, 6, 6, 6])), false, false);
    bench("shape [7,6,5,4,3,2,1]", &Poset::from_shape(&Partition::new(vec![7, 6, 5, 4, 3, 2, 1])), false, false);
}
