/// Check real-rootedness of h*-vectors of GT polytopes.
/// Reads lines of the form "lambda=(3,2,1) h*=[1, 8, 35, 32, 9, 0, 0, 0]" from stdin.
use polynomial_tools::is_real_rooted;

fn main() {
    let stdin = std::io::stdin();
    let mut failures = Vec::new();
    let mut count = 0;

    for line in stdin.lines() {
        let line = line.unwrap();
        if line.trim().is_empty() {
            continue;
        }

        // Parse lambda
        let lambda_str = line
            .split("lambda=(")
            .nth(1)
            .and_then(|s| s.split(')').next())
            .unwrap_or("?");

        // Parse h*-vector, stripping trailing zeros
        let hvec_str = line
            .split("h*=[")
            .nth(1)
            .and_then(|s| s.split(']').next())
            .unwrap_or("");

        let mut coeffs: Vec<i64> = hvec_str
            .split(',')
            .filter_map(|s| s.trim().parse().ok())
            .collect();

        // Strip trailing zeros
        while coeffs.last() == Some(&0) {
            coeffs.pop();
        }

        if coeffs.is_empty() {
            continue;
        }

        count += 1;
        let rr = is_real_rooted(&coeffs);
        if !rr {
            println!("NOT real-rooted: lambda=({}) h*={:?}", lambda_str, coeffs);
            failures.push(lambda_str.to_string());
        } else {
            println!(
                "real-rooted:     lambda=({}) deg={}",
                lambda_str,
                coeffs.len() - 1
            );
        }
    }

    println!(
        "\n--- Summary: {}/{} real-rooted ---",
        count - failures.len(),
        count
    );
    if !failures.is_empty() {
        println!("Failures: {:?}", failures);
    }
}
