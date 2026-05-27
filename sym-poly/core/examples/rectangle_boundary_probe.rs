use std::collections::BTreeMap;

use sym_poly_core::Ssaf;

#[derive(Debug, Clone)]
struct Config {
    height: usize,
    width: usize,
    basement: Vec<u32>,
    first: Vec<u32>,
    last: Vec<u32>,
    show: bool,
}

fn main() {
    let config = parse_args();
    validate_config(&config);

    let ssaf_counts = ssaf_rectangle_counts(&config);
    let direct_counts = direct_rectangle_counts(&config);

    println!("height={}, width={}", config.height, config.width);
    println!("basement={:?}", config.basement);
    println!("first={:?}", config.first);
    println!("last={:?}", config.last);
    println!("ssaf weights: {}", total_count(&ssaf_counts));
    println!(
        "direct flagged-rectangle weights: {}",
        total_count(&direct_counts)
    );
    println!("weight enumerators match: {}", ssaf_counts == direct_counts);

    if ssaf_counts != direct_counts {
        print_difference("in SSAF only", &ssaf_counts, &direct_counts);
        print_difference("in direct only", &direct_counts, &ssaf_counts);
    }

    let bad_columns = filtered_ssafs(&config)
        .into_iter()
        .filter(|f| !filling_columns_strictly_decrease(f.rows(), config.width))
        .collect::<Vec<_>>();
    println!(
        "filtered SSAFs with a non-decreasing filling column: {}",
        bad_columns.len()
    );

    if config.show {
        println!("\nFiltered SSAFs:");
        for f in filtered_ssafs(&config) {
            println!("{f}\n");
        }

        println!("Direct flagged rectangles:");
        for grid in direct_rectangles(&config) {
            print_grid(&grid);
            println!();
        }
    }
}

fn parse_args() -> Config {
    let mut config = Config {
        height: 3,
        width: 3,
        basement: vec![5, 4, 3],
        first: vec![5, 4, 3],
        last: vec![3, 2, 1],
        show: false,
    };

    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--height" => {
                i += 1;
                config.height = args[i].parse().expect("--height expects an integer");
            }
            "--width" => {
                i += 1;
                config.width = args[i].parse().expect("--width expects an integer");
            }
            "--basement" => {
                i += 1;
                config.basement = parse_vec(&args[i]);
            }
            "--first" => {
                i += 1;
                config.first = parse_vec(&args[i]);
            }
            "--last" => {
                i += 1;
                config.last = parse_vec(&args[i]);
            }
            "--show" => {
                config.show = true;
            }
            "--help" | "-h" => {
                print_help_and_exit();
            }
            flag => panic!("unknown argument: {flag}"),
        }
        i += 1;
    }

    config
}

fn parse_vec(s: &str) -> Vec<u32> {
    if s.trim().is_empty() {
        return Vec::new();
    }
    s.split(',')
        .map(|part| part.parse().expect("vector entries must be integers"))
        .collect()
}

fn print_help_and_exit() -> ! {
    println!(
        "Usage: cargo run -p sym-poly-core --example rectangle_boundary_probe -- \\
         [--height H] [--width W] [--basement a,b,c] [--first a,b,c] \\
         [--last a,b,c] [--show]"
    );
    std::process::exit(0);
}

fn validate_config(config: &Config) {
    assert!(config.height > 0, "height must be positive");
    assert!(config.width > 0, "width must be positive");
    assert_eq!(config.basement.len(), config.height);
    assert_eq!(config.first.len(), config.height);
    assert_eq!(config.last.len(), config.height);

    assert!(
        column_strictly_decreases(&config.first),
        "first column must strictly decrease"
    );
    assert!(
        column_strictly_decreases(&config.last),
        "last column must strictly decrease"
    );

    for r in 0..config.height {
        assert!(
            config.basement[r] >= config.first[r],
            "basement must weakly dominate the first column rowwise"
        );
        assert!(
            config.first[r] >= config.last[r],
            "rows must be able to weakly decrease from first to last"
        );
    }
}

fn ssaf_rectangle_counts(config: &Config) -> BTreeMap<Vec<u32>, usize> {
    let alphabet = alphabet_size(config);
    let mut counts = BTreeMap::new();
    for filling in filtered_ssafs(config) {
        let weight = weight_from_rows(filling.rows(), alphabet);
        *counts.entry(weight).or_insert(0) += 1;
    }
    counts
}

fn filtered_ssafs(config: &Config) -> Vec<Ssaf> {
    let shape = vec![config.width as u32; config.height];
    Ssaf::generate(&shape, &config.basement)
        .into_iter()
        .filter(|f| {
            for r in 0..config.height {
                if f.rows()[r][1] != config.first[r] {
                    return false;
                }
                if f.rows()[r][config.width] != config.last[r] {
                    return false;
                }
            }
            true
        })
        .collect()
}

fn direct_rectangle_counts(config: &Config) -> BTreeMap<Vec<u32>, usize> {
    let alphabet = alphabet_size(config);
    let mut counts = BTreeMap::new();
    for grid in direct_rectangles(config) {
        let weight = weight_from_grid(&grid, alphabet);
        *counts.entry(weight).or_insert(0) += 1;
    }
    counts
}

fn direct_rectangles(config: &Config) -> Vec<Vec<Vec<u32>>> {
    let mut grid = vec![vec![0; config.width]; config.height];
    for r in 0..config.height {
        grid[r][0] = config.first[r];
        grid[r][config.width - 1] = config.last[r];
    }

    let mut result = Vec::new();
    fill_direct_cell(config, &mut grid, 0, 0, &mut result);
    result
}

fn fill_direct_cell(
    config: &Config,
    grid: &mut [Vec<u32>],
    r: usize,
    c: usize,
    result: &mut Vec<Vec<Vec<u32>>>,
) {
    if r == config.height {
        result.push(grid.to_vec());
        return;
    }

    let (next_r, next_c) = if c + 1 == config.width {
        (r + 1, 0)
    } else {
        (r, c + 1)
    };

    if c == 0 || c + 1 == config.width {
        if local_direct_ok(grid, r, c) {
            fill_direct_cell(config, grid, next_r, next_c, result);
        }
        return;
    }

    for value in 1..=alphabet_size(config) as u32 {
        grid[r][c] = value;
        if local_direct_ok(grid, r, c) {
            fill_direct_cell(config, grid, next_r, next_c, result);
        }
    }
    grid[r][c] = 0;
}

fn local_direct_ok(grid: &[Vec<u32>], r: usize, c: usize) -> bool {
    let value = grid[r][c];
    if value == 0 {
        return false;
    }
    if c > 0 && grid[r][c - 1] != 0 && grid[r][c - 1] < value {
        return false;
    }
    if r > 0 && grid[r - 1][c] != 0 && grid[r - 1][c] <= value {
        return false;
    }
    true
}

fn filling_columns_strictly_decrease(rows: &[Vec<u32>], width: usize) -> bool {
    for c in 1..=width {
        let col = rows.iter().map(|row| row[c]).collect::<Vec<_>>();
        if !column_strictly_decreases(&col) {
            return false;
        }
    }
    true
}

fn column_strictly_decreases(col: &[u32]) -> bool {
    col.windows(2).all(|w| w[0] > w[1])
}

fn alphabet_size(config: &Config) -> usize {
    *config
        .basement
        .iter()
        .chain(config.first.iter())
        .chain(config.last.iter())
        .max()
        .expect("nonempty data") as usize
}

fn weight_from_rows(rows: &[Vec<u32>], alphabet: usize) -> Vec<u32> {
    let mut weight = vec![0; alphabet];
    for row in rows {
        for &value in row.iter().skip(1) {
            weight[value as usize - 1] += 1;
        }
    }
    weight
}

fn weight_from_grid(grid: &[Vec<u32>], alphabet: usize) -> Vec<u32> {
    let mut weight = vec![0; alphabet];
    for row in grid {
        for &value in row {
            weight[value as usize - 1] += 1;
        }
    }
    weight
}

fn total_count(counts: &BTreeMap<Vec<u32>, usize>) -> usize {
    counts.values().sum()
}

fn print_difference(label: &str, lhs: &BTreeMap<Vec<u32>, usize>, rhs: &BTreeMap<Vec<u32>, usize>) {
    let diff = lhs
        .iter()
        .filter_map(|(weight, &count)| {
            let other = rhs.get(weight).copied().unwrap_or(0);
            (count > other).then_some((weight, count - other))
        })
        .collect::<Vec<_>>();

    if diff.is_empty() {
        return;
    }

    println!("{label}:");
    for (weight, count) in diff.into_iter().take(10) {
        println!("  {count} * x^{weight:?}");
    }
}

fn print_grid(grid: &[Vec<u32>]) {
    for row in grid {
        println!("{row:?}");
    }
}
