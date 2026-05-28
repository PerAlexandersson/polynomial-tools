use std::env;
use std::time::Instant;

use sym_poly_core::packed_sparse_linear_algebra::{packed_sparse_rref, PackedSparseRow};
use sym_poly_core::sparse_linear_algebra::{sparse_rref, sparse_vector, SparseVector};
use sym_poly_core::{PrimeField, Ring};

type F251 = PrimeField<251>;

fn main() {
    let args = env::args().skip(1).collect::<Vec<_>>();
    let config = Config::parse(&args).unwrap_or_else(|message| {
        eprintln!("{message}");
        eprintln!("usage:");
        eprintln!("  cargo run -p sym-poly-core --example packed_sparse_benchmark");
        eprintln!("  cargo run -p sym-poly-core --example packed_sparse_benchmark -- 2000 1000 8");
        std::process::exit(2);
    });

    println!(
        "rows={}, cols={}, nnz_per_row={}, prime=251",
        config.rows, config.cols, config.nnz_per_row
    );

    let raw_rows = generate_rows(config.rows, config.cols, config.nnz_per_row);
    let generic_rows = raw_rows
        .iter()
        .map(|row| {
            sparse_vector(
                config.cols,
                row.iter()
                    .copied()
                    .map(|(col, value)| (col, F251::from_i64(value as i64))),
            )
        })
        .collect::<Vec<SparseVector<F251>>>();
    let packed_rows = raw_rows
        .iter()
        .map(|row| PackedSparseRow::new::<251, _>(config.cols, row.iter().copied()))
        .collect::<Vec<_>>();

    let start = Instant::now();
    let generic = sparse_rref(config.cols, &generic_rows);
    let generic_elapsed = start.elapsed();

    let start = Instant::now();
    let packed = packed_sparse_rref::<251>(config.cols, &packed_rows);
    let packed_elapsed = start.elapsed();

    println!(
        "generic rank={} elapsed={:.3}s",
        generic.rank,
        generic_elapsed.as_secs_f64()
    );
    println!(
        "packed  rank={} elapsed={:.3}s",
        packed.rank,
        packed_elapsed.as_secs_f64()
    );
    println!(
        "same pivot profile: {}",
        generic.pivot_columns == packed.pivot_columns
    );
}

#[derive(Debug, Clone, Copy)]
struct Config {
    rows: usize,
    cols: usize,
    nnz_per_row: usize,
}

impl Config {
    fn parse(args: &[String]) -> Result<Self, String> {
        if args.is_empty() {
            return Ok(Self {
                rows: 1200,
                cols: 700,
                nnz_per_row: 8,
            });
        }
        if args.len() != 3 {
            return Err("expected either no arguments or rows cols nnz_per_row".to_string());
        }
        let rows = parse_positive_usize(&args[0], "rows")?;
        let cols = parse_positive_usize(&args[1], "cols")?;
        let nnz_per_row = parse_positive_usize(&args[2], "nnz_per_row")?;
        Ok(Self {
            rows,
            cols,
            nnz_per_row,
        })
    }
}

fn parse_positive_usize(text: &str, name: &str) -> Result<usize, String> {
    let value = text
        .parse::<usize>()
        .map_err(|_| format!("{name} must be a positive integer"))?;
    if value == 0 {
        return Err(format!("{name} must be positive"));
    }
    Ok(value)
}

fn generate_rows(rows: usize, cols: usize, nnz_per_row: usize) -> Vec<Vec<(usize, u8)>> {
    let mut rng = Lcg::new(0x9e3779b97f4a7c15);
    let mut result = Vec::with_capacity(rows);
    for row_index in 0..rows {
        let mut row = Vec::with_capacity(nnz_per_row + 1);
        let diagonal_col = row_index % cols;
        row.push((diagonal_col, 1 + (rng.next_usize(250) as u8)));
        for _ in 0..nnz_per_row {
            let col = rng.next_usize(cols);
            let value = 1 + (rng.next_usize(250) as u8);
            row.push((col, value));
        }
        result.push(row);
    }
    result
}

#[derive(Debug, Clone, Copy)]
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    fn next_usize(&mut self, modulus: usize) -> usize {
        (self.next_u64() as usize) % modulus
    }
}
