use clap::Parser;
use combinatoric_core::key_polynomial::key_ehrhart_polynomial;
use combinatoric_core::Partition;
use dotenv::from_path;
use mysql::prelude::*;
use mysql::*;
use polynomial_tools::ehrhart_to_hstar;
use std::collections::HashSet;
use std::path::Path;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(
    about = "Batch-compute key Ehrhart data and store it in MariaDB",
    long_about = "Compute key Ehrhart data for pairs (lambda, sigma), where sigma is a permutation of [1..n] and n is the chosen ambient rank. Results are stored in a MariaDB table with crash-safe resume behavior."
)]
struct Args {
    /// Compute only this partition lambda, e.g. 3,2,1,0.
    #[arg(long)]
    lambda: Option<String>,

    /// Maximum partition size |lambda| to compute.
    #[arg(long)]
    max_n: Option<u32>,

    /// Minimum partition size |lambda| to compute.
    #[arg(long)]
    min_n: Option<u32>,

    /// Ambient rank n for sigma. If omitted, use lambda.num_parts().
    #[arg(long)]
    rank: Option<usize>,

    /// MariaDB connection URL. If omitted, build it from env vars in kostka/.env.
    #[arg(long, env = "KOSTKA_DB_URL")]
    db_url: Option<String>,

    /// Only print the work items without computing or writing.
    #[arg(long)]
    dry_run: bool,

    /// Limit to the first N permutations for each lambda, in lexicographic order.
    #[arg(long)]
    limit_perms: Option<usize>,
}

const CREATE_TABLE: &str = r"
CREATE TABLE IF NOT EXISTS key_ehrhart (
    id              BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    lambda          VARCHAR(255) NOT NULL,
    sigma           VARCHAR(255) NOT NULL,
    dimension       INT NOT NULL,
    polynomial      TEXT NOT NULL,
    hstar           TEXT NOT NULL,
    elapsed_ms      BIGINT UNSIGNED,
    created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE KEY uq_key_shape (lambda, sigma)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
";

fn parts_str(p: &Partition) -> String {
    p.parts()
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn sigma_str(sigma: &[usize]) -> String {
    sigma
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn hstar_json(h: &[i64]) -> String {
    let elems: Vec<String> = h.iter().map(|x| x.to_string()).collect();
    format!("[{}]", elems.join(","))
}

fn parse_partition(s: &str) -> Result<Partition, String> {
    let parts: Result<Vec<u32>, _> = s.split(',').map(|x| x.trim().parse::<u32>()).collect();
    parts
        .map(Partition::new)
        .map_err(|e| format!("Invalid partition: {}", e))
}

fn load_db_url(cli_url: Option<String>) -> String {
    if let Some(url) = cli_url {
        return url;
    }

    let env_path = Path::new("/home/paxinum/Dropbox/projects/rust/kostka/.env");
    let _ = from_path(env_path);

    let host = match std::env::var("DB_HOST") {
        Ok(host) if host == "localhost" => "127.0.0.1".to_string(),
        Ok(host) => host,
        Err(_) => "127.0.0.1".to_string(),
    };
    let user = std::env::var("DB_USER").expect("DB_USER missing and no --db-url supplied");
    let password =
        std::env::var("DB_PASSWORD").expect("DB_PASSWORD missing and no --db-url supplied");
    let db = std::env::var("DB_NAME").expect("DB_NAME missing and no --db-url supplied");
    let _charset = std::env::var("DB_CHARSET").unwrap_or_else(|_| "utf8mb4".to_string());

    format!("mysql://{}:{}@{}/{}", user, password, host, db)
}

fn all_permutations(n: usize) -> Vec<Vec<usize>> {
    let mut a: Vec<usize> = (1..=n).collect();
    let mut out = vec![a.clone()];
    while next_permutation(&mut a) {
        out.push(a.clone());
    }
    out
}

fn next_permutation(a: &mut [usize]) -> bool {
    if a.len() < 2 {
        return false;
    }

    let mut i = a.len() - 2;
    loop {
        if a[i] < a[i + 1] {
            break;
        }
        if i == 0 {
            a.reverse();
            return false;
        }
        i -= 1;
    }

    let mut j = a.len() - 1;
    while a[j] <= a[i] {
        j -= 1;
    }
    a.swap(i, j);
    a[i + 1..].reverse();
    true
}

fn already_computed(conn: &mut PooledConn, lambda_str: &str) -> HashSet<String> {
    let rows: Vec<String> = conn
        .exec(
            "SELECT sigma FROM key_ehrhart WHERE lambda = ?",
            (lambda_str,),
        )
        .unwrap_or_default();
    rows.into_iter().collect()
}

fn main() {
    let args = Args::parse();
    let db_url = load_db_url(args.db_url);

    let pool = if args.dry_run {
        None
    } else {
        let opts = Opts::from_url(&db_url).expect("Invalid database URL");
        let pool = Pool::new(opts).expect("Failed to connect to database");
        let mut conn = pool.get_conn().expect("Failed to get connection");
        conn.query_drop(CREATE_TABLE)
            .expect("Failed to create key_ehrhart table");
        Some(pool)
    };

    let mut total_computed = 0u64;
    let mut total_skipped = 0u64;

    let lambdas: Vec<Partition> = if let Some(lambda_str) = args.lambda.as_deref() {
        vec![parse_partition(lambda_str).expect("invalid --lambda")]
    } else {
        let min_n = args.min_n.unwrap_or(1);
        let max_n = args.max_n.expect("either --lambda or --max-n is required");
        let mut all = Vec::new();
        for size in min_n..=max_n {
            all.extend(Partition::all_of_size(size));
        }
        all
    };

    for lambda in lambdas {
        let lambda_str = parts_str(&lambda);
        let rank = args.rank.unwrap_or(lambda.num_parts());
        if rank < lambda.num_parts() {
            eprintln!(
                "Skipping lambda=({}) because rank {} < num_parts {}",
                lambda_str,
                rank,
                lambda.num_parts()
            );
            continue;
        }

        let mut sigmas = all_permutations(rank);
        if let Some(limit) = args.limit_perms {
            sigmas.truncate(limit);
        }

        let done = if let Some(ref pool) = pool {
            let mut conn = pool.get_conn().expect("DB connection failed");
            already_computed(&mut conn, &lambda_str)
        } else {
            HashSet::new()
        };

        for sigma in sigmas {
            let sigma_text = sigma_str(&sigma);
            if done.contains(&sigma_text) {
                total_skipped += 1;
                continue;
            }

            if args.dry_run {
                println!("lambda=({}) sigma=({})", lambda_str, sigma_text);
                total_computed += 1;
                continue;
            }

            let start = Instant::now();
            let poly = key_ehrhart_polynomial(&lambda, &sigma, None);
            let elapsed_ms = start.elapsed().as_millis() as u64;
            let hstar = ehrhart_to_hstar(&poly.coeffs);

            if let Some(ref pool) = pool {
                let mut conn = pool.get_conn().expect("DB connection failed");
                conn.exec_drop(
                    r"INSERT IGNORE INTO key_ehrhart
                          (lambda, sigma, dimension, polynomial, hstar, elapsed_ms)
                          VALUES (?, ?, ?, ?, ?, ?)",
                    (
                        &lambda_str,
                        &sigma_text,
                        poly.degree as i32,
                        poly.display(),
                        hstar_json(&hstar),
                        elapsed_ms,
                    ),
                )
                .expect("DB insert failed");
            }

            total_computed += 1;
        }

        eprintln!(
            "lambda=({}) rank={} -> {} sigmas processed ({} skipped so far)",
            lambda_str, rank, total_computed, total_skipped
        );
    }

    eprintln!(
        "Done: {} computed, {} skipped.",
        total_computed, total_skipped
    );
}
