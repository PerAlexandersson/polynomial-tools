use clap::{Args, Parser, Subcommand, ValueEnum};
use std::io::{self, Write};

use combpoly::{order, parking, permutation, polynomial_builder, statistics, word};
use polynomial_tools::real_rootedness::{format_poly, is_log_concave, is_real_rooted};
use polynomial_tools::recurrence;
use statistics::Stat;

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OrderType {
    Bruhat,
    Weak,
}

#[derive(Parser)]
#[command(name = "combpoly")]
#[command(about = "Explore combinatorial polynomials from permutations and words")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Compute a generating polynomial over a set of objects
    Poly {
        #[command(flatten)]
        source: Source,

        /// Statistic to use
        #[arg(long, value_enum, default_value_t = Stat::Des)]
        stat: Stat,

        /// Check real-rootedness
        #[arg(long)]
        real_rooted: bool,

        /// Check log-concavity
        #[arg(long)]
        log_concave: bool,
    },

    /// For each permutation, compute polynomial over its order ideal
    Scan {
        /// Size of permutations (generates S_n)
        #[arg(long)]
        size: u8,

        #[command(flatten)]
        filters: PermFilters,

        /// Order type for ideal computation
        #[arg(long, value_enum, default_value_t = OrderType::Bruhat)]
        ideal: OrderType,

        /// Statistic to use
        #[arg(long, value_enum, default_value_t = Stat::Des)]
        stat: Stat,

        /// Check real-rootedness of each polynomial
        #[arg(long)]
        real_rooted: bool,

        /// Stop on first non-real-rooted polynomial
        #[arg(long)]
        halt: bool,
    },

    /// List combinatorial objects
    List {
        #[command(flatten)]
        source: Source,
    },

    /// Find a polynomial recurrence P_n(t) for n = min..=max
    Recurrence {
        /// Minimum n (1-based)
        #[arg(long, default_value_t = 1)]
        min_n: u8,

        /// Maximum n
        #[arg(long)]
        max_n: u8,

        /// Statistic to use
        #[arg(long, value_enum, default_value_t = Stat::Des)]
        stat: Stat,

        #[command(flatten)]
        filters: PermFilters,

        /// Automatically search for the simplest recurrence
        #[arg(long)]
        auto: bool,

        /// Show search progress (with --auto)
        #[arg(long)]
        verbose: bool,

        // --- Manual parameters (ignored when --auto is set) ---
        /// Max degree of coefficients in t
        #[arg(long, default_value_t = 1)]
        var_deg: usize,

        /// Max degree of coefficients in n
        #[arg(long, default_value_t = 1)]
        idx_deg: usize,

        /// Max derivative order
        #[arg(long, default_value_t = 0)]
        diff_deg: usize,

        /// Recurrence length (how many previous terms)
        #[arg(long, default_value_t = 2)]
        rec_len: usize,

        /// Allow inhomogeneous recurrence
        #[arg(long)]
        inhomogeneous: bool,

        /// Degree in t of LHS denominator
        #[arg(long, default_value_t = 0)]
        denom_var_deg: usize,

        /// Degree in n of LHS denominator
        #[arg(long, default_value_t = 0)]
        denom_idx_deg: usize,

        // --- Auto-mode upper bounds ---
        /// (--auto) Max recurrence length to search
        #[arg(long, default_value_t = 5)]
        max_rec_len: usize,

        /// (--auto) Max var_deg to search
        #[arg(long, default_value_t = 3)]
        max_var_deg: usize,

        /// (--auto) Max idx_deg to search
        #[arg(long, default_value_t = 3)]
        max_idx_deg: usize,

        /// (--auto) Max diff_deg to search
        #[arg(long, default_value_t = 2)]
        max_diff_deg: usize,

        /// (--auto) Minimum equation surplus
        #[arg(long, default_value_t = 1)]
        min_margin: usize,
    },
}

#[derive(Args)]
struct Source {
    /// Generate all permutations of S_n
    #[arg(long)]
    perms: Option<u8>,

    /// Generate all parking functions of size n
    #[arg(long)]
    parking: Option<u8>,

    /// Bruhat lower ideal of a permutation (e.g., 321 or 3,2,1)
    #[arg(long, value_name = "PERM")]
    bruhat_ideal: Option<String>,

    /// Weak lower ideal of a permutation (e.g., 321 or 3,2,1)
    #[arg(long, value_name = "PERM")]
    weak_ideal: Option<String>,

    /// Generate words on a board (use with --content and --board)
    #[arg(long)]
    words: bool,

    /// Content vector for words (comma-separated, e.g., 2,2,1)
    #[arg(long, value_name = "ALPHA")]
    content: Option<String>,

    /// Board/partition for words (e.g., 22233 or 2,2,2,3,3)
    #[arg(long, value_name = "LAMBDA")]
    board: Option<String>,

    /// Skew lower board (e.g., 11 or 1,1)
    #[arg(long, value_name = "MU")]
    skew: Option<String>,

    #[command(flatten)]
    filters: PermFilters,
}

#[derive(Args)]
struct PermFilters {
    /// Pattern(s) to avoid (e.g., 312). Repeatable.
    #[arg(long, value_name = "PATTERN")]
    avoiding: Vec<String>,

    /// Only alternating permutations (up-down: p1 < p2 > p3 < ...)
    #[arg(long)]
    alternating: bool,

    /// Only derangements (no fixed points)
    #[arg(long)]
    derangement: bool,

    /// Only permutations starting with this value
    #[arg(long, value_name = "VAL")]
    starts_with: Option<u8>,

    /// Only permutations ending with this value
    #[arg(long, value_name = "VAL")]
    ends_with: Option<u8>,
}

fn build_constraints(filters: &PermFilters) -> permutation::PermConstraints {
    permutation::PermConstraints {
        avoiding: filters
            .avoiding
            .iter()
            .map(|s| permutation::parse_sequence(s))
            .collect(),
        derangement: filters.derangement,
        alternating: filters.alternating,
        involution: false, // TODO: add --involution CLI flag
        starts_with: filters.starts_with,
        ends_with: filters.ends_with,
    }
}

fn get_objects(source: &Source) -> Vec<Vec<u8>> {
    if let Some(n) = source.perms {
        let constraints = build_constraints(&source.filters);
        let has_constraints = !constraints.avoiding.is_empty()
            || constraints.derangement
            || constraints.alternating
            || constraints.starts_with.is_some()
            || constraints.ends_with.is_some();
        let perms = if has_constraints {
            permutation::filtered_permutations(n, &constraints)
        } else {
            permutation::all_permutations(n)
        };
        perms
    } else if let Some(n) = source.parking {
        parking::all_parking_functions(n)
    } else if let Some(ref s) = source.bruhat_ideal {
        let perm = permutation::parse_sequence(s);
        order::bruhat_lower_ideal(&perm)
    } else if let Some(ref s) = source.weak_ideal {
        let perm = permutation::parse_sequence(s);
        order::weak_lower_ideal(&perm)
    } else if source.words || source.content.is_some() {
        let content_str = source
            .content
            .as_ref()
            .expect("--content required for words");
        let content = word::parse_content(content_str);
        let n: usize = content.iter().sum();
        let k = content.len() as u8;
        let board = if let Some(ref b) = source.board {
            word::parse_board(b)
        } else {
            vec![k; n]
        };
        let skew = source.skew.as_ref().map(|s| word::parse_board(s));
        word::words_on_board(&content, &board, skew.as_deref())
    } else {
        eprintln!(
            "No source specified. Use --perms, --parking, --bruhat-ideal, --weak-ideal, or --content."
        );
        std::process::exit(1);
    }
}

fn format_obj(p: &[u8]) -> String {
    if p.iter().all(|&x| x < 10) {
        p.iter().map(|x| x.to_string()).collect::<String>()
    } else {
        let strs: Vec<String> = p.iter().map(|x| x.to_string()).collect();
        format!("[{}]", strs.join(","))
    }
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Poly {
            source,
            stat,
            real_rooted,
            log_concave,
        } => {
            let objects = get_objects(&source);
            let coeffs = polynomial_builder::build_generating_polynomial(&objects, stat);

            println!("Objects: {}", objects.len());
            println!("Stat: {}", stat);
            println!("Polynomial: {}", format_poly(&coeffs));
            println!("Coefficients: {:?}", coeffs);

            if real_rooted {
                let rr = is_real_rooted(&coeffs);
                println!("Real-rooted: {}", if rr { "yes" } else { "NO" });
            }
            if log_concave {
                let lc = is_log_concave(&coeffs);
                println!("Log-concave: {}", if lc { "yes" } else { "NO" });
            }
        }

        Commands::Scan {
            size,
            filters,
            ideal,
            stat,
            real_rooted,
            halt,
        } => {
            let constraints = build_constraints(&filters);
            let perms = permutation::filtered_permutations(size, &constraints);

            let total = perms.len();
            eprintln!(
                "Scanning {} permutations (ideal: {:?}, stat: {})",
                total, ideal, stat
            );

            let mut count = 0;
            let mut counterexamples = 0;

            for pi in &perms {
                let ideal_set = match ideal {
                    OrderType::Bruhat => order::bruhat_lower_ideal(pi),
                    OrderType::Weak => order::weak_lower_ideal(pi),
                };
                let coeffs = polynomial_builder::build_generating_polynomial(&ideal_set, stat);

                let mut line = format!("{}\t{:?}", format_obj(pi), coeffs);

                if real_rooted {
                    let rr = is_real_rooted(&coeffs);
                    line.push_str(if rr { "\tRR" } else { "\tNOT-RR" });
                    if !rr {
                        counterexamples += 1;
                        if halt {
                            println!("{}", line);
                            eprintln!("COUNTEREXAMPLE found! Stopping.");
                            std::process::exit(0);
                        }
                    }
                }

                println!("{}", line);
                count += 1;

                if total >= 500 && count % 500 == 0 {
                    eprint!("\r[{}/{}]", count, total);
                    io::stderr().flush().ok();
                }
            }

            if total >= 500 {
                eprintln!();
            }

            if real_rooted {
                eprintln!("---");
                if counterexamples == 0 {
                    eprintln!("All {} polynomials are real-rooted.", total);
                } else {
                    eprintln!(
                        "{} of {} polynomials are NOT real-rooted.",
                        counterexamples, total
                    );
                }
            }
        }

        Commands::List { source } => {
            let objects = get_objects(&source);
            eprintln!("Count: {}", objects.len());
            for obj in &objects {
                println!("{}", format_obj(obj));
            }
        }

        Commands::Recurrence {
            min_n,
            max_n,
            stat,
            filters,
            auto,
            verbose,
            var_deg,
            idx_deg,
            diff_deg,
            rec_len,
            inhomogeneous,
            denom_var_deg,
            denom_idx_deg,
            max_rec_len,
            max_var_deg,
            max_idx_deg,
            max_diff_deg,
            min_margin,
        } => {
            // Compute P_n(t) for each n in the range.
            let mut polys: Vec<Vec<i64>> = Vec::new();
            for n in min_n..=max_n {
                let constraints = build_constraints(&filters);
                let perms = permutation::filtered_permutations(n, &constraints);
                let coeffs = polynomial_builder::build_generating_polynomial(&perms, stat);
                eprintln!(
                    "P_{}(t) = {}  (coeffs: {:?})",
                    n,
                    format_poly(&coeffs),
                    coeffs
                );
                polys.push(coeffs);
            }

            eprintln!("---");

            if auto {
                let search = recurrence::AdaptiveSearchOptions {
                    max_rec_len,
                    max_var_deg,
                    max_idx_deg,
                    max_diff_deg,
                    try_inhomogeneous: inhomogeneous,
                    try_denominator: denom_var_deg > 0 || denom_idx_deg > 0,
                    min_margin,
                    verbose,
                    ..Default::default()
                };

                match recurrence::find_recurrence_adaptive(&polys, &search) {
                    Some(result) => {
                        eprintln!(
                            "Found after {} candidates (unknowns={}, equations={}, \
                             rec_len={}, var_deg={}, idx_deg={}, diff_deg={})",
                            result.candidates_tried,
                            result.num_unknowns,
                            result.num_equations,
                            result.opts.rec_len,
                            result.opts.var_deg,
                            result.opts.idx_deg,
                            result.opts.diff_deg,
                        );
                        println!("{}", result.recurrence);
                    }
                    None => println!("No recurrence found within search bounds."),
                }
            } else {
                let opts = recurrence::RecurrenceOptions {
                    var_deg,
                    idx_deg,
                    diff_deg,
                    rec_len,
                    homogeneous: !inhomogeneous,
                    denom_var_deg,
                    denom_idx_deg,
                };

                match recurrence::find_polynomial_recurrence(&polys, &opts) {
                    Some(rec) => println!("{}", rec),
                    None => println!("No recurrence found with the given parameters."),
                }
            }
        }
    }
}
