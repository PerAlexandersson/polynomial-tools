use std::env;
use std::process;

use combinatoric_core::optimist_sort_step_distribution_via_derangements;

fn parse_max_n() -> usize {
    let mut args = env::args().skip(1);
    let Some(raw) = args.next() else {
        return 10;
    };
    if args.next().is_some() {
        eprintln!("usage: optimist_sorting [max_n]");
        process::exit(2);
    }
    raw.parse::<usize>().unwrap_or_else(|error| {
        eprintln!("invalid max_n `{raw}`: {error}");
        process::exit(2);
    })
}

fn main() {
    let max_n = parse_max_n();
    for n in 1..=max_n {
        let row = optimist_sort_step_distribution_via_derangements(n);
        let row_text = row
            .iter()
            .map(u128::to_string)
            .collect::<Vec<_>>()
            .join(", ");
        println!("{n}: {row_text}");
    }
}
