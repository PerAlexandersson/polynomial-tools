/// Quick recurrence search for H_n(t), using fewer terms to keep it fast.
use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};

fn main() {
    // Use first 14 terms (enough to find, rest to verify)
    let polys: Vec<Vec<i64>> = vec![
        vec![1],
        vec![1],
        vec![1, 1],
        vec![1, 3, 1],
        vec![1, 7, 7, 1],
        vec![1, 14, 31, 14, 1],
        vec![1, 26, 109, 109, 26, 1],
        vec![1, 46, 334, 623, 334, 46, 1],
        vec![1, 79, 937, 2951, 2951, 937, 79, 1],
        vec![1, 133, 2475, 12331, 20641, 12331, 2475, 133, 1],
        vec![1, 221, 6267, 47191, 123216, 123216, 47191, 6267, 221, 1],
        vec![
            1, 364, 15393, 169416, 656683, 1019051, 656683, 169416, 15393, 364, 1,
        ],
        vec![
            1, 596, 36976, 579889, 3217526, 7349140, 7349140, 3217526, 579889, 36976, 596, 1,
        ],
        vec![
            1, 972, 87369, 1914226, 14786816, 47816612, 70148989, 47816612, 14786816, 1914226,
            87369, 972, 1,
        ],
    ];

    let arg = std::env::args().nth(1).unwrap_or("1".to_string());
    let search_id: usize = arg.parse().unwrap_or(1);

    let (desc, opts) = match search_id {
        1 => (
            "len=2 vd=2 id=2 d=0",
            AdaptiveSearchOptions {
                max_rec_len: 2,
                max_var_deg: 2,
                max_idx_deg: 2,
                max_diff_deg: 0,
                try_inhomogeneous: false,
                try_denominator: false,
                max_denom_var_deg: 0,
                max_denom_idx_deg: 0,
                min_margin: 2,
                verbose: false,
            },
        ),
        2 => (
            "len=2 vd=2 id=2 d=1",
            AdaptiveSearchOptions {
                max_rec_len: 2,
                max_var_deg: 2,
                max_idx_deg: 2,
                max_diff_deg: 1,
                try_inhomogeneous: false,
                try_denominator: false,
                max_denom_var_deg: 0,
                max_denom_idx_deg: 0,
                min_margin: 2,
                verbose: false,
            },
        ),
        3 => (
            "len=3 vd=1 id=1 d=0",
            AdaptiveSearchOptions {
                max_rec_len: 3,
                max_var_deg: 1,
                max_idx_deg: 1,
                max_diff_deg: 0,
                try_inhomogeneous: false,
                try_denominator: false,
                max_denom_var_deg: 0,
                max_denom_idx_deg: 0,
                min_margin: 2,
                verbose: false,
            },
        ),
        4 => (
            "len=3 vd=2 id=2 d=0",
            AdaptiveSearchOptions {
                max_rec_len: 3,
                max_var_deg: 2,
                max_idx_deg: 2,
                max_diff_deg: 0,
                try_inhomogeneous: false,
                try_denominator: false,
                max_denom_var_deg: 0,
                max_denom_idx_deg: 0,
                min_margin: 2,
                verbose: false,
            },
        ),
        5 => (
            "len=2 vd=2 id=2 d=2",
            AdaptiveSearchOptions {
                max_rec_len: 2,
                max_var_deg: 2,
                max_idx_deg: 2,
                max_diff_deg: 2,
                try_inhomogeneous: false,
                try_denominator: false,
                max_denom_var_deg: 0,
                max_denom_idx_deg: 0,
                min_margin: 2,
                verbose: false,
            },
        ),
        6 => (
            "len=3 vd=1 id=1 d=1",
            AdaptiveSearchOptions {
                max_rec_len: 3,
                max_var_deg: 1,
                max_idx_deg: 1,
                max_diff_deg: 1,
                try_inhomogeneous: false,
                try_denominator: false,
                max_denom_var_deg: 0,
                max_denom_idx_deg: 0,
                min_margin: 2,
                verbose: false,
            },
        ),
        7 => (
            "len=2 vd=3 id=2 d=1",
            AdaptiveSearchOptions {
                max_rec_len: 2,
                max_var_deg: 3,
                max_idx_deg: 2,
                max_diff_deg: 1,
                try_inhomogeneous: false,
                try_denominator: false,
                max_denom_var_deg: 0,
                max_denom_idx_deg: 0,
                min_margin: 2,
                verbose: false,
            },
        ),
        8 => (
            "len=3 vd=2 id=1 d=1",
            AdaptiveSearchOptions {
                max_rec_len: 3,
                max_var_deg: 2,
                max_idx_deg: 1,
                max_diff_deg: 1,
                try_inhomogeneous: false,
                try_denominator: false,
                max_denom_var_deg: 0,
                max_denom_idx_deg: 0,
                min_margin: 2,
                verbose: false,
            },
        ),
        _ => return,
    };

    println!("Search {}: {}", search_id, desc);
    let t = std::time::Instant::now();
    match find_recurrence_adaptive(&polys, &opts) {
        Some(res) => {
            let s = format!("{}", res.recurrence);
            if s.len() > 500 {
                println!("  overfitted ({} chars) {:?}", s.len(), t.elapsed());
            } else {
                println!("  FOUND {:?}\n  {}", t.elapsed(), s);
            }
        }
        None => println!("  none {:?}", t.elapsed()),
    }
}
