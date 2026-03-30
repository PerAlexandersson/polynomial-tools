use polynomial_tools::recurrence::{find_recurrence_adaptive, AdaptiveSearchOptions};
use polynomial_tools::sequences::narayana_polynomials;

fn main() {
    let np = narayana_polynomials(14);
    let polys: Vec<Vec<i64>> = np.into_iter().map(|mut p| {
        while p.len() > 1 && *p.last().unwrap() == 0 { p.pop(); }
        p
    }).collect();

    for (i, p) in polys.iter().enumerate() {
        eprintln!("N_{}: {:?}", i+1, p);
    }

    let search = AdaptiveSearchOptions {
        max_rec_len: 3,
        max_var_deg: 2,
        max_idx_deg: 2,
        max_diff_deg: 1,
        verbose: true,
        ..Default::default()
    };

    eprintln!("\nSearching...");
    match find_recurrence_adaptive(&polys, &search) {
        Some(res) => {
            println!("\nFound: {}", res.recurrence);
            println!("Tried: {} candidates", res.candidates_tried);
        }
        None => println!("\nNo recurrence found"),
    }
}
