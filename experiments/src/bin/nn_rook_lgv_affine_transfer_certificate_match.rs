//! Check the factorized affine transfer target.
//!
//! The raw affine branch model is too coarse for the LGV involution.  The
//! replacement target is a path-word/transfer stack whose path set is the
//! disjoint union
//!
//!   A paths, t-marked Q paths, s-marked shifted Q paths.
//!
//! This binary checks, with independent row-weight variables and independent
//! affine parameters `s,t`, that this factorized transfer construction matches
//! the determinant packet formula for
//!
//!   A(x) + (s*x + t) Q(x)
//!
//! for both `Q=U` and `Q=L`.

use std::collections::{BTreeMap, HashMap};
use std::env;

use experiments::nn_rook_utils::{fixed_row_count_lgv, partitions, subset_size};

type Exponents = Vec<u8>;
type MultiPoly = BTreeMap<Exponents, i128>;

#[derive(Clone)]
struct Packet {
    c: usize,
    p: usize,
    h: usize,
    a: Vec<MultiPoly>,
    u: Vec<MultiPoly>,
    l: Vec<MultiPoly>,
}

#[derive(Clone, Copy)]
enum Param {
    S,
    T,
}

fn zero() -> MultiPoly {
    BTreeMap::new()
}

fn one(num_rows: usize) -> MultiPoly {
    BTreeMap::from([(vec![0; num_rows + 2], 1)])
}

fn monom(num_rows: usize, coeff: i128, row_mask: u64, s_exp: u8, t_exp: u8) -> MultiPoly {
    if coeff == 0 {
        return zero();
    }
    let mut exps = vec![0; num_rows + 2];
    for row in 0..num_rows {
        if (row_mask & (1u64 << row)) != 0 {
            exps[row] = 1;
        }
    }
    exps[num_rows] = s_exp;
    exps[num_rows + 1] = t_exp;
    BTreeMap::from([(exps, coeff)])
}

fn var(num_rows: usize, row: Option<usize>, param: Option<Param>) -> MultiPoly {
    let mut exps = vec![0; num_rows + 2];
    if let Some(row) = row {
        exps[row] = 1;
    }
    match param {
        Some(Param::S) => exps[num_rows] = 1,
        Some(Param::T) => exps[num_rows + 1] = 1,
        None => {}
    }
    BTreeMap::from([(exps, 1)])
}

fn add_to(dst: &mut MultiPoly, src: &MultiPoly) {
    for (key, value) in src {
        let entry = dst.entry(key.clone()).or_insert(0);
        *entry += value;
        if *entry == 0 {
            dst.remove(key);
        }
    }
}

fn mul(a: &MultiPoly, b: &MultiPoly) -> MultiPoly {
    let mut out = zero();
    for (ea, ca) in a {
        for (eb, cb) in b {
            let mut exps = ea.clone();
            for (i, e) in eb.iter().enumerate() {
                exps[i] = exps[i].checked_add(*e).expect("exponent overflow");
            }
            *out.entry(exps).or_insert(0) += ca * cb;
        }
    }
    out.retain(|_, c| *c != 0);
    out
}

fn scale_var(
    p: &MultiPoly,
    num_rows: usize,
    row: Option<usize>,
    param: Option<Param>,
) -> MultiPoly {
    mul(p, &var(num_rows, row, param))
}

fn trim_sequence(mut seq: Vec<MultiPoly>) -> Vec<MultiPoly> {
    while seq.len() > 1 && seq.last().is_some_and(|p| p.is_empty()) {
        seq.pop();
    }
    seq
}

fn add_coeff(seq: &mut Vec<MultiPoly>, degree: usize, term: MultiPoly) {
    if degree >= seq.len() {
        seq.resize_with(degree + 1, zero);
    }
    add_to(&mut seq[degree], &term);
}

fn seq_add(a: &[MultiPoly], b: &[MultiPoly]) -> Vec<MultiPoly> {
    let len = a.len().max(b.len());
    let mut out = vec![zero(); len];
    for i in 0..len {
        if let Some(p) = a.get(i) {
            add_to(&mut out[i], p);
        }
        if let Some(p) = b.get(i) {
            add_to(&mut out[i], p);
        }
    }
    trim_sequence(out)
}

fn seq_shift_degree(seq: &[MultiPoly]) -> Vec<MultiPoly> {
    let mut out = vec![zero(); seq.len() + 1];
    for (i, p) in seq.iter().enumerate() {
        add_to(&mut out[i + 1], p);
    }
    trim_sequence(out)
}

fn seq_param(seq: &[MultiPoly], num_rows: usize, param: Param) -> Vec<MultiPoly> {
    trim_sequence(
        seq.iter()
            .map(|p| scale_var(p, num_rows, None, Some(param)))
            .collect(),
    )
}

fn fixed_count_cached(
    cache: &mut HashMap<(usize, u64), i64>,
    eta: &[usize],
    strip: usize,
    mask: u64,
) -> i64 {
    if let Some(&value) = cache.get(&(strip, mask)) {
        return value;
    }
    let value = fixed_row_count_lgv(eta, strip, mask);
    cache.insert((strip, mask), value);
    value
}

fn bottom_row_bucket(mask: u64) -> usize {
    if mask == 0 {
        0
    } else {
        63 - mask.leading_zeros() as usize
    }
}

fn fixed_remaining_packet(
    cache: &mut HashMap<(usize, u64), i64>,
    eta: &[usize],
    remaining_mask: u64,
    c: usize,
    p: usize,
) -> (Vec<MultiPoly>, Vec<MultiPoly>, Vec<MultiPoly>) {
    let num_rows = eta.len();
    let k = subset_size(remaining_mask);
    let d_c1 = fixed_count_cached(cache, eta, c + 1, remaining_mask) as i128;
    let d_p1 = fixed_count_cached(cache, eta, p + 1, remaining_mask) as i128;
    let tail_sum: i128 = ((c + 2)..=(p + 1))
        .map(|j| fixed_count_cached(cache, eta, j, remaining_mask) as i128)
        .sum();

    let mut a = vec![zero(); k + 2];
    add_coeff(&mut a, k, monom(num_rows, d_c1, remaining_mask, 0, 0));
    add_coeff(
        &mut a,
        k + 1,
        monom(num_rows, tail_sum, remaining_mask, 0, 0),
    );

    let mut u = vec![zero(); k + 1];
    let mut l = vec![zero(); k + 1];
    add_coeff(
        &mut u,
        k,
        monom(num_rows, d_c1 + d_p1, remaining_mask, 0, 0),
    );
    add_coeff(
        &mut l,
        k,
        monom(num_rows, d_c1 - d_p1, remaining_mask, 0, 0),
    );

    let highest_remaining_row = (0..eta.len())
        .filter(|&r| (remaining_mask & (1u64 << r)) != 0)
        .max();
    let first_forced_row = highest_remaining_row.map_or(0, |r| r + 1);
    for r in first_forced_row..eta.len() {
        if eta[r] > c {
            let forced_mask = remaining_mask | (1u64 << r);
            let term = monom(num_rows, d_c1, forced_mask, 0, 0);
            add_coeff(&mut u, k, term.clone());
            add_coeff(&mut l, k, term);
        }
    }

    (trim_sequence(a), trim_sequence(u), trim_sequence(l))
}

fn determinant_cutoff_packets(eta: &[usize]) -> Vec<Packet> {
    let row_add_width = eta[eta.len() - 1];
    let all_masks: Vec<u64> = (0..(1u64 << eta.len())).collect();
    let mut count_cache: HashMap<(usize, u64), i64> = HashMap::new();
    let mut out = Vec::new();

    for c in 0..row_add_width {
        for p in (c + 1)..row_add_width {
            let mut by_bottom = vec![(vec![zero()], vec![zero()], vec![zero()]); eta.len()];
            for &mask in &all_masks {
                let h = bottom_row_bucket(mask);
                let (a_j, u_j, l_j) = fixed_remaining_packet(&mut count_cache, eta, mask, c, p);
                by_bottom[h].0 = seq_add(&by_bottom[h].0, &a_j);
                by_bottom[h].1 = seq_add(&by_bottom[h].1, &u_j);
                by_bottom[h].2 = seq_add(&by_bottom[h].2, &l_j);
            }

            let mut a_sum = vec![zero()];
            let mut u_sum = vec![zero()];
            let mut l_sum = vec![zero()];
            for (h, (a_h, u_h, l_h)) in by_bottom.iter().enumerate() {
                a_sum = seq_add(&a_sum, a_h);
                u_sum = seq_add(&u_sum, u_h);
                l_sum = seq_add(&l_sum, l_h);
                out.push(Packet {
                    c,
                    p,
                    h,
                    a: a_sum.clone(),
                    u: u_sum.clone(),
                    l: l_sum.clone(),
                });
            }
        }
    }
    out
}

fn strip_transfer(eta: &[usize], h: usize, strip: usize) -> Vec<MultiPoly> {
    let num_rows = eta.len();
    let max_width = eta
        .iter()
        .map(|&width| width.saturating_sub(strip))
        .max()
        .unwrap_or(0);
    let start_bound = max_width + 1;
    let mut states: HashMap<(usize, usize), MultiPoly> = HashMap::new();
    states.insert((start_bound, 0), one(num_rows));

    for (row, &row_width) in eta.iter().enumerate() {
        let width = row_width.saturating_sub(strip);
        let mut next: HashMap<(usize, usize), MultiPoly> = HashMap::new();
        for (&(bound, degree), weight) in &states {
            add_to(next.entry((bound, degree)).or_insert_with(zero), weight);
            if row <= h && width > 0 {
                for q in 1..=width.min(bound.saturating_sub(1)) {
                    let selected = scale_var(weight, num_rows, Some(row), None);
                    add_to(next.entry((q, degree + 1)).or_insert_with(zero), &selected);
                }
            }
        }
        states = next;
    }

    let mut out = vec![zero()];
    for ((_, degree), weight) in states {
        add_coeff(&mut out, degree, weight);
    }
    trim_sequence(out)
}

fn reservoir_transfer(eta: &[usize], h: usize, c: usize) -> Vec<MultiPoly> {
    let num_rows = eta.len();
    let strip = c + 1;
    let max_width = eta
        .iter()
        .map(|&width| width.saturating_sub(strip))
        .max()
        .unwrap_or(0);
    let start_bound = max_width + 1;
    let mut active: HashMap<(usize, usize), MultiPoly> = HashMap::new();
    let mut terminal: HashMap<(usize, usize), MultiPoly> = HashMap::new();
    active.insert((start_bound, 0), one(num_rows));

    for (row, &row_width) in eta.iter().enumerate() {
        let width = row_width.saturating_sub(strip);
        let mut next_active: HashMap<(usize, usize), MultiPoly> = HashMap::new();
        let mut next_terminal = terminal.clone();

        for (&(bound, degree), weight) in &active {
            add_to(
                next_active.entry((bound, degree)).or_insert_with(zero),
                weight,
            );
            if row_width > c {
                let reservoir = scale_var(weight, num_rows, Some(row), None);
                add_to(
                    next_terminal.entry((bound, degree)).or_insert_with(zero),
                    &reservoir,
                );
            }
            if row <= h && width > 0 {
                for q in 1..=width.min(bound.saturating_sub(1)) {
                    let selected = scale_var(weight, num_rows, Some(row), None);
                    add_to(
                        next_active.entry((q, degree + 1)).or_insert_with(zero),
                        &selected,
                    );
                }
            }
        }

        active = next_active;
        terminal = next_terminal;
    }

    let mut out = vec![zero()];
    for ((_, degree), weight) in terminal {
        add_coeff(&mut out, degree, weight);
    }
    trim_sequence(out)
}

fn window_transfer(eta: &[usize], h: usize, c: usize, p: usize) -> Vec<MultiPoly> {
    let num_rows = eta.len();
    let strip = c + 1;
    let window_width = p - c;
    let max_width = eta
        .iter()
        .map(|&width| width.saturating_sub(strip))
        .max()
        .unwrap_or(0);
    let start_bound = max_width + 1;
    let mut unseen: HashMap<(usize, usize), MultiPoly> = HashMap::new();
    let mut seen: HashMap<(usize, usize), MultiPoly> = HashMap::new();
    unseen.insert((start_bound, 0), one(num_rows));

    for (row, &row_width) in eta.iter().enumerate() {
        let width = row_width.saturating_sub(strip);
        let mut next_unseen: HashMap<(usize, usize), MultiPoly> = HashMap::new();
        let mut next_seen: HashMap<(usize, usize), MultiPoly> = HashMap::new();

        for (&(bound, degree), weight) in &unseen {
            add_to(
                next_unseen.entry((bound, degree)).or_insert_with(zero),
                weight,
            );
            if row <= h && width > 0 {
                for q in 1..=width.min(bound.saturating_sub(1)) {
                    let selected = scale_var(weight, num_rows, Some(row), None);
                    if q <= window_width {
                        add_to(
                            next_seen.entry((q, degree + 1)).or_insert_with(zero),
                            &selected,
                        );
                    } else {
                        add_to(
                            next_unseen.entry((q, degree + 1)).or_insert_with(zero),
                            &selected,
                        );
                    }
                }
            }
        }
        for (&(bound, degree), weight) in &seen {
            add_to(
                next_seen.entry((bound, degree)).or_insert_with(zero),
                weight,
            );
            if row <= h && width > 0 {
                for q in 1..=width.min(bound.saturating_sub(1)) {
                    let selected = scale_var(weight, num_rows, Some(row), None);
                    add_to(
                        next_seen.entry((q, degree + 1)).or_insert_with(zero),
                        &selected,
                    );
                }
            }
        }

        unseen = next_unseen;
        seen = next_seen;
    }

    let mut out = vec![zero()];
    for ((_, degree), weight) in seen {
        add_coeff(&mut out, degree, weight);
    }
    trim_sequence(out)
}

fn transfer_components(
    eta: &[usize],
    c: usize,
    p: usize,
    h: usize,
) -> (Vec<MultiPoly>, Vec<MultiPoly>, Vec<MultiPoly>) {
    let base_c1 = strip_transfer(eta, h, c + 1);
    let endpoint = strip_transfer(eta, h, p + 1);
    let mut a = base_c1.clone();
    for strip in (c + 2)..=(p + 1) {
        a = seq_add(&a, &seq_shift_degree(&strip_transfer(eta, h, strip)));
    }

    let reservoir = reservoir_transfer(eta, h, c);
    let u = seq_add(&seq_add(&base_c1, &endpoint), &reservoir);
    let l = seq_add(&reservoir, &window_transfer(eta, h, c, p));
    (a, u, l)
}

fn affine_pencil(a: &[MultiPoly], q: &[MultiPoly], num_rows: usize) -> Vec<MultiPoly> {
    let t_q = seq_param(q, num_rows, Param::T);
    let s_x_q = seq_shift_degree(&seq_param(q, num_rows, Param::S));
    seq_add(a, &seq_add(&t_q, &s_x_q))
}

fn format_poly(p: &[MultiPoly]) -> String {
    p.iter()
        .enumerate()
        .filter(|(_, coeff)| !coeff.is_empty())
        .map(|(i, coeff)| format!("deg {i}: {} terms", coeff.len()))
        .collect::<Vec<_>>()
        .join(", ")
}

fn fail(
    label: &str,
    eta: &[usize],
    packet: &Packet,
    direct: &[MultiPoly],
    transfer: &[MultiPoly],
) -> ! {
    println!("FAIL affine transfer certificate match:");
    println!("  label={label}");
    println!("  eta={eta:?}");
    println!("  c={} p={} h={}", packet.c, packet.p, packet.h);
    println!("  direct={}", format_poly(direct));
    println!("  transfer={}", format_poly(transfer));
    for i in 0..direct.len().max(transfer.len()) {
        let d = direct.get(i).cloned().unwrap_or_default();
        let t = transfer.get(i).cloned().unwrap_or_default();
        if d != t {
            println!("  first differing degree={i}");
            println!("  direct terms={d:?}");
            println!("  transfer terms={t:?}");
            break;
        }
    }
    std::process::exit(1);
}

fn main() {
    let max_n = env::args()
        .nth(1)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(10);

    println!("=== Affine LGV transfer-certificate algebra check ===");
    println!("Checking |eta| <= {max_n}.\n");

    let mut packets_checked = 0usize;
    for n in 1..=max_n {
        for eta in partitions(n) {
            if eta.len() >= 20 || eta[eta.len() - 1] < 2 {
                continue;
            }
            for packet in determinant_cutoff_packets(&eta) {
                packets_checked += 1;
                let (transfer_a, transfer_u, transfer_l) =
                    transfer_components(&eta, packet.c, packet.p, packet.h);

                let direct_u = affine_pencil(&packet.a, &packet.u, eta.len());
                let direct_l = affine_pencil(&packet.a, &packet.l, eta.len());
                let factorized_u = affine_pencil(&transfer_a, &transfer_u, eta.len());
                let factorized_l = affine_pencil(&transfer_a, &transfer_l, eta.len());

                if direct_u != factorized_u {
                    fail("U", &eta, &packet, &direct_u, &factorized_u);
                }
                if direct_l != factorized_l {
                    fail("L", &eta, &packet, &direct_l, &factorized_l);
                }
            }
        }
    }

    println!("packets checked: {packets_checked}");
    println!("All affine factorized transfer weights matched the determinant packet.");
}
