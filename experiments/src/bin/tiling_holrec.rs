//! Find holonomic (P-recursive) recurrence for tiling polynomials:
//!   sum_{j=0}^{order} a_j(n,t) * P_{n-j}(t) = 0
//! where a_j(n,t) = sum c_{j,p,q} n^p t^q
//!
//! Usage: cargo run --release --bin tiling_holrec -- 3,2,2,1 4

use std::env;

type Poly = Vec<i128>;

fn poly_zero() -> Poly { vec![0] }
fn poly_one() -> Poly { vec![1] }
fn poly_is_zero(p: &Poly) -> bool { p.iter().all(|&c| c == 0) }

fn poly_add(a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let mut r = vec![0i128; n];
    for (i, &c) in a.iter().enumerate() { r[i] += c; }
    for (i, &c) in b.iter().enumerate() { r[i] += c; }
    while r.len() > 1 && *r.last().unwrap() == 0 { r.pop(); }
    r
}

fn poly_mul(a: &Poly, b: &Poly) -> Poly {
    if poly_is_zero(a) || poly_is_zero(b) { return poly_zero(); }
    let n = a.len() + b.len() - 1;
    let mut r = vec![0i128; n];
    for (i, &ca) in a.iter().enumerate() {
        if ca == 0 { continue; }
        for (j, &cb) in b.iter().enumerate() { r[i + j] += ca * cb; }
    }
    while r.len() > 1 && *r.last().unwrap() == 0 { r.pop(); }
    r
}

fn deg(p: &Poly) -> usize {
    let mut d = p.len() - 1;
    while d > 0 && p[d] == 0 { d -= 1; }
    d
}

fn encode(state: &[usize], d: usize) -> usize {
    let base = d + 1; let mut idx = 0;
    for &s in state { idx = idx * base + s; } idx
}
fn decode(mut idx: usize, ell: usize, d: usize) -> Vec<usize> {
    let base = d + 1; let mut s = vec![0usize; ell];
    for i in (0..ell).rev() { s[i] = idx % base; idx /= base; } s
}

fn compute_polys(mu: &[usize], d: usize, max_n: usize) -> Vec<Poly> {
    let ell = mu.len() - 1;
    let ns = (d + 1).pow(ell as u32);
    let mut t_mat = vec![vec![poly_zero(); ns]; ns];
    for j in 0..ns {
        let old = decode(j, ell, d);
        for v in 0..=d {
            let mut ok = true;
            if v > 0 {
                for m in 1..=ell {
                    let prev = old[m - 1];
                    if prev > 0 && v < prev + mu[m] { ok = false; break; }
                }
            }
            if ok {
                let mut new_s = vec![0usize; ell];
                new_s[0] = v;
                for k in 1..ell { new_s[k] = old[k - 1]; }
                let i = encode(&new_s, d);
                let w: Poly = if v > 0 { vec![0, 1] } else { vec![1] };
                t_mat[i][j] = poly_add(&t_mat[i][j], &w);
            }
        }
    }
    let mut vec_s: Vec<Poly> = vec![poly_zero(); ns];
    vec_s[0] = poly_one();
    let mut polys = vec![];
    for _ in 0..=max_n {
        polys.push(vec_s[0].clone());
        let mut nv: Vec<Poly> = vec![poly_zero(); ns];
        for i in 0..ns { for j in 0..ns {
            if poly_is_zero(&t_mat[i][j]) || poly_is_zero(&vec_s[j]) { continue; }
            nv[i] = poly_add(&nv[i], &poly_mul(&t_mat[i][j], &vec_s[j]));
        }}
        vec_s = nv;
    }
    polys
}

/// Search for holonomic recurrence of given order, with a_j(n,t) having
/// degree <= max_n_deg in n and degree <= max_t_deg in t.
fn find_holonomic(polys: &[Poly], order: usize, max_n_deg: usize, max_t_deg: usize) -> bool {
    let n_coeffs_per_aj = (max_n_deg + 1) * (max_t_deg + 1);
    let n_unknowns = (order + 1) * n_coeffs_per_aj;

    // Collect equations
    let max_poly_t_deg = polys.iter().map(|p| deg(p)).max().unwrap_or(0);
    let max_r = max_poly_t_deg + max_t_deg + 2;

    let mut rows: Vec<Vec<f64>> = vec![];

    for n in order..polys.len() {
        let nn = n as f64;
        // n_powers[p] = n^p
        let mut n_powers = vec![1.0f64; max_n_deg + 1];
        for p in 1..=max_n_deg { n_powers[p] = n_powers[p-1] * nn; }

        for r in 0..=max_r {
            let mut row = vec![0.0f64; n_unknowns];
            let mut any_nonzero = false;

            for j in 0..=order {
                let pnj = &polys[n - j];
                for p in 0..=max_n_deg {
                    for q in 0..=max_t_deg {
                        if r >= q {
                            let t_idx = r - q;
                            let pval = if t_idx < pnj.len() { pnj[t_idx] as f64 } else { 0.0 };
                            if pval.abs() > 0.5 { any_nonzero = true; }
                            let var_idx = j * n_coeffs_per_aj + p * (max_t_deg + 1) + q;
                            row[var_idx] += n_powers[p] * pval;
                        }
                    }
                }
            }
            if any_nonzero { rows.push(row); }
        }
    }

    // Gaussian elimination
    let nr = rows.len();
    let nc = n_unknowns;
    let mut mat = rows;
    let mut pivot_row = 0;
    let mut pivot_col = 0;
    let mut pivots = vec![];

    while pivot_row < nr && pivot_col < nc {
        let mut max_val = 0.0f64;
        let mut max_r = pivot_row;
        for i in pivot_row..nr {
            if mat[i][pivot_col].abs() > max_val { max_val = mat[i][pivot_col].abs(); max_r = i; }
        }
        if max_val < 1e-8 { pivot_col += 1; continue; }
        mat.swap(pivot_row, max_r);
        let scale = mat[pivot_row][pivot_col];
        for c in 0..nc { mat[pivot_row][c] /= scale; }
        for i in 0..nr {
            if i == pivot_row { continue; }
            let f = mat[i][pivot_col];
            if f.abs() < 1e-12 { continue; }
            for c in 0..nc { mat[i][c] -= f * mat[pivot_row][c]; }
        }
        pivots.push(pivot_col);
        pivot_row += 1;
        pivot_col += 1;
    }

    let rank = pivots.len();
    let null_dim = nc.saturating_sub(rank);

    if null_dim == 0 { return false; }

    // Extract null vector
    let free_cols: Vec<usize> = (0..nc).filter(|c| !pivots.contains(c)).collect();

    for &fc in &free_cols {
        let mut sol = vec![0.0f64; nc];
        sol[fc] = 1.0;
        for (&pc, pr) in pivots.iter().zip(0..rank).rev() {
            let mut val = 0.0;
            for c in (pc+1)..nc { val += mat[pr][c] * sol[c]; }
            sol[pc] = -val;
        }

        // Scale to integers
        let max_abs = sol.iter().map(|&x| x.abs()).fold(0.0f64, f64::max);
        if max_abs < 1e-10 { continue; }
        let min_nz = sol.iter().map(|&x| x.abs()).filter(|&x| x > 1e-8).fold(f64::MAX, f64::min);
        let scaled: Vec<f64> = sol.iter().map(|&x| x / min_nz).collect();

        // Try to round to small integers
        let int_sol: Vec<i64> = scaled.iter().map(|&x| x.round() as i64).collect();

        // Check residual
        let mut max_resid = 0.0f64;
        for s in &scaled {
            let r = (s - s.round()).abs();
            if r > max_resid { max_resid = r; }
        }

        if max_resid > 0.01 { continue; } // not integer solution

        // Verify exactly using i128
        let mut verified = true;
        for n in order..polys.len() {
            let mut n_powers = vec![1i128; max_n_deg + 1];
            for p in 1..=max_n_deg { n_powers[p] = n_powers[p-1] * (n as i128); }

            // For each t-power, check the sum is 0
            let max_r_check = polys.iter().take(n+1).map(|p| deg(p)).max().unwrap_or(0) + max_t_deg;
            for r in 0..=max_r_check {
                let mut total = 0i128;
                for j in 0..=order {
                    let pnj = &polys[n - j];
                    for p in 0..=max_n_deg {
                        for q in 0..=max_t_deg {
                            if r >= q {
                                let t_idx = r - q;
                                let pval = if t_idx < pnj.len() { pnj[t_idx] } else { 0 };
                                let var_idx = j * n_coeffs_per_aj + p * (max_t_deg + 1) + q;
                                total += int_sol[var_idx] as i128 * n_powers[p] * pval;
                            }
                        }
                    }
                }
                if total != 0 { verified = false; break; }
            }
            if !verified { break; }
        }

        if verified {
            println!("  FOUND! order={}, n-deg={}, t-deg={}", order, max_n_deg, max_t_deg);
            for j in 0..=order {
                let mut terms = vec![];
                for p in 0..=max_n_deg {
                    for q in 0..=max_t_deg {
                        let idx = j * n_coeffs_per_aj + p * (max_t_deg + 1) + q;
                        let c = int_sol[idx];
                        if c == 0 { continue; }
                        let npart = match p { 0 => "".into(), 1 => "n".into(), _ => format!("n^{}", p) };
                        let tpart = match q { 0 => "".into(), 1 => "t".into(), _ => format!("t^{}", q) };
                        let var = match (p, q) {
                            (0, 0) => format!("{}", c),
                            (0, _) => if c == 1 { tpart } else if c == -1 { format!("-{}", tpart) } else { format!("{}{}", c, tpart) },
                            (_, 0) => if c == 1 { npart } else if c == -1 { format!("-{}", npart) } else { format!("{}{}", c, npart) },
                            _ => if c == 1 { format!("{}{}", npart, tpart) } else if c == -1 { format!("-{}{}", npart, tpart) } else { format!("{}{}{}",c, npart, tpart) },
                        };
                        terms.push(var);
                    }
                }
                if !terms.is_empty() {
                    println!("    a_{}(n,t) = {}", j, terms.join(" + "));
                }
            }
            // Print nicely
            print!("  Recurrence: ");
            let mut first = true;
            for j in 0..=order {
                let mut any = false;
                for p in 0..=max_n_deg { for q in 0..=max_t_deg {
                    let idx = j * n_coeffs_per_aj + p * (max_t_deg + 1) + q;
                    if int_sol[idx] != 0 { any = true; }
                }}
                if !any { continue; }
                if !first { print!(" + "); }
                print!("a_{}(n,t)*P(n-{})", j, j);
                first = false;
            }
            println!(" = 0");
            return true;
        }
    }
    false
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mu: Vec<usize> = args[0].split(',').map(|s| s.trim().parse().unwrap()).collect();
    let d: usize = args[1].parse().unwrap();
    let max_n: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(40);

    println!("mu = {:?}, d = {}", mu, d);
    let polys = compute_polys(&mu, d, max_n);

    println!("Generated {} polynomials, max deg = {}", polys.len(), polys.iter().map(|p| deg(p)).max().unwrap_or(0));

    println!("\nSearching for holonomic recurrence...");
    // Try increasing order and coefficient complexity
    for order in 1..=6 {
        for max_n_deg in 0..=3 {
            for max_t_deg in 0..=4 {
                let n_unk = (order + 1) * (max_n_deg + 1) * (max_t_deg + 1);
                if n_unk > 200 { continue; } // too many unknowns
                if find_holonomic(&polys, order, max_n_deg, max_t_deg) {
                    return;
                }
            }
        }
        println!("  No recurrence of order {} found", order);
    }
    println!("No short holonomic recurrence found.");
}
