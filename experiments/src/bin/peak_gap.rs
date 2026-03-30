//! Targeted investigation of the DU proof gap.
//!
//! Tests potential strengthenings of the induction hypothesis:
//!   (e1) A_j ≼ W_l for ALL j, l  (unrestricted AW)
//!   (e2) D_j ≼ W_l for ALL j, l  (DW universal)
//!   (e3) D_j ≼ D_l for j ≤ l     (DD forward)
//!   (e4) D_l ≼ D_j for j < l     (DD reversed)
//!   (e5) D_j ≼ tD_l for j < l    (equivalent to e4 by Wagner)
//!   (e6) U_j ≼ U_l for j > l     (UU reversed — expected to FAIL)
//!
//! And at the next level:
//!   (f1) A_j^+ ≼ W_l^+ for ALL j, l  (unrestricted AW at λ+)
//!   (f2) D_j^+ ≼ U_l^+ for ALL j, l  (DU at λ+, i.e. S_j ≼ T_l)

use std::collections::BTreeSet;

fn exact_interlaces(fc: &[i64], gc: &[i64]) -> bool {
    polynomial_tools::check_weak_interlacing(fc, gc).unwrap_or(false)
}

fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> {
    let n = perm.len();
    let mut vis: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut q: BTreeSet<Vec<u8>> = BTreeSet::new();
    q.insert(perm.to_vec());
    while let Some(cur) = q.pop_last() {
        for i in 0..n { for j in i+1..n { if cur[i]>cur[j] {
            let mut c=cur.clone(); c.swap(i,j);
            if !vis.contains(&c) { q.insert(c); }
        }}}
        vis.insert(cur);
    }
    vis.into_iter().collect()
}

fn board_to_perm(b: &[u8]) -> Vec<u8> {
    let n=b.len(); let mut p=vec![0u8;n]; let mut u=vec![false;n+1];
    for i in 0..n { for c in (1..=(b[i] as usize).min(n)).rev() {
        if !u[c] { p[i]=c as u8; u[c]=true; break; }
    }} p
}

fn peaks(w: &[u8]) -> usize {
    if w.len()<3 { return 0; }
    (1..w.len()-1).filter(|&i| w[i-1]<w[i] && w[i]>w[i+1]).count()
}

fn pt(p:&[i64])->Vec<i64>{ let mut v=p.to_vec(); while v.len()>1&&*v.last().unwrap()==0{v.pop();} v }
fn pz(p:&[i64])->bool{ p.iter().all(|&c|c==0) }
fn pa(a:&[i64],b:&[i64])->Vec<i64>{
    let l=a.len().max(b.len()); let mut r=vec![0i64;l];
    for(i,&v)in a.iter().enumerate(){r[i]+=v;}
    for(i,&v)in b.iter().enumerate(){r[i]+=v;} pt(&r)
}
fn pmt(p:&[i64])->Vec<i64>{
    let mut r=vec![0i64;p.len()+1];
    for(i,&v)in p.iter().enumerate(){r[i+1]=v;} pt(&r)
}

fn gen_boards(n:usize)->Vec<Vec<u8>>{
    let mut r=vec![]; let mut c=vec![];
    gb(n,n,0,&mut c,&mut r); r
}
fn gb(n:usize,mx:usize,d:usize,c:&mut Vec<u8>,r:&mut Vec<Vec<u8>>){
    if d==n{r.push(c.clone());return;}
    for v in (d+1).max(if d>0{c[d-1]as usize}else{1})..=mx {
        c.push(v as u8); gb(n,mx,d+1,c,r); c.pop();
    }
}

struct BD {
    m: usize,
    d: Vec<Vec<i64>>,
    u: Vec<Vec<i64>>,
    a: Vec<Vec<i64>>,
    w: Vec<Vec<i64>>,
}

fn compute(board: &[u8]) -> BD {
    let n = board.len();
    let m = (board[0] as usize).min(n);
    let perm = board_to_perm(board);
    let ideal = bruhat_lower_ideal(&perm);
    let mut d = vec![vec![0i64]; m+1];
    let mut u = vec![vec![0i64]; m+1];
    for pi in &ideal {
        if pi.len()<2 { let k=pi[0] as usize; if k<=m { while u[k].len()<1{u[k].push(0);} u[k][0]+=1; } continue; }
        let k = pi[0] as usize; if k>m { continue; }
        let pk = peaks(pi);
        let poly = if pi[0]>pi[1] { &mut d[k] } else { &mut u[k] };
        while poly.len()<=pk { poly.push(0); } poly[pk]+=1;
    }
    let mut a = vec![vec![0i64]; m+1];
    let mut w = vec![vec![0i64]; m+1];
    for k in 1..=m { a[k]=pa(&d[k],&u[k]); w[k]=pa(&pmt(&d[k]),&u[k]); }
    BD { m, d, u, a, w }
}

fn main() {
    let max_n: usize = std::env::args().nth(1).and_then(|s|s.parse().ok()).unwrap_or(6);
    println!("=== DU proof gap investigation: n ≤ {} ===\n", max_n);

    // Counters: [tests, fails]
    let mut aw_all = [0u64; 2];   // A_j ≼ W_l for ALL j,l
    let mut dw_all = [0u64; 2];   // D_j ≼ W_l for ALL j,l
    let mut dd_fwd = [0u64; 2];   // D_j ≼ D_l for j ≤ l
    let mut dd_rev = [0u64; 2];   // D_l ≼ D_j for j < l
    let mut uu_rev = [0u64; 2];   // U_j ≼ U_l for j > l (REVERSE)
    let mut aw_plus = [0u64; 2];  // A_j^+ ≼ W_l^+ for ALL j,l
    let mut st_all = [0u64; 2];   // S_j ≼ T_l for ALL j,l

    for n in 1..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let bd = compute(board);
            let m = bd.m;

            // Build S_k, T_k
            let mut s = vec![vec![0i64]; m+2];
            let mut t = vec![vec![0i64]; m+2];
            for k in 2..=m+1 { s[k] = pa(&s[k-1], &bd.a[k-1]); }
            for k in (1..=m).rev() { t[k] = pa(&t[k+1], &bd.w[k]); }

            // (e1) A_j ≼ W_l for ALL j, l
            for j in 1..=m { for l in 1..=m {
                if pz(&bd.a[j]) || pz(&bd.w[l]) { continue; }
                aw_all[0] += 1;
                if !exact_interlaces(&bd.a[j], &bd.w[l]) { aw_all[1] += 1; }
            }}

            // (e2) D_j ≼ W_l for ALL j, l
            for j in 1..=m { for l in 1..=m {
                if pz(&bd.d[j]) || pz(&bd.w[l]) { continue; }
                dw_all[0] += 1;
                if !exact_interlaces(&bd.d[j], &bd.w[l]) { dw_all[1] += 1; }
            }}

            // (e3) D_j ≼ D_l for j ≤ l (forward)
            for j in 1..=m { for l in j+1..=m {
                if pz(&bd.d[j]) || pz(&bd.d[l]) { continue; }
                dd_fwd[0] += 1;
                if !exact_interlaces(&bd.d[j], &bd.d[l]) { dd_fwd[1] += 1; }
            }}

            // (e4) D_l ≼ D_j for j < l (reversed)
            for j in 1..=m { for l in j+1..=m {
                if pz(&bd.d[l]) || pz(&bd.d[j]) { continue; }
                dd_rev[0] += 1;
                if !exact_interlaces(&bd.d[l], &bd.d[j]) { dd_rev[1] += 1; }
            }}

            // (e6) U_j ≼ U_l for j > l (REVERSE UU)
            for j in 1..=m { for l in 1..j {
                if pz(&bd.u[j]) || pz(&bd.u[l]) { continue; }
                uu_rev[0] += 1;
                if !exact_interlaces(&bd.u[j], &bd.u[l]) { uu_rev[1] += 1; }
            }}

            // (f1) A_j^+ ≼ W_l^+ for ALL j, l
            for mp in 1..=m {
                for j in 1..=mp { for l in 1..=mp {
                    let a_plus = pa(&s[j], &t[j]);
                    let w_plus = pa(&pmt(&s[l]), &t[l]);
                    if pz(&a_plus) || pz(&w_plus) { continue; }
                    aw_plus[0] += 1;
                    if !exact_interlaces(&a_plus, &w_plus) { aw_plus[1] += 1; }
                }}
            }

            // (f2) S_j ≼ T_l for ALL j, l
            for j in 1..=m { for l in 1..=m {
                if pz(&s[j]) || pz(&t[l]) { continue; }
                st_all[0] += 1;
                if !exact_interlaces(&s[j], &t[l]) { st_all[1] += 1; }
            }}
        }
        println!("n={} done", n);
    }

    println!("\n=== RESULTS ===\n");
    println!("At level λ:");
    let show = |name: &str, c: [u64; 2]| {
        if c[0] == 0 { println!("  {}: (no tests)", name); }
        else { println!("  {}: {}/{} pass {}", name, c[0]-c[1], c[0], if c[1]==0 {"✓"} else {"✗ FAIL"}); }
    };
    show("A_j ≼ W_l (ALL j,l)", aw_all);
    show("D_j ≼ W_l (ALL j,l)", dw_all);
    show("D_j ≼ D_l (j≤l, fwd)", dd_fwd);
    show("D_l ≼ D_j (j<l, rev)", dd_rev);
    show("U_j ≼ U_l (j>l, rev)", uu_rev);

    println!("\nAt level λ+:");
    show("A_j^+ ≼ W_l^+ (ALL j,l)", aw_plus);
    show("S_j ≼ T_l (ALL j,l)", st_all);

    if aw_all[1]==0 && st_all[1]==0 && aw_plus[1]==0 {
        println!("\n>>> Key finding: AW universal, DU universal, and AW+ universal ALL hold.");
        println!(">>> Strengthening AW to all j,l appears viable as induction hypothesis.");
    }
}
