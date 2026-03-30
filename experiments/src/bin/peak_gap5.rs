//! Test whether strengthened AW at λ+ for j > l can be decomposed.
//! Specifically: does T_j ≪ W_l^+ hold for j > l? And S_l ≪ S_j for l < j?
//! These are needed to verify A_j^+ = S_j + T_j ≪ tS_l + T_l = W_l^+ via cone.

use num::{BigInt, BigRational};
use num_rational::Ratio;
use polynomial_tools::polynomial::{Polynomial, FieldRing};
use std::collections::BTreeSet;

type Poly = Polynomial<Ratio<BigInt>>;
type BR = Ratio<BigInt>;
fn br(n: i64) -> BR { BR::from_integer(BigInt::from(n)) }
fn to_poly(c: &[i64]) -> Poly {
    Polynomial::new(c.iter().map(|&x| BR::from_integer(BigInt::from(x))).collect())
}

fn poly_long_div(f: &Poly, g: &Poly) -> Poly { let (q, _) = f.div_rem(g); q }

fn polynomial_gcd(f: &Poly, g: &Poly) -> Poly { f.gcd(g) }

fn lagrange(pts: &[BR], vals: &[BR]) -> Vec<BR> {
    let p = Polynomial::lagrange_interpolation(pts, vals);
    let d = p.degree().unwrap_or(0);
    (0..=d).map(|i| p.coeff(i)).collect()
}

fn exact_interlaces(fc: &[i64], gc: &[i64]) -> bool {
    if let Some(result) = polynomial_tools::check_weak_interlacing(fc, gc) {
        return result;
    }
    polynomial_tools::check_interlacing_sturm(fc, gc).unwrap_or(false)
}

fn bruhat_lower_ideal(perm: &[u8]) -> Vec<Vec<u8>> {
    let n=perm.len(); let mut vis:BTreeSet<Vec<u8>>=BTreeSet::new(); let mut q:BTreeSet<Vec<u8>>=BTreeSet::new();
    q.insert(perm.to_vec());
    while let Some(cur)=q.pop_last(){for i in 0..n{for j in i+1..n{if cur[i]>cur[j]{
        let mut c=cur.clone();c.swap(i,j);if !vis.contains(&c){q.insert(c);}}}} vis.insert(cur);}
    vis.into_iter().collect()
}
fn board_to_perm(b:&[u8])->Vec<u8>{let n=b.len();let mut p=vec![0u8;n];let mut u=vec![false;n+1];
    for i in 0..n{for c in(1..=(b[i]as usize).min(n)).rev(){if !u[c]{p[i]=c as u8;u[c]=true;break;}}}p}
fn is_312_avoiding(p:&[u8])->bool{let n=p.len();for i in 0..n{for j in i+1..n{for k in j+1..n{
    if p[k]<p[i]&&p[i]<p[j]{return false;}}}}true}
fn peaks(w:&[u8])->usize{if w.len()<3{return 0;}(1..w.len()-1).filter(|&i|w[i-1]<w[i]&&w[i]>w[i+1]).count()}
fn pt(p:&[i64])->Vec<i64>{let mut v=p.to_vec();while v.len()>1&&*v.last().unwrap()==0{v.pop();}v}
fn pz(p:&[i64])->bool{p.iter().all(|&c|c==0)}
fn pa(a:&[i64],b:&[i64])->Vec<i64>{let l=a.len().max(b.len());let mut r=vec![0i64;l];
    for(i,&v)in a.iter().enumerate(){r[i]+=v;}for(i,&v)in b.iter().enumerate(){r[i]+=v;}pt(&r)}
fn pmt(p:&[i64])->Vec<i64>{let mut r=vec![0i64;p.len()+1];for(i,&v)in p.iter().enumerate(){r[i+1]=v;}pt(&r)}
fn gen_boards(n:usize)->Vec<Vec<u8>>{let mut r=vec![];let mut c=vec![];gb(n,n,0,&mut c,&mut r);r}
fn gb(n:usize,mx:usize,d:usize,c:&mut Vec<u8>,r:&mut Vec<Vec<u8>>){
    if d==n{r.push(c.clone());return;}for v in(d+1).max(if d>0{c[d-1]as usize}else{1})..=mx{c.push(v as u8);gb(n,mx,d+1,c,r);c.pop();}}

fn main() {
    let max_n:usize = std::env::args().nth(1).and_then(|s|s.parse().ok()).unwrap_or(7);
    println!("=== Decomposition tests for AW at λ+ (j>l), n ≤ {} ===\n", max_n);

    let mut ss_fwd=[0u64;2]; // S_l ≪ S_j for l < j
    let mut tj_wl=[0u64;2]; // T_j ≪ W_l^+ for j > l
    let mut sj_wl=[0u64;2]; // S_j ≪ W_l^+ for j > l
    let mut sl_ai=[0u64;2]; // S_l ≪ A_i for i ≥ l (needed for S_l ≪ S_j)

    for n in 1..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let perm = board_to_perm(board);
            if !is_312_avoiding(&perm) { continue; }
            let m = (board[0] as usize).min(n);
            let ideal = bruhat_lower_ideal(&perm);
            let mut d=vec![vec![0i64];m+1]; let mut u=vec![vec![0i64];m+1];
            for pi in &ideal {
                if pi.len()<2{let k=pi[0]as usize;if k<=m{while u[k].len()<1{u[k].push(0);}u[k][0]+=1;}continue;}
                let k=pi[0]as usize;if k>m{continue;} let pk=peaks(pi);
                let poly=if pi[0]>pi[1]{&mut d[k]}else{&mut u[k]};while poly.len()<=pk{poly.push(0);}poly[pk]+=1;
            }
            let mut a=vec![vec![0i64];m+1]; let mut w=vec![vec![0i64];m+1];
            for k in 1..=m{a[k]=pa(&d[k],&u[k]);w[k]=pa(&pmt(&d[k]),&u[k]);}
            let mut s=vec![vec![0i64];m+2]; let mut t_a=vec![vec![0i64];m+2];
            for k in 2..=m+1{s[k]=pa(&s[k-1],&a[k-1]);}
            for k in (1..=m).rev(){t_a[k]=pa(&t_a[k+1],&w[k]);}

            for j in 1..=m { for l in 1..j {
                // S_l ≪ S_j (l < j)
                if !pz(&s[l]) && !pz(&s[j]) {
                    ss_fwd[0]+=1; if !exact_interlaces(&s[l],&s[j]) { ss_fwd[1]+=1; }
                }
                // T_j ≪ W_l^+ for j > l (testing with all valid m')
                let w_l_plus = pa(&pmt(&s[l]), &t_a[l]);
                if !pz(&t_a[j]) && !pz(&w_l_plus) {
                    tj_wl[0]+=1; if !exact_interlaces(&t_a[j],&w_l_plus) { tj_wl[1]+=1; }
                }
                // S_j ≪ W_l^+
                if !pz(&s[j]) && !pz(&w_l_plus) {
                    sj_wl[0]+=1; if !exact_interlaces(&s[j],&w_l_plus) { sj_wl[1]+=1; }
                }
            }}
            // S_l ≪ A_i for i ≥ l
            for l in 1..=m { for i in l..=m {
                if !pz(&s[l]) && !pz(&a[i]) {
                    sl_ai[0]+=1; if !exact_interlaces(&s[l],&a[i]) { sl_ai[1]+=1; }
                }
            }}
        }
        println!("n={} done", n);
    }
    println!("\n=== RESULTS ===\n");
    let show=|n:&str,c:[u64;2]|{if c[0]==0{println!("  {}: (none)",n);}
        else{println!("  {}: {}/{} {}",n,c[0]-c[1],c[0],if c[1]==0{"✓"}else{"✗"});}};
    show("S_l ≪ S_j (l<j)", ss_fwd);
    show("T_j ≪ W_l^+ (j>l)", tj_wl);
    show("S_j ≪ W_l^+ (j>l)", sj_wl);
    show("S_l ≪ A_i (i≥l)", sl_ai);
    println!();
    if ss_fwd[1]==0 && tj_wl[1]==0 {
        println!("Both S_l≪S_j and T_j≪W_l^+ hold!");
        println!("Proof route: A_j^+ = S_j + T_j, each ≪ W_l^+ by cone.");
        println!("  S_j ≪ W_l^+: from S_j ≪ T_l (DU fix) + S_j ≪ tS_l (from S_l ≪ S_j + Wagner)");
        println!("  T_j ≪ W_l^+: directly verified");
    }
}
