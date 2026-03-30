//! Exhaustive test of the INDUCTION STEP for peak polynomial real-rootedness.
//!
//! For every valid (λ, λ⁺) pair, verify:
//!   IH(λ) = {A_j ≼ W_l, W_j ≼ W_l for j ≤ l} at level λ
//!   ⟹ IH(λ⁺) at level λ⁺.
//!
//! This directly tests the gap: W_j⁺ ≼ W_l⁺ for j < l.

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

struct BD { m:usize, a:Vec<Vec<i64>>, w:Vec<Vec<i64>> }

fn compute_board(board:&[u8])->BD {
    let n=board.len(); let m=(board[0] as usize).min(n);
    let perm=board_to_perm(board); let ideal=bruhat_lower_ideal(&perm);
    let mut d=vec![vec![0i64];m+1]; let mut u=vec![vec![0i64];m+1];
    for pi in &ideal {
        if pi.len()<2 { let k=pi[0] as usize; if k<=m { while u[k].len()<1{u[k].push(0);} u[k][0]+=1; } continue; }
        let k=pi[0] as usize; if k>m { continue; }
        let pk=peaks(pi);
        let poly=if pi[0]>pi[1]{&mut d[k]}else{&mut u[k]};
        while poly.len()<=pk{poly.push(0);} poly[pk]+=1;
    }
    let mut a=vec![vec![0i64];m+1]; let mut w=vec![vec![0i64];m+1];
    for k in 1..=m { a[k]=pa(&d[k],&u[k]); w[k]=pa(&pmt(&d[k]),&u[k]); }
    BD{m,a,w}
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

fn main() {
    let max_n:usize = std::env::args().nth(1).and_then(|s|s.parse().ok()).unwrap_or(7);
    println!("INDUCTION STEP test: IH(λ) ⟹ IH(λ⁺) for all (λ,λ⁺) pairs, n ≤ {}\n", max_n);

    let mut total_pairs = 0u64;
    let mut total_aw = 0u64;  // A_j^+ ≼ W_l^+ tests
    let mut total_ww = 0u64;  // W_j^+ ≼ W_l^+ tests (THE GAP)
    let mut aw_fails = 0u64;
    let mut ww_fails = 0u64;
    let mut ih_fails = 0u64;  // IH failing at some level

    for n in 2..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let bd = compute_board(board);
            let m = bd.m;

            // Check IH holds at this level
            let mut ih_ok = true;
            for j in 1..=m { for l in j..=m {
                if !pz(&bd.a[j]) && !pz(&bd.w[l]) && !exact_interlaces(&bd.a[j],&bd.w[l]) { ih_ok=false; }
                if j<l && !pz(&bd.w[j]) && !pz(&bd.w[l]) && !exact_interlaces(&bd.w[j],&bd.w[l]) { ih_ok=false; }
            }}
            if !ih_ok { ih_fails+=1; continue; }

            // Build S_k, T_k from this level
            let mut s = vec![vec![0i64]; m+2];
            let mut t = vec![vec![0i64]; m+2];
            for k in 2..=m+1 { s[k] = pa(&s[k-1], &bd.a[k-1]); }
            for k in (1..=m).rev() { t[k] = pa(&t[k+1], &bd.w[k]); }

            // For each valid m' (width of new first row), test IH at λ⁺
            for mp in 1..=m {
                total_pairs += 1;

                // Compute A_k^+ = S_k + T_k, W_k^+ = tS_k + T_k for k=1..mp
                let mut a_plus: Vec<Vec<i64>> = vec![vec![0i64]; mp+1];
                let mut w_plus: Vec<Vec<i64>> = vec![vec![0i64]; mp+1];
                for k in 1..=mp {
                    a_plus[k] = pa(&s[k], &t[k]);
                    w_plus[k] = pa(&pmt(&s[k]), &t[k]);
                }

                // Test A_j^+ ≼ W_l^+ for j ≤ l
                for j in 1..=mp { for l in j..=mp {
                    if pz(&a_plus[j]) || pz(&w_plus[l]) { continue; }
                    total_aw += 1;
                    if !exact_interlaces(&a_plus[j], &w_plus[l]) {
                        aw_fails += 1;
                        if aw_fails <= 3 {
                            println!("FAIL A_{}^+≼W_{}^+ for λ={:?} m'={}", j, l, board, mp);
                        }
                    }
                }}

                // Test W_j^+ ≼ W_l^+ for j < l  (THE GAP)
                for j in 1..mp { for l in j+1..=mp {
                    if pz(&w_plus[j]) || pz(&w_plus[l]) { continue; }
                    total_ww += 1;
                    if !exact_interlaces(&w_plus[j], &w_plus[l]) {
                        ww_fails += 1;
                        if ww_fails <= 3 {
                            println!("FAIL W_{}^+≼W_{}^+ for λ={:?} m'={}", j, l, board, mp);
                        }
                    }
                }}
            }
        }
        println!("n={}: pairs={}, A≼W tests={} (fail={}), W≼W tests={} (fail={}), IH_fail={}",
            n, total_pairs, total_aw, aw_fails, total_ww, ww_fails, ih_fails);
    }

    println!("\n=== SUMMARY ===");
    println!("Total (λ,λ⁺) pairs: {}", total_pairs);
    println!("A_j^+≼W_l^+ tests: {} (failures: {})", total_aw, aw_fails);
    println!("W_j^+≼W_l^+ tests: {} (failures: {})  ← THE GAP", total_ww, ww_fails);
    println!("IH failures at base level: {}", ih_fails);
    if aw_fails==0 && ww_fails==0 && ih_fails==0 {
        println!("\nALL INDUCTION STEPS VERIFIED.");
    }
}
