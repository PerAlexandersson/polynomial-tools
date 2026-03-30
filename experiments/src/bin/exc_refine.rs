//! Explore excedance polynomial refinements for Ferrers boards.
//!
//! Try different refinements and check interlacing:
//! 1. F_k = exc poly for perms with first entry k (natural from permanent expansion)
//! 2. Check if {F_k} interlaces in some order
//! 3. Check the recursion structure when peeling off first row
//!
//! The exc polynomial = perm(M) where M_{i,j} = t if j>i, 1 if j≤i, 0 if j>λ_i.

use std::collections::BTreeSet;

fn exact_interlaces(fc: &[i64], gc: &[i64]) -> bool {
    polynomial_tools::check_weak_interlacing(fc, gc).unwrap_or(false)
}
fn bruhat_lower_ideal(perm:&[u8])->Vec<Vec<u8>>{
    let n=perm.len(); let mut vis:BTreeSet<Vec<u8>>=BTreeSet::new();
    let mut q:BTreeSet<Vec<u8>>=BTreeSet::new(); q.insert(perm.to_vec());
    while let Some(cur)=q.pop_last(){
        for i in 0..n{for j in i+1..n{if cur[i]>cur[j]{
            let mut c=cur.clone();c.swap(i,j);if!vis.contains(&c){q.insert(c);}
        }}} vis.insert(cur);
    } vis.into_iter().collect()
}
fn board_to_perm(b:&[u8])->Vec<u8>{
    let n=b.len();let mut p=vec![0u8;n];let mut u=vec![false;n+1];
    for i in 0..n{for c in(1..=(b[i]as usize).min(n)).rev(){
        if!u[c]{p[i]=c as u8;u[c]=true;break;}
    }} p
}
fn excedances(w:&[u8])->usize{(0..w.len()).filter(|&i|w[i]as usize>i+1).count()}
fn pt(p:&[i64])->Vec<i64>{let mut v=p.to_vec();while v.len()>1&&*v.last().unwrap()==0{v.pop();}v}
fn pz(p:&[i64])->bool{p.iter().all(|&c|c==0)}
fn pa(a:&[i64],b:&[i64])->Vec<i64>{
    let l=a.len().max(b.len());let mut r=vec![0i64;l];
    for(i,&v)in a.iter().enumerate(){r[i]+=v;}
    for(i,&v)in b.iter().enumerate(){r[i]+=v;}pt(&r)
}
fn pf(p:&[i64])->String{
    let p=pt(p);if pz(&p){"0".into()}else{
    let mut t=vec![];for(i,&c)in p.iter().enumerate(){if c==0{continue;}match(c,i){
        (c,0)=>t.push(format!("{}",c)),(1,1)=>t.push("t".into()),
        (c,1)=>t.push(format!("{}t",c)),(1,e)=>t.push(format!("t^{}",e)),
        (c,e)=>t.push(format!("{}t^{}",c,e)),
    }}t.join(" + ")}
}
fn gen_boards(n:usize)->Vec<Vec<u8>>{let mut r=vec![];let mut c=vec![];
    gb(n,n,0,&mut c,&mut r);r}
fn gb(n:usize,mx:usize,d:usize,c:&mut Vec<u8>,r:&mut Vec<Vec<u8>>){
    if d==n{r.push(c.clone());return;}
    for v in(d+1).max(if d>0{c[d-1]as usize}else{1})..=mx{
        c.push(v as u8);gb(n,mx,d+1,c,r);c.pop();
    }
}

fn main() {
    let max_n:usize=std::env::args().nth(1).and_then(|s|s.parse().ok()).unwrap_or(7);
    println!("=== Excedance refinement: F_k = exc poly for σ_1=k ===\n");

    let mut total = 0;
    let mut fk_fwd = [0u64;2]; // F_j ≼ F_l for j < l
    let mut fk_rev = [0u64;2]; // F_l ≼ F_j for j < l (reversed)
    let mut fk_rr = [0u64;2];  // F_k real-rooted
    let mut sk_tk = [0u64;2];  // S_k ≼ T_k where S_k=Σ_{j<k}F_j, T_k=Σ_{j≥k}tF_j... or similar

    // Also try: G_k = F_k/t for k≥2 (remove the t from exc at pos 1), G_1 = F_1
    // Then P = F_1 + tΣ_{k≥2} G_k
    // Check if {G_k} interlaces
    let mut gk_fwd = [0u64;2];
    let mut gk_rev = [0u64;2];

    for n in 2..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let nn=board.len(); let m=(board[0]as usize).min(nn);
            let perm=board_to_perm(board); let ideal=bruhat_lower_ideal(&perm);
            total += 1;

            // Compute F_k (exc poly conditioned on first entry k)
            let mut fk = vec![vec![0i64]; m+1];
            for pi in &ideal {
                let k = pi[0] as usize;
                if k > m { continue; }
                let exc = excedances(pi);
                while fk[k].len() <= exc { fk[k].push(0); }
                fk[k][exc] += 1;
            }
            for k in 1..=m { fk[k] = pt(&fk[k]); }

            // Check F_k real-rooted
            for k in 1..=m {
                if pz(&fk[k]) { continue; }
                fk_rr[0] += 1;
                if !polynomial_tools::is_real_rooted(&fk[k]) {
                    fk_rr[1] += 1;
                }
            }

            // Check F_j ≼ F_l for j < l (forward)
            for j in 1..m { for l in j+1..=m {
                if pz(&fk[j]) || pz(&fk[l]) { continue; }
                fk_fwd[0]+=1;
                if !exact_interlaces(&fk[j],&fk[l]) { fk_fwd[1]+=1; }
            }}
            // Check F_l ≼ F_j for j < l (reversed)
            for j in 1..m { for l in j+1..=m {
                if pz(&fk[l]) || pz(&fk[j]) { continue; }
                fk_rev[0]+=1;
                if !exact_interlaces(&fk[l],&fk[j]) { fk_rev[1]+=1; }
            }}

            // Compute G_k = F_k / t for k≥2 (divide out the t from exc at pos 1)
            // F_k for k≥2 should have a factor of t (since σ_1=k>1 means exc at pos 1)
            let mut gk = vec![vec![0i64]; m+1];
            gk[1] = fk[1].clone();
            for k in 2..=m {
                // Divide by t: shift coefficients down
                if fk[k].len() > 1 && fk[k][0] == 0 {
                    gk[k] = fk[k][1..].to_vec();
                } else {
                    gk[k] = fk[k].clone(); // can't divide by t
                }
                gk[k] = pt(&gk[k]);
            }

            // Check G_j ≼ G_l forward and reversed
            for j in 1..m { for l in j+1..=m {
                if pz(&gk[j]) || pz(&gk[l]) { continue; }
                gk_fwd[0]+=1;
                if !exact_interlaces(&gk[j],&gk[l]) { gk_fwd[1]+=1; }
            }}
            for j in 1..m { for l in j+1..=m {
                if pz(&gk[l]) || pz(&gk[j]) { continue; }
                gk_rev[0]+=1;
                if !exact_interlaces(&gk[l],&gk[j]) { gk_rev[1]+=1; }
            }}
        }
        println!("n={}: boards={}", n, total);
    }

    println!("\n=== RESULTS ===");
    println!("F_k real-rooted: {}/{} {}", fk_rr[0]-fk_rr[1], fk_rr[0], if fk_rr[1]==0{"✓"}else{"FAIL"});
    println!("F_j≼F_l (j<l, fwd): {}/{} {}", fk_fwd[0]-fk_fwd[1], fk_fwd[0], if fk_fwd[1]==0{"✓"}else{"FAIL"});
    println!("F_l≼F_j (j<l, rev): {}/{} {}", fk_rev[0]-fk_rev[1], fk_rev[0], if fk_rev[1]==0{"✓"}else{"FAIL"});
    println!("G_j≼G_l (j<l, fwd): {}/{} {}", gk_fwd[0]-gk_fwd[1], gk_fwd[0], if gk_fwd[1]==0{"✓"}else{"FAIL"});
    println!("G_l≼G_j (j<l, rev): {}/{} {}", gk_rev[0]-gk_rev[1], gk_rev[0], if gk_rev[1]==0{"✓"}else{"FAIL"});

    // Print example
    println!("\n=== S_5 detail ===");
    let board=vec![5u8,5,5,5,5]; let m=5;
    let perm=board_to_perm(&board); let ideal=bruhat_lower_ideal(&perm);
    let mut fk=vec![vec![0i64];m+1];
    for pi in &ideal { let k=pi[0]as usize; if k>m{continue;} let exc=excedances(pi);
        while fk[k].len()<=exc{fk[k].push(0);} fk[k][exc]+=1; }
    for k in 1..=m { fk[k]=pt(&fk[k]); println!("F_{} = {}", k, pf(&fk[k])); }
}
