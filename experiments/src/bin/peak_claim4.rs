//! Deep exploration of Claim 4: W_k^+ ≼ W_{k+1}^+ = W_k^+ + (t-1)U_k.
//!
//! For each board and each step k, prints:
//!   - W_k^+, U_k, W_{k+1}^+
//!   - Q_k = W_k^+ - U_k, tU_k (the interlacing pair)
//!   - Whether U_k divides W_k^+ or shares factors
//!   - The relationship between U_k and the partial/tail sums
//!
//! Also tests if D_k ≼ U_k (which would give W_k = tD_k+U_k ≼ tU_k trivially).

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
fn ps(a:&[i64],b:&[i64])->Vec<i64>{
    let l=a.len().max(b.len()); let mut r=vec![0i64;l];
    for(i,&v)in a.iter().enumerate(){r[i]+=v;}
    for(i,&v)in b.iter().enumerate(){r[i]-=v;} pt(&r)
}
fn pmt(p:&[i64])->Vec<i64>{
    let mut r=vec![0i64;p.len()+1];
    for(i,&v)in p.iter().enumerate(){r[i+1]=v;} pt(&r)
}
fn psc(p:&[i64],c:i64)->Vec<i64>{ pt(&p.iter().map(|&x|x*c).collect::<Vec<_>>()) }
fn pf(p:&[i64])->String{
    let p=pt(p); if pz(&p){"0".into()} else {
    let mut t=vec![]; for(i,&c)in p.iter().enumerate(){ if c==0{continue;} match(c,i){
        (c,0)=>t.push(format!("{}",c)),(1,1)=>t.push("t".into()),
        (c,1)=>t.push(format!("{}t",c)),(1,e)=>t.push(format!("t^{}",e)),
        (c,e)=>t.push(format!("{}t^{}",c,e)),
    }} t.join(" + ")}
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
    let max_n: usize = std::env::args().nth(1).and_then(|s|s.parse().ok()).unwrap_or(6);

    println!("=== Claim 4 exploration: W_k^+ ≼ W_{{k+1}}^+ = W_k^+ + (t-1)U_k ===\n");

    // Test 1: Is D_k ≼ U_k always? (This is known from earlier.)
    // Test 2: Is U_k ≼ T_k? (U_k ≼ tail sum of W's starting at k)
    // Test 3: Is U_k ≼ S_{k+1}? (U_k ≼ partial sum of A's up to k)
    // Test 4: Is V_k^tail ≼ tU_k? where V_k^tail = Σ_{j>k} U_j
    // Test 5: Is (P_D + V_k) ≼ U_k? where P_D = Σ D_j, V_k = Σ_{j<k} U_j
    //         (Since W_k^+-U_k = t(P_D+V_k) + V_k^tail, and condition is W_k^+-U_k ≼ tU_k)
    // Test 6: Is S_k ≼ U_k? (S_k = Σ_{j<k} A_j ≼ U_k?)
    // Test 7: consecutive W_k^+ - U_k = tS_k + T_k - U_k = tS_k + tD_k + T_{k+1}
    //         = t(S_k + D_k) + T_{k+1}. So condition is t(S_k+D_k)+T_{k+1} ≼ tU_k.
    //         Equivalently after GCD: (S_k+D_k) and T_{k+1} relate to U_k somehow.
    //         Note S_k + D_k = Σ_{j<k} A_j + D_k = Σ_{j<k}(D_j+U_j) + D_k = Σ_{j≤k}D_j + V_k.
    // Test 8: Direct: is (S_k + D_k) ≼ U_k?
    //         And: is T_{k+1} ≼ tU_k?

    let mut tests = [0u64; 10]; // counts
    let mut fails = [0u64; 10]; // fails

    for n in 2..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let nn = board.len();
            let m = (board[0] as usize).min(nn);
            let perm = board_to_perm(board);
            let ideal = bruhat_lower_ideal(&perm);
            let mut d = vec![vec![0i64]; m+1];
            let mut u = vec![vec![0i64]; m+1];
            for pi in &ideal {
                if pi.len()<2 { let k=pi[0] as usize; if k<=m { while u[k].len()<1{u[k].push(0);} u[k][0]+=1; } continue; }
                let k=pi[0] as usize; if k>m{continue;}
                let pk=peaks(pi);
                let poly=if pi[0]>pi[1]{&mut d[k]}else{&mut u[k]};
                while poly.len()<=pk{poly.push(0);} poly[pk]+=1;
            }
            let mut a=vec![vec![0i64];m+1]; let mut w=vec![vec![0i64];m+1];
            for k in 1..=m { a[k]=pa(&d[k],&u[k]); w[k]=pa(&pmt(&d[k]),&u[k]); }

            let mut s=vec![vec![0i64];m+2]; let mut t_arr=vec![vec![0i64];m+2];
            for k in 2..=m+1 { s[k]=pa(&s[k-1],&a[k-1]); }
            for k in (1..=m).rev() { t_arr[k]=pa(&t_arr[k+1],&w[k]); }

            for k in 1..=m {
                if pz(&u[k]) { continue; }

                // Test 1: D_k ≼ U_k
                if !pz(&d[k]) { tests[1]+=1; if !exact_interlaces(&d[k],&u[k]) { fails[1]+=1; } }
                // Test 2: U_k ≼ T_k
                if !pz(&t_arr[k]) { tests[2]+=1; if !exact_interlaces(&u[k],&t_arr[k]) { fails[2]+=1; } }
                // Test 3: U_k ≼ S_{k+1} (if exists)
                if k<m && !pz(&s[k+1]) { tests[3]+=1; if !exact_interlaces(&u[k],&s[k+1]) { fails[3]+=1; } }
                // Test 6: S_k ≼ U_k
                if !pz(&s[k]) { tests[6]+=1; if !exact_interlaces(&s[k],&u[k]) { fails[6]+=1; } }
                // Test 8: (S_k + D_k) ≼ U_k
                let skdk = pa(&s[k], &d[k]);
                if !pz(&skdk) { tests[8]+=1; if !exact_interlaces(&skdk,&u[k]) { fails[8]+=1; } }

                // Test T_{k+1} ≼ tU_k
                let tuk = pmt(&u[k]);
                if k<m && !pz(&t_arr[k+1]) { tests[9]+=1; if !exact_interlaces(&t_arr[k+1],&tuk) { fails[9]+=1; } }
            }
        }
    }

    println!("Results (across all boards n ≤ {}):", max_n);
    let labels = ["", "D_k≼U_k", "U_k≼T_k", "U_k≼S_{k+1}", "", "",
                   "S_k≼U_k", "", "(S_k+D_k)≼U_k", "T_{k+1}≼tU_k"];
    for i in [1,2,3,6,8,9] {
        if tests[i] > 0 {
            println!("  {}: {}/{} pass ({})", labels[i], tests[i]-fails[i], tests[i],
                if fails[i]==0 {"ALL PASS"} else {"FAILURES"});
        }
    }

    // === Additional tests for the complete proof ===
    println!("\n--- Additional conditions ---\n");
    let mut du_all = [0u64; 2]; // D_j ≼ U_l for ALL j,l
    let mut uu_fwd = [0u64; 2]; // U_j ≼ U_l for j ≤ l (forward)
    let mut dd_rev = [0u64; 2]; // D_l ≼ D_j for j ≤ l (reversed: bigger index more LEFT)
    let mut st_all = [0u64; 2]; // S_j ≼ T_l for ALL j,l

    for n in 2..=max_n {
        let boards = gen_boards(n);
        for board in &boards {
            let nn = board.len();
            let m = (board[0] as usize).min(nn);
            let perm = board_to_perm(board);
            let ideal = bruhat_lower_ideal(&perm);
            let mut d = vec![vec![0i64]; m+1];
            let mut u_p = vec![vec![0i64]; m+1];
            for pi in &ideal {
                if pi.len()<2 { let k=pi[0] as usize; if k<=m { while u_p[k].len()<1{u_p[k].push(0);} u_p[k][0]+=1; } continue; }
                let k=pi[0] as usize; if k>m{continue;}
                let pk=peaks(pi);
                let poly=if pi[0]>pi[1]{&mut d[k]}else{&mut u_p[k]};
                while poly.len()<=pk{poly.push(0);} poly[pk]+=1;
            }
            let mut a=vec![vec![0i64];m+1]; let mut w=vec![vec![0i64];m+1];
            for k in 1..=m { a[k]=pa(&d[k],&u_p[k]); w[k]=pa(&pmt(&d[k]),&u_p[k]); }
            let mut s=vec![vec![0i64];m+2]; let mut t_a=vec![vec![0i64];m+2];
            for k in 2..=m+1 { s[k]=pa(&s[k-1],&a[k-1]); }
            for k in (1..=m).rev() { t_a[k]=pa(&t_a[k+1],&w[k]); }

            // D_j ≼ U_l for ALL j,l (both nonzero)
            for j in 1..=m { for l in 1..=m {
                if pz(&d[j]) || pz(&u_p[l]) { continue; }
                du_all[0]+=1;
                if !exact_interlaces(&d[j],&u_p[l]) { du_all[1]+=1; }
            }}
            // U_j ≼ U_l for j ≤ l
            for j in 1..=m { for l in j+1..=m {
                if pz(&u_p[j]) || pz(&u_p[l]) { continue; }
                uu_fwd[0]+=1;
                if !exact_interlaces(&u_p[j],&u_p[l]) { uu_fwd[1]+=1; }
            }}
            // D_l ≼ D_j for j < l (reversed)
            for j in 1..=m { for l in j+1..=m {
                if pz(&d[l]) || pz(&d[j]) { continue; }
                dd_rev[0]+=1;
                if !exact_interlaces(&d[l],&d[j]) { dd_rev[1]+=1; }
            }}
            // S_j ≼ T_l for ALL j,l
            for j in 1..=m { for l in 1..=m {
                if pz(&s[j]) || pz(&t_a[l]) { continue; }
                st_all[0]+=1;
                if !exact_interlaces(&s[j],&t_a[l]) { st_all[1]+=1; }
            }}
        }
    }
    println!("D_j≼U_l (ALL j,l): {}/{} pass {}", du_all[0]-du_all[1], du_all[0], if du_all[1]==0{"✓"}else{"FAIL"});
    println!("U_j≼U_l (j≤l, fwd): {}/{} pass {}", uu_fwd[0]-uu_fwd[1], uu_fwd[0], if uu_fwd[1]==0{"✓"}else{"FAIL"});
    println!("D_l≼D_j (j<l, rev): {}/{} pass {}", dd_rev[0]-dd_rev[1], dd_rev[0], if dd_rev[1]==0{"✓"}else{"FAIL"});
    println!("S_j≼T_l (ALL j,l):  {}/{} pass {}", st_all[0]-st_all[1], st_all[0], if st_all[1]==0{"✓"}else{"FAIL"});

    println!("\n=== Detailed output for S_5 ===\n");
    let board = vec![5u8,5,5,5,5];
    let nn = 5usize; let m = 5;
    let perm = board_to_perm(&board);
    let ideal = bruhat_lower_ideal(&perm);
    let mut d = vec![vec![0i64]; m+1];
    let mut u_arr = vec![vec![0i64]; m+1];
    for pi in &ideal {
        if pi.len()<2 { continue; }
        let k=pi[0] as usize; if k>m{continue;}
        let pk=peaks(pi);
        let poly=if pi[0]>pi[1]{&mut d[k]}else{&mut u_arr[k]};
        while poly.len()<=pk{poly.push(0);} poly[pk]+=1;
    }
    let mut a=vec![vec![0i64];m+1]; let mut w=vec![vec![0i64];m+1];
    for k in 1..=m { a[k]=pa(&d[k],&u_arr[k]); w[k]=pa(&pmt(&d[k]),&u_arr[k]); }
    let mut s=vec![vec![0i64];m+2]; let mut t_arr=vec![vec![0i64];m+2];
    for k in 2..=m+1 { s[k]=pa(&s[k-1],&a[k-1]); }
    for k in (1..=m).rev() { t_arr[k]=pa(&t_arr[k+1],&w[k]); }

    for k in 1..=m {
        println!("k={}:", k);
        println!("  D_k = {}", pf(&d[k]));
        println!("  U_k = {}", pf(&u_arr[k]));
        println!("  A_k = {}", pf(&a[k]));
        println!("  W_k = {}", pf(&w[k]));
        println!("  S_k = {}", pf(&s[k]));
        println!("  T_k = {}", pf(&t_arr[k]));
        let w_plus_k = pa(&pmt(&s[k]), &t_arr[k]);
        let q_k = ps(&w_plus_k, &u_arr[k]); // W_k^+ - U_k
        let tu_k = pmt(&u_arr[k]);
        println!("  W_k^+ = {}", pf(&w_plus_k));
        println!("  W_k^+ - U_k = {}", pf(&q_k));
        println!("  tU_k = {}", pf(&tu_k));
        let skdk = pa(&s[k], &d[k]);
        println!("  S_k+D_k = {}", pf(&skdk));
        if k <= m-1 {
            println!("  T_{{k+1}} = {}", pf(&t_arr[k+1]));
        }
        println!("  (S_k+D_k)≼U_k: {}", if pz(&skdk)||pz(&u_arr[k]){"trivial"} else if exact_interlaces(&skdk,&u_arr[k]){"YES"} else {"NO"});
        println!("  (W_k^+-U_k)≼tU_k: {}", if pz(&q_k)||pz(&tu_k){"trivial"} else if exact_interlaces(&q_k,&tu_k){"YES"} else {"NO"});
        println!();
    }
}
