//! Refined stability test: only test (s,t) with Im(s) >= 0.5, Im(t) >= 0.5.
//! This avoids boundary effects.
#[derive(Clone, Copy)] struct C64 { re: f64, im: f64 }
impl C64 { fn new(re: f64, im: f64) -> Self { Self{re,im} }
    fn mul(self, o: Self) -> Self { C64::new(self.re*o.re-self.im*o.im, self.re*o.im+self.im*o.re) }
    fn add(self, o: Self) -> Self { C64::new(self.re+o.re, self.im+o.im) }
    fn norm2(self) -> f64 { self.re*self.re + self.im*self.im } }
fn perm_c(mat: &[Vec<C64>]) -> C64 { let n = mat.len(); if n == 0 { return C64::new(1.0,0.0); }
    let mut r = C64::new(0.0,0.0); let mut ix: Vec<usize> = (0..n).collect();
    loop { let mut p = C64::new(1.0,0.0); for i in 0..n { p = p.mul(mat[i][ix[i]]); } r = r.add(p); if !np(&mut ix) { break; } } r }
fn np(a: &mut [usize]) -> bool { let n=a.len(); if n<=1 {return false;} let mut i=n-2;
    while a[i]>=a[i+1] { if i==0 {return false;} i-=1; } let mut j=n-1;
    while a[j]<=a[i] { j-=1; } a.swap(i,j); a[i+1..].reverse(); true }
fn board_to_perm(b: &[u8]) -> Vec<u8> { let n=b.len(); let mut p=vec![0u8;n]; let mut u=vec![false;n+1];
    for i in 0..n { for c in (1..=(b[i] as usize).min(n)).rev() { if !u[c] { p[i]=c as u8; u[c]=true; break; } } } p }
fn is_312(perm: &[u8]) -> bool { let n=perm.len(); for i in 0..n { for j in i+1..n { for k in j+1..n {
    if perm[k]<perm[i] && perm[i]<perm[j] { return false; } } } } true }
fn gen_boards(n: usize) -> Vec<Vec<u8>> { let mut r=vec![]; let mut c=vec![]; gb(n,n,0,&mut c,&mut r); r }
fn gb(n: usize, mx: usize, d: usize, c: &mut Vec<u8>, r: &mut Vec<Vec<u8>>) {
    if d==n { r.push(c.clone()); return; }
    for v in (d+1).max(if d>0 {c[d-1] as usize} else {1})..=mx { c.push(v as u8); gb(n,mx,d+1,c,r); c.pop(); } }
fn sub_partitions(lam: &[u8]) -> Vec<Vec<u8>> { let n=lam.len(); let mut res=Vec::new(); let mut mu=vec![0u8;n];
    fn gen(l: &[u8], m: &mut Vec<u8>, p: usize, mx: u8, r: &mut Vec<Vec<u8>>) { if p==l.len() { r.push(m.clone()); return; }
        let u=l[p].min(mx); for v in 0..=u { m[p]=v; gen(l,m,p+1,v,r); } }
    gen(lam,&mut mu,0,lam[0],&mut res); res }
fn main() {
    // Test points with Im >= 0.5
    let pts: Vec<C64> = vec![
        C64::new(0.0, 0.5), C64::new(0.0, 1.0), C64::new(0.0, 2.0), C64::new(0.0, 5.0),
        C64::new(1.0, 0.5), C64::new(1.0, 1.0), C64::new(1.0, 3.0),
        C64::new(-1.0, 0.5), C64::new(-1.0, 1.0), C64::new(-1.0, 2.0),
        C64::new(-2.0, 0.5), C64::new(-2.0, 1.0), C64::new(-3.0, 0.5),
        C64::new(-5.0, 1.0), C64::new(-10.0, 1.0),
        C64::new(0.5, 0.5), C64::new(2.0, 1.0), C64::new(-0.5, 0.5),
        C64::new(3.0, 0.5), C64::new(-4.0, 2.0),
    ];
    let mut total = 0u64;
    let mut min_norm = f64::MAX;
    let mut worst = String::new();
    for n in 2..=6usize {
        for board in &gen_boards(n) {
            let perm = board_to_perm(board); if !is_312(&perm) { continue; }
            let m = *board.iter().max().unwrap() as usize;
            for mu in &sub_partitions(board) {
                for j in 0..n {
                    if mu[j] >= board[j] { continue; }
                    let mut mp = mu.clone(); mp[j] += 1;
                    if j > 0 && mp[j] > mp[j-1] { continue; }
                    let cj = mu[j] as usize + 1;
                    for &sv in &pts { for &tv in &pts {
                        let mut mat = Vec::new();
                        for i in 0..n { let mut row = Vec::new();
                            for k in 1..=m {
                                if k > board[i] as usize { row.push(C64::new(0.0,0.0)); }
                                else if i == j && k == cj { row.push(sv); }
                                else if k <= mp[i] as usize { row.push(C64::new(1.0,0.0)); }
                                else { row.push(tv); }
                            } mat.push(row); }
                        let p = perm_c(&mat);
                        total += 1;
                        if p.norm2() < min_norm { min_norm = p.norm2();
                            worst = format!("b={:?} m'={:?} j={} s={:.1}+{:.1}i t={:.1}+{:.1}i |F|²={:.2e}",
                                board, mp, j, sv.re, sv.im, tv.re, tv.im, p.norm2()); }
                    }}
                }
            }
        }
    }
    println!("=== Stability (Im >= 0.5), n <= 6 ===");
    println!("  {} evaluations", total);
    println!("  Min |F|² = {:.2e}", min_norm);
    println!("  Worst: {}", worst);
    if min_norm > 1e-6 { println!("  STABLE: no near-zeros in upper half-plane <<<"); }
    else { println!("  UNSTABLE: near-zero detected!"); }
}
