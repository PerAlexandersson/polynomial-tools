use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::env;

#[derive(Debug)]
struct Args {
    max_ordinary_moves: usize,
    max_loop_len: usize,
    bad_grid: usize,
    max_extended_len: usize,
    residual_report: bool,
    picture_report: bool,
    stacked_report: bool,
    extended_report: bool,
    verbose: bool,
}

impl Args {
    fn parse() -> Self {
        let mut args = Self {
            max_ordinary_moves: 8,
            max_loop_len: 8,
            bad_grid: 3,
            max_extended_len: 12,
            residual_report: false,
            picture_report: false,
            stacked_report: false,
            extended_report: false,
            verbose: false,
        };

        let mut iter = env::args().skip(1);
        while let Some(arg) = iter.next() {
            match arg.as_str() {
                "--max-ordinary-moves" => {
                    args.max_ordinary_moves = iter
                        .next()
                        .expect("--max-ordinary-moves needs a value")
                        .parse()
                        .expect("invalid --max-ordinary-moves value");
                }
                "--max-loop-len" => {
                    args.max_loop_len = iter
                        .next()
                        .expect("--max-loop-len needs a value")
                        .parse()
                        .expect("invalid --max-loop-len value");
                }
                "--bad-grid" => {
                    args.bad_grid = iter
                        .next()
                        .expect("--bad-grid needs a value")
                        .parse()
                        .expect("invalid --bad-grid value");
                }
                "--max-extended-len" => {
                    args.max_extended_len = iter
                        .next()
                        .expect("--max-extended-len needs a value")
                        .parse()
                        .expect("invalid --max-extended-len value");
                }
                "--residual-report" => args.residual_report = true,
                "--picture-report" => args.picture_report = true,
                "--stacked-report" => args.stacked_report = true,
                "--extended-report" => args.extended_report = true,
                "--verbose" => args.verbose = true,
                "--help" | "-h" => {
                    print_help();
                    std::process::exit(0);
                }
                _ => panic!("unknown argument {arg}; use --help"),
            }
        }
        args
    }
}

fn print_help() {
    println!(
        "\
Brute-force checks for the cylindric MN inclusion-exclusion sum.

Options:
  --max-ordinary-moves N   largest ordinary ribbon move word length [default: 8]
  --max-loop-len N         largest cylinder circumference x+y [default: 8]
  --bad-grid N             scan non-ribbon subsets in an N x N ordinary grid [default: 3]
  --max-extended-len N     largest extended path word length for --extended-report [default: 12]
  --residual-report        print anchored residual checks for the post-peeling cases
  --picture-report         print diagnostics for shapes encoded from toggleProofPics.tex
  --stacked-report         print representative shifted-loop band diagnostics
  --extended-report        run exploratory extended-path diagnostics; not part of pass/fail
  --verbose                print every checked base case
"
    );
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Cell {
    col: i32,
    row: i32,
}

#[derive(Clone, Debug)]
struct Cylinder {
    x: i32,
    y: i32,
}

impl Cylinder {
    fn new(x: i32, y: i32) -> Self {
        assert!(x > 0 && y > 0);
        Self { x, y }
    }

    fn norm(&self, col: i32, row: i32) -> Cell {
        let q = row.div_euclid(self.y);
        Cell {
            col: col - q * self.x,
            row: row - q * self.y,
        }
    }

    fn right(&self, c: Cell) -> Cell {
        self.norm(c.col + 1, c.row)
    }

    fn up(&self, c: Cell) -> Cell {
        self.norm(c.col, c.row + 1)
    }

    fn down(&self, c: Cell) -> Cell {
        self.norm(c.col, c.row - 1)
    }
}

#[derive(Clone, Debug)]
struct Shape {
    cyl: Cylinder,
    cells: BTreeSet<Cell>,
}

#[derive(Clone, Debug)]
struct Poset {
    n: usize,
    cover_edges: Vec<(usize, usize)>,
    strict_edges: Vec<(usize, usize)>,
}

#[derive(Default)]
struct CheckStats {
    checked: usize,
    failed: usize,
}

fn shape_from_path(cyl: Cylinder, word: &[char], closed_loop: bool, shift: i32) -> Shape {
    let mut cells = BTreeSet::new();
    let mut cur = cyl.norm(shift, 0);
    let take = if closed_loop {
        word.len()
    } else {
        word.len() + 1
    };
    for step_index in 0..take {
        cells.insert(cur);
        if step_index == word.len() {
            break;
        }
        cur = match word[step_index] {
            'R' => cyl.right(cur),
            'U' => cyl.up(cur),
            other => panic!("unexpected step {other}"),
        };
    }
    Shape { cyl, cells }
}

fn shape_from_intervals(cyl: Cylinder, intervals: &[(i32, i32, i32)]) -> Shape {
    let mut cells = BTreeSet::new();
    for &(row, start_col, end_col) in intervals {
        for col in start_col..end_col {
            cells.insert(cyl.norm(col, row));
        }
    }
    Shape { cyl, cells }
}

fn union_shifted_loops(cyl: Cylinder, word: &[char], layers: i32) -> Shape {
    let mut cells = BTreeSet::new();
    for shift in 0..layers {
        let layer = shape_from_path(cyl.clone(), word, true, shift);
        cells.extend(layer.cells);
    }
    Shape { cyl, cells }
}

fn poset(shape: &Shape) -> Poset {
    let index: HashMap<Cell, usize> = shape
        .cells
        .iter()
        .enumerate()
        .map(|(i, c)| (*c, i))
        .collect();
    let mut cover_edges = BTreeSet::new();
    let mut strict_edges = BTreeSet::new();

    for &cell in &shape.cells {
        let u = index[&cell];

        let right = shape.cyl.right(cell);
        if let Some(&v) = index.get(&right) {
            cover_edges.insert((u, v));
        }

        let up = shape.cyl.up(cell);
        if let Some(&v_above) = index.get(&up) {
            // The poset relation is "above < below".
            cover_edges.insert((v_above, u));
            strict_edges.insert((v_above, u));
        }
    }

    Poset {
        n: shape.cells.len(),
        cover_edges: cover_edges.into_iter().collect(),
        strict_edges: strict_edges.into_iter().collect(),
    }
}

fn source_vertices(p: &Poset) -> Vec<usize> {
    let mut indeg = vec![0usize; p.n];
    for &(_, v) in &p.cover_edges {
        indeg[v] += 1;
    }
    indeg
        .into_iter()
        .enumerate()
        .filter_map(|(v, d)| (d == 0).then_some(v))
        .collect()
}

fn anchored_residual_poset(shape: &Shape, anchor_weight: i64) -> (Poset, Vec<i64>, usize) {
    let base = poset(shape);
    let sources = source_vertices(&base);
    let mut cover_edges = Vec::new();
    let mut strict_edges = Vec::new();

    for &(u, v) in &base.cover_edges {
        cover_edges.push((u + 1, v + 1));
    }
    for &(u, v) in &base.strict_edges {
        strict_edges.push((u + 1, v + 1));
    }
    for &v in &sources {
        cover_edges.push((0, v + 1));
        strict_edges.push((0, v + 1));
    }

    let mut weights = vec![1i64; base.n + 1];
    weights[0] = anchor_weight;
    (
        Poset {
            n: base.n + 1,
            cover_edges,
            strict_edges,
        },
        weights,
        sources.len(),
    )
}

fn anchored_residual_sum(shape: &Shape, anchor_weight: i64) -> (i64, usize) {
    let (p, weights, source_count) = anchored_residual_poset(shape, anchor_weight);
    if p.strict_edges.len() >= 63 {
        return (i64::MIN, source_count);
    }
    (alternating_sum_weighted(&p, &weights), source_count)
}

fn m_value(p: &Poset, mask: u64) -> usize {
    let weights = vec![1i64; p.n];
    m_value_weighted(p, mask, &weights) as usize
}

fn m_value_weighted(p: &Poset, mask: u64, weights: &[i64]) -> i64 {
    assert_eq!(p.n, weights.len());
    let mut graph = vec![Vec::<usize>::new(); p.n];
    for &(u, v) in &p.cover_edges {
        graph[u].push(v);
    }
    for (i, &(u, v)) in p.strict_edges.iter().enumerate() {
        if ((mask >> i) & 1) == 1 {
            graph[v].push(u);
        }
    }

    let comp = strongly_connected_components(&graph);
    let comp_count = comp.iter().copied().max().map_or(0, |m| m + 1);
    let mut size = vec![0i64; comp_count];
    let mut indeg = vec![0usize; comp_count];
    for (v, &c) in comp.iter().enumerate() {
        size[c] += weights[v];
    }
    for u in 0..p.n {
        for &v in &graph[u] {
            let cu = comp[u];
            let cv = comp[v];
            if cu != cv {
                indeg[cv] += 1;
            }
        }
    }
    let sources: Vec<_> = (0..comp_count).filter(|&c| indeg[c] == 0).collect();
    if sources.len() == 1 {
        size[sources[0]]
    } else {
        0
    }
}

fn alternating_sum(p: &Poset) -> i64 {
    let weights = vec![1i64; p.n];
    alternating_sum_weighted(p, &weights)
}

fn alternating_sum_weighted(p: &Poset, weights: &[i64]) -> i64 {
    assert_eq!(p.n, weights.len());
    assert!(
        p.strict_edges.len() < 63,
        "too many strict edges for u64 mask"
    );
    let mut total = 0i64;
    for mask in 0..(1u64 << p.strict_edges.len()) {
        let sign = if mask.count_ones() % 2 == 0 { 1 } else { -1 };
        total += sign * m_value_weighted(p, mask, weights);
    }
    total
}

fn strongly_connected_components(graph: &[Vec<usize>]) -> Vec<usize> {
    struct Tarjan<'a> {
        graph: &'a [Vec<usize>],
        next_index: usize,
        stack: Vec<usize>,
        on_stack: Vec<bool>,
        index: Vec<Option<usize>>,
        low: Vec<usize>,
        comp: Vec<usize>,
        comp_count: usize,
    }

    impl<'a> Tarjan<'a> {
        fn visit(&mut self, v: usize) {
            self.index[v] = Some(self.next_index);
            self.low[v] = self.next_index;
            self.next_index += 1;
            self.stack.push(v);
            self.on_stack[v] = true;

            for &w in &self.graph[v] {
                if self.index[w].is_none() {
                    self.visit(w);
                    self.low[v] = self.low[v].min(self.low[w]);
                } else if self.on_stack[w] {
                    self.low[v] = self.low[v].min(self.index[w].unwrap());
                }
            }

            if self.low[v] == self.index[v].unwrap() {
                loop {
                    let w = self.stack.pop().unwrap();
                    self.on_stack[w] = false;
                    self.comp[w] = self.comp_count;
                    if w == v {
                        break;
                    }
                }
                self.comp_count += 1;
            }
        }
    }

    let n = graph.len();
    let mut tarjan = Tarjan {
        graph,
        next_index: 0,
        stack: Vec::new(),
        on_stack: vec![false; n],
        index: vec![None; n],
        low: vec![0; n],
        comp: vec![usize::MAX; n],
        comp_count: 0,
    };

    for v in 0..n {
        if tarjan.index[v].is_none() {
            tarjan.visit(v);
        }
    }
    tarjan.comp
}

fn has_2x2(shape: &Shape) -> bool {
    let cells: HashSet<Cell> = shape.cells.iter().copied().collect();
    for &c in &shape.cells {
        let r = shape.cyl.right(c);
        let u = shape.cyl.up(c);
        let ur = shape.cyl.right(u);
        if cells.contains(&r) && cells.contains(&u) && cells.contains(&ur) {
            return true;
        }
    }
    false
}

fn connected(shape: &Shape) -> bool {
    if shape.cells.is_empty() {
        return false;
    }
    let cells: HashSet<Cell> = shape.cells.iter().copied().collect();
    let start = *shape.cells.iter().next().unwrap();
    let mut seen = HashSet::new();
    let mut stack = vec![start];
    seen.insert(start);
    while let Some(c) = stack.pop() {
        let neigh = [
            shape.cyl.right(c),
            shape.cyl.norm(c.col - 1, c.row),
            shape.cyl.up(c),
            shape.cyl.down(c),
        ];
        for v in neigh {
            if cells.contains(&v) && seen.insert(v) {
                stack.push(v);
            }
        }
    }
    seen.len() == shape.cells.len()
}

fn ordinary_height(word: &[char]) -> i64 {
    word.iter().filter(|&&c| c == 'U').count() as i64
}

fn signed_height(height: i64) -> i64 {
    if height % 2 == 0 {
        1
    } else {
        -1
    }
}

fn words_of_len(n: usize, alphabet: &[char]) -> Vec<Vec<char>> {
    let mut out = Vec::new();
    let total = alphabet.len().pow(n as u32);
    for mut k in 0..total {
        let mut w = Vec::with_capacity(n);
        for _ in 0..n {
            w.push(alphabet[k % alphabet.len()]);
            k /= alphabet.len();
        }
        out.push(w);
    }
    out
}

fn loop_words(x: usize, y: usize) -> Vec<Vec<char>> {
    let n = x + y;
    let mut out = Vec::new();
    for mask in 0usize..(1usize << n) {
        if mask.count_ones() as usize != y {
            continue;
        }
        let mut w = Vec::with_capacity(n);
        for i in 0..n {
            w.push(if ((mask >> i) & 1) == 1 { 'U' } else { 'R' });
        }
        out.push(w);
    }
    out
}

fn word_string(w: &[char]) -> String {
    w.iter().collect()
}

fn scan_ordinary(args: &Args) -> CheckStats {
    let mut stats = CheckStats::default();
    for len in 0..=args.max_ordinary_moves {
        for word in words_of_len(len, &['R', 'U']) {
            let cyl = Cylinder::new(100, 100);
            let shape = shape_from_path(cyl, &word, false, 0);
            if !connected(&shape) || has_2x2(&shape) {
                continue;
            }
            let p = poset(&shape);
            let actual = alternating_sum(&p);
            let expected = signed_height(ordinary_height(&word));
            stats.checked += 1;
            if actual != expected {
                stats.failed += 1;
                println!(
                    "ordinary FAIL word={} cells={} strict={} actual={} expected={}",
                    word_string(&word),
                    shape.cells.len(),
                    p.strict_edges.len(),
                    actual,
                    expected
                );
            } else if args.verbose {
                println!(
                    "ordinary ok word={} cells={} strict={} value={}",
                    word_string(&word),
                    shape.cells.len(),
                    p.strict_edges.len(),
                    actual
                );
            }
        }
    }
    stats
}

fn scan_loops(args: &Args) -> CheckStats {
    let mut stats = CheckStats::default();
    for n in 2..=args.max_loop_len {
        for x in 1..n {
            let y = n - x;
            for word in loop_words(x, y) {
                let cyl = Cylinder::new(x as i32, y as i32);
                let shape = shape_from_path(cyl, &word, true, 0);
                if shape.cells.len() != n || !connected(&shape) || has_2x2(&shape) {
                    continue;
                }
                let p = poset(&shape);
                let width = shape.cells.len() as i64 - p.strict_edges.len() as i64;
                let actual = alternating_sum(&p);
                let expected = signed_height(y as i64) * width;
                stats.checked += 1;
                if actual != expected {
                    stats.failed += 1;
                    println!(
                        "loop FAIL x={} y={} word={} cells={} strict={} width={} actual={} expected={}",
                        x,
                        y,
                        word_string(&word),
                        shape.cells.len(),
                        p.strict_edges.len(),
                        width,
                        actual,
                        expected
                    );
                } else if args.verbose {
                    println!(
                        "loop ok x={} y={} word={} strict={} width={} value={}",
                        x,
                        y,
                        word_string(&word),
                        p.strict_edges.len(),
                        width,
                        actual
                    );
                }
            }
        }
    }
    stats
}

fn scan_residual_ordinary(args: &Args) -> CheckStats {
    let mut stats = CheckStats::default();
    let anchor_weight = 17;
    for len in 0..=args.max_ordinary_moves {
        for word in words_of_len(len, &['R', 'U']) {
            let cyl = Cylinder::new(100, 100);
            let shape = shape_from_path(cyl, &word, false, 0);
            if !connected(&shape) || has_2x2(&shape) {
                continue;
            }
            let (actual, source_count) = anchored_residual_sum(&shape, anchor_weight);
            let expected = -signed_height(ordinary_height(&word));
            stats.checked += 1;
            if actual != expected {
                stats.failed += 1;
                println!(
                    "residual ordinary FAIL word={} cells={} sources={} actual={} expected={}",
                    word_string(&word),
                    shape.cells.len(),
                    source_count,
                    actual,
                    expected
                );
            } else if args.residual_report {
                println!(
                    "residual ordinary ok word={} cells={} sources={} value={}",
                    word_string(&word),
                    shape.cells.len(),
                    source_count,
                    actual
                );
            }
        }
    }
    stats
}

fn scan_residual_loops(args: &Args) -> CheckStats {
    let mut stats = CheckStats::default();
    let anchor_weight = 17;
    for n in 2..=args.max_loop_len {
        for x in 1..n {
            let y = n - x;
            for word in loop_words(x, y) {
                let cyl = Cylinder::new(x as i32, y as i32);
                let shape = shape_from_path(cyl, &word, true, 0);
                if shape.cells.len() != n || !connected(&shape) || has_2x2(&shape) {
                    continue;
                }
                let p = poset(&shape);
                let width = shape.cells.len() as i64 - p.strict_edges.len() as i64;
                let (actual, source_count) = anchored_residual_sum(&shape, anchor_weight);
                let expected = -signed_height(y as i64) * width;
                stats.checked += 1;
                if actual != expected {
                    stats.failed += 1;
                    println!(
                        "residual loop FAIL x={} y={} word={} cells={} sources={} width={} actual={} expected={}",
                        x,
                        y,
                        word_string(&word),
                        shape.cells.len(),
                        source_count,
                        width,
                        actual,
                        expected
                    );
                } else if args.residual_report {
                    println!(
                        "residual loop ok x={} y={} word={} sources={} width={} value={}",
                        x,
                        y,
                        word_string(&word),
                        source_count,
                        width,
                        actual
                    );
                }
            }
        }
    }
    stats
}

fn scan_residual_bad_grid(args: &Args) -> CheckStats {
    let mut stats = CheckStats::default();
    let anchor_weight = 17;
    let side = args.bad_grid.min(3);
    if side == 0 {
        return stats;
    }

    let bits = side * side;
    for mask in 1u64..(1u64 << bits) {
        let shape = shape_from_grid_subset(side, mask);
        let is_bad_base_case = !connected(&shape) || has_2x2(&shape);
        if !is_bad_base_case {
            continue;
        }
        let (actual, source_count) = anchored_residual_sum(&shape, anchor_weight);
        stats.checked += 1;
        if actual != 0 {
            stats.failed += 1;
            println!(
                "residual bad-grid FAIL side={} mask={:#x} cells={} sources={} connected={} 2x2={} actual={} expected=0",
                side,
                mask,
                shape.cells.len(),
                source_count,
                connected(&shape),
                has_2x2(&shape),
                actual
            );
        } else if args.residual_report {
            println!(
                "residual bad-grid ok side={} mask={:#x} cells={} sources={}",
                side,
                mask,
                shape.cells.len(),
                source_count
            );
        }
    }
    stats
}

fn shape_from_grid_subset(side: usize, mask: u64) -> Shape {
    let cyl = Cylinder::new(100, 100);
    let mut cells = BTreeSet::new();
    for row in 0..side {
        for col in 0..side {
            let bit = row * side + col;
            if ((mask >> bit) & 1) == 1 {
                cells.insert(cyl.norm(col as i32, row as i32));
            }
        }
    }
    Shape { cyl, cells }
}

fn scan_bad_grid(args: &Args) -> CheckStats {
    let mut stats = CheckStats::default();
    if args.bad_grid == 0 || args.bad_grid > 4 {
        println!("bad-grid scan skipped: use 1 <= --bad-grid <= 4");
        return stats;
    }

    let bits = args.bad_grid * args.bad_grid;
    for mask in 1u64..(1u64 << bits) {
        let shape = shape_from_grid_subset(args.bad_grid, mask);
        let is_bad_base_case = !connected(&shape) || has_2x2(&shape);
        if !is_bad_base_case {
            continue;
        }
        let p = poset(&shape);
        let actual = alternating_sum(&p);
        stats.checked += 1;
        if actual != 0 {
            stats.failed += 1;
            println!(
                "bad-grid FAIL side={} mask={:#x} cells={} connected={} 2x2={} strict={} actual={} expected=0",
                args.bad_grid,
                mask,
                shape.cells.len(),
                connected(&shape),
                has_2x2(&shape),
                p.strict_edges.len(),
                actual
            );
        }
    }
    stats
}

fn scan_extended_open(args: &Args) -> CheckStats {
    let mut stats = CheckStats::default();
    for n in 2..=args.max_loop_len.min(7) {
        for x in 1..n {
            let y = n - x;
            let cyl = Cylinder::new(x as i32, y as i32);
            for moves in n..=args.max_extended_len {
                if moves + 1 <= n || (moves + 1) % n == 0 {
                    continue;
                }
                for word in words_of_len(moves, &['R', 'U']) {
                    let shape = shape_from_path(cyl.clone(), &word, false, 0);
                    if shape.cells.len() != moves + 1 || !connected(&shape) || has_2x2(&shape) {
                        continue;
                    }
                    let u_count = word.iter().filter(|&&c| c == 'U').count();
                    if u_count < y {
                        continue;
                    }
                    let p = poset(&shape);
                    if p.strict_edges.len() >= 31 {
                        continue;
                    }
                    let actual = alternating_sum(&p);
                    let expected = signed_height(u_count as i64);
                    stats.checked += 1;
                    if actual != expected {
                        stats.failed += 1;
                        println!(
                            "extended-open FAIL x={} y={} word={} cells={} strict={} U={} actual={} expected={}",
                            x,
                            y,
                            word_string(&word),
                            shape.cells.len(),
                            p.strict_edges.len(),
                            u_count,
                            actual,
                            expected
                        );
                    } else if args.extended_report {
                        println!(
                            "extended-open ok x={} y={} word={} cells={} strict={} U={} value={}",
                            x,
                            y,
                            word_string(&word),
                            shape.cells.len(),
                            p.strict_edges.len(),
                            u_count,
                            actual
                        );
                    }
                }
            }
        }
    }
    stats
}

fn scan_extended_pure(args: &Args) -> CheckStats {
    let mut stats = CheckStats::default();
    for n in 2..=args.max_loop_len.min(7) {
        for x in 1..n {
            let y = n - x;
            let cyl = Cylinder::new(x as i32, y as i32);
            for periods in 2..=3 {
                let moves = periods * n;
                if moves > args.max_extended_len {
                    continue;
                }
                for word in words_of_len(moves, &['R', 'U']) {
                    let r_count = word.iter().filter(|&&c| c == 'R').count();
                    let u_count = word.len() - r_count;
                    if r_count != periods * x || u_count != periods * y {
                        continue;
                    }
                    let shape = shape_from_path(cyl.clone(), &word, true, 0);
                    if shape.cells.len() != moves || !connected(&shape) || has_2x2(&shape) {
                        continue;
                    }
                    let p = poset(&shape);
                    if p.strict_edges.len() >= 31 {
                        continue;
                    }
                    let actual = alternating_sum(&p);
                    let width_proxy = shape.cells.len() as i64 - p.strict_edges.len() as i64;
                    let expected_total = signed_height(u_count as i64) * width_proxy;
                    let expected_one_loop = signed_height(y as i64) * width_proxy;
                    stats.checked += 1;
                    if actual != expected_total {
                        stats.failed += 1;
                        println!(
                            "extended-pure FAIL x={} y={} periods={} word={} cells={} strict={} width_proxy={} U={} actual={} expected_total={} expected_one_loop={}",
                            x,
                            y,
                            periods,
                            word_string(&word),
                            shape.cells.len(),
                            p.strict_edges.len(),
                            width_proxy,
                            u_count,
                            actual,
                            expected_total,
                            expected_one_loop
                        );
                    } else if args.extended_report {
                        println!(
                            "extended-pure ok x={} y={} periods={} word={} cells={} strict={} width_proxy={} U={} value={} one_loop_sign_value={}",
                            x,
                            y,
                            periods,
                            word_string(&word),
                            shape.cells.len(),
                            p.strict_edges.len(),
                            width_proxy,
                            u_count,
                            actual,
                            expected_one_loop
                        );
                    }
                }
            }
        }
    }
    stats
}

fn stacked_report() {
    println!("\nstacked-band report (parallel shifted loop ribbons; diagnostic only)");
    for (x, y, word) in [
        (3, 2, "RRURU"),
        (4, 2, "RRURRU"),
        (4, 3, "RRURURU"),
        (5, 3, "RRURRURU"),
    ] {
        let word: Vec<char> = word.chars().collect();
        let cyl = Cylinder::new(x, y);
        for layers in 1..=4 {
            let shape = union_shifted_loops(cyl.clone(), &word, layers);
            let p = poset(&shape);
            if p.strict_edges.len() >= 31 {
                continue;
            }
            let actual = alternating_sum(&p);
            let width_proxy = shape.cells.len() as i64 - p.strict_edges.len() as i64;
            let sign_total_height = signed_height((layers * y) as i64);
            let sign_one_loop = signed_height(y as i64);
            println!(
                "x={x} y={y} word={} layers={layers} cells={} strict={} 2x2={} A={} width_proxy={} A/width_proxy={:?} signs(total={}, one_loop={})",
                word_string(&word),
                shape.cells.len(),
                p.strict_edges.len(),
                has_2x2(&shape),
                actual,
                width_proxy,
                if width_proxy != 0 && actual % width_proxy == 0 {
                    Some(actual / width_proxy)
                } else {
                    None
                },
                sign_total_height,
                sign_one_loop
            );
        }
    }
}

fn picture_shape_report() {
    println!("\ntoggleProofPics encoded-shape report (diagnostic only)");
    let l1 = [(3, 0, 3), (4, 2, 5), (5, 4, 5)];
    let l2 = [(1, 0, 1), (2, 0, 4), (3, 3, 5)];
    let l3 = [(0, 2, 4), (1, 1, 5), (2, 4, 5)];
    let f = [(0, 2, 4)];

    for (x, y) in [(5, 2), (4, 3)] {
        let cyl = Cylinder::new(x, y);
        println!("candidate cylinder C_{{{x},{y}}}");
        for (name, intervals, expected_height) in [
            ("L1", l1.as_slice(), y),
            (
                "L1+L2",
                [l1.as_slice(), l2.as_slice()].concat().as_slice(),
                2 * y,
            ),
            (
                "L1+L2+L3",
                [l1.as_slice(), l2.as_slice(), l3.as_slice()]
                    .concat()
                    .as_slice(),
                3 * y,
            ),
            (
                "L1+L2+L3+F",
                [l1.as_slice(), l2.as_slice(), l3.as_slice(), f.as_slice()]
                    .concat()
                    .as_slice(),
                3 * y,
            ),
        ] {
            let shape = shape_from_intervals(cyl.clone(), intervals);
            let p = poset(&shape);
            if p.strict_edges.len() >= 31 {
                println!(
                    "{name}: skipped; strict edge count {}",
                    p.strict_edges.len()
                );
                continue;
            }
            let actual = alternating_sum(&p);
            println!(
                "{name}: cells={} strict={} connected={} 2x2={} A={} sign(ht={})={}",
                shape.cells.len(),
                p.strict_edges.len(),
                connected(&shape),
                has_2x2(&shape),
                actual,
                expected_height,
                signed_height(expected_height as i64)
            );
        }
    }
}

fn distribution_for_example() {
    let cyl = Cylinder::new(3, 2);
    let word: Vec<char> = "RRURU".chars().collect();
    let shape = shape_from_path(cyl, &word, true, 0);
    let p = poset(&shape);
    let mut dist: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for mask in 0..(1u64 << p.strict_edges.len()) {
        let m = m_value(&p, mask);
        let chosen = mask.count_ones() as usize;
        *dist.entry((chosen, m)).or_default() += 1;
    }
    println!("\nexample distribution for loop x=3 y=2 word=RRURU");
    println!("strict_edges={:?}", p.strict_edges);
    for ((chosen, m), count) in dist {
        println!("  |S|={chosen}, m_S={m}: {count}");
    }
}

fn main() {
    let args = Args::parse();

    let ordinary = scan_ordinary(&args);
    println!(
        "ordinary ribbons checked: {}, failures: {}",
        ordinary.checked, ordinary.failed
    );

    let loops = scan_loops(&args);
    println!(
        "loop ribbons checked: {}, failures: {}",
        loops.checked, loops.failed
    );

    let residual_ordinary = scan_residual_ordinary(&args);
    println!(
        "anchored residual ordinary ribbons checked: {}, failures: {}",
        residual_ordinary.checked, residual_ordinary.failed
    );

    let residual_loops = scan_residual_loops(&args);
    println!(
        "anchored residual loop ribbons checked: {}, failures: {}",
        residual_loops.checked, residual_loops.failed
    );

    let residual_bad_grid = scan_residual_bad_grid(&args);
    println!(
        "anchored residual non-ribbon grid subsets checked: {}, failures: {}",
        residual_bad_grid.checked, residual_bad_grid.failed
    );

    let bad_grid = scan_bad_grid(&args);
    println!(
        "ordinary non-ribbon grid subsets checked: {}, failures: {}",
        bad_grid.checked, bad_grid.failed
    );

    if args.extended_report {
        let extended_open = scan_extended_open(&args);
        println!(
            "exploratory extended-open candidates checked: {}, mismatches: {}",
            extended_open.checked, extended_open.failed
        );

        let extended_pure = scan_extended_pure(&args);
        println!(
            "exploratory extended-pure candidates checked: {}, mismatches: {}",
            extended_pure.checked, extended_pure.failed
        );
    }

    distribution_for_example();

    if args.stacked_report {
        stacked_report();
    }

    if args.picture_report {
        picture_shape_report();
    }

    if ordinary.failed > 0
        || loops.failed > 0
        || residual_ordinary.failed > 0
        || residual_loops.failed > 0
        || residual_bad_grid.failed > 0
        || bad_grid.failed > 0
    {
        std::process::exit(1);
    }
}
