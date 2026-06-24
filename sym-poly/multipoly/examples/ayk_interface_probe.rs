use sym_poly_multipoly::{
    canonical_labeling, diagram_weight, is_yamanouchi, kohnert_diagrams, rothe_diagram, Diagram,
    Labeling,
};

type ColumnType = Vec<usize>;
type LabeledState = Vec<(usize, usize)>;

fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    let Some(perm_arg) = args.get(1) else {
        panic!(
            "usage: ayk_interface_probe <comma-separated permutation> [dilation] [max diagrams]"
        );
    };
    let perm = parse_perm(perm_arg);
    let dilation = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4);
    let max_diagrams = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(500_000);

    let stretched = stretch_perm(&perm, dilation);
    let initial = rothe_diagram(&stretched);
    let blocks = repeated_blocks(&initial);
    let diagrams = match kohnert_diagrams(&initial, max_diagrams) {
        Ok(diagrams) => diagrams,
        Err(err) => {
            println!("{err}");
            return;
        }
    };

    println!(
        "u={perm:?}, N={dilation}, N*u={stretched:?}, cells={}, KD={}",
        initial.len(),
        diagrams.len()
    );
    println!("repeated source-column blocks:");
    for (start, len, column_type) in &blocks {
        if *len >= 2 {
            println!(
                "  columns {start}..{} type {:?}",
                start + len - 1,
                column_type
            );
        }
    }

    let mut ayk = Vec::new();
    for diagram in diagrams {
        if !is_yamanouchi(&initial, &diagram) {
            continue;
        }
        let Some(labeling) = canonical_labeling(&initial, &diagram) else {
            continue;
        };
        ayk.push((diagram_weight(&diagram), labeling));
    }
    ayk.sort_by(|a, b| a.0.cmp(&b.0));
    println!("AYK count: {}", ayk.len());

    for (idx, (weight, labeling)) in ayk.iter().enumerate() {
        println!("Y{} wt={:?}", idx + 1, weight);
        for (start, len, column_type) in &blocks {
            if *len < 2 {
                continue;
            }
            let states = block_states(labeling, *start, *len);
            let interfaces = interfaces(&states);
            println!(
                "  block columns {start}..{} type {:?}: {}",
                start + len - 1,
                column_type,
                format_states(&states)
            );
            if interfaces.is_empty() {
                println!("    interfaces: none");
            } else {
                println!("    interfaces after local block columns {:?}", interfaces);
            }
        }
    }
}

fn parse_perm(text: &str) -> Vec<usize> {
    text.split(',')
        .map(|s| {
            s.parse::<usize>()
                .expect("permutation entries are integers")
        })
        .collect()
}

fn stretch_perm(perm: &[usize], dilation: usize) -> Vec<usize> {
    let code = lehmer_code(perm);
    let needed_len = code
        .iter()
        .enumerate()
        .filter(|(_, c)| **c > 0)
        .map(|(i, c)| i + 1 + dilation * c)
        .max()
        .unwrap_or(perm.len())
        .max(perm.len());
    let mut stretched_code = vec![0; needed_len];
    for (i, c) in code.iter().enumerate() {
        stretched_code[i] = dilation * c;
    }
    from_lehmer_code(&stretched_code)
}

fn lehmer_code(perm: &[usize]) -> Vec<usize> {
    (0..perm.len())
        .map(|i| perm[i + 1..].iter().filter(|&&v| v < perm[i]).count())
        .collect()
}

fn from_lehmer_code(code: &[usize]) -> Vec<usize> {
    let n = code.len();
    let mut available = (1..=n).collect::<Vec<_>>();
    let mut result = Vec::with_capacity(n);
    for &c in code {
        result.push(available.remove(c));
    }
    result
}

fn repeated_blocks(diagram: &Diagram) -> Vec<(usize, usize, ColumnType)> {
    let max_col = diagram.iter().map(|cell| cell.col).max().unwrap_or(0);
    let types = (1..=max_col)
        .map(|col| {
            let mut rows = diagram
                .iter()
                .filter_map(|cell| {
                    if cell.col == col {
                        Some(cell.row)
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();
            rows.sort();
            rows
        })
        .collect::<Vec<_>>();

    let mut blocks = Vec::new();
    let mut start = 1usize;
    let mut current = types.first().cloned().unwrap_or_default();
    for (idx, column_type) in types.iter().enumerate().skip(1) {
        let col = idx + 1;
        if *column_type != current {
            blocks.push((start, col - start, current));
            start = col;
            current = column_type.clone();
        }
    }
    if max_col > 0 {
        blocks.push((start, max_col + 1 - start, current));
    }
    blocks
}

fn block_states(labeling: &Labeling, block_start: usize, block_len: usize) -> Vec<LabeledState> {
    (block_start..block_start + block_len)
        .map(|col| {
            let mut state = labeling
                .iter()
                .filter(|(cell, _)| cell.col == col)
                .map(|(cell, label)| (cell.row, *label))
                .collect::<Vec<_>>();
            state.sort();
            state
        })
        .collect()
}

fn interfaces(states: &[LabeledState]) -> Vec<usize> {
    states
        .windows(2)
        .enumerate()
        .filter_map(|(idx, pair)| {
            if pair[0] == pair[1] {
                None
            } else {
                Some(idx + 1)
            }
        })
        .collect()
}

fn format_states(states: &[LabeledState]) -> String {
    states
        .iter()
        .map(|state| {
            let entries = state
                .iter()
                .map(|(row, label)| format!("({row},{label})"))
                .collect::<Vec<_>>()
                .join(",");
            format!("{{{entries}}}")
        })
        .collect::<Vec<_>>()
        .join(" ")
}
