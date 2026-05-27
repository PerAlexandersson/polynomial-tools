use std::collections::BTreeMap;

use flagged_lorentzian::{
    crystal_f, enumerate_tableaux, is_gt_array, RowFlaggedSkewShape, SkewGtPattern, SkewShape,
    TableauRecord,
};

type GtRows = Vec<Vec<u32>>;

#[derive(Debug, Clone)]
struct Args {
    lambda: Vec<u32>,
    mu: Vec<u32>,
    row_flags: Option<Vec<u32>>,
    alphabet: usize,
    lower_label: usize,
    beta: Vec<u32>,
    left_pos: usize,
    right_pos: usize,
}

impl Default for Args {
    fn default() -> Self {
        Self {
            lambda: vec![4, 2],
            mu: vec![2],
            row_flags: None,
            alphabet: 5,
            lower_label: 4,
            beta: vec![0, 1, 1, 0, 0],
            left_pos: 1,
            right_pos: 2,
        }
    }
}

#[derive(Debug, Clone)]
struct Split {
    pos: usize,
    complement_pos: usize,
    content: Vec<u32>,
    rows: GtRows,
}

fn main() {
    let args = match parse_args() {
        Ok(args) => args,
        Err(message) => {
            eprintln!("{message}");
            eprintln!();
            eprintln!("{}", help_text());
            std::process::exit(2);
        }
    };

    let shape = SkewShape::from_parts(args.lambda.clone(), args.mu.clone());
    let flagged = match &args.row_flags {
        Some(row_flags) => {
            RowFlaggedSkewShape::new(shape.clone(), row_flags.clone(), args.alphabet)
        }
        None => RowFlaggedSkewShape::ordinary(shape.clone(), args.alphabet),
    };
    let tableaux = enumerate_tableaux(&flagged, None).unwrap();
    let gt_rows: Vec<_> = tableaux
        .iter()
        .map(|tableau| {
            SkewGtPattern::from_tableau(&shape, &tableau.values, args.alphabet)
                .rows()
                .to_vec()
        })
        .collect();
    let gt_to_pos: BTreeMap<_, _> = gt_rows
        .iter()
        .enumerate()
        .map(|(pos, rows)| (rows.clone(), pos))
        .collect();
    let values_to_pos: BTreeMap<_, _> = tableaux
        .iter()
        .enumerate()
        .map(|(pos, tableau)| (tableau.values.clone(), pos))
        .collect();

    let lower = args.lower_label - 1;
    let upper = lower + 1;
    let ii = add_units(&args.beta, lower, 2);
    let mixed = add_units(&add_units(&args.beta, lower, 1), upper, 1);
    let jj = add_units(&args.beta, upper, 2);
    let ii_indices = indices_with_content(&tableaux, &ii);
    let mixed_indices = indices_with_content(&tableaux, &mixed);
    let jj_indices = indices_with_content(&tableaux, &jj);

    let left_idx = ii_indices[args.left_pos];
    let right_idx = jj_indices[args.right_pos];
    let z = add_gt_rows(&gt_rows[left_idx], &gt_rows[right_idx]);

    println!("two-row split poset dump");
    println!("outer={:?}", shape.outer().parts());
    println!("inner={:?}", shape.inner().parts());
    println!("flags={:?}", flagged.row_flags());
    println!(
        "beta={:?}, active labels={},{}",
        args.beta,
        lower + 1,
        upper + 1
    );
    println!(
        "source=(A{},C{}) = {} / {}",
        args.left_pos,
        args.right_pos,
        format_tableau(&shape, &tableaux[left_idx].values),
        format_tableau(&shape, &tableaux[right_idx].values)
    );
    println!("Z={}", summarize_gt(&z));
    println!();

    let splits = split_poset(&z, &gt_rows, &gt_to_pos, &tableaux);
    println!("splits in P_Z: {}", splits.len());
    for (rank, split) in splits.iter().enumerate() {
        let complement = &tableaux[split.complement_pos];
        println!(
            "T{rank}=all{} comp{} content={:?} active_sum={} y={} :: {} / {}",
            split.pos,
            split.complement_pos,
            split.content,
            split.rows[args.lower_label][0] + split.rows[args.lower_label][1],
            bottom_word(&split.rows),
            format_tableau(&shape, &tableaux[split.pos].values),
            format_tableau(&shape, &complement.values)
        );
    }
    println!();

    let split_pos_to_rank: BTreeMap<_, _> = splits
        .iter()
        .enumerate()
        .map(|(rank, split)| (split.pos, rank))
        .collect();
    println!("source slice beta+2e_i:");
    for (pos, &idx) in ii_indices.iter().enumerate() {
        if let Some(&rank) = split_pos_to_rank.get(&idx) {
            println!("  A{pos}: T{rank}");
        }
    }
    println!("mixed slice beta+e_i+e_j:");
    for (pos, &idx) in mixed_indices.iter().enumerate() {
        if let Some(&rank) = split_pos_to_rank.get(&idx) {
            println!("  M{pos}: T{rank}");
        }
    }
    println!();

    println!("fixed-outside-content beta-line in P_Z:");
    let beta_line: Vec<_> = splits
        .iter()
        .enumerate()
        .filter(|(_, split)| outside_content_matches(&split.content, &args.beta, lower, upper))
        .collect();
    for (rank, split) in &beta_line {
        println!(
            "  T{rank}: active=({},{}), rows={}, y={}, gap={} :: {} / {}",
            split.content[lower],
            split.content[upper],
            format_rows(&split.rows),
            bottom_word(&split.rows),
            gap_word(&split.rows),
            format_tableau(&shape, &tableaux[split.pos].values),
            format_tableau(&shape, &tableaux[split.complement_pos].values)
        );
    }
    println!();

    println!("beta-line grouped by lower prefix and active row:");
    let mut by_lower_prefix = BTreeMap::<Vec<Vec<u32>>, Vec<String>>::new();
    for (rank, split) in &beta_line {
        by_lower_prefix
            .entry(split.rows[..args.lower_label].to_vec())
            .or_default()
            .push(format!(
                "T{rank}:r{}={:?},active=({},{}),x={}",
                args.lower_label,
                split.rows[args.lower_label],
                split.content[lower],
                split.content[upper],
                split.rows[args.lower_label][1]
            ));
    }
    for (lower_prefix, mut items) in by_lower_prefix {
        items.sort();
        let prefix = lower_prefix
            .iter()
            .enumerate()
            .map(|(level, row)| format!("r{level}={row:?}"))
            .collect::<Vec<_>>()
            .join(";");
        println!("  {prefix}: {}", items.join("; "));
    }
    println!();

    println!("carrier edges on beta-line, from beta+2e_i to beta+e_i+e_j:");
    for (source_rank, source) in beta_line
        .iter()
        .copied()
        .filter(|(_, split)| split.content.as_slice() == ii.as_slice())
    {
        let edges: Vec<_> = beta_line
            .iter()
            .copied()
            .filter(|(_, split)| split.content.as_slice() == mixed.as_slice())
            .filter_map(|(target_rank, target)| {
                carrier_delta(&source.rows, &target.rows, args.lower_label)
                    .map(|carrier| (target_rank, carrier))
            })
            .collect();
        let formatted = edges
            .iter()
            .map(|(target_rank, carrier)| format!("T{target_rank}:{}", format_carrier(carrier)))
            .collect::<Vec<_>>()
            .join(", ");
        println!("  T{source_rank} -> {{{formatted}}}");
        let standard_target = crystal_f(
            &shape,
            &tableaux[source.pos].values,
            args.lower_label as u32,
        )
        .and_then(|values| values_to_pos.get(&values).copied())
        .and_then(|pos| split_pos_to_rank.get(&pos).copied());
        match standard_target {
            Some(target_rank) => println!("    ordinary tableau f_i gives T{target_rank}"),
            None => println!("    ordinary tableau f_i leaves the beta-line"),
        }
    }
    println!();

    println!("all adjacent-rank carrier edges on beta-line:");
    let mut undirected_edges = Vec::<(usize, usize)>::new();
    for (source_rank, source) in beta_line.iter().copied() {
        if source.content[lower] == 0 {
            continue;
        }
        let target_content_i = source.content[lower] - 1;
        let edges: Vec<_> = beta_line
            .iter()
            .copied()
            .filter(|(_, target)| target.content[lower] == target_content_i)
            .filter_map(|(target_rank, target)| {
                carrier_delta(&source.rows, &target.rows, args.lower_label)
                    .map(|carrier| (target_rank, carrier))
            })
            .collect();
        if edges.is_empty() {
            continue;
        }
        let formatted = edges
            .iter()
            .map(|(target_rank, carrier)| format!("T{target_rank}:{}", format_carrier(carrier)))
            .collect::<Vec<_>>()
            .join(", ");
        println!(
            "  T{source_rank} rank {} -> rank {} {{{formatted}}}",
            source.content[lower], target_content_i
        );
        for (target_rank, _) in edges {
            undirected_edges.push((source_rank, target_rank));
        }
    }
    println!();

    println!("carrier graph connected components on beta-line:");
    for component in carrier_components(&beta_line, &undirected_edges, lower) {
        println!("  {}", component);
    }
}

fn split_poset(
    z: &[Vec<u32>],
    gt_rows: &[GtRows],
    gt_to_pos: &BTreeMap<GtRows, usize>,
    tableaux: &[TableauRecord],
) -> Vec<Split> {
    let mut splits = Vec::new();
    for (pos, rows) in gt_rows.iter().enumerate() {
        let Some(complement) = subtract_gt_rows(z, rows) else {
            continue;
        };
        if !is_gt_array(&complement) {
            continue;
        }
        let Some(&complement_pos) = gt_to_pos.get(&complement) else {
            continue;
        };
        splits.push(Split {
            pos,
            complement_pos,
            content: tableaux[pos].content.clone(),
            rows: rows.clone(),
        });
    }
    splits.sort_by_key(|split| {
        (
            split.rows.clone(),
            split.content.clone(),
            split.pos,
            split.complement_pos,
        )
    });
    splits
}

fn add_gt_rows(left: &[Vec<u32>], right: &[Vec<u32>]) -> GtRows {
    left.iter()
        .zip(right)
        .map(|(left, right)| left.iter().zip(right).map(|(&a, &b)| a + b).collect())
        .collect()
}

fn subtract_gt_rows(left: &[Vec<u32>], right: &[Vec<u32>]) -> Option<GtRows> {
    left.iter()
        .zip(right)
        .map(|(left, right)| {
            left.iter()
                .zip(right)
                .map(|(&a, &b)| a.checked_sub(b))
                .collect::<Option<Vec<_>>>()
        })
        .collect()
}

fn bottom_word(rows: &[Vec<u32>]) -> String {
    rows.iter()
        .map(|row| row[1].to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn gap_word(rows: &[Vec<u32>]) -> String {
    rows.iter()
        .map(|row| (row[0] as i64 - row[1] as i64).to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn format_rows(rows: &[Vec<u32>]) -> String {
    rows.iter()
        .enumerate()
        .map(|(level, row)| format!("r{level}={row:?}"))
        .collect::<Vec<_>>()
        .join(";")
}

fn summarize_gt(rows: &[Vec<u32>]) -> String {
    rows.iter()
        .enumerate()
        .filter(|(_, row)| row.iter().any(|&entry| entry != 0))
        .map(|(level, row)| format!("r{level}={row:?}"))
        .collect::<Vec<_>>()
        .join(";")
}

fn indices_with_content(tableaux: &[TableauRecord], content: &[u32]) -> Vec<usize> {
    tableaux
        .iter()
        .enumerate()
        .filter_map(|(idx, tableau)| (tableau.content.as_slice() == content).then_some(idx))
        .collect()
}

fn add_units(content: &[u32], index: usize, amount: u32) -> Vec<u32> {
    let mut out = content.to_vec();
    out[index] += amount;
    out
}

fn outside_content_matches(content: &[u32], beta: &[u32], lower: usize, upper: usize) -> bool {
    content
        .iter()
        .enumerate()
        .all(|(idx, &count)| idx == lower || idx == upper || beta.get(idx).copied() == Some(count))
}

fn carrier_delta(
    source: &[Vec<u32>],
    target: &[Vec<u32>],
    lower_label: usize,
) -> Option<Vec<[i32; 2]>> {
    if source.len() != target.len() || source.len() <= lower_label {
        return None;
    }
    let d: Vec<[i32; 2]> = source
        .iter()
        .zip(target)
        .map(|(source_row, target_row)| {
            [
                source_row[0] as i32 - target_row[0] as i32,
                source_row[1] as i32 - target_row[1] as i32,
            ]
        })
        .collect();
    if d.iter().skip(lower_label + 1).any(|&row| row != [0, 0]) {
        return None;
    }
    if d[lower_label] != [1, 0] && d[lower_label] != [0, 1] {
        return None;
    }
    let start = d.iter().position(|&row| row != [0, 0])?;
    if start > lower_label {
        return None;
    }
    if d.iter().take(start).any(|&row| row != [0, 0]) {
        return None;
    }
    if d.iter()
        .take(lower_label)
        .skip(start)
        .any(|&row| row != [1, -1] && row != [-1, 1])
    {
        return None;
    }
    Some(d)
}

fn format_carrier(carrier: &[[i32; 2]]) -> String {
    carrier
        .iter()
        .enumerate()
        .filter(|(_, delta)| **delta != [0, 0])
        .map(|(level, delta)| format!("r{level}{delta:?}"))
        .collect::<Vec<_>>()
        .join("/")
}

fn carrier_components(
    beta_line: &[(usize, &Split)],
    edges: &[(usize, usize)],
    active_index: usize,
) -> Vec<String> {
    let rank_to_line_pos: BTreeMap<_, _> = beta_line
        .iter()
        .enumerate()
        .map(|(line_pos, (rank, _))| (*rank, line_pos))
        .collect();
    let mut graph = vec![Vec::new(); beta_line.len()];
    for &(left_rank, right_rank) in edges {
        let Some(&left) = rank_to_line_pos.get(&left_rank) else {
            continue;
        };
        let Some(&right) = rank_to_line_pos.get(&right_rank) else {
            continue;
        };
        graph[left].push(right);
        graph[right].push(left);
    }

    let mut seen = vec![false; beta_line.len()];
    let mut out = Vec::new();
    for start in 0..beta_line.len() {
        if seen[start] {
            continue;
        }
        let mut stack = vec![start];
        seen[start] = true;
        let mut vertices = Vec::new();
        while let Some(vertex) = stack.pop() {
            vertices.push(vertex);
            for &neighbor in &graph[vertex] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        vertices.sort_unstable_by_key(|&line_pos| beta_line[line_pos].0);
        let mut profile = BTreeMap::<u32, Vec<usize>>::new();
        for line_pos in vertices {
            let (rank, split) = beta_line[line_pos];
            profile
                .entry(split.content[active_index])
                .or_default()
                .push(rank);
        }
        out.push(
            profile
                .into_iter()
                .map(|(rank, vertices)| format!("rank {rank}: T{vertices:?}"))
                .collect::<Vec<_>>()
                .join("; "),
        );
    }
    out.sort();
    out
}

fn format_tableau(shape: &SkewShape, values: &[u32]) -> String {
    let rows = shape
        .cells()
        .iter()
        .map(|cell| cell.row)
        .max()
        .map_or(0, |row| row + 1);
    let cols = shape
        .cells()
        .iter()
        .map(|cell| cell.col)
        .max()
        .map_or(0, |col| col + 1);

    let mut out = Vec::new();
    for row in 0..rows {
        let mut row_values = Vec::new();
        for col in 0..cols {
            let value = shape
                .cell_index(flagged_lorentzian::Cell { row, col })
                .map(|idx| values[idx].to_string())
                .unwrap_or_else(|| ".".to_string());
            row_values.push(value);
        }
        out.push(row_values.join(" "));
    }
    format!("[{}]", out.join(" / "))
}

fn parse_args() -> Result<Args, String> {
    let mut args = Args::default();
    let mut iter = std::env::args().skip(1);
    while let Some(flag) = iter.next() {
        match flag.as_str() {
            "--help" | "-h" => {
                println!("{}", help_text());
                std::process::exit(0);
            }
            "--lambda" => args.lambda = parse_u32_vec(&take_value(&mut iter, &flag)?)?,
            "--mu" => args.mu = parse_u32_vec(&take_value(&mut iter, &flag)?)?,
            "--row-flags" => args.row_flags = Some(parse_u32_vec(&take_value(&mut iter, &flag)?)?),
            "--alphabet" => {
                args.alphabet = take_value(&mut iter, &flag)?
                    .parse()
                    .map_err(|err| format!("invalid --alphabet: {err}"))?
            }
            "--lower-label" => {
                args.lower_label = take_value(&mut iter, &flag)?
                    .parse()
                    .map_err(|err| format!("invalid --lower-label: {err}"))?
            }
            "--beta" => args.beta = parse_u32_vec(&take_value(&mut iter, &flag)?)?,
            "--left-pos" => {
                args.left_pos = take_value(&mut iter, &flag)?
                    .parse()
                    .map_err(|err| format!("invalid --left-pos: {err}"))?
            }
            "--right-pos" => {
                args.right_pos = take_value(&mut iter, &flag)?
                    .parse()
                    .map_err(|err| format!("invalid --right-pos: {err}"))?
            }
            other => return Err(format!("unknown argument `{other}`")),
        }
    }
    if args.lower_label == 0 || args.lower_label >= args.alphabet {
        return Err("--lower-label must be between 1 and alphabet-1".to_string());
    }
    if args.beta.len() != args.alphabet {
        return Err("--beta must have one entry per alphabet letter".to_string());
    }
    Ok(args)
}

fn take_value(iter: &mut impl Iterator<Item = String>, flag: &str) -> Result<String, String> {
    iter.next()
        .ok_or_else(|| format!("missing value after {flag}"))
}

fn parse_u32_vec(input: &str) -> Result<Vec<u32>, String> {
    if input.trim().is_empty() {
        return Ok(Vec::new());
    }
    input
        .split(',')
        .map(|part| {
            part.trim()
                .parse::<u32>()
                .map_err(|err| format!("invalid integer `{part}`: {err}"))
        })
        .collect()
}

fn help_text() -> &'static str {
    "Dump the fixed pair-sum split poset P_Z for one two-row source pair.

USAGE:
  two_row_split_poset_dump [OPTIONS]

OPTIONS:
  --lambda PARTS      Outer partition. Default: 4,2.
  --mu PARTS          Inner partition. Default: 2.
  --row-flags PARTS   Row upper flags. Omit for ordinary skew tableaux.
  --alphabet N        Alphabet size. Default: 5.
  --lower-label I     Active labels are I,I+1. Default: 4.
  --beta CONTENT      Base content. Default: 0,1,1,0,0.
  --left-pos N        Source A position in ii fiber. Default: 1.
  --right-pos N       Source C position in jj fiber. Default: 2.
  --help              Print this help."
}
