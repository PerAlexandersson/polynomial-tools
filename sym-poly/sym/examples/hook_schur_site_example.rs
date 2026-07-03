use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum HookEntry {
    Unprimed(u32),
    Primed(u32),
}

impl HookEntry {
    fn is_unprimed(self) -> bool {
        matches!(self, HookEntry::Unprimed(_))
    }

    fn is_primed(self) -> bool {
        matches!(self, HookEntry::Primed(_))
    }

    fn index(self) -> u32 {
        match self {
            HookEntry::Unprimed(index) | HookEntry::Primed(index) => index,
        }
    }

    fn monomial_index(self) -> usize {
        match self {
            HookEntry::Unprimed(index) => index as usize - 1,
            HookEntry::Primed(index) => index as usize + 1,
        }
    }

    fn display(self) -> String {
        match self {
            HookEntry::Unprimed(index) => index.to_string(),
            HookEntry::Primed(index) => format!("{index}'"),
        }
    }
}

type HookTableau = Vec<Vec<HookEntry>>;
type Expansion = BTreeMap<Vec<u32>, u32>;

fn is_hook_tableau(tableau: &HookTableau) -> bool {
    for row in tableau {
        for pair in row.windows(2) {
            if pair[0] > pair[1] {
                return false;
            }
            if pair[0].is_primed() && pair[1].is_primed() && pair[0].index() >= pair[1].index() {
                return false;
            }
        }
    }

    for row in 0..tableau.len().saturating_sub(1) {
        for col in 0..tableau[row + 1].len() {
            let top = tableau[row][col];
            let bottom = tableau[row + 1][col];
            if top > bottom {
                return false;
            }
            if top.is_unprimed() && bottom.is_unprimed() && top.index() >= bottom.index() {
                return false;
            }
        }
    }

    true
}

fn enumerate_hook_tableaux(shape: &[usize], alphabet: &[HookEntry]) -> Vec<HookTableau> {
    fn backtrack(
        shape: &[usize],
        alphabet: &[HookEntry],
        row: usize,
        col: usize,
        current: &mut HookTableau,
        results: &mut Vec<HookTableau>,
    ) {
        if row == shape.len() {
            if is_hook_tableau(current) {
                results.push(current.clone());
            }
            return;
        }

        let (next_row, next_col) = if col + 1 == shape[row] {
            (row + 1, 0)
        } else {
            (row, col + 1)
        };

        for &entry in alphabet {
            current[row][col] = entry;
            backtrack(shape, alphabet, next_row, next_col, current, results);
        }
    }

    let mut current = shape
        .iter()
        .map(|&length| vec![alphabet[0]; length])
        .collect::<HookTableau>();
    let mut results = Vec::new();
    backtrack(shape, alphabet, 0, 0, &mut current, &mut results);
    results
}

fn weight(tableau: &HookTableau) -> Vec<u32> {
    let mut weight = vec![0; 3];
    for row in tableau {
        for &entry in row {
            weight[entry.monomial_index()] += 1;
        }
    }
    weight
}

fn hook_schur_31_expansion() -> Expansion {
    let alphabet = [
        HookEntry::Unprimed(1),
        HookEntry::Unprimed(2),
        HookEntry::Primed(1),
    ];
    let mut expansion = Expansion::new();

    for tableau in enumerate_hook_tableaux(&[3, 1], &alphabet) {
        *expansion.entry(weight(&tableau)).or_insert(0) += 1;
    }

    expansion
}

fn format_tableau(tableau: &HookTableau) -> String {
    tableau
        .iter()
        .map(|row| {
            row.iter()
                .map(|entry| entry.display())
                .collect::<Vec<_>>()
                .join(" & ")
        })
        .collect::<Vec<_>>()
        .join(" \\\\ ")
}

fn format_monomial(weight: &[u32], coefficient: u32) -> String {
    let variables = ["x_1", "x_2", "y_1"];
    let monomial = variables
        .iter()
        .zip(weight)
        .filter(|(_, &exponent)| exponent > 0)
        .map(|(variable, &exponent)| {
            if exponent == 1 {
                variable.to_string()
            } else {
                format!("{variable}^{exponent}")
            }
        })
        .collect::<Vec<_>>()
        .join("");

    if coefficient == 1 {
        monomial
    } else {
        format!("{coefficient}{monomial}")
    }
}

fn main() {
    let expansion = hook_schur_31_expansion();
    let expected = BTreeMap::from([
        (vec![0, 2, 2], 1),
        (vec![0, 3, 1], 1),
        (vec![1, 1, 2], 1),
        (vec![1, 2, 1], 2),
        (vec![1, 3, 0], 1),
        (vec![2, 0, 2], 1),
        (vec![2, 1, 1], 2),
        (vec![2, 2, 0], 1),
        (vec![3, 0, 1], 1),
        (vec![3, 1, 0], 1),
    ]);
    assert_eq!(expansion, expected);

    let first = vec![
        vec![
            HookEntry::Unprimed(1),
            HookEntry::Unprimed(1),
            HookEntry::Primed(1),
        ],
        vec![HookEntry::Unprimed(2)],
    ];
    let second = vec![
        vec![
            HookEntry::Unprimed(1),
            HookEntry::Unprimed(1),
            HookEntry::Unprimed(2),
        ],
        vec![HookEntry::Primed(1)],
    ];
    assert!(is_hook_tableau(&first));
    assert!(is_hook_tableau(&second));
    assert_eq!(weight(&first), vec![2, 1, 1]);
    assert_eq!(weight(&second), vec![2, 1, 1]);

    println!("hook tableaux with weight x_1^2 x_2 y_1:");
    println!("{}", format_tableau(&first));
    println!("{}", format_tableau(&second));
    println!("s_(3,1)(x_1,x_2 / y_1) =");
    println!(
        "{}",
        expansion
            .iter()
            .map(|(weight, &coefficient)| format_monomial(weight, coefficient))
            .collect::<Vec<_>>()
            .join(" + ")
    );
}
