use sym_poly_sym::{hook_schur_expansion, HookEntry, HookSchurExpansion, HookTableau, HookWeight};

fn format_monomial(weight: &HookWeight, coefficient: u64) -> String {
    let variables = ["x_1", "x_2", "y_1"];
    let exponents = weight
        .x
        .iter()
        .chain(weight.y.iter())
        .copied()
        .collect::<Vec<_>>();
    let monomial = variables
        .iter()
        .zip(exponents)
        .filter(|(_, exponent)| *exponent > 0)
        .map(|(variable, exponent)| {
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
    let expansion = hook_schur_expansion(&[3, 1], 2, 1);
    let expected = HookSchurExpansion::from([
        (
            HookWeight {
                x: vec![0, 2],
                y: vec![2],
            },
            1,
        ),
        (
            HookWeight {
                x: vec![0, 3],
                y: vec![1],
            },
            1,
        ),
        (
            HookWeight {
                x: vec![1, 1],
                y: vec![2],
            },
            1,
        ),
        (
            HookWeight {
                x: vec![1, 2],
                y: vec![1],
            },
            2,
        ),
        (
            HookWeight {
                x: vec![1, 3],
                y: vec![0],
            },
            1,
        ),
        (
            HookWeight {
                x: vec![2, 0],
                y: vec![2],
            },
            1,
        ),
        (
            HookWeight {
                x: vec![2, 1],
                y: vec![1],
            },
            2,
        ),
        (
            HookWeight {
                x: vec![2, 2],
                y: vec![0],
            },
            1,
        ),
        (
            HookWeight {
                x: vec![3, 0],
                y: vec![1],
            },
            1,
        ),
        (
            HookWeight {
                x: vec![3, 1],
                y: vec![0],
            },
            1,
        ),
    ]);
    assert_eq!(expansion, expected);

    let first = HookTableau::new(vec![
        vec![
            HookEntry::Unprimed(1),
            HookEntry::Unprimed(1),
            HookEntry::Primed(1),
        ],
        vec![HookEntry::Unprimed(2)],
    ]);
    let second = HookTableau::new(vec![
        vec![
            HookEntry::Unprimed(1),
            HookEntry::Unprimed(1),
            HookEntry::Unprimed(2),
        ],
        vec![HookEntry::Primed(1)],
    ]);
    assert!(first.is_hook_tableau());
    assert!(second.is_hook_tableau());
    assert_eq!(
        first.weight(2, 1),
        HookWeight {
            x: vec![2, 1],
            y: vec![1]
        }
    );
    assert_eq!(
        second.weight(2, 1),
        HookWeight {
            x: vec![2, 1],
            y: vec![1]
        }
    );

    println!("hook tableaux with weight x_1^2 x_2 y_1:");
    println!("{}", first.display_rows());
    println!("{}", second.display_rows());
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
